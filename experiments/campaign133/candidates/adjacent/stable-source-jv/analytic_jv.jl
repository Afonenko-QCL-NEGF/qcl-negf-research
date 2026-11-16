# WRITE-ONLY prototype: inclusion never builds a model, reads H5 or runs maps.
# Caller supplies pinned EquilibriumGammaPilot GP and StreamingSpectral SS.
module StableSourceJv
using QCLNEGF, LinearAlgebra
const N=QCLNEGF.QCLNumerics
const Field=Array{ComplexF64,4}
mutable struct Limits
    started::Float64
    seconds::Float64
    pi::Int
    jv::Int
    pi_completed::Int
    jv_completed::Int
end
Limits(;seconds=600.0)=Limits(time(),seconds,0,0,0,0)
function entered!(b::Limits,kind)
    time()-b.started<b.seconds || error("diagnostic wall cap")
    if kind===:pi
        b.pi<5 || error("five fresh-Phi cap"); b.pi+=1
    else
        b.jv<4 || error("four analytic-Jv cap"); b.jv+=1
    end
end
function policy(c)
    p=c.problem; a=c.options.algorithms
    a.embedding===:finite_chain && a.embedding_periods==4 || error("finite4 only")
    a.hilbert===:product_integration && a.energy_shift===:sparse_plan || error("PI policy")
    iszero(c.cavity.drop) && iszero(N._scaled_physics(p.physical,p.scales).Eᵖ) || error("V0 only")
    p.models.lo_population===:thermal || error("fixed thermal Bose only")
    c.cache.product_hilbert_operator isa N._ProductionPIFFTHilbertPlan || error("canonical uniform PI")
    for rows in c.cavity.remaps, e in eachindex(rows)
        all(pair->first(pair)==e || iszero(last(pair)),rows[e]) || error("nonidentity V0 remap")
        sum(last(pair) for pair in rows[e] if first(pair)==e)==1 || error("nonunit V0 remap")
    end
    return nothing
end
"One mathematical raw Phi; streaming primal already independently qualified."
function preserve_raw_before_guard(GP,raw,stats,budget::Limits,persist_raw)
    budget.pi_completed+=1 # completed raw contraction, not a cone/return/admission claim
    persist_raw===nothing || persist_raw(raw,stats)
    GP.cone(raw)
    return nothing
end
function fresh_phi(GP,SS,c,gamma::Field,budget::Limits;keep_base=false,persist_raw=nothing)
    policy(c); GP.cone(gamma); entered!(budget,:pi)
    s=SS.stream_spectrum(GP,c,gamma)
    GP.cone(s.A;require_nonzero=false)
    fit=GP.fit_mu(c.problem.grids.ε,s.q,c.target,c.kBT)
    green,pairerror=GP.fd_green(s.GR,s.A,s.condition,s.scale,c.problem.grids.ε,fit.mu,c.kBT)
    raw,channels,sizes=GP.phi(c,green);empty!(channels)
    stats=(mu=fit,pairerror=pairerror,channel_norms=sizes,spectral=s.ledger)
    preserve_raw_before_guard(GP,raw,stats,budget,persist_raw)
    keep_base || return raw,stats
    R,_=GP.local_retarded(c,gamma) # same Herm(PI Gamma)-i Gamma/2; extra linear PI only
    return (;gamma,R,GR=s.GR,A=s.A,condition=s.condition,scale=s.scale,raw,stats)
end

# All matrices and inverse buffers are owned by one logical worker.
struct Work
    m::Vector{Matrix{ComplexF64}}
    inverse::N._ProductionInverseWorkspace
end
Work(n)=Work([zeros(ComplexF64,n,n) for _ in 1:50],N._ProductionInverseWorkspace(n))
function triple!(out,A,B,C,tmp)
    mul!(tmp,A,B); mul!(out,tmp,C); return out
end
function side!(w::Work,E,h,R,gamma,dR,v,C,outputs)
    g,L,dg,dL,gn,Ln,dgn,dLn,D,Lin,dLin,V,tmp,tmp2,term,product=ntuple(i->w.m[i],Val(16))
    fill!(g,0);fill!(L,0);fill!(dg,0);fill!(dL,0)
    for depth in 1:4
        N._cavity_sandwich!(product,C,g,tmp)
        @. D=-h-R-product
        for a in axes(D,1); D[a,a]+=E; end
        copyto!(gn,D);N._cavity_inverse!(gn,w.inverse)
        N._cavity_sandwich!(V,C,dg,tmp); @. V+=dR
        triple!(dgn,gn,V,gn,tmp)
        N._cavity_sandwich!(Lin,C,L,tmp); @. Lin+=im*gamma
        N._cavity_sandwich!(dLin,C,dL,tmp); @. dLin+=im*v
        N._cavity_sandwich!(Ln,gn,Lin,tmp)
        triple!(dLn,dgn,Lin,adjoint(gn),tmp)
        triple!(term,gn,dLin,adjoint(gn),tmp);dLn .+= term
        triple!(term,gn,Lin,adjoint(dgn),tmp);dLn .+= term
        copyto!(g,gn);copyto!(L,Ln);copyto!(dg,dgn);copyto!(dL,dLn)
    end
    for (out,M) in zip(outputs,(g,L,dg,dL))
        N._cavity_sandwich!(out,C,M,tmp)
    end
    return nothing
end
"Actual covariance tangent; no i(deltaG-deltaG') replacement."
function block!(w::Work,E,h,R,gamma,dR,v,T)
    pr,pl,dpr,dpl=ntuple(i->w.m[24+i],Val(4)); mr,ml,dmr,dml=ntuple(i->w.m[28+i],Val(4))
    side!(w,E,h,R,gamma,dR,v,adjoint(T),(mr,ml,dmr,dml))
    side!(w,E,h,R,gamma,dR,v,T,(pr,pl,dpr,dpl))
    D,G,K,dK,V,dG,A,dA,tmp,term=ntuple(i->w.m[32+i],Val(10))
    @. D=-h-R-(pr+mr)
    for a in axes(D,1);D[a,a]+=E;end
    ds=0.0
    for a in axes(D,1)
        row=0.0
        for b in axes(D,2);row+=abs(D[a,b]);end
        ds=max(ds,row)
    end
    ds>0 || throw(SingularException(0))
    @. G=D/ds; N._production_inverse!(G,w.inverse); G ./= ds
    @. K=gamma-im*(pl+ml)
    @. dK=v-im*(dpl+dml)
    @. V=dR+dpr+dmr
    triple!(dG,G,V,G,tmp)
    N._cavity_sandwich!(A,G,K,tmp)
    triple!(dA,dG,K,adjoint(G),tmp)
    triple!(term,G,dK,adjoint(G),tmp);dA .+= term
    triple!(term,G,K,adjoint(dG),tmp);dA .+= term
    return G,A,dA
end
function spectral_tangent(GP,c,base,v::Field)
    size(v)==size(base.gamma) || throw(DimensionMismatch("tangent field"))
    # Signed Hermitian direction is never passed through a PSD projector.
    norm(v-GP.hermitian_copy(v))<=1e-12*norm(v) || error("non-Hermitian direction")
    dR,roundoff=GP.local_retarded(c,v)
    out=similar(v);ne,nk,nb,_=size(v);hH=c.cavity.hlocal
    workers=N._production_worker_count(c.options,ne*nk);sums=zeros(4,workers)
    N._production_worker_ranges!(ne*nk,workers) do worker,first,last
        w=Work(nb);h,R,gamma,dRt,vt=ntuple(i->w.m[42+i],Val(5))
        for linear in first:last
            e=(linear-1)%ne+1;m=(linear-1)÷ne+1
            for b in 1:nb,a in 1:nb
                h[a,b]=hH[m,a,b];R[a,b]=base.R[e,m,a,b];gamma[a,b]=base.gamma[e,m,a,b]
                dRt[a,b]=dR[e,m,a,b];vt[a,b]=v[e,m,a,b]
            end
            G,A,dA=block!(w,c.problem.grids.ε[e],h,R,gamma,dRt,vt,c.problem.basis.T₊)
            weight=c.problem.grids.wᴱ[e]*c.problem.grids.wᵏ[m]
            for b in 1:nb,a in 1:nb
                out[e,m,a,b]=dA[a,b]
                sums[1,worker]+=weight*abs2(G[a,b]-base.GR[e,m,a,b])
                sums[2,worker]+=weight*abs2(base.GR[e,m,a,b])
                sums[3,worker]+=weight*abs2(A[a,b]-base.A[e,m,a,b])
                sums[4,worker]+=weight*abs2(base.A[e,m,a,b])
            end
        end
    end
    all(isfinite,out) || error("nonfinite spectral tangent")
    s=vec(sum(sums;dims=2))
    return out,(primal_GR_relative=GP.relative(sqrt(s[1]),sqrt(s[2])),
        primal_A_relative=GP.relative(sqrt(s[3]),sqrt(s[4])),PI_roundoff=roundoff)
end
function number_tangent(GP,c,base,dA::Field)
    g=c.problem.grids;mu=base.stats.mu.mu;pref=c.problem.physical.g_s/(2pi)
    numerator=0.0;chi=0.0
    for e in axes(dA,1)
        f,hole=GP.fermi_pair((g.ε[e]-mu)/c.kBT)
        for m in axes(dA,2),a in axes(dA,3)
            w=pref*g.wᴱ[e]*g.wᵏ[m]
            numerator+=w*f*real(dA[e,m,a,a])
            chi+=w*f*hole/c.kBT*real(base.A[e,m,a,a])
        end
    end
    isfinite(chi) && chi>0 || error("nonpositive/nonfinite susceptibility")
    dmu=-numerator/chi;GL=similar(dA);GG=similar(dA);dN=0.0;scale=0.0
    for e in axes(dA,1)
        f,hole=GP.fermi_pair((g.ε[e]-mu)/c.kBT);q=f*hole/c.kBT
        for m in axes(dA,2),b in axes(dA,4),a in axes(dA,3)
            dF=f*dA[e,m,a,b]+q*dmu*base.A[e,m,a,b]
            GL[e,m,a,b]=im*dF;GG[e,m,a,b]=-im*(dA[e,m,a,b]-dF)
            if a==b
                w=pref*g.wᴱ[e]*g.wᵏ[m];dN+=w*real(dF)
                scale+=abs(w*f*real(dA[e,m,a,a]))+abs(w*q*dmu*real(base.A[e,m,a,a]))
            end
        end
    end
    return GL,GG,(chi=chi,dmu=dmu,number_numerator=numerator,
        dN=dN,dN_scale=scale,dN_relative=GP.relative(abs(dN),scale))
end
function jvp(GP,c,base,v::Field,budget::Limits)
    entered!(budget,:jv)
    dA,spectral=spectral_tangent(GP,c,base,v)
    GL,GG,number=number_tangent(GP,c,base,dA)
    # phi reads only A/GL/GG. GR is an unused storage alias, not deltaGR.
    derivative_green=GreenState(dA,GL,GG,dA,base.condition,base.scale)
    raw,channels,sizes=GP.phi(c,derivative_green);empty!(channels)
    all(isfinite,raw) || error("nonfinite raw analytic Jv")
    herm=GP.hermitian_copy(raw)
    representation=GP.norms(raw,herm,c.problem.grids.wᴱ,c.problem.grids.wᵏ)
    budget.jv_completed+=1
    return raw,(spectral=spectral,number=number,channels=sizes,representation=representation)
end
function central_control(GP,SS,c,base,budget)
    pairs=Field[];records=Any[]
    for amplitude in (1e-3,5e-4)
        plus,sp=fresh_phi(GP,SS,c,(1+amplitude).*base.gamma,budget)
        minus,sm=fresh_phi(GP,SS,c,(1-amplitude).*base.gamma,budget)
        @. plus=(plus-minus)/(2amplitude)
        push!(pairs,plus);push!(records,(amplitude=amplitude,plus=sp,minus=sm))
    end
    analytic,info=jvp(GP,c,base,base.gamma,budget)
    g=c.problem.grids
    comparisons=[GP.norms(analytic,D,g.wᴱ,g.wᵏ) for D in pairs]
    refinement=GP.norms(pairs[1],pairs[2],g.wᴱ,g.wᵏ)
    return (;analytic,pairs,records,info,comparisons,refinement)
end
"At most three Arnoldi columns: one earlier validation Jv leaves total cap4."
function arnoldi3(GP,c,base,budget;qualify)
    g=c.problem.grids;rawR=base.raw.-base.gamma;seed=GP.hermitian_copy(rawR)
    seed_conversion=GP.norms(rawR,seed,g.wᴱ,g.wᵏ)
    beta=GP.wnorm(seed,g.wᴱ,g.wᵏ);beta>0 || error("zero Arnoldi seed")
    seed ./= beta;Q=Field[seed];H=zeros(Float64,4,3);records=Any[];columns=0
    for j in 1:3
        raw,record=jvp(GP,c,base,Q[j],budget)
        qualify(record) || error("Jv qualification failed; do not report Ritz spectrum")
        w=GP.hermitian_copy(raw);push!(records,record)
        # Real Hermitian vector space, two-pass weighted reorthogonalization.
        for pass in 1:2,i in 1:j
            h=GP.inner(Q[i],w,g.wᴱ,g.wᵏ);H[i,j]+=h;@. w-=h*Q[i]
        end
        H[j+1,j]=GP.wnorm(w,g.wᴱ,g.wᵏ);columns=j
        H[j+1,j]>0 || break # exact breakdown; no silent norm threshold
        w ./= H[j+1,j];push!(Q,w)
    end
    m=columns;Hb=H[1:m+1,1:m];Hs=Hb[1:m,:]
    eig=eigen(Hs);ritz_residual=abs.(Hb[m+1,m].*eig.vectors[m,:])
    B=-copy(Hb);for j in 1:m;B[j,j]+=1;end
    sv=svd(B);z=sv.V[:,end]
    rhs=zeros(m+1);rhs[1]=beta;gmres_coordinates=B\rhs
    # This is a restricted near-null witness. No Gamma update is formed.
    return (;Q,Hbar=Hb,Ritz=eig.values,Ritz_residual=ritz_residual,
        restricted_B_singular_values=sv.S,near_null_coordinates=z,
        gmres_coordinates,gmres_predicted_residual=norm(B*gmres_coordinates-rhs),
        seed_conversion,seed_norm=beta,records,
        qualification="local projected evidence only; no spectral-radius or full-sigma-min certificate")
end
end # module
