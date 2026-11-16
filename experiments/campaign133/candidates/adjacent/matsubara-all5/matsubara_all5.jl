# WRITE-ONLY prototype; no loader/CLI/autostart. All energies MUST share native units.
module MatsubaraAll5
using LinearAlgebra, QCLNEGF
const QN = QCLNEGF.QCLNumerics
const STATIC = (:acoustic, :impurity, :IFR, :alloy)
const Field = Array{ComplexF64,4}
fermi(x) = x >= 0 ? exp(-x)/(1+exp(-x)) : 1/(1+exp(x))
herm(X) = (X+X')/2
relative(a,b) = iszero(b) ? (iszero(a) ? 0.0 : Inf) : a/b

"ROOT supplies source arrays; this constructor performs no file access."
function context(hH, Tplus, kernels, wk, qK, beta, mu, Omega; spin=2)
    beta>0 && Omega>0 && all(>(0),wk) || error("invalid beta/Omega/radial measure")
    nk,nb,nb2=size(hH); nb==nb2 && length(wk)==nk || error("shape")
    Set(keys(kernels))==Set((STATIC...,:LO)) || error("exactly ALL5 required")
    ops=Dict(s=>QN._production_kernel(kernels[s],AlgorithmOptions(contraction=:dense_blas)) for s in keys(kernels))
    full=Matrix{ComplexF64}[]; atoms=[]; centre=(4nb+1):(5nb)
    for k in 1:nk
        h=Matrix(hH[k,:,:]); ishermitian(h) || error("supply exact qualified hH")
        H=zeros(ComplexF64,9nb,9nb)
        for j in 0:8
            r=(j*nb+1):((j+1)*nb); H[r,r]=h
            if j<8
                t=((j+1)*nb+1):((j+2)*nb)
                H[r,t]=Tplus; H[t,r]=Tplus'
            end
        end
        es=eigen(Hermitian(H)); V=Matrix(es.vectors[centre,:])
        push!(full,H); push!(atoms,(eps=es.values.-mu,V=V))
    end
    return (;full,atoms,ops,kernels,wk=Float64.(wk),qK,beta=Float64(beta),mu=Float64(mu),
        Omega=Float64(Omega),NB=inv(expm1(beta*Omega)),nk,nb,centre,spin)
end

function applyK(c,s,X)
    out=zeros(ComplexF64,size(X))
    QN._production_static_contraction!(out,c.ops[s],X,c.wk,c.qK[s],32;
        parallel_backend=:blas,worker_count=1)
    return out
end

"Independent six-index loop; does not call flattening/native contraction helpers."
function literalK(c,s,X)
    ne,nk,nb,_=size(X); Y=zeros(ComplexF64,size(X)); K=c.kernels[s]
    for e in 1:ne,m in 1:nk,a in 1:nb,b in 1:nb
        v=0.0im
        for mp in 1:nk,cc in 1:nb,d in 1:nb
            v+=c.wk[mp]*K[m,mp,a,cc,d,b]*X[e,mp,cc,d]
        end
        Y[e,m,a,b]=c.qK[s]*v
    end
    return Y
end

"Prospective qualification:17 non-Hermitian slices, all five raw source kernels."
function qualify_native_ops(c; rtol=2e-12)
    X=Array{ComplexF64}(undef,17,c.nk,c.nb,c.nb)
    for i in eachindex(X); X[i]=complex(sin(i/7),cos(i/11)); end
    errors=Dict{Symbol,Float64}()
    for s in keys(c.ops)
        a=applyK(c,s,X); b=literalK(c,s,X)
        errors[s]=relative(norm(a-b),norm(b))
        errors[s]<=rtol || error("native/literal contraction mismatch: $s")
    end
    return errors
end

"Positive n=0:N-1 reference data; z physical=mu+iω, eps atoms relative mu."
function reference(c,N)
    w=[(2n+1)*pi/c.beta for n in 0:N-1]
    G=zeros(ComplexF64,N,c.nk,c.nb,c.nb); C=similar(G)
    rho=zeros(ComplexF64,c.nk,c.nb,c.nb)
    for k in 1:c.nk
        (;eps,V)=c.atoms[k]; f=fermi.(c.beta.*eps)
        rho[k,:,:]=V*Diagonal(f)*V'
        for n in 1:N
            z=im*w[n]
            G[n,k,:,:]=V*Diagonal(inv.(z.-eps))*V'
            coef=(c.NB+1 .- f)./(z.-eps.-c.Omega)+(c.NB .+ f)./(z.-eps.+c.Omega)
            C[n,k,:,:]=V*Diagonal(coef)*V'
        end
    end
    return (;w,G,lo=applyK(c,:LO,C),rho)
end

"Exact symmetric finite Mats convolution of ΔG; omitted ΔG tail remains an error."
function lo_delta(c,delta,w)
    N=length(w); both=zeros(ComplexF64,2N,c.nk,c.nb,c.nb)
    for j in -N:N-1,k in 1:c.nk
        both[j+N+1,k,:,:]=j>=0 ? delta[j+1,k,:,:] : Matrix(delta[-j,k,:,:])'
    end
    D=[2*c.Omega/((2pi*(n-j)/c.beta)^2+c.Omega^2) for n in 0:N-1,j in -N:N-1]
    conv=reshape((D*reshape(both,2N,:))/c.beta,size(delta))
    return applyK(c,:LO,conv)
end

"Same central-driven Σ repeated in all nine cells, not a full cellwise SCBA."
function greens(c,S,w)
    G=similar(S); dyson=0.0; eye=Matrix{ComplexF64}(I,9*c.nb,9*c.nb)
    for n in eachindex(w),k in 1:c.nk
        D=(c.mu+im*w[n])*eye-c.full[k]
        for j in 0:8
            r=(j*c.nb+1):((j+1)*c.nb); D[r,r]-=S[n,k,:,:]
        end
        invD=D\eye
        dyson=max(dyson,norm(D*invD-eye)/(norm(D)*norm(invD)+norm(eye)))
        G[n,k,:,:]=invD[c.centre,c.centre]
    end
    return G,dyson
end

function density(c,G,ref)
    rho=copy(ref.rho)
    for k in 1:c.nk,n in eachindex(ref.w)
        d=Matrix(G[n,k,:,:]-ref.G[n,k,:,:]); rho[k,:,:].+=(2/c.beta).*herm(d)
    end
    number=c.spin*sum(c.wk[k]*real(tr(rho[k,:,:])) for k in 1:c.nk)
    pauli_min=minimum(eigmin(Hermitian(rho[k,:,:])) for k in 1:c.nk)
    pauli_max=maximum(eigmax(Hermitian(rho[k,:,:])) for k in 1:c.nk)
    return rho,(;number,pauli_min,pauli_max)
end

function causal_max(X)
    maximum(eigmax(Hermitian((Matrix(X[n,k,:,:])-Matrix(X[n,k,:,:])')/(2im)))
        for n in axes(X,1),k in axes(X,2))
end

function cp_bound(c)
    X=zeros(ComplexF64,1,c.nk,c.nb,c.nb)
    for k in 1:c.nk; X[1,k,:,:]=Matrix{ComplexF64}(I,c.nb,c.nb); end
    B=sum(applyK(c,s,X) for s in STATIC)+(2*c.NB+1).*applyK(c,:LO,X)
    L=maximum(opnorm(Hermitian(Matrix(B[1,k,:,:]))) for k in 1:c.nk)
    return (;L,omega0=pi/c.beta,conditional_L_over_omega0_squared=L/(pi/c.beta)^2,
        scope="infinite Mats causal map only; not a finite-cutoff uniqueness assertion")
end

"N128/256, max16 each intended. No density/μ normalization or automated μ fit."
function pilot(c,N;steps=16,alpha=0.3,S0=nothing,soft_seconds=540.0,callback=identity)
    N in (128,256) && 1<=steps<=16 && 0<alpha<=1 || error("bounded pilot options")
    ref=reference(c,N); S=S0===nothing ? zeros(ComplexF64,size(ref.G)) : ComplexF64.(S0)
    size(S)==size(ref.G) || error("initial Σ dimensions")
    hist=[]; began=time(); latest=nothing
    for it in 1:steps
        time()-began<soft_seconds || return latest===nothing ? (;status="soft_budget_unaligned",S,history=hist) :
            (;status="soft_budget_aligned",latest...,history=hist)
        G,rD=greens(c,S,ref.w)
        Phi=sum(applyK(c,s,G) for s in STATIC)+ref.lo+lo_delta(c,G-ref.G,ref.w)
        rho,d= density(c,G,ref)
        rec=(;iteration=it,raw_sigma_absolute=norm(Phi-S),
            raw_sigma_relative=relative(norm(Phi-S),norm(Phi)),dyson_backward=rD,
            sigma_causality_max=causal_max(Phi),G_causality_max=causal_max(G),d...)
        push!(hist,rec); latest=(;S,G,Phi,rho,w=ref.w,record=rec); callback(latest)
        st=512eps(Float64)*max(maximum(abs,Phi),eps(Float64))
        gt=512eps(Float64)*max(maximum(abs,G),eps(Float64))
        if rec.sigma_causality_max>st || rec.G_causality_max>gt
            return (;status="noncausal_candidate",S,G,Phi,rho,w=ref.w,history=hist,
                sigma_roundoff_allowance=st,G_roundoff_allowance=gt)
        end
        if it==steps
            return (;status="iteration_budget",S,G,Phi,rho,w=ref.w,history=hist,
                truncation="ΔG outside ±N omitted; compare128/256, no certified tail bound")
        end
        S=(1-alpha).*S.+alpha.*Phi # preserve aligned latest.S
    end
end

"Independent one-level raw Mats sum vs exact rational Fock; no reference subtraction."
function onelevel_oracle(beta,Omega;N=256,g2=0.001,epsilon=0.02)
    NB=inv(expm1(beta*Omega)); f=fermi(beta*epsilon)
    rows=[]
    for n in (0,1,7),M in (N,2N)
        z=im*(2n+1)*pi/beta
        exact=g2*((NB+1-f)/(z-epsilon-Omega)+(NB+f)/(z-epsilon+Omega))
        raw=0.0im
        for j in -M:M-1
            wj=(2j+1)*pi/beta; nu=2pi*(n-j)/beta
            raw+=g2*2Omega/(nu^2+Omega^2)/(im*wj-epsilon)/beta
        end
        push!(rows,(;n,M,raw,exact,absolute=abs(raw-exact),relative=relative(abs(raw-exact),abs(exact))))
    end
    return rows
end

"Shared low-frequency comparison, fixed μ; both arrays retained only at this small N."
function compare_cutoffs(a,b,c;low=32)
    n=min(low,length(a.w),length(b.w)); ar=@view a.G[1:n,:,:,:]; br=@view b.G[1:n,:,:,:]
    sa=@view a.S[1:n,:,:,:]; sb=@view b.S[1:n,:,:,:]
    return (;low=n,G_absolute=norm(ar-br),G_relative=relative(norm(ar-br),norm(br)),
        S_absolute=norm(sa-sb),S_relative=relative(norm(sa-sb),norm(sb)),
        rho_absolute=norm(a.rho-b.rho),number_difference=
            c.spin*sum(c.wk[k]*real(tr(a.rho[k,:,:]-b.rho[k,:,:])) for k in 1:c.nk))
end
"Actual reference/lo_delta/greens/density wiring: two-pole G from known auxiliary level."
function wired_onelevel(beta,Omega;g2=Omega^2/5)
    e=Omega/4; eb=-3Omega/4; a=Omega^2/4
    split=sqrt((e-eb)^2+4a); lp=(e+eb+split)/2; lm=(e+eb-split)/2
    wp=(lp-eb)/split; wm=1-wp
    raw=Dict(s=>fill(s===:LO ? g2 : 0.0,1,1,1,1,1,1) for s in (STATIC...,:LO))
    c=context(reshape(ComplexF64[e],1,1,1),zeros(ComplexF64,1,1),raw,[1.0],
        Dict(s=>1.0 for s in keys(raw)),beta,0.0,Omega;spin=1)
    records=[]
    exactF(z,x)=g2*((c.NB+1-fermi(beta*x))/(z-x-Omega)+(c.NB+fermi(beta*x))/(z-x+Omega))
    for N in (128,256)
        ref=reference(c,N); S=reshape([a/(im*w-eb) for w in ref.w],N,1,1,1)
        G,rD=greens(c,S,ref.w); rho,_=density(c,G,ref)
        target=reshape([wp/(im*w-lp)+wm/(im*w-lm) for w in ref.w],N,1,1,1)
        exact=reshape([wp*exactF(im*w,lp)+wm*exactF(im*w,lm) for w in ref.w],N,1,1,1)
        expectedref=reshape([exactF(im*w,e) for w in ref.w],N,1,1,1)
        out=ref.lo+lo_delta(c,G-ref.G,ref.w)
        push!(records,(;N,G_relative=relative(norm(G-target),norm(target)),
            reference_LO_relative=relative(norm(ref.lo-expectedref),norm(expectedref)),
            LO_relative=relative(norm(out-exact),norm(exact)),LO_absolute=norm(out-exact),
            density_absolute=abs(real(rho[1,1,1])-(wp*fermi(beta*lp)+wm*fermi(beta*lm))),
            dyson_backward=rD))
    end
    return records
end
end # module
