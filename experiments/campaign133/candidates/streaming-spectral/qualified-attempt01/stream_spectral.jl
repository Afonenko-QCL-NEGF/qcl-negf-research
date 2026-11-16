# Experimental operator-only helper. Frozen Core e171; no execution on include.
module StreamingSpectral
using QCLNEGF, LinearAlgebra
const N=QCLNEGF.QCLNumerics
struct LazyL{A}; gamma::A; end
Base.@propagate_inbounds Base.getindex(x::LazyL,e,m,a,b)=im*x.gamma[e,m,a,b]
struct LazyZero end
Base.getindex(::LazyZero,e,m,a,b)=zero(ComplexF64)
# Per-worker destination: native _cavity_side! writes only its current (e,m).
struct BlockSink; matrix::Matrix{ComplexF64}; end
Base.@propagate_inbounds Base.setindex!(x::BlockSink,z,e,m,a,b)=(x.matrix[a,b]=z)
block_destination(n)=(Σᴿ=BlockSink(zeros(ComplexF64,n,n)),Σˡ=BlockSink(zeros(ComplexF64,n,n)),Σᵍ=BlockSink(zeros(ComplexF64,n,n)))
function number_weights(A,g,spin)
    q=zeros(Float64,size(A,1));pref=Float64(spin)/(2pi)
    @inbounds for e in axes(A,1)
        momentum=0.0
        for m in axes(A,2)
            trace=0.0
            for a in axes(A,3);trace+=real(A[e,m,a,a]);end
            momentum+=g.wᵏ[m]*trace
        end
        q[e]=pref*g.wᴱ[e]*momentum
    end
    return q
end
function stream_spectrum(GP,c,gamma::Array{ComplexF64,4};selected_indices=Int[])
    p=c.problem;o=c.options;g=p.grids;cp=c.cavity;NE,Nk,Nb,Nb2=size(gamma)
    (NE,Nk,Nb)==(p.numerical.N_E,p.numerical.N_k,p.numerical.N_b) && Nb==Nb2 || throw(DimensionMismatch("Gamma axes"))
    o.algorithms.embedding===:finite_chain && o.algorithms.embedding_periods==4 || error("finite_chain4 only")
    o.algorithms.energy_shift===:sparse_plan && o.algorithms.hilbert===:product_integration || error("frozen operator policy")
    iszero(N._scaled_physics(p.physical,p.scales).Eᵖ) || error("V0 only")
    c.cache.product_hilbert_operator isa N._ProductionPIFFTHilbertPlan || error("uniform PI FFT only")
    all(isfinite,gamma) || error("nonfinite Gamma")
    N._cache_source_contract(c.cache,p);N._cache_shift_contract(c.cache,p)
    cp.source_grids===g && cp.source_h==c.h && cp.source_energy==g.ε && cp.source_weights==g.wᴱ &&
        cp.periods==4 && cp.energy_shift===:sparse_plan && iszero(cp.drop) || error("stale original cavity plan")
    hH=copy(cp.hlocal);plan=N._cavity_production_plan(p,hH,o)
    isequal(plan.hlocal,hH) || error("all9 hH repeat differs")
    dh=norm(hH-c.h);hs=norm(c.h)
    GP.relative(dh,hs)<=1e-10 || error("original h correction target")
    chosen=sort!(unique!(Int.(selected_indices)))
    length(chosen)<=128 && all(e->1<=e<=NE,chosen) || error("selected-block envelope")
    position=zeros(Int,NE);for (i,e) in enumerate(chosen);position[e]=i;end
    names=(:Gamma,:Lambda_raw,:Lambda_H_preassembly,:localR,:totalR,:GR,:A_cov,:A_sub,:Gamma_physical_total,:Gamma_eff)
    samples=Dict(name=>zeros(ComplexF64,length(chosen),Nk,Nb,Nb) for name in names)
    raw,roundoff=N._production_pi_fft_hilbert_transform(gamma,g.ε,g.wᴱ,c.cache.product_hilbert_operator)
    LH=GP.hermitian_copy(raw);correction=GP.norms(raw,LH,g.wᴱ,g.wᵏ)
    correction.relative<=1e-10 && correction.weighted_relative<=1e-10 || error("original PI symmetry target")
    @inbounds for (i,e) in enumerate(chosen),m in 1:Nk,b in 1:Nb,a in 1:Nb
        samples[:Gamma][i,m,a,b]=gamma[e,m,a,b]
        samples[:Lambda_raw][i,m,a,b]=raw[e,m,a,b]
        samples[:Lambda_H_preassembly][i,m,a,b]=LH[e,m,a,b]
    end
    raw=nothing;GC.gc(true) # release global diagnostic field before the sweep
    S=LH;LH=nothing
    @inbounds @simd for i in eachindex(S);S[i]-=im*gamma[i]/2;end
    family=(Σᴿ=S,Σˡ=LazyL(gamma),Σᵍ=LazyZero())
    GR=similar(gamma);A=similar(gamma);condition=zeros(Float64,NE,Nk);scale=similar(condition)
    workers=N._production_worker_count(o,NE*Nk);partials=zeros(Float64,8,workers)
    Tplus=p.basis.T₊;Tminus=adjoint(Tplus)
    N._production_worker_ranges!(NE*Nk,workers) do worker,first,last
        w=N._CavityWorkspace(Nb);plus=block_destination(Nb);minus=block_destination(Nb)
        D=zeros(ComplexF64,Nb,Nb);gr=similar(D);K=similar(D);tmp=similar(D);ac=similar(D);totalR=similar(D)
        invwork=N._ProductionInverseWorkspace(Nb)
        @inbounds for linear in first:last
            e=(linear-1)%NE+1;m=(linear-1)÷NE+1
            N._cavity_side!(w,family,hH,plan.remaps,e,m,plan.source_energy,plan.drop,4,Tminus,:minus,minus)
            N._cavity_side!(w,family,hH,plan.remaps,e,m,plan.source_energy,plan.drop,4,Tplus,:plus,plus)
            for b in 1:Nb,a in 1:Nb
                totalR[a,b]=S[e,m,a,b]+(plus.Σᴿ.matrix[a,b]+minus.Σᴿ.matrix[a,b])
                K[a,b]=gamma[e,m,a,b]-im*(plus.Σˡ.matrix[a,b]+minus.Σˡ.matrix[a,b])
                D[a,b]=(a==b ? g.ε[e]+im*0.0 : 0.0+0.0im)-hH[m,a,b]-totalR[a,b]
            end
            # Same owning inverse, scale arithmetic and condition as native central Dyson.
            ds=0.0
            for a in 1:Nb
                rowsum=0.0
                for b in 1:Nb;rowsum+=abs(D[a,b]);end
                ds=max(ds,rowsum)
            end
            ds>0 || throw(SingularException(0))
            for b in 1:Nb,a in 1:Nb;gr[a,b]=D[a,b]/ds;end
            N._production_inverse!(gr,invwork)
            for b in 1:Nb,a in 1:Nb;gr[a,b]/=ds;GR[e,m,a,b]=gr[a,b];end
            condition[e,m]=cond(D);scale[e,m]=ds
            mul!(tmp,gr,K);mul!(ac,tmp,adjoint(gr)) # qualified covariance arithmetic
            sample=position[e];weight=g.wᴱ[e]*g.wᵏ[m]
            for b in 1:Nb,a in 1:Nb
                A[e,m,a,b]=ac[a,b]
                sub=im*(gr[a,b]-conj(gr[b,a]))
                eff=im*(hH[m,a,b]+totalR[a,b]-conj(hH[m,b,a]+totalR[b,a]))
                d=abs2(ac[a,b]-sub);n=abs2(ac[a,b]);v=abs2(eff-K[a,b]);z=abs2(eff)
                partials[1,worker]+=d;partials[2,worker]+=n;partials[3,worker]+=weight*d;partials[4,worker]+=weight*n
                partials[5,worker]+=v;partials[6,worker]+=z;partials[7,worker]+=weight*v;partials[8,worker]+=weight*z
                if sample>0
                    samples[:localR][sample,m,a,b]=S[e,m,a,b];samples[:totalR][sample,m,a,b]=totalR[a,b]
                    samples[:GR][sample,m,a,b]=gr[a,b];samples[:A_cov][sample,m,a,b]=ac[a,b];samples[:A_sub][sample,m,a,b]=sub
                    samples[:Gamma_physical_total][sample,m,a,b]=K[a,b];samples[:Gamma_eff][sample,m,a,b]=eff
                end
            end
        end
    end
    all(isfinite,GR) && all(isfinite,A) && all(isfinite,condition) && all(isfinite,scale) || error("nonfinite spectral result")
    q=number_weights(A,g,p.physical.g_s);all(isfinite,q) || error("nonfinite q")
    sums=vec(sum(partials;dims=2))
    ledger=Dict("operator"=>"stream_native_covariance_exact_plan_hlocal_all9","native_acceptance_override"=>false,
        "h_correction"=>Dict("absolute"=>dh,"scale"=>hs,"relative"=>GP.relative(dh,hs)),
        "Lambda_projection"=>Dict(string(k)=>v for (k,v) in pairs(correction)),"roundoff"=>roundoff,
        "spectral_difference_absolute"=>sqrt(sums[1]),"spectral_difference_relative"=>GP.relative(sqrt(sums[1]),sqrt(sums[2])),
        "spectral_difference_weighted_absolute"=>sqrt(sums[3]),"spectral_difference_weighted_relative"=>GP.relative(sqrt(sums[3]),sqrt(sums[4])),
        "source_difference_absolute"=>sqrt(sums[5]),"source_difference_relative"=>GP.relative(sqrt(sums[5]),sqrt(sums[6])),
        "source_difference_weighted_absolute"=>sqrt(sums[7]),"source_difference_weighted_relative"=>GP.relative(sqrt(sums[7]),sqrt(sums[8])),
        "scalar_reduction"=>"per-worker ordered blocks; may differ from parent full-field norm reduction in last bits",
        "q_minimum"=>minimum(q),"q_negative_count"=>count(<(0.0),q),"q_capacity"=>sum(q),
        "cavity_plan_bytes"=>Base.summarysize(plan),"workers"=>workers,"source_greater_storage"=>"exact-zero scalar adapter")
    return (;GR,A,condition,scale,roundoff,q,ledger,selected_indices=chosen,samples)
end
end # module
