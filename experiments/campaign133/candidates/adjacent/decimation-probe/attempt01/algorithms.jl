module FrozenBoundaryAlgorithms
using LinearAlgebra
export BoundaryWorkspace, finite_centres!, decimate!, diagnostics, bloch_integral,
       scalar_surface, full_finite_centre, selfcheck, spectral!

# Independently written implementation of Schur elimination / Lopez Sancho
# doubling. Mathematical attribution and the exact recurrences: algorithm.md.
# No SIESTA/Kwant source is included, imported, or executed.
mutable struct BoundaryWorkspace
    D0::Matrix{ComplexF64}
    bulk::Matrix{ComplexF64}
    surfaceR::Matrix{ComplexF64}
    surfaceL::Matrix{ComplexF64}
    alpha::Matrix{ComplexF64}
    beta::Matrix{ComplexF64}
    factor::Matrix{ComplexF64}
    rhs::Matrix{ComplexF64}
    cR::Matrix{ComplexF64}
    cL::Matrix{ComplexF64}
    aNew::Matrix{ComplexF64}
    bNew::Matrix{ComplexF64}
    gR::Matrix{ComplexF64}
    gL::Matrix{ComplexF64}
    temp::Matrix{ComplexF64}
    denominator::Matrix{ComplexF64}
    G::Matrix{ComplexF64}
    eye::Matrix{ComplexF64}
end
function BoundaryWorkspace(n::Int)
    arrays = [zeros(ComplexF64,n,n) for _ in 1:16]
    BoundaryWorkspace(arrays[1:7]...,zeros(ComplexF64,n,2n),arrays[8:16]...,
                      Matrix{ComplexF64}(I,n,n))
end

function inverse_into!(out, factor, denominator, eye)
    copyto!(factor,denominator)
    F = lu!(factor;check=true) # bounded LU pivot metadata, no custom inverse
    copyto!(out,eye)
    ldiv!(F,out)
    out
end
function sandwich!(out,T,g,temp)
    mul!(temp,T,g)
    mul!(out,temp,adjoint(T))
    out
end
function spectral!(A,G)
    @inbounds for b in axes(G,2),a in axes(G,1)
        A[a,b]=im*(G[a,b]-conj(G[b,a]))
    end
    A
end

# Record centres at several depths in a single ordinary finite recursion.
# D0 is the physical local inverse resolvent, with no extra eta.
function finite_centres!(dest,w::BoundaryWorkspace,T,depths;conditions=nothing,residuals=nothing)
    fill!(w.gR,0);fill!(w.gL,0)
    target=1
    for depth in 1:last(depths)
        sandwich!(w.cR,T,w.gR,w.temp)
        @. w.denominator=w.D0-w.cR
        inverse_into!(w.gR,w.factor,w.denominator,w.eye)
        sandwich!(w.cL,adjoint(T),w.gL,w.temp)
        @. w.denominator=w.D0-w.cL
        inverse_into!(w.gL,w.factor,w.denominator,w.eye)
        if depth==depths[target]
            sandwich!(w.cR,T,w.gR,w.temp)
            sandwich!(w.cL,adjoint(T),w.gL,w.temp)
            @. w.denominator=w.D0-w.cR-w.cL
            inverse_into!(dest[target],w.factor,w.denominator,w.eye)
            isnothing(conditions) || (conditions[target]=cond(w.denominator))
            isnothing(residuals) || (residuals[target]=backward_residual(w.denominator,dest[target]))
            target+=1
            target>length(depths) && break
        end
    end
    dest
end

function backward_residual(D,G)
    n=size(D,1)
    norm(D*G-I)/(norm(D)*norm(G)+sqrt(n))
end
function diagnostics(D0,T,gR,gL,G)
    dR=D0-T*gR*T';dL=D0-T'*gL*T
    dC=D0-T*gR*T'-T'*gL*T
    psdmin(g)=eigmin(Hermitian((im*(g-g'))))
    scaleR=max(norm(gR),floatmin(Float64))
    scaleL=max(norm(gL),floatmin(Float64))
    (;surface_residual_R=backward_residual(dR,gR),
      surface_residual_L=backward_residual(dL,gL),
      centre_residual=backward_residual(dC,G),
      surface_min_spectral_R=psdmin(gR),surface_min_spectral_L=psdmin(gL),
      surface_min_spectral_scaled_R=psdmin(gR)/scaleR,
      surface_min_spectral_scaled_L=psdmin(gL)/scaleL,
      centre_condition=cond(dC),surface_condition_R=cond(dR),surface_condition_L=cond(dL))
end

# One factorization and 2n RHS at each doubling iteration. Both retarded
# orientations evolve independently: bNew is NOT aNew' for complex D0.
function decimate!(w::BoundaryWorkspace,T;atol_energy,rtol,maxiter=64)
    n=size(T,1)
    copyto!(w.bulk,w.D0);copyto!(w.surfaceR,w.D0);copyto!(w.surfaceL,w.D0)
    copyto!(w.alpha,T);copyto!(w.beta,adjoint(T))
    correction=Inf;stopscale=NaN;iteration=0;converged=false
    max_condition_sample=0.0
    for step in 1:maxiter
        iteration=step
        copyto!(w.factor,w.bulk)
        F=lu!(w.factor;check=true)
        @views copyto!(w.rhs[:,1:n],w.alpha)
        @views copyto!(w.rhs[:,n+1:2n],w.beta)
        ldiv!(F,w.rhs)
        @views mul!(w.cR,w.alpha,w.rhs[:,n+1:2n])
        @views mul!(w.cL,w.beta,w.rhs[:,1:n])
        @views mul!(w.aNew,w.alpha,w.rhs[:,1:n])
        @views mul!(w.bNew,w.beta,w.rhs[:,n+1:2n])
        w.surfaceR .-= w.cR;w.surfaceL .-= w.cL
        w.bulk .-= w.cR;w.bulk .-= w.cL
        copyto!(w.alpha,w.aNew);copyto!(w.beta,w.bNew)
        correction=max(norm(w.cR),norm(w.cL),norm(w.alpha),norm(w.beta))
        stopscale=max(norm(w.D0),norm(w.surfaceR),norm(w.surfaceL))
        all(isfinite,w.bulk) && isfinite(correction) || error("nonfinite doubling state")
        # Sample the initial and near-stop renormalized operator; final
        # physical surface/centre conditioning is recorded independently.
        if step==1 || correction<=atol_energy+rtol*stopscale || step==maxiter
            max_condition_sample=max(max_condition_sample,cond(w.bulk))
        end
        if correction<=atol_energy+rtol*stopscale
            converged=true;break
        end
    end
    inverse_into!(w.gR,w.factor,w.surfaceR,w.eye)
    inverse_into!(w.gL,w.factor,w.surfaceL,w.eye)
    sandwich!(w.cR,T,w.gR,w.temp)
    sandwich!(w.cL,adjoint(T),w.gL,w.temp)
    @. w.denominator=w.D0-w.cR-w.cL
    inverse_into!(w.G,w.factor,w.denominator,w.eye)
    (;converged,iterations=iteration,correction_energy=correction,
      correction_scale_energy=stopscale,bulk_condition_sample=max_condition_sample)
end

# Periodic trapezoid quadrature for the infinite identical-cell bulk.
# Only selected actual (E,k) blocks are sent here, never the full energy grid.
function bloch_integral(D0,T,ntheta::Int)
    n=size(T,1);D=similar(D0);factor=similar(D0);G=similar(D0)
    result=zeros(ComplexF64,n,n);eye=Matrix{ComplexF64}(I,n,n)
    max_condition=0.0
    for j in 0:ntheta-1
        phase=cis(2π*j/ntheta)
        @inbounds for b in 1:n,a in 1:n
            D[a,b]=D0[a,b]-T[a,b]*phase-conj(T[b,a])*conj(phase)
        end
        inverse_into!(G,factor,D,eye)
        result .+= G
        max_condition=max(max_condition,cond(D))
    end
    result ./= ntheta
    (;G=result,max_condition)
end

function scalar_surface(z::Complex,t)
    imag(z)>0 || error("analytic branch requires upper-half-plane z")
    iszero(t) && return inv(z)
    2/(z+sqrt(z-2abs(t))*sqrt(z+2abs(t)))
end
function full_finite_centre(D0,T,depth::Int)
    n=size(T,1);count=2depth+1;D=zeros(ComplexF64,count*n,count*n)
    for cell in 1:count
        rows=(cell-1)*n+1:cell*n
        D[rows,rows].=D0
        if cell<count
            nxt=cell*n+1:(cell+1)*n
            D[rows,nxt].=-T;D[nxt,rows].=-T'
        end
    end
    G=lu(D)\Matrix{ComplexF64}(I,count*n,count*n)
    rows=depth*n+1:(depth+1)*n
    (;centre=G[rows,rows],full_condition=cond(D),backward=backward_residual(D,G))
end

function selfcheck()
    results=Dict{String,Any}()
    # Independent scalar analytic branch, including energies outside band.
    for (j,z) in enumerate((0.37+0.08im,2.4+0.02im,-2.4+0.03im,
                            2.4+1.3e-15im,-2.4+1.3e-15im))
        w=BoundaryWorkspace(1);w.D0[1,1]=z;T=reshape(ComplexF64[0.61],1,1)
        state=decimate!(w,T;atol_energy=1e-14,rtol=1e-13)
        exact=scalar_surface(z,T[1,1]);err=abs(w.gR[1,1]-exact)/abs(exact)
        diag=diagnostics(w.D0,T,w.gR,w.gL,w.G)
        finite=all(isfinite,w.gR) && all(isfinite,w.gL) && all(isfinite,w.G) &&
               all(isfinite,values(diag)) && isfinite(err)
        finite && state.converged && err<1e-10 &&
            max(diag.surface_residual_R,diag.surface_residual_L)<1e-10 || error("scalar branch selfcheck failed")
        results["scalar_"*string(j)]=Dict("relative_error"=>err,"iterations"=>state.iterations,
            "imaginary_energy"=>imag(z),"all_finite"=>finite,
            "surface_residual_R"=>diag.surface_residual_R,"surface_residual_L"=>diag.surface_residual_L)
    end
    # Noncommuting complex 2x2 hopping and dissipative local inverse operator.
    h=ComplexF64[0.13 0.07+0.02im;0.07-0.02im -0.21]
    T=ComplexF64[0.17+0.03im 0.05-0.04im;-0.02+0.06im -0.11+0.01im]
    gamma=ComplexF64[0.12 0.02+0.01im;0.02-0.01im 0.09]
    D0=0.08I-h+0.5im*gamma
    w=BoundaryWorkspace(2);copyto!(w.D0,D0)
    destinations=[zeros(ComplexF64,2,2) for _ in 1:3]
    finite_centres!(destinations,w,T,[1,2,4])
    for (j,depth) in enumerate((1,2,4))
        oracle=full_finite_centre(D0,T,depth)
        err=norm(destinations[j]-oracle.centre)/norm(oracle.centre)
        err<1e-11 || error("complex 2x2 finite block selfcheck failed")
        results["complex_finite_"*string(depth)]=Dict("relative_error"=>err,"full_backward"=>oracle.backward)
    end
    state=decimate!(w,T;atol_energy=1e-14,rtol=1e-13)
    diag=diagnostics(D0,T,w.gR,w.gL,w.G)
    B=bloch_integral(D0,T,2048)
    err=norm(w.G-B.G)/norm(B.G)
    state.converged && err<1e-10 && max(diag.surface_residual_R,diag.surface_residual_L)<1e-10 ||
        error("complex 2x2 surface/Bloch selfcheck failed")
    results["complex_infinite"]=Dict("Bloch_relative_error"=>err,
        "iterations"=>state.iterations,"surface_residual_R"=>diag.surface_residual_R,
        "surface_residual_L"=>diag.surface_residual_L,
        "beta_not_adjoint_alpha_defect"=>norm(w.beta-w.alpha'))
    # Weak, PSD rank-one damping with an exactly undamped orthogonal
    # direction. The real energy lies outside all bands (bounded by
    # ||h||+2||T||), so no artificial eta or pole regularization is needed.
    v=ComplexF64[1,0.4+0.2im]
    weak_gamma=3e-15*(v*v')
    energy=1.7
    energy>opnorm(h)+2opnorm(T) || error("weak oracle not outside spectrum bound")
    weakD=energy*I-h+0.5im*weak_gamma
    weak=BoundaryWorkspace(2);copyto!(weak.D0,weakD)
    weakfinite=[zeros(ComplexF64,2,2)]
    finite_centres!(weakfinite,weak,T,[4])
    full=full_finite_centre(weakD,T,4)
    finite_error=norm(weakfinite[1]-full.centre)/norm(full.centre)
    weakstate=decimate!(weak,T;atol_energy=1e-14,rtol=1e-13)
    weakdiag=diagnostics(weakD,T,weak.gR,weak.gL,weak.G)
    B1=bloch_integral(weakD,T,2048);B2=bloch_integral(weakD,T,4096)
    refinement=norm(B2.G-B1.G)/norm(B2.G)
    agreement=norm(B2.G-weak.G)/norm(B2.G)
    finite=all(isfinite,weak.G) && all(isfinite,weak.gR) && all(isfinite,weak.gL) &&
        all(isfinite,B1.G) && all(isfinite,B2.G) && all(isfinite,values(weakdiag)) &&
        all(isfinite,(finite_error,full.backward,refinement,agreement))
    finite && weakstate.converged && max(finite_error,refinement,agreement)<1e-10 &&
        max(weakdiag.surface_residual_R,weakdiag.surface_residual_L)<1e-10 &&
        min(weakdiag.surface_min_spectral_scaled_R,weakdiag.surface_min_spectral_scaled_L)>-1e-10 ||
        error("weak rank-deficient damping selfcheck failed")
    results["weak_rank_one_damping"]=Dict("all_finite"=>finite,
        "gamma_scale"=>norm(weak_gamma),"outside_spectrum_margin"=>energy-opnorm(h)-2opnorm(T),
        "finite4_full_block_relative"=>finite_error,"full_backward"=>full.backward,
        "iterations"=>weakstate.iterations,"surface_residual_R"=>weakdiag.surface_residual_R,
        "surface_residual_L"=>weakdiag.surface_residual_L,
        "Bloch_2048_to_4096_relative"=>refinement,"Bloch_decimation_relative"=>agreement)
    results
end
end # module
