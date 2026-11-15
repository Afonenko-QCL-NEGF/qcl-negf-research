module EquilibriumGammaPilot

using QCLNEGF, Unitful, LinearAlgebra, SparseArrays, TOML, SHA, HDF5, Dates
const N = QCLNEGF.QCLNumerics
const Field = Array{ComplexF64,4}
const STATIC_IDS = (:acoustic, :impurity, :IFR, :alloy)

"A failed physical domain check; programming/API errors are never a mixing fallback."
struct Infeasible <: Exception
    reason::String
end
Base.showerror(io::IO, e::Infeasible) = print(io, e.reason)
struct NativeStop <: Exception
    state::SCBAResult
end

"Only known inverse failures; EVERY nested cause must be numerical."
function known_inverse_failure(error)
    error isa LinearAlgebra.SingularException && return error.info >= 0
    # Negative LAPACK INFO means an invalid argument, not numerical domain.
    error isa LinearAlgebra.LAPACKException && return error.info > 0
    if error isa TaskFailedException
        istaskfailed(error.task) || return false
        stack = Base.current_exceptions(error.task;backtrace=false)
        return !isempty(stack) && all(entry->known_inverse_failure(entry.exception),stack)
    elseif error isa CompositeException
        return !isempty(error.exceptions) && all(known_inverse_failure,error.exceptions)
    elseif error isa CapturedException
        return known_inverse_failure(error.ex)
    end
    return false
end
function inverse_domain(f::F) where F
    try
        return f()
    catch error
        known_inverse_failure(error) || rethrow()
        throw(Infeasible("known numerical inverse failure: "*sprint(showerror,error)))
    end
end

mutable struct Budget
    started::Float64
    stage_started::Float64
    stage::Symbol
    pi::Int
    precheck::Int
    native::Int
    native_events::Dict{String,Int}
    cfg::Dict{String,Any}
end
Budget(cfg) = Budget(time(), time(), :precheck, 0, 0, 0,
    Dict("materialization_started"=>0,"materialization_completed"=>0,
        "solver_entered"=>0,"solver_observed_completed"=>0,
        "audit_started"=>0,"audit_map_completed"=>0,"audit_report_written"=>0),cfg)
function native_accounting(b::Budget)
    ev = copy(b.native_events)
    upper = ev["materialization_started"]+30*ev["solver_entered"]+ev["audit_started"]
    lower = ev["materialization_started"]+ev["solver_observed_completed"]+ev["audit_started"]
    completed = ev["materialization_completed"]+ev["solver_observed_completed"]+ev["audit_map_completed"]
    return Dict("events"=>ev,"started_lower_bound"=>lower,"started_upper_bound"=>upper,
        "completed_lower_bound"=>completed,"structural_cap"=>32,"solver_reserved_cap"=>30,
        "meaning"=>"materialization/audit counted before call; solver observations count completed maps only; an interrupted solver map is bounded by reserved max30, not observed as completed")
end
function check_time(b::Budget)
    time()-b.started < b.cfg["overall_seconds"] || error("overall wall budget exhausted")
    key = b.stage === :precheck ? "precheck_seconds" :
          b.stage === :pilot ? "pilot_seconds" : "native_seconds"
    time()-b.stage_started < b.cfg[key] || error("$(b.stage) wall budget exhausted")
end
stage_expired(b::Budget) = time()-b.started >= b.cfg["overall_seconds"] ||
    time()-b.stage_started >= b.cfg[b.stage === :precheck ? "precheck_seconds" :
        b.stage === :pilot ? "pilot_seconds" : "native_seconds"]
function actual_call!(b::Budget; pi=false)
    check_time(b)
    if pi
        b.pi < b.cfg["pi_calls_total"] || error("total Pi budget exhausted")
        b.pi += 1                         # failures/rejections count too
    end
    if b.stage === :precheck
        b.precheck < b.cfg["precheck_actual_calls"] || error("precheck call cap")
        b.precheck += 1
    end
end

"Concrete type parameters form the barrier around all native array operators."
struct Context{P,C,O,H,CP}
    problem::P
    cache::C
    options::O
    h::H
    cavity::CP
    kBT::Float64
    target::Float64
end

relative(n::Real, d::Real) = d == 0 ? (n == 0 ? 0.0 : Inf) : n/d
function inner(X::Field, Y::Field, we::Vector{Float64}, wk::Vector{Float64})
    s = 0.0
    @inbounds for b in axes(X,4), a in axes(X,3), m in axes(X,2), e in axes(X,1)
        s += we[e]*wk[m]*real(conj(X[e,m,a,b])*Y[e,m,a,b])
    end
    return s
end
wnorm(X::Field, we, wk) = sqrt(max(0.0, inner(X,X,we,wk)))
function norms(X::Field, Y::Field, we, wk)
    # Numerator is a raw map difference, never a mixed step. No additive floor.
    d2, n2 = 0.0, 0.0
    dw2, nw2 = 0.0, 0.0
    @inbounds for b in axes(X,4), a in axes(X,3), m in axes(X,2), e in axes(X,1)
        q = abs2(X[e,m,a,b]-Y[e,m,a,b]); s = abs2(X[e,m,a,b])
        w = we[e]*wk[m]
        d2 += q; n2 += s; dw2 += w*q; nw2 += w*s
    end
    return (absolute=sqrt(d2), scale=sqrt(n2), relative=relative(sqrt(d2),sqrt(n2)),
            weighted_absolute=sqrt(dw2), weighted_scale=sqrt(nw2),
            weighted_relative=relative(sqrt(dw2),sqrt(nw2)))
end

"Native backward-error allowance (Frobenius block scale), no clipping."
function cone(X::Field; require_nonzero=true)
    all(isfinite, X) || throw(Infeasible("nonfinite field"))
    require_nonzero && iszero(norm(X)) && throw(Infeasible("zero Gamma/domain not admitted"))
    nb = size(X,3)
    block = zeros(ComplexF64, nb, nb)
    values = zeros(Float64,nb)
    work = zeros(ComplexF64,max(1,3nb))
    rwork = zeros(Float64,max(1,5nb))
    info = zeros(LinearAlgebra.BlasInt,1)
    minvalue, maxdefect, worst_e, worst_m = Inf, 0.0, 0, 0
    worst_scale,worst_allowance = NaN,NaN
    minmargin = Inf
    @inbounds for m in axes(X,2), e in axes(X,1)
        sq, defect = 0.0, 0.0
        for b in 1:nb, a in 1:nb
            z = X[e,m,a,b]
            sq += abs2(z)
            defect += abs2(z-conj(X[e,m,b,a]))
            block[a,b] = (z+conj(X[e,m,b,a]))/2
        end
        scale = sqrt(sq)
        allowance = N.psd_error_budget(scale,nb;construction_scale=1.0)
        herm = sqrt(defect)
        herm <= allowance || throw(Infeasible("non-Hermitian field at ($e,$m): $herm > $allowance"))
        # Reuse the owning production LAPACK workspace, avoiding millions of
        # per-block eigvals! workspace allocations on the actual S01 grid.
        N._production_hermitian_eigvals!(block,values,work,rwork,info) ||
            throw(Infeasible("PSD eigensolve failure at ($e,$m)"))
        low = first(values)
        margin = low+allowance
        margin >= 0 || throw(Infeasible("negative PSD field at ($e,$m): $low; budget=$allowance"))
        if low < minvalue
            minvalue, worst_e, worst_m = low, e, m
            worst_scale,worst_allowance = scale,allowance
        end
        minmargin = min(minmargin,margin)
        maxdefect = max(maxdefect,herm)
    end
    return (minimum_eigenvalue=minvalue, minimum_margin=minmargin,
            hermiticity_absolute=maxdefect, energy_index=worst_e, k_index=worst_m,
            block_frobenius_scale=worst_scale,backward_error=worst_allowance)
end
function hermitian_copy(X::Field)
    Y = similar(X)
    @inbounds for b in axes(X,4), a in axes(X,3), m in axes(X,2), e in axes(X,1)
        Y[e,m,a,b] = (X[e,m,a,b]+conj(X[e,m,b,a]))/2
    end
    return Y
end
function gamma_residuals(rawphi::Field,gamma::Field,we,wk)
    mapped = hermitian_copy(rawphi)
    raw = norms(rawphi,gamma,we,wk)
    working = norms(mapped,gamma,we,wk)
    representation = norms(rawphi,mapped,we,wk)
    return mapped,raw,working,representation
end

@inline function fermi_pair(x::Float64)
    if x >= 0
        z = exp(-x); d = 1+z
        return z/d, 1/d
    else
        z = exp(x); d = 1+z
        return 1/d, z/d
    end
end
function number_mu(energy, q, mu, kBT)
    n, derivative = 0.0, 0.0
    @inbounds for e in eachindex(energy,q)
        f,hole = fermi_pair((energy[e]-mu)/kBT)
        n += q[e]*f
        derivative += q[e]*f*hole/kBT
    end
    return n,derivative
end
function fit_mu(energy::Vector{Float64}, q::Vector{Float64}, target::Float64, kBT::Float64)
    all(isfinite,q) && all(>=(0.0),q) || throw(Infeasible("nonpositive/nonfinite raw spectral weights"))
    capacity = sum(q)
    0 < target < capacity || throw(Infeasible("mu capacity failure target=$target capacity=$capacity"))
    lo,hi = first(energy)-50kBT,last(energy)+50kBT
    nlo,_ = number_mu(energy,q,lo,kBT); nhi,_ = number_mu(energy,q,hi,kBT)
    nlo < target < nhi || throw(Infeasible("mu root not bracketed"))
    mu,number,derivative,iteration = NaN,NaN,NaN,0
    for i in 1:160
        iteration = i
        mu = (lo+hi)/2
        number,derivative = number_mu(energy,q,mu,kBT)
        abs(number-target) <= 2e-13*target && break
        number < target ? (lo=mu) : (hi=mu)
        nextfloat(lo) >= hi && break
    end
    isfinite(mu) && isfinite(derivative) && derivative > 0 ||
        throw(Infeasible("ill-conditioned or nonfinite mu fit"))
    abs(number-target) <= 2e-13*target || throw(Infeasible("mu fit precision not attained"))
    return (mu=mu,number=number,target=target,capacity=capacity,
            number_relative_error=(number-target)/target,derivative=derivative,
            mu_condition=target/derivative,bracket_width=hi-lo,iterations=iteration)
end
function fd_green(GR::Field, A::Field, condition, scale, energy, mu, kBT)
    GL,GG = similar(A),similar(A)
    maxpair = 0.0
    @inbounds for e in axes(A,1)
        f,hole = fermi_pair((energy[e]-mu)/kBT)
        maxpair = max(maxpair,abs(f+hole-1))
        for b in axes(A,4), a in axes(A,3), m in axes(A,2)
            GL[e,m,a,b] = im*f*A[e,m,a,b]
            GG[e,m,a,b] = -im*hole*A[e,m,a,b]
        end
    end
    return GreenState(GR,GL,GG,A,condition,scale),maxpair
end
function pair_A_error(green::GreenState)
    value = 0.0
    @inbounds @simd for i in eachindex(green.A)
        value += abs2(im*(green.Gᵍ[i]-green.Gˡ[i])-green.A[i])
    end
    return sqrt(value)
end

function local_retarded(c::Context, gamma::Field)
    g = c.problem.grids
    rawpart,roundoff = N._production_pi_fft_hilbert_transform(
        gamma,g.ε,g.wᴱ,c.cache.product_hilbert_operator)
    realpart = hermitian_copy(rawpart)
    correction = norms(rawpart,realpart,g.wᴱ,g.wᵏ)
    correction.relative <= 1e-10 && correction.weighted_relative <= 1e-10 ||
        throw(Infeasible("PI real-part symmetry correction exceeds original arithmetic target"))
    println("PI_HERMITIAN_CORRECTION ",correction); flush(stdout)
    @inbounds @simd for i in eachindex(gamma)
        realpart[i] -= im*gamma[i]/2
    end
    return realpart,roundoff
end
"R-only barrier: all dummy lesser/greater and cavity-side arrays die on return."
function _retarded_state(c::Context, gamma::Field)
    R,roundoff = local_retarded(c,gamma)
    z = zeros(ComplexF64,size(gamma))
    dummy = SelfEnergyFamily(R,z,z) # native cavity only reads these fields
    embedding,plus,minus = N._cavity_embedding_production(
        c.problem,dummy,c.h,c.options;plan=c.cavity)
    R .+= embedding.Σᴿ
    GR,condition,scale = N._retarded_green_production(
        c.problem.grids.ε,c.h,R;η=0.0,options=c.options)
    all(isfinite,GR) && all(isfinite,condition) && all(isfinite,scale) ||
        throw(Infeasible("nonfinite Dyson inverse/condition"))
    return GR,condition,scale,roundoff
end
retarded_state(c::Context,gamma::Field) = inverse_domain(()->_retarded_state(c,gamma))
function phi(c::Context, green::GreenState; keep_channels=false)
    p = c.problem
    output = zeros(ComplexF64,size(green.A))
    channels = Dict{Symbol,Field}()
    sizes = Dict{String,Float64}()
    for id in STATIC_IDS
        field = N._selected_static_contraction(p,c.cache,id,green.A,p.kernels.qᴷ[id],c.options)::Field
        output .+= field
        sizes[string(id)] = norm(field)
        keep_channels && (channels[id] = field)
    end
    lesser,greater = N._production_lo_contraction(
        c.cache.kernels[:LO],green,p,c.cache,p.kernels.qᴷ[:LO],c.options)
    @inbounds @simd for i in eachindex(lesser)
        lesser[i] = im*(greater[i]-lesser[i])
    end
    output .+= lesser
    sizes["LO"] = norm(lesser)
    keep_channels && (channels[:LO] = lesser)
    return output,channels,sizes
end

# At most 3 QR columns. All coefficients and weighted scalar products are REAL.
function aa_trial(states::Vector{Field}, residuals::Vector{Field}, we, wk, cfg)
    x,f = last(states),last(residuals)
    beta,alpha = Float64(cfg["aa_beta"]),Float64(cfg["linear_alpha"])
    linear() = x .+ alpha.*f
    length(states) < 2 && return linear(),(kind="linear",rank=0,condition=1.0,l1=1.0,coefficients=Float64[],firstpair=0,R_diagonal=Float64[],initial_scales=Float64[])
    for firstpair in 1:length(states)-1
        m = length(states)-firstpair
        columns = Field[]
        R = zeros(Float64,m,m)
        initial_scales = Float64[]
        valid = true
        for j in 1:m
            k = firstpair+j-1
            v = residuals[k+1] .- residuals[k]
            initial = wnorm(v,we,wk)
            push!(initial_scales,initial)
            for pass in 1:2, i in 1:j-1
                projection = inner(columns[i],v,we,wk)
                R[i,j] += projection
                v .-= projection.*columns[i]
            end
            R[j,j] = wnorm(v,we,wk)
            if initial == 0 || R[j,j] <= 64eps(Float64)*initial
                valid = false; break
            end
            v ./= R[j,j]
            push!(columns,v)
        end
        valid || continue
        condition = cond(R)
        isfinite(condition) && condition <= cfg["qr_condition_max"] || continue
        rhs = [inner(q,f,we,wk) for q in columns]
        coeff = UpperTriangular(R)\rhs
        # Expansion in past (x_i+beta*f_i); constrained affine coefficients.
        affine = zeros(Float64,m+1); affine[end] = 1.0
        for j in 1:m
            affine[j] += coeff[j]; affine[j+1] -= coeff[j]
        end
        l1 = sum(abs,affine)
        isfinite(l1) && l1 <= cfg["affine_l1_max"] || continue
        trial = x .+ beta.*f
        for j in 1:m
            k = firstpair+j-1
            # Fused broadcast; no retained DeltaGamma/DeltaF arrays here.
            d = coeff[j]
            xnew,xold,fnew,fold = states[k+1],states[k],residuals[k+1],residuals[k]
            @. trial -= d*((xnew-xold)+beta*(fnew-fold))
        end
        return trial,(kind="AAII",rank=m,condition=condition,l1=l1,coefficients=coeff,
            firstpair=firstpair,R_diagonal=diag(R),initial_scales=initial_scales)
    end
    return linear(),(kind="linear_rank_reset",rank=0,condition=Inf,l1=1.0,coefficients=Float64[],firstpair=0,R_diagonal=Float64[],initial_scales=Float64[])
end

function data_bytes(output)
    sum(filesize(joinpath(root,f)) for (root,_,files) in walkdir(output) for f in files;init=0)
end
function atomic_h5(f::Function, path::String, cfg; expected_bytes::Int)
    ispath(path) && error("Refusing to overwrite scientific artifact $path")
    expected_bytes <= cfg["single_file_max_bytes"] || error("file admission exceeds cap")
    data_bytes(dirname(path))+expected_bytes <= cfg["data_max_bytes"] || error("data admission exceeds cap; earlier raw preserved")
    HDF5.h5open(path*".tmp","w") do file
        f(file)
    end
    filesize(path*".tmp") <= cfg["single_file_max_bytes"] || error("file cap exceeded; raw tmp preserved")
    data_bytes(dirname(path)) <= cfg["data_max_bytes"] || error("data cap exceeded; raw tmp preserved")
    mv(path*".tmp",path)
end
function write_complex(group, name::String, array::Array{ComplexF64,D}) where D
    # HDF5 native hyperslab API, <=~8 MiB per real/imag writer buffer.
    shape = size(array)
    rows = min(shape[1],max(1,fld(8*1024^2,8*prod(shape[2:end]))))
    chunk = (rows,shape[2:end]...)
    for (suffix,part) in (("_real",real),("_imag",imag))
        dataset = HDF5.create_dataset(group,name*suffix,Float64,shape;chunk)
        try
            for firstrow in 1:rows:shape[1]
                indices = (firstrow:min(firstrow+rows-1,shape[1]),ntuple(_->Colon(),D-1)...)
                dataset[indices...] = part.(@view array[indices...])
            end
        finally
            close(dataset)
        end
    end
end
function gamma_checkpoint(output,name,gamma::Field,cfg,call)
    atomic_h5(joinpath(output,name*".h5"),cfg;expected_bytes=sizeof(gamma)+1024^2) do file
        file["schema"] = "campaign133-gamma-checkpoint-v1"
        file["identity"] = cfg["identity"]
        file["pi_calls_total"] = call
        file["restart_kind"] = "Gamma field only; fresh AA history, not native exact restart"
        write_complex(file,"Gamma",gamma)
    end
end
function write_green(file,green::GreenState)
    group = HDF5.create_group(file,"green")
    for (name,array) in (("GR",green.Gᴿ),("GL",green.Gˡ),("GG",green.Gᵍ),("A",green.A))
        write_complex(group,name,array)
    end
    group["condition_number"] = green.condition_number
    group["dyson_scale"] = green.dyson_scale
end
function write_metadata(group,name::String,value)
    # BitArray is packed and has no dense pointer/stride contract for HDF5.
    group[name] = value isa BitArray ? Int8.(value) : value isa Symbol ? string(value) : value
    return nothing
end
function save_model(c::Context,output,cfg)
    model_estimate = sum(sizeof(K) for K in values(c.problem.kernels.K))+32*1024^2
    atomic_h5(joinpath(output,"model.h5"),cfg;expected_bytes=model_estimate) do file
        p = c.problem
        file["schema"] = "campaign133-gamma-model-v2"
        file["mask_storage"] = "BitArray metadata materialized as dense Int8 {0,1}"
        file["algorithms_toml"] = sprint(io->TOML.print(io,Dict(string(k)=>
            (getfield(c.options.algorithms,k) isa Symbol ? string(getfield(c.options.algorithms,k)) : getfield(c.options.algorithms,k))
            for k in fieldnames(typeof(c.options.algorithms)))))
        for name in fieldnames(typeof(p.grids))
            write_metadata(file,"grid_"*string(name),getfield(p.grids,name))
        end
        for name in fieldnames(typeof(p.profiles))
            write_metadata(file,"profile_"*string(name),getfield(p.profiles,name))
        end
        for name in fieldnames(typeof(p.basis))
            value = getfield(p.basis,name)
            if value isa Array{ComplexF64}
                write_complex(file,"basis_"*string(name),value)
            else
                write_metadata(file,"basis_"*string(name),value)
            end
        end
        write_complex(file,"h_fixedU0",c.h)
        kernels = HDF5.create_group(file,"kernels")
        for id in p.kernels.enabled
            group = HDF5.create_group(kernels,string(id))
            group["qK"] = p.kernels.qᴷ[id]
            write_complex(group,"K",p.kernels.K[id])
        end
        # Actual sparse stencils are saved as their owning field values.
        shifts = HDF5.create_group(file,"shift_plans")
        for name in (:W̃₊ᴸᴼ,:W̃₋ᴸᴼ)
            plan = getfield(c.cache,name)
            group = HDF5.create_group(shifts,string(name))
            for key in fieldnames(typeof(plan))
                v = getfield(plan,key)
                write_metadata(group,string(key),v)
            end
        end
    end
end

function sparse_shift(plan)
    rows,columns,values = Int[],Int[],Float64[]
    for e in eachindex(plan.valid)
        plan.valid[e] || continue
        push!(rows,e); push!(columns,plan.left[e]); push!(values,1-plan.weight_right[e])
        push!(rows,e); push!(columns,plan.right[e]); push!(values,plan.weight_right[e])
    end
    return sparse(rows,columns,values,length(plan.valid),length(plan.valid))
end
function stencil_witness(c::Context,mu)
    g = c.problem.grids
    plus,minus = c.cache.W̃₊ᴸᴼ,c.cache.W̃₋ᴸᴼ
    quantum = Float64(ustrip(u"eV",c.problem.physical.ħωᴸᴼ))/c.problem.scales.E₀_eV
    bose = inv(expm1(quantum/c.kBT))
    maximum_defect,leakage = 0.0,0.0
    for (plan,direction) in ((plus,1),(minus,-1)),e in eachindex(g.ε)
        plan.valid[e] || continue
        fe,_ = fermi_pair((g.ε[e]-mu)/c.kBT)
        wr = plan.weight_right[e]
        leakage = max(leakage,min(abs(wr),abs(1-wr)))
        defect = 0.0
        for (j,w) in ((plan.left[e],1-wr),(plan.right[e],wr))
            fj,hj = fermi_pair((g.ε[j]-mu)/c.kBT)
            incoming = direction == 1 ? (bose+1)*fj : bose*fj
            outgoing = direction == 1 ? bose*hj : (bose+1)*hj
            defect += abs(w*(incoming-fe*(incoming+outgoing)))
        end
        maximum_defect = max(maximum_defect,defect)
    end
    wp,wm = sparse_shift(plus),sparse_shift(minus)
    W = spdiagm(0=>g.wᴱ)
    d = W*wp-transpose(wm)*W
    return Dict("phonon_scaled"=>quantum,"bose"=>bose,
        "shift_plus_scaled"=>plus.δ,"shift_minus_scaled"=>minus.δ,
        "maximum_minor_interpolation_weight"=>leakage,
        "independent_KMS_coefficient_absolute_bound"=>maximum_defect,
        "weighted_shift_adjoint_absolute_F"=>norm(d),
        "weighted_shift_adjoint_scale_F"=>norm(W*wp),
        "weighted_shift_adjoint_relative_F"=>relative(norm(d),norm(W*wp)),
        "meaning"=>"coefficient witness on actual nodes; KMS map measured separately; trapezoid endpoints may break weighted adjoint")
end

function kms_metrics(family::SelfEnergyFamily,gammaK::Field,g,mu::Float64,kBT::Float64)
    n2,s2,nw2,sw2 = 0.0,0.0,0.0,0.0
    maxabsolute,maxscale,maxrelative = 0.0,0.0,0.0
    worst_e,worst_m = 0,0
    @inbounds for m in axes(gammaK,2),e in axes(gammaK,1)
        f,hole = fermi_pair((g.ε[e]-mu)/kBT)
        bn,bs = 0.0,0.0
        for b in axes(gammaK,4),a in axes(gammaK,3)
            bn += abs2(family.Σˡ[e,m,a,b]-im*f*gammaK[e,m,a,b])
            bn += abs2(family.Σᵍ[e,m,a,b]+im*hole*gammaK[e,m,a,b])
            bs += abs2(family.Σˡ[e,m,a,b])+abs2(family.Σᵍ[e,m,a,b])
        end
        w = g.wᴱ[e]*g.wᵏ[m]
        n2 += bn; s2 += bs; nw2 += w*bn; sw2 += w*bs
        if sqrt(bn) > maxabsolute
            maxabsolute,maxscale,worst_e,worst_m = sqrt(bn),sqrt(bs),e,m
        end
        maxrelative = max(maxrelative,relative(sqrt(bn),sqrt(bs)))
    end
    return Dict("absolute"=>sqrt(n2),"scale"=>sqrt(s2),"relative"=>relative(sqrt(n2),sqrt(s2)),
        "weighted_absolute"=>sqrt(nw2),"weighted_scale"=>sqrt(sw2),
        "weighted_relative"=>relative(sqrt(nw2),sqrt(sw2)),
        "maximum_block_absolute"=>maxabsolute,"scale_at_maximum_block"=>maxscale,
        "energy_index"=>worst_e,"k_index"=>worst_m,"maximum_block_relative"=>maxrelative,
        "measure"=>"full-domain joint lesser/greater Frobenius; weighted uses wE*wk; absolute fields in scaled self-energy E0",
        "scope"=>"initial seed-spectrum consumer control at current mu, not a uniform KMS certificate over future Gamma/mu")
end
function channel_check_record(witness,kms::AbstractDict)
    # Scalar norm fields and a nested KMS table share this persisted record.
    record = Dict{String,Any}(string(k)=>v for (k,v) in pairs(witness))
    record["KMS"] = kms
    record["kms_relative"] = kms["relative"]
    return record
end
function compare_actual(c::Context,gamma::Field,green::GreenState,rawphi::Field,channels,
                        mu,output,cfg,budget)
    actual_call!(budget)
    sc,roundoff = N._scattering_candidate_production(c.problem,green,c.cache,c.options;return_roundoff=true)
    g = c.problem.grids
    broad = zeros(ComplexF64,size(gamma)); ret = similar(broad); fill!(ret,0)
    checks = Dict{String,Any}("native_roundoff"=>roundoff)
    checks["stencil"] = stencil_witness(c,mu)
    worst = 0.0
    for id in c.problem.kernels.enabled
        family = sc[id]
        field = im.*(family.Σᵍ.-family.Σˡ)
        witness = norms(channels[id],field,g.wᴱ,g.wᵏ)
        worst = max(worst,witness.relative)
        kms = kms_metrics(family,field,g,mu,c.kBT)
        checks[string(id)] = channel_check_record(witness,kms)
        worst = max(worst,kms["relative"])
        broad .+= field; ret .+= family.Σᴿ
    end
    reconstructed,_ = local_retarded(c,rawphi)
    checks["retarded_sum"] = Dict(string(k)=>v for (k,v) in pairs(norms(ret,reconstructed,g.wᴱ,g.wᵏ)))
    worst = max(worst,checks["retarded_sum"]["relative"])
    # Compare SAME input Gamma cavity R: dummy kinetic vs physical current-FD.
    inputR,_ = local_retarded(c,gamma)
    inputL,inputG = similar(gamma),similar(gamma)
    @inbounds for e in axes(gamma,1)
        f,hole = fermi_pair((g.ε[e]-mu)/c.kBT)
        for b in axes(gamma,4),a in axes(gamma,3),m in axes(gamma,2)
            inputL[e,m,a,b] = im*f*gamma[e,m,a,b]
            inputG[e,m,a,b] = -im*hole*gamma[e,m,a,b]
        end
    end
    actual_call!(budget)
    emb,ep,em = N._cavity_embedding_production(c.problem,SelfEnergyFamily(inputR,inputL,inputG),c.h,c.options;plan=c.cavity)
    inputR .+= emb.Σᴿ
    checkGR,_,_ = N._retarded_green_production(g.ε,c.h,inputR;η=0.0,options=c.options)
    checks["same_input_cavity_GR"] = Dict(string(k)=>v for (k,v) in pairs(norms(green.Gᴿ,checkGR,g.wᴱ,g.wᵏ)))
    worst = max(worst,checks["same_input_cavity_GR"]["relative"])
    checks["worst_relative"] = worst
    checks["passed"] = worst <= cfg["arithmetic_relative"]
    atomic_h5(joinpath(output,"actual-one-map-raw.h5"),cfg;expected_bytes=3*sizeof(gamma)+16*1024^2) do file
        file["schema"] = "campaign133-gamma-native-map-control-v1"
        write_complex(file,"raw_Phi",rawphi)
        write_complex(file,"native_Gamma_sum",broad)
        write_complex(file,"native_retarded_sum",ret)
    end
    open(joinpath(output,"actual-one-map.toml"),"w") do io; TOML.print(io,checks); end
    checks["passed"] || throw(Infeasible("actual native map equivalence/KMS failed: $worst"))
    return checks
end

"One Pi call; full Green/raw-Phi die at this barrier after optional persistence."
function evaluate(c::Context,gamma::Field,budget,output,cfg,streak; precheck=false)
    check_time(budget)
    t = time()
    inputcone = cone(gamma)
    actual_call!(budget;pi=true) # cone-admitted inputs only; failed map counts
    call = budget.pi
    GR,condition,scale,roundoff = retarded_state(c,gamma)
    A = spectral_function(GR)
    acone = cone(A;require_nonzero=false)
    q = N._seed_number_weights_production(GR,c.problem.grids,c.problem.physical.g_s;options=c.options)
    fitted = fit_mu(c.problem.grids.ε,q,c.target,c.kBT)
    green,pairerror = fd_green(GR,A,condition,scale,c.problem.grids.ε,fitted.mu,c.kBT)
    rawphi,channels,sizes = phi(c,green;keep_channels=precheck)
    # Save raw operator output before representation conversion if a cone fails.
    mapcone = try
        cone(rawphi)
    catch e
        e isa Infeasible || rethrow()
        gamma_checkpoint(output,"failed_raw_Phi_call$(call)",rawphi,cfg,call)
        rethrow()
    end
    precheck && compare_actual(c,gamma,green,rawphi,channels,fitted.mu,output,cfg,budget)
    empty!(channels)
    g = c.problem.grids
    mapped,residual,working_residual,representation = gamma_residuals(rawphi,gamma,g.wᴱ,g.wᵏ)
    newstreak = residual.relative < cfg["trigger_relative"] &&
                residual.weighted_relative < cfg["trigger_relative"] ? streak+1 : 0
    trigger = newstreak >= cfg["trigger_passes"]
    terminal = trigger || call == cfg["pi_calls_total"] || stage_expired(budget)
    stats = Dict{String,Any}("call"=>call,"seconds"=>time()-t,
        "mu"=>Dict(string(k)=>v for (k,v) in pairs(fitted)),
        "residual"=>Dict(string(k)=>v for (k,v) in pairs(residual)),
        "residual_definition"=>"RAW Phi-Gamma; raw Phi norm denominator; handoff trigger and best score",
        "AA_working_residual"=>Dict(string(k)=>v for (k,v) in pairs(working_residual)),
        "representation_conversion"=>Dict(string(k)=>v for (k,v) in pairs(representation)),
        "absolute_E0_eV"=>c.problem.scales.E₀_eV,
        "input_cone"=>Dict(string(k)=>v for (k,v) in pairs(inputcone)),
        "map_cone"=>Dict(string(k)=>v for (k,v) in pairs(mapcone)),
        "A_cone"=>Dict(string(k)=>v for (k,v) in pairs(acone)),
        "raw_phi_hermitization_absolute"=>representation.absolute,
        "fermi_pair_sum_error"=>pairerror,"green_pair_A_identity_absolute"=>pair_A_error(green),
        "maximum_Dyson_condition"=>maximum(condition),"minimum_Dyson_scale"=>minimum(scale),
        "roundoff"=>roundoff,
        "channel_Gamma_norms"=>sizes,"trigger"=>trigger)
    open(joinpath(output,"map-$(lpad(call,3,'0')).toml"),"w") do io; TOML.print(io,stats); end
    if terminal
        terminal_bytes = 6*sizeof(gamma)+sizeof(condition)+sizeof(scale)+16*1024^2
        atomic_h5(joinpath(output,"reduced-terminal.h5"),cfg;expected_bytes=terminal_bytes) do file
            file["schema"] = "campaign133-gamma-reduced-terminal-v1"
            file["mu"] = fitted.mu
            file["native_converged"] = false
            file["handoff_trigger"] = trigger
            write_complex(file,"Gamma",gamma)
            write_complex(file,"raw_Phi",rawphi)
            write_green(file,green)
        end
    end
    return mapped,mapped.-gamma,stats,newstreak,trigger
end

function record_attempt(output,call,kind,reason,budget;info=nothing)
    row = Dict{String,Any}("call"=>call,"kind"=>kind,"reason"=>reason,
        "elapsed_seconds"=>time()-budget.started)
    info !== nothing && (row["AA"] = Dict(string(k)=>v for (k,v) in pairs(info)))
    # One immutable TOML per attempt avoids nested-table append ambiguity.
    attempt_index = count(name->startswith(name,"attempt-"),readdir(output))+1
    open(joinpath(output,"attempt-$(lpad(attempt_index,3,'0')).toml"),"w") do io
        TOML.print(io,row); flush(io)
    end
end

"Returns only Gamma/small scalars. Histories, QR and all old Greens cannot escape."
function reduced_stage(c::Context,budget,output,cfg)
    p = c.problem
    gamma = load_warm_seed(c,cfg,output)
    gamma_checkpoint(output,"Gamma-first",gamma,cfg,0)
    mapped,F,stats,streak,trigger = evaluate(c,gamma,budget,output,cfg,0;precheck=true)
    record_attempt(output,budget.pi,"precheck_seed","accepted",budget)
    if stage_expired(budget)
        return gamma,false,Dict("status"=>"precheck_wall_cap","pi_calls"=>budget.pi,
            "precheck_actual_calls"=>budget.precheck,"last_map"=>stats)
    end
    budget.stage = :pilot; budget.stage_started = time()
    states,residuals = Field[gamma],Field[F]
    best = copy(gamma); bestnorm = stats["residual"]["relative"]
    status = "pilot_call_cap"
    try
    while budget.pi < cfg["pi_calls_total"] && !trigger
        if stage_expired(budget)
            status = "pilot_wall_cap"; break
        end
        trial,info = aa_trial(states,residuals,p.grids.wᴱ,p.grids.wᵏ,cfg)
        # AA rejection before Pi is recorded but consumes no Pi call. A map-domain
        # failure DOES consume a call. Exactly one convex fallback is permitted.
        fallback = false
        next = try
            evaluate(c,trial,budget,output,cfg,streak)
        catch e
            e isa Infeasible || rethrow()
            record_attempt(output,budget.pi,info.kind,sprint(showerror,e),budget;info)
            gamma_checkpoint(output,"rejected-trial-$(budget.pi)",trial,cfg,budget.pi)
            if info.kind != "AAII" || budget.pi >= cfg["pi_calls_total"]
                status = "infeasible_no_fallback"; break
            end
            empty!(states); empty!(residuals)
            alpha = Float64(cfg["linear_alpha"])
            trial = gamma .+ alpha.*F
            fallback = true
            try
                evaluate(c,trial,budget,output,cfg,0)
            catch fallback_error
                fallback_error isa Infeasible || rethrow()
                record_attempt(output,budget.pi,"convex_fallback",sprint(showerror,fallback_error),budget)
                gamma_checkpoint(output,"failed-fallback-$(budget.pi)",trial,cfg,budget.pi)
                status = "infeasible_fallback"; break
            end
        end
        gamma = trial
        mapped,F,stats,streak,trigger = next
        record_attempt(output,budget.pi,fallback ? "convex_fallback" : info.kind,"accepted",budget;info)
        push!(states,gamma); push!(residuals,F)
        while length(states) > cfg["history_differences"]+1
            popfirst!(states); popfirst!(residuals)
        end
        r = stats["residual"]["relative"]
        if r < bestnorm
            best = copy(gamma); bestnorm = r
        end
        if budget.pi == cfg["checkpoint_call"]
            gamma_checkpoint(output,"Gamma-intermediate",gamma,cfg,budget.pi)
        end
        println("Pi=",budget.pi," rGamma=",r," weighted=",stats["residual"]["weighted_relative"]," kind=",fallback ? "convex_fallback" : info.kind)
        flush(stdout)
    end
    finally
        # Also preserve the last valid field if a wall/API error interrupts AA
        # or the next map. A failed trial never replaces this accepted Gamma.
        gamma_checkpoint(output,"Gamma-best",best,cfg,budget.pi)
        gamma_checkpoint(output,"Gamma-last",gamma,cfg,budget.pi)
    end
    trigger && (status="reduced_trigger_native_pending")
    return gamma,trigger,Dict("status"=>status,"best_relative"=>bestnorm,
        "pi_calls"=>budget.pi,"precheck_actual_calls"=>budget.precheck,
        "last_map"=>stats)
end

function native_warm(c::Context,gamma::Field,strict)
    # Stage barrier after reduced history release; regenerate current physical FD.
    GR,condition,scale,_ = retarded_state(c,gamma)
    A = spectral_function(GR)
    q = N._seed_number_weights_production(GR,c.problem.grids,c.problem.physical.g_s;options=c.options)
    fit = fit_mu(c.problem.grids.ε,q,c.target,c.kBT)
    green,_ = fd_green(GR,A,condition,scale,c.problem.grids.ε,fit.mu,c.kBT)
    sc,_ = N._scattering_candidate_production(c.problem,green,c.cache,c.options;return_roundoff=true)
    emb,plus,minus = N._embedding_self_energy_selected(c.problem,green,c.cache,c.options;scattering=sc,h=c.h,cavity_plan=c.cavity)
    physical,lambda,total = N._current_green_production(c.problem,c.h,sc,emb,c.options)
    contract = N._solver_restart_contract(strict,c.options.algorithms)
    return SCBAResult(physical,sc,emb,plus,minus,SCBAIteration[],false,
        :equilibrium_gamma_warm_start,:unresolved,contract),lambda
end
function save_native(c::Context,result,output,cfg)
    native_bytes = 28*sizeof(result.green.A)+sizeof(result.green.condition_number)+
        sizeof(result.green.dyson_scale)+32*1024^2
    atomic_h5(joinpath(output,"native-terminal.h5"),cfg;expected_bytes=native_bytes) do file
        file["schema"] = "campaign133-equilibrium-gamma-native-v1"
        file["converged"] = result.converged
        file["status"] = string(result.status)
        file["quality"] = string(result.quality)
        write_green(file,result.green)
        families = copy(result.scattering)
        merge!(families,Dict(:embedding=>result.embedding,:embedding_plus=>result.embedding_plus,:embedding_minus=>result.embedding_minus))
        group = HDF5.create_group(file,"selfenergy")
        for (id,family) in families
            sub = HDF5.create_group(group,string(id))
            for (name,array) in (("retarded",family.Σᴿ),("lesser",family.Σˡ),("greater",family.Σᵍ))
                write_complex(sub,name,array)
            end
        end
        file["restart_contract_toml"] = sprint(io->TOML.print(io,result.restart_contract))
    end
end

"Signed collision quadrature. No absolute value is taken before summation."
function collision_metrics(family::SelfEnergyFamily,green::GreenState,g)
    particle,power,absolute,imaginary = 0.0,0.0,0.0,0.0
    @inbounds for m in axes(green.A,2),e in axes(green.A,1)
        incoming,outgoing = 0.0+0im,0.0+0im
        for b in axes(green.A,4),a in axes(green.A,3)
            incoming += family.Σˡ[e,m,a,b]*green.Gᵍ[e,m,b,a]
            outgoing += family.Σᵍ[e,m,a,b]*green.Gˡ[e,m,b,a]
        end
        w = g.wᴱ[e]*g.wᵏ[m]/(2π)
        delta = incoming-outgoing
        particle += w*real(delta)
        power += w*g.ε[e]*real(delta)
        absolute += w*(abs(incoming)+abs(outgoing))
        imaginary += w*imag(delta)
    end
    return Dict("signed_particle"=>particle,"signed_energy"=>power,
        "incoming_outgoing_absolute_scale"=>absolute,"signed_imaginary"=>imaginary,
        "particle_relative"=>relative(abs(particle),absolute),
        "measure"=>"wE*wk/(2pi); no spin factor, no hbar/physical rate conversion")
end
function retarded_identity(family::SelfEnergyFamily)
    numerator,scale = 0.0,0.0
    @inbounds for b in axes(family.Σᴿ,4),a in axes(family.Σᴿ,3),m in axes(family.Σᴿ,2),e in axes(family.Σᴿ,1)
        R = im*(family.Σᴿ[e,m,a,b]-conj(family.Σᴿ[e,m,b,a]))
        K = im*(family.Σᵍ[e,m,a,b]-family.Σˡ[e,m,a,b])
        numerator += abs2(R-K); scale += abs2(K)
    end
    return Dict("absolute"=>sqrt(numerator),"scale"=>sqrt(scale),
        "relative"=>relative(sqrt(numerator),sqrt(scale)))
end
function fdr_native(c::Context,green::GreenState,total::SelfEnergyFamily)
    g = c.problem.grids
    q = N._seed_number_weights_production(green.Gᴿ,g,c.problem.physical.g_s;options=c.options)
    fitted = fit_mu(g.ε,q,c.target,c.kBT) # diagnostic ONLY; native arrays unchanged
    nb = size(green.A,3)
    gr,sl,sg,tmp,rl,rg = (zeros(ComplexF64,nb,nb) for _ in 1:6)
    n2raw,n2pair,scale2,ident2 = 0.0,0.0,0.0,0.0
    nwraw,nwpair,sw,A2,Aw2 = 0.0,0.0,0.0,0.0,0.0
    occraw,occpair,occscale = 0.0,0.0,0.0
    occwraw,occwpair,occwscale = 0.0,0.0,0.0
    traw,tpair,tscale = 0.0,0.0,0.0
    @inbounds for m in axes(green.A,2),e in axes(green.A,1)
        f,hole = fermi_pair((g.ε[e]-fitted.mu)/c.kBT)
        w = g.wᴱ[e]*g.wᵏ[m]
        for b in 1:nb,a in 1:nb
            gr[a,b] = green.Gᴿ[e,m,a,b]
            sl[a,b] = total.Σˡ[e,m,a,b]; sg[a,b] = total.Σᵍ[e,m,a,b]
        end
        mul!(tmp,gr,sl); mul!(rl,tmp,adjoint(gr))
        mul!(tmp,gr,sg); mul!(rg,tmp,adjoint(gr))
        for b in 1:nb,a in 1:nb
            lf = im*f*green.A[e,m,a,b]; gf = -im*hole*green.A[e,m,a,b]
            olraw = abs2(rl[a,b]-lf); olpair = abs2(green.Gˡ[e,m,a,b]-lf)
            eraw = olraw+abs2(rg[a,b]-gf)
            epair = olpair+abs2(green.Gᵍ[e,m,a,b]-gf)
            s = abs2(lf)+abs2(gf)
            n2raw += eraw; n2pair += epair; scale2 += s
            nwraw += w*eraw; nwpair += w*epair; sw += w*s
            occraw += olraw; occpair += olpair; occscale += abs2(lf)
            occwraw += w*olraw; occwpair += w*olpair; occwscale += w*abs2(lf)
            A2 += abs2(green.A[e,m,a,b]); Aw2 += w*abs2(green.A[e,m,a,b])
            if g.trusted_energy[e]
                traw += w*eraw; tpair += w*epair; tscale += w*s
            end
            ident2 += abs2(im*(green.Gᵍ[e,m,a,b]-green.Gˡ[e,m,a,b])-green.A[e,m,a,b])
        end
    end
    return Dict("mu_fit"=>Dict(string(k)=>v for (k,v) in pairs(fitted)),
        "pair_sources"=>"raw=GR*stored_total_lesser/greater*GA before paired constraint; normalized=published native Green; fresh-candidate kinetic consistency is measured separately by native raw rK/lambda",
        "measure"=>"full-domain unweighted joint lesser/greater Frobenius; denominator sqrt(norm(ifA)^2+norm(-ihA)^2); absolute Green fields in scaled 1/E0",
        "weighted_measure"=>"same matrix Frobenius sums with wE*wk (dimensionless solver grid), no spin/2pi factor",
        "raw_pair_absolute"=>sqrt(n2raw),"normalized_pair_absolute"=>sqrt(n2pair),
        "fd_scale"=>sqrt(scale2),"raw_pair_relative"=>relative(sqrt(n2raw),sqrt(scale2)),
        "normalized_pair_relative"=>relative(sqrt(n2pair),sqrt(scale2)),
        "weighted"=>Dict("raw_absolute"=>sqrt(nwraw),"normalized_absolute"=>sqrt(nwpair),
            "fd_scale"=>sqrt(sw),"raw_relative"=>relative(sqrt(nwraw),sqrt(sw)),
            "normalized_relative"=>relative(sqrt(nwpair),sqrt(sw))),
        "occupied_lesser"=>Dict("raw_absolute"=>sqrt(occraw),"normalized_absolute"=>sqrt(occpair),
            "fd_scale"=>sqrt(occscale),"raw_relative"=>relative(sqrt(occraw),sqrt(occscale)),
            "normalized_relative"=>relative(sqrt(occpair),sqrt(occscale)),
            "weighted_raw_absolute"=>sqrt(occwraw),"weighted_normalized_absolute"=>sqrt(occwpair),
            "weighted_fd_scale"=>sqrt(occwscale),"weighted_raw_relative"=>relative(sqrt(occwraw),sqrt(occwscale)),
            "weighted_normalized_relative"=>relative(sqrt(occwpair),sqrt(occwscale))),
        "trusted_weighted_pair"=>Dict("raw_absolute"=>sqrt(traw),"normalized_absolute"=>sqrt(tpair),
            "fd_scale"=>sqrt(tscale),"raw_relative"=>relative(sqrt(traw),sqrt(tscale)),
            "normalized_relative"=>relative(sqrt(tpair),sqrt(tscale))),
        "A_scale"=>sqrt(A2),"normalized_pair_over_A"=>relative(sqrt(n2pair),sqrt(A2)),
        "weighted_A_scale"=>sqrt(Aw2),"weighted_normalized_pair_over_A"=>relative(sqrt(nwpair),sqrt(Aw2)),
        "normalized_pair_A_identity_absolute"=>sqrt(ident2),"projection_applied"=>false)
end

function component_metrics(stored::SelfEnergyFamily,fresh::SelfEnergyFamily,g)
    result = Dict{String,Any}()
    for (name,fresh_array,old) in (("retarded",fresh.Σᴿ,stored.Σᴿ),
                          ("lesser",fresh.Σˡ,stored.Σˡ),("greater",fresh.Σᵍ,stored.Σᵍ))
        # norms is a typed streaming Field barrier; no full difference temporary.
        result[name] = Dict(string(k)=>v for (k,v) in pairs(norms(fresh_array,old,g.wᴱ,g.wᵏ)))
    end
    return result
end

function raw_charge(green::GreenState,family::SelfEnergyFamily,g,spin)
    nb = size(green.A,3)
    gr,sl,tmp,out = (zeros(ComplexF64,nb,nb) for _ in 1:4)
    number = 0.0
    @inbounds for m in axes(green.A,2),e in axes(green.A,1)
        for b in 1:nb,a in 1:nb
            gr[a,b] = green.Gᴿ[e,m,a,b]; sl[a,b] = family.Σˡ[e,m,a,b]
        end
        mul!(tmp,gr,sl); mul!(out,tmp,adjoint(gr))
        for a in 1:nb
            number += spin*g.wᴱ[e]*g.wᵏ[m]*real(-im*out[a,a])/(2π)
        end
    end
    return number
end
function spectral_sum(green::GreenState,g)
    nb = size(green.A,3)
    block = zeros(ComplexF64,nb,nb)
    value = 0.0
    @inbounds for m in axes(green.A,2)
        fill!(block,0)
        for b in 1:nb,a in 1:nb,e in axes(green.A,1)
            block[a,b] += g.wᴱ[e]*green.A[e,m,a,b]/(2π)
        end
        for a in 1:nb; block[a,a] -= 1; end
        value = max(value,norm(block)/sqrt(nb))
    end
    return value
end
function audit_native(c::Context,result,output,cfg,budget)
    check_time(budget)
    budget.native < cfg["native_max_maps"] || error("native map cap exhausted before audit")
    budget.native += 1
    budget.native_events["audit_started"] += 1
    sc,roundoff = N._scattering_candidate_production(c.problem,result.green,c.cache,c.options;return_roundoff=true)
    emb,plus,minus = N._embedding_self_energy_selected(c.problem,result.green,c.cache,c.options;
        scattering=sc,h=c.h,cavity_plan=c.cavity)
    budget.native_events["audit_map_completed"] += 1
    total = N._total_family_production(result.scattering,result.embedding,c.problem.kernels.enabled,c.options)
    candidate_total = N._total_family_production(sc,emb,c.problem.kernels.enabled,c.options)
    charge = raw_charge(result.green,total,c.problem.grids,c.problem.physical.g_s)
    charge_candidate = raw_charge(result.green,candidate_total,c.problem.grids,c.problem.physical.g_s)
    workspace = ProductionResidualWorkspace(c.problem;worker_count=c.options.worker_count)
    measured = N._production_fixed_point_measurement(c.problem,c.h,result.green,total,candidate_total,
        result.scattering,result.embedding,result.embedding_plus,result.embedding_minus,
        sc,emb,plus,minus,charge/c.target,workspace;worker_count=c.options.worker_count)
    fresh = Dict{String,Any}(string(k)=>v for (k,v) in pairs(measured))
    fresh["raw_charge"] = charge; fresh["target_charge"] = c.target
    fresh["lambda"] = charge/c.target; fresh["candidate_lambda"] = charge_candidate/c.target
    fresh["r_lambda_candidate"] = abs(charge_candidate/c.target-1)
    fresh["normalized_charge"] = N._number_functional_production(result.green.Gˡ,c.problem.grids,c.problem.physical.g_s,c.options)
    fresh["roundoff"] = roundoff
    # Candidate Dyson/spectral consistency at the SAME published native Green.
    candidate_metrics = production_residual_suite(c.problem.grids.ε,c.h,result.green,candidate_total,candidate_total;
        workspace,worker_count=c.options.worker_count)
    fresh["r_D_candidate"] = candidate_metrics.r_D
    fresh["r_A_candidate"] = candidate_metrics.r_A
    report = Dict{String,Any}("fresh_native"=>fresh,"per_channel"=>Dict{String,Any}(),
        "not_a_Poisson_solution"=>true,"native_residual_definition"=>"unchanged production raw native gates, including their existing normalizations")
    families = Dict{Symbol,Tuple{SelfEnergyFamily,SelfEnergyFamily}}(
        id=>(result.scattering[id],sc[id]) for id in c.problem.kernels.enabled)
    families[:embedding] = (result.embedding,emb)
    families[:embedding_plus] = (result.embedding_plus,plus)
    families[:embedding_minus] = (result.embedding_minus,minus)
    for (id,(stored,fresh_family)) in families
        report["per_channel"][string(id)] = Dict("stored_collision"=>collision_metrics(stored,result.green,c.problem.grids),
            "fresh_collision"=>collision_metrics(fresh_family,result.green,c.problem.grids),
            "stored_retarded_kinetic"=>retarded_identity(stored),
            "fresh_retarded_kinetic"=>retarded_identity(fresh_family),
            "fresh_component_residuals"=>component_metrics(stored,fresh_family,c.problem.grids),
            "component_residual_definition"=>"fresh-stored absolute/relative; scale fresh; weighted wE*wk; zero/zero=0, nonzero/zero=Inf; diagnostic only, does not replace native max rSigma")
    end
    report["spectral_sum_native_residual"] = spectral_sum(result.green,c.problem.grids)
    try
        report["FDR"] = fdr_native(c,result.green,total)
    catch e
        e isa Infeasible || rethrow()
        report["FDR"] = Dict("not_measured"=>true,"reason"=>sprint(showerror,e))
    end
    open(joinpath(output,"native-fresh-audit.toml"),"w") do io; TOML.print(io,report); end
    budget.native_events["audit_report_written"] += 1
    return nothing
end

function native_stage(c::Context,gamma::Field,budget,output,cfg)
    cfg["native_max_maps"] == 32 || error("native structural cap must be1+30+1")
    budget.stage = :native; budget.stage_started = time()
    # Many rejected raw Gamma fields can consume the output budget. Preserve
    # them and reject handoff admission rather than deleting evidence or
    # allocating a full state which cannot be saved under the agreed cap.
    native_bytes = 28*sizeof(gamma)+2*8*size(gamma,1)*size(gamma,2)+48*1024^2
    if data_bytes(output)+native_bytes > cfg["data_max_bytes"]
        return Dict("not_measured"=>true,"reason"=>"native full-raw output admission exceeds6GiB")
    end
    # Two fresh five-channel maps are reserved: materialization and final audit.
    strict = SolverOptions(α_Σ=.2,max_scba=30,max_poisson=1,
        convergence=ConvergencePolicy(mode=:strict_fail_fast,minimum_scba_iterations=2,
            required_consecutive_scba_passes=3,stagnation_window=0,
            stagnation_relative_improvement=0.0,diagnostic_quality=DiagnosticQualityPolicy(enabled=false)))
    check_time(budget)
    budget.native += 1
    budget.native_events["materialization_started"] += 1
    warm,lambda = native_warm(c,gamma,strict)
    budget.native_events["materialization_completed"] += 1
    if stage_expired(budget)
        save_native(c,warm,output,cfg)
        return Dict("status"=>"native_materialization_wall_cap","converged"=>false,
            "quality"=>"unresolved","warm_lambda"=>lambda,"accounting"=>native_accounting(budget),
            "audit_not_measured"=>true)
    end
    fields = [f for f in fieldnames(SCBAIteration) if fieldtype(SCBAIteration,f)<:Real]
    native_options = N.with_production_options(c.options;
        checkpoint_request=(_->stage_expired(budget)))
    result = open(joinpath(output,"native-history.csv"),"w") do io
        println(io,join(["elapsed_s";string.(fields)],","))
        observer = function(row)
            budget.native += 1
            budget.native_events["solver_observed_completed"] += 1
            println(io,join([string(time()-budget.started);[string(getfield(row,f)) for f in fields]],",")); flush(io)
        end
        stop_at_aligned = aligned -> throw(NativeStop(aligned))
        try
            budget.native_events["solver_entered"] += 1
            solve_scba_production(c.problem,zeros(c.problem.numerical.N_z);options=strict,
                production_options=native_options,cache=c.cache,initial=warm,
                consume_initial=true,history_observer=observer,iteration_callback=stop_at_aligned)
        catch e
            e isa NativeStop || rethrow()
            # The callback does not return: no subsequent in-place mixing occurs.
            e.state
        end
    end
    save_native(c,result,output,cfg) # Preserve native state BEFORE fresh audit allocation.
    warm = nothing
    GC.gc(true) # release consumed warm Green/scratch before fresh full map
    if !stage_expired(budget)
        audit_native(c,result,output,cfg,budget)
    else
        open(joinpath(output,"native-audit-not-measured.toml"),"w") do io
            TOML.print(io,Dict("reason"=>"native or overall wall cap","not_measured"=>true))
        end
    end
    return Dict("status"=>string(result.status),"converged"=>result.converged,
        "quality"=>string(result.quality),"warm_lambda"=>lambda,"accounting"=>native_accounting(budget))
end

function build_context(cfg)
    physical = reference_parameters(V_period=0u"mV",Tᴸ=cfg["temperature_K"]*u"K",
        Tᴸᴼ=cfg["temperature_K"]*u"K",ΔV_alloy=.60u"eV",Ω₀=4.5168e-29u"m^3")
    numerical = NumericalParameters(N_z=cfg["Nz"],N_b=cfg["Nb"],P_basis=cfg["basis_periods"],
        E_min=cfg["energy_min_eV"]*u"eV",E_max=cfg["energy_max_eV"]*u"eV",N_E=cfg["NE"],
        M_E=cfg["trust_margin_meV"]*u"meV",k_max=cfg["kmax_per_nm"]*u"nm^-1",N_k=cfg["Nk"],
        N_φ=cfg["Nphi"],qz_max=cfg["qzmax_per_nm"]*u"nm^-1",N_qz=cfg["Nqz"],
        η_seed=cfg["seed_halfwidth_meV"]*u"meV")
    algorithms = AlgorithmOptions(kernel_build=:direct,embedding=:finite_chain,embedding_periods=4,
        mixing=:linear,energy_shift=:sparse_plan,hilbert=:product_integration,elastic_retarded=:hilbert)
    options = ProductionOptions(algorithms=algorithms,worker_count=cfg["workers"],hilbert_columns=1,
        memory_budget_bytes=cfg["memory_max_bytes"],parallel_backend=:threads,
        verify_fft_roundoff=true,checkpoint_every_scba=0,checkpoint_every_outer=0,progress_every_scba=1)
    built = build_configured_scattering_problem(physical=physical,numerical=numerical,scales=ScaleSystem(),
        scattering=ScatteringOptions(alloy=true),algorithms=algorithms,kernel_options=ProductionKernelOptions())
    p = built.problem
    Set(p.kernels.enabled) == Set((STATIC_IDS...,:LO)) || error("five-channel source identity mismatch")
    p.models.lo_population === :thermal || error("thermal LO required")
    h = project_hamiltonians(p,zeros(numerical.N_z))
    cache = build_production_cache(p;options)
    cache.product_hilbert_operator isa N._ProductionPIFFTHilbertPlan || error("approved uniform PI FFT plan required")
    cavity = N._cavity_production_plan(p,h,options)
    target = N._scaled_physics(p.physical,p.scales).Nᴰ²ᴰ
    kBT = Float64(ustrip(u"eV",CODATA.kᴮₑᵥ*p.physical.Tᴸ))/p.scales.E₀_eV
    return Context(p,cache,options,h,cavity,kBT,target)
end

include("warm_seed.jl")
end # module
