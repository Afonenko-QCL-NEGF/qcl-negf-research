using QCLNEGF, Unitful, LinearAlgebra, HDF5, TOML, SHA, Dates
const N = QCLNEGF.QCLNumerics
const SOURCE_SHA = "e6fffb997f36397f877aee1e2c98b1dbad7d97417d45458a042b5f074ce57a6d"
const INPUT_SHA = "0fa0d0bc2fc27602fa99b8fe0797071eae99c66ab2661dfde05bd00ebe19b86e"
const CORE_SHA = "2a500740a82dfd5b44f6be966ef63dd54010d097"
const LOCAL_NAMES = ["IFR", "LO", "acoustic", "alloy", "impurity"]
const MU_EV = 0.008495393439818052
sha_file(path) = bytes2hex(open(sha256, path))
sha_array(x) = bytes2hex(sha256(reinterpret(UInt8, vec(x))))
read_complex(group, name) = complex.(read(group[name*"_real"]), read(group[name*"_imag"]))
read_k(group, name, k) = complex.(group[name*"_real"][:,k,:,:], group[name*"_imag"][:,k,:,:])
function write_complex(group, name, x)
    group[name*"_real"] = real.(x)
    group[name*"_imag"] = imag.(x)
end

function parameters(c)
    p = reference_parameters(V_period=c["voltage_mV"]*u"mV", Tᴸ=c["temperature_K"]*u"K",
        Tᴸᴼ=c["temperature_K"]*u"K", ΔV_alloy=.60u"eV", Ω₀=4.5168e-29u"m^3")
    n = NumericalParameters(N_z=c["Nz"], N_b=c["Nb"], P_basis=c["basis_periods"],
        E_min=c["energy_min_eV"]*u"eV", E_max=c["energy_max_eV"]*u"eV", N_E=c["NE"],
        M_E=c["trust_margin_meV"]*u"meV", k_max=c["kmax_per_nm"]*u"nm^-1", N_k=c["Nk"],
        N_φ=c["Nphi"], qz_max=c["qzmax_per_nm"]*u"nm^-1", N_qz=c["Nqz"], η_seed=c["seed_meV"]*u"meV")
    return p, n, ScaleSystem()
end

function source_frame(file, c)
    p, n, s = parameters(c)
    g = build_grids(p,n,s)
    for (name, actual) in (("energy",g.ε),("weight_energy",g.wᴱ),("k",g.κ),
                           ("weight_k",g.wᵏ),("x",g.x),("weight_x",g.wˣ))
        read(file[name]) == actual || error("source axis differs: $name")
    end
    read(file["schema"]) == "campaign133-inner-raw-v3" || error("source is not v3")
    read(file["energy_scale_eV"]) == s.E₀_eV || error("source E0 differs")
    read(file["length_scale_m"]) == s.L₀_m || error("source L0 differs")
    read(file["spin_degeneracy"]) == p.g_s || error("source spin differs")
    all(iszero,read(file["hartree"])) || error("source Hartree differs")
    Set(keys(file["selfenergy"])) == Set([LOCAL_NAMES;"embedding";"embedding_plus";"embedding_minus"]) || error("source families differ")
    profiles = build_profiles(p,n,s,g)
    rebuilt = build_basis(p,n,s,g,profiles; localization=:pzp)
    mapping = Dict(:Φ=>"Phi",:χ=>"chi",:H₀=>"H0",:Z=>"Z",:M⁻¹=>"Minv",:T₊=>"Tplus",:T₋=>"Tminus")
    arrays = Dict(name=>read_complex(file["basis"], stored) for (name,stored) in mapping)
    comparison = Dict(string(name)=>Dict("source_sha256"=>sha_array(array),
        "rebuild_absolute_Frobenius"=>norm(array-getfield(rebuilt,name)),
        "rebuild_scale_Frobenius"=>norm(array)) for (name,array) in arrays)
    # Numerical source operators are bound exactly to v3. Only unpersisted
    # diagnostic metadata comes from the matching rebuilt basis container.
    basis = BasisData((get(arrays,name,getfield(rebuilt,name)) for name in fieldnames(BasisData))...)
    empty_kernels = KernelSet(Dict{Symbol,Array{ComplexF64,6}}(),Dict{Symbol,Float64}(),zeros(ComplexF64,0,0,0),Symbol[])
    sp = QCLNEGF._scaled_physics(p,s)
    shifts = map(δ->build_shift_plan(g.ε,δ),(sp.Eᵖ,-sp.Eᵖ,sp.ħωᴸᴼ,-sp.ħωᴸᴼ))
    container = NEGFProblem(p,n,ScatteringOptions(LO=false,acoustic=false,impurity=false,IFR=false,alloy=false),
        s,g,profiles,basis,empty_kernels,shifts...,:nodal_linear)
    h = project_hamiltonians(container,zeros(n.N_z))
    return container, h, comparison
end

function psd_witness(field, options)
    ne,nk,nb,_ = size(field)
    workers = N._production_worker_count(options,ne*nk)
    records = [Dict{String,Any}() for _ in 1:workers]
    N._production_worker_ranges!(ne*nk,workers) do worker, first_linear, last_linear
        minimum_eigenvalue, minimum_scaled, max_hermitian = Inf, Inf, 0.0
        worst = (1,1); block = zeros(ComplexF64,nb,nb); vector = zeros(ComplexF64,nb)
        residual, blockscale = NaN, NaN
        for linear in first_linear:last_linear
            e,k = (linear-1)%ne+1,(linear-1)÷ne+1
            B = Matrix(view(field,e,k,:,:)); scale = norm(B)
            max_hermitian = max(max_hermitian,norm(B-B')/max(scale,eps(Float64)))
            ev = eigen(Hermitian((B+B')/2))
            λ = ev.values[1]; scaled = λ/max(scale,eps(Float64))
            minimum_eigenvalue = min(minimum_eigenvalue,λ)
            if scaled < minimum_scaled
                minimum_scaled = scaled; worst = (e,k); block = B
                vector = ev.vectors[:,1]; blockscale = scale
                residual = norm(B*vector-λ*vector)/max(scale,eps(Float64))
            end
        end
        records[worker] = Dict("minimum_eigenvalue"=>minimum_eigenvalue,"minimum_scaled_eigenvalue"=>minimum_scaled,
            "max_relative_antihermitian"=>max_hermitian,"energy_index"=>worst[1],"k_index"=>worst[2],
            "block"=>block,"vector"=>vector,"witness_eigenvalue"=>real(dot(vector,block*vector)),
            "eigenvector_residual"=>residual,"block_scale"=>blockscale)
    end
    worst = records[argmin([x["minimum_scaled_eigenvalue"] for x in records])]
    worst["minimum_eigenvalue"] = minimum(x["minimum_eigenvalue"] for x in records)
    worst["max_relative_antihermitian"] = maximum(x["max_relative_antihermitian"] for x in records)
    return worst
end

function norm_comparison(delta, reference, g, f)
    # Fields are small enough to scan without constructing another delta array.
    ne,nk,nb,_ = size(reference)
    raw, weighted, occupied, scale, weighted_scale, occupied_scale = zeros(6)
    for k in 1:nk, e in 1:ne, b in 1:nb, a in 1:nb
        z = delta(e,k,a,b); r = reference[e,k,a,b]; w = g.wᴱ[e]*g.wᵏ[k]
        raw += abs2(z); scale += abs2(r)
        weighted += w*abs2(z); weighted_scale += w*abs2(r)
        occupied += w*f[e]^2*abs2(z); occupied_scale += w*f[e]^2*abs2(r)
    end
    return Dict("absolute_Frobenius"=>sqrt(raw),"reference_Frobenius"=>sqrt(scale),
        "relative_Frobenius"=>sqrt(raw)/max(sqrt(scale),eps(Float64)),
        "weighted_L2_absolute"=>sqrt(weighted),"weighted_L2_reference"=>sqrt(weighted_scale),
        "weighted_L2_relative"=>sqrt(weighted)/max(sqrt(weighted_scale),eps(Float64)),
        "occupied_L2_absolute"=>sqrt(occupied),"occupied_L2_reference"=>sqrt(occupied_scale),
        "occupied_L2_relative"=>sqrt(occupied)/max(sqrt(occupied_scale),eps(Float64)))
end

function integrals(A, g, f)
    ne,nk,nb,_ = size(A)
    M,F = zeros(ComplexF64,nk,nb,nb),zeros(ComplexF64,nk,nb,nb)
    for k in 1:nk,e in 1:ne,b in 1:nb,a in 1:nb
        value = g.wᴱ[e]/(2π)*A[e,k,a,b]
        M[k,a,b] += value; F[k,a,b] += f[e]*value
    end
    return M,F
end
charge(F,p,g) = p.g_s*sum(g.wᵏ[k]*real(tr(Matrix(view(F,k,:,:)))) for k in axes(F,1))
rsum(M) = maximum(norm(Matrix(view(M,k,:,:))-I)/sqrt(size(M,2)) for k in axes(M,1))

function fresh_response(localR,problem,h,options)
    zeros_K = zeros(ComplexF64,size(localR))
    local_family = SelfEnergyFamily(localR,zeros_K,zeros_K)
    plan = N._cavity_production_plan(problem,h,options)
    total,plus,minus = N._cavity_embedding_production(problem,local_family,h,options;plan)
    embeddingR = total.Σᴿ
    total = plus = minus = local_family = zeros_K = plan = nothing
    GC.gc()
    sigma = localR + embeddingR
    GR,condition,scale = N._retarded_green_production(problem.grids.ε,h,sigma;η=0.0,options)
    sigma = nothing
    A = similar(GR)
    N._production_parallel_linear!(size(GR,1)*size(GR,2),options) do linear
        e,k = (linear-1)%size(GR,1)+1,(linear-1)÷size(GR,1)+1
        for b in axes(GR,4), a in axes(GR,3)
            A[e,k,a,b] = im*(GR[e,k,a,b]-conj(GR[e,k,b,a]))
        end
    end
    return (;GR,A,embeddingR,max_condition=maximum(condition),max_scale=maximum(scale))
end

function edges(Γ,h,basis,geometry)
    ne,nk,nb,_ = size(Γ); periods = 4; count = 2periods+1
    left,right = zeros(ComplexF64,nk,nb,nb),zeros(ComplexF64,nk,nb,nb)
    coefficients = [N._product_integration_lag_coefficient(j) for j in 1:ne]
    for k in 1:nk,e in 1:ne,b in 1:nb,a in 1:nb
        left[k,a,b] -= coefficients[e]*Γ[e,k,a,b]
        right[k,a,b] += coefficients[ne+1-e]*Γ[e,k,a,b]
    end
    result = Dict{String,Any}()
    a = geometry.origin; step = geometry.spacing
    for (label,en,sigma) in (("left",a-step,left),("right",a+ne*step,right)),k in 1:nk
        D = zeros(ComplexF64,count*nb,count*nb)
        for cell in 1:count
            rows = (cell-1)*nb+1:cell*nb
            D[rows,rows] .= en*I-Matrix(view(h,k,:,:))-Matrix(view(sigma,k,:,:))
            if cell < count
                next = cell*nb+1:(cell+1)*nb
                D[rows,next] .= -basis.T₊; D[next,rows] .= -basis.T₊'
            end
        end
        hermitian = (D+D')/2
        ev = eigen(Hermitian(hermitian))
        result[label*"_k"*string(k)] = Dict("energy_dimensionless"=>en,"D"=>D,
            "eigenvalues"=>ev.values,"lambda_min"=>first(ev.values),"lambda_max"=>last(ev.values),
            "eigen_residual"=>norm(hermitian*ev.vectors-ev.vectors*Diagonal(ev.values))/max(norm(hermitian),eps(Float64)),
            "eigenvectors_orthogonality"=>norm(ev.vectors'*ev.vectors-I),
            "relative_antihermitian"=>norm(D-D')/max(norm(D),eps(Float64)))
    end
    return result
end

function main()
    BLAS.set_num_threads(1)
    length(ARGS)==1 || error("usage: frozen_control.jl NEW_OUTPUT_DIRECTORY")
    output = abspath(ARGS[1]); ispath(output) && error("output exists")
    mkpath(output); cp(@__FILE__,joinpath(output,"frozen_control.jl"))
    source = normpath(joinpath(@__DIR__,"../raw/S01-warm-FD64/state.h5"))
    input = normpath(joinpath(dirname(source),"input.toml"))
    identity = normpath(joinpath(dirname(source),"identity.toml"))
    marker = normpath(joinpath(dirname(source),"markers/final_markers.toml"))
    started = time(); stage = Ref("identity")
    checkstage(name) = begin
        stage[]=name
        time()-started <= 285 || error("phase-boundary time stop")
        bytes = sum(filesize(joinpath(root,file)) for (root,_,files) in walkdir(output) for file in files;init=0)
        bytes <= 32*1024^2 || error("output cap exceeded")
        println("stage=",name," elapsed_s=",time()-started); flush(stdout)
    end
    summary = Dict{String,Any}("status"=>"running","scope"=>"frozen retarded sensitivity; not SCBA, no Keldysh claim",
        "started_utc"=>string(Dates.now(Dates.UTC)),"source"=>source,"source_sha256"=>SOURCE_SHA,
        "script_sha256"=>sha_file(@__FILE__),"mu_eV"=>MU_EV,"Julia_threads"=>Threads.nthreads(),
        "BLAS_threads"=>BLAS.get_num_threads(),"current_stage"=>"identity")
    try
        sha_file(source)==SOURCE_SHA || error("source SHA mismatch")
        sha_file(input)==INPUT_SHA || error("source input SHA mismatch")
        cfg = TOML.parsefile(input); ident = TOML.parsefile(identity); markers = TOML.parsefile(marker)
        ident["input_sha256"]==INPUT_SHA || error("source identity differs")
        ident["core_head"]=="f50fa9bdfc04e299f82b6260f70fd8f829d838f3" || error("source Core identity differs")
        ident["runner_head"]=="18f45edb3bd086ebc7d5c984d68266a937358d4e" || error("source Runner identity differs")
        source_script = joinpath(dirname(source),"run_inner.jl")
        sha_file(source_script)==ident["script_sha256"] || error("source writer identity differs")
        summary["source_writer_sha256"] = sha_file(source_script)
        markers["equilibrium_mu_eV"]==MU_EV || error("source mu differs")
        cfg["embedding"]=="finite_chain" && cfg["embedding_periods"]==4 || error("embedding identity differs")
        cfg["voltage_mV"]==0 && cfg["temperature_K"]==70 && cfg["alloy"] || error("source physical identity differs")
        root = dirname(dirname(pathof(QCLNEGF)))
        head = readchomp(`git -C $root rev-parse HEAD`)
        head==CORE_SHA || error("Core candidate differs: $head")
        isempty(readchomp(`git -C $root status --porcelain`)) || error("Core candidate dirty")
        project = Base.active_project(); manifest=joinpath(dirname(project),"Manifest.toml")
        merge!(summary,Dict("input_sha256"=>sha_file(input),"identity_sha256"=>sha_file(identity),
            "marker_sha256"=>sha_file(marker),"core_head"=>head,"core_source"=>pathof(QCLNEGF),
            "core_operators_sha256"=>sha_file(joinpath(root,"src/numerics/optimized/production_operators.jl")),
            "core_fft_sha256"=>sha_file(joinpath(root,"src/numerics/optimized/fft_parallel.jl")),
            "project_sha256"=>sha_file(project),"manifest_sha256"=>sha_file(manifest),"julia_version"=>string(VERSION)))
        checkstage("source_frame")
        problem,h,basis_comparison = h5open(source,"r") do file; source_frame(file,cfg); end
        g,p,s,n = problem.grids,problem.physical,problem.scales,problem.numerical
        iszero(QCLNEGF._scaled_physics(p,s).Eᵖ) || error("source bias is nonzero")
        summary["basis_source_and_rebuild_comparison"] = basis_comparison
        summary["basis_metadata_note"] = "7 numerical operators from source v3; unpersisted diagnostic metadata from rebuilt basis; no full static/scattering acceptance"
        summary["h_sha256"] = sha_array(h)
        summary["Tminus_adjoint_defect"] = norm(problem.basis.T₋-problem.basis.T₊')
        options = ProductionOptions(worker_count=4,hilbert_columns=1,memory_budget_bytes=cfg["memory_budget_bytes"],
            algorithms=AlgorithmOptions(embedding=:finite_chain,embedding_periods=4,hilbert=:fft))
        pi_options = QCLNEGF.with_production_options(options;algorithms=AlgorithmOptions(embedding=:finite_chain,embedding_periods=4,hilbert=:product_integration))
        # Retain the source full-five-channel numerical budget guard; the
        # harness does not allocate microscopic kernels or ProductionCache.
        for (label,opt) in (("nodal",options),("PI",pi_options))
            estimate=estimate_production_memory(n,5;dense_mechanism_count=3,lo_kernel=:dense,options=opt,energy_nodes=g.ε)
            estimate.peak_bytes<=opt.memory_budget_bytes || error("source numerical budget rejected $label")
            summary[label*"_full_scientific_peak_estimate_bytes"] = estimate.peak_bytes
        end
        summary["source_numerical_budget_bytes"] = cfg["memory_budget_bytes"]
        target = QCLNEGF._scaled_physics(p,s).Nᴰ²ᴰ
        kBT = QCLNEGF._electronvolts(CODATA.kᴮₑᵥ*p.Tᴸ)/s.E₀_eV
        f = [QCLNEGF._fermi(e,MU_EV/s.E₀_eV,kBT) for e in g.ε]
        shape=(n.N_E,n.N_k,n.N_b,n.N_b); sourceR=zeros(ComplexF64,shape)
        h5open(source,"r") do file
            for k in 1:n.N_k,name in LOCAL_NAMES
                sourceR[:,k,:,:] .+= read_k(file["selfenergy/"*name],"retarded",k)
            end
        end
        Γ=similar(sourceR)
        for k in 1:n.N_k,e in 1:n.N_E,b in 1:n.N_b,a in 1:n.N_b
            Γ[e,k,a,b]=im*(sourceR[e,k,a,b]-conj(sourceR[e,k,b,a]))
        end
        summary["frozen_local_retarded_sha256"] = sha_array(sourceR)
        summary["frozen_local_gamma_sha256"] = sha_array(Γ)
        summary["families_summed_once"] = LOCAL_NAMES
        summary["families_excluded"] = ["embedding","embedding_plus","embedding_minus"]
        checkstage("Hilbert_reconstruction")
        nodal_plan=ProductionFFTHilbertPlan(g.ε,n.N_k*n.N_b^2;column_chunk=1,worker_count=4,verify_roundoff=true)
        pi_plan=N._ProductionPIFFTHilbertPlan(g.ε,n.N_k*n.N_b^2;column_chunk=1,worker_count=4)
        Rnodal,nodal_roundoff=production_fft_hilbert_transform(Γ,g.ε,g.wᴱ,nodal_plan;return_roundoff=true)
        Rpi,pi_roundoff=N._production_pi_fft_hilbert_transform(Γ,g.ε,g.wᴱ,pi_plan)
        summary["PI_geometry"] = Dict(string(k)=>v for (k,v) in pairs(pi_plan.geometry))
        summary["PI_minus_nodal_Lambda"] = norm_comparison((e,k,a,b)->Rpi[e,k,a,b]-Rnodal[e,k,a,b],Rnodal,g,f)
        summary["nodal_roundoff"] = nodal_roundoff; summary["PI_roundoff"] = pi_roundoff
        summary["nodal_plan_explicit_bytes"] = production_fft_workspace_bytes(nodal_plan)+production_fft_workspace_bytes(nodal_plan;verification=true)
        summary["PI_plan_explicit_bytes"] = production_fft_workspace_bytes(pi_plan)+production_fft_workspace_bytes(pi_plan;verification=true)
        function hermitian_defect(X)
            sqrt(sum(abs2(X[e,k,a,b]-conj(X[e,k,b,a])) for e in axes(X,1),k in axes(X,2),a in axes(X,3),b in axes(X,4)))/max(norm(X),eps(Float64))
        end
        summary["nodal_Lambda_Hermiticity"] = hermitian_defect(Rnodal)
        summary["PI_Lambda_Hermiticity"] = hermitian_defect(Rpi)
        Rnodal .-= 0.5im .* Γ; Rpi .-= 0.5im .* Γ
        summary["nodal_rebuild_minus_source_R"] = norm_comparison((e,k,a,b)->Rnodal[e,k,a,b]-sourceR[e,k,a,b],sourceR,g,f)
        for (label,R) in (("nodal",Rnodal),("PI",Rpi))
            summary[label*"_reconstructed_gamma_minus_frozen"] = norm_comparison(
                (e,k,a,b)->im*(R[e,k,a,b]-conj(R[e,k,b,a]))-Γ[e,k,a,b],Γ,g,f)
        end
        sourceR=nodal_plan=pi_plan=nothing; GC.gc()
        checkstage("Gamma_PSD_and_ideal_edges")
        gamma_psd=psd_witness(Γ,options)
        edge_witness=edges(Γ,h,problem.basis,N._uniform_product_integration_geometry(g.ε))
        checkstage("fresh_nodal_cavity_and_Dyson")
        nodal=fresh_response(Rnodal,problem,h,options)
        checkstage("fresh_PI_cavity_and_Dyson")
        pi=fresh_response(Rpi,problem,h,pi_options)
        checkstage("spectral_integrals_and_witnesses")
        MN,FN=integrals(nodal.A,g,f); MP,FP=integrals(pi.A,g,f)
        MS,FS=zeros(ComplexF64,n.N_k,n.N_b,n.N_b),zeros(ComplexF64,n.N_k,n.N_b,n.N_b)
        archived_raw,archived_scale=0.0,0.0
        h5open(source,"r") do file
            for k in 1:n.N_k
                A=read_k(file["green"],"A",k)
                for e in 1:n.N_E,b in 1:n.N_b,a in 1:n.N_b
                    z=A[e,a,b]; w=g.wᴱ[e]/(2π)
                    MS[k,a,b]+=w*z; FS[k,a,b]+=w*f[e]*z
                    archived_raw+=g.wᴱ[e]*g.wᵏ[k]*abs2(nodal.A[e,k,a,b]-z)
                    archived_scale+=g.wᴱ[e]*g.wᵏ[k]*abs2(z)
                end
            end
        end
        summary["fresh_nodal_minus_archived_A_weighted_L2_numerator"] = sqrt(archived_raw)
        summary["archived_A_weighted_L2_scale"] = sqrt(archived_scale)
        summary["PI_minus_nodal_A"] = norm_comparison((e,k,a,b)->pi.A[e,k,a,b]-nodal.A[e,k,a,b],nodal.A,g,f)
        summary["PI_minus_nodal_GR"] = norm_comparison((e,k,a,b)->pi.GR[e,k,a,b]-nodal.GR[e,k,a,b],nodal.GR,g,f)
        boundsM,boundsF=zeros(n.N_k),zeros(n.N_k)
        maximum_difference,maximum_occupied_difference = -Inf,-Inf
        worst_difference,worst_occupied_difference = (1,1),(1,1)
        for k in 1:n.N_k,e in 1:n.N_E
            difference=norm(Matrix(view(pi.A,e,k,:,:))-Matrix(view(nodal.A,e,k,:,:)))
            boundsM[k]+=g.wᴱ[e]/(2π)*difference; boundsF[k]+=g.wᴱ[e]/(2π)*f[e]*difference
            if difference > maximum_difference
                maximum_difference=difference;worst_difference=(e,k)
            end
            if f[e]*difference > maximum_occupied_difference
                maximum_occupied_difference=f[e]*difference;worst_occupied_difference=(e,k)
            end
        end
        summary["PI_minus_nodal_discrete_N_absolute_bound_over_target"] = p.g_s*sqrt(n.N_b)*sum(g.wᵏ.*boundsF)/target
        summary["PI_minus_nodal_M_absolute_Frobenius"] = norm(MP-MN)
        summary["PI_minus_nodal_N_over_target"] = (charge(FP,p,g)-charge(FN,p,g))/target
        summary["target_dimensionless"] = target; summary["E0_eV"] = s.E₀_eV
        summary["integral_measure"] = "source positive trapezoid wE/(2pi), source wk; no number normalization or mu refit"
        summary["norm_units"] = "all arrays/weights dimensionless with source E0,L0; unweighted and positive wE*wk L2 sums; occupied norm uses f^2; relative floors eps(Float64) in these scaled norms"
        for (label,response,M,F) in (("fresh_nodal",nodal,MN,FN),("fresh_PI",pi,MP,FP))
            summary[label] = Dict("r_sum"=>rsum(M),"N_fA_over_target"=>charge(F,p,g)/target,
                "max_condition_number"=>response.max_condition,"max_Dyson_scale"=>response.max_scale)
        end
        summary["archived_source_A"] = Dict("r_sum"=>rsum(MS),"N_fA_over_target"=>charge(FS,p,g)/target)
        nodal_psd=psd_witness(nodal.A,options); pi_psd=psd_witness(pi.A,options)
        h5open(joinpath(output,"integrals_and_witnesses.h5"),"w") do file
            file["schema"]="campaign133-frozen-local-gamma-retarded-control-v1"
            file["energy_dimensionless"]=g.ε; file["weight_energy_dimensionless"]=g.wᴱ
            file["k_dimensionless"]=g.κ; file["weight_k_dimensionless"]=g.wᵏ; write_complex(file,"h",h)
            file["per_k_PI_minus_nodal_M_triangle_bound"]=boundsM
            file["per_k_PI_minus_nodal_F_triangle_bound"]=boundsF
            for (label,M,F) in (("fresh_nodal",MN,FN),("fresh_PI",MP,FP),("archived_source_A",MS,FS))
                group=create_group(file,label);write_complex(group,"M",M);write_complex(group,"F",F)
            end
            for (label,witness) in (("Gamma_PSD",gamma_psd),("nodal_A_PSD",nodal_psd),("PI_A_PSD",pi_psd))
                group=create_group(file,label)
                write_complex(group,"block",pop!(witness,"block"));write_complex(group,"vector",pop!(witness,"vector"))
                summary[label]=witness
            end
            blockgroup=create_group(file,"response_blocks")
            for (label,(e,k)) in (("max_A_difference",worst_difference),
                                  ("max_fA_difference",worst_occupied_difference))
                block=create_group(blockgroup,label)
                block["energy_index"]=e;block["k_index"]=k;block["energy_dimensionless"]=g.ε[e]
                block["energy_eV_relative_to_Eref"]=g.ε[e]*s.E₀_eV;block["k_dimensionless"]=g.κ[k]
                for (name,array) in (("Gamma",Γ),("local_R_nodal",Rnodal),("local_R_PI",Rpi),
                    ("embedding_R_nodal",nodal.embeddingR),("embedding_R_PI",pi.embeddingR),
                    ("GR_nodal",nodal.GR),("GR_PI",pi.GR),("A_nodal",nodal.A),("A_PI",pi.A))
                    write_complex(block,name,Matrix(view(array,e,k,:,:)))
                end
            end
            group=create_group(file,"ideal_ghost_edge_D45_witnesses")
            for (label,witness) in edge_witness
                edge=create_group(group,label);write_complex(edge,"D",pop!(witness,"D"));edge["eigenvalues"]=pop!(witness,"eigenvalues")
                summary["edge_"*label]=witness
            end
        end
        summary["edge_scope"]="floating full9cell45block witnesses; no interval error bounds; original window excludes ghost strips; dark interior atoms not excluded; no pole attribution from failed sign"
        summary["source_scope"]="unconverged FD64 iterate; fresh nodal minus archived includes embedding/history alignment, PI minus fresh nodal is frozen operator sensitivity"
        summary["Keldysh_SCBA_Poisson_status"]="not_measured"
        checkstage("complete")
        summary["status"]="completed"
    catch error
        summary["status"]="failed";summary["error"]=sprint(showerror,error)
        rethrow()
    finally
        summary["current_stage"]=stage[];summary["elapsed_seconds"]=time()-started
        summary["finish_utc"]=string(Dates.now(Dates.UTC))
        open(joinpath(output,"summary.toml"),"w") do io;TOML.print(io,summary);end
    end
    println(summary)
end
main()
