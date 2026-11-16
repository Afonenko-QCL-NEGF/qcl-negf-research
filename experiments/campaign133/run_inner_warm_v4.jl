using QCLNEGF, QCLNEGFRunner, Unitful, LinearAlgebra, TOML, SHA, HDF5, Dates

length(ARGS) == 2 || error("usage: run_inner.jl INPUT.toml NEW_OUTPUT_DIRECTORY")
input_path, output = abspath.(ARGS)
ispath(output) && error("Refusing to overwrite existing run: $output")
mkpath(output)
cp(input_path, joinpath(output, "input.toml"))
cp(@__FILE__, joinpath(output, "run_inner.jl"))
cfg = TOML.parsefile(input_path)
BLAS.set_num_threads(1)
started = time()
stamp() = string(Dates.now(Dates.UTC)) * "Z"
function checkpoint_inner(path, problem, result)
    HDF5.h5open(path * ".tmp", "w") do file
        file["schema"] = "campaign133-inner-raw-v3"
        file["status"] = string(result.status)
        file["converged"] = result.converged
        file["quality"] = string(result.quality)
        for (name, array) in (("energy", problem.grids.ε), ("k", problem.grids.κ), ("x", problem.grids.x), ("weight_energy", problem.grids.wᴱ), ("weight_k", problem.grids.wᵏ), ("weight_x", problem.grids.wˣ))
            file[name] = array
        end
        file["hartree"] = zeros(problem.numerical.N_z)
        QCLNEGFRunner._write_physical_markers!(file, result.history)
        function complex_array(group, name, array)
            group[name * "_real"] = real.(array)
            group[name * "_imag"] = imag.(array)
        end
        file["energy_scale_eV"] = problem.scales.E₀_eV
        file["length_scale_m"] = problem.scales.L₀_m
        file["spin_degeneracy"] = problem.physical.g_s
        file["trusted_energy"] = Int8.(problem.grids.trusted_energy)
        basis = HDF5.create_group(file,"basis")
        for (name,array) in (("Phi",problem.basis.Φ),("chi",problem.basis.χ),("H0",problem.basis.H₀),("Z",problem.basis.Z),("Minv",problem.basis.M⁻¹),("Tplus",problem.basis.T₊),("Tminus",problem.basis.T₋))
            complex_array(basis,name,array)
        end
        kernels=HDF5.create_group(file,"kernels")
        for name in problem.kernels.enabled
            group=HDF5.create_group(kernels,string(name))
            complex_array(group,"normalized_K",problem.kernels.K[name])
            group["qK"] = problem.kernels.qᴷ[name]
        end
        green = HDF5.create_group(file, "green")
        for (name, array) in (("GR", result.green.Gᴿ), ("GL", result.green.Gˡ), ("GG", result.green.Gᵍ), ("A", result.green.A))
            complex_array(green, name, array)
        end
        green["condition_number"] = result.green.condition_number
        green["dyson_scale"] = result.green.dyson_scale
        families = copy(result.scattering)
        merge!(families, Dict(:embedding=>result.embedding, :embedding_plus=>result.embedding_plus, :embedding_minus=>result.embedding_minus))
        selfenergy = HDF5.create_group(file, "selfenergy")
        for (name, family) in families
            group = HDF5.create_group(selfenergy, string(name))
            for (component, array) in (("retarded", family.Σᴿ), ("lesser", family.Σˡ), ("greater", family.Σᵍ))
                complex_array(group, component, array)
            end
        end
        witnesses = HDF5.create_group(file, "witnesses")
        for row in result.history
            witness = row.witness
            witness === nothing && continue
            group = HDF5.create_group(witnesses, string(row.ν))
            for field in fieldnames(typeof(witness))
                value = getfield(witness, field)
                value isa Symbol && (value = string(value))
                value isa Dict && continue
                if value isa AbstractArray{<:Complex}
                    complex_array(group, string(field), value)
                else
                    group[string(field)] = value
                end
            end
            blocks = HDF5.create_group(group, "matrices")
            for (name, array) in witness.matrices
                complex_array(blocks, name, array)
            end
        end
    end
    mv(path * ".tmp", path)
end

function load_initial_guess(problem,production,cfg)
    haskey(cfg,"initial_state") || return nothing
    path=normpath(joinpath(dirname(input_path),cfg["initial_state"]))
    digest=bytes2hex(open(sha256,path))
    digest==cfg["initial_state_sha256"] || error("initial source hash differs")
    mode=cfg["initialization_mode"]
    mode in ("unmodified","FD_repair_from_retarded") || error("unknown initial mode")
    metadata=Dict{String,Any}("source_state_sha256"=>digest,"mode"=>mode,
        "resume_iterations"=>false,"history_reset"=>true,"consume_private_arrays"=>true,
        "scope"=>"new warm initial guess, not exact restart; first unmodified boundary repeats source iteration; FD first lambda is constructive identity")
    result=HDF5.h5open(path,"r") do file
        read(file["schema"]) in ("campaign133-inner-raw-v2","campaign133-inner-raw-v3") || error("source schema differs")
        read(file["energy"])==problem.grids.ε || error("initial energies differ")
        read(file["k"])==problem.grids.κ || error("initial momentum differs")
        read(file["weight_energy"])==problem.grids.wᴱ || error("initial energy weights differ")
        read(file["weight_k"])==problem.grids.wᵏ || error("initial momentum weights differ")
        all(iszero,read(file["hartree"])) || error("initial Hartree is not zero")
        readcomplex(group,name)=complex.(read(group[name*"_real"]),read(group[name*"_imag"]))
        group=file["green"]
        green=GreenState(readcomplex(group,"GR"),readcomplex(group,"GL"),readcomplex(group,"GG"),readcomplex(group,"A"),read(group["condition_number"]),read(group["dyson_scale"]))
        expected=Set(string.([problem.kernels.enabled;:embedding;:embedding_plus;:embedding_minus]))
        Set(keys(file["selfenergy"]))==expected || error("initial mechanisms differ")
        families=Dict{Symbol,SelfEnergyFamily}()
        for name in keys(file["selfenergy"])
            g=file["selfenergy/"*name]
            families[Symbol(name)]=SelfEnergyFamily(readcomplex(g,"retarded"),readcomplex(g,"lesser"),readcomplex(g,"greater"))
        end
        if mode=="FD_repair_from_retarded"
            iszero(problem.physical.F_bias) && problem.physical.Tᴸ==problem.physical.Tᴸᴼ || error("thermal init requires zero bias and common bath temperature")
            g=problem.grids;n=problem.numerical;s=problem.scales
            target=QCLNEGF._scaled_physics(problem.physical,s).Nᴰ²ᴰ
            q=[problem.physical.g_s/(2π)*g.wᴱ[e]*sum(g.wᵏ[m]*sum(real(green.A[e,m,a,a]) for a=1:n.N_b) for m=1:n.N_k) for e=1:n.N_E]
            minimum(q)>=0 && sum(q)>=target || error("invalid initial spectral capacity")
            kBT=QCLNEGF._electronvolts(CODATA.kᴮₑᵥ*problem.physical.Tᴸ)/s.E₀_eV
            number(mu)=sum(q[e]*QCLNEGF._fermi(g.ε[e],mu,kBT) for e=1:n.N_E)
            lo=first(g.ε)-50kBT;hi=last(g.ε)+50kBT
            number(lo)<=target<=number(hi) || error("initial mu not bracketed")
            for j=1:100
                mid=(lo+hi)/2
                number(mid)<target ? (lo=mid) : (hi=mid)
            end
            mu=(lo+hi)/2
            f=[QCLNEGF._fermi(e,mu,kBT) for e in g.ε]
            for family in values(families)
                for m=1:n.N_k,e=1:n.N_E,b=1:n.N_b,a=1:n.N_b
                    gamma=im*(family.Σᴿ[e,m,a,b]-conj(family.Σᴿ[e,m,b,a]))
                    family.Σˡ[e,m,a,b]=im*f[e]*gamma
                    family.Σᵍ[e,m,a,b]=-im*(1-f[e])*gamma
                end
            end
            for m=1:n.N_k,e=1:n.N_E,b=1:n.N_b,a=1:n.N_b
                green.Gˡ[e,m,a,b]=im*f[e]*green.A[e,m,a,b]
                green.Gᵍ[e,m,a,b]=-im*(1-f[e])*green.A[e,m,a,b]
            end
            metadata["common_mu_eV"]=mu*s.E₀_eV
            metadata["FD_number_over_target"]=number(mu)/target
            metadata["retarded_arrays_changed"]=false
        end
        scattering=Dict(name=>families[name] for name in problem.kernels.enabled)
        SCBAResult(green,scattering,families[:embedding],families[:embedding_plus],families[:embedding_minus],SCBAIteration[],false,:max_iterations,:unresolved)
    end
    open(joinpath(output,"initialization.toml"),"w") do io
        TOML.print(io,metadata)
    end
    println("Initialization=",metadata)
    flush(stdout)
    return result
end

alloy = get(cfg, "alloy", false)
physical = reference_parameters(V_period=cfg["voltage_mV"]*u"mV", Tᴸ=cfg["temperature_K"]*u"K", Tᴸᴼ=cfg["temperature_K"]*u"K", ΔV_alloy=alloy ? .60u"eV" : nothing, Ω₀=alloy ? 4.5168e-29u"m^3" : nothing)
numerical = NumericalParameters(N_z=cfg["Nz"], N_b=cfg["Nb"], P_basis=cfg["basis_periods"], E_min=cfg["energy_min_eV"]*u"eV", E_max=cfg["energy_max_eV"]*u"eV", N_E=cfg["NE"], M_E=cfg["trust_margin_meV"]*u"meV", k_max=cfg["kmax_per_nm"]*u"nm^-1", N_k=cfg["Nk"], N_φ=cfg["Nphi"], qz_max=cfg["qzmax_per_nm"]*u"nm^-1", N_qz=cfg["Nqz"], η_seed=cfg["seed_meV"]*u"meV")
algorithms = AlgorithmOptions(kernel_build=:direct, embedding=Symbol(cfg["embedding"]), embedding_periods=cfg["embedding_periods"], mixing=Symbol(get(cfg, "mixing", "linear")), energy_shift=Symbol(get(cfg, "energy_shift", "sparse_plan")))
first_request = request -> get(cfg,"save_first_state",false) && request.stage === :scba && request.iteration == 1
production = ProductionOptions(algorithms=algorithms, worker_count=cfg["workers"], hilbert_columns=get(cfg,"hilbert_columns",32), memory_budget_bytes=cfg["memory_budget_bytes"], phase_timing=true, checkpoint_every_scba=0, checkpoint_every_outer=0, checkpoint_request=first_request, progress_every_scba=1)
options = SolverOptions(α_Σ=cfg["alpha_sigma"], α_P=.15, max_scba=cfg["max_scba"], max_poisson=1, convergence=ConvergencePolicy(stagnation_window=0, stagnation_relative_improvement=0.0))
println("Start=", stamp(), " Julia=", VERSION, " Threads=", Threads.nthreads(), " BLAS=", BLAS.get_num_threads())
println("Core=", pathof(QCLNEGF), " Runner=", pathof(QCLNEGFRunner))
println("Input=", cfg)
flush(stdout)
identity = Dict{String,Any}("start_utc"=>stamp(), "julia_version"=>string(VERSION), "core_source"=>pathof(QCLNEGF), "runner_source"=>pathof(QCLNEGFRunner), "input_sha256"=>bytes2hex(sha256(read(input_path))), "script_sha256"=>bytes2hex(sha256(read(@__FILE__))), "core_head"=>readchomp(`git -C $(dirname(dirname(pathof(QCLNEGF)))) rev-parse HEAD`), "runner_head"=>readchomp(`git -C $(dirname(dirname(pathof(QCLNEGFRunner)))) rev-parse HEAD`))
open(joinpath(output, "identity.toml"), "w") do io
    TOML.print(io, identity)
end
built = build_configured_scattering_problem(physical=physical, numerical=numerical, scales=ScaleSystem(), scattering=ScatteringOptions(alloy=alloy), algorithms=algorithms, kernel_options=ProductionKernelOptions())
problem = built.problem
println("BuiltSeconds=", time()-started, " Problem=", problem)
flush(stdout)
history_path = joinpath(output, "history.csv")
fields = [name for name in fieldnames(SCBAIteration) if fieldtype(SCBAIteration, name) <: Real]
result = open(history_path, "w") do io
    println(io, join(["elapsed_s"; string.(fields)], ","))
    observer = function(row)
        println(io, join([string(time()-started); [string(getfield(row, field)) for field in fields]], ","))
        flush(io)
        println("Iteration=", row.ν, " rK=", row.r_K, " rSigma=", row.r_Σ, " rLambda=", row.r_λ, " PSD=", row.r_PSD, " elapsed_s=", time()-started)
        flush(stdout)
    end
    first_saved = Ref(false)
    callback = function(aligned)
        if get(cfg,"save_first_state",false) && !first_saved[]
            checkpoint_inner(joinpath(output,"first_state.h5"),problem,aligned)
            first_saved[] = true
            println("FirstStateSavedSeconds=",time()-started)
            flush(stdout)
        end
    end
    solve_scba_production(problem, zeros(numerical.N_z); options=options, production_options=production, history_observer=observer, iteration_callback=callback, initial=load_initial_guess(problem,production,cfg),consume_initial=true,resume_iterations=false)
end
checkpoint_inner(joinpath(output, "state.h5"), problem, result)
summary = Dict{String,Any}("finish_utc"=>stamp(), "elapsed_seconds"=>time()-started, "status"=>string(result.status), "converged"=>result.converged, "quality"=>string(result.quality), "iterations"=>length(result.history), "state_sha256"=>bytes2hex(open(sha256, joinpath(output,"state.h5"))))
for field in fields
    summary[string(field)] = getfield(last(result.history), field)
end
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, summary)
end
println("Finished=", summary)
flush(stdout)
