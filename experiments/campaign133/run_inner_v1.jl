using QCLNEGF, QCLNEGFRunner, Unitful, LinearAlgebra, TOML, SHA, HDF5, Dates

length(ARGS) == 2 || error("usage: run_inner.jl INPUT.toml NEW_OUTPUT_DIRECTORY")
input_path, output = abspath.(ARGS)
ispath(output) && error("Refusing to overwrite existing run: $output")
mkpath(output)
cp(input_path, joinpath(output, "input.toml"))
cfg = TOML.parsefile(input_path)
BLAS.set_num_threads(1)
started = time()
stamp() = string(Dates.now(Dates.UTC)) * "Z"
function checkpoint_inner(path, problem, result)
    HDF5.h5open(path * ".tmp", "w") do file
        file["schema"] = "campaign133-inner-raw-v1"
        file["status"] = string(result.status)
        file["converged"] = result.converged
        file["quality"] = string(result.quality)
        for (name, array) in (("energy", problem.grids.ε), ("k", problem.grids.κ), ("x", problem.grids.x), ("weight_energy", problem.grids.wᴱ), ("weight_k", problem.grids.wᵏ), ("weight_x", problem.grids.wˣ))
            file[name] = array
        end
        file["hartree"] = zeros(problem.numerical.N_z)
        function complex_array(group, name, array)
            group[name * "_real"] = real.(array)
            group[name * "_imag"] = imag.(array)
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

physical = reference_parameters(V_period=cfg["voltage_mV"]*u"mV", Tᴸ=cfg["temperature_K"]*u"K", Tᴸᴼ=cfg["temperature_K"]*u"K")
numerical = NumericalParameters(N_z=cfg["Nz"], N_b=cfg["Nb"], P_basis=cfg["basis_periods"], E_min=cfg["energy_min_eV"]*u"eV", E_max=cfg["energy_max_eV"]*u"eV", N_E=cfg["NE"], M_E=cfg["trust_margin_meV"]*u"meV", k_max=cfg["kmax_per_nm"]*u"nm^-1", N_k=cfg["Nk"], N_φ=cfg["Nphi"], qz_max=cfg["qzmax_per_nm"]*u"nm^-1", N_qz=cfg["Nqz"], η_seed=cfg["seed_meV"]*u"meV")
algorithms = AlgorithmOptions(kernel_build=:direct, embedding=Symbol(cfg["embedding"]), embedding_periods=cfg["embedding_periods"], mixing=Symbol(get(cfg, "mixing", "linear")), energy_shift=Symbol(get(cfg, "energy_shift", "sparse_plan")))
production = ProductionOptions(algorithms=algorithms, worker_count=cfg["workers"], memory_budget_bytes=cfg["memory_budget_bytes"], phase_timing=true, checkpoint_every_scba=0, checkpoint_every_outer=0, progress_every_scba=1)
options = SolverOptions(α_Σ=cfg["alpha_sigma"], α_P=.15, max_scba=cfg["max_scba"], max_poisson=1, convergence=ConvergencePolicy(stagnation_window=0, stagnation_relative_improvement=0.0))
println("Start=", stamp(), " Julia=", VERSION, " Threads=", Threads.nthreads(), " BLAS=", BLAS.get_num_threads())
println("Core=", pathof(QCLNEGF), " Runner=", pathof(QCLNEGFRunner))
println("Input=", cfg)
flush(stdout)
identity = Dict{String,Any}("start_utc"=>stamp(), "julia_version"=>string(VERSION), "core_source"=>pathof(QCLNEGF), "runner_source"=>pathof(QCLNEGFRunner), "input_sha256"=>bytes2hex(sha256(read(input_path))), "script_sha256"=>bytes2hex(sha256(read(@__FILE__))), "core_head"=>readchomp(`git -C $(dirname(dirname(pathof(QCLNEGF)))) rev-parse HEAD`), "runner_head"=>readchomp(`git -C $(dirname(dirname(pathof(QCLNEGFRunner)))) rev-parse HEAD`))
open(joinpath(output, "identity.toml"), "w") do io
    TOML.print(io, identity)
end
built = build_configured_scattering_problem(physical=physical, numerical=numerical, scales=ScaleSystem(), scattering=default_scattering(), algorithms=algorithms, kernel_options=ProductionKernelOptions())
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
    solve_scba_production(problem, zeros(numerical.N_z); options=options, production_options=production, history_observer=observer)
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
