using QCLNEGF, QCLNEGFRunner, Unitful, LinearAlgebra, HDF5, Test
BLAS.set_num_threads(1)
length(ARGS) == 3 || error("usage: check_pair.jl CORE RUNNER NEW_OUTPUT.h5")
core, runner, output = abspath.(ARGS)
@test VERSION == v"1.13.3"
@test realpath(pkgdir(QCLNEGF)) == realpath(core)
@test realpath(pkgdir(QCLNEGFRunner)) == realpath(runner)
@test readchomp(`git -C $core rev-parse HEAD`) == "db26a0ef49673d8246ae336cbe73bafe6d7482ac"
@test readchomp(`git -C $runner rev-parse HEAD`) == "53b3a280109228a980df7ddea59e4f1b75b836be"
ispath(output) && error("output already exists")
println("Julia=", VERSION, " Core=", pathof(QCLNEGF), " Runner=", pathof(QCLNEGFRunner), " project=", Base.active_project())
physical = reference_parameters()
numerical = tutorial_numerics()
algorithms = AlgorithmOptions(kernel_build=:direct)
production = ProductionOptions(algorithms=algorithms, worker_count=2, memory_budget_bytes=2*1024^3, checkpoint_every_scba=0, checkpoint_every_outer=0)
built = build_configured_scattering_problem(physical=physical, numerical=numerical, scattering=ScatteringOptions(LO=false, acoustic=true, impurity=false, IFR=false, alloy=false), scales=ScaleSystem(), algorithms=algorithms, kernel_options=ProductionKernelOptions())
options = SolverOptions(max_scba=2, max_poisson=1, α_Σ=.2, α_P=.15)
result = solve_scba_production(built.problem, zeros(numerical.N_z); options=options, production_options=production)
@test !result.converged
@test length(result.history) == 2
@test all(isfinite(getfield(row, field)) for row in result.history for field in (:r_D,:r_A,:r_K,:r_Σ,:r_λ,:λ))
green = result.green
h5open(output, "w") do f
    write(f, "GR_real", real.(green.Gᴿ)); write(f, "GR_imag", imag.(green.Gᴿ))
    write(f, "Gl_real", real.(green.Gˡ)); write(f, "Gl_imag", imag.(green.Gˡ))
end
h5open(output, "r") do f
    @test read(f["GR_real"]) == real.(green.Gᴿ)
    @test read(f["GR_imag"]) == imag.(green.Gᴿ)
    @test read(f["Gl_real"]) == real.(green.Gˡ)
    @test read(f["Gl_imag"]) == imag.(green.Gˡ)
end
println("FixedUStatus=", result.status, " converged=", result.converged, " iterations=", length(result.history))
println("FinalPairChecksOK")
