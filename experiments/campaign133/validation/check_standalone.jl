using QCLNEGF, QCLNEGFRunner, Unitful, LinearAlgebra, Test, Pkg
BLAS.set_num_threads(1)
println("Julia=", VERSION, " Core=", pathof(QCLNEGF), " Runner=", pathof(QCLNEGFRunner))
compat = Pkg.Types.semver_spec("~1.13.0")
@test v"1.13.0" in compat
@test VERSION in compat
@test !(v"1.12.9" in compat)
@test !(v"1.14.0" in compat)
physical = reference_parameters()
numerical = tutorial_numerics()
algorithms = AlgorithmOptions(kernel_build=:direct)
production = ProductionOptions(algorithms=algorithms, worker_count=2, memory_budget_bytes=2*1024^3, checkpoint_every_scba=0, checkpoint_every_outer=0)
built = build_configured_scattering_problem(physical=physical, numerical=numerical, scattering=ScatteringOptions(LO=false, acoustic=true, impurity=false, IFR=false, alloy=false), scales=ScaleSystem(), algorithms=algorithms, kernel_options=ProductionKernelOptions())
options = SolverOptions(max_scba=2, max_poisson=1, α_Σ=.2, α_P=.15)
solution = QCLNEGFRunner.solve_production(built.problem; options=options, production_options=production)
@test solution.converged == false
@test !isempty(solution.scba.history)
@test all(isfinite(getfield(row, field)) for row in solution.scba.history for field in (:r_D,:r_A,:r_K,:r_Σ,:r_λ,:λ))
checkpoint = joinpath(@__DIR__, "raw", "standalone-smoke.h5")
QCLNEGFRunner.save_checkpoint(checkpoint, solution)
restored = QCLNEGFRunner.load_checkpoint(checkpoint)
@test restored.Gᴿ == solution.scba.green.Gᴿ
@test restored.Gˡ == solution.scba.green.Gˡ
@test restored.Uᴴ == solution.Uᴴ
println("SmokeStatus=", solution.status, " SCBA=", solution.scba.status, " report=", solution.report.metrics)
println("DirectChecksOK")
