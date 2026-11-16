using QCLNEGF, QCLNEGFRunner, Unitful, LinearAlgebra, Test
BLAS.set_num_threads(1)
root = normpath(joinpath(@__DIR__, "..", ".."))
core = joinpath(root, "components", "QCLNEGF.jl")
println("Julia=", VERSION, " Core=", pathof(QCLNEGF), " Runner=", pathof(QCLNEGFRunner))
for name in ("scalar_dyson_and_spectral_identity", "reference2019_units_geometry_and_quadratures", "manufactured_periodic_poisson_solution", "threaded_production_seed_equals_the_educational_discretization", "production_lo_equilibrium_discretization")
    include(joinpath(core, "test", "physics", name * ".jl"))
end
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
checkpoint = joinpath(@__DIR__, "raw", "compat-smoke.h5")
QCLNEGFRunner.save_checkpoint(checkpoint, solution)
restored = QCLNEGFRunner.load_checkpoint(checkpoint)
@test restored.Gᴿ == solution.scba.green.Gᴿ
@test restored.Gˡ == solution.scba.green.Gˡ
@test restored.Uᴴ == solution.Uᴴ
println("SmokeStatus=", solution.status, " SCBA=", solution.scba.status, " report=", solution.report.metrics)
println("DirectChecksOK")
