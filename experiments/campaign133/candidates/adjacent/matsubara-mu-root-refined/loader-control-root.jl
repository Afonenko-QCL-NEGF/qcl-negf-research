using Test,TOML,LinearAlgebra
BLAS.set_num_threads(1)
include("run_mu_root.jl")
function exercise_loader(c,sconf)
    Base.include(Main,BoundedMuRootRefined.path(c,"source_adapter"))
    Base.invokelatest(() -> Main.load_modules(sconf))
    Base.include(Main,BoundedMuRootRefined.path(c,"refined_module"))
    Base.include(Main,BoundedMuRootRefined.path(c,"M3_module"))
end
c=TOML.parsefile(joinpath(@__DIR__,"input.toml"))
s=TOML.parsefile(BoundedMuRootRefined.path(c,"source_input"))
@testset "latest-world source loading without model or HDF arrays" begin
    @test c["runtime_authorized"]===false
    exercise_loader(c,s)
    @test all(isdefined(Main,n) for n in (:EquilibriumGammaPilot,:MatsubaraAll5,
        :GammaDensityHermitianFreezeIO,:MatsubaraAll5Refined,:MatsubaraMoments))
    @test isdefined(Main,:source) && isdefined(Main,:save_state)
    @test BLAS.get_num_threads()==1 && Threads.nthreads()==1
end
