using Test,LinearAlgebra
BLAS.set_num_threads(1)
include("matsubara_all5_refined.jl")
include("run_mu_root.jl")
const MR=MatsubaraAll5Refined
@testset "refined cutoff domain and bounded charge controls" begin
    @test BLAS.get_num_threads()==1
    @test BoundedMuRootRefined.sourcefree_controls().pass
    rows=MR.wired_onelevel(16.578,0.367)
    @test all(r.G_relative<1e-12 && r.reference_LO_relative<1e-12 for r in rows)
    @test rows[2].LO_absolute<rows[1].LO_absolute
    @test rows[2].density_absolute<rows[1].density_absolute
    raw=Dict(s=>fill(s===:LO ? 0.001 : 0.0,1,1,1,1,1,1) for s in (MR.STATIC...,:LO))
    c=MR.context(reshape(ComplexF64[0.1],1,1,1),zeros(ComplexF64,1,1),raw,[1.0],
        Dict(s=>1.0 for s in keys(raw)),16.578,0.05,0.367;spin=1)
    for N in (512,1024)
        snapshots=Any[]
        r=MR.pilot(c,N;steps=2,alpha=1.0,callback=s->push!(snapshots,deepcopy(s)))
        @test r.status=="iteration_budget" && length(r.history)==2
        @test size(r.S)==(N,1,1,1) && all(isfinite,r.Phi)
        @test all(iszero,snapshots[1].S) && snapshots[2].S==snapshots[1].Phi
        @test r.S==snapshots[2].S && r.G==snapshots[2].G && r.rho==snapshots[2].rho
    end
end
