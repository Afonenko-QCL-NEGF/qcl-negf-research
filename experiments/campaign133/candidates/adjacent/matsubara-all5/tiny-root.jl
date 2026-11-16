using Test, LinearAlgebra
include("matsubara_all5.jl")
const M=MatsubaraAll5
@testset "Mats ALL5 actual primitive wiring and aligned snapshots" begin
    beta=inv(0.060321332835);Omega=0.367
    rows=M.wired_onelevel(beta,Omega)
    @test length(rows)==2
    for r in rows
        @test r.G_relative<1e-12
        @test r.reference_LO_relative<1e-12
        @test r.dyson_backward<1e-12
        @test isfinite(r.LO_absolute) && isfinite(r.density_absolute)
    end
    @test rows[2].LO_absolute<=rows[1].LO_absolute+1e-12
    @test rows[2].density_absolute<=rows[1].density_absolute+1e-12
    pole=M.onelevel_oracle(beta,Omega;N=256,g2=0.001,epsilon=0.02)
    @test length(pole)==6
    for i in (1,3,5)
        @test pole[i+1].absolute < pole[i].absolute
    end
    raw=Dict(s=>fill(s===:LO ? 0.001 : 0.0,1,1,1,1,1,1) for s in (M.STATIC...,:LO))
    c=M.context(reshape(ComplexF64[0.1],1,1,1),zeros(ComplexF64,1,1),raw,[1.0],
        Dict(s=>1.0 for s in keys(raw)),beta,0.05,Omega;spin=1)
    @test maximum(values(M.qualify_native_ops(c)))<2e-12
    snapshots=Any[]
    state=M.pilot(c,128;steps=2,alpha=1.0,callback=s->push!(snapshots,deepcopy(s)))
    @test length(snapshots)==2 && state.status=="iteration_budget"
    @test all(iszero,snapshots[1].S)
    @test snapshots[2].S==snapshots[1].Phi
    @test state.S==snapshots[2].S && state.G==snapshots[2].G && state.Phi==snapshots[2].Phi
    @test state.rho==snapshots[2].rho && state.w==snapshots[2].w
    @test M.relative(0.0,0.0)==0.0 && isinf(M.relative(1.0,0.0))
    println("WIRED_ONELEVEL_RECORDS ", rows)
    println("ONELEVEL_RAW_SUM ", pole)
end
