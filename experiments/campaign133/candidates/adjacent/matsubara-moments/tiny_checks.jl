# ROOT-only source-free controls; author did not execute.
using Test
include("../matsubara-all5/matsubara_all5.jl")
include("moment_tail.jl")
function tiny_controls()
    M=MatsubaraMoments
    @test M.relative(0.0,0.0)==0 && M.relative(1.0,0.0)==Inf
    for N in (128,256)
        t=M.tails(N,10.0);@test 0<t.t4_integral_lower<=t.t4<=t.t4_integral_upper
    end
    o=M.scalar_oracle(10.0);@test o.moment_absolute_error<=1e-14
    @test all(r.M3_error<=r.bare_error+1e-13 for r in o.rows)
    ids=(:acoustic,:impurity,:IFR,:alloy,:LO)
    K=Dict(s=>fill(s==:LO ? .03 : .005,1,1,1,1,1,1) for s in ids)
    c=MatsubaraAll5.context(reshape(ComplexF64[.1],1,1,1),zeros(ComplexF64,1,1),K,[2.0],Dict(s=>1.0 for s in ids),10.0,0.0,.2;spin=1)
    r=M.correct(c,fill(.2+0im,1,1,1),128);t=r.tails.t4
    S0=.04+(2*c.NB+1)*.06;B0=.2*S0+.004+.06*((2*c.NB+1)*.1+.2)
    @test r.rho[1,1,1]≈(.2+t*B0)/(1+2*.2*t*.06)
    @test r.oneshot[1,1,1]≈.2+t*(B0-2*.2*.06*.2)
    @test r.linear_weighted_absolute≈sqrt(2)*r.linear_residual_absolute
    @test r.eta<1 && r.linear_residual_relative<1e-12
    return "moment controls complete; no actual source/model/inner Φ used"
end
