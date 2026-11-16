using Test, LinearAlgebra, TOML
include("aa_reference.jl")
function affine_small_checks(cfg)
    E=EquilibriumGammaPilot
    scalar(x)=fill(ComplexF64(x),1,1,1,1)
    @testset "AA telemetry and bounded affine4096 variant" begin
        c10=copy(cfg);c10["affine_l1_max"]=10.0
        phi(x)=1.0+0.9*(x-1.0)
        x0=2.0;f0=phi(x0)-x0;x1=x0+0.2*f0;f1=phi(x1)-x1
        xs,fs=E.Field[scalar(x0),scalar(x1)],E.Field[scalar(f0),scalar(f1)]
        rejected,info10=E.aa_trial(xs,fs,[1.0],[1.0],c10;cached_norms=[abs(f0),abs(f1)])
        c512=copy(cfg);c512["affine_l1_max"]=512.0
        admitted,info512=E.aa_trial(xs,fs,[1.0],[1.0],c512;cached_norms=[abs(f0),abs(f1)])
        bold,info4096=E.aa_trial(xs,fs,[1.0],[1.0],cfg;cached_norms=[abs(f0),abs(f1)])
        @test info10.kind=="linear_rank_reset"
        @test info10.suffix_decisions[end]["reason"]=="affine_l1_limit"
        @test info512.kind=="AAII" && info4096.kind=="AAII"
        @test real(only(rejected))>1.9
        @test real(only(admitted))≈1.0 atol=1e-12
        @test reinterpret(UInt64,vec(admitted))==reinterpret(UInt64,vec(bold))
        @test info512.l1≈99.0 atol=1e-10
        @test info10.suffix_decisions[end]["m1_a_cached_working_norm"]==abs(f1)
        @test info10.suffix_decisions[end]["m1_b_cached_working_norm"]==abs(f0)
        # Include accepted multidimensional and rank-deficient suffix cases.
        x2=x1+0.2*f1;f2=phi(x2)-x2;x3=x2+0.2*f2;f3=phi(x3)-x3
        vector(v)=reshape(ComplexF64.(v),2,2,1,1)
        vx=E.Field[vector([.1,.3,.4,.8]),vector([.2,.5,.7,1.1]),vector([.4,.6,.9,1.4])]
        vf=E.Field[vector([.3,.7,.1,.5]),vector([.4,.2,.3,.4]),vector([.2,.4,.6,.1])]
        fixtures=((xs[1:1],fs[1:1],[1.0],[1.0]),(xs,fs,[1.0],[1.0]),
            (E.Field[scalar(x0),scalar(x1),scalar(x2),scalar(x3)],E.Field[scalar(f0),scalar(f1),scalar(f2),scalar(f3)],[1.0],[1.0]),
            (vx,vf,[.3,.7],[.2,.8]))
        for (states,residuals,we,wk) in fixtures
            old,oldinfo=ParentAAOracle.parent_aa_trial(states,residuals,we,wk,c10)
            new,newinfo=E.aa_trial(states,residuals,we,wk,c10)
            @test reinterpret(UInt64,vec(old))==reinterpret(UInt64,vec(new))
            @test all(isequal(getproperty(oldinfo,k),getproperty(newinfo,k)) for k in propertynames(oldinfo))
        end
        row=Dict{String,Any}("AA"=>Dict(string(k)=>v for (k,v) in pairs(info10)))
        roundtrip=TOML.parse(sprint(io->TOML.print(io,row)))
        @test roundtrip["AA"]["suffix_decisions"][end]["reason"]=="affine_l1_limit"
        @test roundtrip["AA"]["suffix_decisions"][end]["actual_affine_l1"]==info10.suffix_decisions[end]["actual_affine_l1"]
        parent=copy(cfg);parent["pi_calls_total"]=128;parent["affine_l1_max"]=10.0
        parent["pilot_seconds"]=1200.0;parent["overall_seconds"]=2400.0
        identity=Dict("core_head"=>cfg["expected_core_head"],"core_tree"=>cfg["expected_core_tree"],
            "source_sha256"=>Dict("gamma_functions.jl"=>E.PARENT_GAMMA_OPERATOR_SHA))
        @test E.restart_contract(cfg,parent,identity)===nothing
        wrong=copy(cfg);wrong["qr_condition_max"]=1e9
        @test_throws ErrorException E.restart_contract(wrong,parent,identity)
        wrong=copy(cfg);wrong["affine_l1_max"]=512.0
        @test_throws ErrorException E.restart_contract(wrong,parent,identity)
        @test E.restart_map_name(cfg)=="map-128.toml"
        budget=E.Budget(cfg);budget.stage=:pilot
        for _ in 1:32;E.actual_call!(budget;pi=true);end
        @test_throws ErrorException E.actual_call!(budget;pi=true)
    end
end
