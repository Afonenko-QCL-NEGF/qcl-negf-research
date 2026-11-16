using Test, TOML
function restart_small_checks(cfg)
    E=EquilibriumGammaPilot
    @testset "stable-A GammaLast16 strict continuation metadata and budget" begin
        parent=copy(cfg);parent["identity"]="parent";parent["runtime_authorized"]=true
        parent["pi_calls_total"]=16;parent["project_AA_trial"]=false;parent["pilot_seconds"]=600.0;parent["overall_seconds"]=1800.0;parent["initialization"]="warm_source_GammaK_ghost0"
        identity=Dict("core_head"=>cfg["expected_core_head"],"core_tree"=>cfg["expected_core_tree"],
            "source_sha256"=>Dict("gamma_functions.jl"=>E.PARENT_GAMMA_OPERATOR_SHA,
                "stable_spectral.jl"=>E.PARENT_SPECTRAL_SHA))
        # Parent is a fresh-mode input; remove child-only source bindings.
        for key in collect(keys(parent));startswith(key,"parent_") && delete!(parent,key);end
        @test E.restart_contract(cfg,parent,identity)===nothing
        for (key,value) in (("NE",30721),("affine_l1_max",512.0),("history_differences",4),
                           ("spectral_operator","other"),("project_AA_trial",false),("pi_calls_total",65))
            changed=copy(cfg);changed[key]=value
            @test_throws ErrorException E.restart_contract(changed,parent,identity)
        end
        changed=deepcopy(identity);changed["source_sha256"]["stable_spectral.jl"]="wrong"
        @test_throws ErrorException E.restart_contract(cfg,parent,changed)
        changed=deepcopy(identity);changed["source_sha256"]["gamma_functions.jl"]="wrong"
        @test_throws ErrorException E.restart_contract(cfg,parent,changed)
        changed=deepcopy(identity);changed["core_head"]="wrong"
        @test_throws ErrorException E.restart_contract(cfg,parent,changed)
        @test E.restart_map_name(cfg)=="map-016.toml"
        current=Dict("residual"=>Dict("absolute"=>1.0,"scale"=>2.0,"relative"=>0.5,
            "weighted_absolute"=>0.3,"weighted_scale"=>0.6,"weighted_relative"=>0.5),
            "mu"=>Dict("mu"=>0.01,"number"=>0.4,"derivative"=>0.2))
        prior=deepcopy(current);metric=E.scalar_boundary_metrics(current,prior)
        @test metric["passed"] && length(metric["fields"])==9
        @test all(row["bitwise_equal"] for row in values(metric["fields"]))
        prior["residual"]["absolute"]=1.01
        @test !E.scalar_boundary_metrics(current,prior)["passed"]
        @test E.machine_same(0.0,0.0) && !E.machine_same(0.0,1e-100)
        @test E.machine_same(1.0,nextfloat(1.0))
        parsed=TOML.parse(sprint(io->TOML.print(io,metric)))
        @test parsed["fields"]["mu.mu"]["bitwise_equal"]===true
        budget=E.Budget(cfg);budget.stage=:pilot
        for _ in 1:64;E.actual_call!(budget;pi=true);end
        @test budget.pi==64
        @test_throws ErrorException E.actual_call!(budget;pi=true)
    end
end
