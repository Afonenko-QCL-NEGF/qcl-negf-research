using Test, HDF5, LinearAlgebra
function restart_small_checks(cfg)
    E=EquilibriumGammaPilot
    @testset "Gamma-only restart source and HDF admission" begin
        parent=Dict{String,Any}("identity"=>"parent","runtime_authorized"=>true,"initialization"=>"warm_source_GammaK_ghost0",
            "retarded_real_symmetry"=>"hermitian_projection","pi_calls_total"=>16,"NE"=>30704)
        identity=Dict("core_head"=>cfg["expected_core_head"],"core_tree"=>cfg["expected_core_tree"],
            "source_sha256"=>Dict("gamma_functions.jl"=>E.PARENT_GAMMA_OPERATOR_SHA))
        candidate=copy(cfg)
        @test E.restart_contract(candidate,parent,identity)===nothing
        changed=copy(candidate);changed["NE"]=30721
        @test_throws ErrorException E.restart_contract(changed,parent,identity)
        changed=deepcopy(identity);changed["source_sha256"]["gamma_functions.jl"]="wrong"
        @test_throws ErrorException E.restart_contract(candidate,parent,changed)
        @test E.machine_same(0.0,0.0)
        @test !E.machine_same(0.0,1e-100)
        @test E.machine_same(1.0,nextfloat(1.0))
        @test !E.machine_same(1.0,1.001)
        budget=E.Budget(cfg);budget.stage=:pilot
        for _ in 1:128;E.actual_call!(budget;pi=true);end
        @test budget.pi==128
        @test_throws ErrorException E.actual_call!(budget;pi=true)
        mktempdir() do directory
            # Existing policy yields row chunk3 and pads logical4 rows to6.
            # Only16MiB logical/24MiB allocated; old sizeof admission is unsafe.
            X=zeros(ComplexF64,4,262145);X[end,end]=1+2im
            path=joinpath(directory,"padded.h5")
            declared=E.expected_padded_complex_bytes(size(X))+1024^2
            HDF5.h5open(path,"w") do file;E.write_complex(file,"probe",X);end
            @test filesize(path)>sizeof(X)
            @test declared>=filesize(path)
            HDF5.h5open(path,"r") do file
                @test file["probe_real"][4,262145]==1.0
                @test file["probe_imag"][4,262145]==2.0
            end
            checkpoint=joinpath(directory,"Gamma.h5")
            gamma=reshape(ComplexF64[2 0.2im;-0.2im 1],1,1,2,2)
            E.gamma_checkpoint(directory,"Gamma",gamma,cfg,16)
            HDF5.h5open(checkpoint,"r") do file
                @test read(file["schema"])=="campaign133-gamma-checkpoint-v1"
                @test E.read_complex(file,"Gamma")==gamma
            end
        end
    end
end
