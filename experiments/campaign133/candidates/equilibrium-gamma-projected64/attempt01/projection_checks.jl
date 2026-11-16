using Test, LinearAlgebra, TOML
function projection_small_checks(cfg)
    E=EquilibriumGammaPilot
    field(M)=reshape(ComplexF64.(M),1,1,size(M,1),size(M,2))
    matrix(X)=X[1,1,:,:]
    @testset "AA-only complex nearest-PSD projection and raw acceptance" begin
        H=ComplexF64[0.5 2.5im;-2.5im 0.5]
        expected=ComplexF64[1.5 1.5im;-1.5im 1.5] # independent analytic eigensystem {-2,3}
        X=field(H);before=copy(X);P,record,_=E.project_psd_trial(X,[1.0],[1.0])
        @test X==before
        @test matrix(P)≈expected rtol=1e-14 atol=1e-14
        @test record["absolute"]≈2.0 rtol=1e-14
        @test record["expected_orthogonal_distance"]≈2.0 rtol=1e-14
        @test record["changed_blocks"]==1 && record["negative_input_eigenvalues"]==1
        @test record["maximum_eigen_backward_relative"]<1e-13
        @test record["maximum_eigenvector_orthogonality_defect"]<1e-13
        twice,_,_=E.project_psd_trial(P,[1.0],[1.0])
        @test norm(twice-P)<1e-13
        # The analytic nearest point beats other independently specified PSD matrices.
        for candidate in (zeros(ComplexF64,2,2),0.5expected,expected+0.2I,ComplexF64[2 im;-im 2])
            @test eigmin(Hermitian(candidate))>=-1e-14
            @test norm(matrix(P)-H)<=norm(candidate-H)+1e-13
        end
        asymmetric=field(H);asymmetric[1,1,1,1]+=0.2im
        hermP,hermrecord,_=E.project_psd_trial(asymmetric,[0.2],[0.3])
        @test matrix(hermP)≈expected rtol=1e-14 atol=1e-14
        @test hermrecord["absolute"]≈sqrt(4.04) rtol=1e-14
        @test hermrecord["weighted_absolute"]≈sqrt(0.06)*hermrecord["absolute"] rtol=1e-14
        feasible=field(ComplexF64[2 0.1im;-0.1im 1]);same,same_record,_=E.project_psd_trial(feasible,[1.0],[1.0])
        @test same===feasible && reinterpret(UInt64,vec(same))==reinterpret(UInt64,vec(feasible))
        @test same_record["changed_blocks"]==0 && same_record["unchanged_field_identity"]===true
        raw=fill(-1.0+0im,1,1,1,1);zero,_,_=E.project_psd_trial(raw,[1.0],[1.0])
        @test only(zero)==0 # zero eigenvalue allowed; no positive floor
        _,rawres,_,_=E.gamma_residuals(raw,zero,[1.0],[1.0])
        @test rawres.relative==1 && rawres.relative>cfg["trigger_relative"] # false projected FP is not a root
        _,rootres,_,_=E.gamma_residuals(feasible,feasible,[1.0],[1.0])
        @test rootres.absolute==0 && rootres.relative==0
        @test raw==fill(-1.0+0im,1,1,1,1) # raw Phi not altered by projection
        parsed=TOML.parse(sprint(io->TOML.print(io,record)))
        @test parsed["minimum_input_location"]==[1,1]
        @test parsed["output_positive_rank_counts"]==[0,1,0]
    end
end
