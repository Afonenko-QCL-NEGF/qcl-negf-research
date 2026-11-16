using Test, LinearAlgebra, HDF5

function warm_small_checks(cfg)
    E=EquilibriumGammaPilot
    @testset "Warm Gamma import controls: no S01 map" begin
        old=[0.0,1.0,2.0];new=[-0.5,0.0,0.5,1.5,2.0,2.25]
        p=E.warm_interpolation_plan(old,new)
        @test all(w->0<=w<=1,p.weight_right)
        @test p.left[1]==0 && p.right[end]==4
        @test p.weight_right[end]==0.25
        @test_throws ErrorException E.warm_interpolation_plan(old,[-1.1,0.0])
        V=ComplexF64[1 im;0 2];Q=ComplexF64[.2 1+im;1-im -.3]
        A=V*adjoint(V);B=Q*adjoint(Q)
        X=zeros(ComplexF64,3,1,2,2);X[1,1,:,:]=A;X[2,1,:,:]=B;X[3,1,:,:]=A
        Y=E.interpolate_gamma(X,1,p,1:length(new),3)
        @test Y[3,1,:,:] ≈ (A+B)/2
        @test Y[1,1,:,:] ≈ A/2
        @test Y[end,1,:,:] ≈ 0.75A
        @test E.cone(Y).minimum_eigenvalue>=-1e-14
        constants=similar(X);constants[1,1,:,:]=A;constants[2,1,:,:]=A;constants[3,1,:,:]=A
        @test E.interpolate_gamma(constants,1,p,2:5,3)[2,1,:,:] ≈ A
        R=zeros(ComplexF64,1,1,2,2);L=similar(R);G=similar(R)
        gamma=ComplexF64[2 1+im;1-im 3];H=ComplexF64[1 0.3im;-0.3im 2]
        R[1,1,:,:]=H-im*gamma/2;L[1,1,:,:]=im*gamma/3;G[1,1,:,:]=-2im*gamma/3
        K,T=E.warm_broadenings(R,L,G)
        @test K[1,1,:,:] ≈ gamma
        @test T[1,1,:,:] ≈ gamma
        @test norm(T[1,1,:,:]-(-2 .*imag.(R[1,1,:,:])))>0.1
        cold=Dict{String,Any}("identity"=>"cold","runtime_authorized"=>true,"NE"=>30704,"temperature_K"=>70.0)
        candidate=copy(cold);candidate["identity"]="warm";candidate["runtime_authorized"]=false
        @test E.same_pilot_options(candidate,cold)===nothing
        candidate["NE"]=30721
        @test_throws ErrorException E.same_pilot_options(candidate,cold)
        mktempdir() do directory
            path=joinpath(directory,"orientation.h5")
            probe=zeros(ComplexF64,3,2,2,2)
            for e in 1:3,m in 1:2,a in 1:2,b in 1:2
                probe[e,m,a,b]=1000*e+100*m+10*a+b+im*(e-m+a-b)
            end
            HDF5.h5open(path,"w") do file
                datasets=E.create_complex_chunks(file,"probe",size(probe),2)
                try
                    E.write_complex_rows(datasets,probe[1:2,:,:,:],1:2)
                    E.write_complex_rows(datasets,probe[3:3,:,:,:],3:3)
                finally
                    E.close_complex(datasets)
                end
            end
            HDF5.h5open(path,"r") do file
                @test size(file["probe_real"])==(3,2,2,2)
                @test E.read_complex_rows(file,"probe",2:3)==probe[2:3,:,:,:]
            end
            # Independent C/HDF order is (b,a,k,E), not Julia (E,k,a,b).
            raw=joinpath(directory,"external.bin")
            read(`h5dump -d /probe_real -s 1,0,1,2 -c 1,1,1,1 -b NATIVE -o $raw $path`,String)
            @test open(io->read(io,Float64),raw)==real(probe[3,2,1,2])
            record=Dict{String,String}()
            E.require_equal("probe",probe,copy(probe),record)
            @test haskey(record,"probe")
            bad=copy(probe);bad[1]+=1
            @test_throws ErrorException E.require_equal("probe",probe,bad,record)
        end
    end
end
