# This file defines checks. Inclusion does not execute them; run_gamma invokes
# these only AFTER the explicit runtime/source gate.
using Test, LinearAlgebra

function small_checks(cfg)
    E = EquilibriumGammaPilot
    @testset "Gamma pilot: independent small controls" begin
        for x in (-1000.0,-100.0,-1.0,0.0,1.0,100.0,1000.0)
            f,h = E.fermi_pair(x)
            @test 0 <= f <= 1 && 0 <= h <= 1
            @test abs(f+h-1) <= eps(Float64)
        end
        @test E.fermi_pair(-100.0)[2] > 0 # subtraction 1-f would lose it
        @test E.relative(0.0,0.0) == 0
        @test E.relative(1.0,0.0) == Inf
        fit = E.fit_mu([.25],[2.0],.3,.1)
        @test isapprox(fit.mu,.25-.1*log((2-.3)/.3);atol=3e-13,rtol=0)
        @test fit.derivative > 0
        @test_throws E.Infeasible E.fit_mu([0.0],[0.0],.3,.1)
        @test_throws E.Infeasible E.fit_mu([0.0],[.2],.3,.1)
        @test_throws E.Infeasible E.fit_mu([0.0],[-.2],.3,.1)
        @test_throws E.Infeasible E.cone(zeros(ComplexF64,1,1,1,1))
        @test_throws E.Infeasible E.cone(fill(-1.0+0im,1,1,1,1))
        @test_throws E.Infeasible E.cone(fill(1.0+im,1,1,1,1))
        X = reshape(ComplexF64[1 2+3im;2-3im 4],1,1,2,2)
        Y = reshape(ComplexF64[2 -1+.4im;-1-.4im .5],1,1,2,2)
        @test E.inner(X,Y,[.3],[.7]) ≈ .3*.7*real(tr(adjoint(X[1,1,:,:])*Y[1,1,:,:])) atol=1e-15

        # CP complex-matrix oracle: independent explicit V X Vdagger loops
        # against actual literal Core contraction; non-Hermitian rational GR.
        V = ComplexF64[1 im;-im 2]
        W = ComplexF64[.2 1+im;1-im -.3]
        z = .7+.31im
        G = ComplexF64[z-2 1+im;1-im z]/(z*z-2z-2)
        K = zeros(ComplexF64,1,1,2,2,2,2)
        expected = zeros(ComplexF64,2,2)
        for (weight,vertex) in ((.3,V),(.7,W)),a in 1:2,b in 1:2,c in 1:2,d in 1:2
            K[1,1,a,c,d,b] += weight*vertex[a,c]*conj(vertex[b,d])
            expected[a,b] += weight*vertex[a,c]*G[c,d]*conj(vertex[b,d])
        end
        field = zeros(ComplexF64,1,1,2,2); field[1,1,:,:] = G
        measured = E.N._static_contraction(K,field,[1.0])
        @test measured[1,1,:,:] ≈ expected rtol=2e-14 atol=1e-15
        production = E.N.production_static_contraction(K,field,[1.0];energy_chunk=1)
        @test production[1,1,:,:] ≈ expected rtol=2e-14 atol=1e-15
        A = im*(G-adjoint(G)); field[1,1,:,:] = A
        gamma = E.N._static_contraction(K,field,[1.0])[1,1,:,:]
        @test im*(expected-adjoint(expected)) ≈ gamma rtol=2e-14 atol=1e-15

        # Weighted QR/type-II step compared with independent dense real LS.
        field1 = reshape(ComplexF64[.1,.3,.4,.8],2,2,1,1)
        field2 = reshape(ComplexF64[.2,.5,.7,1.1],2,2,1,1)
        field3 = reshape(ComplexF64[.4,.6,.9,1.4],2,2,1,1)
        f1 = reshape(ComplexF64[.3,.7,.1,.5],2,2,1,1)
        f2 = reshape(ComplexF64[.4,.2,.3,.4],2,2,1,1)
        f3 = reshape(ComplexF64[.2,.4,.6,.1],2,2,1,1)
        xs,fs = E.Field[field1,field2,field3],E.Field[f1,f2,f3]
        we,wk = [.3,.7],[.2,.8]
        trial,info = E.aa_trial(xs,fs,we,wk,cfg)
        weights = sqrt.(vec(we*transpose(wk)))
        DF = hcat(real.(vec(f2-f1)),real.(vec(f3-f2)))
        DX = hcat(vec(field2-field1),vec(field3-field2))
        coeff = (weights.*DF)\(weights.*real.(vec(f3)))
        oracle = vec(field3)+cfg["aa_beta"]*vec(f3)-(DX+cfg["aa_beta"]*DF)*coeff
        @test info.kind == "AAII"
        @test eltype(info.coefficients) <: Real
        @test vec(trial) ≈ oracle rtol=3e-13 atol=3e-14

        # Nonlinear positive scalar fixed point x=.5+.2*x^2, starting x=.1.
        xs,fs = E.Field[],E.Field[]
        x = fill(.1+0im,1,1,1,1)
        for iteration in 1:20
            mapped = @. .5+.2*x*x
            f = mapped-x
            push!(xs,x); push!(fs,f)
            length(xs)>4 && (popfirst!(xs);popfirst!(fs))
            trial,_ = E.aa_trial(xs,fs,[1.0],[1.0],cfg)
            x = try
                E.cone(trial); trial
            catch error
                error isa E.Infeasible || rethrow()
                x+.2*f
            end
        end
        @test real(only(x)) ≈ (1-sqrt(.6))/.4 atol=2e-11
        @test abs(only(.5 .+ .2.*x.*x .- x)) < 1e-11
        @test E.cone(x).minimum_eigenvalue > 0

        # Check the REAL artifact consumer path before expensive native work.
        mktempdir() do directory
            for shape in ((3,2),(3,2,2),(3,2,2,2),(3,1,1,2,2,2))
                source = reshape(ComplexF64.(1:prod(shape)).*(1+.2im),shape)
                path = joinpath(directory,"writer$(length(shape)).h5")
                E.atomic_h5(path,cfg;expected_bytes=1024^2) do file
                    E.write_complex(file,"probe",source)
                end
                E.HDF5.h5open(path,"r") do file
                    @test read(file["probe_real"])+im*read(file["probe_imag"]) == source
                end
            end
        end
        budget = E.Budget(cfg)
        budget.stage = :pilot
        for _ in 1:16; E.actual_call!(budget;pi=true); end
        @test budget.pi == 16
        @test_throws ErrorException E.actual_call!(budget;pi=true)
    end
    return nothing
end
