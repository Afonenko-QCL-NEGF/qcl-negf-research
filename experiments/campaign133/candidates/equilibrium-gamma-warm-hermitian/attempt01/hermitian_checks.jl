function hermitian_small_checks()
    E=EquilibriumGammaPilot
    @testset "PI symmetry restoration and tiny broadening" begin
        gamma=reshape(ComplexF64[2e-27 0.3e-27im; -0.3e-27im 1e-27],1,1,2,2)
        H=reshape(ComplexF64[0.2 0.1+0.3im; 0.1-0.3im -0.4],1,1,2,2)
        raw=copy(H);raw[1,1,1,1]+=1e-18im
        restored=E.hermitian_copy(raw)
        @test restored==H
        R=restored.-0.5im.*gamma
        @test imag(R[1,1,1,1]) == -real(gamma[1,1,1,1])/2
        @test imag(R[1,1,2,2]) == -real(gamma[1,1,2,2])/2
        @test gamma[1,1,1,1]==2e-27
        @test E.norms(raw,restored,[1.0],[1.0]).absolute ≈ 1e-18 rtol=1e-14
        # Independent real-weight dense linear map commutes with the symmetry.
        weights=[0.1 -0.4; 0.7 0.3]
        fields=cat(H,2H;dims=1)
        exact=reshape(weights*reshape(fields,2,:),2,1,2,2)
        @test E.hermitian_copy(exact) ≈ exact rtol=2e-15
    end
end
