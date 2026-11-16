using Test, LinearAlgebra, TOML
function stable_small_checks(cfg)
    E=EquilibriumGammaPilot
    @testset "stable-A source and native-exit graph integration" begin
        toy=copy(cfg)
        for (key,value) in (("NE",65),("Nz",32),("Nb",2),("Nk",3),("Nphi",8),("Nqz",9));toy[key]=value;end
        c0=E.build_context(toy)
        rawh=copy(c0.h);rawh[1,1,1]+=1e-16im
        nativeplan=E.N._cavity_production_plan(c0.problem,rawh,c0.options)
        c=E.Context(c0.problem,c0.cache,c0.options,rawh,nativeplan,c0.kBT,c0.target)
        hbefore=copy(c.h);oldplan=c.cavity
        gamma=zeros(ComplexF64,65,3,2,2);block=ComplexF64[.02 .003im;-.003im .01]
        for e in 1:65,k in 1:3;gamma[e,k,:,:]=block;end
        GR,A,condition,scale,_,q,ledger=E.stable_spectral_state(c,gamma)
        @test isequal(c.h,hbefore) && c.cavity===oldplan
        @test ledger["h_correction"]["absolute"]>0
        @test ledger["q_cov"]["negative_count"]==0
        hH=c.cavity.hlocal;e=33;k=2;nb=2;central=9:10
        localR,_=E.local_retarded(c,gamma);D=zeros(ComplexF64,18,18);Kall=similar(D);fill!(Kall,0)
        T=c.problem.basis.T₊
        for cell in 1:9
            index=(2cell-1):(2cell)
            D[index,index]=c.problem.grids.ε[e]*I-hH[k,:,:]-localR[e,k,:,:]
            Kall[index,index]=block
            if cell<9
                next=(2cell+1):(2cell+2)
                D[index,next]=-T;D[next,index]=-T'
            end
        end
        full=inv(D);expectedGR=full[central,central]
        # Independent full-chain source sum, no production cavity recursion.
        expectedA=zeros(ComplexF64,2,2)
        for a in 1:2,b in 1:2,i in 1:18,j in 1:18
            expectedA[a,b]+=full[central[a],i]*Kall[i,j]*conj(full[central[b],j])
        end
        @test norm(GR[e,k,:,:]-expectedGR)/norm(expectedGR)<1e-9
        @test norm(A[e,k,:,:]-expectedA)/norm(expectedA)<1e-9
        target=0.3sum(q);fit=E.fit_mu(c.problem.grids.ε,q,target,c.kBT)
        green,_=E.fd_green(GR,A,condition,scale,c.problem.grids.ε,fit.mu,c.kBT)
        @test E.N._number_functional_production(green.Gˡ,c.problem.grids,c.problem.physical.g_s,c.options)≈target rtol=1e-11
        @test E.pair_A_error(green)/norm(A)<1e-14
        phi,channels,_=E.phi(c,green;keep_channels=true)
        native,_=E.N._scattering_candidate_production(c.problem,green,c.cache,c.options;return_roundoff=true)
        expected=zeros(ComplexF64,size(phi))
        for id in c.problem.kernels.enabled;expected .+= im.*(native[id].Σᵍ.-native[id].Σˡ);end
        @test E.norms(phi,expected,c.problem.grids.wᴱ,c.problem.grids.wᵏ).relative<1e-12
        @test isequal(c.h,hbefore) && c.cavity===oldplan # all consumers left native graph intact
        parsed=TOML.parse(sprint(io->TOML.print(io,ledger)))
        @test parsed["q_cov"]["negative_count"]==0
        @test parsed["native_exit_h_scope"]==ledger["native_exit_h_scope"]
    end
end
