# Definitions only. ROOT may run --tiny-mu separately; no preserved H5 inputs.
function tiny_mu_checks()
    @testset "native thermal LO affine FD source, complex matrices" begin
        cfg=Dict{String,Any}("temperature_K"=>70.0,"Nz"=>24,"Nb"=>2,"basis_periods"=>1,
            "NE"=>17,"energy_min_eV"=>-.3,"energy_max_eV"=>.9000001064962726,
            "trust_margin_meV"=>80.0,"kmax_per_nm"=>.6,"Nk"=>3,"Nphi"=>4,"Nqz"=>5,
            "qzmax_per_nm"=>10.0,"seed_halfwidth_meV"=>1.25,"workers"=>1,"memory_max_bytes"=>13019119616)
        c=GP.build_context(cfg); p=c.problem; energies=p.grids.ε
        A=zeros(ComplexF64,17,3,2,2)
        for e in 1:17,m in 1:3
            V=ComplexF64[1 .1im*e; .2+.1im*m 1+.03*e]
            A[e,m,:,:].=V*V'
        end
        f0=[.2+.1*sin(e) for e in 1:17]; f1=[.6+.1*cos(e) for e in 1:17]
        function native(f)
            GL=similar(A); GG=similar(A)
            for e in 1:17; GL[e,:,:,:].=im*f[e].*A[e,:,:,:]; GG[e,:,:,:].=-im*(1-f[e]).*A[e,:,:,:]; end
            green=GreenState(zeros(ComplexF64,size(A)),GL,GG,A,zeros(17,3),ones(17,3))
            L,G=N._production_lo_contraction(c.cache.kernels[:LO],green,p,c.cache,p.kernels.qᴷ[:LO],c.options)
            return im.*(G.-L)
        end
        # Independent actual-coordinate interpolation, not the native shift plan.
        function shifted(X,delta)
            Y=zeros(ComplexF64,size(X))
            for e in eachindex(energies)
                x=energies[e]+delta
                first(energies)<=x<=last(energies) || continue
                j=searchsortedlast(energies,x)
                if j==length(energies); Y[e,:,:,:].=X[j,:,:,:]
                else
                    t=(x-energies[j])/(energies[j+1]-energies[j])
                    Y[e,:,:,:].=(1-t).*X[j,:,:,:].+t.*X[j+1,:,:,:]
                end
            end
            return Y
        end
        X=similar(A)
        for e in 1:17; X[e,:,:,:].=(f1[e]-f0[e]).*A[e,:,:,:]; end
        quantum=N._scaled_physics(p.physical,p.scales).ħωᴸᴼ
        # Bose-independent exact delta: K[ shift+(Δf A) - shift-(Δf A) ].
        oracle=QCLNEGF._static_contraction(p.kernels.K[:LO],shifted(X,quantum)-shifted(X,-quantum),p.grids.wᵏ,p.kernels.qᴷ[:LO])
        B0=native(f0); B1=native(f1); mix=.37
        @test norm(B1-B0-oracle)<=3e-13*max(norm(B0),norm(B1))
        @test norm(native(mix*f1+(1-mix)*f0)-mix*B1-(1-mix)*B0)<=3e-13*max(norm(B0),norm(B1))
        @test norm(native(f0)-B0)==0
    end
    @testset "actual prolongation and streamed complex norm ledger" begin
        coarse_energy=[0.0,.4,2.0]; fine_energy=[0.0,.2,.4,1.2,2.0]
        X=reshape(ComplexF64[complex(i,sin(i)) for i in 1:24],3,2,2,2)
        P=zeros(ComplexF64,5,2,2,2)
        P[1:2:end,:,:,:].=X
        P[2,:,:,:].=(X[1,:,:,:]+X[2,:,:,:])/2
        P[4,:,:,:].=(X[2,:,:,:]+X[3,:,:,:])/2
        actual=mktempdir() do dir
            path=joinpath(dir,"tiny.h5")
            h5open(path,"w") do file; GP.write_complex(file,"Phi",X); end
            h5open(path,"r") do file; prolong_field(file,"Phi",coarse_energy,fine_energy); end
        end
        @test actual≈P atol=4e-15 rtol=4e-16
        D=(.03+.02im).*P; M=(-.01+.04im).*P
        we=[.1,.2,.3,.4,.5]; wk=[.2,.7]; s=empty_mu_moments()
        for b in 1:2,a in 1:2,m in 1:2,e in 1:5
            add_mu_moments!(s,1,D[e,m,a,b],M[e,m,a,b],P[e,m,a,b],2D[e,m,a,b],we[e]*wk[m])
        end
        # V1 dictionary reference preserves old update order; independent of slots.
        old=Dict{String,Float64}()
        for name in ("D","mu","fixed_mu","reference","fine_residual","fine_residual_fixed_mu"), suffix in ("_squared","_weighted_squared"); old[name*suffix]=0.0; end
        for name in ("D_mu","D_fixed","mu_fixed"),suffix in ("_inner","_weighted_inner"); old[name*suffix]=0.0; end
        for b in 1:2,a in 1:2,m in 1:2,e in 1:5
            d=D[e,m,a,b];mu=M[e,m,a,b];p=P[e,m,a,b];rr=2*d;w=we[e]*wk[m];f=d-mu;rf=rr-mu
            for (name,x) in (("D",d),("mu",mu),("fixed_mu",f),("reference",p),("fine_residual",rr),("fine_residual_fixed_mu",rf))
                q=abs2(x); old[name*"_squared"]+=q; old[name*"_weighted_squared"]+=w*q
            end
            for (name,x,y) in (("D_mu",d,mu),("D_fixed",d,f),("mu_fixed",mu,f))
                q=real(conj(x)*y);old[name*"_inner"]+=q;old[name*"_weighted_inner"]+=w*q
            end
        end
        r=finish_mu_moments(s)
        @test all(isequal(r[k],v) for (k,v) in old)
        @test r["D_absolute"]≈norm(D)
        @test r["fixed_mu_absolute"]≈norm(D-M)
        @test r["mu_fixed_inner"]≈real(dot(vec(M),vec(D-M)))
        independent=sum(we[e]*wk[m]*abs2(D[e,m,a,b]-M[e,m,a,b]) for b in 1:2,a in 1:2,m in 1:2,e in 1:5)
        @test r["fixed_mu_weighted_squared"]≈independent
        @test abs(r["squared_identity_error"])<=32eps()*r["reference_squared"]
        @test abs(r["weighted_squared_identity_error"])<=32eps()*r["reference_weighted_squared"]
        @test GP.relative(0.0,0.0)==0.0
        @test GP.relative(1.0,0.0)==Inf
        for bad in (:coarse,:fine,:gamma)
            mktempdir() do dir
                arrays=Dict(:coarse=>copy(X),:fine=>copy(P),:gamma=>copy(P))
                arrays[bad][1]=ComplexF64(NaN,0)
                for name in (:coarse,:fine,:gamma)
                    h5open(joinpath(dir,string(name)*".h5"),"w") do file
                        GP.write_complex(file,name===:gamma ? "Gamma" : "Phi",arrays[name])
                        name===:coarse && (file["energy"]=coarse_energy)
                    end
                end
                c=(problem=(grids=(ε=fine_energy,wᴱ=we,wᵏ=wk),scales=(E₀_eV=.1,)),)
                h5open(joinpath(dir,"fine.h5"),"r") do fine
                    h5open(joinpath(dir,"coarse.h5"),"r") do coarse
                        h5open(joinpath(dir,"gamma.h5"),"r") do gamma
                            @test_throws ErrorException stream_mu_decomposition(c,zero(P),fine,coarse,gamma,dir,Dict("anchor_indices"=>[1]))
                        end
                    end
                end
            end
        end
    end
end
