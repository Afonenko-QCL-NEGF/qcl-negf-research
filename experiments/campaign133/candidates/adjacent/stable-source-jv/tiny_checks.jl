# Definitions only. Caller must explicitly invoke tiny_checks(GP) after source gate.
using Test, LinearAlgebra, QCLNEGF
function tiny_checks(GP)
    J=StableSourceJv
    @testset "raw-Phi callback precedes rejecting cone and later work" begin
        raw=fill(-1.0+0im,1,1,1,1);before=copy(raw);events=Symbol[];saved=Ref{Any}()
        gp=(cone=x->begin push!(events,:guard);error("manufactured cone failure") end,)
        budget=J.Limits();budget.pi=1
        persist=(x,stats)->begin push!(events,:persist);saved[]=(copy(x),stats) end
        @test_throws ErrorException J.preserve_raw_before_guard(gp,raw,(mu=.25,),budget,persist)
        @test events==[:persist,:guard] && saved[][1]==before && raw==before
        @test saved[][2].mu==.25 && budget.pi_completed==1 && budget.pi==1
        # If persistence itself fails, the guard/later PI cannot consume the raw.
        empty!(events);budget=J.Limits();budget.pi=1
        failpersist=(x,stats)->begin push!(events,:persist);error("manufactured IO failure") end
        @test_throws ErrorException J.preserve_raw_before_guard(gp,raw,(mu=.25,),budget,failpersist)
        @test events==[:persist] && raw==before && budget.pi_completed==1
    end
    @testset "noncommuting complex finite9 covariance derivative" begin
        for n in (2,5)
            X=ComplexF64[sin(a+2b)+im*cos(2a-b) for a in 1:n,b in 1:n]
            Y=ComplexF64[cos(2a+b)+im*sin(a-2b) for a in 1:n,b in 1:n]
            h=(X+X')/12;Lambda=(Y+Y')/17
            gamma=Matrix{ComplexF64}(I,n,n)/5+X*X'/30
            v=(Y+Y')/40;dLambda=(X+X')/23
            R=Lambda-im*gamma/2;dR=dLambda-im*v/2;T=(X+im*Y)/25;E=.71
            w=J.Work(n);G,A,dA=J.block!(w,E,h,R,gamma,dR,v,T)
            D=zeros(ComplexF64,9n,9n)
            for cell in 1:9
                i=(cell-1)*n+1:cell*n;D[i,i]=E*I-h-R
                if cell<9
                    j=cell*n+1:(cell+1)*n;D[i,j]=-T;D[j,i]=-T'
                end
            end
            F=inv(D);deltaF=F*kron(Matrix{ComplexF64}(I,9,9),dR)*F
            K=kron(Matrix{ComplexF64}(I,9,9),gamma)
            deltaK=kron(Matrix{ComplexF64}(I,9,9),v)
            fullA=F*K*F';fulldA=deltaF*K*F'+F*deltaK*F'+F*K*deltaF'
            cc=4n+1:5n
            @test norm(G-F[cc,cc])<2e-12*norm(F[cc,cc])
            @test norm(A-fullA[cc,cc])<3e-12*norm(fullA[cc,cc])
            @test norm(dA-fulldA[cc,cc])<5e-12*norm(fulldA[cc,cc])
            # Independent full-chain central difference, no cavity/tangent reuse.
            function full_at(t)
                F=inv(D-t*kron(Matrix{ComplexF64}(I,9,9),dR))
                return (F*(K+t*deltaK)*F')[cc,cc]
            end
            fd=(full_at(1e-5)-full_at(-1e-5))/(2e-5)
            @test norm(fd-dA)<2e-8*norm(dA)
        end
    end
    cfg=Dict{String,Any}("temperature_K"=>70.0,"Nz"=>24,"Nb"=>2,"basis_periods"=>1,
        "NE"=>17,"energy_min_eV"=>-.3,"energy_max_eV"=>.9000001064962726,
        "trust_margin_meV"=>80.0,"kmax_per_nm"=>.6,"Nk"=>3,"Nphi"=>4,"Nqz"=>5,
        "qzmax_per_nm"=>10.0,"seed_halfwidth_meV"=>1.25,"workers"=>1,"memory_max_bytes"=>13019119616)
    c=GP.build_context(cfg);shape=(17,3,2,2);A=zeros(ComplexF64,shape);dA=similar(A)
    for e in 1:17,m in 1:3
        A[e,m,:,:]=ComplexF64[1+.01e .1im;-.1im 2+.02m]
        dA[e,m,:,:]=ComplexF64[sin(e)/10 .03+.02im;.03-.02im cos(m)/7]
    end
    base=(A=A,stats=(mu=(mu=.12,),))
    @testset "implicit-current-mu and actual five native linear channels" begin
        GL,GG,number=J.number_tangent(GP,c,base,dA)
        @test number.chi>0 && number.dN_relative<3e-14
        @test norm(im*(GG-GL)-dA)<4e-16*norm(dA)
        one=GreenState(dA,GL,GG,dA,zeros(17,3),ones(17,3))
        GL2=reverse(GL;dims=1);GG2=reverse(GG;dims=1);A2=reverse(dA;dims=1)
        two=GreenState(A2,GL2,GG2,A2,zeros(17,3),ones(17,3))
        a,b=.7,-.3
        mix=GreenState(A2,a*GL+b*GL2,a*GG+b*GG2,a*dA+b*A2,zeros(17,3),ones(17,3))
        p1,k1,_=GP.phi(c,one;keep_channels=true);p2,k2,_=GP.phi(c,two;keep_channels=true)
        pm,km,_=GP.phi(c,mix;keep_channels=true)
        @test Set(keys(km))==Set((GP.STATIC_IDS...,:LO))
        for id in keys(km)
            scale=norm(a*k1[id])+norm(b*k2[id]);err=norm(km[id]-a*k1[id]-b*k2[id])
            @test err<=2e-12*scale
        end
        @test norm(pm-a*p1-b*p2)<=2e-12*(norm(a*p1)+norm(b*p2))
    end
    return nothing
end
