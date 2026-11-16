# WRITE-ONLY until root explicitly runs this; uses Julia stdlibs only.
using Test, LinearAlgebra
include("CauchyPanelTree.jl")
using .CauchyPanelTree
function stackmat(matrices)
    n=length(matrices); b=size(first(matrices),1)
    out=zeros(ComplexF64,n,b,b)
    for j in 1:n; out[j,:,:] .= matrices[j]; end
    return out
end
# Independent positive-hat divided difference of F(t)=t*log(t).
# No panel_integral, tree, quadrature, or translated moment helper is reused.
function hat_reference(E,G,z,ghosts; real_pv=false)
    setprecision(BigFloat,256) do
        x=BigFloat[ghosts[1]; E; ghosts[2]]
        zz=real_pv ? BigFloat(z) : Complex{BigFloat}(z)
        F(t)=iszero(t) ? zero(t) : t*(real_pv ? log(abs(t)) : log(t))
        out=zeros(Complex{BigFloat},size(G,2),size(G,3))
        for j in eachindex(E)
            a,c,b=x[j],x[j+1],x[j+2]
            w=(F(zz-a)-F(zz-c))/(c-a)-(F(zz-c)-F(zz-b))/(b-c)
            out .+= (w/(2big(pi))).*Complex{BigFloat}.(G[j,:,:])
        end
        if real_pv
            # Independent scalar interpolation for the single boundary imaginary term.
            if first(x)<=z<=last(x)
                j=clamp(searchsortedlast(x,zz),1,length(x)-1)
                t=(zz-x[j])/(x[j+1]-x[j])
                gl=j==1 ? zeros(Complex{BigFloat},size(G,2),size(G,3)) : G[j-1,:,:]
                gr=j+1==length(x) ? zero(gl) : G[j,:,:]
                out .-= (big(0.5)*im).*((1-t).*gl.+t.*gr)
            end
        end
        return ComplexF64.(out)
    end
end
@testset "Asymmetric diagonal, endpoint limits, and zero mass" begin
    A=Matrix(Diagonal([1.0,2.0])); E=[0.0,2.0]; G=stackmat([A,zero(A)])
    T=build_tree(E,G;ghosts=(-1.0,3.0),leaf_panels=1)
    R,L=real_boundary(T,0.0;tau_far=1e-12)
    @test R ≈ A*log(0.5)/(2pi)-0.5im*A atol=2e-14
    @test L.covered_panels==3 && !L.formal_whole_certificate
    for e in (-1.0,-0.5,0.0,2.0,3.0,-1.0+eps(),nextfloat(0.0),prevfloat(2.0))
        r,_=real_boundary(T,e;tau_far=1e-12)
        # Ghost endpoints have zero Gamma, so hat-reference imaginary part is zero.
        @test r ≈ hat_reference(E,G,e,(-1.0,3.0);real_pv=true) atol=1e-12
        @test all(isfinite,r)
    end
    Z=build_tree(E,zero(G);ghosts=(-1.0,3.0),leaf_panels=1)
    for z in (0.0,-1.0,4.0)
        r,l=real_boundary(Z,z); @test iszero(norm(r)); @test l.covered_panels==3
    end
    r,l=complex_cauchy(Z,0.1+0.3im); @test iszero(norm(r)) && l.far_truncation_estimate==0
end
@testset "Noncommuting PSD hats, translations, and phase" begin
    A=ComplexF64[2 0.3+0.2im;0.3-0.2im 1]
    B=ComplexF64[1 -0.2+0.1im;-0.2-0.1im 2]
    @test eigmin(Hermitian(A))>0 && eigmin(Hermitian(B))>0 && norm(A*B-B*A)>0.1
    E=collect(range(0.0,1.0;length=65)).^1.7 .* 3 .- 1.2
    G=stackmat([(1+0.2sin(3e))*A+(0.5+0.1cos(5e))*B for e in E])
    ghosts=(E[1]-0.07,E[end]+0.31)
    T=build_tree(E,G;ghosts,leaf_panels=2,pmax=24)
    @test T.geometry_witness.max_depth>1
    @test T.geometry_witness.max_adjacent_gap_ratio>1
    root=T.nodes[T.root]; S0,S1=raw_moments(T)
    @test reshape(root.moments[1,:],2,2) ≈ S0 atol=2e-13
    @test root.center*S0+root.radius*reshape(root.moments[2,:],2,2) ≈ S1 atol=3e-13
    # Direct Gauss leaf integration across the root bypasses ALL translations.
    tq,wq=CauchyPanelTree.gauss_rule(13)
    direct=CauchyPanelTree.leaf_moments(T.x,T.gamma,1,length(T.x)-1,
        root.center,root.radius,24,tq,wq)
    @test root.moments ≈ direct atol=2e-12
    @test maximum(norm(root.moments[q,:]) for q in 1:24)<=root.mass+1e-12
    for z in (-0.1+0.03im,0.3+0.4im,8.0+7im)
        r,l=complex_cauchy(T,z;tau_far=1e-9)
        ref=hat_reference(E,G,z,ghosts)
        @test norm(r-ref)<=l.far_truncation_estimate+3e-11
        @test eigmin(Hermitian(-(r-r')/(2im))) >= -l.far_truncation_estimate-3e-11
        @test l.covered_panels==66 && l.max_far_radius_over_distance<=0.5
        @test l.far_truncation_estimate<=1e-9*(1+1e-12)
        rr,_=complex_cauchy(T,z;tau_far=1e-9,reverse_order=true)
        @test norm(rr-r)<=2e-9+1e-11
    end
    r,l=complex_cauchy(T,8.0+7im;tau_far=1e-9)
    @test l.far_panels>0  # Exercise multipoles, not merely a direct reference path.
    for e in (E[20],(E[20]+E[21])/2,E[20]+1e-10,E[20]-1e-10)
        r,l=real_boundary(T,e;tau_far=1e-9)
        @test norm(r-hat_reference(E,G,e,ghosts;real_pv=true))<=l.far_truncation_estimate+3e-10
    end
    # Affine energy covariance: Gamma values unchanged, Sigma unchanged at z'=a*z+b.
    a,b=0.037,1.13
    U=build_tree(a.*E.+b,G;ghosts=(a*ghosts[1]+b,a*ghosts[2]+b),leaf_panels=2)
    z=0.2+0.17im; r,_=complex_cauchy(T,z); s,_=complex_cauchy(U,a*z+b)
    @test r ≈ s atol=3e-9
    z=1e5+2e5im; r,_=complex_cauchy(T,z;tau_far=1e-13)
    @test norm(z*r-S0-S1/z)<1e-8
end
@testset "Five-basis column layout and raw-input preservation" begin
    P=reshape(ComplexF64.(1:15),5,3) .+ im.*reshape(ComplexF64.(15:-1:1),5,3)
    A=P*P'/100+I; D=Diagonal(cis.(collect(1:5))); B=D*A*D'+I
    E=[-0.8,-0.2,1.0]; G=stackmat([A,B,A+B])
    T=build_tree(E,G;ghosts=(-1.2,1.7),leaf_panels=1,hermitian_atol=1e-14)
    @test T.gamma[2:end-1,:]==reshape(G,3,25)
    r,l=complex_cauchy(T,0.17+0.29im;tau_far=1e-10)
    @test size(r)==(5,5)
    @test norm(r-hat_reference(E,G,0.17+0.29im,(-1.2,1.7)))<=l.far_truncation_estimate+1e-11
    @test !T.input_witness.representation_changed
end
@testset "Input and target contracts" begin
    G=ones(ComplexF64,2,1,1); E=[0.0,1.0]
    @test_throws Exception build_tree(E,G;ghosts=(-1.0,2.0),pmax=25)
    @test_throws Exception build_tree([0.0,0.0],G;ghosts=(-1.0,2.0))
    @test_throws Exception build_tree(E,-G;ghosts=(-1.0,2.0))
    @test_throws Exception build_tree(E,im.*G;ghosts=(-1.0,2.0))
    @test_throws Exception build_tree(E,fill(NaN+0im,2,1,1);ghosts=(-1.0,2.0))
    T=build_tree(E,G;ghosts=(-1.0,2.0))
    @test_throws Exception real_boundary(T,NaN)
    @test_throws Exception real_boundary(T,0.2;tau_far=0.0)
    @test_throws Exception complex_cauchy(T,0.2-0.1im)
    @test_throws Exception complex_cauchy(T,Inf+0.1im)
end
