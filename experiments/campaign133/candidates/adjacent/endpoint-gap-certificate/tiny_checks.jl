# Prospective tiny only; no actual Gamma/H5 reads or run on include.
using Test,LinearAlgebra
include("endpoint_certificate.jl")
function tiny_checks()
    before=BLAS.get_num_threads();BLAS.set_num_threads(1)
    println("tiny_BLAS_threads_before_setter=",before," after_setter=",BLAS.get_num_threads());@test BLAS.get_num_threads()==1
    @testset "outward basic operations and exact rational constants" begin
        for a in (0.0,-0.0,1.2,-2.7,1e-300),b in (.3,-1.7,0.0)
            x=iadd(point(a),point(b));@test Q(x[1])<=Q(a)+Q(b)<=Q(x[2])
            x=imul(point(a),point(b));@test Q(x[1])<=Q(a)*Q(b)<=Q(x[2])
        end
        pi,k1=constants()
        @test 333//106<Q(pi[1])<Q(pi[2])<355//113
        setprecision(256) do
            @test BigFloat(pi[1])<big(π)<BigFloat(pi[2])
            for n in (1,2,8,30704)
                box=kernel(point(n),pi,k1;one=n==1)
                F(x)=iszero(x) ? big(0) : x*log(x)
                t=BigFloat(n);exact=(F(t+1)-2F(t)+F(t-1))/(2big(π))
                @test BigFloat(box[1])<=exact<=BigFloat(box[2])
            end
        end
        @test Q(pi[1])>3 && k1[1]>0
    end
    @testset "complex Gershgorin sign and ghost area/factor" begin
        H=ComplexF64[3 .5im;-.5im 4];S=ComplexF64[.2 0;0 .2]
        d=gap(H,S,S,point(0),1);@test 2<d<eigmin(Hermitian(H-S))
        @test gap(zeros(ComplexF64,2,2),zeros(ComplexF64,2,2),zeros(ComplexF64,2,2),point(0),1)<0
        h=1//8;c=3//7;@test area(c,h,c)==h/2 && area(c,h,c+h)==h
        delta=1//big(2)^52;@test 0<area(c+h-delta,h,c)==delta^2/(2h)
        pi,_=constants();b=ghostbound(.5,0.0,2.0,3.0,pi)
        @test b>=Float64(inv(16big(π))) && b<Float64(inv(16big(π)))*(1+1e-13)
        @test_throws ErrorException ghostbound(.5,0.0,-2.0,3.0,pi)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    ARGS==["--tiny"] || error("usage: --tiny");tiny_checks()
end
