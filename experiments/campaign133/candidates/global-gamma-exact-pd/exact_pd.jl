# Exact certificate of stored binary Float64 values, not an error-bound oracle.
using SHA, TOML, HDF5, Dates, Test, LinearAlgebra
const GI = Complex{BigInt}
const QR = Rational{BigInt}
const QC = Complex{QR}
samebits(x::Float64,y::Float64) = (iszero(x)&&iszero(y)) || reinterpret(UInt64,x)==reinterpret(UInt64,y)

"IEEE binary64 = signed integer mantissa * 2^exponent; signed zeros are exact zero."
function dyadic(x::Float64)
    u=reinterpret(UInt64,x);ef=Int((u>>52)&0x7ff);f=u&0x000fffffffffffff
    ef==0x7ff && throw(ArgumentError("nonfinite IEEE input"))
    m=ef==0 ? Int64(f) : Int64(f|0x0010000000000000)
    (u>>63)!=0 && (m=-m)
    return m,ef==0 ? -1074 : ef-1075
end

mutable struct Scratch
    X::Matrix{ComplexF64}
    B::Matrix{GI}
    pivots::Vector{BigInt}
    exponent::Int
    divisions::Int
end
Scratch()=Scratch(zeros(ComplexF64,5,5),zeros(GI,5,5),zeros(BigInt,5),0,0)

"No projection. Exact binary equality, except equivalent signed zeros."
function exact_hermitian(X)
    n=size(X,1)
    for j in 1:n,i in 1:n
        samebits(real(X[i,j]),real(X[j,i])) || return false
        samebits(imag(X[i,j]),-imag(X[j,i])) || return false
    end
    return true
end

function integerize!(s::Scratch)
    emin=typemax(Int)
    for z in s.X,x in (real(z),imag(z))
        m,e=dyadic(x);iszero(m) || (emin=min(emin,e))
    end
    emin==typemax(Int) && (emin=0)
    s.exponent=emin
    for j in 1:5,i in 1:5
        mr,er=dyadic(real(s.X[i,j]));mi,ei=dyadic(imag(s.X[i,j]))
        r=iszero(mr) ? BigInt(0) : BigInt(mr)<<(er-emin)
        q=iszero(mi) ? BigInt(0) : BigInt(mi)<<(ei-emin)
        s.B[i,j]=GI(r,q)
    end
    return s
end

"Bareiss: diagonal stage k is det(M[1:k,1:k]); all divisions are checked."
function certify!(s::Scratch)
    fill!(s.pivots,0);s.divisions=0;s.exponent=0
    all(isfinite,s.X) || return (:nonfinite,0)
    exact_hermitian(s.X) || return (:nonhermitian,0)
    integerize!(s);previous=BigInt(1)
    for k in 1:5
        p=s.B[k,k]
        iszero(imag(p)) || error("INTERNAL: nonreal exact Bareiss pivot")
        pivot=real(p);s.pivots[k]=pivot
        pivot>0 || return (:nonpositive,k)
        for j in k+1:5,i in k+1:5
            v=pivot*s.B[i,j]-s.B[i,k]*s.B[k,j]
            rr,remr=divrem(real(v),previous);ii,remi=divrem(imag(v),previous)
            s.divisions+=2
            iszero(remr)&&iszero(remi) || error("INTERNAL: nonexact Bareiss division")
            s.B[i,j]=GI(rr,ii)
        end
        previous=pivot
    end
    return (:pd,5)
end

"Compare a*2^ea and b*2^eb exactly (positive determinant witnesses)."
function dyadic_less(a::BigInt,ea::Int,b::BigInt,eb::Int)
    e=min(ea,eb);return (a<<(ea-e)) < (b<<(eb-e))
end
function witness(s::Scratch,e,k,status,stage)
    Dict{String,Any}("energy_index"=>e,"k_index"=>k,"status"=>string(status),
        "stage"=>stage,"integer_scale_exponent"=>s.exponent,
        "leading_minors_integer"=>string.(s.pivots[1:stage]),
        "physical_minor_rule"=>"Delta_k = integer_minor_k * 2^(k*integer_scale_exponent)",
        "Gamma_real"=>[real(s.X[i,j]) for i in 1:5 for j in 1:5],
        "Gamma_imag"=>[imag(s.X[i,j]) for i in 1:5 for j in 1:5],
        "matrix_serialization"=>"row-major, five rows, no symmetrization")
end

# Tiny oracle only: a different algorithm with exact rational LDL reconstruction.
function rational_ldl(X::Matrix{ComplexF64})
    Q=QC.(QR.(real.(X)),QR.(imag.(X)));L=zeros(QC,5,5);d=zeros(QR,5)
    Q==Q' || return false,d
    for j in 1:5
        L[j,j]=1;p=Q[j,j]
        for k in 1:j-1;p-=L[j,k]*d[k]*conj(L[j,k]);end
        iszero(imag(p)) || error("nonreal exact rational pivot")
        d[j]=real(p);d[j]>0 || return false,d
        for i in j+1:5
            z=Q[i,j]
            for k in 1:j-1;z-=L[i,k]*d[k]*conj(L[j,k]);end
            L[i,j]=z/d[j]
        end
    end
    @assert L*Diagonal(d)*L'==Q
    return true,d
end
function compare_reference(X)
    s=Scratch();s.X.=X;status,_=certify!(s);pd,d=rational_ldl(X)
    @test (status==:pd)==pd
    if pd
        determinant=QR(1)
        for k in 1:5
            determinant*=d[k];exponent=k*s.exponent
            exact=exponent>=0 ? QR(s.pivots[k]<<exponent) : s.pivots[k]//(BigInt(1)<<(-exponent))
            @test exact==determinant
        end
        @test s.divisions==60
    end
    return status
end

function tiny_checks()
    @testset "exact IEEE / Bareiss versus rational LDL" begin
        # Fixed complex manufactured positive Gram; no numerical eigensolver oracle.
        L=ComplexF64[1 0 0 0 0; .5im 1 0 0 0; .25 -.5im 1 0 0;
            .125 .25im -.5 1 0; .25im -.125 .5im .25 1]
        X=L*Diagonal([1.,2.,4.,8.,16.])*L'
        @test compare_reference(X)==:pd
        @test compare_reference(2.0^(-200)*X)==:pd
        @test compare_reference(2.0^200*X)==:pd
        Y=copy(X);Y[1,1]=-1;@test compare_reference(Y)==:nonpositive
        @test compare_reference(zeros(ComplexF64,5,5))==:nonpositive
        Y=copy(X);Y[1,2]+=2.0^-40;s=Scratch();s.X.=Y;@test certify!(s)[1]==:nonhermitian
        Y=copy(X);Y[1,1]=Inf;s.X.=Y;@test certify!(s)[1]==:nonfinite
        Y=Matrix{ComplexF64}(I,5,5)*nextfloat(0.0)
        @test compare_reference(Y)==:pd
        Y=zeros(ComplexF64,5,5);for i in 1:5;Y[i,i]=1;end
        Y[1,2]=complex(0.,-0.);Y[2,1]=complex(-0.,0.);@test compare_reference(Y)==:pd
        @test dyadic(nextfloat(0.0))==(Int64(1),-1074)
        @test dyadic(-0.0)[1]==0
        @test_throws ArgumentError dyadic(NaN)
        @test_throws ArgumentError dyadic(Inf)
        @test dyadic_less(BigInt(1),-10,BigInt(3),-11)
    end
end

digest(path)=bytes2hex(open(sha256,path)) # streaming whole-file SHA, before HDF reading
function checked(path,sha)
    occursin(r"^[0-9a-f]{64}$",sha) || error("unfrozen source SHA: $path")
    digest(path)==sha || error("source SHA mismatch: $path")
    return path
end
function selected_fixture_checks(cfg,base)
    path=checked(normpath(joinpath(base,cfg["fixture_path"])),cfg["fixture_sha256"])
    @testset "frozen selected Gamma fixtures vs rational LDL" begin
        h5open(path,"r") do f
            r=read(f["Gamma_real"]);q=read(f["Gamma_imag"])
            @test size(r)==size(q)==(4,11,5,5)
            s=Scratch()
            for e in 1:4,k in 1:11
                for j in 1:5,i in 1:5;s.X[i,j]=complex(r[e,k,i,j],q[e,k,i,j]);end
                compare_reference(s.X)
            end
        end
    end
end
