# Selected 5x5 diagnostic. Include defines functions; execution requires Root's
# authorized input and external 2GiB/120s scope. No Core/Runner/solver imports.
using LinearAlgebra, HDF5, TOML, SHA, Dates, Test
const RQ = Rational{BigInt}
const CQ = Complex{RQ}
digest(path) = bytes2hex(open(sha256,path))
exact_sign(x) = x < 0 ? -1 : x > 0 ? 1 : 0

"Exact binary Float64 value, without decimal rationalization or rounding."
function binary_rational(x::Float64)::RQ
    isfinite(x) || throw(ArgumentError("nonfinite Float64"))
    bits = reinterpret(UInt64,x)
    exponent = Int((bits >> 52) & 0x7ff)
    fraction = bits & 0x000fffffffffffff
    significand = exponent == 0 ? BigInt(fraction) : BigInt(fraction | 0x0010000000000000)
    (bits >> 63) != 0 && (significand = -significand)
    shift = exponent == 0 ? -1074 : exponent - 1023 - 52
    return shift >= 0 ? (significand << shift)//BigInt(1) : significand//(BigInt(1) << (-shift))
end
exact_complex(z::ComplexF64) = CQ(binary_rational(real(z)),binary_rational(imag(z)))
function guard_hermitian(X::Matrix{ComplexF64})
    n=size(X,1);size(X,2)==n || error("square block required")
    B=similar(X);sq=0.0;defect=0.0
    # Same b,a order and expression as frozen cone: block=(X+Xdagger)/2.
    for b in 1:n,a in 1:n
        z=X[a,b];sq+=abs2(z);defect+=abs2(z-conj(X[b,a]))
        B[a,b]=(z+conj(X[b,a]))/2
    end
    return B,sqrt(sq),sqrt(defect)
end
function exact_hermitian_raw(X::Matrix{ComplexF64})
    Q=exact_complex.(X);B=similar(Q);half=BigInt(1)//BigInt(2)
    for b in axes(Q,2),a in axes(Q,1)
        B[a,b]=(Q[a,b]+conj(Q[b,a]))*half
    end
    return B
end
function exact_rayleigh(B::Matrix{CQ},v::Vector{ComplexF64})
    qv=exact_complex.(v);n=length(v);z=zero(CQ);vnorm=zero(RQ)
    for a in 1:n
        vnorm+=real(conj(qv[a])*qv[a])
        for b in 1:n;z+=conj(qv[a])*B[a,b]*qv[b];end
    end
    vnorm>0 || error("zero eigenvector")
    iszero(imag(z)) || error("exact Hermitian Rayleigh value has imaginary part")
    return real(z)/vnorm,vnorm
end
function rational_record(q::RQ,vnorm::RQ,allowance::Float64)
    decimal256=setprecision(BigFloat,256) do
        string(BigFloat(numerator(q))/BigFloat(denominator(q)))
    end
    return Dict{String,Any}("sign"=>exact_sign(q),"numerator"=>string(numerator(q)),
        "denominator"=>string(denominator(q)),"vector_norm2_numerator"=>string(numerator(vnorm)),
        "vector_norm2_denominator"=>string(denominator(vnorm)),
        "decimal256"=>decimal256,
        "Float64_display"=>Float64(q),"below_negative_declared_allowance"=>q < -binary_rational(allowance),
        "claim"=>"exact sign of this rational quadratic form; no eigenvalue enclosure or whole-Pi accuracy certificate")
end
function read_block(file,e,k,expected_shape)
    re,impart=file["Gamma_real"],file["Gamma_imag"]
    size(re)==size(impart)==expected_shape || error("checkpoint logical axes differ")
    1<=e<=expected_shape[1] && 1<=k<=expected_shape[2] || error("selected index outside axes")
    # Julia HDF API logical(E,k,a,b); the disk C order is(b,a,k,E).
    return Matrix(reshape(complex.(re[e:e,k:k,:,:],impart[e:e,k:k,:,:]),expected_shape[3],expected_shape[4]))
end
function check_digest(path,expected)
    occursin(r"^[0-9a-f]{64}$",expected) || error("Root must freeze a SHA256 before runtime: $path")
    observed=digest(path);observed==expected || error("source SHA mismatch: $path")
    return observed
end
function analyze_block(X,declared_allowance)
    all(isfinite,X) || error("nonfinite selected block")
    B,scale,defect=guard_hermitian(X)
    # Accepted Core e171, positivity_diagnostics.jl: PSD metric v3. The cone
    # call uses Frobenius block scale and construction_scale=1, error=0.
    allowance=64.0*size(B,1)*eps(Float64)*max(scale,1.0)
    allowance==declared_allowance || error("declared allowance differs from frozen cone formula")
    ev=eigen(Hermitian(B,:U));j=argmin(ev.values);lambda=ev.values[j];v=Vector(ev.vectors[:,j])
    guardq,vn=exact_rayleigh(exact_complex.(B),v)
    rawq,_=exact_rayleigh(exact_hermitian_raw(X),v)
    residual=norm(B*v-lambda*v);residual_scale=(norm(B)+abs(lambda))*norm(v)
    r=Dict{String,Any}("raw_Frobenius_scale"=>scale,"raw_antiHermitian_F"=>defect,
        "guard_allowance"=>allowance,"guard_allowance_formula"=>"64*n*eps(Float64)*max(Frobenius(X),1), accepted Core e171 metricv3",
        "float_min_eigenvalue"=>lambda,"float_Rayleigh"=>real(dot(v,B*v))/real(dot(v,v)),
        "eigenvector_norm"=>norm(v),"eigenpair_residual_absolute"=>residual,
        "eigenpair_residual_scale"=>residual_scale,
        "eigenpair_residual_relative"=>(residual_scale==0 ? (residual==0 ? 0.0 : Inf) : residual/residual_scale),
        "exact_guard_averaged_matrix"=>rational_record(guardq,vn,allowance),
        "exact_Hermitian_part_raw_entries"=>rational_record(rawq,vn,allowance),
        "matrix_averaging_difference_on_vector"=>rational_record(rawq-guardq,vn,allowance))
    return r,B,v
end
function tiny_checks()
    @testset "exact binary Rayleigh and selected HDF" begin
        @test binary_rational(0.5)==BigInt(1)//BigInt(2)
        @test binary_rational(-0.0)==0
        @test binary_rational(nextfloat(1.0))==(BigInt(2)^52+1)//(BigInt(2)^52)
        @test binary_rational(nextfloat(0.0))==BigInt(1)//(BigInt(2)^1074)
        @test_throws ArgumentError binary_rational(Inf)
        B=ComplexF64[1 im;-im 0];v=ComplexF64[1,im]
        q,vn=exact_rayleigh(exact_complex.(B),v)
        @test q==(-BigInt(1)//BigInt(2)) && vn==2
        @test exact_sign(q)==-1
        X=copy(B);X[1,2]+=1e-15
        F,_,_=guard_hermitian(X)
        @test F==F'
        qraw,_=exact_rayleigh(exact_hermitian_raw(X),v)
        @test qraw==q
        rec,_,_=analyze_block(B,64*2*eps(Float64)*max(norm(B),1))
        @test rec["exact_guard_averaged_matrix"]["sign"]==-1
        roundtrip=TOML.parse(sprint(io->TOML.print(io,rec)))["exact_guard_averaged_matrix"]
        original=rec["exact_guard_averaged_matrix"]
        @test roundtrip["numerator"]==original["numerator"]
        @test roundtrip["denominator"]==original["denominator"]
        mktempdir() do directory
            path=joinpath(directory,"axes.h5");A=zeros(ComplexF64,3,2,2,2)
            A[2,1,:,:]=B
            HDF5.h5open(path,"w") do f;f["Gamma_real"]=real.(A);f["Gamma_imag"]=imag.(A);end
            HDF5.h5open(path,"r") do f;@test read_block(f,2,1,size(A))==B;end
        end
    end
end
function run_diagnostic(input,output)
    cfg=TOML.parsefile(input);cfg["runtime_authorized"]===true || error("write-only freeze: runtime disabled")
    string(VERSION)==cfg["expected_julia_version"] && Threads.nthreads()==1 || error("Julia/thread identity differs")
    BLAS.set_num_threads(1)
    length(cfg["blocks"]) in 1:7 || error("only1..7 selected blocks admitted")
    ispath(output) && error("refusing to overwrite output")
    base=dirname(abspath(input));mkpath(output);cp(input,joinpath(output,"input.toml"));cp(@__FILE__,joinpath(output,"rayleigh.jl"))
    started=time();summary=Dict{String,Any}("execution"=>"running","script_sha256"=>digest(@__FILE__),
        "input_sha256"=>digest(input),"scope"=>"selected exact-rational Rayleigh signs only; no Pi/solver/acceptance",
        "units"=>"source scaled Gamma energy units; no spin or quadrature weights in a local matrix Rayleigh quotient",
        "started_utc"=>string(Dates.now(Dates.UTC)),"blocks"=>Dict{String,Any}())
    try
        tiny_checks();summary["tiny_checks"]="pass"
        sourceinput=normpath(joinpath(base,cfg["parent_input_path"]))
        check_digest(sourceinput,cfg["parent_input_sha256"])
        parent=TOML.parsefile(sourceinput);parent["identity"]==cfg["expected_parent_identity"] || error("parent identity differs")
        shape=(30704,11,5,5)
        for block in cfg["blocks"]
            time()-started<110 || error("internal110s diagnostic wall stop")
            occursin(r"^[A-Za-z0-9_-]+$",block["name"]) || error("invalid output block name")
            path=normpath(joinpath(base,block["source_path"]));attemptpath=normpath(joinpath(base,block["attempt_path"]))
            # Whole-file SHA streaming ALWAYS precedes HDF5 field access.
            check_digest(path,block["source_sha256"]);check_digest(attemptpath,block["attempt_sha256"])
            attempt=TOML.parsefile(attemptpath);attempt["call"]==block["call"] && attempt["kind"]=="AAII" || error("attempt role differs")
            match_reason=match(r"^negative PSD field at \((\d+),(\d+)\): (.*); budget=(.*)$",attempt["reason"])
            match_reason===nothing && error("attempt does not record a negative-cone rejection")
            e,k=block["energy_index"],block["k_index"]
            parse(Int,match_reason[1])==e && parse(Int,match_reason[2])==k || error("selected reason indices differ")
            declared=parse(Float64,match_reason[4])
            X=HDF5.h5open(path,"r") do f
                read(f["schema"])=="campaign133-gamma-checkpoint-v1" || error("checkpoint schema differs")
                read(f["identity"])==parent["identity"] && read(f["pi_calls_total"])==block["call"] || error("checkpoint provenance differs")
                read_block(f,e,k,shape)
            end
            result,B,v=analyze_block(X,declared)
            result["source_sha256"]=block["source_sha256"];result["attempt_sha256"]=block["attempt_sha256"]
            result["energy_index"]=e;result["k_index"]=k;result["declared_float_eigenvalue"]=parse(Float64,match_reason[3])
            result["AA_l1"]=attempt["AA"]["l1"];result["source_call"]=block["call"]
            summary["blocks"][block["name"]]=result
            HDF5.h5open(joinpath(output,block["name"]*".h5"),"w") do f
                for (name,A) in (("raw",X),("guard_Herm",B),("float_eigenvector",v))
                    f[name*"_real"]=real.(A);f[name*"_imag"]=imag.(A)
                end
            end
        end
        summary["execution"]="complete"
    catch error
        summary["execution"]="stopped";summary["reason"]=sprint(showerror,error,catch_backtrace());rethrow()
    finally
        summary["elapsed_seconds"]=time()-started
        open(joinpath(output,"summary.toml"),"w") do io;TOML.print(io,summary);end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    ARGS==["--tiny"] ? tiny_checks() : (length(ARGS)==2 ? run_diagnostic(abspath(ARGS[1]),abspath(ARGS[2])) : error("--tiny OR INPUT NEW_OUTPUT"))
end
