# Read-only snapshot diagnostic. No Core dependency or new Phi evaluation.
using HDF5, LinearAlgebra, TOML, SHA, Dates, Test
sha(p)=bytes2hex(open(sha256,p))
finite(xs...)=all(isfinite,xs)||error("nonfinite input/intermediate")
ratio(a,b)=iszero(b) ? (iszero(a) ? 0.0 : Inf) : a/b
fraction(a,b)=iszero(b) ? NaN : a/b
function replay_ok(value,ref,tol)
    isfinite(tol)&&tol>=0&&!isnan(value)&&!isnan(ref)||return false
    (isequal(value,ref)||iszero(value)&&iszero(ref))&&return true
    isfinite(value)&&isfinite(ref)||return false
    delta=abs(value-ref); limit=tol*abs(ref)
    return isfinite(delta)&&isfinite(limit)&&delta<=limit
end
fermi(x)=x>=0 ? exp(-x)/(1+exp(-x)) : inv(1+exp(x))
readrows(f,n,r)=complex.(f[n*"_real"][r,:,:,:],f[n*"_imag"][r,:,:,:])
blockstats(P,G)=(sum(abs2,P-G),sum(abs2,P))
function measures(r2,s2,E,we,wk,mu,kBT)
    W=we*transpose(wk); R=W.*r2; S=W.*s2; f=fermi.((E.-mu)./kBT)
    zones=[E.<mu-5kBT,abs.(E.-mu).<=5kBT,E.>mu+5kBT]
    norms=Dict("absolute"=>sqrt(sum(r2)),"scale"=>sqrt(sum(s2)),
        "weighted_absolute"=>sqrt(sum(R)),"weighted_scale"=>sqrt(sum(S)))
    norms["relative"]=ratio(norms["absolute"],norms["scale"])
    norms["weighted_relative"]=ratio(norms["weighted_absolute"],norms["weighted_scale"])
    occ=Dict("absolute"=>sqrt(sum(R.*f.^2)),"scale"=>sqrt(sum(S.*f.^2)))
    occ["relative"]=ratio(occ["absolute"],occ["scale"])
    fractions=[fraction(sum(R[z,:]),sum(R)) for z in zones]
    return (;R,S,f,norms,occ,fractions)
end
function fixture()
    @testset "source-free complex residual support" begin
        P=ComplexF64[1 im;-im 2]; r=zeros(3,2); s=similar(r)
        for e in 1:3,k in 1:2
            d=e+2(k-1); G=P-ComplexF64[d im*d;0 0]; r[e,k],s[e,k]=blockstats(P,G)
        end
        m=measures(r,s,[-1.,0.,1.],[1.,2.,4.],[.25,.75],0.,.1)
        @test sum(r)==128 && sum(s)==42 && sum(m.R)==234 && sum(m.S)==49
        @test m.fractions≈[14/234,52/234,168/234]
        @test m.occ["absolute"]^2≈14*fermi(-10.)^2+13+168*fermi(10.)^2
        @test ratio(0.,0.)==0 && ratio(1.,0.)==Inf && isnan(fraction(0.,0.))
        @test replay_ok(1.,1.,1e-10) && replay_ok(1+5e-11,1.,1e-10)
        @test !replay_ok(1+2e-10,1.,1e-10) && !replay_ok(1e-12,0.,1e-10)
        @test replay_ok(0.,0.,1e-10) && replay_ok(-0.,0.,1e-10)
        @test replay_ok(Inf,Inf,1e-10) && replay_ok(-Inf,-Inf,1e-10)
        @test !replay_ok(Inf,-Inf,1e-10) && !replay_ok(1.,Inf,1e-10) && !replay_ok(Inf,1.,1e-10)
        @test !replay_ok(NaN,NaN,1e-10) && !replay_ok(NaN,1.,1e-10) && !replay_ok(1.,NaN,1e-10)
        @test !replay_ok(1.,1.,Inf) && !replay_ok(1.,1.,NaN) && !replay_ok(1.,1.,-1.)
        @test !replay_ok(-floatmax(Float64),floatmax(Float64),1e-10)
    end
end
function main(input,output)
    c=TOML.parsefile(input); base=dirname(abspath(input)); paths=Dict{String,String}(); hashstart=time()
    finite(c["replay_relative"],c["seconds"],c["output_bytes"])
    for key in ("terminal","model","map","source_input","identity","writer")
        p=normpath(joinpath(base,c[key])); sha(p)==c[key*"_sha256"]||error("SHA mismatch: $key"); paths[key]=p
    end
    hashseconds=time()-hashstart; started=time(); BLAS.set_num_threads(1)
    Threads.nthreads()==1||error("requires Julia1"); ispath(output)&&error("output exists"); mkpath(output)
    cp(input,joinpath(output,"input.toml")); cp(@__FILE__,joinpath(output,"probe.jl"))
    checktime()=(time()-started<=c["seconds"]||error("phase/slab time cap"))
    src=TOML.parsefile(paths["source_input"]); id=TOML.parsefile(paths["identity"]); mp=TOML.parsefile(paths["map"])
    all(dirname(realpath(p))==dirname(realpath(paths["terminal"])) for p in values(paths))||error("attempt ownership")
    id["core_head"]==c["core_head"] && id["core_tree"]==c["core_tree"] && isempty(id["core_dirty"])||error("Core identity")
    id["input_sha256"]==c["source_input_sha256"] && id["source_sha256"]["gamma_functions.jl"]==c["writer_sha256"]||error("writer binding")
    src["retarded_real_symmetry"]=="hermitian_projection" && mp["call"]==128||error("wrong closure/map")
    ne,nk,nb=src["NE"],src["Nk"],src["Nb"]; E0=mp["absolute_E0_eV"]; L0=c["L0_m"]
    finite(E0,L0,src["temperature_K"]); E0>0 && L0>0 && src["temperature_K"]==70||error("units/bath")
    E,we,k,wk=h5open(paths["model"],"r") do f
        read(f["schema"])=="campaign133-gamma-model-v2"||error("model schema")
        alg=TOML.parse(read(f["algorithms_toml"]))
        alg["hilbert"]=="product_integration" && alg["embedding"]=="finite_chain" && alg["embedding_periods"]==4||error("model algorithms")
        read(f["grid_ε"]),read(f["grid_wᴱ"]),read(f["grid_κ"]),read(f["grid_wᵏ"])
    end
    length(E)==length(we)==ne && length(k)==length(wk)==nk||error("frame dimensions")
    all(isfinite,E)&&all(isfinite,we)&&all(isfinite,k)&&all(isfinite,wk)&&all(>(0),we)&&all(>=(0),wk)||error("frame values")
    expected=src["kmax_per_nm"]*1e9*L0; finite(expected); isapprox(maximum(k),expected;rtol=128eps(),atol=0)||error("L0/kappa")
    h=(last(E)-first(E))/(ne-1); finite(h); h>0||error("energy spacing")
    maximum(abs.(E.-[first(E)+(j-1)*h for j in 1:ne]))<=64eps()*max(1.,maximum(abs,E))||error("nonuniform E")
    kBT=(1.380649e-23/1.602176634e-19)*src["temperature_K"]/E0
    r2=zeros(ne,nk); s2=similar(r2); widths=similar(r2); summary=Dict{String,Any}()
    h5open(paths["terminal"],"r") do f
        read(f["schema"])=="campaign133-gamma-reduced-terminal-v1"||error("terminal schema")
        for name in ("Gamma_real","Gamma_imag","raw_Phi_real","raw_Phi_imag")
            size(f[name])==(ne,nk,nb,nb)||error("logical E,k,a,b axes")
        end
        mu=read(f["mu"]); finite(mu,kBT); mu==mp["mu"]["mu"]&&kBT>0||error("mu/bath")
        for firste in 1:1024:ne
            checktime(); rows=firste:min(firste+1023,ne); G=readrows(f,"Gamma",rows); P=readrows(f,"raw_Phi",rows)
            all(isfinite,G)&&all(isfinite,P)||error("nonfinite fields")
            for j in eachindex(rows),q in 1:nk
                e=rows[j]; r2[e,q],s2[e,q]=blockstats(Matrix(view(P,j,q,:,:)),Matrix(view(G,j,q,:,:)))
                widths[e,q]=real(tr(Matrix(view(G,j,q,:,:))))/nb*E0
            end
        end
        m=measures(r2,s2,E,we,wk,mu,kBT)
        finite(sum(r2),sum(s2),sum(m.R),sum(m.S))
        for (name,value) in m.norms
            ref=mp["residual"][name]; replay_ok(value,ref,c["replay_relative"])||error("raw norm replay: $name")
        end
        summary=Dict("status"=>"completed","raw_map_norms"=>m.norms,"f2_weighted_diagnostic"=>m.occ,
            "weighted_numerator_zone_fractions"=>m.fractions,"zone_fraction_status"=>iszero(sum(m.R)) ? "unavailable_zero_numerator" : "measured",
            "zones"=>["E<mu-5kBT","abs(E-mu)<=5kBT","E>mu+5kBT"],"mu_dimensionless"=>mu,"kBT_dimensionless"=>kBT,
            "E0_eV"=>E0,"L0_m"=>L0,"raw_squared_numerator"=>sum(r2),"raw_squared_scale"=>sum(s2),
            "weighted_squared_numerator"=>sum(m.R),"weighted_squared_scale"=>sum(m.S))
        summary["optional_saved_condition_available"]=haskey(f,"green/condition_number")
        summary["optional_saved_A_available"]=haskey(f,"green/A_real")&&haskey(f,"green/A_imag")
        summary["map_raw_norms"]=mp["residual"]
        order=partialsortperm(vec(m.R),1:min(20,length(m.R));rev=true)
        open(joinpath(output,"top20.csv"),"w") do io
            println(io,"rank,e,k,E_eV,kappa,f,weighted_r2,weighted_fraction,raw_r2,raw_ref2,Gamma_trace_perNb_eV,width_over_dE,HermGamma_min_eV,Gamma_antiHerm_F,condition,Atrace_dimensionless_real,Atrace_dimensionless_imag")
            for (rank,linear) in enumerate(order)
                e,q=Tuple(CartesianIndices(m.R)[linear]); G=complex.(f["Gamma_real"][e,q,:,:],f["Gamma_imag"][e,q,:,:])
                condition=haskey(f,"green/condition_number") ? f["green/condition_number"][e,q] : NaN
                at=haskey(f,"green/A_real")&&haskey(f,"green/A_imag") ? tr(complex.(f["green/A_real"][e,q,:,:],f["green/A_imag"][e,q,:,:])) : complex(NaN,NaN)
                println(io,join((rank,e,q,E[e]*E0,k[q],m.f[e],m.R[e,q],fraction(m.R[e,q],sum(m.R)),r2[e,q],s2[e,q],widths[e,q],widths[e,q]/(h*E0),eigmin(Hermitian((G+G')/2))*E0,norm(G-G'),condition,real(at),imag(at)),','))
            end
        end
        for (name,axis,reduce) in (("perE.csv",E*E0,X->vec(sum(X;dims=2))),("perK.csv",k,X->vec(sum(X;dims=1))))
            open(joinpath(output,name),"w") do io
                println(io,"coordinate,raw_r2,raw_ref2,weighted_r2,weighted_ref2,weighted_fraction,f2_weighted_r2,f2_weighted_ref2")
                cols=(axis,reduce(r2),reduce(s2),reduce(m.R),reduce(m.S),fraction.(reduce(m.R),sum(m.R)),reduce(m.R.*m.f.^2),reduce(m.S.*m.f.^2))
                for row in zip(cols...); println(io,join(row,',')); end
            end
        end
    end
    checktime(); summary["script_sha256"]=sha(@__FILE__); summary["input_sha256"]=sha(input)
    summary["source_hashes"]=Dict(key=>c[key*"_sha256"] for key in keys(paths)); summary["sha_seconds_separate"]=hashseconds
    summary["elapsed_seconds_excluding_SHA"]=time()-started; summary["scientific_scope"]="snapshot support only; no Jacobian/linewidth/SCBA convergence inference"
    open(joinpath(output,"summary.toml"),"w") do io; TOML.print(io,summary); end
    sum(filesize(joinpath(output,n)) for n in readdir(output))<=c["output_bytes"]||error("output cap")
end
if abspath(PROGRAM_FILE)==@__FILE__
    ARGS==["--fixture"] ? fixture() : (length(ARGS)==2 ? main(abspath(ARGS[1]),abspath(ARGS[2])) : error("--fixture or INPUT NEW_OUTPUT"))
end
