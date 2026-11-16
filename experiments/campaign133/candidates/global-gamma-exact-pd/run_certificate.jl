# Standalone, bounded, read-only HDF consumer. No Core/Runner/Pkg imports.
include("exact_pd.jl")

mutable struct Worker
    scratch::Scratch
    counts::Dict{Symbol,Int} # one table per worker, never per block
    divisions::Int
    max_integer_bits::Int
    first::Union{Nothing,Dict{String,Any}}
    failures::Vector{Dict{String,Any}}
    worst_integer::BigInt
    worst_exponent::Int
    worst::Union{Nothing,Scratch} # owning fixed-size cache; serialize only at chunk boundary
    worst_energy::Int
    worst_k::Int
end
Worker()=Worker(Scratch(),Dict(:pd=>0,:nonpositive=>0,:nonhermitian=>0,:nonfinite=>0),
    0,0,nothing,Dict{String,Any}[],BigInt(0),0,nothing,0,0)

function consume!(w::Worker,r,q,local_e,k,global_e)
    s=w.scratch
    for j in 1:5,i in 1:5;s.X[i,j]=complex(r[local_e,k,i,j],q[local_e,k,i,j]);end
    status,stage=certify!(s);w.counts[status]+=1;w.divisions+=s.divisions
    global_e==1&&k==1 && (w.first=witness(s,global_e,k,status,stage))
    if status==:pd
        w.max_integer_bits=max(w.max_integer_bits,ndigits(s.pivots[5],base=2))
        exponent=5*s.exponent
        if isnothing(w.worst)||dyadic_less(s.pivots[5],exponent,w.worst_integer,w.worst_exponent)
            w.worst_integer=s.pivots[5];w.worst_exponent=exponent
            isnothing(w.worst)&&(w.worst=Scratch())
            copyto!(w.worst.X,s.X);copyto!(w.worst.pivots,s.pivots)
            w.worst.exponent=s.exponent;w.worst_energy=global_e;w.worst_k=k
        end
    elseif length(w.failures)<8
        push!(w.failures,witness(s,global_e,k,status,stage))
    end
    return nothing
end

function save_summary(path,record,cap)
    tmp=path*".tmp"
    open(tmp,"w") do io;TOML.print(io,record);end
    filesize(tmp)<=cap || error("output cap; oversized tmp retained")
    mv(tmp,path;force=true)
end
function aggregate(workers)
    counts=Dict(string(s)=>sum(w.counts[s] for w in workers) for s in keys(workers[1].counts))
    witnesses=Dict{String,Any}[];first=Dict{String,Any}[]
    for w in workers
        append!(witnesses,w.failures);isnothing(w.first)||push!(first,w.first)
    end
    sort!(witnesses;by=r->(r["energy_index"],r["k_index"]))
    resize!(witnesses,min(8,length(witnesses)))
    worst=nothing
    for w in workers
        isnothing(w.worst)&&continue
        if isnothing(worst)||dyadic_less(w.worst_integer,w.worst_exponent,worst.worst_integer,worst.worst_exponent)
            worst=w
        end
    end
    Dict{String,Any}("counts"=>counts,"checked_blocks"=>sum(values(counts)),
        "checked_exact_divrem_components"=>sum(w.divisions for w in workers),
        "max_integer_determinant_bits"=>maximum(w.max_integer_bits for w in workers),
        "first_block"=>first,"failure_witnesses_first8"=>witnesses,
        "smallest_physical_determinant_positive_block"=>(isnothing(worst) ? Dict{String,Any}() :
            witness(worst.worst,worst.worst_energy,worst.worst_k,:pd,5)),
        "worst_definition"=>"minimum exact positive Delta5 in ORIGINAL stored units; not a minimum eigenvalue or condition estimate")
end

function source_gate(cfg,base)
    expected=["gamma","model","input","identity","map","operator","spectral","ledger","summary"]
    sort(collect(keys(cfg["bindings"])))==sort(expected) || error("all source bindings required")
    paths=Dict{String,String}();hashes=Dict{String,String}()
    # ALL whole-file SHA checks precede first HDF open. No source arrays imported here.
    for label in expected
        b=cfg["bindings"][label];path=normpath(joinpath(base,b["path"]))
        paths[label]=checked(path,b["sha256"]);hashes[label]=b["sha256"]
    end
    input=TOML.parsefile(paths["input"]);identity=TOML.parsefile(paths["identity"])
    map=TOML.parsefile(paths["map"]);summary=TOML.parsefile(paths["summary"])
    ledger=TOML.parsefile(paths["ledger"])
    input["identity"]==cfg["source_identity"] || error("source input identity")
    input["pi_calls_total"]==cfg["source_call"]==64 || error("source call cap")
    input["project_AA_trial"]===false || error("source was projected")
    input["spectral_operator"]==cfg["source_spectral_operator"] || error("source operator policy")
    identity["input_sha256"]==hashes["input"] || error("identity/input binding")
    identity["source_sha256"]["gamma_functions.jl"]==hashes["operator"] || error("operator/identity binding")
    identity["source_sha256"]["stable_spectral.jl"]==hashes["spectral"] || error("spectral/identity binding")
    map["call"]==cfg["source_call"] || error("map call")
    map["spectral_operator_ledger"]["operator"]==cfg["source_spectral_operator"] || error("map operator")
    ledger==map["spectral_operator_ledger"] || error("standalone spectral ledger differs from map")
    summary["execution"]=="complete"&&summary["pi_calls"]==cfg["source_call"] || error("source incomplete")
    summary["reduced"]["last_map"]==map || error("source last accepted map differs")
    (input["NE"],input["Nk"],input["Nb"],input["Nb"])==Tuple(cfg["shape"]) || error("source dimensions")
    return paths,hashes
end

function tiny_runner_checks()
    @testset "owning worst cache / bounded record TOML" begin
        w=Worker();r=zeros(Float64,1,1,5,5);q=similar(r);fill!(q,0)
        for a in 1:5;r[1,1,a,a]=2;end
        consume!(w,r,q,1,1,1)
        for a in 1:5;r[1,1,a,a]=1;end
        consume!(w,r,q,1,1,2)
        fill!(r,0);consume!(w,r,q,1,1,3) # overwrites working scratch, not cached best
        rec=aggregate([w])
        @test rec["counts"]["pd"]==2&&rec["counts"]["nonpositive"]==1
        @test rec["smallest_physical_determinant_positive_block"]["energy_index"]==2
        @test rec["smallest_physical_determinant_positive_block"]["leading_minors_integer"][5]==string(BigInt(1)<<260)
        io=IOBuffer();TOML.print(io,rec);@test TOML.parse(String(take!(io)))==rec
    end
end

function run_certificate(input_path,output)
    started=time();cfg=TOML.parsefile(input_path);base=dirname(abspath(input_path))
    cfg["runtime_authorized"]===true || error("runtime not authorized")
    string(VERSION)==cfg["expected_julia_version"] || error("Julia version")
    Threads.nthreads()<=4 || error("maximum four workers")
    mode=cfg["mode"];mode in ("benchmark5000","all") || error("unknown mode")
    shape=Tuple(Int.(cfg["shape"]));shape==(30704,11,5,5) || error("fixed source shape")
    cfg["chunk_energy"] in 1:1024 || error("chunk bound")
    0<cfg["max_seconds"]<=280 || error("internal wall bound")
    0<cfg["output_cap_bytes"]<=32*1024^2 || error("output bound")
    limit=mode=="all" ? shape[1]*shape[2] : 5000
    paths,hashes=source_gate(cfg,base);workers=[Worker() for _ in 1:Threads.nthreads()]
    mkpath(output);isempty(readdir(output)) || error("new empty output directory required")
    summary=Dict{String,Any}("identity"=>cfg["identity"],"execution"=>"running",
        "mode"=>mode,"expected_all_blocks"=>shape[1]*shape[2],"requested_blocks"=>limit,
        "input_sha256"=>digest(input_path),"script_sha256"=>digest(@__FILE__),
        "math_sha256"=>digest(joinpath(@__DIR__,"exact_pd.jl")),"source_sha256"=>hashes,
        "julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),
        "start_utc"=>string(now(UTC)),"all_blocks_strict_PD_certified"=>false,
        "SCBA_convergence"=>"not_measured","native_acceptance"=>"not_measured",
        "causal_Mats_density_acceptance"=>"not_measured; PD alone is only one prerequisite",
        "scope"=>"exact stored UNPROJECTED GammaLast64 only; intended Hermitian hH(all9) future variant; not Gamma128, no forward/grid accuracy certificate")
    complete=false;stop_reason="requested block limit"
    try
        h5open(paths["gamma"],"r") do f
            read(f["schema"])=="campaign133-gamma-checkpoint-v1" || error("checkpoint schema")
            read(f["identity"])==cfg["source_identity"] || error("checkpoint identity")
            read(f["pi_calls_total"])==cfg["source_call"] || error("checkpoint call")
            rds=f["Gamma_real"];qds=f["Gamma_imag"]
            size(rds)==size(qds)==shape || error("checkpoint axis shape")
            eltype(rds)==eltype(qds)==Float64 || error("requires original Float64")
            for first in 1:cfg["chunk_energy"]:shape[1]
                if time()-started>cfg["max_seconds"]
                    stop_reason="internal wall cap; incomplete is not ALL certificate";break
                end
                last=min(first+cfg["chunk_energy"]-1,shape[1],cld(limit,shape[2]))
                first>last&&break
                r=rds[first:last,:,:,:];q=qds[first:last,:,:,:]
                count=min((last-first+1)*shape[2],limit-(first-1)*shape[2])
                Threads.@threads :static for worker in eachindex(workers)
                    w=workers[worker]
                    for index in worker:length(workers):count
                        local_e=fld(index-1,shape[2])+1;k=mod(index-1,shape[2])+1
                        consume!(w,r,q,local_e,k,first+local_e-1)
                    end
                end
                summary["results"]=aggregate(workers);summary["elapsed_seconds"]=time()-started
                save_summary(joinpath(output,"summary.toml"),summary,cfg["output_cap_bytes"])
                if summary["results"]["checked_blocks"]==limit;complete=true;break;end
            end
        end
        summary["execution"]=complete ? "complete" : "bounded_incomplete"
    catch e
        summary["execution"]="failed";summary["failure"]=sprint(showerror,e)
        rethrow()
    finally
        summary["results"]=aggregate(workers);summary["elapsed_seconds"]=time()-started
        summary["stop_reason"]=stop_reason
        allpd=summary["results"]["counts"]["pd"]==shape[1]*shape[2]
        summary["all_blocks_strict_PD_certified"]=complete&&mode=="all"&&allpd&&summary["execution"]=="complete"
        summary["requested_subset_strict_PD_certified"]=complete&&summary["results"]["counts"]["pd"]==limit&&summary["execution"]=="complete"
        save_summary(joinpath(output,"summary.toml"),summary,cfg["output_cap_bytes"])
    end
    return summary
end

if abspath(PROGRAM_FILE)==@__FILE__
    if length(ARGS)==1&&ARGS[1]=="--tiny"
        tiny_checks();tiny_runner_checks()
    elseif length(ARGS)==2&&ARGS[1]=="--fixtures"
        cfg=TOML.parsefile(ARGS[2]);tiny_checks();tiny_runner_checks();selected_fixture_checks(cfg,dirname(abspath(ARGS[2])))
    elseif length(ARGS)==2
        run_certificate(abspath(ARGS[1]),abspath(ARGS[2]))
    else
        error("usage: run_certificate.jl --tiny | --fixtures input.toml | input.toml new-output-dir")
    end
end
