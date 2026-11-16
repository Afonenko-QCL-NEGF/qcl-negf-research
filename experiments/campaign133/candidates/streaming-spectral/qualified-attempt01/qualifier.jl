# WRITE ONLY; ROOT owns the external scope and runtime input.
using QCLNEGF, HDF5, LinearAlgebra, TOML, SHA, Test
const QSOURCE=normpath(joinpath(@__DIR__,"../equilibrium-gamma-stableA-cont64/attempt01"))
qsha(path)=bytes2hex(open(sha256,path))
const QDEFS=Dict("gamma_functions.jl"=>"50bcfaaeb3c399baaaa2a4b9b9e023167e97b4d87c6af95d83353cb453072211",
    "warm_seed.jl"=>"68fd8582334a1ed3f257efcc0b1096a429361dd73aa5be3be415dd79ddcbc157",
    "stable_spectral.jl"=>"eeca1f29710d9bb60bf6ce70ce0ebd1620fb86d4c77314b4bd8ffa9a5bcd5aa5",
    "gamma_restart.jl"=>"2a805be29aeef33d7a31356850b7fdaa80da7e68506d8f96b980a8108bd56a00")
for (name,sha) in QDEFS;qsha(joinpath(QSOURCE,name))==sha || error("old definition identity: "*name);end
qsha(joinpath(@__DIR__,"stream_spectral.jl"))=="8df7885cbf547ffd4077e611f6a88b933dc75c154053f29e420230de538db898" || error("frozen stream helper changed")
include(joinpath(QSOURCE,"gamma_functions.jl"));const QGP=EquilibriumGammaPilot
include("stream_spectral.jl")
function response_metrics(X,Y,g)
    size(X)==size(Y) || throw(DimensionMismatch("response comparator"))
    all(isfinite,X) && all(isfinite,Y) || error("nonfinite comparison input")
    ne=size(X,1);nk=ndims(X)==1 ? 1 : size(X,2)
    xr=reshape(X,ne,nk,:);yr=reshape(Y,ne,nk,:)
    d2=n2=dw2=nw2=maxabs=maxrel=maxscale=0.0;worst=(0,0);worstrel=(0,0)
    @inbounds for m in 1:nk,e in 1:ne
        bd=bn=0.0
        for a in axes(xr,3);bd+=abs2(yr[e,m,a]-xr[e,m,a]);bn+=abs2(xr[e,m,a]);end
        w=g.wᴱ[e]*(ndims(X)==1 ? 1.0 : g.wᵏ[m]);d2+=bd;n2+=bn;dw2+=w*bd;nw2+=w*bn
        da=sqrt(bd);scale=sqrt(bn);rel=QGP.relative(da,scale)
        if da>maxabs;maxabs=da;worst=(e,m);end
        if rel>maxrel;maxrel=rel;worstrel=(e,m);end
        maxscale=max(maxscale,scale)
    end
    r=QGP.relative(sqrt(d2),sqrt(n2));rw=QGP.relative(sqrt(dw2),sqrt(nw2));rb=QGP.relative(maxabs,maxscale)
    return Dict("absolute"=>sqrt(d2),"scale"=>sqrt(n2),"relative"=>r,"weighted_absolute"=>sqrt(dw2),"weighted_scale"=>sqrt(nw2),"weighted_relative"=>rw,
        "maximum_block_absolute"=>maxabs,"maximum_reference_block_scale"=>maxscale,"maximum_block_scaled"=>rb,"maximum_block_relative"=>maxrel,
        "maximum_absolute_energy_index"=>worst[1],"maximum_absolute_k_index"=>worst[2],"maximum_relative_energy_index"=>worstrel[1],"maximum_relative_k_index"=>worstrel[2],
        "bitwise_equal"=>reinterpret(UInt64,vec(X))==reinterpret(UInt64,vec(Y)),
        "original_arithmetic_pass"=>max(r,rw,rb)<=1e-10,"planned_tight_pass"=>max(r,rw,rb)<=1e-12,
        "weights"=>(ndims(X)==1 ? "wE diagnostic (q already contains wE,wk,spin/2pi)" : "wE*wk diagnostic"))
end
function qualifier_gate(cfg)
    cfg["runtime_authorized"]===true || error("write-only qualifier")
    cfg["memory_max_bytes"]==13019119616 && cfg["memory_high_bytes"]==12348030976 && cfg["swap_max_bytes"]==0 || error("source numerical/OS cap identity")
    cfg["workers"]==4 && cfg["BLAS_threads"]==1 && cfg["responses_max"]==2 && cfg["Phi_max"]==0 || error("operator call policy")
    Threads.nthreads()==4 || error("Julia4 required");BLAS.set_num_threads(1)
    cfg["external_seconds"]==300 && 0<cfg["soft_seconds"]<=290 || error("wall cap")
    cfg["single_file_max_bytes"]==cfg["data_max_bytes"]==32*1024^2 || error("output cap")
    root=dirname(dirname(pathof(QCLNEGF)))
    readchomp(`git -C $root rev-parse HEAD`)==cfg["core_head"] || error("Core HEAD")
    readchomp(`git -C $root rev-parse $("HEAD^{tree}")`)==cfg["core_tree"] || error("Core tree")
    isempty(readchomp(`git -C $root status --porcelain --untracked-files=no`)) || error("dirty Core")
    entry=only(filter(x->startswith(x,"0::"),readlines("/proc/self/cgroup")))
    cg=joinpath("/sys/fs/cgroup",lstrip(split(entry,':';limit=3)[3],'/'))
    for (file,key) in (("memory.max","memory_max_bytes"),("memory.high","memory_high_bytes"),("memory.swap.max","swap_max_bytes"))
        strip(read(joinpath(cg,file),String))==string(cfg[key]) || error("external cap "*file)
    end
    cpu=split(read(joinpath(cg,"cpu.max"),String));cpu[1]!="max" && parse(Int,cpu[1])<=4*parse(Int,cpu[2]) || error("CPU exceeds400%")
    limits=split(only(filter(x->startswith(x,"Max file size"),readlines("/proc/self/limits"))))
    limits[4]!="unlimited" && parse(Int,limits[4])<=cfg["single_file_max_bytes"] || error("external file cap")
    return Dict("cgroup"=>cg,"Core_import"=>pathof(QCLNEGF),"Core_HEAD"=>cfg["core_head"],"Julia"=>string(VERSION),"Julia_threads"=>Threads.nthreads(),"BLAS_threads"=>BLAS.get_num_threads())
end
function run_qualifier(input,out)
    cfg=TOML.parsefile(input);identity=qualifier_gate(cfg)
    ispath(out) && error("fresh output required");mkpath(out);begun=time();responses=Ref(0)
    report=Dict{String,Any}("identity"=>identity,"execution"=>"running","Phi_calls"=>0,"mu_fits"=>0,"mixing_steps"=>0,"scientific_acceptance"=>false,"input_sha256"=>qsha(input))
    function checkpoint(stage)
        report["stage"]=stage;report["responses_started"]=responses[];report["seconds"]=time()-begun;report["maxrss_bytes"]=Sys.maxrss()
        open(joinpath(out,"qualifier-summary.toml"),"w") do f;TOML.print(f,report);end
        QGP.data_bytes(out)<=cfg["data_max_bytes"] || error("data cap")
        time()-begun<=cfg["soft_seconds"] || error("soft wall cap")
        println("stage=",stage," responses=",responses[]);flush(stdout)
    end
    try
        for name in ("qualifier.jl","stream_spectral.jl","qualifier-protocol.md");cp(joinpath(@__DIR__,name),joinpath(out,name));end
        cp(input,joinpath(out,"input.toml"))
        hashes=Dict{String,String}()
        for (name,expected) in cfg["source_sha256"]
            length(expected)==64 && expected!="UNSET" || error("ROOT must bind whole source SHA: "*name)
            hashes[name]=qsha(joinpath(QSOURCE,name));hashes[name]==expected || error("source whole SHA: "*name)
        end
        Set(keys(hashes))==Set(("Gamma-last.h5","model.h5","input.toml","identity.toml","map-064.toml","summary.toml")) || error("source set")
        report["source_sha256"]=hashes;report["source_definition_sha256"]=QDEFS;report["hash_seconds"]=time()-begun
        report["code_sha256"]=Dict(name=>qsha(joinpath(@__DIR__,name)) for name in ("qualifier.jl","stream_spectral.jl","qualifier-protocol.md"))
        source=TOML.parsefile(joinpath(QSOURCE,"input.toml"));sid=TOML.parsefile(joinpath(QSOURCE,"identity.toml"));completed=TOML.parsefile(joinpath(QSOURCE,"summary.toml"))
        completed["execution"]=="complete" && completed["pi_calls"]==64 || error("source completion")
        source["memory_max_bytes"]==cfg["memory_max_bytes"] && source["workers"]==4 || error("source numeric budget/workers")
        sid["source_sha256"]["stable_spectral.jl"]==QDEFS["stable_spectral.jl"] && sid["input_sha256"]==hashes["input.toml"] || error("source identity")
        string(VERSION)==source["expected_julia_version"] || error("source Julia version")
        report["manifest_sha256"]=qsha(joinpath(dirname(Base.active_project()),"Manifest.toml"))
        report["manifest_sha256"]==sid["manifest_sha256"] || error("source manifest")
        checkpoint("all_sources_hashed")
        c=QGP.build_context(source)
        c.options.algorithms.embedding===:finite_chain && c.options.algorithms.embedding_periods==4 &&
            iszero(QGP.N._scaled_physics(c.problem.physical,c.problem.scales).Eᵖ) || error("V0 finite_chain4 qualifier only")
        TOML.parsefile(joinpath(QSOURCE,"map-064.toml"))["call"]==64 || error("source map locator")
        h5open(joinpath(QSOURCE,"model.h5"),"r") do f;report["model_binding"]=QGP.verify_problem_model(c,f);end
        G=h5open(joinpath(QSOURCE,"Gamma-last.h5"),"r") do f
            read(f["schema"])=="campaign133-gamma-checkpoint-v1" && read(f["pi_calls_total"])==64 && read(f["identity"])==source["identity"] || error("checkpoint metadata")
            QGP.read_complex(f,"Gamma")
        end
        size(G)==(30704,11,5,5) || error("source Gamma dimensions")
        report["Gamma_array_sha256"]=QGP.array_digest(G);report["input_cone"]=Dict(string(k)=>v for (k,v) in pairs(QGP.cone(G))) # no strict-PD fiction
        F=sizeof(ComplexF64)*length(G)
        report["memory_array_ledger"]=Dict("full_field_bytes"=>F,"old_conservative_spectral_peak_21F"=>21F,
            "stream_own_sweep_4F"=>4F,"retained_old_GR_A_2F"=>2F,"paired_stream_sweep_6F"=>6F,
            "retained_old_small_arrays_bytes"=>2*sizeof(Float64)*size(G,1)*size(G,2)+sizeof(Float64)*size(G,1),
            "cache_model_plan_PI_allocator_RSS_extra"=>"not included; real OS cap retained","native_full_guard_unchanged"=>true,
            "native_full_estimate_bytes"=>c.cache.estimate.peak_bytes,"numeric_budget_bytes"=>c.options.memory_budget_bytes)
        responses[]+=1;checkpoint("old_response_started");t=time()
        old=QGP.stable_spectral_state(c,G;output=out,cfg=cfg,call=0)
        report["old_seconds_observed"]=time()-t;report["old_response_completed"]=true
        oldGR,oldA,oldcondition,oldscale,oldroundoff,oldq,oldledger=old;old=nothing;GC.gc(true)
        selectedpath=joinpath(out,"spectral-selected-000.h5")
        chosen=h5open(selectedpath,"r") do f;Int.(read(f["selected_energy_indices"]));end
        1<=length(chosen)<=8 || error("selected response cap")
        responses[]<2 || error("response cap");responses[]+=1;checkpoint("stream_response_started");t=time()
        fresh=StreamingSpectral.stream_spectrum(QGP,c,G;selected_indices=chosen)
        report["stream_seconds_observed"]=time()-t;report["stream_response_completed"]=true
        report["timing_scope"]="one first-call pair; old includes full diagnostics/selected H5 I/O; stream includes reduced diagnostics; JIT/order effects; no speed ratio or warmed benchmark claim"
        comparisons=Dict{String,Any}()
        for (name,X,Y) in (("GR",oldGR,fresh.GR),("A",oldA,fresh.A),("q",oldq,fresh.q),("condition",oldcondition,fresh.condition),("scale",oldscale,fresh.scale))
            comparisons[name]=response_metrics(X,Y,c.problem.grids)
        end
        report["all_entry_comparisons"]=comparisons;report["old_PI_roundoff"]=oldroundoff;report["stream_PI_roundoff"]=fresh.roundoff
        report["old_ledger"]=oldledger;report["stream_ledger"]=fresh.ledger
        selected=Dict{String,Any}()
        h5open(selectedpath,"r") do f
            for (key,name) in ((:Gamma,"Gamma"),(:Lambda_raw,"Lambda_raw"),(:Lambda_H_preassembly,"Lambda_H_preassembly"),(:GR,"GR"),(:A_cov,"A_cov"),(:A_sub,"A_sub"),(:localR,"localR"),(:totalR,"totalR"),(:Gamma_physical_total,"Gamma_physical_total"),(:Gamma_eff,"Gamma_eff"))
                X=QGP.read_complex(f,name);Y=fresh.samples[key]
                metric=response_metrics(X,Y,(wᴱ=c.problem.grids.wᴱ[chosen],wᵏ=c.problem.grids.wᵏ))
                selected[name]=metric
            end
        end
        report["selected_comparisons"]=selected;report["selected_energy_indices"]=chosen
        QGP.atomic_h5(joinpath(out,"stream-selected.h5"),cfg;expected_bytes=2*1024^2) do f
            f["selected_energy_indices"]=chosen;f["energy"]=c.problem.grids.ε[chosen];f["weight_k"]=c.problem.grids.wᵏ
            for name in (:Gamma,:Lambda_raw,:Lambda_H_preassembly,:GR,:A_cov);QGP.write_complex(f,string(name),fresh.samples[name]);end
        end
        report["original_arithmetic_pass"]=all(v["original_arithmetic_pass"] for v in values(comparisons)) && all(v["original_arithmetic_pass"] for v in values(selected))
        report["planned_tight_pass"]=all(v["planned_tight_pass"] for v in values(comparisons)) && all(v["planned_tight_pass"] for v in values(selected))
        report["Gamma_input_unchanged_after"]=QGP.array_digest(G)==report["Gamma_array_sha256"]
        report["Gamma_input_unchanged_after"] || error("borrowed Gamma mutated")
        report["sources_unchanged_after"]=all(qsha(joinpath(QSOURCE,name))==sha for (name,sha) in hashes)
        report["sources_unchanged_after"] || error("source changed during qualifier")
        report["definitions_unchanged_after"]=all(qsha(joinpath(QSOURCE,name))==sha for (name,sha) in QDEFS) &&
            qsha(joinpath(@__DIR__,"stream_spectral.jl"))=="8df7885cbf547ffd4077e611f6a88b933dc75c154053f29e420230de538db898"
        report["definitions_unchanged_after"] || error("included code changed during qualifier")
        report["original_arithmetic_pass"] || error("original arithmetic equivalence target failed")
        report["execution"]="complete_operator_equivalence_only";checkpoint("complete")
    catch err
        report["execution"]="stopped";report["reason"]=sprint(showerror,err);rethrow()
    finally
        report["seconds"]=time()-begun;report["maxrss_bytes"]=Sys.maxrss();report["responses_started"]=responses[]
        open(joinpath(out,"qualifier-summary.toml"),"w") do f;TOML.print(f,report);end
    end
end
function tiny_qualifier_checks()
    g=(wᴱ=[.2,.7],wᵏ=[.3,.8]);X=reshape(ComplexF64[complex(i,.1*i) for i in 1:8],2,2,2,1);Y=(1+.01im)*X
    @testset "all-entry comparator independent norm/weight and near-zero convention" begin
        r=response_metrics(X,Y,g)
        @test r["absolute"]≈norm(Y-X)
        @test r["weighted_absolute"]≈sqrt(sum(g.wᴱ[e]*g.wᵏ[m]*abs2(Y[e,m,a,1]-X[e,m,a,1]) for a in 1:2,m in 1:2,e in 1:2))
        @test response_metrics(X,X,g)["bitwise_equal"]
        @test response_metrics(zeros(2),zeros(2),g)["relative"]==0
        @test response_metrics(zeros(2),ones(2),g)["relative"]==Inf
        @test_throws ErrorException response_metrics([NaN,0.0],zeros(2),g)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS==["--tiny-qualifier"];tiny_qualifier_checks()
    else
        length(ARGS)==2 || error("usage: qualifier.jl --tiny-qualifier | INPUT NEW_OUTPUT")
        run_qualifier(ARGS[1],abspath(ARGS[2]))
    end
end
