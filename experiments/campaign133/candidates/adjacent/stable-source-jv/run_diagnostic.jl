# Finite adapter only. No Newton update, native solve, Pkg operation or run on include.
const PROCESS_STARTED=time()
using TOML, SHA, LinearAlgebra, QCLNEGF, HDF5, Dates
const PROTOTYPE_SHA="5a56b333d40a71794de69d1d95483e0ef5ff532d5f3e4506b8e78aed3ac22615"
const SERIALIZER_SHA="6157cdb7dba7c923f1e06ea73bc43d9010f50d6b59dc9a387b813c0de682eea4"
const STREAM_SHA="8df7885cbf547ffd4077e611f6a88b933dc75c154053f29e420230de538db898"
const PARENT_SHA="50bcfaaeb3c399baaaa2a4b9b9e023167e97b4d87c6af95d83353cb453072211"
const SPECTRAL_SHA="eeca1f29710d9bb60bf6ce70ce0ebd1620fb86d4c77314b4bd8ffa9a5bcd5aa5"
digest(path)=bytes2hex(open(sha256,path))
checksha(path,expected)=(digest(path)==expected || error("SHA mismatch: $path"))
resolve(cfg,key)=normpath(joinpath(@__DIR__,cfg[key])) # frozen frame, independent of authorized-copy location
include("serialization.jl")
using .DiagnosticSerialization:plain,small_toml
function time_gate(cfg)
    time()-PROCESS_STARTED<cfg["wall_seconds"] || error("whole-process600s cap")
end
function source_gate(cfg;tiny=false)
    frozen=TOML.parsefile(joinpath(@__DIR__,"adapter-input.toml"))
    for (k,v) in frozen
        k=="runtime_authorized" && continue
        haskey(cfg,k) && isequal(cfg[k],v) || error("frozen option changed: $k")
    end
    tiny || cfg["runtime_authorized"]===true || error("runtime_authorized=false")
    core=dirname(dirname(pathof(QCLNEGF)))
    readchomp(`git -C $core rev-parse HEAD`)==cfg["expected_core_head"] || error("Core HEAD")
    readchomp(`git -C $core rev-parse $("HEAD^{tree}")`)==cfg["expected_core_tree"] || error("Core tree")
    isempty(readchomp(`git -C $core status --porcelain --untracked-files=no`)) || error("dirty Core")
    string(VERSION)==cfg["expected_julia_version"] || error("Julia version")
    project=Base.active_project();project===nothing && error("no project")
    manifest=joinpath(dirname(project),"Manifest.toml");checksha(manifest,cfg["manifest_sha256"])
    parent=resolve(cfg,"source_directory");stream=resolve(cfg,"streaming_helper")
    checksha(joinpath(@__DIR__,"analytic_jv.jl"),PROTOTYPE_SHA);checksha(stream,STREAM_SHA)
    checksha(joinpath(@__DIR__,"serialization.jl"),SERIALIZER_SHA)
    checksha(joinpath(parent,"gamma_functions.jl"),PARENT_SHA)
    checksha(joinpath(parent,"stable_spectral.jl"),SPECTRAL_SHA)
    checksha(joinpath(parent,"identity.toml"),cfg["source_identity_sha256"])
    identity=TOML.parsefile(joinpath(parent,"identity.toml"))
    identity["core_head"]==cfg["expected_core_head"] && identity["core_tree"]==cfg["expected_core_tree"] || error("parent Core")
    identity["manifest_sha256"]==cfg["manifest_sha256"] || error("parent manifest")
    for (name,hash) in identity["source_sha256"];checksha(joinpath(parent,name),hash);end
    # Verify this adapter's complete text frame before model arrays.
    workspace=normpath(joinpath(@__DIR__,"../../../.."))
    for line in eachline(joinpath(@__DIR__,"adapter-SHA256SUMS"))
        hash,path=split(line;limit=2);checksha(joinpath(workspace,strip(path)),hash)
    end
    if !tiny
        Threads.nthreads()==cfg["workers"] || error("J4 required")
        member=only(filter(s->startswith(s,"0::"),readlines("/proc/self/cgroup")))
        group=joinpath("/sys/fs/cgroup",lstrip(split(member,':';limit=3)[3],'/'))
        for (file,key) in (("memory.max","memory_max_bytes"),("memory.high","memory_high_bytes"),("memory.swap.max","swap_max_bytes"))
            strip(read(joinpath(group,file),String))==string(cfg[key]) || error("cgroup $file")
        end
        cpu=split(strip(read(joinpath(group,"cpu.max"),String)))
        cpu[1]!="max" && parse(Int,cpu[1])/parse(Int,cpu[2])<=4 || error("CPU cap")
        limit=split(only(filter(s->startswith(s,"Max file size"),readlines("/proc/self/limits"))))
        limit[4]!="unlimited" && parse(Int,limit[4])<=cfg["single_file_max_bytes"] || error("LimitFSIZE2GiB")
    end
    return parent,stream,project,manifest,identity
end
function import_frozen(parent,stream)
    Base.include(Main,joinpath(parent,"gamma_functions.jl"))
    Base.include(Main,stream)
    Base.include(Main,joinpath(@__DIR__,"analytic_jv.jl"))
end
function run_actual(input,output)
    cfg=TOML.parsefile(input);cfg["_input_directory"]=dirname(abspath(input))
    parent,stream,project,manifest,parentid=source_gate(cfg)
    if !isdefined(Main,:StableSourceJv)
        import_frozen(parent,stream)
        return Base.invokelatest(run_actual,input,output) # newly included methods, no old-world calls
    end
    for (name,key) in (("Gamma-last.h5","source_Gamma_sha256"),("model.h5","source_model_sha256"),
        ("reduced-terminal.h5","source_terminal_sha256"),("input.toml","source_input_sha256"),
        ("map-064.toml","source_map_sha256"),("summary.toml","source_summary_sha256"))
        checksha(joinpath(parent,name),cfg[key]);time_gate(cfg)
    end
    pcfg=TOML.parsefile(joinpath(parent,"input.toml"));summary=TOML.parsefile(joinpath(parent,"summary.toml"))
    pcfg["project_AA_trial"]===false && summary["execution"]=="complete" && summary["pi_calls"]==64 || error("parent state")
    parentid["input_sha256"]==cfg["source_input_sha256"] || error("parent input identity")
    HDF5.h5open(joinpath(parent,"Gamma-last.h5"),"r") do f
        read(f["schema"])=="campaign133-gamma-checkpoint-v1" && read(f["identity"])==pcfg["identity"] && read(f["pi_calls_total"])==64 || error("checkpoint identity")
        expected=(pcfg["NE"],pcfg["Nk"],pcfg["Nb"],pcfg["Nb"])
        size(f["Gamma_real"])==size(f["Gamma_imag"])==expected || error("checkpoint axes")
    end
    HDF5.h5open(joinpath(parent,"model.h5"),"r") do f
        read(f["schema"])=="campaign133-gamma-model-v2" || error("model schema")
    end
    ispath(output) && error("output exists");mkpath(output)
    # Archive first, then hash the archived sources for identity.
    for (dir,names) in (("prototype",readdir(@__DIR__)),("parent",collect(keys(parentid["source_sha256"]))))
        dest=joinpath(output,"source",dir);mkpath(dest)
        for name in names
            src=dir=="prototype" ? joinpath(@__DIR__,name) : joinpath(parent,name)
            isfile(src) && cp(src,joinpath(dest,name))
        end
    end
    cp(stream,joinpath(output,"source","stream_spectral.jl"));cp(input,joinpath(output,"input.toml"))
    cp(project,joinpath(output,"environment-Project.toml"));cp(manifest,joinpath(output,"environment-Manifest.toml"))
    sources=Dict(relpath(joinpath(r,f),output)=>digest(joinpath(r,f)) for (r,_,fs) in walkdir(joinpath(output,"source")) for f in fs)
    small_toml(joinpath(output,"identity.toml"),(input_sha256=digest(input),sources=sources,
        source_Gamma_sha256=cfg["source_Gamma_sha256"],manifest_sha256=digest(manifest),
        core_head=cfg["expected_core_head"],core_tree=cfg["expected_core_tree"],start_utc=string(Dates.now(Dates.UTC)),
        method="analytic covariance Jv; diagnostic only; source exact-PD proof pending"))
    GP=Main.EquilibriumGammaPilot;SS=Main.StreamingSpectral;J=Main.StableSourceJv
    BLAS.set_num_threads(1);budget=J.Limits(PROCESS_STARTED,cfg["wall_seconds"],0,0,0,0)
    report=Dict{String,Any}("scientific_acceptance"=>false,"native_calls"=>0,"Gamma_updates"=>0)
    io_cfg=copy(cfg)
    function request(label;phi_more=0,jv_more=0)
        time_gate(cfg)
        small_toml(joinpath(output,"request-"*label*".toml"),(Phi_entered_lower=budget.pi,
            Phi_entered_upper=budget.pi+phi_more,Jv_entered_lower=budget.jv,
            Jv_entered_upper=budget.jv+jv_more,Phi_completed=budget.pi_completed,
            Jv_completed=budget.jv_completed,elapsed_seconds=time()-PROCESS_STARTED,
            meaning="pre-call structural bounds; hard-kill need not execute finally; not a completed map claim"))
    end
    save_raw(name,X)=GP.atomic_h5(joinpath(output,name*".h5"),io_cfg;expected_bytes=GP.expected_padded_complex_bytes(size(X))+1024^2) do f
        f["schema"]="campaign133-analytic-Jv-raw-v1";GP.write_complex(f,name,X)
    end
    function preserver(name,metadata)
        return function(raw,stats)
            save_raw(name,raw)
            small_toml(joinpath(output,metadata*".toml"),merge(stats,
                (Phi_entered=budget.pi,Phi_completed=budget.pi_completed,
                 completed_scope="raw contraction produced; cone/return/extra-PI not yet admitted")))
        end
    end
    g=nothing
    try
        time_gate(cfg);c=GP.build_context(pcfg)
        HDF5.h5open(joinpath(parent,"model.h5"),"r") do f;GP.verify_problem_model(c,f);end
        gamma=HDF5.h5open(joinpath(parent,"Gamma-last.h5"),"r") do f
            read(f["schema"])=="campaign133-gamma-checkpoint-v1" && read(f["identity"])==pcfg["identity"] && read(f["pi_calls_total"])==64 || error("checkpoint identity")
            expected=(pcfg["NE"],pcfg["Nk"],pcfg["Nb"],pcfg["Nb"])
            size(f["Gamma_real"])==size(f["Gamma_imag"])==expected || error("checkpoint axes")
            GP.read_complex(f,"Gamma")
        end
        g=c.problem.grids;time_gate(cfg)
        request("base";phi_more=1)
        base=J.fresh_phi(GP,SS,c,gamma,budget;keep_base=true,
            persist_raw=preserver("base_raw_Phi","base-raw-produced"))
        small_toml(joinpath(output,"base.toml"),merge(base.stats,
            (maximum_Dyson_condition=maximum(base.condition),minimum_Dyson_scale=minimum(base.scale))))
        replaycfg=Dict{String,Any}("_input_directory"=>dirname(parent),"parent_directory"=>parent,"parent_call"=>64)
        replaystats=Dict("residual"=>plain(GP.norms(base.raw,gamma,g.wᴱ,g.wᵏ)),"mu"=>plain(base.stats.mu))
        GP.restart_boundary_audit(c,base.raw,replaystats,replaycfg,output)
        pairs=Array{ComplexF64,4}[]
        for (i,eps) in enumerate(cfg["probe_relative_amplitudes"])
            request("probe$(i)-plus";phi_more=1)
            plus,_=J.fresh_phi(GP,SS,c,(1+eps).*gamma,budget;
                persist_raw=preserver("probe$(i)_plus","probe$(i)-plus"))
            request("probe$(i)-minus";phi_more=1)
            minus,_=J.fresh_phi(GP,SS,c,(1-eps).*gamma,budget;
                persist_raw=preserver("probe$(i)_minus","probe$(i)-minus"))
            @. plus=(plus-minus)/(2eps)
            save_raw("central$(i)",plus);push!(pairs,plus)
        end
        request("validation-Jv";jv_more=1);analytic,meta=J.jvp(GP,c,base,gamma,budget)
        save_raw("analytic_Gamma_direction",analytic);small_toml(joinpath(output,"validation-Jv.toml"),meta)
        comparisons=[GP.norms(analytic,D,g.wᴱ,g.wᵏ) for D in pairs]
        diff=GP.norms(pairs[1],pairs[2],g.wᴱ,g.wᵏ)
        arith=cfg["arithmetic_relative"]
        qualifies(r)=isfinite(r.spectral.primal_GR_relative) && r.spectral.primal_GR_relative<=arith &&
            isfinite(r.spectral.primal_A_relative) && r.spectral.primal_A_relative<=arith &&
            isfinite(r.spectral.PI_roundoff) && r.spectral.PI_roundoff<=arith &&
            isfinite(r.number.dN_relative) && r.number.dN_relative<=arith &&
            isfinite(r.representation.relative) && r.representation.relative<=arith &&
            isfinite(r.representation.weighted_relative) && r.representation.weighted_relative<=arith
        # Empirical Richardson consistency only: leading central errors have ratio4.
        fine,coarse=comparisons[2],comparisons[1]
        fdok=(fine.absolute<=diff.absolute || fine.relative<=arith) &&
             (fine.weighted_absolute<=diff.weighted_absolute || fine.weighted_relative<=arith) &&
             (fine.absolute<=coarse.absolute || coarse.relative<=arith)
        small_toml(joinpath(output,"finite-difference-validation.toml"),(comparisons=comparisons,
            central_difference=diff,empirical_consistency=fdok,
            qualification="empirical epsilon-halving consistency, not rigorous derivative-error enclosure"))
        qualifies(meta) && fdok || error("derivative evidence not qualified; no Arnoldi")
        request("Arnoldi3";jv_more=3)
        function qualify_and_save(r)
            small_toml(joinpath(output,"krylov-Jv-$(budget.jv).toml"),merge(r,
                (Jv_entered=budget.jv,Jv_completed=budget.jv_completed,)))
            return qualifies(r)
        end
        k=J.arnoldi3(GP,c,base,budget;qualify=qualify_and_save)
        small_toml(joinpath(output,"krylov.toml"),(Hbar=k.Hbar,Ritz_real=real.(k.Ritz),Ritz_imag=imag.(k.Ritz),
            Ritz_residual=k.Ritz_residual,restricted_B_singular_values=k.restricted_B_singular_values,
            near_null_coordinates=k.near_null_coordinates,gmres_coordinates=k.gmres_coordinates,
            gmres_predicted_residual=k.gmres_predicted_residual,seed_conversion=k.seed_conversion,
            seed_norm=k.seed_norm,records=k.records,qualification=k.qualification))
        for (i,Q) in enumerate(k.Q);save_raw("Arnoldi_Q$(i)",Q);end
        report["execution"]="complete"
    catch err
        report["execution"]="failed";report["reason"]=sprint(showerror,err)
        rethrow()
    finally
        merge!(report,Dict("Phi_entered"=>budget.pi,"Phi_completed"=>budget.pi_completed,
            "Phi_completed_scope"=>"raw contraction produced, not cone/admission/extra-PI success",
            "Jv_entered"=>budget.jv,"Jv_completed"=>budget.jv_completed,
            "elapsed_seconds"=>time()-PROCESS_STARTED,"data_bytes"=>GP.data_bytes(output)))
        small_toml(joinpath(output,"summary.toml"),report)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS==["--tiny"]
        cfg=TOML.parsefile(joinpath(@__DIR__,"adapter-input.toml"));cfg["_input_directory"]=@__DIR__
        parent,stream,_,_,_=source_gate(cfg;tiny=true);import_frozen(parent,stream)
        BLAS.set_num_threads(1);include("serialization_checks.jl");serialization_checks()
        include("tiny_checks.jl");tiny_checks(Main.EquilibriumGammaPilot)
    else
        length(ARGS)==2 || error("usage: --tiny OR AUTHORIZED_INPUT NEW_OUTPUT")
        run_actual(abspath(ARGS[1]),abspath(ARGS[2]))
    end
end
