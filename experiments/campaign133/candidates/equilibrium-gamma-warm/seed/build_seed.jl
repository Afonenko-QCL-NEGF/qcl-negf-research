using QCLNEGF, HDF5, LinearAlgebra, TOML, SHA, Dates
include("gamma_functions.jl")
include("warm_checks.jl")
const W=EquilibriumGammaPilot

function builder_gate(cfg)
    cfg["runtime_authorized"]===true || error("Builder implementation-only freeze; runtime not authorized")
    string(VERSION)==cfg["expected_julia_version"] && Threads.nthreads()==1 || error("builder Julia/thread identity differs")
    cfg["source_state_sha256"]==W.WARM_SOURCE_SHA && cfg["target_model_sha256"]==W.WARM_MODEL_SHA || error("builder source contract differs")
    0<cfg["memory_max_bytes"]<=4*1024^3 && 0<cfg["memory_high_bytes"]<=min(3*1024^3,cfg["memory_max_bytes"]) || error("builder memory envelope exceeds4/3GiB")
    0<cfg["builder_seconds"]<=300 && 1<=cfg["chunk_energies"]<=1024 || error("builder time/chunk envelope differs")
    cfg["single_file_max_bytes"]<=512*1024^2 && cfg["data_max_bytes"]<=512*1024^2 || error("builder output envelope exceeds512MiB")
    core=dirname(dirname(pathof(QCLNEGF)))
    head=readchomp(`git -C $core rev-parse HEAD`)
    tree=readchomp(`git -C $core rev-parse $("HEAD^{tree}")`)
    head==cfg["expected_core_head"] && tree==cfg["expected_core_tree"] || error("builder imported source differs")
    isempty(readchomp(`git -C $core status --porcelain --untracked-files=no`)) || error("builder Core tracked source is dirty")
    membership=only(filter(line->startswith(line,"0::"),readlines("/proc/self/cgroup")))
    cgroup=joinpath("/sys/fs/cgroup",lstrip(split(membership,':';limit=3)[3],'/'))
    for (field,key) in (("memory.max","memory_max_bytes"),("memory.high","memory_high_bytes"),("memory.swap.max","swap_max_bytes"))
        strip(read(joinpath(cgroup,field),String))==string(cfg[key]) || error("builder cgroup $field differs")
    end
    cpu=split(strip(read(joinpath(cgroup,"cpu.max"),String)))
    cpu[1]!="max" && parse(Int,cpu[1])/parse(Int,cpu[2])<=1 || error("builder CPU exceeds100%")
    limit=split(only(filter(line->startswith(line,"Max file size"),readlines("/proc/self/limits"))))
    limit[4]!="unlimited" && parse(Int,limit[4])<=cfg["single_file_max_bytes"] || error("builder needs LimitFSIZE<=512MiB")
    return Dict("core_head"=>head,"core_tree"=>tree,"core_source"=>pathof(QCLNEGF),"cgroup"=>cgroup)
end

function build_seed(input,output)
    cfg=TOML.parsefile(input);cfg["_input_directory"]=dirname(abspath(input))
    report=Dict{String,Any}("source_identity"=>builder_gate(cfg))
    ispath(output) && error("Refusing to overwrite builder output")
    mkpath(output);cp(input,joinpath(output,"input.toml"))
    files=("build_seed.jl","gamma_functions.jl","warm_seed.jl","warm_checks.jl","protocol.md","algorithm.md")
    for name in files;cp(joinpath(@__DIR__,name),joinpath(output,name));end
    started=time();checkwall()=time()-started<cfg["builder_seconds"] || error("builder wall cap exhausted")
    BLAS.set_num_threads(1)
    try
        warm_small_checks(cfg)
        report["small_checks"]="pass"
        for (pathkey,shakey) in (("source_state_path","source_state_sha256"),
                              ("source_input_path","source_input_sha256"),
                              ("target_model_path","target_model_sha256"),
                              ("target_input_path","target_input_sha256"))
            checkwall();W.checked_digest(W.warm_path(cfg,pathkey),cfg[shakey])
        end
        W.verify_model_inputs(cfg)
        sourcepath,targetpath=W.warm_path(cfg,"source_state_path"),W.warm_path(cfg,"target_model_path")
        fingerprints=HDF5.h5open(sourcepath,"r") do source
            read(source["schema"])=="campaign133-inner-raw-v3" || error("source schema differs")
            HDF5.h5open(targetpath,"r") do target;W.verify_source_model(source,target);end
        end
        report["source_target_fingerprints"]=fingerprints
        for (key,name) in (("source_input_path","source-input.toml"),("source_identity_path","source-identity.toml"),
                           ("target_input_path","target-input.toml"),("target_identity_path","target-identity.toml"))
            cp(W.warm_path(cfg,key),joinpath(output,name))
        end
        W.build_source_gamma(sourcepath,output,cfg,checkwall,report)
        W.build_target_seed(targetpath,output,cfg,checkwall,report)
        checkwall()
        manifest=Dict{String,Any}("schema"=>"campaign133-warm-seed-manifest-v1","execution"=>"complete",
            "identity"=>cfg["identity"],"source_state_path"=>sourcepath,"source_state_sha256"=>cfg["source_state_sha256"],
            "source_state_bytes"=>filesize(sourcepath),"target_model_path"=>targetpath,
            "source_input_sha256"=>cfg["source_input_sha256"],"target_input_sha256"=>cfg["target_input_sha256"],
            "target_model_sha256"=>cfg["target_model_sha256"],"seed_file"=>"gamma-seed.h5",
            "seed_sha256"=>W.file_digest(joinpath(output,"gamma-seed.h5")),
            "source_gamma_sha256"=>W.file_digest(joinpath(output,"source-gamma.h5")),
            "builder_input_sha256"=>W.file_digest(input),"core_head"=>cfg["expected_core_head"],
            "core_tree"=>cfg["expected_core_tree"],"builder_source_sha256"=>Dict(name=>W.file_digest(joinpath(@__DIR__,name)) for name in files),
            "metric_definition"=>"source GammaK-GammaR: scale GammaK, weighted source wE*wk; representation: raw-Herm with scale raw and target wE*wk; zero/zero0, nonzero/zeroInf; no eigen-clipping",
            "report"=>report,"scope"=>"Gamma-only warm initialization; raw and explicit Hermitian representation both retained; no eigen-clipping or source occupations imported")
        open(joinpath(output,"seed-manifest.toml"),"w") do io;TOML.print(io,manifest);end
        report["execution"]="complete"
        println("BuilderComplete seedSHA=",manifest["seed_sha256"]," manifestSHA=",W.file_digest(joinpath(output,"seed-manifest.toml")))
    catch error
        report["execution"]="stopped";report["reason"]=sprint(showerror,error,catch_backtrace());rethrow()
    finally
        report["elapsed_seconds"]=time()-started;report["data_bytes"]=W.data_bytes(output)
        open(joinpath(output,"builder-summary.toml"),"w") do io;TOML.print(io,report);end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("usage: build_seed.jl BUILDER_INPUT.toml NEW_OUTPUT_DIRECTORY")
    build_seed(abspath(ARGS[1]),abspath(ARGS[2]))
end
