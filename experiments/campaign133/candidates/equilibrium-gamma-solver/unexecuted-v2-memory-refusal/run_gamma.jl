# Experimental campaign harness. No Pkg operations, no Runner dependency.
using QCLNEGF, TOML, SHA, LinearAlgebra, Dates
include("gamma_functions.jl")
include("checks.jl")

function source_gate(cfg)
    cfg["runtime_authorized"] === true || error("Implementation frozen for review; runtime_authorized=false")
    string(VERSION) == cfg["expected_julia_version"] || error("Julia version identity mismatch")
    cfg["expected_core_head"] != "UNSET" || error("Accepted Core source is not frozen")
    core = dirname(dirname(pathof(QCLNEGF)))
    head = readchomp(`git -C $core rev-parse HEAD`)
    tree = readchomp(`git -C $core rev-parse $("HEAD^{tree}")`)
    dirty = readchomp(`git -C $core status --porcelain --untracked-files=no`)
    head == cfg["expected_core_head"] || error("Imported Core HEAD mismatch: $head")
    tree == cfg["expected_core_tree"] || error("Imported Core tree mismatch: $tree")
    isempty(dirty) || error("Imported Core contains tracked changes")
    Threads.nthreads() == cfg["workers"] || error("Julia thread count mismatch")
    cfg["pi_calls_total"] == 16 || error("This frozen harness admits only pilot16")
    cfg["native_max_maps"] == 32 || error("Frozen native structural cap is1+30+1")
    cfg["memory_max_bytes"] == 13287555072 || error("OS/Core memory budget mismatch")
    # OS controls are imposed by Root's existing systemd-run invocation; read the
    # active cgroup rather than pretending ProductionOptions enforces RSS/CPU.
    membership = only(filter(line->startswith(line,"0::"),readlines("/proc/self/cgroup")))
    cgroup = joinpath("/sys/fs/cgroup",lstrip(split(membership,':';limit=3)[3],'/'))
    controls = Dict{String,Any}()
    for (name,key) in (("memory.max","memory_max_bytes"),("memory.high","memory_high_bytes"),("memory.swap.max","swap_max_bytes"))
        raw = strip(read(joinpath(cgroup,name),String))
        raw == string(cfg[key]) || error("cgroup $name=$raw, expected $(cfg[key])")
        controls[name] = raw
    end
    cpu = split(strip(read(joinpath(cgroup,"cpu.max"),String)))
    cpu[1] != "max" && parse(Int,cpu[1])/parse(Int,cpu[2]) <= 4 || error("CPU cgroup exceeds400%")
    controls["cpu.max"] = join(cpu," ")
    limits = readlines("/proc/self/limits")
    file_limit = split(only(filter(line->startswith(line,"Max file size"),limits)))
    file_limit[4] != "unlimited" && parse(Int,file_limit[4]) <= cfg["single_file_max_bytes"] ||
        error("Set external LimitFSIZE<=5GiB before runtime")
    controls["file_size_soft_limit"] = file_limit[4]
    return Dict("core_source"=>pathof(QCLNEGF),"core_head"=>head,"core_tree"=>tree,
        "core_dirty"=>dirty,"julia_version"=>string(VERSION),"cgroup"=>cgroup,"controls"=>controls)
end
function run_pilot(input,output)
    cfg = TOML.parsefile(input)
    identity = source_gate(cfg) # NO model allocations before explicit gate
    ispath(output) && error("Refusing to overwrite output $output")
    mkpath(output)
    for file in ("run_gamma.jl","gamma_functions.jl","checks.jl","protocol.md","algorithm.md")
        cp(joinpath(@__DIR__,file),joinpath(output,file))
    end
    cp(input,joinpath(output,"input.toml"))
    project = Base.active_project()
    project === nothing && error("A native frozen Julia project is required")
    manifest = joinpath(dirname(project),"Manifest.toml")
    isfile(manifest) || error("Native Pkg Manifest is required for provenance")
    cp(project,joinpath(output,"environment-Project.toml"))
    cp(manifest,joinpath(output,"environment-Manifest.toml"))
    identity["active_project"] = project
    identity["manifest_sha256"] = bytes2hex(sha256(read(manifest)))
    identity["start_utc"] = string(Dates.now(Dates.UTC))*"Z"
    identity["input_sha256"] = bytes2hex(sha256(read(input)))
    identity["source_sha256"] = Dict(file=>bytes2hex(sha256(read(joinpath(@__DIR__,file))))
        for file in ("run_gamma.jl","gamma_functions.jl","checks.jl","protocol.md","algorithm.md"))
    identity["method"] = "custom Gamma AAII m3; native operator policy linear; explicit PI/Hilbert; current common-mu FD restriction"
    open(joinpath(output,"identity.toml"),"w") do io; TOML.print(io,identity); end
    BLAS.set_num_threads(1)
    budget = EquilibriumGammaPilot.Budget(cfg)
    summary = Dict{String,Any}("scientific_acceptance"=>false)
    try
        small_checks(cfg)
        summary["small_checks"] = "pass"
        context = EquilibriumGammaPilot.build_context(cfg)
        EquilibriumGammaPilot.check_time(budget)
        EquilibriumGammaPilot.save_model(context,output,cfg)
        gamma,trigger,reduced = EquilibriumGammaPilot.reduced_stage(context,budget,output,cfg)
        summary["reduced"] = reduced
        # Hard lifetime barrier: reduced_stage returned ONLY Gamma and scalars;
        # histories, QR columns and its old Green arrays are unreachable.
        GC.gc(true)
        if trigger && !EquilibriumGammaPilot.stage_expired(budget)
            summary["native"] = EquilibriumGammaPilot.native_stage(context,gamma,budget,output,cfg)
        else
            summary["native"] = Dict("not_measured"=>true,"reason"=>trigger ? "wall budget" : "no reduced trigger")
        end
        summary["execution"] = "complete"
    catch error
        summary["execution"] = "stopped"
        summary["reason"] = sprint(showerror,error,catch_backtrace())
        rethrow()
    finally
        summary["pi_calls"] = budget.pi
        summary["precheck_actual_calls"] = budget.precheck
        summary["native_accounting"] = EquilibriumGammaPilot.native_accounting(budget)
        summary["elapsed_seconds"] = time()-budget.started
        summary["data_bytes"] = EquilibriumGammaPilot.data_bytes(output)
        open(joinpath(output,"summary.toml"),"w") do io; TOML.print(io,summary); end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS)==2 || error("usage: run_gamma.jl INPUT.toml NEW_OUTPUT_DIRECTORY")
    run_pilot(abspath(ARGS[1]),abspath(ARGS[2]))
end
