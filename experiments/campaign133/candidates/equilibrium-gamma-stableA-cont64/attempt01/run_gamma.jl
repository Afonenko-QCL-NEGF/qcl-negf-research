# Experimental campaign harness. No Pkg operations, no Runner dependency.
using QCLNEGF, TOML, SHA, LinearAlgebra, Dates
include("gamma_functions.jl")
include("checks.jl")
include("warm_checks.jl")
include("hermitian_checks.jl")
include("stable_checks.jl")
include("restart_checks.jl")

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
    cfg["pi_calls_total"] == 64 || error("This frozen harness admits only64TOTAL Pi calls")
    cfg["pilot_seconds"] == 600.0 && cfg["overall_seconds"] == 1800.0 && cfg["data_max_bytes"] == 8*1024^3 || error("fixed affine experiment envelope differs")
    cfg["affine_l1_max"] == 4096.0 && cfg["history_differences"]==3 && cfg["linear_alpha"]==0.2 && cfg["aa_beta"]==1.0 && cfg["qr_condition_max"]==1e8 || error("AA controls differ")
    cfg["spectral_operator"]=="native_covariance_exact_plan_hlocal_all9" && cfg["project_AA_trial"]===false || error("stable-A only identity differs")
    cfg["qualifier_summary_sha256"]!="UNSET" || error("qualified operator result is not frozen")
    qualifier=EquilibriumGammaPilot.warm_path(cfg,"qualifier_summary_path")
    EquilibriumGammaPilot.checked_digest(qualifier,cfg["qualifier_summary_sha256"])
    q=TOML.parsefile(qualifier)
    q["source_sha256"]=="a34b1da110293ba289a78a437354436970ed8362a184779eb99e1f87690ade18" && q["execution"]=="complete" && q["maps_completed"]==3 || error("qualifier source/completion differs")
    q["q_cov"]["negative_count"]==0 && q["native_comparability"]["kinetic_KMS_control_pass"]===true || error("operator qualification did not pass")
    cfg["native_max_maps"] == 32 || error("Frozen native structural cap is1+30+1")
    cfg["initialization"] == "Gamma_checkpoint_emptyAA" && cfg["parent_call"]==16 || error("GammaLast16 continuation identity mismatch")
    cfg["retarded_real_symmetry"] == "hermitian_projection" || error("Hermitian-PI variant identity mismatch")
    0 < cfg["memory_max_bytes"] <= 13287555072 || error("OS/Core memory budget exceeds the approved envelope")
    0 < cfg["memory_high_bytes"] <= cfg["memory_max_bytes"] || error("memory.high must fit the hard cap")
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
    cfg["_input_directory"] = dirname(abspath(input))
    identity = source_gate(cfg) # NO model allocations before explicit gate
    ispath(output) && error("Refusing to overwrite output $output")
    mkpath(output)
    for file in ("run_gamma.jl","gamma_functions.jl","checks.jl","warm_seed.jl","warm_checks.jl","hermitian_checks.jl","stable_spectral.jl","stable_checks.jl","gamma_restart.jl","restart_checks.jl","aa_reference.jl","protocol.md","algorithm.md")
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
        for file in ("run_gamma.jl","gamma_functions.jl","checks.jl","warm_seed.jl","warm_checks.jl","hermitian_checks.jl","stable_spectral.jl","stable_checks.jl","gamma_restart.jl","restart_checks.jl","aa_reference.jl","protocol.md","algorithm.md"))
    identity["method"] = "custom Gamma AAII m3 affine4096; stable-A covariance exact plan.hlocal all9; GammaLast16 restart emptyAA/freshmu; original raw-h native exit; no PSD projector"
    open(joinpath(output,"identity.toml"),"w") do io; TOML.print(io,identity); end
    BLAS.set_num_threads(1)
    budget = EquilibriumGammaPilot.Budget(cfg)
    summary = Dict{String,Any}("scientific_acceptance"=>false)
    try
        restart_small_checks(cfg)
        summary["restart_metadata_budget_checks"] = "pass"
        summary["physical_104_checks"] = "reused unchanged physical104, qualified primitive tiny8 and stable-A16 toy integration/source hashes; no scientific tests repeated"
        summary["native_exit_scope"]="original raw-h Context/native public solver; reduced stable-A hH map acceptance alone cannot pass native gates"
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
        summary["domain_failure_records"] = sort(filter(name->startswith(name,"domain-failure-"),readdir(output)))
        open(joinpath(output,"summary.toml"),"w") do io; TOML.print(io,summary); end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if ARGS==["--tiny"]
        BLAS.set_num_threads(1)
        restart_small_checks(TOML.parsefile(joinpath(@__DIR__,"cont64.toml")))
    else
        length(ARGS)==2 || error("usage: run_gamma.jl --tiny OR INPUT.toml NEW_OUTPUT_DIRECTORY")
        run_pilot(abspath(ARGS[1]),abspath(ARGS[2]))
    end
end
