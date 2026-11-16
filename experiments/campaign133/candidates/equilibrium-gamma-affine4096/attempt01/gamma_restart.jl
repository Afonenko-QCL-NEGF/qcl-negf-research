# Γ-only checkpoint loader. Existing warm_seed.jl supplies only digest/model
# utilities here; load_warm_seed/builder/interpolation are never called.
const PARENT_GAMMA_OPERATOR_SHA="57f2c277255ae1f8265659e65176cc88f97d5a40cbf8d6e32b7ea9257ccfb032"
const REPLAY_RELATIVE=256eps(Float64) # machine replay guard, NOT solver tolerance
restart_parent(cfg)=warm_path(cfg,"parent_directory")
restart_map_name(cfg)="map-"*lpad(cfg["parent_call"],3,'0')*".toml"
function restart_contract(cfg,parent,identity)
    cfg["initialization"]=="Gamma_checkpoint_emptyAA" || error("restart mode differs")
    cfg["parent_operator_sha256"]==PARENT_GAMMA_OPERATOR_SHA &&
        identity["source_sha256"]["gamma_functions.jl"]==PARENT_GAMMA_OPERATOR_SHA || error("parent operator differs")
    identity["core_head"]==cfg["expected_core_head"] && identity["core_tree"]==cfg["expected_core_tree"] || error("parent Core differs")
    parent["retarded_real_symmetry"]==cfg["retarded_real_symmetry"]=="hermitian_projection" || error("PI symmetry scheme differs")
    allowed=("identity","runtime_authorized","parent_directory","parent_input_sha256","parent_identity_sha256",
        "parent_checkpoint_sha256","parent_terminal_sha256","parent_map_sha256","parent_operator_sha256",
        "target_model_path","target_model_sha256","pi_calls_total","pilot_seconds","overall_seconds","affine_l1_max")
    for (key,value) in parent
        key in allowed && continue
        haskey(cfg,key) && isequal(cfg[key],value) || error("restart changed inherited option: $key")
    end
    cfg["parent_call"]==128 && parent["pi_calls_total"]==128 && cfg["pi_calls_total"]==32 || error("fixed call budgets differ")
    parent["affine_l1_max"]==10.0 && cfg["affine_l1_max"]==4096.0 || error("single affine parameter variant differs")
    cfg["pilot_seconds"]==600.0 && cfg["overall_seconds"]==1800.0 && cfg["data_max_bytes"]==8*1024^3 || error("finite affine experiment envelope differs")
    return nothing
end
function load_gamma_restart(c::Context,cfg,output)
    parent=restart_parent(cfg)
    for (name,key) in (("input.toml","parent_input_sha256"),("identity.toml","parent_identity_sha256"),
                       ("Gamma-last.h5","parent_checkpoint_sha256"),(restart_map_name(cfg),"parent_map_sha256"),
                       ("reduced-terminal.h5","parent_terminal_sha256"),("gamma_functions.jl","parent_operator_sha256"))
        checked_digest(joinpath(parent,name),cfg[key])
    end
    options=TOML.parsefile(joinpath(parent,"input.toml"));identity=TOML.parsefile(joinpath(parent,"identity.toml"))
    identity["input_sha256"]==cfg["parent_input_sha256"] || error("parent input/identity mismatch")
    restart_contract(cfg,options,identity)
    modelpath=warm_path(cfg,"target_model_path");checked_digest(modelpath,cfg["target_model_sha256"])
    normpath(modelpath)==normpath(joinpath(parent,"model.h5")) || error("model is not the pinned parent frame")
    fingerprints=HDF5.h5open(modelpath,"r") do file;verify_problem_model(c,file);end
    gamma=HDF5.h5open(joinpath(parent,"Gamma-last.h5"),"r") do file
        read(file["schema"])=="campaign133-gamma-checkpoint-v1" || error("checkpoint schema differs")
        read(file["identity"])==options["identity"] && read(file["pi_calls_total"])==cfg["parent_call"] || error("checkpoint parent identity/call differs")
        expected=(c.problem.numerical.N_E,c.problem.numerical.N_k,c.problem.numerical.N_b,c.problem.numerical.N_b)
        size(file["Gamma_real"])==expected && size(file["Gamma_imag"])==expected || error("checkpoint axes/shape differs")
        read_complex(file,"Gamma")::Field # no symmetrization/rescaling/interpolation
    end
    witness=cone(gamma)
    record=Dict{String,Any}("parent_checkpoint_sha256"=>cfg["parent_checkpoint_sha256"],
        "parent_input_sha256"=>cfg["parent_input_sha256"],"parent_identity_sha256"=>cfg["parent_identity_sha256"],
        "parent_model_sha256"=>cfg["target_model_sha256"],"parent_operator_sha256"=>PARENT_GAMMA_OPERATOR_SHA,
        "actual_problem_fingerprints"=>fingerprints,"input_cone"=>Dict(string(k)=>v for (k,v) in pairs(witness)),
        "scope"=>"Gamma only, fresh mu, empty AA/streak; no R/L/G/embedding/history imported; checkpoint axes bound by parent input+identity+model SHA")
    open(joinpath(output,"restart-initialization.toml"),"w") do io;TOML.print(io,record);end
    return gamma
end
machine_same(a::Real,b::Real)=isfinite(a) && isfinite(b) && (a==b || abs(a-b)<=REPLAY_RELATIVE*max(abs(a),abs(b)))
function scalar_boundary_metrics(current,parent)
    fields=(("residual","absolute"),("residual","scale"),("residual","relative"),
        ("residual","weighted_absolute"),("residual","weighted_scale"),("residual","weighted_relative"),
        ("mu","mu"),("mu","number"),("mu","derivative"))
    rows=Dict{String,Any}();passed=true
    for (group,key) in fields
        a,b=current[group][key],parent[group][key];ok=machine_same(a,b);passed &= ok
        rows[group*"."*key]=Dict("current"=>a,"parent"=>b,"absolute_difference"=>abs(a-b),"passed"=>ok)
    end
    return Dict{String,Any}("passed"=>passed,"fields"=>rows,"relative_allowance"=>REPLAY_RELATIVE,
        "zero_policy"=>"zero/zero exact; no additive norm floor")
end
function restart_boundary_audit(c::Context,rawphi::Field,stats,cfg,output)
    parent=restart_parent(cfg);prior=TOML.parsefile(joinpath(parent,restart_map_name(cfg)))
    report=scalar_boundary_metrics(stats,prior);metric=WarmNorm();g=c.problem.grids
    HDF5.h5open(joinpath(parent,"reduced-terminal.h5"),"r") do file
        size(file["raw_Phi_real"])==size(rawphi) || error("parent raw Phi axes differ")
        for start in 1:1024:size(rawphi,1)
            rows=start:min(start+1023,size(rawphi,1))
            reference=read_complex_rows(file,"raw_Phi",rows)
            # Small independent row buffers only; no parent Green/state loading.
            accumulate_norm!(metric,reference,rawphi[rows,:,:,:],g.wᴱ,g.wᵏ,start-1)
        end
    end
    norms=warm_norm_dict(metric);report["raw_Phi_field_replay"]=norms
    report["passed"] &= norms["relative"]<=REPLAY_RELATIVE && norms["weighted_relative"]<=REPLAY_RELATIVE
    report["parent_call"],report["new_call"] =cfg["parent_call"],1
    open(joinpath(output,"boundary-replay.toml"),"w") do io;TOML.print(io,report);end
    report["passed"] || error("fresh first map did not repeat parent frozen call within machine allowance")
    return nothing
end
