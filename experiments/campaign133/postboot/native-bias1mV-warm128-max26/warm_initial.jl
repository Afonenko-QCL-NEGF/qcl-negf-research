# Private unmodified full-state reader, adapted from archived run_inner_warm_v4.jl.
# No interpolation, retargeting, FD/number modification, or source array writes.
const COLD_DRIVER_SHA="1a25007ef353aabb7ee41de1474bf9043be45f0302eea6609943de42c6068726"
parent_info(c)=c["parents"][string(Int(c["voltage_mV"]))]
parent_path(c)=normpath(joinpath(@__DIR__,parent_info(c)["directory"]))
function parent_contract(c,pc,id,summary; tiny=false)
    c["warm_mode"]=="unmodified_private_fullstate" || error("warm mode differs")
    pc["max_scba"]==(tiny ? 128 : 32) && c["max_scba"]==128 || error("parent/new native budgets differ")
    expected_parent_caps=tiny ? (c["memory_max_bytes"],c["memory_high_bytes"]) : (12616466432,12348030976)
    (pc["memory_max_bytes"],pc["memory_high_bytes"])==expected_parent_caps || error("declared parent resource delta differs")
    (c["memory_max_bytes"],c["memory_high_bytes"])==(23085449216,22548578304) || error("declared child resource delta differs")
    allowed=("memory_max_bytes","memory_high_bytes","identity","runtime_authorized","max_scba","solver_seconds","overall_soft_seconds","data_max_bytes","tiny_cold_fixture")
    for (key,value) in pc
        key in allowed && continue
        haskey(c,key) && isequal(c[key],value) || error("warm changed inherited frame/control: $key")
    end
    id["core_head"]==c["core_head"] || error("parent Core differs")
    summary["execution"]=="complete" && 1<=summary["maps_completed"]<=(tiny ? 2 : 32) || error("parent has no complete aligned state")
    id["manifest_sha256"]==sha(joinpath(dirname(Base.active_project()),"Manifest.toml")) || error("parent Manifest differs")
    return nothing
end
function read_warm_field(group,name,shape)
    # Dataset metadata must agree BEFORE broadcast/materialization can hide a singleton axis.
    size(group[name*"_real"])==shape && size(group[name*"_imag"])==shape ||
        throw(DimensionMismatch("warm real/imag dataset shapes differ: $name"))
    X=GP.read_complex(group,name)
    X isa Array{ComplexF64,4} && size(X)==shape && all(isfinite,X) || error("warm field shape/type/finite: $name")
    return X
end
function load_initial(ctx,c,out;tiny=false)
    tiny && get(c,"tiny_cold_fixture",false) && return nothing
    info=parent_info(c);dir=parent_path(c);expected_driver=tiny ? sha(joinpath(@__DIR__,"run_native.jl")) : COLD_DRIVER_SHA
    info["driver_sha256"]==expected_driver || error("parent driver binding differs")
    for (file,key) in (("state.h5","state_sha256"),("model.h5","model_sha256"),("input.toml","input_sha256"),
                       ("identity.toml","identity_sha256"),("summary.toml","summary_sha256"),("run_native.jl","driver_sha256"))
        info[key]!="UNSET" && sha(joinpath(dir,file))==info[key] || error("parent whole SHA differs: $file")
    end
    pc=TOML.parsefile(joinpath(dir,"input.toml"));id=TOML.parsefile(joinpath(dir,"identity.toml"));sm=TOML.parsefile(joinpath(dir,"summary.toml"))
    id["input_sha256"]==info["input_sha256"] && id["script_sha256"]==info["driver_sha256"] || error("parent identity lineage")
    sm["state_sha256"]==info["state_sha256"] && sm["model_sha256"]==info["model_sha256"] || error("parent terminal ownership")
    parent_contract(c,pc,id,sm;tiny)
    proof=h5open(joinpath(dir,"model.h5"),"r") do f
        read(f["schema"])=="campaign133-gamma-model-v2" || error("parent model schema")
        GP.verify_problem_model(ctx,f) # all grids/profiles/basis/h/raw kernels/qK, exact before loading large state
    end
    p=ctx.problem;shape=(p.numerical.N_E,p.numerical.N_k,p.numerical.N_b,p.numerical.N_b)
    state=h5open(joinpath(dir,"state.h5"),"r") do f
        read(f["schema"])=="campaign133-native-bias-fullstate-v1" || error("parent state schema")
        read(f["voltage_per_period_mV"])==c["voltage_mV"] || error("parent voltage differs")
        isequal(read(f["energy"]),p.grids.ε) && isequal(read(f["weight_energy"]),p.grids.wᴱ) && isequal(read(f["weight_k"]),p.grids.wᵏ) || error("parent axes differ")
        U=read(f["hartree"]);length(U)==p.numerical.N_z && all(iszero,U) || error("parent U0 differs")
        read(f["energy_scale_eV"])==p.scales.E₀_eV && read(f["length_scale_m"])==p.scales.L₀_m && read(f["spin_degeneracy"])==p.physical.g_s || error("parent units differ")
        g=f["green"];GR,GL,GG,A=(read_warm_field(g,n,shape) for n in ("GR","GL","GG","A"))
        condition=read(g["condition_number"]);scale=read(g["dyson_scale"])
        size(condition)==size(scale)==shape[1:2] && all(isfinite,condition) && all(isfinite,scale) || error("parent Green diagnostics differ")
        green=GreenState(GR,GL,GG,A,condition,scale)
        expected=Set(string.([p.kernels.enabled;:embedding;:embedding_plus;:embedding_minus]))
        Set(keys(f["selfenergy"]))==expected || error("parent self-energy channels differ")
        families=Dict{Symbol,SelfEnergyFamily}()
        for name in expected
            group=f["selfenergy/"*name];R,L,G=(read_warm_field(group,n,shape) for n in ("retarded","lesser","greater"))
            families[Symbol(name)]=SelfEnergyFamily(R,L,G)
        end
        contract=TOML.parse(read(f["restart_contract_toml"]))
        SCBAResult(green,Dict(name=>families[name] for name in p.kernels.enabled),families[:embedding],families[:embedding_plus],families[:embedding_minus],
            SCBAIteration[],false,:warm_initial_guess,:unresolved,contract)
    end
    write_toml(joinpath(out,"warm-initialization.toml"),Dict("mode"=>"unmodified private full28fields; no FD/interpolation/rescale",
        "parent"=>dir,"source_sha256"=>info,"actual_frame_digests"=>proof,"parent_completed"=>sm["maps_completed"],
        "parent_converged"=>get(sm,"converged",false),"resume_iterations"=>false,"consume_initial"=>true,"history_reset"=>true,
        "scope"=>"source is initialization only; all new native maps/gates fresh, no parent residual accepted as convergence"))
    return state
end
function boundary_audit(ctx,state,c,out)
    dir=parent_path(c);prior=TOML.parsefile(joinpath(dir,"summary.toml"))["last_row"];row=last(state.history)
    fields=Dict{String,Any}();greenpass=true;tol=256eps(Float64)
    h5open(joinpath(dir,"state.h5"),"r") do f
        for (name,X) in (("GR",state.green.Gᴿ),("GL",state.green.Gˡ),("GG",state.green.Gᵍ),("A",state.green.A))
            d2=0.0;s2=0.0;dw=0.0;sw=0.0;exact=true
            for start in 1:1024:size(X,1)
                rows=start:min(start+1023,size(X,1));Y=GP.read_complex_rows(f["green"],name,rows)
                @inbounds for b in axes(Y,4),a in axes(Y,3),m in axes(Y,2),e in axes(Y,1)
                    global_e=start+e-1;v=X[global_e,m,a,b];ref=Y[e,m,a,b];q=abs2(v-ref);r=abs2(ref)
                    weight=ctx.problem.grids.wᴱ[global_e]*ctx.problem.grids.wᵏ[m]
                    d2+=q;s2+=r;dw+=weight*q;sw+=weight*r;exact &= isequal(v,ref)
                end
            end
            relative=GP.relative(sqrt(d2),sqrt(s2));weighted=GP.relative(sqrt(dw),sqrt(sw));greenpass &= relative<=tol && weighted<=tol
            fields[name]=Dict("absolute"=>sqrt(d2),"scale"=>sqrt(s2),"relative"=>relative,"weighted_absolute"=>sqrt(dw),"weighted_scale"=>sqrt(sw),"weighted_relative"=>weighted,"bitwise_equal"=>exact)
        end
    end
    scalars=Dict{String,Any}()
    for name in (:r_D,:r_A,:r_K,:r_Σ,:r_λ,:λ,:r_PSD,:r_caus,:r_roundoff,:J,:raw_charge,:target_charge,:normalized_charge)
        value=getfield(row,name);old=prior[string(name)]
        scalars[string(name)]=Dict("fresh"=>value,"parent"=>old,"absolute_difference"=>abs(value-old),"bitwise_equal"=>isequal(value,old))
    end
    write_toml(joinpath(out,"boundary-replay.toml"),Dict("green_replay_passed"=>greenpass,"green_relative_allowance"=>tol,
        "green_fields"=>fields,"fresh_equation_rows"=>scalars,"new_iteration"=>row.ν,"parent_iteration"=>prior["ν"],
        "reset_fields"=>"nu/history/streak/Jchange/population/seed_mu intentionally reset; not equivalence gates",
        "scope"=>"same full state must reproduce Green; scalar raw equation differences recorded, not convergence acceptance"))
    greenpass || error("warm first-boundary Green did not replay parent within machine allowance")
    return nothing
end
