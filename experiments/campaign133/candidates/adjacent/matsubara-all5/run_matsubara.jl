# ROOT-owned runtime only. No Pkg, Γ loading, μ fitting or automatic extension.
const PROCESS_STARTED=time()
using TOML, SHA, LinearAlgebra, QCLNEGF, HDF5, Dates
sha(p)=bytes2hex(open(sha256,p))
checksha(p,d)=(sha(p)==d || error("SHA mismatch: $p"))
resolve(c,k)=normpath(joinpath(@__DIR__,c[k]))
plain(x::NamedTuple)=Dict(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::AbstractDict)=Dict(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::Complex)=Dict("real"=>real(x),"imag"=>imag(x))
plain(x::AbstractVector)=[plain(v) for v in x]
plain(x)=x
bytesout(out)=sum(filesize(joinpath(r,f)) for (r,_,fs) in walkdir(out) for f in fs;init=0)
function toml_atomic(path,x,c)
    io=IOBuffer();TOML.print(io,plain(x));b=take!(io)
    length(b)<=c["single_file_max_bytes"] || error("TOML file cap")
    bytesout(dirname(path))+length(b)<=c["data_max_bytes"] || error("TOML data cap")
    write(path*".tmp",b);mv(path*".tmp",path;force=true)
end
checktime(c)=(time()-PROCESS_STARTED<c["wall_seconds"] || error("global600s budget"))
function gate(c)
    frozen=TOML.parsefile(joinpath(@__DIR__,"runtime-input.toml"))
    for (k,v) in frozen
        k=="runtime_authorized" && continue
        haskey(c,k) && isequal(c[k],v) || error("frozen setting changed: $k")
    end
    c["runtime_authorized"]===true || error("runtime_authorized=false")
    checksha(@__FILE__,c["adapter_sha256"])
    core=dirname(dirname(pathof(QCLNEGF)))
    readchomp(`git -C $core rev-parse HEAD`)==c["core_head"] || error("Core HEAD")
    readchomp(`git -C $core rev-parse $("HEAD^{tree}")`)==c["core_tree"] || error("Core tree")
    isempty(readchomp(`git -C $core status --porcelain --untracked-files=no`)) || error("dirty Core")
    string(VERSION)==c["julia_version"] && Threads.nthreads()==1 || error("Julia identity")
    project=Base.active_project();project===nothing && error("missing accepted project")
    checksha(project,c["project_sha256"]);manifest=joinpath(dirname(project),"Manifest.toml")
    checksha(manifest,c["manifest_sha256"])
    member=only(filter(s->startswith(s,"0::"),readlines("/proc/self/cgroup")))
    group=joinpath("/sys/fs/cgroup",lstrip(split(member,':';limit=3)[3],'/'))
    for (file,key) in (("memory.max","memory_max_bytes"),("memory.swap.max","swap_max_bytes"))
        strip(read(joinpath(group,file),String))==string(c[key]) || error("cgroup $file")
    end
    cpu=split(strip(read(joinpath(group,"cpu.max"),String)))
    cpu[1]!="max" && 0<parse(Int,cpu[1])/parse(Int,cpu[2])<=2 || error("CPU≤200% required")
    limit=split(only(filter(s->startswith(s,"Max file size"),readlines("/proc/self/limits"))))
    limit[4]!="unlimited" && parse(Int,limit[4])<=c["single_file_max_bytes"] || error("LimitFSIZE32MiB")
    for (name,d) in c["definition_sha256"];checksha(joinpath(resolve(c,"source_directory"),name),d);end
    for (key,d) in c["text_sha256"];checksha(resolve(c,key),d);end
    return (;project,manifest,group,core_source=pathof(QCLNEGF))
end
function load_modules(c)
    Base.include(Main,joinpath(resolve(c,"source_directory"),"gamma_functions.jl"))
    Base.include(Main,resolve(c,"digest_helper"))
    Base.include(Main,resolve(c,"prototype"))
end
function source(c)
    parent=resolve(c,"source_directory")
    for (name,d) in c["source_file_sha256"]
        checksha(joinpath(parent,name),d);checktime(c)
    end
    pc=TOML.parsefile(joinpath(parent,"input.toml"));id=TOML.parsefile(joinpath(parent,"identity.toml"))
    mp=TOML.parsefile(joinpath(parent,"map-064.toml"));ledger=TOML.parsefile(joinpath(parent,"spectral-ledger-064.toml"))
    sm=TOML.parsefile(joinpath(parent,"summary.toml"));dg=TOML.parsefile(resolve(c,"digest_packet"))
    id["core_head"]==c["core_head"] && id["core_tree"]==c["core_tree"] && isempty(id["core_dirty"]) || error("parent Core")
    id["input_sha256"]==c["source_file_sha256"]["input.toml"] && id["manifest_sha256"]==c["manifest_sha256"] || error("parent frame")
    pc["project_AA_trial"]===false && pc["temperature_K"]==70 || error("model physical frame")
    (pc["NE"],pc["Nk"],pc["Nb"],pc["Nz"],pc["basis_periods"])==(30704,11,5,96,3) || error("axes/basis")
    sm["execution"]=="complete" && sm["pi_calls"]==64 && sm["reduced"]["last_map"]==mp || error("parent completion")
    ledger==mp["spectral_operator_ledger"] && mp["call"]==64 || error("spectral ownership")
    dg["status"]=="completed" && dg["mu_dimensionless"]==mp["mu"]["mu"]==c["mu_native"] || error("digest/μ binding")
    dg["core_head"]==c["core_head"] && dg["core_tree"]==c["core_tree"] || error("digest Core")
    for (key,name) in (("model_h5","model.h5"),("spectral_selected","spectral-selected-064.h5"),
        ("source_input","input.toml"),("source_identity","identity.toml"),("map_summary","map-064.toml"),
        ("source_ledger","spectral-ledger-064.toml"),("source_summary","summary.toml"))
        dg["source_file_sha256"][key]==c["source_file_sha256"][name] || error("digest ownership: $key")
    end
    GP=Main.EquilibriumGammaPilot;D=Main.GammaDensityHermitianFreezeIO
    pc["_input_directory"]=parent;ctx=GP.build_context(pc) # source options preserved, NOT an admitted native run
    checktime(c)
    proof=h5open(joinpath(parent,"model.h5"),"r") do f
        read(f["schema"])=="campaign133-gamma-model-v2" || error("model schema")
        GP.verify_problem_model(ctx,f)
    end
    hH=h5open(joinpath(parent,"spectral-selected-064.h5"),"r") do f
        read(f["schema"])=="campaign133-stableA-spectral-selected-v1" || error("spectral schema")
        isequal(D.read_complex(f,"h_raw"),ctx.h) || error("saved/model raw h")
        isequal(read(f["energy"]),ctx.problem.grids.ε) && isequal(read(f["weight_k"]),ctx.problem.grids.wᵏ) || error("spectral axes")
        read(f["energy_scale_eV"])==ctx.problem.scales.E₀_eV==0.1 || error("E0 units")
        isequal(D.read_complex(f,"Tplus"),ctx.problem.basis.T₊) || error("spectral Tplus")
        TOML.parse(read(f["ledger_toml"]))==ledger || error("saved ledger")
        h=D.read_complex(f,"h_H")
        size(h)==(11,5,5) && all(isfinite,h) && D.sha_array(h)==dg["array_sha256"]["hcell"] || error("saved h_H digest")
        all(ishermitian(Matrix(h[k,:,:])) for k in 1:11) || error("exact Hermitian h_H required")
        h
    end
    p=ctx.problem;all(p.basis.T₋[a,b]==conj(p.basis.T₊[b,a]) for a in 1:5,b in 1:5) || error("exact adjoint coupling")
    D.sha_array(p.grids.wᵏ)==dg["array_sha256"]["weight_momentum"] || error("radial measure digest")
    sp=QCLNEGF.QCLNumerics._scaled_physics(p.physical,p.scales)
    iszero(sp.Eᵖ) && p.scales.L₀_m==1e-8 && p.physical.g_s==2 || error("V0/length/spin units")
    beta=inv(ctx.kBT);Omega=sp.ħωᴸᴼ
    isfinite(beta) && beta>0 && Omega>0 && p.models.lo_population===:thermal || error("bath units")
    m=Main.MatsubaraAll5.context(hH,p.basis.T₊,p.kernels.K,p.grids.wᵏ,p.kernels.qᴷ,beta,c["mu_native"],Omega;spin=p.physical.g_s)
    binding=(;model_proof=proof,hH_sha256=D.sha_array(hH),Tplus_sha256=D.sha_array(p.basis.T₊),
        beta_native=beta,Omega_native=Omega,mu_native=c["mu_native"],mu_eV=c["mu_native"]*p.scales.E₀_eV,
        E0_eV=p.scales.E₀_eV,L0_m=p.scales.L₀_m,spin=p.physical.g_s,
        source_cache_full_native_estimate=ctx.cache.estimate.peak_bytes,
        note="source model reconstruction only; no Γ arrays, native solve or inherited PD acceptance")
    ctx=nothing;GC.gc(true);return m,binding
end
function save_state(path,snapshot,c;full=false,replace=false)
    GP=Main.EquilibriumGammaPilot
    ispath(path) && !replace && error("state exists")
    fields=full ? (("Sigma",snapshot.S),("G",snapshot.G),("Phi",snapshot.Phi),("rho",snapshot.rho)) :
        (("Sigma",snapshot.S),("G",snapshot.G),("rho",snapshot.rho))
    expected=sum(GP.expected_padded_complex_bytes(size(X)) for (_,X) in fields)+65536
    expected<=c["single_file_max_bytes"] && bytesout(dirname(path))+expected<=c["data_max_bytes"] || error("raw-state admission cap")
    h5open(path*".tmp","w") do f
        f["schema"]="campaign133-Mats-ALL5-aligned-raw-v1";f["omega_native"]=snapshot.w
        f["iteration"]=snapshot.record.iteration;f["mu_native"]=c["mu_native"]
        f["record_toml"]=sprint(io->TOML.print(io,plain(snapshot.record)))
        for (name,X) in fields;GP.write_complex(f,name,X);end
    end
    filesize(path*".tmp")<=c["single_file_max_bytes"] && bytesout(dirname(path))<=c["data_max_bytes"] || error("raw tmp cap; preserved")
    mv(path*".tmp",path;force=replace)
end
function run_actual(input,out)
    c=TOML.parsefile(input);environment=gate(c)
    if !isdefined(Main,:MatsubaraAll5)
        load_modules(c);return Base.invokelatest(run_actual,input,out)
    end
    ispath(out) && error("new output required");mkpath(out)
    cp(input,joinpath(out,"input.toml"));cp(environment.project,joinpath(out,"Project.toml"))
    cp(environment.manifest,joinpath(out,"Manifest.toml"))
    for key in ("prototype","digest_helper");cp(resolve(c,key),joinpath(out,basename(resolve(c,key))));end
    cp(@__FILE__,joinpath(out,"run_matsubara.jl"))
    result=Dict{String,Any}("execution"=>"started","scientific_acceptance"=>false,"native_pass"=>false,
        "variant"=>"fixedμ/all5/direct-static-C/full-line-LO/reference-tail/finite45central-repeat/hH-all9",
        "source_Gamma_used"=>false,"inner_requested"=>0,"inner_completed"=>0,"completed_N"=>Int[])
    toml_atomic(joinpath(out,"identity.toml"),(;environment...,input_sha256=sha(input),
        definition_sha256=c["definition_sha256"],text_sha256=c["text_sha256"],source_file_sha256=c["source_file_sha256"]),c)
    BLAS.set_num_threads(1)
    try
        checktime(c);m,binding=source(c);toml_atomic(joinpath(out,"source-binding.toml"),binding,c)
        M=Main.MatsubaraAll5;checktime(c)
        opcheck=M.qualify_native_ops(m);pole=M.onelevel_oracle(m.beta,m.Omega;N=256,g2=m.Omega^2/5,epsilon=m.Omega/4)
        wired=M.wired_onelevel(m.beta,m.Omega)
        toml_atomic(joinpath(out,"controls.toml"),(;opcheck,pole,wired,cp=M.cp_bound(m)),c)
        all(r.G_relative<=2e-12 && r.reference_LO_relative<=2e-12 for r in wired) || error("actual onelevel wiring")
        all(pole[i+1].absolute<=pole[i].absolute+1e-12 for i in (1,3,5)) || error("independent onepole refinement")
        wired[2].LO_absolute<=wired[1].LO_absolute+1e-12 && wired[2].density_absolute<=wired[1].density_absolute+1e-12 || error("onelevel tail refinement")
        states=Dict{Int,Any}()
        for N in c["N_values"]
            checktime(c);result["inner_requested"]+=16
            toml_atomic(joinpath(out,"request-N$(N).toml"),(;N,requested_inner=16,total_requested=result["inner_requested"],
                completed=result["inner_completed"],entered_lower=result["inner_completed"],entered_upper=result["inner_requested"]),c)
            callback=function(snap)
                tag="N$(N)-i$(lpad(snap.record.iteration,2,'0'))"
                save_state(joinpath(out,tag*".h5"),snap,c)
                save_state(joinpath(out,"N$(N)-latest.h5"),snap,c;full=true,replace=true)
                result["inner_completed"]+=1
                toml_atomic(joinpath(out,tag*".toml"),snap.record,c)
                result["current_N"]=N;result["last_record"]=plain(snap.record)
                toml_atomic(joinpath(out,"summary.toml"),result,c);checktime(c)
            end
            remain=c["wall_seconds"]-(time()-PROCESS_STARTED)-2
            remain>0 || error("global wall before pilot")
            state=M.pilot(m,N;steps=16,alpha=c["alpha"],soft_seconds=remain,callback)
            states[N]=state;result["status_N$(N)"]=state.status
            push!(result["completed_N"],N)
            state.status=="iteration_budget" || break
        end
        if all(haskey(states,N) && hasproperty(states[N],:rho) for N in (128,256))
            result["cutoff_comparison"]=plain(M.compare_cutoffs(states[128],states[256],m))
        end
        result["execution"]="complete_bounded_diagnostic"
    catch err
        result["execution"]="stopped";result["reason"]=sprint(showerror,err);rethrow()
    finally
        result["elapsed_seconds_including_imports"]=time()-PROCESS_STARTED
        result["maxrss_bytes"]=Sys.maxrss();result["finish_utc"]=string(Dates.now(Dates.UTC))
        result["native_pass"]=false;result["density_normalized"]=false
        toml_atomic(joinpath(out,"summary.toml"),result,c)
    end
    return result
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("INPUT_TOML NEW_OUTPUT")
    run_actual(abspath(ARGS[1]),abspath(ARGS[2]))
end
