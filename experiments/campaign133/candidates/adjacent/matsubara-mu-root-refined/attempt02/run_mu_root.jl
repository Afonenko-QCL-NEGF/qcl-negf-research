# WRITE-ONLY preparation. ROOT owns authorization, M3 qualification and execution.
module BoundedMuRootRefined
using TOML,SHA,LinearAlgebra,HDF5,QCLNEGF
const STARTED=time()
blas1()=(BLAS.set_num_threads(1);BLAS.get_num_threads()==1 || error("BLAS getter is not1");1)
blas1() # OpenMP OpenBLAS may ignore OPENBLAS_NUM_THREADS; explicit setting is required.
sha(p)=bytes2hex(open(sha256,p))
path(c,k)=normpath(joinpath(@__DIR__,c[k]))
needsha(p,h)=(sha(p)==h || error("SHA mismatch: $p"))
function preflight(c)
    frozen=TOML.parsefile(joinpath(@__DIR__,"input.toml"))
    for (k,v) in frozen
        k=="runtime_authorized" && continue
        haskey(c,k) && isequal(c[k],v) || error("frozen wrapper setting: $k")
    end
    c["runtime_authorized"]===true && c["M3_qualified"]===true || error("HOLD: authorization/M3 qualification")
    needsha(@__FILE__,c["wrapper_sha256"])
    for k in ("source_adapter","source_module","source_input","refined_module","M3_module","M3_review","M3_qualification","M3_qualification_input")
        needsha(path(c,k),c[k*"_sha256"])
    end
    all(isfinite,(c["source_mu_native"],c["bracket_halfwidth_native"],c["target_native"],c["wall_seconds"])) || error("finite scalar inputs")
    c["target_native"]>0 && c["bracket_halfwidth_native"]>0 || error("target/bracket domain")
    s=TOML.parsefile(path(c,"source_input"));s["mu_native"]==c["source_mu_native"] || error("source μ")
    core=dirname(dirname(pathof(QCLNEGF)))
    readchomp(`git -C $core rev-parse HEAD`)==s["core_head"] || error("Core HEAD")
    readchomp(`git -C $core rev-parse $("HEAD^{tree}")`)==s["core_tree"] || error("Core tree")
    isempty(readchomp(`git -C $core status --porcelain --untracked-files=no`)) || error("dirty Core")
    string(VERSION)==s["julia_version"] && Threads.nthreads()==1 || error("Julia identity")
    BLAS.get_num_threads()==1 || error("BLAS1 required")
    project=Base.active_project();project===nothing && error("missing project")
    needsha(project,s["project_sha256"]);needsha(joinpath(dirname(project),"Manifest.toml"),s["manifest_sha256"])
    member=only(filter(x->startswith(x,"0::"),readlines("/proc/self/cgroup")))
    group=joinpath("/sys/fs/cgroup",lstrip(split(member,':';limit=3)[3],'/'))
    for (f,k) in (("memory.max","memory_max_bytes"),("memory.swap.max","swap_max_bytes"))
        strip(read(joinpath(group,f),String))==string(c[k]) || error("cgroup $f")
    end
    cpu=split(strip(read(joinpath(group,"cpu.max"),String)))
    cpu[1]!="max" && 0<parse(Int,cpu[1])/parse(Int,cpu[2])<=2 || error("CPU cap")
    fs=split(only(filter(x->startswith(x,"Max file size"),readlines("/proc/self/limits"))))
    fs[4]!="unlimited" && parse(Int,fs[4])<=c["single_file_max_bytes"] || error("file cap")
    base=dirname(path(c,"source_input"))
    for (name,h) in s["definition_sha256"];needsha(normpath(joinpath(base,s["source_directory"],name)),h);end
    for (k,h) in s["text_sha256"];needsha(normpath(joinpath(base,s[k])),h);end
    qualify_M3(c,s)
    return s,(;project,group,BLAS_threads=BLAS.get_num_threads(),core_head=s["core_head"],core_tree=s["core_tree"])
end
function qualify_M3(c,s)
    q=TOML.parsefile(path(c,"M3_qualification"));qi=TOML.parsefile(path(c,"M3_qualification_input"))
    q["execution"]=="completed_postprocess_only" && q["inner_maps"]==q["mu_fits"]==0 || error("M3 incomplete/wrong scope")
    q["BLAS_threads"]==q["BLAS_threads_final"]==1 || error("M3 actual BLAS qualification")
    q["text_sha256"]==qi["text_sha256"] && q["source_sha256"]==qi["source_sha256"] || error("M3 input/report binding")
    q["text_sha256"]["moment_module"]==c["M3_module_sha256"] && q["text_sha256"]["old_adapter"]==c["source_adapter_sha256"] || error("M3 module/API adapter binding")
    q["text_sha256"]["adapter"]==c["M3_postprocess_adapter_sha256"] || error("M3 qualified adapter")
    qi["mu_native"]==c["source_mu_native"] && qi["runtime_authorized"]===true || error("M3 source μ/auth")
    qi["actual_directory"]==c["M3_source_directory"] || error("M3 source locator")
    Set(keys(q["source_sha256"]))==Set(("N128-latest.h5","N256-latest.h5","summary.toml","input.toml",
        "source-binding.toml","N128-i16.toml","N256-i16.toml","identity.toml")) || error("M3 complete source packet")
    for (file,h) in q["source_sha256"];needsha(joinpath(path(c,"M3_source_directory"),file),h);end
    old=TOML.parsefile(joinpath(path(c,"M3_source_directory"),"input.toml"))
    all(k=="runtime_authorized" || isequal(old[k],v) for (k,v) in s) || error("M3 source scientific config")
    for N in (128,256)
        r=q["N$N"];all(isfinite,(r["M3_number"],r["eta"],r["linear_residual_relative"])) || error("M3 qualification finite")
        r["N"]==N && r["eta"]<1 && r["linear_residual_relative"]<=c["M3_linear_relative_tolerance"] || error("M3 qualification linear solve")
    end
end
function check(c,out;extra=0,final=false)
    final || time()-STARTED<c["wall_seconds"] || error("global wall budget")
    reserve=final ? 0 : 1048576
    Main.bytesout(out)+extra+reserve<=c["data_max_bytes"] || error("global data budget")
end
function emit(c,out,file,x;final=false)
    check(c,out;final);io=IOBuffer();TOML.print(io,Main.plain(x));b=take!(io)
    length(b)<=c["single_file_max_bytes"] || error("TOML file budget")
    check(c,out;extra=length(b),final);p=joinpath(out,file);write(p*".tmp",b);mv(p*".tmp",p;force=true)
end
function proposal(a,b)
    fa=a.target_relative;fb=b.target_relative
    mu=b.mu_native-fb*(b.mu_native-a.mu_native)/(fb-fa)
    (!isfinite(mu) || !(a.mu_native<mu<b.mu_native)) ? (a.mu_native+b.mu_native)/2 : mu
end
"Source-free toy checks; ROOT may call after static review. No solver/file IO."
function sourcefree_controls()
    c=Dict("steps_per_call"=>16,"inner_required_consecutive"=>3,
        "inner_raw_relative_tolerance"=>1e-9,"inner_raw_absolute_tolerance_native"=>1e-12)
    row(a,r)=(;raw_sigma_absolute=a,raw_sigma_relative=r)
    good=fill(row(1e-13,1e-10),16);@assert qualified(good,c)
    for bad in (row(2e-12,1e-10),row(1e-13,2e-9),row(Inf,Inf),row(NaN,0.0))
        h=copy(good);h[end-2]=bad;@assert !qualified(h,c)
    end
    @assert !qualified(good[1:15],c)
    V=ones(ComplexF64,1,1);base=(;mu=0.1,atoms=[(;eps=[0.2],V)],full=:unchanged)
    shifted=at_mu(base,0.11);@assert shifted.full===base.full && shifted.atoms[1].V===V
    @assert isapprox(shifted.atoms[1].eps[1],0.19;atol=1e-15)
    a=(;mu_native=0.0,target_relative=-2.0);b=(;mu_native=3.0,target_relative=1.0)
    @assert proposal(a,b)==2.0
    @assert proposal(a,(;mu_native=3.0,target_relative=-2.0))==1.5
    return (;pass=true,scope="pure qualification/context-shift/secant metadata toys; no scientific run")
end
"Only μ/relative pole energies change; all physical Hamiltonians/vertices stay fixed."
function at_mu(c,mu)
    delta=mu-c.mu
    atoms=[(;eps=a.eps.-delta,V=a.V) for a in c.atoms]
    merge(c,(;mu=Float64(mu),atoms))
end
function qualified(h,c)
    length(h)==c["steps_per_call"] || return false
    all(h[end-j].raw_sigma_relative<=c["inner_raw_relative_tolerance"] &&
        h[end-j].raw_sigma_absolute<=c["inner_raw_absolute_tolerance_native"] &&
        isfinite(h[end-j].raw_sigma_relative) && isfinite(h[end-j].raw_sigma_absolute)
        for j in 0:c["inner_required_consecutive"]-1)
end
function weighted_raw(s,c)
    num=sum(c.wk[k]*sum(abs2,s.Phi[:,k,:,:]-s.S[:,k,:,:]) for k in 1:c.nk)
    den=sum(c.wk[k]*sum(abs2,s.Phi[:,k,:,:]) for k in 1:c.nk)
    (;absolute=sqrt(num),scale=sqrt(den),relative=Main.MatsubaraAll5.relative(sqrt(num),sqrt(den)),
       scope="wk-weighted diagnostic only; qualification uses unweighted Frobenius")
end
function terminal(c,out,dir,m,s,cor)
    p=joinpath(dir,"density-M3.h5");check(c,out;extra=65536)
    h5open(p,"w") do f
        f["schema"]="campaign133-Mats-ALL5-mu-root-density-v1";f["mu_native"]=m.mu;f["candidate_accepted"]=false
        for (name,a) in (("raw_rho",s.rho),("oneshot_rho",cor.oneshot),("linear_M3_rho",cor.rho),
            ("S0",cor.S0),("B0",cor.B0))
            Main.EquilibriumGammaPilot.write_complex(f,name,a)
        end
    end
    filesize(p)<=c["single_file_max_bytes"] || error("density file cap");check(c,out)
end
function run_loaded(c,sconf,env,out)
    M=Main.MatsubaraAll5Refined;MM=Main.MatsubaraMoments
    blas1();base,binding=Main.source(sconf);check(c,out);blas1()
    Main.plain(binding)==TOML.parsefile(joinpath(path(c,"M3_source_directory"),"source-binding.toml")) || error("M3 actual/current model binding")
    emit(c,out,"source-binding.toml",binding);emit(c,out,"environment.toml",env)
    ops=M.qualify_native_ops(base);emit(c,out,"native-operator-controls.toml",(;ops,cp=M.cp_bound(base)))
    calls=Any[];maps=Ref(0);outer=Ref(0);final=Ref(0)
    function evaluate(mu,N)
        check(c,out);N==c["outer_N"] ? (outer[]+=1) : (final[]+=1)
        outer[]<=c["max_outer_calls"] && final[]<=1 || error("call budget")
        m=at_mu(base,mu);dir=joinpath(out,"call-$(length(calls)+1)-N$N");mkdir(dir)
        saveconf=copy(sconf);saveconf["mu_native"]=mu
        saveconf["data_max_bytes"]=c["data_max_bytes"]
        request=(;mu_native=mu,N,status="started_reserved",reserved_maps=16,alpha=1.0,S0="zero")
        push!(calls,request);emit(c,out,relpath(joinpath(dir,"request.toml"),out),request)
        emit(c,out,"calls.toml",(;calls,maps_completed_callbacks=maps[],reserved_maps=16*length(calls)))
        function callback(snap)
            check(c,out;extra=5*sizeof(snap.S)+262144);maps[]+=1
            maps[]<=c["max_inner_maps"] || error("map budget")
            it=snap.record.iteration
            Main.save_state(joinpath(dir,"i$(lpad(it,2,'0')).h5"),snap,saveconf)
            Main.save_state(joinpath(dir,"latest.h5"),snap,saveconf;full=true,replace=true)
            emit(c,out,relpath(joinpath(dir,"i$(lpad(it,2,'0')).toml"),out),
                (;raw=snap.record,raw_unweighted_scale=norm(snap.Phi),
                    weighted_raw=weighted_raw(snap,m),mu_native=mu,N))
        end
        remaining=c["wall_seconds"]-(time()-STARTED)
        blas1();result=M.pilot(m,N;steps=16,alpha=1.0,S0=nothing,soft_seconds=remaining,callback)
        ok=result.status=="iteration_budget" && qualified(result.history,c)
        calls[end]=(;request...,status=result.status,inner_qualified=ok,history=result.history)
        emit(c,out,"calls.toml",(;calls,maps_completed_callbacks=maps[],reserved_maps=16*length(calls)))
        ok || error("unqualified inner call; no density root value")
        emit(c,out,relpath(joinpath(dir,"M3-request.toml"),out),(;status="started_candidate_not_measured",N,mu_native=mu))
        cor=MM.correct(m,result.rho,N;check=()->check(c,out))
        terminal(c,out,dir,m,result,cor)
        rec=(;candidate_accepted=false,mu_native=mu,N,number=cor.number,raw_number=result.history[end].number,
            oneshot_number=cor.oneshot_number,linear_minus_oneshot=cor.number-cor.oneshot_number,
            linear_minus_oneshot_Frobenius=norm(cor.feedback),
            target_relative=(cor.number-c["target_native"])/c["target_native"],
            pauli_final=cor.pauli,oneshot_pauli_final=cor.oneshot_pauli,eta=cor.eta,
            linear_condition=cor.linear_condition,linear_residual_absolute=cor.linear_residual_absolute,
            linear_residual_scale=cor.linear_residual_scale,
            linear_residual_relative=cor.linear_residual_relative,
            linear_weighted_absolute=cor.linear_weighted_absolute,linear_weighted_scale=cor.linear_weighted_scale,
            linear_weighted_relative=cor.linear_weighted_relative,
            conditional_linear_error_bound=cor.conditional_linear_error_bound,tails=cor.tails,
            scope="asymptotic linear M3 closure; no infinite-N error certificate")
        emit(c,out,relpath(joinpath(dir,"density.toml"),out),rec)
        all(isfinite,cor.rho) && all(isfinite,cor.oneshot) || error("nonfinite M3 density")
        all(isfinite,(cor.number,cor.oneshot_number,cor.eta,cor.linear_condition,
            cor.linear_residual_absolute,cor.linear_residual_relative)) || error("nonfinite M3 diagnostics")
        cor.eta<1 && cor.linear_residual_relative<=c["M3_linear_relative_tolerance"] || error("M3 linear qualification")
        emit(c,out,relpath(joinpath(dir,"density-qualification.toml"),out),(;candidate_accepted=true,scope="scalar root proposal only"))
        return rec
    end
    outcome="failed";detail="";selected=nothing;witness=nothing
    try
        a=evaluate(base.mu-c["bracket_halfwidth_native"],c["outer_N"])
        b=evaluate(base.mu+c["bracket_halfwidth_native"],c["outer_N"])
        tol=c["charge_relative_tolerance"]
        if abs(a.target_relative)<=tol;selected=a
        elseif abs(b.target_relative)<=tol;selected=b
        else
            signbit(a.target_relative)!=signbit(b.target_relative) || error("no sign bracket; no expansion")
            while outer[]<c["max_outer_calls"]
                mu=proposal(a,b)
                a.mu_native<mu<b.mu_native || error("bracket exhausted floating spacing")
                r=evaluate(mu,c["outer_N"])
                if abs(r.target_relative)<=tol;selected=r;break;end
                signbit(r.target_relative)==signbit(a.target_relative) ? (a=r) : (b=r)
            end
        end
        selected===nothing && error("outer call budget; charge target not reached")
        witness=evaluate(selected.mu_native,c["check_N"])
        abs(witness.target_relative)<=c["witness_charge_relative_tolerance"] || error("independent witness charge target failed")
        abs(witness.number-selected.number)/c["target_native"]<=c["cutoff_relative_tolerance"] || error("outer/witness number cutoff mismatch")
        outcome="qualified_fixedcharge_M3_variant";detail="measured finite-N number cutoff witness only"
    catch err
        detail=sprint(showerror,err)
    end
    selected_record=selected===nothing ? Dict{String,Any}() : selected
    witness_record=witness===nothing ? Dict{String,Any}() : witness
    emit(c,out,"summary.toml",(;identity=c["identity"],outcome,detail,selected=selected_record,witness=witness_record,calls,
        maps_completed_callbacks=maps[],reserved_maps=16*length(calls),structural_map_bound=144,
        outer_calls=outer[],final_calls=final[],elapsed_seconds=time()-STARTED,BLAS_threads_final=BLAS.get_num_threads(),
        native_pass=false,infinite_N_certified=false,scientific_acceptance=false);final=true)
end
function run(input,out)
    blas1();c=TOML.parsefile(input);sconf,env=preflight(c)
    ispath(out) && error("output exists");mkdir(out)
    Base.include(Main,path(c,"source_adapter"));Base.invokelatest(() -> Main.load_modules(sconf))
    Base.include(Main,path(c,"refined_module"));Base.include(Main,path(c,"M3_module"))
    for k in ("source_input","source_adapter","source_module","refined_module","M3_module","M3_review","M3_qualification","M3_qualification_input")
        cp(path(c,k),joinpath(out,k*"-"*basename(path(c,k))))
    end
    cp(input,joinpath(out,"input.toml"));cp(@__FILE__,joinpath(out,"run_mu_root.jl"))
    Base.invokelatest(run_loaded,c,sconf,env,out)
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("usage: run_mu_root.jl authorized-input.toml new-output")
    BoundedMuRootRefined.run(ARGS[1],abspath(ARGS[2]))
end
