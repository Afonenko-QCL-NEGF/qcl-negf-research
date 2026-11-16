# ROOT-only isolated kernel benchmark. No new Phi, solver/model build or Core edit.
module InfoPaddingBench
const STARTED=time()
using QCLNEGF,LinearAlgebra,HDF5,TOML,SHA
const N=QCLNEGF.QCLNumerics
const IDS=(:LO,:acoustic,:impurity,:IFR,:alloy,:embedding,:embedding_plus,:embedding_minus)
sha(p)=bytes2hex(open(sha256,p))
loc(c,k)=normpath(joinpath(@__DIR__,c[k]))
needsha(p,h)=(sha(p)==h || error("SHA mismatch: $p"))
plain(x::NamedTuple)=Dict(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::AbstractDict)=Dict(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::AbstractVector)=plain.(x)
plain(x::Symbol)=string(x)
plain(x::Nothing)="not_measured"
plain(x)=x
function write_report(out,r,c)
    b=sprint(io->TOML.print(io,plain(r)))
    other=sum(filesize(joinpath(out,f)) for f in readdir(out) if !startswith(f,"summary.toml");init=0)
    other+2sizeof(b)<=c["output_max_bytes"] || error("metadata output cap including atomic replace")
    write(joinpath(out,"summary.toml.tmp"),b)
    mv(joinpath(out,"summary.toml.tmp"),joinpath(out,"summary.toml");force=true)
end
function gate(c)
    f=TOML.parsefile(joinpath(@__DIR__,"input.toml"))
    all(k=="runtime_authorized" || isequal(c[k],v) for (k,v) in f) || error("frozen settings")
    c["runtime_authorized"]===true || error("runtime_authorized=false")
    needsha(@__FILE__,c["script_sha256"])
    BLAS.set_num_threads(1);BLAS.get_num_threads()==1 || error("BLAS1")
    string(VERSION)==c["julia_version"] && Threads.nthreads(:default)==26 || error("Julia/J26")
    get(ENV,"OMP_NUM_THREADS","")=="1" || error("OMP1")
    core=dirname(dirname(pathof(QCLNEGF)))
    realpath(core)==realpath(loc(c,"core")) || error("Core path")
    readchomp(`git -C $core rev-parse HEAD`)==c["core_head"] || error("Core HEAD")
    readchomp(`git -C $core rev-parse $("HEAD^{tree}")`)==c["core_tree"] || error("Core tree")
    isempty(readchomp(`git -C $core status --porcelain --untracked-files=no`)) || error("dirty Core")
    needsha(Base.active_project(),c["project_sha256"])
    needsha(joinpath(dirname(Base.active_project()),"Manifest.toml"),c["manifest_sha256"])
    member=only(filter(x->startswith(x,"0::"),readlines("/proc/self/cgroup")))
    cg=joinpath("/sys/fs/cgroup",lstrip(split(member,':';limit=3)[3],'/'))
    for (file,key) in (("memory.max","memory_max_bytes"),("memory.swap.max","swap_max_bytes"))
        strip(read(joinpath(cg,file),String))==string(c[key]) || error("OS $file")
    end
    cpu=split(strip(read(joinpath(cg,"cpu.max"),String)))
    cpu[1]!="max" && parse(Int,cpu[1])==26*parse(Int,cpu[2]) || error("CPU2600 required")
    lim=split(only(filter(x->startswith(x,"Max file size"),readlines("/proc/self/limits"))))
    lim[4]!="unlimited" && parse(Int,lim[4])<=c["output_max_bytes"] || error("OS file cap")
    for (file,h) in c["helper_sha256"];needsha(joinpath(loc(c,"helper_directory"),file),h);end
    for (file,h) in c["source_sha256"];needsha(joinpath(loc(c,"parent"),file),h);end
    return (;core,project=Base.active_project(),cgroup=cg,BLAS_threads=BLAS.get_num_threads())
end
function readfield(GP,g,name,shape)
    size(g[name*"_real"])==size(g[name*"_imag"])==shape || error("field axes: $name")
    a=GP.read_complex(g,name)
    a isa Array{ComplexF64} && size(a)==shape && all(isfinite,a) || error("field finite/type")
    GC.gc(false);return a
end
function load_fixture(c,GP)
    parent=loc(c,"parent");id=TOML.parsefile(joinpath(parent,"identity.toml"))
    sm=TOML.parsefile(joinpath(parent,"summary.toml"));pc=TOML.parsefile(joinpath(parent,"input.toml"))
    id["core_head"]==c["core_head"] && id["manifest_sha256"]==c["manifest_sha256"] || error("parent environment")
    id["input_sha256"]==c["source_sha256"]["input.toml"] && id["script_sha256"]==c["source_sha256"]["run_native.jl"] || error("parent lineage")
    sm["execution"]=="complete" && sm["maps_completed"]==128 && pc["voltage_mV"]==0 || error("parent completed case0")
    sm["state_sha256"]==c["source_sha256"]["state.h5"] && sm["model_sha256"]==c["source_sha256"]["model.h5"] || error("terminal ownership")
    shape=(24001,11,5,5)
    energy,h=h5open(joinpath(parent,"model.h5"),"r") do f
        read(f["schema"])=="campaign133-gamma-model-v2" || error("model schema")
        read(f["grid_ε"]),readfield(GP,f,"h_fixedU0",(11,5,5))
    end
    length(energy)==shape[1] && all(isfinite,energy) && issorted(energy) || error("energy")
    green,families=h5open(joinpath(parent,"state.h5"),"r") do f
        read(f["schema"])=="campaign133-native-bias-fullstate-v1" && read(f["voltage_per_period_mV"])==0 || error("state schema/bias")
        isequal(read(f["energy"]),energy) && all(iszero,read(f["hartree"])) || error("state energy/U0")
        read(f["energy_scale_eV"])==0.1 && read(f["length_scale_m"])==1e-8 && read(f["spin_degeneracy"])==2 || error("native units")
        g=f["green"];a=[readfield(GP,g,s,shape) for s in ("GR","GL","GG","A")]
        condition=read(g["condition_number"]);scale=read(g["dyson_scale"])
        size(condition)==size(scale)==shape[1:2] && all(isfinite,condition) && all(isfinite,scale) || error("Green diagnostics")
        gs=GreenState(a...,condition,scale)
        Set(keys(f["selfenergy"]))==Set(string.(IDS)) || error("eight source families")
        fs=[SelfEnergyFamily((readfield(GP,f["selfenergy/"*string(s)],n,shape)
            for n in ("retarded","lesser","greater"))...) for s in IDS]
        gs,fs
    end
    # Deliberate SURROGATE: all8 includes embedding sum AND both sides; not physical total.
    total=SelfEnergyFamily(zeros(ComplexF64,shape),zeros(ComplexF64,shape),zeros(ComplexF64,shape))
    for family in families, name in (:Σᴿ,:Σˡ,:Σᵍ);getfield(total,name).+=getfield(family,name);end
    return (;energy,h,green,families,total,candidate=total)
end
function line_size()
    base="/sys/devices/system/cpu/cpu0/cache"
    for d in readdir(base;join=true)
        isfile(joinpath(d,"type")) || continue
        strip(read(joinpath(d,"type"),String))=="Data" || continue
        strip(read(joinpath(d,"level"),String))=="1" || continue
        return parse(Int,strip(read(joinpath(d,"coherency_line_size"),String)))
    end
    error("L1 data coherency line size unavailable")
end
addr(a)=UInt(pointer(a)) # rooted vectors remain owned by workspace/layout tables.
function layout_metadata(info,line)
    [(;worker=w,length=length(a),pointer=string(addr(a)),line=string(addr(a)÷UInt(line)),
        offset=Int(addr(a)%UInt(line))) for (w,a) in enumerate(info)]
end
const FIXED_BUFFERS=(:spectral_numerator,:spectral_denominator,:keldysh_numerator,
    :keldysh_denominator,:dyson_partial,:positivity_partial,:positivity_witness,:causality_partial)
function fixed_addresses(ws)
    v=UInt[addr(getfield(ws,f)) for f in FIXED_BUFFERS]
    for f in (:D,:Γ,:temporary_1,:temporary_2,:hermitian,:lapack_values,:lapack_work,:lapack_rwork)
        append!(v,addr.(getfield(ws,f)))
    end
    return v
end
swap_info!(ws,layout)=(foreach(w->ws.lapack_info[w]=layout[w],eachindex(layout));nothing)
function snapshot(metrics,ws)
    w=N._production_positivity_witness(ws)
    (;metrics,witness=NamedTuple{fieldnames(typeof(w))}(Tuple(getfield(w,f) for f in fieldnames(typeof(w)))))
end
function fingerprints(GP,f)
    a=[f.energy,f.h,f.green.Gᴿ,f.green.Gˡ,f.green.Gᵍ,f.green.A,f.green.condition_number,f.green.dyson_scale]
    for s in [f.families;[f.total]], n in (:Σᴿ,:Σˡ,:Σᵍ);push!(a,getfield(s,n));end
    GP.array_digest.(a)
end
function run_loaded(c,out,env)
    GP=Main.EquilibriumGammaPilot;rows=Any[]
    report=Dict{String,Any}("execution"=>"started","scope"=>"residual suite ONLY; sum8 surrogate and candidate===total; NO Phi/SCBA",
        "scientific_acceptance"=>false,"false_sharing_diagnosed"=>false,"environment"=>env,"rows"=>rows)
    try
        f=load_fixture(c,GP);ws=ProductionResidualWorkspace(24001,11,5;chunk_size=64,worker_count=26)
        original=copy(ws.lapack_info);all(length(a)==1 for a in original) || error("baseline INFO size")
        sentinel=LinearAlgebra.BlasInt(0x13579);padded=[fill(sentinel,8) for _ in 1:26]
        foreach(a->a[1]=0,padded);line=line_size();fixed=fixed_addresses(ws);before=fingerprints(GP,f)
        layouts=Dict("A"=>original,"B"=>padded)
        addresses=Dict(k=>layout_metadata(v,line) for (k,v) in layouts)
        report["layout"]=Dict("BlasInt_bytes"=>sizeof(LinearAlgebra.BlasInt),"cache_line_bytes"=>line,"info"=>addresses,
            "unchanged_buffers_pointer_digest"=>bytes2hex(sha256(reinterpret(UInt8,fixed))))
        report["fixture"]=Dict("shape"=>[24001,11,5,5],"families_order"=>collect(string.(IDS)),"candidate_same_object"=>f.candidate===f.total,
            "parent_converged"=>false,"initial_logical_digests"=>before,"workspace_capacity"=>26,"chunk_size"=>64)
        report["setup_seconds"]=time()-STARTED;write_report(out,report,c)
        reference=Ref{Any}(nothing);computed=time();worst=Dict(4=>0.0,26=>0.0)
        function call(W,key,repeat,position,warmup)
            time()-computed+max(5.0,1.5*worst[W])<c["compute_soft_seconds"] || error("compute soft budget: no next kernel admitted")
            swap_info!(ws,layouts[key]);fixed_addresses(ws)==fixed || error("nonINFO buffer moved")
            a=N._production_counter_sample()
            metrics=production_residual_suite(f.energy,f.h,f.green,f.total,f.candidate;workspace=ws,worker_count=W)
            b=N._production_counter_sample();wall=(b.monotonic_ns-a.monotonic_ns)*1e-9;worst[W]=max(worst[W],wall)
            answer=snapshot(metrics,ws);reference[]===nothing && (reference[]=answer)
            equal=isequal(answer,reference[])
            tails=all(all(==(sentinel),@view v[2:8]) for v in padded)
            stable=fixed_addresses(ws)==fixed && all(layout_metadata(layouts[k],line)==addresses[k] for k in keys(layouts))
            row=(;workers=W,layout=key,repeat,position,warmup,wall_seconds=wall,cpu_seconds=b.cpu_seconds-a.cpu_seconds,
                allocated_bytes=b.allocated_bytes-a.allocated_bytes,gc_seconds=b.gc_seconds-a.gc_seconds,gc_count=b.gc_count-a.gc_count,
                outputs=answer,exact_reference_equal=equal,padding_tail_unchanged=tails,buffer_addresses_unchanged=stable)
            push!(rows,row);report["compute_elapsed_seconds"]=time()-computed;write_report(out,report,c)
            equal && tails && stable || error("strict output/witness/padding/address qualification failed")
        end
        # Same specialization; both widths warm up, excluded from timing comparisons.
        call(26,"A",0,0,true);call(4,"A",0,0,true)
        for W in (4,26), repeat in 1:2, (position,key) in enumerate(("A","B","A"))
            call(W,key,repeat,position,false)
        end
        report["compute_elapsed_seconds"]=time()-computed
        report["inputs_bitwise_unchanged"]=fingerprints(GP,f)==before
        report["inputs_bitwise_unchanged"] || error("fixture mutated")
        report["execution"]="completed_benchmark_only"
    catch err
        report["execution"]="stopped";report["reason"]=sprint(showerror,err);rethrow()
    finally
        report["elapsed_including_imports_seconds"]=time()-STARTED;report["maxrss_bytes"]=Sys.maxrss()
        report["BLAS_threads_final"]=BLAS.get_num_threads();write_report(out,report,c)
    end
end
function run(input,out)
    c=TOML.parsefile(input);env=gate(c);ispath(out) && error("new output required");mkdir(out)
    cp(input,joinpath(out,"input.toml"));cp(@__FILE__,joinpath(out,"bench.jl"))
    Base.include(Main,joinpath(loc(c,"helper_directory"),"gamma_functions.jl"))
    Base.invokelatest(run_loaded,c,out,env)
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("bench.jl AUTHORIZED_INPUT NEW_OUTPUT")
    InfoPaddingBench.run(abspath(ARGS[1]),abspath(ARGS[2]))
end
