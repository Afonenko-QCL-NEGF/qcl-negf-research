# WRITE ONLY. ROOT owns bounded runtime. No full Phi/Pi call in this driver.
using SHA, TOML
const REFINE=normpath(joinpath(@__DIR__,"../refinement-control"))
mu_sha(p)=bytes2hex(open(sha256,p))
const DEFINITION_SHA=Dict(
    "run_control.jl"=>"5084c1f1f97e6b90bb8a369e3a212d6e68f94f9930a2aaa467406728d7d484cb",
    "operators.jl"=>"7d71dadce03edc17bcd7c68bae2fd5d1057e17cf95aae66432859ce4caa4a4dd",
    "tiny_checks.jl"=>"565a8d182372c0d5cfab45f7598df10a0f593c14cfd39445140aa69438e33120")
for (name,digest) in DEFINITION_SHA
    mu_sha(joinpath(REFINE,name))==digest || error("frozen refinement definition changed: "*name)
end
const GP_PARENT=normpath(joinpath(REFINE,"../equilibrium-gamma-hermitian-cont128/attempt01"))
for (name,digest) in (("gamma_functions.jl","57f2c277255ae1f8265659e65176cc88f97d5a40cbf8d6e32b7ea9257ccfb032"),
                     ("warm_seed.jl","d786df94a73aaeebb1bda1a2633806f3a5fef578e389bd1524d464fd9a048644"),
                     ("gamma_restart.jl","b8d1f638c3bb30b22798a79767ed66e957dd2f27a01031aaffb29080ba7affb7"))
    mu_sha(joinpath(GP_PARENT,name))==digest || error("frozen Gamma definition changed: "*name)
end
include(joinpath(REFINE,"run_control.jl")) # definitions only; PROGRAM_FILE differs

function mu_resource_gate(cfg)
    cfg["runtime_authorized"]===true || error("write-only; ROOT runtime authorization required")
    cfg["memory_max_bytes"]==13019119616 && cfg["memory_high_bytes"]==12348030976 && cfg["swap_max_bytes"]==0 || error("frozen admission/cap identity")
    cfg["workers"]==cfg["BLAS_threads"]==1 && cfg["retarded_cap"]==1 && cfg["LO_cap"]==2 && cfg["Pi_cap"]==0 || error("operator caps")
    Threads.nthreads()==1 || error("Julia1 required"); BLAS.set_num_threads(1)
    cfg["external_seconds"]==300 && 0<cfg["soft_seconds"]<=290 || error("wall envelope")
    cfg["single_file_max_bytes"]==32*1024^2 && cfg["data_max_bytes"]==64*1024^2 || error("output envelope")
    root=dirname(dirname(pathof(QCLNEGF)))
    readchomp(`git -C $root rev-parse HEAD`)==cfg["core_head"] || error("Core HEAD")
    isempty(readchomp(`git -C $root status --porcelain --untracked-files=no`)) || error("dirty Core")
    entry=only(filter(s->startswith(s,"0::"),readlines("/proc/self/cgroup")))
    cg=joinpath("/sys/fs/cgroup",lstrip(split(entry,':';limit=3)[3],'/'))
    for (field,key) in (("memory.max","memory_max_bytes"),("memory.high","memory_high_bytes"),("memory.swap.max","swap_max_bytes"))
        strip(read(joinpath(cg,field),String))==string(cfg[key]) || error("external cap: "*field)
    end
    cpu=split(read(joinpath(cg,"cpu.max"),String)); cpu[1]!="max" && parse(Int,cpu[1])<=parse(Int,cpu[2]) || error("CPU exceeds100%")
    limits=split(only(filter(s->startswith(s,"Max file size"),readlines("/proc/self/limits"))))
    limits[4]!="unlimited" && parse(Int,limits[4])<=cfg["single_file_max_bytes"] || error("external file cap")
    return Dict("cgroup"=>cg,"Core_source"=>pathof(QCLNEGF),"core_head"=>cfg["core_head"],"Julia"=>string(VERSION),
        "Julia_threads"=>Threads.nthreads(),"BLAS_threads"=>BLAS.get_num_threads(),"native_calls"=>0)
end

# Fixed slots: D,mu,fixed,reference,R,Rfixed squared+weighted (1:12),
# then D⋅mu,D⋅fixed,mu⋅fixed inner+weighted (13:18). No hot string keys.
function add_mu_moments!(s::Matrix{Float64},col::Int,D,M,P,R,w::Float64)
    F=D-M; RF=R-M
    d=abs2(D); m=abs2(M); f=abs2(F); p=abs2(P); r=abs2(R); rf=abs2(RF)
    dm=real(conj(D)*M); df=real(conj(D)*F); mf=real(conj(M)*F)
    @inbounds begin
        s[1,col]+=d; s[2,col]+=w*d
        s[3,col]+=m; s[4,col]+=w*m
        s[5,col]+=f; s[6,col]+=w*f
        s[7,col]+=p; s[8,col]+=w*p
        s[9,col]+=r; s[10,col]+=w*r
        s[11,col]+=rf; s[12,col]+=w*rf
        s[13,col]+=dm; s[14,col]+=w*dm
        s[15,col]+=df; s[16,col]+=w*df
        s[17,col]+=mf; s[18,col]+=w*mf
    end
    return nothing
end
empty_mu_moments()=zeros(Float64,18,3) # full, old nodes, midpoints
function finish_mu_moments(s::Matrix{Float64},col::Int=1)
    d=Dict{String,Any}()
    names=("D","mu","fixed_mu","reference","fine_residual","fine_residual_fixed_mu")
    for (i,name) in enumerate(names)
        d[name*"_squared"]=s[2*i-1,col]; d[name*"_weighted_squared"]=s[2*i,col]
    end
    for (i,name) in enumerate(("D_mu","D_fixed","mu_fixed"))
        d[name*"_inner"]=s[11+2*i,col]; d[name*"_weighted_inner"]=s[12+2*i,col]
    end
    for name in ("D","mu","fixed_mu","fine_residual","fine_residual_fixed_mu")
        d[name*"_absolute"]=sqrt(d[name*"_squared"])
        d[name*"_weighted_absolute"]=sqrt(d[name*"_weighted_squared"])
        d[name*"_relative_to_PPhi"]=GP.relative(d[name*"_absolute"],sqrt(d["reference_squared"]))
        d[name*"_weighted_relative_to_PPhi"]=GP.relative(d[name*"_weighted_absolute"],sqrt(d["reference_weighted_squared"]))
    end
    d["squared_identity_error"]=d["D_squared"]-d["mu_squared"]-d["fixed_mu_squared"]-2*d["mu_fixed_inner"]
    d["weighted_squared_identity_error"]=d["D_weighted_squared"]-d["mu_weighted_squared"]-d["fixed_mu_weighted_squared"]-2*d["mu_fixed_weighted_inner"]
    d["mu_projection_coefficient"]=GP.relative(d["D_mu_inner"],d["D_squared"])
    d["mu_projection_coefficient_weighted"]=GP.relative(d["D_mu_weighted_inner"],d["D_weighted_squared"])
    return d
end
function lo_at_mu(c,GR,A,condition,scale,mu,calls,checkpoint)
    calls[]<2 || error("two LO cap"); calls[]+=1; checkpoint("LO_started")
    green,pair=GP.fd_green(GR,A,condition,scale,c.problem.grids.ε,mu,c.kBT)
    lesser,greater=N._production_lo_contraction(c.cache.kernels[:LO],green,c.problem,c.cache,c.problem.kernels.qᴷ[:LO],c.options)
    @inbounds @simd for i in eachindex(lesser); lesser[i]=im*(greater[i]-lesser[i]); end
    all(isfinite,lesser) || error("nonfinite LO Gamma")
    return lesser,Dict("mu"=>mu,"Gamma_norm"=>norm(lesser),"FD_pair_error"=>pair,"pair_A_error"=>GP.pair_A_error(green))
end
function stream_mu_decomposition(c,delta,finefile,coarsefile,gammafile,out,cfg)
    g=c.problem.grids; n=length(g.ε); source_energy=read(coarsefile["energy"])
    length(source_energy)*2-1==n && g.ε[1:2:end]==source_energy || error("prolong geometry")
    stats=empty_mu_moments()
    selected=Set{Int}()
    for anchor in cfg["anchor_indices"], sideband in -40:40, offset in -3:3, nearby in -1:1
        old=anchor+939*sideband+offset; fineidx=2*old-1+nearby
        1<=old<=length(source_energy) && 1<=fineidx<=n && push!(selected,fineidx)
    end
    open(joinpath(out,"selected-LO-mu.csv"),"w") do io
        println(io,"e,k,E_eV,D_squared,mu_squared,fixed_mu_squared,D_mu_real_inner,mu_fixed_real_inner,mu_trace_real,mu_trace_imag")
        for start in 1:512:n
            rows=start:min(start+511,n)
            fine=GP.read_complex_rows(finefile,"Phi",rows)
            P=prolong_rows(coarsefile,"Phi",rows,source_energy,g.ε)
            gamma=GP.read_complex_rows(gammafile,"Gamma",rows)
            all(isfinite,fine) && all(isfinite,P) && all(isfinite,gamma) || error("nonfinite saved fine/coarse/Gamma chunk")
            @inbounds for (r,e) in enumerate(rows),m in axes(fine,2)
                ds=ms=fs=dm=mf=0.0; trmu=0.0+0.0im
                for b in axes(fine,4),a in axes(fine,3)
                    D=fine[r,m,a,b]-P[r,m,a,b]; M=delta[e,m,a,b]; F=D-M; R=fine[r,m,a,b]-gamma[r,m,a,b]
                    w=g.wᴱ[e]*g.wᵏ[m]
                    add_mu_moments!(stats,1,D,M,P[r,m,a,b],R,w)
                    add_mu_moments!(stats,isodd(e) ? 2 : 3,D,M,P[r,m,a,b],R,w)
                    ds+=abs2(D); ms+=abs2(M); fs+=abs2(F); dm+=real(conj(D)*M); mf+=real(conj(M)*F)
                    a==b && (trmu+=M)
                end
                e in selected && println(io,join((e,m,g.ε[e]*c.problem.scales.E₀_eV,ds,ms,fs,dm,mf,real(trmu),imag(trmu)),','))
            end
        end
    end
    return Dict("moments"=>Dict(k=>finish_mu_moments(stats,i) for (i,k) in enumerate(("full","old_nodes","midpoints"))),"selected_indices"=>sort!(collect(selected)),
        "definitions"=>"D=Phi_f(saved,mu_f)-P Phi_c(saved,mu_c); M=LO_f(mu_f)-LO_f(mu_c); F=D-M; R=Phi_f-Gamma_f; no floors/clipping; weights=wE_f*wk")
end
function run_mu(input,out)
    cfg=TOML.parsefile(input); id=mu_resource_gate(cfg)
    ispath(out) && error("new output directory required"); mkpath(out)
    begun=time(); Rcalls=Ref(0); LOcalls=Ref(0)
    summary=Dict{String,Any}("identity"=>id,"execution"=>"running","Pi_calls"=>0,"full_Phi_calls"=>0,"input_sha256"=>mu_sha(input))
    function checkpoint(stage)
        summary["stage"]=stage; summary["retarded_started"]=Rcalls[]; summary["LO_started"]=LOcalls[]
        summary["seconds"]=time()-begun; summary["maxrss_bytes"]=Sys.maxrss()
        open(joinpath(out,"summary.toml"),"w") do io; TOML.print(io,summary); end
        GP.data_bytes(out)<=cfg["data_max_bytes"] || error("data cap")
        time()-begun<=cfg["soft_seconds"] || error("soft wall cap")
        println("stage=",stage," R=",Rcalls[]," LO=",LOcalls[]); flush(stdout)
    end
    try
        for file in ("run_mu.jl","tiny_checks.jl","protocol.md"); cp(joinpath(@__DIR__,file),joinpath(out,file)); end
        cp(input,joinpath(out,"input.toml"))
        files=Dict("coarse.h5"=>joinpath(REFINE,"attempt01/coarse.h5"),"fine-input.h5"=>joinpath(REFINE,"attempt01/fine-input.h5"),
            "fine.h5"=>joinpath(REFINE,"attempt01/fine.h5"),"summary.toml"=>joinpath(REFINE,"attempt01/summary.toml"),
            "model.h5"=>joinpath(PARENT,"model.h5"),"parent-input.toml"=>joinpath(PARENT,"input.toml"),"parent-identity.toml"=>joinpath(PARENT,"identity.toml"))
        # Every WHOLE source digest is verified before opening ANY source HDF5.
        hashes=Dict{String,String}()
        for (name,path) in files
            expected=cfg["source_sha256"][name]
            length(expected)==64 && expected!="UNSET" || error("ROOT must bind source SHA: "*name)
            hashes[name]=mu_sha(path); hashes[name]==expected || error("source SHA changed: "*name)
        end
        summary["source_sha256"]=hashes; summary["hash_seconds"]=time()-begun
        summary["code_sha256"]=Dict(name=>mu_sha(joinpath(@__DIR__,name)) for name in ("run_mu.jl","tiny_checks.jl","protocol.md"))
        summary["definition_sha256"]=DEFINITION_SHA
        saved=TOML.parsefile(files["summary.toml"])
        saved["Pi_started"]==saved["Pi_completed"]==2 && saved["stage"]=="Pi_complete" || error("source maps incomplete")
        saved["source_code_sha256"]==DEFINITION_SHA || error("source definition identity")
        summary["manifest_sha256"]=mu_sha(joinpath(dirname(Base.active_project()),"Manifest.toml"))
        summary["manifest_sha256"]==saved["manifest_sha256"] || error("source environment manifest")
        parent=TOML.parsefile(files["parent-input.toml"])
        parent["memory_max_bytes"]==cfg["memory_max_bytes"] || error("source numerical/OS budget differs")
        string(VERSION)==parent["expected_julia_version"] || error("source Julia")
        checkpoint("all_sources_hashed")
        coarse=GP.build_context(parent)
        h5open(files["model.h5"],"r") do f; summary["model_binding"]=GP.verify_problem_model(coarse,f); end
        h5open(files["coarse.h5"],"r") do f
            read(f["energy"])==coarse.problem.grids.ε && read(f["weight_energy"])==coarse.problem.grids.wᴱ || error("coarse source axes/weights")
        end
        c,ledger=refined_context(coarse); summary["fine_cache_admission"]=ledger
        coarse=nothing; GC.gc(true)
        G=h5open(files["fine-input.h5"],"r") do f
            read(f["energy"])==c.problem.grids.ε && read(f["weight_energy"])==c.problem.grids.wᴱ || error("fine axes")
            GP.read_complex(f,"Gamma")
        end
        size(G)==(61407,11,5,5) && GP.array_digest(G)==saved["fine_Gamma_sha256"] || error("fine Gamma identity")
        for name in ("fine.h5","coarse.h5")
            h5open(files[name],"r") do f
                expected=name=="fine.h5" ? c.problem.grids.ε : c.problem.grids.ε[1:2:end]
                read(f["energy"])==expected || error("saved Phi axes")
                name=="fine.h5" && read(f["weight_energy"])!=c.problem.grids.wᴱ && error("saved fine Phi weights")
            end
        end
        Rcalls[]+=1; checkpoint("retarded_started")
        GR,condition,scale,roundoff=GP.retarded_state(c,G)
        G=nothing; GC.gc(true)
        A=spectral_function(GR)
        q=N._seed_number_weights_production(GR,c.problem.grids,c.problem.physical.g_s;options=c.options)
        muf=saved["fine_map"]["mu"]["mu"]; muc=saved["coarse_map"]["mu"]["mu"]
        summary["mu_fine"]=muf; summary["mu_coarse"]=muc; summary["mu_difference"]=muf-muc
        checks=Dict{String,Any}()
        for (name,mu,expected) in (("fine",muf,saved["fine_map"]["mu"]["number"]),("coarse",muc,saved["fine_map"]["fixed_mu_diagnostic"]["N"]))
            number,derivative=GP.number_mu(c.problem.grids.ε,q,mu,c.kBT)
            err=abs(number-expected)/c.target; err<=1e-10 || error("fresh same-fine q/N replay: "*name)
            checks[name]=Dict("mu"=>mu,"N"=>number,"expected_N"=>expected,"absolute_N_error"=>abs(number-expected),"N_error_over_target"=>err,"dN_dmu"=>derivative)
        end
        summary["fresh_number_replay"]=checks; summary["PI_roundoff"]=roundoff; summary["capacity"]=sum(q)
        checkpoint("retarded_and_number_complete")
        delta,fstats=lo_at_mu(c,GR,A,condition,scale,muf,LOcalls,checkpoint)
        expected=saved["fine_map"]["channel_Gamma_norms"]["LO"]
        relative_error=abs(fstats["Gamma_norm"]-expected)/expected
        relative_error<=1e-10 || error("fresh native fine LO norm replay")
        fstats["saved_Gamma_norm"]=expected; fstats["norm_relative_error"]=relative_error
        other,cstats=lo_at_mu(c,GR,A,condition,scale,muc,LOcalls,checkpoint)
        @inbounds @simd for i in eachindex(delta); delta[i]-=other[i]; end
        other=GR=A=condition=scale=q=nothing; GC.gc(true)
        summary["fine_mu_LO"]=fstats; summary["coarse_mu_LO_on_same_fine_A"]=cstats
        checkpoint("two_LO_complete")
        summary["decomposition"]=h5open(files["fine.h5"],"r") do fine
            h5open(files["coarse.h5"],"r") do old
                h5open(files["fine-input.h5"],"r") do gamma
                    stream_mu_decomposition(c,delta,fine,old,gamma,out,cfg)
                end
            end
        end
        moments=summary["decomposition"]["moments"]["full"]
        for (actual,key) in (("fine_residual_absolute","absolute"),("fine_residual_weighted_absolute","weighted_absolute"))
            expected=saved["fine_map"]["raw_residual"][key]
            abs(moments[actual]-expected)/expected<=1e-10 || error("saved fine residual stream replay")
        end
        # Stream WHOLE hashes again; source outputs stay immutable.
        summary["sources_unchanged_after"]=all(mu_sha(path)==hashes[name] for (name,path) in files)
        summary["sources_unchanged_after"] || error("source mutated during diagnostic")
        summary["execution"]="complete_mu_attribution_only"; checkpoint("complete")
    catch e
        summary["execution"]="stopped"; summary["reason"]=sprint(showerror,e); rethrow()
    finally
        summary["seconds"]=time()-begun; summary["maxrss_bytes"]=Sys.maxrss()
        summary["retarded_started"]=Rcalls[]; summary["LO_started"]=LOcalls[]
        open(joinpath(out,"summary.toml"),"w") do io; TOML.print(io,summary); end
    end
end
include("tiny_checks.jl")
if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS==["--tiny-mu"]; tiny_mu_checks()
    else
        length(ARGS)==2 || error("usage: run_mu.jl --tiny-mu | INPUT NEW_OUTPUT")
        run_mu(ARGS[1],abspath(ARGS[2]))
    end
end
