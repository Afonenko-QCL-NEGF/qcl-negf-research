# WRITE-ONLY; --tiny is source-free. ROOT owns all runtime/admission.
using QCLNEGF, Unitful, LinearAlgebra, HDF5, TOML, SHA, Test
const PARENT=normpath(joinpath(@__DIR__,"../equilibrium-gamma-hermitian-cont128/attempt01"))
include(joinpath(PARENT,"gamma_functions.jl"))
const GP=EquilibriumGammaPilot
const N=QCLNEGF.QCLNumerics
include("operators.jl")
sha_file(p)=bytes2hex(open(sha256,p))
function resource_gate(cfg)
    cfg["runtime_authorized"]===true || error("write-only; ROOT authorization required")
    cfg["history"]==0 && cfg["Pi_start_cap"]==2 && cfg["workers"]==1 && cfg["BLAS_threads"]==1 || error("fixed map envelope")
    cfg["fine_NE"]==61407 && cfg["fine_LO_intervals"]==1878 || error("refinement identity")
    cfg["memory_max_bytes"]==13019119616 && 0<cfg["memory_high_bytes"]<=cfg["memory_max_bytes"] || error("hard memory identity")
    0<cfg["soft_seconds"]<=840 && cfg["external_seconds"]==900 || error("time envelope")
    0<cfg["single_file_max_bytes"]<=2*1024^3 && 0<cfg["data_max_bytes"]<=3*1024^3 || error("data envelope")
    Threads.nthreads()==1 || error("Julia1 required"); BLAS.set_num_threads(1)
    root=dirname(dirname(pathof(QCLNEGF)))
    readchomp(`git -C $root rev-parse HEAD`)==cfg["core_head"] || error("Core HEAD")
    isempty(readchomp(`git -C $root status --porcelain --untracked-files=no`)) || error("dirty Core")
    entry=only(filter(s->startswith(s,"0::"),readlines("/proc/self/cgroup")))
    cg=joinpath("/sys/fs/cgroup",lstrip(split(entry,':';limit=3)[3],'/'))
    for (field,key) in (("memory.max","memory_max_bytes"),("memory.high","memory_high_bytes"),("memory.swap.max","swap_max_bytes"))
        strip(read(joinpath(cg,field),String))==string(cfg[key]) || error("external cap mismatch "*field)
    end
    cpu=split(read(joinpath(cg,"cpu.max"),String)); cpu[1]!="max" && parse(Int,cpu[1])<=parse(Int,cpu[2]) || error("CPU cap exceeds100%")
    limits=split(only(filter(s->startswith(s,"Max file size"),readlines("/proc/self/limits"))))
    limits[4]!="unlimited" && parse(Int,limits[4])<=cfg["single_file_max_bytes"] || error("external file cap")
    return Dict("cgroup"=>cg,"core_import_path"=>pathof(QCLNEGF),"core_head"=>cfg["core_head"],
        "Julia_version"=>string(VERSION),"memory_max"=>cfg["memory_max_bytes"],"native_stage"=>"absent")
end
function defect_and_samples(Phi,G,c,coarse,out,cfg)
    g=c.problem.grids; d=GP.WarmNorm(); n=(size(Phi,1)+1)÷2; source_energy=read(coarse["energy"])
    selected=Int[]
    for a in cfg["anchor_indices"],m in -40:40,s in -3:3,j in -1:1
        e=a+939m+s; q=2e-1+j
        1<=e<=n && 1<=q<=size(Phi,1) && push!(selected,q)
    end
    sort!(unique!(selected))
    values=Set(selected); mid_n=0.0; mid_s=0.0; mid_w=0.0; mid_ws=0.0
    open(joinpath(out,"selected-defect.csv"),"w") do io
        println(io,"e,k,E_eV,raw_r2,collocation_defect_r2,Phi_reference_r2,Gamma_trace_perNb_eV")
        for firste in 1:1024:size(Phi,1)
            rows=firste:min(firste+1023,size(Phi,1)); ref=prolong_rows(coarse,"Phi",rows,source_energy,g.ε)
            GP.accumulate_norm!(d,ref,Phi[rows,:,:,:],g.wᴱ,g.wᵏ,firste-1)
            for (i,e) in enumerate(rows)
                if iseven(e)
                    for k in axes(Phi,2)
                        bn=sum(abs2,Phi[e,k,:,:]-ref[i,k,:,:]); bs=sum(abs2,ref[i,k,:,:]); w=g.wᴱ[e]*g.wᵏ[k]
                        mid_n+=bn; mid_s+=bs; mid_w+=w*bn; mid_ws+=w*bs
                    end
                end
                e in values || continue
                for k in axes(Phi,2)
                    R=Phi[e,k,:,:]-G[e,k,:,:]; D=Phi[e,k,:,:]-ref[i,k,:,:]
                    println(io,join((e,k,g.ε[e]*.1,sum(abs2,R),sum(abs2,D),sum(abs2,ref[i,k,:,:]),real(tr(G[e,k,:,:]))/5*.1),','))
                end
            end
        end
    end
    return Dict("full_fine_defect"=>GP.warm_norm_dict(d),"selected_fine_indices"=>selected,
        "all_new_midpoint_defect"=>Dict("absolute"=>sqrt(mid_n),"scale"=>sqrt(mid_s),"relative"=>GP.relative(sqrt(mid_n),sqrt(mid_s)),
            "weighted_absolute"=>sqrt(mid_w),"weighted_scale"=>sqrt(mid_ws),"weighted_relative"=>GP.relative(sqrt(mid_w),sqrt(mid_ws))),
        "coarse_restriction_defect"=>compare_file(Phi,coarse,"Phi",read(coarse["weight_energy"]),g.wᵏ;restricted=true),
        "meaning"=>"raw Phi_f(PGamma)-P Phi_c; native criteria not changed")
end
function run_control(input,out)
    cfg=TOML.parsefile(input); identity=resource_gate(cfg)
    ispath(out) && error("new output directory required"); mkpath(out)
    for file in ("run_control.jl","operators.jl","tiny_checks.jl","protocol.md")
        cp(joinpath(@__DIR__,file),joinpath(out,file))
    end
    cp(input,joinpath(out,"input.toml")); started=Ref(0); completed=Ref(0); begun=time()
    summary=Dict{String,Any}("identity"=>identity,"execution"=>"running","SCBA_converged"=>false,
        "input_sha256"=>sha_file(input),"history_fields"=>0,"pi_start_cap"=>2,"PI_FFT_controls_outside_Pi"=>3)
    function checkpoint(stage)
        summary["stage"]=stage; summary["Pi_started"]=started[]; summary["Pi_completed"]=completed[]
        summary["seconds"]=time()-begun; summary["maxrss_bytes"]=Sys.maxrss()
        open(joinpath(out,"summary.toml"),"w") do io; TOML.print(io,summary); end
        GP.data_bytes(out)<=cfg["data_max_bytes"] || error("data cap")
        time()-begun<=cfg["soft_seconds"] || error("soft time cap")
        println("stage=",stage," Pi_started=",started[]," elapsed=",time()-begun); flush(stdout)
    end
    try
        parent=TOML.parsefile(joinpath(PARENT,"input.toml")); prior=TOML.parsefile(joinpath(PARENT,"map-128.toml"))
        string(VERSION)==parent["expected_julia_version"] || error("source Julia version")
        for (name,expected) in cfg["source_sha256"]
            sha_file(joinpath(PARENT,name))==expected || error("source SHA differs: "*name)
        end
        sourceid=TOML.parsefile(joinpath(PARENT,"identity.toml"))
        for name in ("gamma_functions.jl","warm_seed.jl","gamma_restart.jl")
            sha_file(joinpath(PARENT,name))==sourceid["source_sha256"][name] || error("included source identity")
        end
        summary["source_code_sha256"]=Dict(name=>sha_file(joinpath(@__DIR__,name)) for name in ("run_control.jl","operators.jl","tiny_checks.jl"))
        summary["active_project"]=Base.active_project(); summary["manifest_sha256"]=sha_file(joinpath(dirname(Base.active_project()),"Manifest.toml"))
        parent["NE"]==30704 && parent["energy_max_eV"]==.9000001064962726 || error("source grid")
        for (k,v) in (("workers",1),("memory_max_bytes",cfg["memory_max_bytes"])); parent[k]=v; end
        checkpoint("source_identity")
        c=GP.build_context(parent) # public native-cache builder ONLY on admitted coarse grid
        h5open(joinpath(PARENT,"model.h5"),"r") do f; summary["coarse_model_binding"]=GP.verify_problem_model(c,f); end
        G=h5open(joinpath(PARENT,"reduced-terminal.h5"),"r") do f
            size(f["Gamma_real"])==(30704,11,5,5) || error("source axes")
            read(f["mu"])==prior["mu"]["mu"] || error("source mu identity")
            GP.read_complex(f,"Gamma")
        end
        mu0=prior["mu"]["mu"]; GP.cone(G)
        private,cacheledger=gamma_cache(c.problem,c.options)
        L,equality=cache_equal(c.cache,private,c.problem,G)
        summary["coarse_cache_admission"]=cacheledger; summary["coarse_cache_equivalence"]=equality
        ghostbound=log(2)/(2pi)*(norm(view(G,1,:,:,:))+norm(view(G,size(G,1),:,:,:)))
        summary["ghost_model_change_bound_native"]=ghostbound
        private=nothing; GC.gc(true); checkpoint("coarse_cache_identical")
        Phi,stats=one_map(c,G,mu0,started,checkpoint); completed[]+=1
        replay=h5open(joinpath(PARENT,"reduced-terminal.h5"),"r") do f
            compare_file(Phi,f,"raw_Phi",c.problem.grids.wᴱ,c.problem.grids.wᵏ)
        end
        summary["coarse_replay"]=replay; summary["coarse_map"]=stats
        replay["relative"]<=1e-10 && replay["weighted_relative"]<=1e-10 || error("fresh coarse Phi replay failed")
        save_fields(joinpath(out,"coarse.h5"),cfg,("Gamma"=>G,"Lambda"=>L,"Phi"=>Phi);energy=c.problem.grids.ε,weights=c.problem.grids.wᴱ)
        G=Phi=L=nothing; GC.gc(true); checkpoint("coarse_saved")
        source_energy=copy(c.problem.grids.ε)
        fine,ledger=refined_context(c); summary["fine_cache_admission"]=ledger
        c=nothing; GC.gc(true)
        G=h5open(joinpath(PARENT,"reduced-terminal.h5"),"r") do f; prolong_field(f,"Gamma",source_energy,fine.problem.grids.ε); end
        GP.cone(G); summary["fine_Gamma_sha256"]=GP.array_digest(G)
        L,roundoff=N._production_pi_fft_hilbert_transform(G,fine.problem.grids.ε,fine.problem.grids.wᴱ,fine.cache.product_hilbert_operator)
        geometry=h5open(joinpath(out,"coarse.h5"),"r") do f
            compare_file(L,f,"Lambda",read(f["weight_energy"]),fine.problem.grids.wᵏ;restricted=true)
        end
        summary["HPI_source_node_difference"]=geometry; summary["fine_HPI_control_roundoff"]=roundoff
        geometry["relative"]<=1e-10 || error("HPI refinement arithmetic control failed")
        L=nothing; GC.gc(true); checkpoint("fine_HPI_control")
        save_fields(joinpath(out,"fine-input.h5"),cfg,("Gamma"=>G,);energy=fine.problem.grids.ε,weights=fine.problem.grids.wᴱ)
        Phi,stats=one_map(fine,G,mu0,started,checkpoint); completed[]+=1; summary["fine_map"]=stats
        save_fields(joinpath(out,"fine.h5"),cfg,("Phi"=>Phi,);energy=fine.problem.grids.ε,weights=fine.problem.grids.wᴱ)
        summary["defect"]=h5open(joinpath(out,"coarse.h5"),"r") do f; defect_and_samples(Phi,G,fine,f,out,cfg); end
        summary["execution"]="complete_frozen_collocation_control"; checkpoint("complete")
    catch e
        summary["execution"]="stopped"; summary["reason"]=sprint(showerror,e); rethrow()
    finally
        summary["Pi_started"]=started[]; summary["Pi_completed"]=completed[]; summary["maxrss_bytes"]=Sys.maxrss()
        summary["seconds"]=time()-begun
        open(joinpath(out,"summary.toml"),"w") do io; TOML.print(io,summary); end
    end
end
include("tiny_checks.jl")
if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS==["--tiny"]; tiny_checks()
    else
        length(ARGS)==2 || error("usage: run_control.jl --tiny | INPUT NEW_OUTPUT")
        run_control(ARGS[1],abspath(ARGS[2]))
    end
end
