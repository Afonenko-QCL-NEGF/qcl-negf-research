# WRITE-ONLY standalone bounded stage; ROOT owns static/tiny/runtime.
using QCLNEGF, LinearAlgebra, Unitful, HDF5, SHA, TOML, Test
const SOURCE=normpath(joinpath(@__DIR__,"../equilibrium-gamma-stableA-cont64/attempt01"))
const STREAM=normpath(joinpath(@__DIR__,"../streaming-spectral/stream_spectral.jl"))
const CACHE=normpath(joinpath(@__DIR__,"../refinement-control/attempt01/operators.jl"))
sha_file(p)=bytes2hex(open(sha256,p))
const DEFS=Dict("gamma_functions.jl"=>"50bcfaaeb3c399baaaa2a4b9b9e023167e97b4d87c6af95d83353cb453072211",
    "warm_seed.jl"=>"68fd8582334a1ed3f257efcc0b1096a429361dd73aa5be3be415dd79ddcbc157",
    "stable_spectral.jl"=>"eeca1f29710d9bb60bf6ce70ce0ebd1620fb86d4c77314b4bd8ffa9a5bcd5aa5",
    "gamma_restart.jl"=>"2a805be29aeef33d7a31356850b7fdaa80da7e68506d8f96b980a8108bd56a00")
for (name,digest) in DEFS;sha_file(joinpath(SOURCE,name))==digest || error("preinclude source SHA: "*name);end
sha_file(STREAM)=="8df7885cbf547ffd4077e611f6a88b933dc75c154053f29e420230de538db898" || error("preinclude stream SHA")
sha_file(CACHE)=="7d71dadce03edc17bcd7c68bae2fd5d1057e17cf95aae66432859ce4caa4a4dd" || error("preinclude cache SHA")
include(joinpath(SOURCE,"gamma_functions.jl")); const GP=EquilibriumGammaPilot; const N=QCLNEGF.QCLNumerics
include(STREAM); const SS=StreamingSpectral
include(CACHE) # exact owning-cache/prolongation helpers; old one_map/refined_context NEVER called
pick(n) = sort!(unique!([q for e in (1,2,382,525,6955,7894,7907,8040,8846,9785,10724,11663,19702,30704)
    for q in (n==30704 ? (e,) : (2 * e-2,2 * e-1,2 * e)) if 1<=q<=n]))
function fine_context(c)
    p=c.problem;a=p.numerical
    n=NumericalParameters(a.N_z,a.N_b,a.P_basis,a.E_min,a.E_max,2*a.N_E-1,a.M_E,a.k_max,a.N_k,a.N_φ,a.qz_max,a.N_qz,a.η_seed)
    g=build_grids(p.physical,n,p.scales);g.ε[1:2:end]==p.grids.ε || error("nonnested canonical grids")
    for key in (:x,:wˣ,:κ,:wᵏ,:qᶻ,:wᑫᶻ,:φ);getfield(g,key)==getfield(p.grids,key) || error("nonenergy axis");end
    sp=N._scaled_physics(p.physical,p.scales);shifts=map(d->build_shift_plan(g.ε,d),(sp.Eᵖ,-sp.Eᵖ,sp.ħωᴸᴼ,-sp.ħωᴸᴼ))
    q=NEGFProblem(p.physical,n,p.scattering,p.scales,g,p.profiles,p.basis,p.kernels,shifts...,:nodal_linear,p.models)
    envelope=31*fieldbytes(q)+3*1024^3;envelope<=c.options.memory_budget_bytes || error("history1 Gamma envelope")
    refused=false
    try;build_production_cache(q;options=c.options)
    catch e;e isa ArgumentError && occursin("estimated peak memory",sprint(showerror,e)) || rethrow();refused=true;end
    refused || error("expected unmodified native full-budget refusal")
    cache,ledger=gamma_cache(q,c.options);ledger["stage_history1_envelope_bytes"]=envelope
    ledger["public_native_cache_builder_refused_before_allocation"]=refused;ledger["PI_workers"]=1
    cp=N._cavity_production_plan(q,c.h,c.options)
    ledger["cavity_plan_bytes"]=Base.summarysize(cp);ledger["nominal_LO_intervals"]=1878
    ledger["source_energy_sha256"]=GP.array_digest(p.grids.ε);ledger["fine_energy_sha256"]=GP.array_digest(g.ε);ledger["fine_weight_sha256"]=GP.array_digest(g.wᴱ)
    return GP.Context(q,cache,c.options,c.h,cp,c.kBT,c.target),ledger
end
function checkpoint_fields(out,name,G,Phi,c,cfg)
    save_fields(joinpath(out,name*".h5"),cfg,("Gamma"=>G,"raw_Phi"=>Phi);energy=c.problem.grids.ε,weights=c.problem.grids.wᴱ)
end
function evaluate_stream(c,G,group,counts,out,cfg,check)
    cap=group=="fine" ? 12 : 1;counts[group*"_attempts"]<cap || error("attempt cap")
    counts[group*"_attempts"]+=1;tag=group*"_"*lpad(counts[group*"_attempts"],3,'0');check("attempt_"*tag)
    input=GP.cone(G);counts[group*"_entered_Pi"]+=1;t=time()
    s=GP.inverse_domain(()->SS.stream_spectrum(GP,c,G;selected_indices=pick(c.problem.numerical.N_E)))
    GP.atomic_h5(joinpath(out,"response_"*tag*".h5"),cfg;expected_bytes=16*1024^2) do f
        f["q"]=s.q;f["mu_record"]="map_"*tag*".toml if fit completed; otherwise not_measured"
        f["selected_indices"]=s.selected_indices;f["energy"]=c.problem.grids.ε;f["weight_energy"]=c.problem.grids.wᴱ
        for (name,X) in s.samples;GP.write_complex(f,string(name),X);end
    end
    fit=GP.fit_mu(c.problem.grids.ε,s.q,c.target,c.kBT);green,pair=GP.fd_green(s.GR,s.A,s.condition,s.scale,c.problem.grids.ε,fit.mu,c.kBT)
    Phi,channels,sizes=GP.phi(c,green;keep_channels=false);counts[group*"_completed_Phi"]+=1
    mapped,raw,working,conversion=GP.gamma_residuals(Phi,G,c.problem.grids.wᴱ,c.problem.grids.wᵏ);F=mapped.-G;mapped=nothing
    stats=Dict{String,Any}("tag"=>tag,"seconds"=>time()-t,"residual"=>record(raw),"working_residual"=>record(working),
        "representation_conversion"=>record(conversion),"mu"=>record(fit),"input_cone"=>record(input),
        "pair_A_error"=>GP.pair_A_error(green),"FD_pair_error"=>pair,
        "maximum_condition"=>maximum(s.condition),"minimum_scale"=>minimum(s.scale),"channel_norms"=>sizes,"spectral_ledger"=>s.ledger)
    try;stats["map_cone"]=record(GP.cone(Phi));stats["A_cone"]=record(GP.cone(s.A;require_nonzero=false))
    catch e
        stats["domain_failure"]=sprint(showerror,e);open(joinpath(out,"map_"*tag*".toml"),"w") do io;TOML.print(io,stats);end
        if e isa GP.Infeasible && !isfile(joinpath(out,"first-domain-failure.h5"));checkpoint_fields(out,"first-domain-failure",G,Phi,c,cfg);end
        rethrow()
    end
    open(joinpath(out,"map_"*tag*".toml"),"w") do io;TOML.print(io,stats);end
    counts[group*"_successful"]+=1;check("complete_"*tag;enforce=false);return Phi,F,stats
end
function runtime_gate(cfg)
    cfg["runtime_authorized"]===true || error("write-only")
    cfg["history_differences"]==1 && cfg["affine_l1_max"]==4096.0 && cfg["linear_alpha"]==0.2 && cfg["aa_beta"]==1.0 && cfg["qr_condition_max"]==1e8 || error("AA identity")
    cfg["project_AA_trial"]===false && cfg["trigger_relative"]==1e-9 && cfg["trigger_passes"]==3 || error("gates/projection identity")
    cfg["memory_max_bytes"]==13019119616 && cfg["memory_high_bytes"]==12348030976 && cfg["swap_max_bytes"]==0 || error("memory identity")
    cfg["single_file_max_bytes"]==2*1024^3 && cfg["data_max_bytes"]==4*1024^3 && cfg["soft_seconds"]==840 && cfg["external_seconds"]==900 || error("resource envelope")
    cfg["workers"]==4 && cfg["BLAS_threads"]==1 && cfg["fine_loop_seconds"]==750 || error("execution identity")
    Threads.nthreads()==4 || error("Julia4");BLAS.set_num_threads(1)
    root=dirname(dirname(pathof(QCLNEGF)));readchomp(`git -C $root rev-parse HEAD`)==cfg["core_head"] || error("Core HEAD")
    isempty(readchomp(`git -C $root status --porcelain --untracked-files=no`)) || error("Core dirty")
    cg=only(filter(s->startswith(s,"0::"),readlines("/proc/self/cgroup")));path=joinpath("/sys/fs/cgroup",lstrip(split(cg,':';limit=3)[3],'/'))
    for (a,b) in (("memory.max","memory_max_bytes"),("memory.high","memory_high_bytes"),("memory.swap.max","swap_max_bytes"))
        strip(read(joinpath(path,a),String))==string(cfg[b]) || error("external cgroup "*a)
    end
    cpu=split(read(joinpath(path,"cpu.max"),String));cpu[1]!="max" && parse(Int,cpu[1])<=4*parse(Int,cpu[2]) || error("CPU cap")
    limits=split(only(filter(s->startswith(s,"Max file size"),readlines("/proc/self/limits"))))
    limits[4]!="unlimited" && parse(Int,limits[4])<=2*1024^3 || error("external file cap")
    qpath=normpath(joinpath(@__DIR__,cfg["qualifier_summary_path"]));sha_file(qpath)==cfg["qualifier_summary_sha256"] || error("qualifier SHA")
    q=TOML.parsefile(qpath);q["execution"]=="complete_operator_equivalence_only" && q["original_arithmetic_pass"] && q["Gamma_input_unchanged_after"] && q["sources_unchanged_after"] && q["definitions_unchanged_after"] || error("operator not qualified")
    for name in ("Gamma-last.h5","model.h5","input.toml","identity.toml","map-064.toml","summary.toml")
        q["source_sha256"][name]==cfg["source_sha256"][name] || error("qualifier bound to another source")
    end
    return path
end
function run_stage(input,out)
    cfg=TOML.parsefile(input);cgroup=runtime_gate(cfg);ispath(out) && error("new output required");mkpath(out)
    for name in ("run_fine.jl","protocol.md");cp(joinpath(@__DIR__,name),joinpath(out,name));end
    cp(input,joinpath(out,"input.toml"));counts=Dict(g*"_"*s=>0 for g in ("coarse_initial","fine","coarse_terminal") for s in ("attempts","entered_Pi","completed_Phi","successful"))
    begun=time();report=Dict{String,Any}("status"=>"running","SCBA_accepted"=>false,"counts"=>counts,"cgroup"=>cgroup,"history"=>1,"AA_trial_PSD_projected"=>false,"PI_FFT_controls_outside_Pi"=>2)
    function check(stage;enforce=true)
        report["stage"]=stage;report["seconds"]=time()-begun;report["maxrss_bytes"]=Sys.maxrss()
        open(joinpath(out,"summary.toml"),"w") do io;TOML.print(io,report);end
        (!enforce || time()-begun<=840) && GP.data_bytes(out)<=4*1024^3 || error("stage envelope exhausted")
        println("stage=",stage," counts=",counts," seconds=",time()-begun);flush(stdout)
    end
    try
        for (name,digest) in cfg["source_sha256"];sha_file(joinpath(SOURCE,name))==digest || error("source SHA "*name);end
        old=TOML.parsefile(joinpath(SOURCE,"input.toml"));done=TOML.parsefile(joinpath(SOURCE,"summary.toml"))
        done["execution"]=="complete" && done["pi_calls"]==64 && old["project_AA_trial"]===false || error("source scope")
        manifest=sha_file(joinpath(dirname(Base.active_project()),"Manifest.toml"));sid=TOML.parsefile(joinpath(SOURCE,"identity.toml"))
        manifest==sid["manifest_sha256"] && string(VERSION)==old["expected_julia_version"] || error("runtime source identity")
        report["input_sha256"]=sha_file(input);report["code_sha256"]=sha_file(@__FILE__);report["source_sha256"]=cfg["source_sha256"];report["manifest_sha256"]=manifest
        c=GP.build_context(old);h5open(joinpath(SOURCE,"model.h5"),"r") do f;report["model_binding"]=GP.verify_problem_model(c,f);end
        G=h5open(joinpath(SOURCE,"Gamma-last.h5"),"r") do f;size(f["Gamma_real"])==(30704,11,5,5) || error("Gamma axes");GP.read_complex(f,"Gamma");end
        qualifier=TOML.parsefile(normpath(joinpath(@__DIR__,cfg["qualifier_summary_path"])));GP.array_digest(G)==qualifier["Gamma_array_sha256"] || error("qualifier Gamma payload")
        private,cacheledger=gamma_cache(c.problem,c.options);L,equal=cache_equal(c.cache,private,c.problem,G)
        report["coarse_cache_equal"]=equal;L=private=nothing;GC.gc(true)
        Phi,F,stats=evaluate_stream(c,G,"coarse_initial",counts,out,cfg,check)
        replay=h5open(joinpath(SOURCE,"reduced-terminal.h5"),"r") do f;compare_file(Phi,f,"raw_Phi",c.problem.grids.wᴱ,c.problem.grids.wᵏ);end
        report["coarse_replay"]=replay;max(replay["relative"],replay["weighted_relative"])<=1e-10 || error("coarse Phi replay")
        fine,ledger=fine_context(c);report["fine_admission"]=ledger;Ec=c.problem.grids.ε
        report["ghost_taper_change_bound_native"]=log(2)/(2pi)*(norm(view(G,1,:,:,:))+norm(view(G,size(G,1),:,:,:)))
        G=Phi=F=nothing;GC.gc(true)
        G=h5open(joinpath(SOURCE,"Gamma-last.h5"),"r") do f;prolong_field(f,"Gamma",Ec,fine.problem.grids.ε);end
        GP.gamma_checkpoint(out,"fine-seed-Gamma",G,cfg,0)
        Phi=nothing;F=nothing;states=GP.Field[];residuals=GP.Field[];streak=0;laststats=Dict{String,Any}();events=Dict{String,Any}[]
        try
            Phi,F,laststats=evaluate_stream(fine,G,"fine",counts,out,cfg,check);push!(states,G);push!(residuals,F)
            r=laststats["residual"];streak=r["relative"]<1e-9 && r["weighted_relative"]<1e-9 ? 1 : 0
            report["events"]=events;report["fine_last"]=laststats
            checkpoint_fields(out,"fine-initial",G,Phi,fine,cfg)
            while counts["fine_attempts"]<12 && time()-begun<750 && streak<3
                trial,info=GP.aa_trial(states,residuals,fine.problem.grids.wᴱ,fine.problem.grids.wᵏ,cfg);fallback=false
                value=try;evaluate_stream(fine,trial,"fine",counts,out,cfg,check)
                catch e
                    e isa GP.Infeasible || rethrow();push!(events,Dict("attempt"=>counts["fine_attempts"],"kind"=>info.kind,"reason"=>sprint(showerror,e)));streak=0
                    if info.kind!="AAII" || counts["fine_attempts"]>=12;break;end
                    empty!(states);empty!(residuals);push!(states,G);push!(residuals,F)
                    trial = G .+ 0.2 .* F;fallback=true
                    try;evaluate_stream(fine,trial,"fine",counts,out,cfg,check)
                    catch failure
                        failure isa GP.Infeasible && push!(events,Dict("attempt"=>counts["fine_attempts"],"kind"=>"convex_fallback","reason"=>sprint(showerror,failure)))
                        rethrow()
                    end
                end
                G=trial;Phi,F,laststats=value;push!(states,G);push!(residuals,F)
                while length(states)>2;popfirst!(states);popfirst!(residuals);end
                r=laststats["residual"];streak=r["relative"]<1e-9 && r["weighted_relative"]<1e-9 ? streak+1 : 0
                push!(events,Dict("attempt"=>counts["fine_attempts"],"kind"=>fallback ? "convex_fallback" : info.kind,"AA"=>record(info),"streak"=>streak))
                report["fine_last"]=laststats;report["events"]=events;check("fine_update")
            end
        finally
            if Phi!==nothing;checkpoint_fields(out,"fine-terminal",G,Phi,fine,cfg);end
        end
        report["fine_last"]=laststats;report["events"]=events;report["strict_reduced_streak"]=streak
        empty!(states);empty!(residuals);F=value=trial=nothing;GC.gc(true)
        restricted=G[1:2:end,:,:,:];coarsePhi,coarseF,cs=evaluate_stream(c,restricted,"coarse_terminal",counts,out,cfg,check)
        checkpoint_fields(out,"coarse-restriction-terminal",restricted,coarsePhi,c,cfg)
        report["terminal_restriction"]=record(GP.norms(coarsePhi,Phi[1:2:end,:,:,:],c.problem.grids.wᴱ,c.problem.grids.wᵏ));report["coarse_terminal"]=cs
        report["status"]="completed_bounded_fine_progress_no_native_exit";check("complete")
    catch e;report["status"]="stopped";report["reason"]=sprint(showerror,e);rethrow()
    finally;report["seconds"]=time()-begun;report["maxrss_bytes"]=Sys.maxrss();open(joinpath(out,"summary.toml"),"w") do io;TOML.print(io,report);end;end
end
function tiny()
    @test pick(61407) isa Vector{Int};@test 19569 in pick(61407);@test length(pick(61407))<=128
    @test 31*(61407*11*25*16)+3*1024^3<13019119616
    cfg=Dict{String,Any}("temperature_K"=>70.0,"Nz"=>24,"Nb"=>2,"basis_periods"=>1,"NE"=>17,"energy_min_eV"=>-0.3,"energy_max_eV"=>0.9000001064962726,
        "trust_margin_meV"=>80.0,"kmax_per_nm"=>0.6,"Nk"=>3,"Nphi"=>4,"Nqz"=>5,"qzmax_per_nm"=>10.0,"seed_halfwidth_meV"=>1.25,"workers"=>4,"memory_max_bytes"=>13019119616)
    c=GP.build_context(cfg);G=fill(0.0+0im,17,3,2,2);for e in 1:17,k in 1:3;G[e,k,:,:] .= 0.03 .* Matrix{ComplexF64}(I,2,2);end
    pv,_=gamma_cache(c.problem,c.options);_,equal=cache_equal(c.cache,pv,c.problem,G);@test equal["exact_value_equality"]
    low=N.with_production_options(c.options;memory_budget_bytes=1);@test_throws ArgumentError build_production_cache(c.problem;options=low)
    s=SS.stream_spectrum(GP,c,G);old=GP.stable_spectral_state(c,G)
    @test s.GR==old[1];@test s.A==old[2];@test s.q==old[6]
    p=copy(G);p[1,1,1,2]+=1e-15im;m,r,w,h=GP.gamma_residuals(p,G,c.problem.grids.wᴱ,c.problem.grids.wᵏ)
    @test r.absolute>0 && h.absolute>0;@test m[1,1,:,:]==m[1,1,:,:]';@test p[1,1,1,2]!=m[1,1,1,2]
    aa=Dict("aa_beta"=>1.0,"linear_alpha"=>0.2,"qr_condition_max"=>1e8,"affine_l1_max"=>4096.0)
    x0=copy(G);f0 = 0.01 .* (G ./ 0.03) .- 0.5 .* x0;x1 = x0 .+ 0.2 .* f0;f1 = 0.01 .* (G ./ 0.03) .- 0.5 .* x1
    trial,info=GP.aa_trial(GP.Field[x0,x1],GP.Field[f0,f1],c.problem.grids.wᴱ,c.problem.grids.wᵏ,aa)
    @test info.kind=="AAII" && info.rank==1 && info.l1<=4096
    @test norm(trial .- (0.02 / 0.03) .* G)<=1e-12*norm(G)
end
if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS==["--tiny"];tiny()
    else;length(ARGS)==2 || error("usage: run_fine.jl --tiny | INPUT NEW_OUTPUT");run_stage(ARGS[1],abspath(ARGS[2]));end
end
