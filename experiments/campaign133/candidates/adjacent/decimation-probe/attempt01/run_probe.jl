using QCLNEGF, Unitful, LinearAlgebra, HDF5, TOML, SHA, Dates
include(joinpath(@__DIR__,"source.jl"))
include(joinpath(@__DIR__,"algorithms.jl"))
using .FrozenBoundaryAlgorithms
const N=QCLNEGF.QCLNumerics
const LABELS=["finite4","finite8","finite16","finite32","infinite_decimation"]
const PAIRS=[(1,2),(1,3),(1,4),(1,5),(2,3),(3,4),(4,5)]
const PAIR_LABELS=[LABELS[b]*"_minus_"*LABELS[a] for (a,b) in PAIRS]

function source_integrals(A,g,f)
    ne,nk,nb,_=size(A);M=zeros(ComplexF64,nk,nb,nb);F=similar(M);fill!(F,0)
    for k in 1:nk,e in 1:ne,b in 1:nb,a in 1:nb
        x=g.wᴱ[e]/(2π)*A[e,k,a,b]
        M[k,a,b]+=x;F[k,a,b]+=f[e]*x
    end
    M,F
end
charge(F,p,g)=p.g_s*sum(g.wᵏ[k]*real(sum(F[k,a,a] for a in axes(F,2))) for k in axes(F,1))
rsum(M)=maximum(norm(Matrix(view(M,k,:,:))-I)/sqrt(size(M,2)) for k in axes(M,1))

function native_c4_A(localR,problem,h,options)
    zeroK=zeros(ComplexF64,size(localR))
    family=SelfEnergyFamily(localR,zeroK,zeroK)
    plan=N._cavity_production_plan(problem,h,options)
    total,plus,minus=N._cavity_embedding_production(problem,family,h,options;plan)
    sigma=localR+total.Σᴿ
    total=plus=minus=family=zeroK=plan=nothing;GC.gc()
    GR,condition,scale=N._retarded_green_production(problem.grids.ε,h,sigma;η=0.0,options)
    sigma=nothing
    A=similar(GR)
    for k in axes(GR,2),e in axes(GR,1),b in axes(GR,4),a in axes(GR,3)
        A[e,k,a,b]=im*(GR[e,k,a,b]-conj(GR[e,k,b,a]))
    end
    (A=A,max_condition=maximum(condition),max_scale=maximum(scale))
end

function align_native!(summary,output,control,A,g,f,p,target,alignment)
    M,F=source_integrals(A,g,f)
    oldM,oldF=h5open(alignment,"r") do file
        read_complex(file["fresh_PI"],"M"),read_complex(file["fresh_PI"],"F")
    end
    mr=norm(M-oldM)/max(norm(oldM),eps(Float64))
    fr=norm(F-oldF)/max(norm(oldF),eps(Float64))
    old=TOML.parsefile(normpath(joinpath(@__DIR__,control["alignment_summary"])))
    sr=abs(rsum(M)-old["fresh_PI"]["r_sum"])
    nr=abs(charge(F,p,g)/target-old["fresh_PI"]["N_fA_over_target"])
    summary["native_c4_alignment"]=Dict("M_relative"=>mr,"F_relative"=>fr,
        "r_sum_absolute_difference"=>sr,"N_over_target_absolute_difference"=>nr,
        "r_sum"=>rsum(M),"N_fA_over_target"=>charge(F,p,g)/target)
    h5open(joinpath(output,"alignment.h5"),"w") do file
        write_complex(file,"M_native_c4",M);write_complex(file,"F_native_c4",F)
        write_complex(file,"M_attempt01_PI",oldM);write_complex(file,"F_attempt01_PI",oldF)
    end
    max(mr,fr)<=control["alignment_MF_relative_maximum"] || error("native c4 M/F alignment failed")
    sr<=control["alignment_rsum_absolute_maximum"] || error("native c4 rsum alignment failed")
    nr<=control["alignment_charge_over_target_absolute_maximum"] || error("native c4 charge alignment failed")
    M,F
end

function sweep!(Ms,Fs,completed,kstats,sensitivity,R,h,T,nativeA,g,f,control,checkpoint)
    ne,nk,nb,_=size(R);depths=Vector{Int}(control["finite_side_depths"])
    atol=Float64(control["correction_absolute_eV"])/Float64(control["E0_eV"])
    rtol=Float64(control["correction_relative"])
    maxiter=Int(control["decimation_max_iterations"])
    equation_tolerance=Float64(control["surface_backward_maximum"])
    causality_tolerance=Float64(control["surface_spectral_negative_relative_maximum"])
    Threads.@threads :static for k in 1:nk
        started=time();w=BoundaryWorkspace(nb)
        centres=[zeros(ComplexF64,nb,nb) for _ in 1:5]
        As=[zeros(ComplexF64,nb,nb) for _ in 1:5]
        localM=[zeros(ComplexF64,nb,nb) for _ in 1:5]
        localF=[zeros(ComplexF64,nb,nb) for _ in 1:5]
        conditions=zeros(5);residuals=zeros(5)
        maxcond=zeros(5);maxback=zeros(5);mineig=fill(Inf,5)
        maxcond_e=zeros(Int,5);mineig_e=zeros(Int,5)
        pairnum=zeros(7);pairden=zeros(7);occnum=zeros(7);occden=zeros(7)
        boundM=zeros(7);boundF=zeros(7);native_num=0.0;native_den=0.0
        surfmax=zeros(2);surfmin=fill(Inf,2);surfmax_e=zeros(Int,2);surfmin_e=zeros(Int,2)
        surfminabs=fill(Inf,2);surfcondmax=zeros(2)
        iterations_max=0;iterations_sum=0;iterations_hist=zeros(Int,maxiter)
        failed_stop=0;failed_equation=0;failed_causality=0;failed_any=0;first_failed_e=0
        nonfinite_nodes=0;finite_failed_nodes=0
        max_correction_ratio=0.0;max_bulk_condition_sample=0.0
        for e in 1:ne
            e%256==1 && checkpoint("sweep k=$k e=$e";announce=false)
            @inbounds for b in 1:nb,a in 1:nb
                w.D0[a,b]=(a==b ? g.ε[e] : 0.0)-h[k,a,b]-R[e,k,a,b]
            end
            finite_centres!(centres,w,T,depths;conditions,residuals)
            state=decimate!(w,T;atol_energy=atol,rtol,maxiter)
            copyto!(centres[5],w.G)
            all(G->all(isfinite,G),centres) || error("nonfinite centre block at k=$k e=$e")
            diag=diagnostics(w.D0,T,w.gR,w.gL,w.G)
            finite_diagnostics=all(isfinite,values(diag)) &&
                all(isfinite,(state.correction_energy,state.correction_scale_energy,state.bulk_condition_sample))
            nonfinite_nodes+=!finite_diagnostics
            finite_checks_pass=all(isfinite,view(conditions,1:4)) &&
                all(isfinite,view(residuals,1:4)) && maximum(view(residuals,1:4))<=equation_tolerance
            finite_failed_nodes+=!finite_checks_pass
            conditions[5]=diag.centre_condition;residuals[5]=diag.centre_residual
            iterations_max=max(iterations_max,state.iterations);iterations_sum+=state.iterations
            iterations_hist[state.iterations]+=1
            max_correction_ratio=max(max_correction_ratio,state.correction_energy/
                (atol+rtol*state.correction_scale_energy))
            max_bulk_condition_sample=max(max_bulk_condition_sample,state.bulk_condition_sample)
            eqfail=!finite_diagnostics || max(diag.surface_residual_R,diag.surface_residual_L,diag.centre_residual)>equation_tolerance
            psdfail=!finite_diagnostics || min(diag.surface_min_spectral_scaled_R,diag.surface_min_spectral_scaled_L)<
                    -causality_tolerance
            failed_stop+=!state.converged;failed_equation+=eqfail;failed_causality+=psdfail
            failed_any+=(!state.converged || eqfail || psdfail)
            if (!state.converged || eqfail || psdfail) && first_failed_e==0;first_failed_e=e;end
            for (j,r) in enumerate((diag.surface_residual_R,diag.surface_residual_L))
                if r>surfmax[j];surfmax[j]=r;surfmax_e[j]=e;end
            end
            for (j,r) in enumerate((diag.surface_min_spectral_scaled_R,diag.surface_min_spectral_scaled_L))
                if r<surfmin[j];surfmin[j]=r;surfmin_e[j]=e;end
            end
            surfminabs[1]=min(surfminabs[1],diag.surface_min_spectral_R)
            surfminabs[2]=min(surfminabs[2],diag.surface_min_spectral_L)
            surfcondmax[1]=max(surfcondmax[1],diag.surface_condition_R)
            surfcondmax[2]=max(surfcondmax[2],diag.surface_condition_L)
            for v in 1:5
                spectral!(As[v],centres[v])
                if conditions[v]>maxcond[v];maxcond[v]=conditions[v];maxcond_e[v]=e;end
                maxback[v]=max(maxback[v],residuals[v])
                ev=eigmin(Hermitian(As[v]))/max(norm(As[v]),floatmin(Float64))
                if ev<mineig[v];mineig[v]=ev;mineig_e[v]=e;end
                factor=g.wᴱ[e]/(2π)
                @inbounds for b in 1:nb,a in 1:nb
                    localM[v][a,b]+=factor*As[v][a,b]
                    localF[v][a,b]+=factor*f[e]*As[v][a,b]
                end
            end
            for (j,(va,vb)) in enumerate(PAIRS)
                d2=0.0;r2=0.0
                @inbounds for b in 1:nb,a in 1:nb
                    d2+=abs2(As[vb][a,b]-As[va][a,b]);r2+=abs2(As[va][a,b])
                end
                weight=g.wᴱ[e]*g.wᵏ[k]
                pairnum[j]+=weight*d2;pairden[j]+=weight*r2
                occnum[j]+=weight*f[e]^2*d2;occden[j]+=weight*f[e]^2*r2
                boundM[j]+=g.wᴱ[e]/(2π)*sqrt(d2)
                boundF[j]+=g.wᴱ[e]/(2π)*f[e]*sqrt(d2)
            end
            @inbounds for b in 1:nb,a in 1:nb
                weight=g.wᴱ[e]*g.wᵏ[k]
                native_num+=weight*abs2(As[1][a,b]-nativeA[e,k,a,b])
                native_den+=weight*abs2(nativeA[e,k,a,b])
            end
        end
        for v in 1:5
            Ms[v][k,:,:].=localM[v];Fs[v][k,:,:].=localF[v]
        end
        sensitivity[k]=(;pairnum,pairden,occnum,occden,boundM,boundF,native_num,native_den)
        kstats[k]=Dict("k_index"=>k,"elapsed_seconds"=>time()-started,
            "max_condition"=>maxcond,"max_condition_energy_indices"=>maxcond_e,
            "max_backward"=>maxback,"minimum_A_eigenvalue_scaled"=>mineig,
            "minimum_A_energy_indices"=>mineig_e,"surface_max_backward"=>surfmax,
            "surface_max_backward_energy_indices"=>surfmax_e,"surface_min_spectral_scaled"=>surfmin,
            "surface_min_spectral_absolute"=>surfminabs,"surface_max_condition"=>surfcondmax,
            "surface_min_spectral_energy_indices"=>surfmin_e,"max_iterations"=>iterations_max,
            "mean_iterations"=>iterations_sum/ne,"iteration_histogram"=>iterations_hist,
            "stop_failed_nodes"=>failed_stop,"equation_failed_nodes"=>failed_equation,
            "causality_failed_nodes"=>failed_causality,"first_failed_energy_index"=>first_failed_e,
            "any_failed_nodes"=>failed_any,
            "nonfinite_diagnostic_nodes"=>nonfinite_nodes,"finite_centre_failed_nodes"=>finite_failed_nodes,
            "maximum_correction_over_threshold"=>max_correction_ratio,
            "max_bulk_condition_sample"=>max_bulk_condition_sample)
        completed[k]=true
        println("completed_k=",k," elapsed_s=",time()-started);flush(stdout)
    end
end

function save_sweep(output,Ms,Fs,completed,g,kstats,sensitivity)
    h5open(joinpath(output,"boundary_integrals.h5"),"w") do file
        file["schema"]="campaign133-frozen-PI-boundary-limit-v1"
        file["variant_labels"]=LABELS;file["sensitivity_pair_labels"]=PAIR_LABELS
        file["completed_k"]=Int.(completed)
        file["energy"]=g.ε;file["weight_energy"]=g.wᴱ;file["k"]=g.κ;file["weight_k"]=g.wᵏ
        for v in 1:5
            group=create_group(file,LABELS[v]);write_complex(group,"M",Ms[v]);write_complex(group,"F",Fs[v])
        end
        for k in eachindex(completed)
            completed[k] || continue
            group=create_group(file,"k"*string(k));rec=sensitivity[k]
            for key in propertynames(rec);group[string(key)]=getproperty(rec,key);end
        end
    end
end

function selected_oracles!(output,summary,R,h,T,nativeA,g,f,control,checkpoint)
    records=Dict{String,Any}();nb=size(T,1)
    summary["selected_oracles"]=records # retain partial witnesses on failure
    h5open(joinpath(output,"selected_oracles.h5"),"w") do file
        for k in control["bloch_k_indices"]
            checkpoint("selected_oracles k=$k")
            e=argmax([f[j]*real(sum(nativeA[j,k,a,a] for a in 1:nb)) for j in eachindex(f)])
            D0=g.ε[e]*I-Matrix(view(h,k,:,:))-Matrix(view(R,e,k,:,:))
            w=BoundaryWorkspace(nb);copyto!(w.D0,D0)
            state=decimate!(w,T;atol_energy=control["correction_absolute_eV"]/control["E0_eV"],
                rtol=control["correction_relative"],maxiter=control["decimation_max_iterations"])
            group=create_group(file,"occupied_peak_k"*string(k))
            group["energy_index"]=e;group["energy_eV"]=g.ε[e]*control["E0_eV"]
            write_complex(group,"D0",D0);write_complex(group,"Tplus",T)
            write_complex(group,"surface_R",w.gR);write_complex(group,"surface_L",w.gL)
            write_complex(group,"bulk_decimation",w.G)
            finite=[zeros(ComplexF64,nb,nb)]
            finite_centres!(finite,w,T,[4])
            oracle=full_finite_centre(D0,T,4)
            finite_error=norm(finite[1]-oracle.centre)/norm(oracle.centre)
            write_complex(group,"finite4_full45_oracle",oracle.centre)
            previous=nothing;levels=Dict{String,Any}();converged=false;agreement=false;all_levels_finite=true
            for level in control["bloch_theta_levels"]
                checkpoint("Bloch k=$k theta=$level")
                B=bloch_integral(D0,T,level)
                refinement=isnothing(previous) ? NaN : norm(B.G-previous)/max(norm(B.G),eps(Float64))
                difference=norm(B.G-w.G)/max(norm(B.G),eps(Float64))
                level_finite=all(isfinite,B.G) && isfinite(B.max_condition) && isfinite(difference) &&
                    (isnothing(previous) || isfinite(refinement))
                all_levels_finite &= level_finite
                levels[string(level)]=Dict("refinement_relative"=>refinement,
                    "decimation_relative"=>difference,"max_condition"=>B.max_condition,
                    "all_finite"=>level_finite,"refinement_measured"=>!isnothing(previous))
                write_complex(group,"Bloch_"*string(level),B.G)
                converged=level_finite && !isnothing(previous) && refinement<=control["bloch_refinement_relative_target"]
                agreement=level_finite && difference<=control["bloch_decimation_relative_target"]
                previous=B.G
                converged && agreement && break
            end
            all_finite=all_levels_finite && all(isfinite,w.G) &&
                all(isfinite,(finite_error,oracle.backward,oracle.full_condition))
            finite_check_pass=all_finite && finite_error<=control["independent_c4_native_A_relative_maximum"] &&
                oracle.backward<=control["surface_backward_maximum"]
            records["k"*string(k)]=Dict("energy_index"=>e,"energy_eV"=>g.ε[e]*control["E0_eV"],
                "selection"=>"argmax_E f(E;fixed source mu,T) Tr A_native_c4(E,k)",
                "finite4_vs_full45_relative"=>finite_error,"full45_backward"=>oracle.backward,
                "full45_condition"=>oracle.full_condition,"decimation_stop_converged"=>state.converged,
                "all_finite"=>all_finite,"finite_check_pass"=>finite_check_pass,
                "Bloch_refinement_converged"=>converged,"Bloch_decimation_agreement"=>agreement,"theta_levels"=>levels)
            finite_check_pass || error("actual finite4 full45/nonfinite witness failed")
        end
    end
    summary["selected_oracles"]=records
end

function main()
    length(ARGS)==1 || error("usage: run_probe.jl NEW_OUTPUT_DIRECTORY")
    control=TOML.parsefile(joinpath(@__DIR__,"input.toml"));BLAS.set_num_threads(1)
    Threads.nthreads()==control["Julia_threads"] || error("Julia thread count differs")
    output=abspath(ARGS[1]);ispath(output) && error("output already exists");mkpath(output)
    for name in ("run_probe.jl","algorithms.jl","source.jl","input.toml","protocol.md","algorithm.md")
        cp(joinpath(@__DIR__,name),joinpath(output,name))
    end
    started=time();stage=Ref("identity");phase_started=Ref(started);timings=Dict{String,Float64}()
    Ms=nothing;Fs=nothing;completed=nothing;g=nothing
    kstats=Dict{String,Any}[];sensitivity=Any[]
    summary=Dict{String,Any}("id"=>control["id"],"status"=>"running",
        "started_utc"=>string(Dates.now(Dates.UTC)),"scientific_scope"=>control["scope"],
        "SCBA_Poisson_discretization_experiment_validation"=>"not_measured",
        "infinite_numerical_status"=>"not_evaluated","surface_numerical_status"=>"not_evaluated",
        "Julia_threads"=>Threads.nthreads(),"BLAS_threads"=>BLAS.get_num_threads())
    function checkpoint(name;announce=true)
        if announce # worker checks do not mutate shared stage/timings
            timings[stage[]]=get(timings,stage[],0.0)+time()-phase_started[]
            stage[]=name;phase_started[]=time()
        end
        time()-started<=control["script_soft_stop_seconds"] || error("soft elapsed-time stop")
        bytes=sum(filesize(joinpath(root,file)) for (root,_,files) in walkdir(output) for file in files;init=0)
        bytes<=control["maximum_output_bytes"] || error("output cap exceeded")
        if announce;println("stage=",name," elapsed_s=",time()-started);flush(stdout);end
    end
    try
        for name in ("source","alignment_h5","alignment_summary")
            path=normpath(joinpath(@__DIR__,control[name]))
            sha_file(path)==control[name*"_sha256"] || error("identity mismatch: $name")
            summary[name*"_sha256"]=control[name*"_sha256"]
        end
        root=dirname(dirname(pathof(QCLNEGF)))
        head=readchomp(`git -C $root rev-parse HEAD`)
        head==control["core_head"] || error("Core HEAD differs: $head")
        isempty(readchomp(`git -C $root status --porcelain`)) || error("accepted Core is dirty")
        summary["core_head"]=head;summary["core_import_path"]=pathof(QCLNEGF)
        project=Base.active_project();manifest=joinpath(dirname(project),"Manifest.toml")
        summary["project_sha256"]=sha_file(project);summary["manifest_sha256"]=sha_file(manifest)
        summary["julia_version"]=string(VERSION)
        summary["code_sha256"]=Dict(name=>sha_file(joinpath(@__DIR__,name)) for name in
            ("run_probe.jl","algorithms.jl","source.jl","input.toml","protocol.md","algorithm.md"))
        checkpoint("small_independent_oracles")
        summary["small_selfchecks"]=selfcheck()
        source=normpath(joinpath(@__DIR__,control["source"]));srcinput=joinpath(dirname(source),"input.toml")
        prior=TOML.parsefile(normpath(joinpath(@__DIR__,control["alignment_summary"])))
        for (relative,key) in (("identity.toml","identity_sha256"),
                               ("markers/final_markers.toml","marker_sha256"),
                               ("run_inner.jl","source_writer_sha256"))
            digest=sha_file(joinpath(dirname(source),relative))
            digest==prior[key] || error("source provenance differs: $relative")
            summary[key]=digest
        end
        sha_file(srcinput)==control["source_input_sha256"] || error("source input SHA differs")
        cfg=TOML.parsefile(srcinput);ident=TOML.parsefile(joinpath(dirname(source),"identity.toml"))
        markers=TOML.parsefile(joinpath(dirname(source),"markers/final_markers.toml"))
        ident["input_sha256"]==control["source_input_sha256"] || error("writer input identity differs")
        sha_file(joinpath(dirname(source),"run_inner.jl"))==ident["script_sha256"] || error("source writer differs")
        markers["equilibrium_mu_eV"]==control["mu_eV"] || error("fixed source mu differs")
        cfg["voltage_mV"]==0 && cfg["temperature_K"]==70 && cfg["embedding_periods"]==4 || error("source case differs")
        cfg["NE"]==30721 && cfg["Nk"]==11 && cfg["Nb"]==5 && cfg["alloy"] &&
            cfg["embedding"]=="finite_chain" || error("source discretization/families differ")
        control["finite_side_depths"]==[4,8,16,32] || error("variant labels require fixed depths")
        checkpoint("source_frame")
        problem,h,basis_comparison=h5open(source,"r") do file;source_frame(file,cfg);end
        g=problem.grids;p=problem.physical;s=problem.scales;n=problem.numerical
        control["E0_eV"]=s.E₀_eV
        iszero(QCLNEGF._scaled_physics(p,s).Eᵖ) || error("not identical zero-field cells")
        sha_array(h)==control["projected_h_sha256"] || error("projected source h differs")
        norm(problem.basis.T₋-problem.basis.T₊')==0 || error("Tminus differs from Tplus adjoint")
        summary["basis_source_and_rebuild_comparison"]=basis_comparison
        summary["projected_h_sha256"]=sha_array(h);summary["Tplus_sha256"]=sha_array(problem.basis.T₊)
        summary["mu_eV"]=control["mu_eV"];summary["E0_eV"]=s.E₀_eV
        target=QCLNEGF._scaled_physics(p,s).Nᴰ²ᴰ
        kBT=QCLNEGF._electronvolts(CODATA.kᴮₑᵥ*p.Tᴸ)/s.E₀_eV
        f=[QCLNEGF._fermi(e,control["mu_eV"]/s.E₀_eV,kBT) for e in g.ε]
        sourceR=zeros(ComplexF64,n.N_E,n.N_k,n.N_b,n.N_b)
        h5open(source,"r") do file
            for k in 1:n.N_k,name in LOCAL_NAMES
                sourceR[:,k,:,:].+=read_k(file["selfenergy/"*name],"retarded",k)
            end
        end
        gamma=similar(sourceR)
        for k in 1:n.N_k,e in 1:n.N_E,b in 1:n.N_b,a in 1:n.N_b
            gamma[e,k,a,b]=im*(sourceR[e,k,a,b]-conj(sourceR[e,k,b,a]))
        end
        sha_array(gamma)==control["frozen_gamma_sha256"] || error("frozen local gamma differs")
        summary["frozen_gamma_sha256"]=sha_array(gamma)
        summary["local_families_summed_once"]=LOCAL_NAMES
        summary["excluded_families"]=["embedding","embedding_plus","embedding_minus"]
        sourceR=nothing;GC.gc();checkpoint("native_uniform_PI")
        plan=N._ProductionPIFFTHilbertPlan(g.ε,n.N_k*n.N_b^2;column_chunk=1,worker_count=4)
        R,roundoff=N._production_pi_fft_hilbert_transform(gamma,g.ε,g.wᴱ,plan)
        summary["PI_roundoff"]=roundoff;summary["PI_geometry"]=Dict(string(k)=>v for (k,v) in pairs(plan.geometry))
        R .-=0.5im.*gamma;gamma=plan=nothing;GC.gc()
        summary["local_PI_R_sha256"]=sha_array(R)
        options=ProductionOptions(worker_count=4,hilbert_columns=1,memory_budget_bytes=cfg["memory_budget_bytes"],
            algorithms=AlgorithmOptions(embedding=:finite_chain,embedding_periods=4,hilbert=:product_integration))
        checkpoint("native_c4_alignment")
        native=native_c4_A(R,problem,h,options)
        align_native!(summary,output,control,native.A,g,f,p,target,
            normpath(joinpath(@__DIR__,control["alignment_h5"])))
        summary["native_c4_max_condition"]=native.max_condition
        summary["native_c4_max_scale"]=native.max_scale
        Ms=[fill(complex(NaN,NaN),n.N_k,n.N_b,n.N_b) for _ in 1:5]
        Fs=deepcopy(Ms);completed=fill(false,n.N_k)
        kstats=[Dict{String,Any}() for _ in 1:n.N_k];sensitivity=Any[nothing for _ in 1:n.N_k]
        checkpoint("independent_finite_and_decimation_sweep")
        sweep!(Ms,Fs,completed,kstats,sensitivity,R,h,problem.basis.T₊,native.A,g,f,control,checkpoint)
        save_sweep(output,Ms,Fs,completed,g,kstats,sensitivity)
        summary["per_k_diagnostics"]=Dict("k"*string(k)=>kstats[k] for k in 1:n.N_k)
        summary["target_dimensionless"]=target;summary["variant_labels"]=LABELS
        summary["measure"]="source dimensionless wE/(2pi), wk, spin degeneracy; source fixed mu,T; no density normalization"
        summary["variants"]=Dict(LABELS[v]=>Dict("r_sum"=>rsum(Ms[v]),
            "N_fA_over_target"=>charge(Fs[v],p,g)/target) for v in 1:5)
        native_relative=sqrt(sum(x.native_num for x in sensitivity)/max(sum(x.native_den for x in sensitivity),eps(Float64)))
        summary["independent_finite4_minus_native_A_weighted_relative"]=native_relative
        native_relative<=control["independent_c4_native_A_relative_maximum"] || error("independent finite4/native spectrum mismatch")
        summary["sensitivity"]=Dict(PAIR_LABELS[j]=>Dict(
            "weighted_L2_absolute"=>sqrt(sum(x.pairnum[j] for x in sensitivity)),
            "weighted_L2_reference"=>sqrt(sum(x.pairden[j] for x in sensitivity)),
            "occupied_L2_absolute"=>sqrt(sum(x.occnum[j] for x in sensitivity)),
            "occupied_L2_reference"=>sqrt(sum(x.occden[j] for x in sensitivity)),
            "delta_M_Frobenius"=>norm(Ms[b]-Ms[a]),"delta_F_Frobenius"=>norm(Fs[b]-Fs[a]),
            "delta_N_over_target"=>(charge(Fs[b],p,g)-charge(Fs[a],p,g))/target,
            "discrete_triangle_N_bound_over_target"=>p.g_s*sqrt(n.N_b)*sum(g.wᵏ[k]*sensitivity[k].boundF[j] for k in 1:n.N_k)/target
        ) for (j,(a,b)) in enumerate(PAIRS))
        summary["surface_failed_nodes"]=sum(x["any_failed_nodes"] for x in kstats)
        summary["surface_numerical_status"]=summary["surface_failed_nodes"]==0 ? "surface_witnesses_pass" : "unresolved_surface_nodes"
        summary["finite_centre_failed_nodes"]=sum(x["finite_centre_failed_nodes"] for x in kstats)
        summary["nonfinite_diagnostic_nodes"]=sum(x["nonfinite_diagnostic_nodes"] for x in kstats)
        integrals_finite=all(M->all(isfinite,M) && isfinite(norm(M)),Ms) &&
            all(F->all(isfinite,F) && isfinite(norm(F)),Fs) &&
            all(x->all(isfinite,x.pairnum) && all(isfinite,x.pairden) &&
                   all(isfinite,x.occnum) && all(isfinite,x.occden) &&
                   all(isfinite,x.boundM) && all(isfinite,x.boundF),sensitivity)
        summary["all_boundary_integrals_and_sensitivities_finite"]=integrals_finite
        checkpoint("selected_Bloch_and_full45_oracles")
        selected_oracles!(output,summary,R,h,problem.basis.T₊,native.A,g,f,control,checkpoint)
        selected=summary["selected_oracles"]
        selected_pass=length(selected)==length(control["bloch_k_indices"]) && all(
            rec["all_finite"] && rec["finite_check_pass"] && rec["decimation_stop_converged"] &&
            rec["Bloch_refinement_converged"] && rec["Bloch_decimation_agreement"] for rec in values(selected))
        combined_pass=all(completed) && summary["surface_failed_nodes"]==0 &&
            summary["finite_centre_failed_nodes"]==0 && summary["nonfinite_diagnostic_nodes"]==0 &&
            integrals_finite && selected_pass
        summary["infinite_numerical_status"]=combined_pass ?
            "surface_finite_and_selected_Bloch_witnesses_pass" : "unresolved_boundary_witnesses"
        checkpoint("complete");summary["status"]="completed_execution_only"
    catch err
        summary["status"]="failed";summary["error"]=sprint(showerror,err);rethrow()
    finally
        if !isnothing(Ms)
            save_sweep(output,Ms,Fs,completed,g,kstats,sensitivity)
            summary["completed_k_indices"]=findall(completed)
            summary["per_k_diagnostics"]=Dict("k"*string(k)=>kstats[k] for k in eachindex(completed) if completed[k])
        end
        summary["last_stage"]=stage[];summary["elapsed_seconds"]=time()-started
        timings[stage[]]=get(timings,stage[],0.0)+time()-phase_started[]
        summary["stage_timings_seconds"]=timings
        summary["finished_utc"]=string(Dates.now(Dates.UTC))
        open(joinpath(output,"summary.toml"),"w") do io;TOML.print(io,summary);end
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
