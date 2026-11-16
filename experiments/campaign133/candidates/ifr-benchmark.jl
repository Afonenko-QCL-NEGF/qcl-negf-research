using QCLNEGF, Unitful, HDF5, LinearAlgebra, TOML, SHA, Dates
const N = QCLNEGF.QCLNumerics
sha_file(p)=bytes2hex(open(sha256,p))
sha_array(A)=bytes2hex(sha256(reinterpret(UInt8,vec(A))))
complex_read(g,name)=complex.(read(g[name*"_real"]),read(g[name*"_imag"]))
function numerical(c,ne)
    NumericalParameters(N_z=c["Nz"],N_b=c["Nb"],P_basis=c["basis_periods"],
        E_min=c["energy_min_eV"]*u"eV",E_max=c["energy_max_eV"]*u"eV",N_E=ne,
        M_E=c["trust_margin_meV"]*u"meV",k_max=c["kmax_per_nm"]*u"nm^-1",N_k=c["Nk"],
        N_φ=c["Nphi"],qz_max=c["qzmax_per_nm"]*u"nm^-1",N_qz=c["Nqz"],η_seed=c["seed_meV"]*u"meV")
end
function error_metrics(A,B)
    difference,reference,maxabs,maxreference,worst=0.0,0.0,0.0,0.0,1
    for i in eachindex(A,B)
        isfinite(A[i]) && isfinite(B[i]) || error("nonfinite contraction output")
        delta=abs(A[i]-B[i]); value=abs(B[i])
        difference+=delta^2;reference+=value^2;maxreference=max(maxreference,value)
        if delta>maxabs;maxabs=delta;worst=i;end
    end
    Dict("absolute_Frobenius"=>sqrt(difference),"reference_Frobenius"=>sqrt(reference),
        "relative_Frobenius"=>sqrt(difference)/max(sqrt(reference),eps(Float64)),
        "maximum_absolute_entry_error"=>maxabs,"maximum_reference_entry"=>maxreference,
        "maximum_entry_scaled_error"=>maxabs/max(maxreference,eps(Float64)),
        "worst_linear_index"=>worst)
end
function main()
    length(ARGS)==2 || error("usage: ifr-benchmark.jl INPUT.toml NEW_OUTPUT")
    input,output=abspath.(ARGS);ispath(output) && error("output exists")
    cfg=TOML.parsefile(input);base=dirname(input)
    BLAS.set_num_threads(1)
    Threads.nthreads(:default)==cfg["workers"] || error("Julia worker launch differs")
    cfg["repetitions"]>=3 && isodd(cfg["repetitions"]) || error("need odd repetitions >=3")
    mkpath(output);cp(input,joinpath(output,"input.toml"));cp(@__FILE__,joinpath(output,"ifr-benchmark.jl"))
    started=time();summary=Dict{String,Any}("status"=>"running","start_utc"=>string(Dates.now(Dates.UTC)),
        "scope"=>"IFR operator-only on preserved real full GR; no SCBA/LO/physics change; no full-field outputs",
        "script_sha256"=>sha_file(@__FILE__),"input_sha256"=>sha_file(input),
        "Julia_threads"=>Threads.nthreads(:default),"BLAS_threads"=>BLAS.get_num_threads())
    function phase(name)
        time()-started<285 || error("phase boundary time stop")
        println("phase=",name," elapsed_s=",time()-started);flush(stdout)
        summary["phase"]=name
    end
    try
        root=dirname(dirname(pathof(QCLNEGF)))
        realpath(root)==realpath(joinpath(base,cfg["core_worktree"])) || error("actual imported Core differs")
        readchomp(`git -C $root rev-parse HEAD`)==cfg["core_base"] || error("Core base differs")
        patch=bytes2hex(sha256(read(`git -C $root diff`)))
        implementation=bytes2hex(sha256(read(`git -C $root diff -- src`)))
        patch==cfg["core_patch_sha256"] || error("Core candidate patch differs")
        implementation==cfg["core_implementation_diff_sha256"] || error("Core implementation differs")
        summary["core_source"]=pathof(QCLNEGF);summary["core_base"]=cfg["core_base"]
        summary["core_patch_sha256"]=patch;summary["core_implementation_diff_sha256"]=implementation
        summary["julia_version"]=string(VERSION)
        summary["project_sha256"]=sha_file(Base.active_project())
        summary["manifest_sha256"]=sha_file(joinpath(dirname(Base.active_project()),"Manifest.toml"))
        source=normpath(joinpath(base,cfg["source_state"]));sourceinput=normpath(joinpath(base,cfg["source_input"]))
        phase("stream_source_hash")
        hashstart=time();digest=sha_file(source);summary["source_file_hash_seconds"]=time()-hashstart
        digest==cfg["source_state_sha256"] || error("preserved source state differs")
        sha_file(sourceinput)==cfg["source_input_sha256"] || error("source input differs")
        c=TOML.parsefile(sourceinput)
        c["NE"]==cfg["full_energy_nodes"] && c["Nb"]==5 && c["Nk"]==11 || error("source dimensions differ")
        summary["source_state_sha256"]=digest;summary["source_input_sha256"]=cfg["source_input_sha256"]
        summary["source_dataset"]=cfg["source_green_dataset"]
        p=reference_parameters(V_period=c["voltage_mV"]*u"mV",Tᴸ=c["temperature_K"]*u"K",
            Tᴸᴼ=c["temperature_K"]*u"K",ΔV_alloy=.60u"eV",Ω₀=4.5168e-29u"m^3")
        s=ScaleSystem();smalln=numerical(c,cfg["model_energy_nodes"]);fulln=numerical(c,cfg["full_energy_nodes"])
        fullgrid=build_grids(p,fulln,s)
        options=ProductionOptions(worker_count=cfg["workers"],energy_chunk=cfg["chunk"],
            memory_budget_bytes=c["memory_budget_bytes"],
            algorithms=AlgorithmOptions(hilbert=:direct,contraction=:dense_blas,kernel_build=:direct))
        phase("canonical_energy_independent_IFR_model")
        # Only the energy-independent IFR operator is constructed. NE17 is
        # explicitly not a replacement full solver grid or memory admission.
        problem=build_problem_production(;physical=p,numerical=smalln,scales=s,
            production_options=options,
            scattering=ScatteringOptions(LO=false,acoustic=false,impurity=false,IFR=true,alloy=false)).problem
        shape=(fulln.N_E,fulln.N_k,fulln.N_b,fulln.N_b)
        G=Array{ComplexF64}(undef,shape)
        K,q,wk=h5open(source,"r") do file
            read(file["schema"])=="campaign133-inner-raw-v3" || error("source schema differs")
            read(file["energy"])==fullgrid.ε && read(file["weight_energy"])==fullgrid.wᴱ || error("source energy axis differs")
            read(file["k"])==problem.grids.κ && read(file["weight_k"])==problem.grids.wᵏ || error("source k measure differs")
            read(file["x"])==problem.grids.x && read(file["weight_x"])==problem.grids.wˣ || error("source spatial measure differs")
            read(file["energy_scale_eV"])==s.E₀_eV && read(file["length_scale_m"])==s.L₀_m || error("source scales differ")
            read(file["spin_degeneracy"])==p.g_s || error("source spin differs")
            for (field,name) in ((:Φ,"Phi"),(:χ,"chi"),(:H₀,"H0"),(:Z,"Z"),(:M⁻¹,"Minv"),(:T₊,"Tplus"),(:T₋,"Tminus"))
                saved=complex_read(file["basis"],name)
                saved==getfield(problem.basis,field) || error("source basis frame differs: $name")
            end
            savedK=complex_read(file["kernels/IFR"],"normalized_K");savedq=read(file["kernels/IFR/qK"])
            savedK==problem.kernels.K[:IFR] && savedq==problem.kernels.qᴷ[:IFR] || error("canonical IFR differs from source kernel")
            name=cfg["source_green_dataset"]
            size(file["green/"*name*"_real"])==shape || error("source Green shape differs")
            for k in 1:fulln.N_k
                re=file["green/"*name*"_real"][:,k,:,:]
                iv=file["green/"*name*"_imag"][:,k,:,:]
                for b in 1:fulln.N_b,a in 1:fulln.N_b,e in 1:fulln.N_E
                    G[e,k,a,b]=complex(re[e,a,b],iv[e,a,b])
                end
            end
            (savedK,savedq,read(file["weight_k"]))
        end
        all(isfinite,G) && all(>(0),wk) || error("invalid source field/weights")
        cache=build_production_cache(problem;options)
        structured=cache.kernels[:IFR]
        structured isa N._SeparableIFRProductionKernel || error("canonical structured IFR was not selected")
        dense=QCLNEGF._production_kernel(K,AlgorithmOptions(contraction=:dense_blas))
        dense isa DenseProductionKernel || error("reference operator is not the original dense path")
        summary["source_field_sha256"]=sha_array(G);summary["source_field_bytes"]=sizeof(G)
        summary["source_kernel_sha256"]=sha_array(K);summary["source_weight_k_sha256"]=sha_array(wk)
        summary["qK"]=q;summary["model_NE"]=smalln.N_E;summary["field_shape"]=collect(shape)
        summary["model_note"]="IFR independent of energy; canonical NE17 builder checked bitwise against source basis/K/qK/wk; full preserved GR has NE30721"
        summary["structured_type"]=string(typeof(structured));summary["dense_type"]=string(typeof(dense))
        summary["structured_operator_bytes"]=N._kernel_operator_bytes(structured)
        summary["dense_operator_bytes"]=N._kernel_operator_bytes(dense)
        summary["workers"]=cfg["workers"];summary["chunk"]=cfg["chunk"]
        output_dense=similar(G);output_structured=similar(G)
        apply!(out,op)=N._production_static_contraction!(out,op,G,wk,q,cfg["chunk"];
            α=1.0,β=0.0,parallel_backend=:threads,worker_count=cfg["workers"])
        phase("warm_each_once")
        warmdense=@timed apply!(output_dense,dense)
        warmstructured=@timed apply!(output_structured,structured)
        summary["warm_dense_seconds"]=warmdense.time;summary["warm_structured_seconds"]=warmstructured.time
        records=Dict{String,Any}[]
        for repeat in 1:cfg["repetitions"]
            phase("measured_pair_"*string(repeat));GC.gc()
            order=isodd(repeat) ? (:dense,:structured) : (:structured,:dense)
            for label in order
                out,op=label===:dense ? (output_dense,dense) : (output_structured,structured)
                measured=@timed apply!(out,op)
                push!(records,Dict("repetition"=>repeat,"operator"=>string(label),
                    "seconds"=>measured.time,"allocated_bytes"=>measured.bytes,"GC_seconds"=>measured.gctime))
            end
        end
        metrics=error_metrics(output_structured,output_dense)
        metrics["relative_Frobenius"]<=cfg["relative_arithmetic_limit"] &&
        metrics["maximum_entry_scaled_error"]<=cfg["relative_arithmetic_limit"] || error("operator arithmetic mismatch: $metrics")
        summary["all_complex_entry_comparison"]=metrics
        summary["source_field_sha256_after"]=sha_array(G)
        summary["source_field_sha256_after"]==summary["source_field_sha256"] || error("source field changed")
        summary["dense_output_sha256"]=sha_array(output_dense)
        summary["structured_output_sha256"]=sha_array(output_structured)
        summary["samples"]=records
        for label in ("dense","structured")
            values=sort([r["seconds"] for r in records if r["operator"]==label])
            allocations=[r["allocated_bytes"] for r in records if r["operator"]==label]
            summary[label]=Dict("median_seconds"=>values[cld(length(values),2)],"minimum_seconds"=>first(values),
                "allocated_bytes_minimum"=>minimum(allocations),"allocated_bytes_maximum"=>maximum(allocations))
        end
        summary["operator_median_speedup"]=summary["dense"]["median_seconds"]/summary["structured"]["median_seconds"]
        summary["timing_scope"]="in-place apply including bounded scratch allocation/packing; excludes source hashing/read/model/cache construction/warmup; alternating pair order, GC before each pair"
        summary["end_to_end_SCBA_speedup"]="not_measured"
        summary["status"]="completed";phase("complete")
    catch error
        summary["status"]="failed";summary["error"]=sprint(showerror,error);rethrow()
    finally
        summary["elapsed_seconds"]=time()-started;summary["finish_utc"]=string(Dates.now(Dates.UTC))
        open(joinpath(output,"summary.toml"),"w") do io;TOML.print(io,summary);end
        bytes=sum(filesize(joinpath(root,file)) for (root,_,files) in walkdir(output) for file in files;init=0)
        bytes<=cfg["output_limit_bytes"] || error("small-output cap exceeded")
    end
    println(summary)
end
main()
