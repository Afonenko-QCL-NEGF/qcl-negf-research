using QCLNEGF, Unitful, LinearAlgebra, HDF5, TOML, SHA, Dates
const BN=QCLNEGF
const QN=QCLNEGF.QCLNumerics
const ComplexField=Array{ComplexF64,4}
digest(path)=bytes2hex(open(sha256,path))
array_digest(A)=bytes2hex(sha256(reinterpret(UInt8,vec(A))))
readcomplex(file,path)=complex.(read(file[path*"_real"]),read(file[path*"_imag"]))
ratio(a,b)=iszero(b) ? (iszero(a) ? 0.0 : Inf) : a/b
write_toml(path,value)=open(io->TOML.print(io,value),path,"w")

"""Raw component norms, with source quadrature as a separate measure.
No dimensionful factor is silently introduced into the dimensionless operator.
Near zero: absolute norm and scale are authoritative; relative=Inf if scale0.
"""
function field_norms(problem::P,X::ComplexField) where {P}
    g=problem.grids;ne,nk,nb,_=size(X)
    full=trusted=unweighted=0.;maximum_block=0.;where=(1,1)
    for m=1:nk,e=1:ne
        sq=0.
        @inbounds for b=1:nb,a=1:nb;sq+=abs2(X[e,m,a,b]);end
        unweighted+=sq
        w=g.wᴱ[e]*g.wᵏ[m]*problem.physical.g_s/(2π)
        full+=w*sq
        physicalE=g.ε[e]*problem.scales.E₀_eV
        -.22<=physicalE<=.82 && (trusted+=w*sq)
        if sq>maximum_block;maximum_block=sq;where=(e,m);end
    end
    Dict("unweighted_Frobenius"=>sqrt(unweighted),"full_weighted_L2"=>sqrt(full),
        "trusted_weighted_L2"=>sqrt(trusted),"maximum_block_Frobenius"=>sqrt(maximum_block),
        "maximum_block_energy_index"=>where[1],"maximum_block_k_index"=>where[2])
end

function hermitian_parts!(H::ComplexField,J::ComplexField,X::ComplexField)
    ne,nk,nb,_=size(X)
    @inbounds for m=1:nk,e=1:ne,b=1:nb,a=1:nb
        H[e,m,a,b]=(X[e,m,a,b]+conj(X[e,m,b,a]))/2
        J[e,m,a,b]=(X[e,m,a,b]-conj(X[e,m,b,a]))/2
    end
    nothing
end


# Keep every element loop behind concrete-array argument barriers. The native
# core fields are already concrete; these boundaries prevent inference loss from
# outer heterogeneous dictionaries and captured locals from reaching the loops.
function primary_H(X::ComplexField,problem::P,cache::C,options::O)::ComplexField where {P,C,O}
    QN.production_fft_hilbert_transform(X,problem.grids.ε,problem.grids.wᴱ,
        cache.fft_hilbert_plan;return_roundoff=false,worker_count=4)
end

function input_identity_defect(green::GreenState)::ComplexField
    defect=similar(green.A)
    @inbounds for i in eachindex(defect)
        defect[i]=green.A[i]-(-im*green.Gˡ[i]+im*green.Gᵍ[i])
    end
    return defect
end

function input_hilbert_defect(green::GreenState,problem::P,cache::C,options::O)::ComplexField where {P,C,O}
    X=primary_H(green.A,problem,cache,options)
    ne,nk,nb,_=size(X)
    @inbounds for m=1:nk,e=1:ne,b=1:nb,a=1:nb
        X[e,m,a,b]=(green.Gᴿ[e,m,a,b]+conj(green.Gᴿ[e,m,b,a]))/2-X[e,m,a,b]
    end
    return X
end

function audit_channel!(f,name::Symbol,problem::P,cacheD::C,poD::O,
    X::ComplexField,defect::ComplexField,Rold::ComplexField,
    R::ComplexField,L::ComplexField,G::ComplexField) where {P,C,O}
    ne,nk,nb,_=size(X)
    q=problem.kernels.qᴷ[name]
    # Three nonzero terms of the exact same-input identity.
    firstterm=QN._selected_static_contraction(problem,cacheD,name,X,q,poD)::ComplexField
    Kdefect=QN._selected_static_contraction(problem,cacheD,name,defect,q,poD)::ComplexField
    Hdefect=primary_H(Kdefect,problem,cacheD,poD)
    delta=similar(X);herm=similar(X);anti=similar(X)
    gamma_mismatch=similar(X);gamma_identity_residual=similar(X)
    identity_residual=similar(X)
    @inbounds for m=1:nk,e=1:ne,b=1:nb,a=1:nb
        delta[e,m,a,b]=R[e,m,a,b]-Rold[e,m,a,b]
        gammaR=im*(R[e,m,a,b]-conj(R[e,m,b,a]))
        gammaK=im*(G[e,m,a,b]-L[e,m,a,b])
        gamma_mismatch[e,m,a,b]=gammaR-gammaK
        gamma_identity_residual[e,m,a,b]=gammaR-gammaK-Kdefect[e,m,a,b]
        identity_residual[e,m,a,b]=delta[e,m,a,b]-firstterm[e,m,a,b]-Hdefect[e,m,a,b]+.5im*Kdefect[e,m,a,b]
    end
    hermitian_parts!(herm,anti,delta)
    rec=Dict{String,Any}("delta_R"=>field_norms(problem,delta),"Hermitian_delta_R"=>field_norms(problem,herm),
        "antiHermitian_delta_R"=>field_norms(problem,anti),"K_HermGR_minus_HA"=>field_norms(problem,firstterm),
        "K_input_identity_defect"=>field_norms(problem,Kdefect),"H_K_input_identity_defect"=>field_norms(problem,Hdefect),
        "GammaR_minus_GammaK"=>field_norms(problem,gamma_mismatch),
        "GammaR_minus_GammaK_minus_K_AminusB"=>field_norms(problem,gamma_identity_residual),
        "full_discrete_delta_identity_residual"=>field_norms(problem,identity_residual),
        "old_R_scale"=>field_norms(problem,Rold),"new_R_scale"=>field_norms(problem,R))
    rec["delta_R_relative_to_old_unweighted"]=ratio(rec["delta_R"]["unweighted_Frobenius"],rec["old_R_scale"]["unweighted_Frobenius"])
    rec["identity_residual_relative_to_delta_unweighted"]=ratio(rec["full_discrete_delta_identity_residual"]["unweighted_Frobenius"],rec["delta_R"]["unweighted_Frobenius"])
    # Preserve the worst delta block and the worst identity-residual block;
    # enough to reconstruct signs/conjugations without full arrays.
    group=create_group(f,string(name))
    for (label,which) in (("worst_delta",rec["delta_R"]),("worst_identity",rec["full_discrete_delta_identity_residual"]))
        e=which["maximum_block_energy_index"];m=which["maximum_block_k_index"]
        subgroup=create_group(group,label)
        subgroup["energy_index"]=e;subgroup["k_index"]=m;subgroup["energy_eV"]=problem.grids.ε[e]*problem.scales.E₀_eV
        for (key,array) in (("old_R",Rold),("new_R",R),("delta",delta),
            ("firstterm",firstterm),("Kdefect",Kdefect),("Hdefect",Hdefect),("identity_residual",identity_residual))
            block=Matrix(view(array,e,m,:,:))
            subgroup[key*"_real"]=real.(block);subgroup[key*"_imag"]=imag.(block)
        end
    end
    return rec
end

function main()
    length(ARGS)==2 || error("usage: compare_static_maps.jl INPUT.toml NEW_OUTPUT_DIRECTORY")
    inputpath,out=abspath.(ARGS)
    ispath(out) && error("refusing overwrite: $out")
    cfg=TOML.parsefile(inputpath)
    Threads.nthreads()==cfg["julia_threads"] || error("Julia thread count differs")
    BLAS.set_num_threads(cfg["blas_threads"])
    BLAS.get_num_threads()==1 || error("BLAS must use1 thread")
    source=normpath(joinpath(dirname(inputpath),cfg["source_state"]))
    sourceinput=normpath(joinpath(dirname(inputpath),cfg["source_input"]))
    sourceidentity=normpath(joinpath(dirname(inputpath),cfg["source_identity"]))
    core=dirname(dirname(pathof(BN)))
    head=readchomp(`git -C $core rev-parse HEAD`)
    head==cfg["expected_core_head"] || error("wrong Core revision: $head")
    tree=readchomp(Cmd(["git","-C",core,"rev-parse","HEAD^{tree}"]))
    tree==cfg["expected_core_tree"] || error("wrong Core tree")
    isempty(readchomp(`git -C $core status --porcelain`)) || error("Core worktree dirty")
    srcsha=digest(source)
    srcsha==cfg["expected_source_sha256"] || error("wrong source SHA")
    c=TOML.parsefile(sourceinput)
    @assert c["alloy"] && c["voltage_mV"]==0 && c["temperature_K"]==70
    @assert (c["NE"],c["Nk"],c["Nb"],c["Nz"])==(30721,11,5,96)
    @assert cfg["trusted_min_eV"]==-.22 && cfg["trusted_max_eV"]==.82
    mkpath(out)
    cp(@__FILE__,joinpath(out,"compare_static_maps.jl"))
    cp(inputpath,joinpath(out,"input.toml"))
    cp(sourceinput,joinpath(out,"source-input.toml"))
    # Original identity stays unchanged; its hash is enough in this output.
    started=time()
    identity=Dict{String,Any}("start_utc"=>string(Dates.now(Dates.UTC))*"Z","julia_version"=>string(VERSION),
        "source_state_sha256"=>srcsha,"source_input_sha256"=>digest(sourceinput),
        "source_identity_sha256"=>digest(sourceidentity),"script_sha256"=>digest(@__FILE__),
        "input_sha256"=>digest(inputpath),"core_head"=>head,"core_tree"=>tree,"core_source"=>pathof(BN),
        "core_source_hashes"=>Dict(name=>digest(joinpath(core,"src/numerics/optimized",name)) for name in
            ("production_operators.jl","production_cache.jl","fft_parallel.jl","scattering_backend.jl")),
        "active_project_sha256"=>digest(Base.active_project()),"manifest_sha256"=>digest(joinpath(dirname(Base.active_project()),"Manifest.toml")),
        "scope"=>cfg["scope"],"acceptance_decision"=>"not_requested_not_evaluated",
        "harness_revision"=>"v2 concrete per-channel audit barriers; same scientific arithmetic")
    p=reference_parameters(V_period=0u"mV",Tᴸ=70u"K",Tᴸᴼ=70u"K",ΔV_alloy=.60u"eV",Ω₀=4.5168e-29u"m^3")
    n=NumericalParameters(N_z=c["Nz"],N_b=c["Nb"],P_basis=c["basis_periods"],
        E_min=c["energy_min_eV"]*u"eV",E_max=c["energy_max_eV"]*u"eV",N_E=c["NE"],
        M_E=c["trust_margin_meV"]*u"meV",k_max=c["kmax_per_nm"]*u"nm^-1",N_k=c["Nk"],
        N_φ=c["Nphi"],qz_max=c["qzmax_per_nm"]*u"nm^-1",N_qz=c["Nqz"],η_seed=c["seed_meV"]*u"meV")
    function algorithm(mode)
        AlgorithmOptions(kernel_build=:direct,embedding=:finite_chain,embedding_periods=c["embedding_periods"],
            mixing=:linear,energy_shift=:sparse_plan,hilbert=:fft,contraction=:dense_blas,
            elastic_retarded=mode,retarded_real_part=:kramers_kronig,self_energy_structure=:full,transverse_momentum=:resolved)
    end
    algH=algorithm(:hilbert);algD=algorithm(:direct_green)
    production(alg)=ProductionOptions(algorithms=alg,worker_count=4,hilbert_columns=cfg["hilbert_columns"],
        memory_budget_bytes=cfg["numerical_memory_budget_bytes"],checkpoint_every_scba=0,checkpoint_every_outer=0,
        verify_fft_roundoff=true,phase_timing=false)
    poH=production(algH);poD=production(algD)
    problem=build_configured_scattering_problem(physical=p,numerical=n,scales=ScaleSystem(),
        scattering=ScatteringOptions(alloy=true),algorithms=algH,kernel_options=ProductionKernelOptions()).problem
    @assert Set(problem.kernels.enabled)==Set((:LO,:acoustic,:impurity,:IFR,:alloy))
    # The honest native full-solver estimate remains its own guard, independent of
    # the externally enforced8GiB scope for this narrower diagnostic trajectory.
    estimate=estimate_production_memory(problem;options=poH)
    identity["native_estimate_peak_bytes"]=estimate.peak_bytes
    identity["native_memory_budget_bytes"]=poH.memory_budget_bytes
    cacheH=build_production_cache(problem;options=poH)
    cacheD=QN._validated_production_cache(cacheH,problem,poD)
    @assert cacheH.kernels===cacheD.kernels && cacheH.fft_hilbert_plan===cacheD.fft_hilbert_plan
    identity["shared_cache"]="native validation accepted both selectors; same kernels/FFT plan; only elastic_retarded differs"
    identity["old_algorithms"]=algorithm_manifest(algH)
    identity["new_algorithms"]=algorithm_manifest(algD)
    write_toml(joinpath(out,"identity.toml"),identity)
    println("cache built elapsed_s=",time()-started);flush(stdout)
    green=h5open(source,"r") do f
        @assert read(f["energy"])==problem.grids.ε && read(f["weight_energy"])==problem.grids.wᴱ
        @assert read(f["k"])==problem.grids.κ && read(f["weight_k"])==problem.grids.wᵏ
        @assert all(iszero,read(f["hartree"]))
        GreenState(readcomplex(f,"green/GR"),readcomplex(f,"green/GL"),readcomplex(f,"green/GG"),readcomplex(f,"green/A"),
            read(f["green/condition_number"]),read(f["green/dyson_scale"]))
    end
    sourcehash=Dict(string(name)=>array_digest(getfield(green,name)) for name in (:Gᴿ,:Gˡ,:Gᵍ,:A))
    @assert all(all(isfinite,getfield(green,name)) for name in (:Gᴿ,:Gˡ,:Gᵍ,:A))
    identity["source_green_hashes"]=sourcehash
    identity["hash_format"]="native little-endian ComplexF64 interleaved real/imag; Julia column-major E,k,a,b"
    write_toml(joinpath(out,"identity.toml"),identity)
    # Evaluate the actual selector twice on the identical unmodified GreenState.
    old,oldround=QN._scattering_candidate_production(problem,green,cacheH,poH;return_roundoff=true)
    old::Dict{Symbol,SelfEnergyFamily}
    oldR=Dict{Symbol,ComplexField}(name=>family.Σᴿ for (name,family) in old)
    oldhash=Dict(string(name)=>Dict(string(f)=>array_digest(getfield(family,f)) for f in (:Σᴿ,:Σˡ,:Σᵍ)) for (name,family) in old)
    old=nothing;GC.gc()
    println("old native map evaluated elapsed_s=",time()-started);flush(stdout)
    new,newround=QN._scattering_candidate_production(problem,green,cacheD,poD;return_roundoff=true)
    new::Dict{Symbol,SelfEnergyFamily}
    newhash=Dict(string(name)=>Dict(string(f)=>array_digest(getfield(family,f)) for f in (:Σᴿ,:Σˡ,:Σᵍ)) for (name,family) in new)
    # The separately evaluated Keldysh components must remain byte-identical.
    unchanged=Dict(string(name)=>Dict(string(f)=>oldhash[string(name)][string(f)]==newhash[string(name)][string(f)] for f in (:Σˡ,:Σᵍ)) for name in problem.kernels.enabled)
    lo_unchanged=Dict(string(f)=>oldhash["LO"][string(f)]==newhash["LO"][string(f)] for f in (:Σᴿ,:Σˡ,:Σᵍ))
    @assert all(all(values(x)) for x in values(unchanged)) && all(values(lo_unchanged))
    write_toml(joinpath(out,"selector-hashes.toml"),Dict("old"=>oldhash,"new"=>newhash,
        "unchanged_lesser_greater"=>unchanged,"unchanged_LO"=>lo_unchanged,"electron_electron"=>"not enabled in frozen all5 model; neither option changes its route",
        "old_fft_roundoff"=>oldround,"new_fft_roundoff"=>newround))
    println("new native map evaluated elapsed_s=",time()-started);flush(stdout)
    # Same input arithmetic in concrete-array helpers; no repair/projection.
    defect=input_identity_defect(green)::ComplexField
    inputnorm=field_norms(problem,defect)
    X=input_hilbert_defect(green,problem,cacheH,poH)::ComplexField
    evidence=Dict{String,Any}()
    h5open(joinpath(out,"witnesses.h5"),"w") do f
        f["schema"]="campaign133-static-map-control-witnesses-v1"
        for name in (:acoustic,:impurity,:IFR,:alloy)
            # Hoist dictionary lookups once per channel, before entering the
            # specialized audit barrier (never once per element).
            Rold=oldR[name]::ComplexField
            family=new[name]::SelfEnergyFamily
            timed=@timed audit_channel!(f,name,problem,cacheD,poD,X,defect,
                Rold,family.Σᴿ,family.Σˡ,family.Σᵍ)
            rec=timed.value
            rec["audit_seconds"]=timed.time
            rec["audit_allocated_bytes"]=timed.bytes
            rec["audit_gc_seconds"]=timed.gctime
            evidence[string(name)]=rec
            write_toml(joinpath(out,"channel-"*string(name)*".toml"),rec)
            println("channel=",name," audit_s=",timed.time," allocated_bytes=",timed.bytes,
                " gc_s=",timed.gctime," elapsed_s=",time()-started);flush(stdout)
            # Per-channel full-array temporaries are now out of scope.
            GC.gc()
        end
    end
    finalhash=Dict(string(name)=>array_digest(getfield(green,name)) for name in (:Gᴿ,:Gˡ,:Gᵍ,:A))
    @assert finalhash==sourcehash
    summary=Dict("status"=>"same_input_operator_diagnostic_completed_without_acceptance",
        "input_A_minus_n_minus_p"=>inputnorm,"source_green_unchanged"=>true,"channels"=>evidence,
        "identity"=>"DeltaR=K[HermGR-H_W^q A]+H_W^q K[A-B]-(i/2)K[A-B], B=-iGL+iGG",
        "weight_note"=>"full/trusted squared component norms use gs/(2pi)*wE*wk; all arrays dimensionless; trusted[-.22,.82]eV; no density normalization",
        "closure_scope"=>"direct static closure is distinct from finite-window Hilbert reconstruction; not baseline equivalence, not SCBA solution, not tail attribution",
        "embedding_control"=>"not performed; no embedding or central Dyson calls",
        "elapsed_s"=>time()-started,"finish_utc"=>string(Dates.now(Dates.UTC))*"Z")
    write_toml(joinpath(out,"summary.toml"),summary)
    println("completed elapsed_s=",time()-started);flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end
