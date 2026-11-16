# WRITE-ONLY: ROOT runs preflight first, then chooses a bounded resource gap.
using QCLNEGF, LinearAlgebra, HDF5, SHA, TOML, Dates
include(joinpath(@__DIR__,"tiny_tests.jl")) # source-free tiny tests BEFORE any large read
const CT=CauchyPanelTree
const PN=QCLNEGF.QCLNumerics
const PIPlan=PN._ProductionPIFFTHilbertPlan
const PIFunction=PN._production_pi_fft_hilbert_transform
sha_file(p)=bytes2hex(open(sha256,p))
sha_array(A)=bytes2hex(sha256(reinterpret(UInt8,vec(A))))
matrix_record(A)=Dict("real"=>vec(real.(A)),"imag"=>vec(imag.(A)),"shape"=>collect(size(A)))
# Same local-R extraction definition as source_selected_k.jl; source SHA is checked ONCE before opening.
function selected_gamma(f,k)
    R=zeros(ComplexF64,30721,5,5)
    for name in ("IFR","LO","acoustic","alloy","impurity")
        group=f["selfenergy/"*name]
        for suffix in ("_real","_imag"); size(group["retarded"*suffix])==(30721,11,5,5) || error("source axes"); end
        R .+= complex.(group["retarded_real"][:,k,:,:],group["retarded_imag"][:,k,:,:])
    end
    G=similar(R)
    for b in 1:5,a in 1:5,e in 1:30721; G[e,a,b]=im*(R[e,a,b]-conj(R[e,b,a])); end
    CT.finiteall(G) && CT.finiteall(R) || error("nonfinite source slice")
    return G
end
# Independent nonuniform positive-hat divided difference; cached F at ALL source+ghost knots.
function independent_hat(E,G,ghosts,z,realpv,cache)
    weights=get!(cache,(real(z),imag(z))) do
      setprecision(BigFloat,192) do
        x=BigFloat[ghosts[1];E;ghosts[2]]; zz=realpv ? BigFloat(z) : Complex{BigFloat}(z)
        F(t)=iszero(t) ? zero(t) : t*(realpv ? log(abs(t)) : log(t))
        values=F.(zz.-x); w=Vector{ComplexF64}(undef,length(E))
        for j in eachindex(E)
            w[j]=ComplexF64(((values[j]-values[j+1])/(x[j+1]-x[j])-
                (values[j+1]-values[j+2])/(x[j+2]-x[j+1]))/(2big(pi)))
        end; w
      end
    end
    flat=reshape(G,length(E),25)
    value=reshape(transpose(weights)*flat,5,5) # transpose, NOT adjoint
    reversevalue=reshape(transpose(view(weights,length(E):-1:1))*view(flat,length(E):-1:1,:),5,5)
    sensitivity=norm(value-reversevalue)
    if realpv
        x=[ghosts[1];E;ghosts[2]]; j=clamp(searchsortedlast(x,z),1,length(x)-1)
        theta=(z-x[j])/(x[j+1]-x[j])
        l=j==1 ? zeros(ComplexF64,5,5) : G[j-1,:,:]
        r=j+1==length(x) ? zeros(ComplexF64,5,5) : G[j,:,:]
        value .-=0.5im.*((1-theta).*l.+theta.*r)
    end
    return value,sensitivity
end
function main()
    if ARGS==["--preflight"]
        println("source_free_preflight=true tree=",CT.PanelTree," plan=",PIPlan," function=",PIFunction); return
    end
    length(ARGS)==1 || error("usage: run_frozen.jl --preflight | NEW_OUTPUT_DIRECTORY")
    cfg=TOML.parsefile(joinpath(@__DIR__,"frozen-input.toml")); BLAS.set_num_threads(1)
    Threads.nthreads()==1 || error("Julia1 required")
    cfg["k_indices"]==[1,6,11] && cfg["tau_far"]==[1e-9,1e-11] || error("probe identity controls")
    length(unique(cfg["source_indices"]))==16 && all(j->1<=j<=30721,cfg["source_indices"]) || error("declared indices")
    isfinite(cfg["soft_stop_seconds"]) && cfg["soft_stop_seconds"]>0 || error("time control")
    out=abspath(ARGS[1]); ispath(out) && error("new output directory required"); mkpath(out)
    started=time(); source=normpath(joinpath(@__DIR__,cfg["source"])); records=Dict{String,Any}[]
    summary=Dict{String,Any}("id"=>cfg["id"],"status"=>"running","records"=>records,
        "formal_whole_certificate"=>false,"density_SCBA_adaptive_economy"=>"not_measured")
    function checkpoint(stage)
        summary["stage"]=stage; summary["elapsed_s"]=time()-started; summary["maxrss_bytes"]=Sys.maxrss()
        open(joinpath(out,"summary.toml"),"w") do io; TOML.print(io,summary); end
        bytes=sum(filesize(joinpath(d,n)) for (d,_,fs) in walkdir(out) for n in fs;init=0)
        bytes<=cfg["maximum_output_bytes"] || error("output cap")
        time()-started<=cfg["soft_stop_seconds"] || error("soft time cap")
        println("stage=",stage," elapsed_s=",time()-started); flush(stdout)
    end
    try
        for name in ("CauchyPanelTree.jl","tiny_tests.jl")
            sha_file(joinpath(@__DIR__,name))==cfg["source_sha256"][name] || error("frozen prototype SHA")
        end
        for name in ("run_frozen.jl","frozen-input.toml","frozen-protocol.md","CauchyPanelTree.jl","tiny_tests.jl")
            cp(joinpath(@__DIR__,name),joinpath(out,name))
        end
        summary["source_code_sha256"]=Dict(n=>sha_file(joinpath(@__DIR__,n)) for n in ("run_frozen.jl","frozen-input.toml","CauchyPanelTree.jl","tiny_tests.jl"))
        summary["julia_version"]=string(VERSION); summary["Julia_threads"]=Threads.nthreads(); summary["BLAS_threads"]=BLAS.get_num_threads()
        summary["resource_enforcement"]="external ROOT scope; script only soft-time/output checks"
        core=dirname(dirname(pathof(QCLNEGF)))
        readchomp(`git -C $core rev-parse HEAD`)==cfg["core_head"] || error("Core HEAD")
        isempty(readchomp(`git -C $core status --porcelain`)) || error("dirty accepted Core")
        summary["core_import_path"]=pathof(QCLNEGF); summary["core_head"]=cfg["core_head"]; summary["project_sha256"]=sha_file(Base.active_project())
        manifest=joinpath(dirname(Base.active_project()),"Manifest.toml"); summary["manifest_sha256"]=sha_file(manifest)
        priorpath=normpath(joinpath(@__DIR__,cfg["prior_summary"])); sha_file(priorpath)==cfg["prior_summary_sha256"] || error("prior SHA")
        prior=TOML.parsefile(priorpath); dir=dirname(source)
        for (file,key) in (("input.toml","input_sha256"),("identity.toml","identity_sha256"),("run_inner.jl","source_writer_sha256"),("markers/final_markers.toml","marker_sha256"))
            expected=file=="input.toml" ? cfg["source_input_sha256"] : prior[key]
            sha_file(joinpath(dir,file))==expected || error("writer provenance: "*file)
        end
        ident=TOML.parsefile(joinpath(dir,"identity.toml")); ident["input_sha256"]==cfg["source_input_sha256"] || error("writer input")
        ident["script_sha256"]==prior["source_writer_sha256"] || error("writer identity")
        TOML.parsefile(joinpath(dir,"markers/final_markers.toml"))["equilibrium_mu_eV"]==cfg["mu_eV"] || error("source mu")
        sha_file(source)==cfg["source_file_sha256"] || error("source SHA before scientific reads")
        prior["frozen_local_gamma_sha256"]==cfg["aggregate_gamma_reference_sha256"] || error("aggregate Gamma reference")
        summary["aggregate_gamma_digest_scope"]="prior reference; not recomputed from selected k"
        checkpoint("identity")
        h5open(source,"r") do f
            read(f["schema"])=="campaign133-inner-raw-v3" && read(f["energy_scale_eV"])==0.1 || error("source schema/scales")
            all(iszero,read(f["hartree"])) || error("Hartree source differs")
            for (stored,key) in (("Phi","Φ"),("chi","χ"),("H0","H₀"),("Z","Z"),("Minv","M⁻¹"),("Tplus","T₊"),("Tminus","T₋"))
                expected=stored in ("Phi","chi") ? (96,5) : (5,5)
                for suffix in ("_real","_imag"); size(f["basis/"*stored*suffix])==expected || error("basis dimensions"); end
                A=complex.(read(f["basis/"*stored*"_real"]),read(f["basis/"*stored*"_imag"]))
                sha_array(A)==prior["basis_source_and_rebuild_comparison"][key]["source_sha256"] || error("basis digest")
            end
            E=read(f["energy"]); w=read(f["weight_energy"]); length(E)==30721 || error("energy count")
            ghosts=(E[1]-(E[2]-E[1]),E[end]+(E[end]-E[end-1]))
            weight_cache=Dict{Tuple{Float64,Float64},Vector{ComplexF64}}()
            for k in cfg["k_indices"]
                checkpoint("load_k"*string(k)); G=selected_gamma(f,k); scale=maximum(norm(view(G,e,:,:)) for e in eachindex(E))
                allowance=64eps(Float64)*5*max(scale,eps(Float64)) # explicit eigensolve/backward witness allowance
                tree=CT.build_tree(E,G;ghosts,pmax=24,leaf_panels=8,hermitian_atol=allowance,psd_atol=allowance)
                plan=PIPlan(E,25;column_chunk=1,worker_count=1); fft,fft_round=PIFunction(reshape(G,length(E),1,5,5),E,w,plan)
                meta=Dict(string(a)=>b for (a,b) in pairs(plan.geometry)); root=tree.nodes[tree.root]
                info=Dict("k"=>k,"gamma_k_sha256"=>sha_array(G),"energy_sha256"=>sha_array(E),"ghosts"=>collect(ghosts),
                    "FFT_geometry"=>meta,"FFT_roundoff"=>fft_round,"input_witness"=>Dict(string(a)=>b for (a,b) in pairs(tree.input_witness)),
                    "moment_reverse_Frobenius"=>norm(root.moments-root.reverse_moments),"geometry"=>Dict(string(a)=>b for (a,b) in pairs(tree.geometry_witness)))
                fg=[plan.geometry.origin-plan.geometry.spacing,muladd(length(E),plan.geometry.spacing,plan.geometry.origin)]
                info["FFT_ideal_ghosts"]=fg; info["ghost_geometry_difference"]=maximum(abs.(fg.-collect(ghosts)))
                summary["k"*string(k)]=info; info["allowance_definition"]="64eps*Nb*max(block_Frobenius_scale,eps); numeric witness only"
                checkpoint("built_tree_fft_k"*string(k))
                targets=[(E[j]+E[j+1])/2 for j in cfg["offgrid_midpoint_indices"]]
                append!(targets,[(ghosts[1]+E[1])/2]); zs=[ComplexF64(z) for z in targets]
                push!(zs,ComplexF64(cfg["mu_eV"]/0.1,cfg["complex_imag_native"]))
                refs=[independent_hat(E,G,ghosts,imag(z)==0 ? real(z) : z,imag(z)==0,weight_cache) for z in zs]
                summary["independent_reference_weight_cache_entries"]=length(weight_cache)
                for tau in cfg["tau_far"]
                    queries=[(ComplexF64(E[j]),Matrix(fft[j,1,:,:])-0.5im*G[j,:,:],"FFT_node_"*string(j),0.0) for j in cfg["source_indices"]]
                    append!(queries,[(zs[j],refs[j][1],"independent_offgrid_"*string(j),refs[j][2]) for j in eachindex(zs)])
                    for (z,ref,label,ref_round) in queries
                        checkpoint("query_k"*string(k)); t=time()
                        value,ledger=imag(z)==0 ? CT.real_boundary(tree,real(z);tau_far=tau) : CT.complex_cauchy(tree,z;tau_far=tau)
                        diff=norm(value-ref); allowed=ledger.far_truncation_estimate+max(1e-11,1e-10*norm(ref)); passed=diff<=allowed
                        push!(records,Dict("k"=>k,"tau"=>tau,"label"=>label,"z_real"=>real(z),"z_imag"=>imag(z),"elapsed_s"=>time()-t,
                            "difference_Frobenius"=>diff,"diagnostic_allowance"=>allowed,"diagnostic_pass"=>passed,"reference_reverse_witness"=>ref_round,
                            "value"=>matrix_record(value),"reference"=>matrix_record(ref),"ledger"=>Dict(string(a)=>b for (a,b) in pairs(ledger))))
                        checkpoint("saved_k"*string(k)); passed || error("diagnostic comparison failed: "*label)
                    end
                end
                tree=plan=G=fft=refs=nothing; GC.gc(); checkpoint("complete_k"*string(k))
            end
        end
        summary["status"]="completed_diagnostic_only"
    catch err
        summary["status"]="failed_or_budget_stop"; summary["error"]=sprint(showerror,err); rethrow()
    finally
        summary["elapsed_s"]=time()-started; summary["maxrss_bytes"]=Sys.maxrss()
        open(joinpath(out,"summary.toml"),"w") do io; TOML.print(io,summary); end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
