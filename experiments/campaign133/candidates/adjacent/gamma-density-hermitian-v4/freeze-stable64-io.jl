# Write-only preparation: root-owned bounded read-only HDF5 digest extraction.
module GammaDensityHermitianFreezeIO
using HDF5, SHA, TOML, Dates
sha_file(p)=bytes2hex(open(sha256,p))
# EXACT V3 logical-array digest expression; no h5py-axis reinterpretation.
sha_array(x)=bytes2hex(sha256(reinterpret(UInt8,vec(x))))
read_complex(g,n)=complex.(read(g[n*"_real"]),read(g[n*"_imag"]))
function pi_geometry(E)
    N=length(E); N>=2 && all(isfinite,E)||error("invalid E")
    a,b=Float64(first(E)),Float64(last(E)); width=b-a
    isfinite(width)&&width>0||error("invalid E span")
    canonical=range(a,b;length=N); h=width/(N-1); nd,sd=0.0,0.0
    for i in eachindex(E)
        E[i]==canonical[i]||error("not the native range geometry")
        if i>1
            step=Float64(E[i]-E[i-1]); step>0||error("nonmonotone E")
            sd=max(sd,abs(step-h))
        end
        nd=max(nd,abs(Float64(E[i])-muladd(i-1,h,a)))
    end
    sd<=8eps(Float64)*width||error("nonuniform native geometry")
    return Dict("origin"=>a,"spacing"=>h,"max_node_deviation"=>nd,"relative_node_deviation"=>nd/h)
end
function digest_one(f,path;complex_array=false)
    x=complex_array ? read_complex(f,path) : read(f[path])
    all(isfinite,x)||error("nonfinite: $path")
    record=Dict("sha256"=>sha_array(x),"logical_shape"=>collect(size(x)),"eltype"=>string(eltype(x)))
    geometry=path=="grid_ε" ? pi_geometry(x) : nothing
    x=nothing; GC.gc(true)
    return record,geometry
end
function main(input,output)
    c=TOML.parsefile(input); base=dirname(abspath(input)); ispath(output)&&error("output exists")
    Threads.nthreads()==1||error("requires Julia1")
    started=time(); paths=Dict{String,String}(); hashes=Dict{String,String}()
    checktime()=(time()-started<=105||error("freeze-IO time boundary"))
    project=Base.active_project(); manifest=joinpath(dirname(project),"Manifest.toml")
    sha_file(project)==c["project_sha256"]&&sha_file(manifest)==c["manifest_sha256"]||error("oracle environment")
    hashes["oracle_project"]=sha_file(project); hashes["oracle_manifest"]=sha_file(manifest)
    for key in ("terminal_h5","model_h5","source_input","source_identity","source_writer","map_summary","numerical_parent","oracle_script",
                "checkpoint","spectral_selected","spectral_writer","source_ledger","source_summary")
        checktime(); p=normpath(joinpath(base,c[key])); d=sha_file(p)
        d==c[key*"_sha256"]||error("whole-file SHA: $key"); paths[key]=p; hashes[key]=d
    end
    owner=dirname(realpath(paths["terminal_h5"]))
    all(dirname(realpath(paths[k]))==owner for k in ("model_h5","source_input","source_identity","source_writer","map_summary",
        "checkpoint","spectral_selected","spectral_writer","source_ledger","source_summary"))||error("attempt ownership")
    ident=TOML.parsefile(paths["source_identity"]); src=TOML.parsefile(paths["source_input"]); mp=TOML.parsefile(paths["map_summary"])
    ident["core_head"]==c["source_core_head"]&&ident["core_tree"]==c["source_core_tree"]&&isempty(ident["core_dirty"])||error("Core identity")
    ident["input_sha256"]==hashes["source_input"]&&ident["source_sha256"]["gamma_functions.jl"]==hashes["source_writer"]||error("writer binding")
    ident["source_sha256"]["stable_spectral.jl"]==hashes["spectral_writer"]||error("stable writer binding")
    src["project_AA_trial"]===false&&src["spectral_operator"]=="native_covariance_exact_plan_hlocal_all9"&&mp["call"]==src["pi_calls_total"]==64||error("wrong snapshot")
    ledger=TOML.parsefile(paths["source_ledger"]); source_summary=TOML.parsefile(paths["source_summary"])
    ledger==mp["spectral_operator_ledger"]&&source_summary["execution"]=="complete"&&source_summary["pi_calls"]==64&&source_summary["reduced"]["last_map"]==mp||error("ledger/terminal source")
    shape=(src["NE"],src["Nk"],src["Nb"],src["Nb"])
    shape==(c["expected_NE"],c["expected_Nk"],c["expected_Nb"],c["expected_Nb"])||error("shape freeze")
    mkpath(output); cp(input,joinpath(output,"input-pre-digest.toml")); cp(@__FILE__,joinpath(output,"freeze-stable64-io.jl"))
    result=Dict{String,Any}("status"=>"running","source_file_sha256"=>hashes,"input_sha256"=>sha_file(input),
        "helper_sha256"=>sha_file(@__FILE__),"core_head"=>ident["core_head"],"core_tree"=>ident["core_tree"],
        "array_sha256"=>Dict{String,String}(),"reference_A_k_sha256"=>Dict{String,String}(),"array_metadata"=>Dict{String,Any}())
    model_raw=nothing
    try
        h5open(paths["model_h5"],"r") do f
            read(f["schema"])=="campaign133-gamma-model-v2"||error("model schema")
            alg=TOML.parse(read(f["algorithms_toml"]))
            alg["hilbert"]=="product_integration"&&alg["embedding"]=="finite_chain"&&alg["embedding_periods"]==4||error("algorithms")
            for (name,path,iscomplex) in (("energy","grid_ε",false),("weight_energy","grid_wᴱ",false),
                ("momentum","grid_κ",false),("weight_momentum","grid_wᵏ",false),("h_raw","h_fixedU0",true),
                ("Phi","basis_Φ",true),("chi","basis_χ",true),("H0","basis_H₀",true),("Z","basis_Z",true),
                ("Minv","basis_M⁻¹",true),("Tplus","basis_T₊",true),("Tminus","basis_T₋",true))
                checktime(); record,geometry=digest_one(f,path;complex_array=iscomplex)
                expected=name in ("energy","weight_energy") ? (shape[1],) : name in ("momentum","weight_momentum") ?
                    (shape[2],) : name=="h_raw" ? shape[2:4] : name in ("H0","Z","Minv","Tplus","Tminus") ? shape[3:4] : nothing
                expected===nothing||record["logical_shape"]==collect(expected)||error("model logical shape: $name")
                if name=="h_raw"; result["h_raw_sha256"]=record["sha256"]
                else; result["array_sha256"][name]=record["sha256"]; end
                result["array_metadata"][name]=record
                if geometry!==nothing; result["PI_geometry"]=geometry; end
            end
            model_raw=read_complex(f,"h_fixedU0") # 11x5x5 only, for exact policy crosscheck.
            tp=read_complex(f,"basis_T₊"); tm=read_complex(f,"basis_T₋")
            all(tm[a,b]==conj(tp[b,a]) for a in 1:5,b in 1:5)||error("couplings not EXACT adjoints")
        end
        h5open(paths["spectral_selected"],"r") do f
            read(f["schema"])=="campaign133-stableA-spectral-selected-v1"||error("spectral selected schema")
            isequal(read_complex(f,"h_raw"),model_raw)||error("spectral/model h_raw differ")
            for (name,path) in (("energy","energy"),("weight_energy","weight_energy"),("weight_momentum","weight_k"))
                sha_array(read(f[path]))==result["array_sha256"][name]||error("spectral grid binding")
            end
            read(f["energy_scale_eV"])==mp["absolute_E0_eV"]&&sha_array(read_complex(f,"Tplus"))==result["array_sha256"]["Tplus"]||error("spectral units/coupling")
            TOML.parse(read(f["ledger_toml"]))==ledger||error("saved ledger")
            hh=read_complex(f,"h_H"); size(hh)==shape[2:4]&&eltype(hh)==ComplexF64&&all(isfinite,hh)||error("h_H shape/type/finite")
            all(hh[k,a,b]==conj(hh[k,b,a]) for k in 1:shape[2],a in 1:5,b in 1:5)||error("h_H not EXACT Hermitian")
            result["array_sha256"]["hcell"]=sha_array(hh)
            result["array_metadata"]["hcell"]=Dict("logical_shape"=>collect(size(hh)),"eltype"=>string(eltype(hh)),"source"=>"actual saved h_H")
        end
        h5open(paths["terminal_h5"],"r") do f
            read(f["schema"])=="campaign133-gamma-reduced-terminal-v1"||error("terminal schema")
            mu=Float64(read(f["mu"])); isfinite(mu)&&mu==c["mu_dimensionless"]==mp["mu"]["mu"]||error("mu binding")
            result["mu_dimensionless"]=mu
            for n in ("Gamma_real","Gamma_imag","green/A_real","green/A_imag"); size(f[n])==shape||error("logical shape: $n"); end
            checktime(); record,_=digest_one(f,"Gamma";complex_array=true)
            result["array_sha256"]["Gamma"]=record["sha256"]; result["array_metadata"]["Gamma"]=record
            h5open(paths["checkpoint"],"r") do checkpoint
                read(checkpoint["schema"])=="campaign133-gamma-checkpoint-v1"&&read(checkpoint["pi_calls_total"])==64&&read(checkpoint["identity"])==src["identity"]||error("checkpoint identity")
                for n in ("Gamma_real","Gamma_imag")
                    size(checkpoint[n])==shape&&eltype(checkpoint[n])==eltype(f[n])==Float64||error("checkpoint axes/type")
                    for first in 1:1024:shape[1]
                        checktime(); rows=first:min(first+1023,shape[1])
                        isequal(f[n][rows,:,:,:],checkpoint[n][rows,:,:,:])||error("terminal/checkpoint Gamma differs")
                    end
                end
            end
            result["checkpoint_terminal_Gamma_exact_equal"]=true; result["compared_complex_entries"]=prod(shape)
            for k in 1:src["Nk"]
                checktime(); A=complex.(f["green/A_real"][:,k,:,:],f["green/A_imag"][:,k,:,:])
                all(isfinite,A)||error("nonfinite A slice"); result["reference_A_k_sha256"][string(k)]=sha_array(A)
                A=nothing; GC.gc(true)
            end
        end
        checktime(); result["status"]="completed"; result["scope"]="digests/geometry/exact source equality only; not PD, density or SCBA acceptance"
        result["global_PD_status"]="not evaluated by freeze-IO; completed ALL certificate required by actual Matsubara gate"
    catch err
        result["status"]="failed"; result["error"]=sprint(showerror,err); rethrow()
    finally
        result["elapsed_seconds_including_whole_file_SHA"]=time()-started
        result["finish_utc"]=string(Dates.now(Dates.UTC))
        open(joinpath(output,"digest.toml"),"w") do io; TOML.print(io,result); end
    end
    sum(filesize(joinpath(output,n)) for n in readdir(output))<=32*1024^2||error("output cap")
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2||error("INPUT_CONT128 NEW_OUTPUT")
    GammaDensityHermitianFreezeIO.main(abspath(ARGS[1]),abspath(ARGS[2]))
end
