# ROOT-owned postprocess only, no new physical fixed-point map or μ update.
const STARTED=time()
using TOML,SHA,HDF5,LinearAlgebra,Dates
include("moment_tail.jl")
module Loader end
resolve(c,k)=normpath(joinpath(@__DIR__,c[k]))
sha(p)=bytes2hex(open(sha256,p))
tick(c)=(time()-STARTED<c["wall_seconds"] || error("whole-process120s cap"))
function postprocess(input,out)
    c=TOML.parsefile(input);frozen=TOML.parsefile(joinpath(@__DIR__,"input.toml"))
    c["runtime_authorized"]===true || error("write-only; runtime_authorized=false")
    for (k,v) in frozen;k=="runtime_authorized" || isequal(c[k],v) || error("frozen option: $k");end
    for (key,d) in c["text_sha256"];sha(resolve(c,key))==d || error("text SHA: $key");end
    source=resolve(c,"actual_directory")
    for (file,d) in c["source_sha256"];tick(c);sha(joinpath(source,file))==d || error("source SHA: $file");end
    old=TOML.parsefile(joinpath(source,"input.toml"))
    if !isdefined(Loader,:source)
        Base.include(Loader,resolve(c,"old_adapter"))
        Base.invokelatest(Loader.gate,old)
        if !all(isdefined(Main,s) for s in (:EquilibriumGammaPilot,:MatsubaraAll5,:GammaDensityHermitianFreezeIO));Base.invokelatest(Loader.load_modules,old);end
        return Base.invokelatest(postprocess,input,out)
    end
    Loader.gate(old);tick(c) # exact Core/env/J1/2GiB/Swap0/CPU≤200/file32MiB gates
    sm=TOML.parsefile(joinpath(source,"summary.toml"));bind=TOML.parsefile(joinpath(source,"source-binding.toml"))
    sm["execution"]=="complete_bounded_diagnostic" && sm["inner_completed"]==32 && sm["completed_N"]==[128,256] || error("incomplete actual source")
    old["alpha"]==1 && old["mu_native"]==c["mu_native"] && sm["density_normalized"]===false || error("actual identity")
    model,newbind=Loader.source(old);tick(c)
    for key in ("model_proof","hH_sha256","Tplus_sha256","beta_native","Omega_native","mu_native","E0_eV","L0_m","spin")
        Loader.plain(getproperty(newbind,Symbol(key)))==bind[key] || error("model/unit binding: $key")
    end
    ispath(out) && error("new output required");mkpath(out)
    cp(input,joinpath(out,"input.toml"));cp(@__FILE__,joinpath(out,"run_moments.jl"))
    cp(resolve(c,"moment_module"),joinpath(out,"moment_tail.jl"))
    report=Dict{String,Any}("execution"=>"started","inner_maps"=>0,"mu_fits"=>0,"native_pass"=>false,
        "identity"=>"savedρbare + leading one-shot and NEW selfconsistent M3 asymptotic closure",
        "finite_N_tail_certificate"=>false,"source_sha256"=>c["source_sha256"],"text_sha256"=>c["text_sha256"])
    function save(name,x);Loader.toml_atomic(joinpath(out,name*".toml"),x,c);end
    try
        oracle=MatsubaraMoments.scalar_oracle(model.beta);save("scalar-oracle",oracle)
        oracle.moment_absolute_error<=1e-12 && all(r.M3_error<=r.bare_error+1e-12 for r in oracle.rows) || error("independent scalar moment oracle")
        computed=Dict{Int,Any}();densities=Dict{Int,Any}()
        for N in (128,256)
            tick(c);last=TOML.parsefile(joinpath(source,"N$(N)-i16.toml"))
            rho,G=h5open(joinpath(source,"N$(N)-latest.h5"),"r") do f
                read(f["schema"])=="campaign133-Mats-ALL5-aligned-raw-v1" && read(f["iteration"])==16 || error("state schema/iteration")
                read(f["mu_native"])==c["mu_native"] && TOML.parse(read(f["record_toml"]))==last || error("state μ/record")
                D=Main.GammaDensityHermitianFreezeIO;S=D.read_complex(f,"Sigma");P=D.read_complex(f,"Phi")
                S==P && last["raw_sigma_absolute"]==last["raw_sigma_relative"]==0 || error("not actual raw finite-N root")
                g=D.read_complex(f,"G");r=D.read_complex(f,"rho");w=read(f["omega_native"])
                size(S)==size(g)==(N,11,5,5) && size(r)==(11,5,5) || error("state axes")
                w==[(2n+1)*pi/model.beta for n in 0:N-1] && all(isfinite,g) && all(isfinite,r) || error("frequencies/finiteness")
                r,g
            end
            M=MatsubaraMoments;M.number(model,rho)==last["number"] || error("saved rho/number mismatch")
            answer=M.correct(model,rho,N;check=()->tick(c));classic=M.classic_M1(model,G,N)
            record=(;N,bare_number=M.number(model,rho),bare_pauli=M.pauli(rho),classic_M1_number=M.number(model,classic),
                oneshot_number=answer.oneshot_number,M3_number=answer.number,
                oneshot_pauli=answer.oneshot_pauli,M3_pauli=answer.pauli,
                feedback_absolute=norm(answer.feedback),correction_absolute=norm(answer.correction),
                t=answer.tails,eta=answer.eta,linear_condition=answer.linear_condition,
                linear_residual_absolute=answer.linear_residual_absolute,linear_residual_scale=answer.linear_residual_scale,linear_residual_relative=answer.linear_residual_relative,
                linear_weighted_absolute=answer.linear_weighted_absolute,linear_weighted_scale=answer.linear_weighted_scale,linear_weighted_relative=answer.linear_weighted_relative,
                conditional_linear_error_bound=answer.conditional_linear_error_bound,note=answer.note)
            save("N$(N)",record);computed[N]=record
            densities[N]=(;bare=rho,classic,oneshot=answer.oneshot,M3=answer.rho)
            Main.EquilibriumGammaPilot.atomic_h5(joinpath(out,"N$(N)-densities.h5"),c;expected_bytes=1024^2) do f
                f["schema"]="campaign133-Mats-density-M3-postprocess-v1";f["mu_native"]=c["mu_native"]
                for (name,X) in (("bare",rho),("classic_M1",classic),("oneshot",answer.oneshot),("M3",answer.rho),
                    ("feedback",answer.feedback),("interaction_S0",reshape(answer.S0,size(rho))),("B0",reshape(answer.B0,size(rho))))
                    Main.EquilibriumGammaPilot.write_complex(f,name,X)
                end
            end
            report["N$(N)"]=Loader.plain(record)
        end
        a,b=computed[128],computed[256]
        report["N256_minus_N128"]=Dict("bare"=>b.bare_number-a.bare_number,"oneshot"=>b.oneshot_number-a.oneshot_number,
            "M3"=>b.M3_number-a.M3_number,"classic_M1"=>b.classic_M1_number-a.classic_M1_number)
        report["matrix_Frobenius_N256_minus_N128"]=Dict(string(k)=>norm(getproperty(densities[256],k)-getproperty(densities[128],k)) for k in (:bare,:classic,:oneshot,:M3))
        report["execution"]="completed_postprocess_only"
    catch err
        report["execution"]="stopped";report["reason"]=sprint(showerror,err);rethrow()
    finally
        report["elapsed_seconds_including_imports"]=time()-STARTED;report["maxrss_bytes"]=Sys.maxrss()
        save("summary",report)
    end
    report
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("INPUT NEW_OUTPUT");postprocess(abspath(ARGS[1]),abspath(ARGS[2]))
end
