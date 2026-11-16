using QCLNEGF, Unitful, LinearAlgebra, HDF5, TOML, SHA, Dates
const BN = QCLNEGF
const started = time()
length(ARGS) == 2 || error("usage: one_fresh_fd_map.jl INPUT.toml NEW_OUTPUT_DIRECTORY")
const config_path = abspath(ARGS[1])
const output = abspath(ARGS[2])
ispath(output) && error("refusing to overwrite output: $output")
mkpath(output)
# Archive exactly the executable source before hashing its archived copy.
cp(@__FILE__, joinpath(output, "one_fresh_fd_map.jl"))
cp(config_path, joinpath(output, "input.toml"))
const cfg = TOML.parsefile(config_path)
BLAS.set_num_threads(cfg["blas_threads"])
Threads.nthreads() == cfg["julia_threads"] == 1 || error("one Julia thread required")
BLAS.get_num_threads() == 1 || error("one BLAS thread required")
write_toml(name, value) = open(io -> TOML.print(io, value), joinpath(output, name), "w")
digest(path) = bytes2hex(open(sha256, path))
path_from_config(name) = normpath(joinpath(dirname(config_path), cfg[name]))
ratio(n,d) = iszero(d) ? (iszero(n) ? 0.0 : Inf) : n/d
stamp() = string(Dates.now(Dates.UTC))*"Z"
function progress(message)
    println(stamp(), " ", message, " elapsed_s=", time()-started)
    flush(stdout)
end
fermi_pair(x) = x >= 0 ? (let z=exp(-x); (z/(1+z),1/(1+z)) end) :
                             (let z=exp(x); (1/(1+z),z/(1+z)) end)
read_complex(file, path) = complex.(read(file[path*"_real"]), read(file[path*"_imag"]))
norm_record(sq, scale_sq) = Dict("absolute"=>sqrt(sq), "scale"=>sqrt(scale_sq),
    "relative"=>sqrt(ratio(sq,scale_sq)))
function sandwich!(result, GR, matrix, temporary)
    mul!(temporary, GR, matrix)
    mul!(result, temporary, adjoint(GR))
    result
end

"""Same component scale as Core markers plus independent quadrature-weighted data.
The reference is source SigmaR for R and FD-repaired source GammaR for </>.
No candidate Keldysh component is normalized or projected.
"""
function family_audit(name, candidate, sourceR, green, f, holes, problem)
    NE,Nk,Nb,_=size(green.A)
    g=problem.grids; gs=problem.physical.g_s
    domains=("full","trusted")
    # component x (difference squared, candidate squared, reference squared)
    unweighted=zeros(3,3)
    weighted=Dict(domain=>zeros(3,3) for domain in domains)
    # KMS sq, KMS scale sq, sn sq, sp sq, Gamma identity sq, Gamma scale sq,
    # C0, CE_mu, CE_zero, activity, local abs collision, Im C0, Im CE_mu.
    moments=Dict(domain=>zeros(13) for domain in domains)
    mu_scaled = CURRENT_MU[]
    for e=1:NE, m=1:Nk
        incoming=0.0im; outgoing=0.0im
        kms2=kms_scale2=sn2=sp2=gid2=gscale2=0.0
        component=zeros(3,3)
        for b=1:Nb,a=1:Nb
            gammaS=im*(sourceR[e,m,a,b]-conj(sourceR[e,m,b,a]))
            oldR=sourceR[e,m,a,b]
            oldL=im*f[e]*gammaS
            oldG=-im*holes[e]*gammaS
            newR=candidate.Σᴿ[e,m,a,b]
            newL=candidate.Σˡ[e,m,a,b]
            newG=candidate.Σᵍ[e,m,a,b]
            for (j,new,old) in ((1,newR,oldR),(2,newL,oldL),(3,newG,oldG))
                component[j,1]+=abs2(new-old)
                component[j,2]+=abs2(new)
                component[j,3]+=abs2(old)
            end
            sn=-im*newL; sp=im*newG
            gammaK=sn+sp
            gammaR=im*(newR-conj(candidate.Σᴿ[e,m,b,a]))
            D=holes[e]*sn-f[e]*sp
            kms2+=abs2(D); kms_scale2+=abs2(sn)+abs2(sp)
            sn2+=abs2(sn);sp2+=abs2(sp)
            gid2+=abs2(gammaK-gammaR);gscale2+=abs2(gammaK)
            incoming+=sn*holes[e]*green.A[e,m,b,a]
            outgoing+=sp*f[e]*green.A[e,m,b,a]
        end
        unweighted .+= component
        collision=incoming-outgoing
        weight=gs/(2pi)*g.wᴱ[e]*g.wᵏ[m]
        trusted=cfg["trusted_min_eV"] <= g.ε[e]*problem.scales.E₀_eV <= cfg["trusted_max_eV"]
        for domain in domains
            domain=="trusted" && !trusted && continue
            weighted[domain] .+= weight .* component
            moments[domain] .+= weight .* [kms2,kms_scale2,sn2,sp2,gid2,gscale2,
                real(collision),(g.ε[e]-mu_scaled)*real(collision),g.ε[e]*real(collision),
                abs(incoming)+abs(outgoing),abs(collision),imag(collision),
                (g.ε[e]-mu_scaled)*imag(collision)]
        end
    end
    components=Dict{String,Any}()
    csvrows=Vector{Vector{Any}}()
    for (j,component) in enumerate(("retarded","lesser","greater"))
        rec=Dict{String,Any}("unweighted"=>norm_record(unweighted[j,1],unweighted[j,2]),
            "unweighted_reference_norm"=>sqrt(unweighted[j,3]),
            "core_marker_relative_with_1em14_floor"=>sqrt(unweighted[j,1])/(sqrt(unweighted[j,2])+1e-14),
            "reference_kind"=>(j==1 ? "unchanged archived SigmaR" : "FD repair from archived GammaR"))
        for domain in domains
            v=weighted[domain]
            rec[domain*"_weighted"]=norm_record(v[j,1],v[j,2])
            rec[domain*"_weighted_reference_norm"]=sqrt(v[j,3])
        end
        components[component]=rec
        push!(csvrows,Any[string(name),component,sqrt(unweighted[j,1]),sqrt(unweighted[j,2]),
            sqrt(unweighted[j,3]),rec["core_marker_relative_with_1em14_floor"],
            rec["full_weighted"]["absolute"],rec["full_weighted"]["scale"],rec["full_weighted"]["relative"]])
    end
    summary=Dict{String,Any}("channel"=>string(name),"components"=>components)
    E0=problem.scales.E₀_eV; L0=problem.scales.L₀_m
    hbar_eVs=ustrip(u"eV*s",BN.CODATA.ħₑᵥ)
    chargeC=ustrip(u"C",BN.CODATA.e)
    rate_scale=E0/(hbar_eVs*L0^2)
    power_scale=E0^2*chargeC/(hbar_eVs*L0^2)
    omega=ustrip(u"eV",problem.physical.ħωᴸᴼ)/E0
    for domain in domains
        v=moments[domain]
        summary[domain]=Dict("KMS"=>norm_record(v[1],v[2]),
            "sn_weighted_norm"=>sqrt(v[3]),"sp_weighted_norm"=>sqrt(v[4]),
            "GammaK_vs_GammaR"=>norm_record(v[5],v[6]),
            "C0_signed_solver"=>v[7],"CE_mu_signed_solver"=>v[8],"CE_zero_signed_solver"=>v[9],
            "activity_solver"=>v[10],"local_collision_absolute_solver"=>v[11],
            "C0_imaginary_solver"=>v[12],"CE_mu_imaginary_solver"=>v[13],
            "C0_over_activity"=>ratio(v[7],v[10]),"CE_mu_over_LO_activity"=>ratio(v[8],omega*v[10]),
            "C0_per_m2_s"=>v[7]*rate_scale,"CE_mu_W_per_m2"=>v[8]*power_scale,
            "activity_per_m2_s"=>v[10]*rate_scale)
    end
    write_toml("channel-"*string(name)*".toml",summary)
    return summary,csvrows
end

const CURRENT_MU = Ref(NaN)
function decomposition(problem, green, source_totalR, candidate, f, holes, target)
    NE,Nk,Nb,_=size(green.A);g=problem.grids;spin=problem.physical.g_s/(2pi)
    labels=("occupied_raw_error","hole_raw_error","occupied_broadening_part","hole_broadening_part",
        "KMS_response","occupied_spectral_identity_correction","hole_spectral_identity_correction",
        "occupied_recombination_error","hole_recombination_error","source_spectral_identity_error",
        "occupied_input","hole_input","spectral_input")
    squares=Dict(domain=>Dict(label=>zeros(2) for label in labels) for domain in ("full","trusted"))
    number_labels=("FD_number","FD_holes","candidate_raw_number","candidate_raw_holes",
        "broadening_number_part","KMS_number_part","spectral_identity_number_correction",
        "current_raw_number","number_recombination_error")
    number=Dict(domain=>Dict(label=>zeros(3) for label in number_labels) for domain in ("full","trusted"))
    # Reuse small blocks only; never materialize full candidate Green products.
    GR,A,gs,gc,sn,sp,D,dg,B,T,As,nc,pc,tmp=(zeros(ComplexF64,Nb,Nb) for i=1:14)
    for e=1:NE,m=1:Nk
        for b=1:Nb,a=1:Nb
            GR[a,b]=green.Gᴿ[e,m,a,b];A[a,b]=green.A[e,m,a,b]
            gs[a,b]=im*(source_totalR[e,m,a,b]-conj(source_totalR[e,m,b,a]))
            sn[a,b]=-im*candidate.Σˡ[e,m,a,b];sp[a,b]=im*candidate.Σᵍ[e,m,a,b]
            gc[a,b]=sn[a,b]+sp[a,b]
            dg[a,b]=gc[a,b]-gs[a,b]
            D[a,b]=holes[e]*sn[a,b]-f[e]*sp[a,b]
        end
        sandwich!(B,GR,dg,tmp);sandwich!(T,GR,D,tmp)
        sandwich!(As,GR,gs,tmp);sandwich!(nc,GR,sn,tmp);sandwich!(pc,GR,sp,tmp)
        block_squares=zeros(length(labels))
        block_numbers=zeros(ComplexF64,length(number_labels))
        for b=1:Nb,a=1:Nb
            n0=f[e]*A[a,b];p0=holes[e]*A[a,b]
            Es=A[a,b]-As[a,b]
            dn=nc[a,b]-n0;dp=pc[a,b]-p0
            bn=f[e]*B[a,b];bp=holes[e]*B[a,b]
            en=-f[e]*Es;ep=-holes[e]*Es
            rn=dn-(bn+T[a,b]+en);rp=dp-(bp-T[a,b]+ep)
            vals=(dn,dp,bn,bp,T[a,b],en,ep,rn,rp,Es,n0,p0,A[a,b])
            for j in eachindex(vals)
                block_squares[j]+=abs2(vals[j])
            end
            if a==b
                block_numbers .+= (n0,p0,nc[a,b],pc[a,b],bn,T[a,b],en,
                    f[e]*As[a,b],nc[a,b]-n0-bn-T[a,b]-en)
            end
        end
        weight=spin*g.wᴱ[e]*g.wᵏ[m]
        trusted=cfg["trusted_min_eV"] <= g.ε[e]*problem.scales.E₀_eV <= cfg["trusted_max_eV"]
        for domain in ("full","trusted")
            domain=="trusted" && !trusted && continue
            for (j,label) in enumerate(labels)
                squares[domain][label][1]+=block_squares[j]
                squares[domain][label][2]+=weight*block_squares[j]
            end
            for (j,label) in enumerate(number_labels)
                z=block_numbers[j]
                number[domain][label] .+= (weight*real(z),weight*abs(z),weight*imag(z))
            end
        end
    end
    result=Dict{String,Any}("definition"=>"candidate raw products on fixed archived GR; no normalization, no inversion",
        "mu_eV"=>CURRENT_MU[]*problem.scales.E₀_eV,"target_number"=>target)
    for domain in ("full","trusted")
        norms=Dict{String,Any}()
        for label in labels
            x=squares[domain][label]
            # All primary absolute norms are reported; input scales remain separate.
            norms[label]=Dict("unweighted_absolute"=>sqrt(x[1]),"weighted_absolute"=>sqrt(x[2]))
        end
        number_out=Dict{String,Any}()
        for label in number_labels
            x=number[domain][label]
            number_out[label]=Dict("signed"=>x[1],"absolute_local_integral"=>x[2],"imaginary"=>x[3],
                "signed_over_target"=>x[1]/target)
        end
        occscale=squares[domain]["occupied_input"]
        holescale=squares[domain]["hole_input"]
        occresponses=Dict(label=>Dict("unweighted_relative"=>sqrt(ratio(squares[domain][label][1],occscale[1])),
            "weighted_relative"=>sqrt(ratio(squares[domain][label][2],occscale[2]))) for label in
            ("occupied_raw_error","occupied_broadening_part","KMS_response","occupied_spectral_identity_correction","occupied_recombination_error"))
        holeresponses=Dict(label=>Dict("unweighted_relative"=>sqrt(ratio(squares[domain][label][1],holescale[1])),
            "weighted_relative"=>sqrt(ratio(squares[domain][label][2],holescale[2]))) for label in
            ("hole_raw_error","hole_broadening_part","KMS_response","hole_spectral_identity_correction","hole_recombination_error"))
        nr=number[domain]["candidate_raw_number"][1]
        ns=number[domain]["current_raw_number"][1]
        nfd=number[domain]["FD_number"][1]
        result[domain]=Dict("norms"=>norms,"numbers"=>number_out,"occupied_responses"=>occresponses,
            "hole_responses"=>holeresponses,"lambda_candidate_raw"=>nr/target,
            "lambda_current_FD_repaired_raw"=>ns/target,"FD_input_number_over_target"=>nfd/target,
            "r_K_core_lesser_reconstructed"=>sqrt(squares[domain]["occupied_raw_error"][1])/(sqrt(occscale[1])+1e-14),
            "candidate_number_minus_FD_input_over_target"=>(nr-nfd)/target)
    end
    result
end

function main()
    source=path_from_config("source_state")
    source_input=path_from_config("source_input")
    source_identity=path_from_config("source_identity")
    copied_script=joinpath(output,"one_fresh_fd_map.jl")
    core=dirname(dirname(pathof(QCLNEGF)))
    head=readchomp(`git -C $core rev-parse HEAD`)
    tree=readchomp(Cmd(["git","-C",core,"rev-parse","HEAD^{tree}"]))
    head==cfg["expected_core_head"] && tree==cfg["expected_core_tree"] || error("frozen Core differs")
    isempty(readchomp(`git -C $core status --porcelain`)) || error("Core checkout has local changes")
    input=TOML.parsefile(source_input)
    input["voltage_mV"]==0.0 && input["temperature_K"]==70.0 && input["alloy"] || error("unexpected physical input")
    input["embedding"]=="finite_chain" && input["embedding_periods"]==4 || error("embedding differs")
    (input["NE"],input["Nk"],input["Nb"])==(30721,11,5) || error("grid differs")
    state_sha=digest(source)
    state_sha==cfg["expected_source_sha256"] || error("source state SHA differs")
    cp(source_input,joinpath(output,"source-input.toml"))
    cp(source_identity,joinpath(output,"source-identity.toml"))
    baseline_channels=path_from_config("source_baseline_channels")
    isfile(baseline_channels) && cp(baseline_channels,joinpath(output,"source-baseline-channels.csv"))
    identity=Dict{String,Any}("id"=>cfg["id"],"start_utc"=>stamp(),"source_state"=>source,
        "source_state_sha256"=>state_sha,"source_input_sha256"=>digest(source_input),
        "source_identity_sha256"=>digest(source_identity),"core_head"=>head,"core_tree"=>tree,
        "script_sha256"=>digest(copied_script),"input_sha256"=>digest(joinpath(output,"input.toml")),
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "scope"=>cfg["scope"],"acceptance_decision"=>"not_requested_not_evaluated",
        "note_limits"=>"requested limits are metadata; external launcher must enforce them")
    write_toml("identity.toml",identity)
    progress("source identity verified")
    physical=reference_parameters(V_period=input["voltage_mV"]*u"mV",Tᴸ=input["temperature_K"]*u"K",Tᴸᴼ=input["temperature_K"]*u"K",
        ΔV_alloy=.60u"eV",Ω₀=4.5168e-29u"m^3")
    numerical=NumericalParameters(N_z=input["Nz"],N_b=input["Nb"],P_basis=input["basis_periods"],
        E_min=input["energy_min_eV"]*u"eV",E_max=input["energy_max_eV"]*u"eV",N_E=input["NE"],
        M_E=input["trust_margin_meV"]*u"meV",k_max=input["kmax_per_nm"]*u"nm^-1",N_k=input["Nk"],
        N_φ=input["Nphi"],qz_max=input["qzmax_per_nm"]*u"nm^-1",N_qz=input["Nqz"],η_seed=input["seed_meV"]*u"meV")
    algorithms=AlgorithmOptions(kernel_build=:direct,embedding=:finite_chain,embedding_periods=4,
        mixing=:linear,energy_shift=:sparse_plan)
    production=ProductionOptions(algorithms=algorithms,worker_count=1,hilbert_columns=cfg["hilbert_columns"],
        memory_budget_bytes=cfg["numerical_memory_budget_bytes"],phase_timing=false,
        checkpoint_every_scba=0,checkpoint_every_outer=0)
    scales=ScaleSystem()
    built=build_configured_scattering_problem(physical=physical,numerical=numerical,scales=scales,
        scattering=ScatteringOptions(alloy=true),algorithms=algorithms,kernel_options=ProductionKernelOptions())
    problem=built.problem
    estimate=BN.estimate_production_memory(problem;options=production)
    identity["native_solver_estimate_peak_bytes"]=estimate.peak_bytes
    identity["numerical_memory_budget_bytes"]=production.memory_budget_bytes
    identity["worker_count"]=production.worker_count
    identity["alloy_deltaV_eV"]=0.60
    identity["alloy_Omega0_m3"]=4.5168e-29
    identity["hilbert_columns"]=production.hilbert_columns
    Set(problem.kernels.enabled)==Set((:LO,:acoustic,:impurity,:IFR,:alloy)) || error("rebuilt all-five channels differ")
    write_toml("identity.toml",identity)
    # Keep the native solver estimate and its guard; do not shrink/bypass it for this diagnostic.
    cache=build_production_cache(problem;options=production)
    h=project_hamiltonians(problem,zeros(numerical.N_z))
    plan=BN.QCLNumerics._cavity_production_plan(problem,h,production)
    progress("same all-five-channel problem and production cache built")
    HDF5.h5open(source,"r") do file
        schema=read(file["schema"])
        schema=="campaign133-inner-raw-v2" || error("unexpected raw source schema")
        read(file["energy"])==problem.grids.ε || error("energy axis differs")
        read(file["k"])==problem.grids.κ || error("k axis differs")
        read(file["weight_energy"])==problem.grids.wᴱ || error("energy weights differ")
        read(file["weight_k"])==problem.grids.wᵏ || error("k weights differ")
        all(iszero,read(file["hartree"])) || error("source Hartree is not fixed zero")
        expected=Set(string.([problem.kernels.enabled;:embedding;:embedding_plus;:embedding_minus]))
        Set(keys(file["selfenergy"]))==expected || error("source mechanisms differ")
        GR=read_complex(file,"green/GR");A=read_complex(file,"green/A")
        shape=(numerical.N_E,numerical.N_k,numerical.N_b,numerical.N_b)
        size(GR)==size(A)==shape || error("Green shape differs")
        all(isfinite,GR) && all(isfinite,A) || error("nonfinite source Green")
        identity["source_A_sha256"]=bytes2hex(sha256(reinterpret(UInt8,vec(A))))
        identity["source_GR_sha256"]=bytes2hex(sha256(reinterpret(UInt8,vec(GR))))
        identity["complex_array_hash_format"]="native little-endian ComplexF64 interleaved real/imag; Julia column-major E,k,a,b"
        energy=problem.grids.ε;g=problem.grids;gs=physical.g_s
        kBT=ustrip(u"eV",BN.CODATA.kᴮₑᵥ*physical.Tᴸ)/scales.E₀_eV
        target=BN._scaled_physics(physical,scales).Nᴰ²ᴰ
        weights=[gs/(2pi)*g.wᴱ[e]*sum(g.wᵏ[m]*sum(real(A[e,m,a,a]) for a=1:numerical.N_b) for m=1:numerical.N_k) for e=1:numerical.N_E]
        all(isfinite,weights) && minimum(weights)>=0 || error("nonfinite/negative source spectral trace weights")
        number(mu)=sum(weights[e]*first(fermi_pair((energy[e]-mu)/kBT)) for e in eachindex(energy))
        lo=first(energy)-50kBT;hi=last(energy)+50kBT
        number(lo)<=target<=number(hi) || error("FD spectral target not bracketed")
        for step=1:100
            mid=(lo+hi)/2
            number(mid)<target ? (lo=mid) : (hi=mid)
        end
        mu=(lo+hi)/2;CURRENT_MU[]=mu
        f=zeros(length(energy));holes=similar(f)
        for e in eachindex(energy)
            f[e],holes[e]=fermi_pair((energy[e]-mu)/kBT)
        end
        lesser=similar(A);greater=similar(A)
        for e in eachindex(energy),m=1:numerical.N_k,b=1:numerical.N_b,a=1:numerical.N_b
            lesser[e,m,a,b]=im*f[e]*A[e,m,a,b]
            greater[e,m,a,b]=-im*holes[e]*A[e,m,a,b]
        end
        green=GreenState(GR,lesser,greater,A,read(file["green/condition_number"]),read(file["green/dyson_scale"]))
        identity["mu_eV"]=mu*scales.E₀_eV;identity["kBT_eV"]=kBT*scales.E₀_eV
        identity["capacity_solver"]=sum(weights);identity["target_number_solver"]=target
        identity["FD_number_solver"]=number(mu)
        identity["E0_eV"]=scales.E₀_eV;identity["L0_m"]=scales.L₀_m
        identity["input_selection"]="archived GR/A unchanged; one fitted mu; only Green Keldysh input is FD; original SigmaR retained as audit reference"
        write_toml("identity.toml",identity)
        progress("FD pair constructed on unchanged archived A")
        fresh_scattering=BN._scattering_candidate_production(problem,green,cache,production)
        progress("one unprojected scattering map completed")
        fresh_embedding,fresh_plus,fresh_minus=BN._embedding_self_energy_selected(problem,green,cache,production;
            scattering=fresh_scattering,h=h,cavity_plan=plan)
        candidate=BN._total_family_production(fresh_scattering,fresh_embedding,problem.kernels.enabled,production)
        progress("fresh finite-chain embedding from fresh local scattering completed")
        families=Tuple{Symbol,SelfEnergyFamily}[(name,fresh_scattering[name]) for name in problem.kernels.enabled]
        append!(families,[(:embedding,fresh_embedding),(:embedding_plus,fresh_plus),(:embedding_minus,fresh_minus)])
        source_totalR=zeros(ComplexF64,shape)
        all_csv=Vector{Vector{Any}}();channel_summaries=Dict{String,Any}()
        for (name,family) in families
            sourceR=read_complex(file,"selfenergy/"*string(name)*"/retarded")
            size(sourceR)==shape && all(isfinite,sourceR) || error("source retarded family differs")
            if name in problem.kernels.enabled || name==:embedding
                source_totalR .+= sourceR
            end
            summary,rows=family_audit(name,family,sourceR,green,f,holes,problem)
            write_toml("channel-$(name).toml",summary)
            channel_summaries[string(name)]=summary
            append!(all_csv,rows)
            sourceR=nothing
        end
        total_summary,rows=family_audit(:total,candidate,source_totalR,green,f,holes,problem)
        write_toml("channel-total.toml",total_summary)
        channel_summaries["total"]=total_summary;append!(all_csv,rows)
        open(joinpath(output,"channels.csv"),"w") do io
            println(io,"channel,component,residual_absolute_unweighted,candidate_scale_unweighted,reference_norm_unweighted,core_marker_relative,weighted_absolute,weighted_candidate_scale,weighted_relative")
            for row in all_csv
                println(io,join(row,","))
            end
        end
        # Release per-channel full arrays; the summed candidate owns independent storage.
        families=nothing;fresh_scattering=nothing;fresh_embedding=nothing;fresh_plus=nothing;fresh_minus=nothing
        GC.gc()
        dec=decomposition(problem,green,source_totalR,candidate,f,holes,target)
        write_toml("decomposition.toml",dec)
        progress("raw Keldysh/number decomposition completed")
        if cfg["compute_core_residual_suites"]
            currentL=similar(source_totalR);currentG=similar(source_totalR)
            for e in eachindex(energy),m=1:numerical.N_k,b=1:numerical.N_b,a=1:numerical.N_b
                gamma=im*(source_totalR[e,m,a,b]-conj(source_totalR[e,m,b,a]))
                currentL[e,m,a,b]=im*f[e]*gamma
                currentG[e,m,a,b]=-im*holes[e]*gamma
            end
            current=SelfEnergyFamily(source_totalR,currentL,currentG)
            suite_current=BN.production_residual_suite(energy,h,green,current,candidate;worker_count=1)
            suite_candidate=BN.production_residual_suite(energy,h,green,candidate,candidate;worker_count=1)
            write_toml("core-raw-audit.toml",Dict("current_FD_repaired"=>Dict(string(k)=>v for (k,v) in pairs(suite_current)),
                "fresh_candidate_on_unchanged_GR"=>Dict(string(k)=>v for (k,v) in pairs(suite_candidate)),
                "note"=>"Fresh rD/rA measure candidate mismatch on old GR; fresh PSD checks raw candidate products. No convergence/acceptance decision."))
        end
        write_toml("completion.toml",Dict("status"=>"one_map_diagnostic_completed_without_acceptance_decision",
            "elapsed_seconds"=>time()-started,"finish_utc"=>stamp(),"map_evaluations"=>1,
            "solver_calls"=>0,"central_Dyson_recomputations"=>0,"mixing_steps"=>0,"Poisson_calls"=>0,
            "cavity_inversions"=>"native finite-chain map retained"))
        progress("all diagnostic outputs written")
    end
end

try
    main()
catch err
    write_toml("failure.toml",Dict("status"=>"diagnostic_failed","error"=>sprint(showerror,err),
        "elapsed_seconds"=>time()-started,"time_utc"=>stamp()))
    rethrow()
end
