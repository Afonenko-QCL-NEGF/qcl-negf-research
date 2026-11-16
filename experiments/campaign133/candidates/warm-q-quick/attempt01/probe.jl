using QCLNEGF, HDF5, LinearAlgebra, TOML, SHA, Dates
include(abspath(".build/campaign133/equilibrium-gamma-warm/run_gamma.jl"))
const E = EquilibriumGammaPilot
const N = QCLNEGF.QCLNumerics

function main()
    input, output = abspath.(ARGS)
    ispath(output) && error("Refusing to overwrite output")
    cfg=TOML.parsefile(input); cfg["_input_directory"]=dirname(input)
    core=dirname(dirname(pathof(QCLNEGF)))
    readchomp(`git -C $core rev-parse HEAD`)==cfg["expected_core_head"] || error("Core identity")
    isempty(readchomp(`git -C $core status --porcelain --untracked-files=no`)) || error("Dirty Core")
    BLAS.set_num_threads(1); mkpath(output)
    cp(@__FILE__,joinpath(output,"probe.jl")); cp(input,joinpath(output,"input.toml"))
    started=time(); report=Dict{String,Any}("scope"=>"one native finite4 R/K spectral probe; no mu fit/scattering/Phi/SCBA",
        "started_utc"=>string(Dates.now(Dates.UTC)),"core_head"=>cfg["expected_core_head"],
        "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"input_sha256"=>bytes2hex(sha256(read(input))))
    try
        c=E.build_context(cfg); g=c.problem.grids
        gamma=E.load_warm_seed(c,cfg,output)
        localR,roundoff=E.local_retarded(c,gamma)
        # Retarded recurrence is independent of kinetic components. The filled
        # f=1 probe carries physical spectral sources, not fitted occupations.
        probeL=im.*gamma; zeroG=zeros(ComplexF64,size(gamma))
        emb,plus,minus=N._cavity_embedding_production(c.problem,SelfEnergyFamily(localR,probeL,zeroG),c.h,c.options;plan=c.cavity)
        totalR=localR+emb.Σᴿ
        GR,condition,scale=N._retarded_green_production(g.ε,c.h,totalR;η=0.0,options=c.options)
        q=N._seed_number_weights_production(GR,g,c.problem.physical.g_s;options=c.options)
        A=E.spectral_function(GR); report["subtraction_A_cone"]=Dict(string(k)=>v for (k,v) in pairs(E.cone(A;require_nonzero=false)))
        bad=findall(<(0.0),q); worst=argmin(q)
        chosen=sort(unique([worst; isempty(bad) ? Int[] : [first(bad),last(bad)]; argmin(abs.(g.ε.-0.085)); length(q)]))
        ne,nk,nb,_=size(GR); selection=Dict(e=>j for (j,e) in enumerate(chosen))
        qeff=zeros(ne); qphysical=zeros(ne)
        selected=Dict(name=>zeros(ComplexF64,length(chosen),nk,nb,nb) for name in ("Gamma_eff","Gamma_physical_total","A_eff","A_physical"))
        G=zeros(ComplexF64,nb,nb); B=similar(G); H=similar(G); K=similar(G); tmp=similar(G); P=similar(G)
        max_h_defect=0.0; max_lambda_defect=0.0; max_product_imaginary_trace=0.0
        pref=c.problem.physical.g_s/(2pi)
        for e in 1:ne, m in 1:nk
            for b in 1:nb, a in 1:nb
                G[a,b]=GR[e,m,a,b]
                H[a,b]=c.h[m,a,b]+totalR[e,m,a,b]
                K[a,b]=gamma[e,m,a,b]-im*emb.Σˡ[e,m,a,b]
                B[a,b]=im*(H[a,b]-conj(c.h[m,b,a]+totalR[e,m,b,a]))
            end
            mul!(tmp,G,B); mul!(P,tmp,adjoint(G))
            qeff[e]+=pref*g.wᴱ[e]*g.wᵏ[m]*real(tr(P))
            if haskey(selection,e)
                j=selection[e]; selected["Gamma_eff"][j,m,:,:].=B; selected["A_eff"][j,m,:,:].=P
            end
            mul!(tmp,G,K); mul!(P,tmp,adjoint(G))
            qphysical[e]+=pref*g.wᴱ[e]*g.wᵏ[m]*real(tr(P))
            max_product_imaginary_trace=max(max_product_imaginary_trace,abs(imag(tr(P))))
            if haskey(selection,e)
                j=selection[e]; selected["Gamma_physical_total"][j,m,:,:].=K; selected["A_physical"][j,m,:,:].=P
                hb=Matrix(view(c.h,m,:,:)); lambda=Matrix(view(localR,e,m,:,:))+0.5im*Matrix(view(gamma,e,m,:,:))
                max_h_defect=max(max_h_defect,norm(hb-hb'))
                max_lambda_defect=max(max_lambda_defect,norm(lambda-lambda'))
            end
        end
        for (name,v) in (("raw",q),("same_G_effective_broadening",qeff),("same_G_physical_sources",qphysical))
            negative=findall(<(0.0),v)
            report[name]=Dict("minimum"=>minimum(v),"maximum"=>maximum(v),"negative_count"=>length(negative),
                "zero_count"=>count(iszero,v),"negative_absolute_mass"=>-sum(v[negative];init=0.0),
                "negative_mass_over_target"=>-sum(v[negative];init=0.0)/c.target,
                "capacity"=>sum(v),"all_finite"=>all(isfinite,v))
        end
        report["worst_energy_index"]=worst; report["worst_energy_eV"]=g.ε[worst]*c.problem.scales.E₀_eV
        report["selected_energy_indices"]=chosen; report["selected_h_antihermiticity"]=max_h_defect
        report["selected_lambda_antihermiticity"]=max_lambda_defect
        report["maximum_product_imaginary_trace"]=max_product_imaginary_trace
        report["raw_minus_effective_L1_over_target"]=sum(abs,q-qeff)/c.target
        report["effective_minus_physical_L1_over_target"]=sum(abs,qeff-qphysical)/c.target
        report["raw_minus_physical_L1_over_target"]=sum(abs,q-qphysical)/c.target
        report["PI_roundoff"]=roundoff; report["maximum_condition"]=maximum(condition)
        h5open(joinpath(output,"selected.h5"),"w") do f
            f["schema"]="campaign133-warm-q-selected-v1"; f["energy"]=g.ε; f["weight_energy"]=g.wᴱ; f["weight_k"]=g.wᵏ
            f["energy_scale_eV"]=c.problem.scales.E₀_eV; f["spin_degeneracy"]=c.problem.physical.g_s; f["target"]=c.target
            f["selected_energy_indices"]=chosen; f["q_raw"]=q; f["q_eff"]=qeff; f["q_physical"]=qphysical
            for (name,field) in (("Gamma",gamma),("localR",localR),("totalR",totalR),("GR",GR),("A_sub",A),("embeddingR",emb.Σᴿ))
                E.write_complex(f,name,field[chosen,:,:,:])
            end
            for (name,field) in selected; E.write_complex(f,name,field); end
            E.write_complex(f,"h",c.h); E.write_complex(f,"Tplus",c.problem.basis.T₊)
        end
        report["execution"]="complete"
        println("QProbeComplete negatives=",length(bad)," minimum=",minimum(q)," physical_negatives=",count(<(0.0),qphysical))
    catch err
        report["execution"]="stopped"; report["reason"]=sprint(showerror,err,catch_backtrace()); rethrow()
    finally
        report["elapsed_seconds"]=time()-started
        open(joinpath(output,"summary.toml"),"w") do io; TOML.print(io,report); end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("probe INPUT NEW_OUTPUT")
    main()
end
