using QCLNEGF, Unitful, LinearAlgebra, HDF5, TOML, SHA, Dates
const BN = QCLNEGF
BLAS.set_num_threads(1)
Threads.nthreads() == 1 || error("one Julia thread required")
const started = time()
const out = @__DIR__
const cfgpath = joinpath(out, "inputs.toml")
const cfg = TOML.parsefile(cfgpath)
const base = dirname(out)
const source = normpath(joinpath(base, cfg["source_state"]))
const source_input = normpath(joinpath(base, cfg["source_input"]))
const input = TOML.parsefile(source_input)
write_toml(name, value) = open(io -> TOML.print(io, value), joinpath(out,name), "w")
ratio(n,d) = iszero(d) ? (iszero(n) ? 0.0 : Inf) : n/d
fermi_pair(x) = x >= 0 ? (let z=exp(-x); (z/(1+z),1/(1+z)) end) :
                             (let z=exp(x); (1/(1+z),z/(1+z)) end)

function rows_for_shift(energy, delta; integer_steps=0)
    count=length(energy)
    rows=Vector{NTuple{4,Float64}}(undef,count)
    for i in eachindex(energy)
        if integer_steps != 0
            j=i+integer_steps
            rows[i] = 1 <= j <= count ? (Float64(j),Float64(j),1.0,0.0) : (1.0,1.0,0.0,0.0)
            continue
        end
        x=energy[i]+delta
        if x < first(energy) || x > last(energy)
            rows[i]=(1.0,1.0,0.0,0.0)
        elseif x == last(energy)
            rows[i]=(Float64(count),Float64(count),1.0,0.0)
        else
            j=clamp(searchsortedlast(energy,x),1,count-1)
            theta=(x-energy[j])/(energy[j+1]-energy[j])
            rows[i]=(Float64(j),Float64(j+1),1-theta,theta)
        end
    end
    rows
end

function resample(X, old, new)
    Y=zeros(ComplexF64,length(new),size(X,2),size(X,3),size(X,4))
    for i in eachindex(new)
        j=clamp(searchsortedlast(old,new[i]),1,length(old)-1)
        theta=(new[i]-old[j])/(old[j+1]-old[j])
        @views Y[i,:,:,:] .= (1-theta).*X[j,:,:,:] .+ theta.*X[j+1,:,:,:]
    end
    Y
end

# Literal compound-index assembly, independent of production packing helpers.
function independent_contraction(K, A, wk, qK)
    NE,Nk,Nb,_=size(A)
    P=zeros(ComplexF64,Nk*Nb^2,Nk*Nb^2)
    X=zeros(ComplexF64,Nk*Nb^2,NE)
    for m=1:Nk, mp=1:Nk, a=1:Nb, b=1:Nb, c=1:Nb, d=1:Nb
        o=(m-1)*Nb^2+(b-1)*Nb+a
        j=(mp-1)*Nb^2+(d-1)*Nb+c
        P[o,j]=qK*wk[mp]*K[m,mp,a,c,d,b]
    end
    for e=1:NE, m=1:Nk, a=1:Nb, b=1:Nb
        X[(m-1)*Nb^2+(b-1)*Nb+a,e]=A[e,m,a,b]
    end
    Y=P*X
    B=similar(A)
    for e=1:NE, m=1:Nk, a=1:Nb, b=1:Nb
        B[e,m,a,b]=Y[(m-1)*Nb^2+(b-1)*Nb+a,e]
    end
    B
end

function probe(name, energy, A, B, wk, mu, kBT, omega, bose, gs, E0, L0;
               integer_steps=0, production_pair=nothing)
    NE,Nk,Nb,_=size(A)
    h=(last(energy)-first(energy))/(NE-1)
    wE=fill(h,NE); wE[[1,end]]./=2
    f=zeros(NE); holes=zeros(NE)
    for e=1:NE
        f[e],holes[e]=fermi_pair((energy[e]-mu)/kBT)
    end
    plus=rows_for_shift(energy,omega;integer_steps)
    minus=rows_for_shift(energy,-omega;integer_steps=-integer_steps)
    # Accumulators: absolute squared defect, self-energy norm squared,
    # source-prediction squared difference, C0, CE, energy moment, activity,
    # absolute local collision, C0 imaginary, CE imaginary, number, capacity.
    values=Dict(domain=>zeros(12) for domain in ("full","common_trusted"))
    prod_diff2=prod_scale2=0.0
    dmatrix=zeros(ComplexF64,Nb,Nb)
    sn=similar(dmatrix); sp=similar(dmatrix); pred=similar(dmatrix)
    for e=1:NE, m=1:Nk
        fill!(sn,0); fill!(sp,0); fill!(pred,0)
        for (rows, sign) in ((plus,1),(minus,-1))
            j1,j2,t1,t2=rows[e]
            for (jf,w) in ((j1,t1),(j2,t2))
                iszero(w) && continue
                j=Int(jf)
                cn=sign==1 ? (bose+1)*f[j] : bose*f[j]
                cp=sign==1 ? bose*holes[j] : (bose+1)*holes[j]
                # Independent thermodynamic source coefficient oracle.
                coefficient=sign==1 ? (bose+1)*holes[e]*f[j]-bose*f[e]*holes[j] :
                    bose*holes[e]*f[j]-(bose+1)*f[e]*holes[j]
                for b=1:Nb,a=1:Nb
                    val=w*B[j,m,a,b]
                    sn[a,b]+=cn*val
                    sp[a,b]+=cp*val
                    pred[a,b]+=coefficient*val
                end
            end
        end
        inc=0.0im; outgoing=0.0im; err2=scale2=prediction_error2=0.0
        for b=1:Nb,a=1:Nb
            D=holes[e]*sn[a,b]-f[e]*sp[a,b]
            dmatrix[a,b]=D
            err2+=abs2(D)
            scale2+=abs2(sn[a,b])+abs2(sp[a,b])
            prediction_error2+=abs2(D-pred[a,b])
            inc+=sn[a,b]*holes[e]*A[e,m,b,a]
            outgoing+=sp[a,b]*f[e]*A[e,m,b,a]
            if production_pair !== nothing
                sl,sg=production_pair
                pn=-im*sl[e,m,a,b]; pp=im*sg[e,m,a,b]
                prod_diff2+=wE[e]*wk[m]*(abs2(sn[a,b]-pn)+abs2(sp[a,b]-pp))
                prod_scale2+=wE[e]*wk[m]*(abs2(sn[a,b])+abs2(sp[a,b]))
            end
        end
        collision=inc-outgoing
        weight=gs/(2pi)*wE[e]*wk[m]
        traceA=sum(real(A[e,m,a,a]) for a=1:Nb)
        for domain in ("full","common_trusted")
            domain=="common_trusted" && !(cfg["trusted_min_eV"]/E0 <= energy[e] <= cfg["trusted_max_eV"]/E0) && continue
            v=values[domain]
            v .+= weight .* [err2,scale2,prediction_error2,real(collision),
                (energy[e]-mu)*real(collision),energy[e]*real(collision),
                abs(inc)+abs(outgoing),abs(collision),imag(collision),
                (energy[e]-mu)*imag(collision),f[e]*traceA,traceA]
        end
    end
    hbar_eVs=6.582119569e-16
    electron_charge=1.602176634e-19
    rate_scale=E0/(hbar_eVs*L0^2)
    power_scale=E0^2*electron_charge/(hbar_eVs*L0^2)
    result=Dict{String,Any}("name"=>name,"NE"=>NE,"energy_min_eV"=>first(energy)*E0,
        "energy_max_eV"=>last(energy)*E0,"dE_eV"=>h*E0,"LO_energy_eV"=>omega*E0,
        "LO_over_dE"=>omega/h,"integer_index_shift"=>integer_steps,"mu_eV"=>mu*E0,
        "A_representation"=>(name=="original_fractional" ? "actual archived A unchanged" :
           "linear interpolation of actual archived A, no renormalization"),
        "B_representation"=>(name=="original_fractional" ? "independent KLO[actual A]" :
           "same linear interpolation of KLO[actual A]; KLO is energy-independent"))
    for (domain,v) in values
        result[domain]=Dict("KMS_absolute_weighted_solver"=>sqrt(v[1]),
          "KMS_scale_weighted_solver"=>sqrt(v[2]),"KMS_relative"=>sqrt(ratio(v[1],v[2])),
          "source_prediction_difference_absolute"=>sqrt(v[3]),
          "source_prediction_difference_relative"=>sqrt(ratio(v[3],v[2])),
          "C0_signed_solver"=>v[4],"CE_mu_signed_solver"=>v[5],"CE_zero_signed_solver"=>v[6],
          "activity_solver"=>v[7],"local_collision_absolute_solver"=>v[8],
          "C0_imaginary_solver"=>v[9],"CE_mu_imaginary_solver"=>v[10],
          "C0_over_activity"=>ratio(v[4],v[7]),"CE_mu_over_LO_activity"=>ratio(v[5],omega*v[7]),
          "number_solver"=>v[11],"capacity_solver"=>v[12],
          "C0_per_m2_s"=>v[4]*rate_scale,"CE_mu_W_per_m2"=>v[5]*power_scale,
          "activity_per_m2_s"=>v[7]*rate_scale)
    end
    if production_pair !== nothing
        result["production_comparison"]=Dict("absolute_weighted_solver"=>sqrt(gs/(2pi)*prod_diff2),
            "relative"=>sqrt(ratio(prod_diff2,prod_scale2)))
    end
    write_toml(name*".toml",result)
    println(name," KMS=",result["full"]["KMS_relative"]," CE=",result["full"]["CE_mu_over_LO_activity"]," time=",time()-started)
    flush(stdout)
    result
end

function main()
    cfg["core_head"]==readchomp(`git -C $(dirname(dirname(pathof(QCLNEGF)))) rev-parse HEAD`) || error("Core HEAD changed")
    physical=reference_parameters(V_period=0u"mV",Tᴸ=input["temperature_K"]*u"K",Tᴸᴼ=input["temperature_K"]*u"K")
    numerical=NumericalParameters(N_z=input["Nz"],N_b=input["Nb"],P_basis=input["basis_periods"],
        E_min=input["energy_min_eV"]*u"eV",E_max=input["energy_max_eV"]*u"eV",N_E=input["NE"],
        M_E=input["trust_margin_meV"]*u"meV",k_max=input["kmax_per_nm"]*u"nm^-1",N_k=input["Nk"],
        N_φ=input["Nphi"],qz_max=input["qzmax_per_nm"]*u"nm^-1",N_qz=input["Nqz"],η_seed=input["seed_meV"]*u"meV")
    scales=ScaleSystem()
    algorithms=AlgorithmOptions(kernel_build=:direct,energy_shift=:sparse_plan,embedding=:finite_chain,embedding_periods=4)
    options=ProductionOptions(algorithms=algorithms,worker_count=1,memory_budget_bytes=2*1024^3)
    built=build_configured_scattering_problem(physical=physical,numerical=numerical,scales=scales,
        scattering=ScatteringOptions(LO=true,acoustic=false,impurity=false,IFR=false,alloy=false),
        algorithms=algorithms,kernel_options=ProductionKernelOptions())
    problem=built.problem
    println("built LO-only kernel, seconds=",time()-started);flush(stdout)
    A,energy,wE,wk=HDF5.h5open(source,"r") do file
        ar=read(file["green/A_real"]);ai=read(file["green/A_imag"])
        (complex.(ar,ai),read(file["energy"]),read(file["weight_energy"]),read(file["weight_k"]))
    end
    size(A)==(input["NE"],input["Nk"],input["Nb"],input["Nb"]) || error("A shape differs")
    energy==problem.grids.ε && wk==problem.grids.wᵏ || error("rebuilt grids differ from archive")
    all(isfinite,A) || error("nonfinite actual A")
    Ahash=bytes2hex(sha256(reinterpret(UInt8,vec(A))))
    E0=scales.E₀_eV;L0=scales.L₀_m;gs=physical.g_s
    kBT=ustrip(u"eV",BN.CODATA.kᴮₑᵥ*physical.Tᴸ)/E0
    omega=ustrip(u"eV",physical.ħωᴸᴼ)/E0
    bose=inv(expm1(omega/kBT))
    weights=[gs/(2pi)*wE[e]*sum(wk[m]*sum(real(A[e,m,a,a]) for a=1:size(A,3)) for m=1:length(wk)) for e=1:length(energy)]
    minimum(weights)>=0 || error("negative spectral trace weight")
    target=BN._scaled_physics(physical,scales).Nᴰ²ᴰ
    charge(mu)=sum(weights[e]*first(fermi_pair((energy[e]-mu)/kBT)) for e in eachindex(energy))
    lo=first(energy)-50kBT;hi=last(energy)+50kBT
    charge(lo)<=target<=charge(hi) || error("spectral target unbracketed")
    for i=1:100
        mid=(lo+hi)/2
        charge(mid)<target ? (lo=mid) : (hi=mid)
    end
    mu=(lo+hi)/2
    spectral_min=Inf; hermitian_max=0.0
    for e in eachindex(energy),m in eachindex(wk)
        block=Matrix(view(A,e,m,:,:))
        hermitian_max=max(hermitian_max,norm(block-block'))
        spectral_min=min(spectral_min,eigmin(Hermitian((block+block')/2)))
    end
    boundary_lo=findall(e->energy[e]<first(energy)+omega || energy[e]>last(energy)-omega,eachindex(energy))
    identity=Dict("source_state"=>source,"source_state_sha256"=>bytes2hex(open(sha256,source)),
        "source_A_sha256"=>Ahash,"source_A_hash_format"=>"ComplexF64 native little-endian interleaved real/imag, Julia column-major E,k,a,b",
        "source_input_sha256"=>bytes2hex(sha256(read(source_input))),"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "input_sha256"=>bytes2hex(sha256(read(cfgpath))),"core_head"=>cfg["core_head"],
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "E0_eV"=>E0,"L0_m"=>L0,"kBT_eV"=>kBT*E0,"LO_energy_eV"=>omega*E0,"NLO"=>bose,
        "mu_eV"=>mu*E0,"target_number_solver"=>target,"FD_number_solver"=>charge(mu),"capacity_solver"=>sum(weights),
        "spectral_min_eigenvalue_solver"=>spectral_min,"spectral_max_antiHermitian_norm_solver"=>hermitian_max,
        "LO_boundary_spectral_fraction"=>sum(weights[boundary_lo])/sum(weights),
        "LO_boundary_FD_occupied_fraction"=>sum(weights[e]*first(fermi_pair((energy[e]-mu)/kBT)) for e in boundary_lo)/target,
        "comparison_scope"=>"frozen B0 spectrum; not SCBA, not solution-grid convergence, not acceptance of previous Core trees")
    write_toml("identity.toml",identity)
    B=independent_contraction(problem.kernels.K[:LO],A,wk,problem.kernels.qᴷ[:LO])
    println("one K[A] complete, seconds=",time()-started);flush(stdout)
    # One production LO call on the original FD pair, solely for algebraic comparison.
    lesser=similar(A);greater=similar(A)
    for e in eachindex(energy),m in eachindex(wk),a in axes(A,3),b in axes(A,4)
        f,h=fermi_pair((energy[e]-mu)/kBT)
        lesser[e,m,a,b]=im*f*A[e,m,a,b];greater[e,m,a,b]=-im*h*A[e,m,a,b]
    end
    green=GreenState(zeros(ComplexF64,size(A)),lesser,greater,A,ones(length(energy),length(wk)),ones(length(energy),length(wk)))
    cache=build_production_cache(problem;options)
    production_pair=BN._production_lo_contraction(cache.kernels[:LO],green,problem,cache,problem.kernels.qᴷ[:LO],options)
    probe("original_fractional",energy,A,B,wk,mu,kBT,omega,bose,gs,E0,L0;production_pair)
    production_pair=nothing;green=nothing;lesser=nothing;greater=nothing;cache=nothing;GC.gc()
    half=collect(range(first(energy),last(energy);length=2length(energy)-1))
    Af=resample(A,energy,half);Bf=resample(B,energy,half)
    probe("half_step_fractional",half,Af,Bf,wk,mu,kBT,omega,bose,gs,E0,L0)
    Af=nothing;Bf=nothing;GC.gc()
    hc=omega/cfg["commensurate_steps"]
    count=floor(Int,(last(energy)-first(energy))/hc)+1
    comm=[first(energy)+(i-1)*hc for i=1:count]
    Ac=resample(A,energy,comm);Bc=resample(B,energy,comm)
    probe("commensurate_128",comm,Ac,Bc,wk,mu,kBT,omega,bose,gs,E0,L0;integer_steps=cfg["commensurate_steps"])
    write_toml("completion.toml",Dict("status"=>"diagnostic_completed_without_acceptance_decision",
        "elapsed_seconds"=>time()-started,"finish_utc"=>string(Dates.now(Dates.UTC))*"Z"))
end
main()
