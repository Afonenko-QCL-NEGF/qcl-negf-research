using QCLNEGF, HDF5, LinearAlgebra, TOML, SHA, Test
const HELPERS=joinpath(@__DIR__,"../equilibrium-gamma-affine4096/attempt01/gamma_functions.jl")
bytes2hex(open(sha256,HELPERS))=="cf923bfb7e5dc56d5aa296c6863e8615a92984f82e125433ea1002bac9c9ee63" || error("frozen helper differs")
include(HELPERS)
const E=EquilibriumGammaPilot
const N=QCLNEGF.QCLNumerics
const Field=E.Field

"Generic native-scale covariance arithmetic; no eigen/spectral clipping."
function covariance_A(GR::Field,K::Field)
    A=similar(GR);ne,nk,nb,_=size(GR)
    N._production_worker_ranges!(ne*nk,min(4,Threads.nthreads())) do _,first,last
        g=zeros(ComplexF64,nb,nb);k=similar(g);tmp=similar(g);out=similar(g)
        @inbounds for linear in first:last
            e=(linear-1)%ne+1;m=(linear-1)÷ne+1
            for b in 1:nb,a in 1:nb;g[a,b]=GR[e,m,a,b];k[a,b]=K[e,m,a,b];end
            mul!(tmp,g,k);mul!(out,tmp,adjoint(g))
            for b in 1:nb,a in 1:nb;A[e,m,a,b]=out[a,b];end
        end
    end
    return A
end
function number_weights_A(A::Field,g,spin)
    q=zeros(Float64,size(A,1));pref=Float64(spin)/(2pi)
    @inbounds for e in axes(A,1)
        momentum=0.0
        for m in axes(A,2)
            trace=0.0
            for a in axes(A,3);trace+=real(A[e,m,a,a]);end
            momentum+=g.wᵏ[m]*trace
        end
        q[e]=pref*g.wᴱ[e]*momentum
    end
    return q
end
function qinfo(q,target)
    bad=findall(<(0.0),q)
    Dict("minimum"=>minimum(q),"negative_count"=>length(bad),"all_finite"=>all(isfinite,q),
        "negative_absolute_mass"=>-sum(q[bad];init=0.0),"negative_over_target"=>-sum(q[bad];init=0.0)/target,
        "capacity"=>sum(q),"absolute_mass"=>sum(abs,q),"minimum_index"=>argmin(q))
end
function spectrum(c,gamma)
    g=c.problem.grids
    Lraw,roundoff=N._production_pi_fft_hilbert_transform(gamma,g.ε,g.wᴱ,c.cache.product_hilbert_operator)
    LH=E.hermitian_copy(Lraw);correction=E.norms(Lraw,LH,g.wᴱ,g.wᵏ)
    correction.relative<=1e-10 && correction.weighted_relative<=1e-10 || error("existing PI symmetry target exceeded")
    S=copy(LH)
    @inbounds @simd for i in eachindex(S);S[i]-=im*gamma[i]/2;end
    probe=SelfEnergyFamily(S,im.*gamma,zeros(ComplexF64,size(gamma)))
    emb,plus,minus=N._cavity_embedding_production(c.problem,probe,c.h,c.options;plan=c.cavity)
    totalR=S+emb.Σᴿ;physical=gamma.-im.*emb.Σˡ
    GR,condition,scale=N._retarded_green_production(g.ε,c.h,totalR;η=0.0,options=c.options)
    A=covariance_A(GR,physical);Asub=spectral_function(GR)
    effective=similar(physical)
    @inbounds for b in axes(GR,4),a in axes(GR,3),m in axes(GR,2),e in axes(GR,1)
        effective[e,m,a,b]=im*(c.h[m,a,b]+totalR[e,m,a,b]-conj(c.h[m,b,a]+totalR[e,m,b,a]))
    end
    return (GR=GR,A=A,Asub=Asub,condition=condition,scale=scale,physical=physical,effective=effective,
        localR=S,totalR=totalR,Lambda_raw=Lraw,Lambda_H_preassembly=LH,roundoff=roundoff,correction=correction)
end
function native_compare(c,green,rawphi,channels,mu)
    p=c.problem;g=p.grids
    families,roundoff=N._scattering_candidate_production(p,green,c.cache,c.options;return_roundoff=true)
    records=Dict{String,Any}();sumR=zeros(ComplexF64,size(rawphi));maximum_error=0.0
    for id in p.kernels.enabled
        family=families[id];K=im.*(family.Σᵍ.-family.Σˡ)
        metric=E.norms(channels[id],K,g.wᴱ,g.wᵏ);kms=E.kms_metrics(family,K,g,mu,c.kBT)
        records[string(id)]=E.channel_check_record(metric,kms)
        maximum_error=max(maximum_error,metric.relative,kms["relative"]);sumR.+=family.Σᴿ
    end
    projected,_=E.local_retarded(c,rawphi)
    retarded=E.norms(sumR,projected,g.wᴱ,g.wᵏ)
    return Dict("per_channel"=>records,"retarded_sum_vs_HermPi"=>Dict(string(k)=>v for (k,v) in pairs(retarded)),
        "roundoff"=>roundoff,"maximum_kinetic_KMS_relative"=>maximum_error,
        "original_kinetic_KMS_target"=>1e-10,"kinetic_KMS_control_pass"=>(maximum_error<=1e-10),
        "scope"=>"same supplied stable-A Green; native kinetic map comparison, NOT unconstrained native fixed-point acceptance")
end
function tiny_checks()
    @testset "same-A positive-source primitive" begin
        gamma=reshape(ComplexF64[1e-27],1,1,1,1);G=reshape(ComplexF64[inv(1.7+0.5e-27im)],1,1,1,1)
        A=covariance_A(G,gamma)
        @test real(only(A))>0
        @test real(only(A))≈1e-27/(1.7^2+0.25e-54) rtol=1e-14
        q=number_weights_A(A,(wᴱ=[.2],wᵏ=[.3]),2)
        @test only(q)≈real(only(A))*0.2*0.3/pi rtol=1e-14
        H=ComplexF64[0 1+im;1-im 2];K=ComplexF64[.04 .01im;-.01im .03]
        gr=inv(3I-H+0.5im*K);gf=zeros(ComplexF64,1,1,2,2);kf=similar(gf)
        gf[1,1,:,:]=gr;kf[1,1,:,:]=K
        @test covariance_A(gf,kf)[1,1,:,:]≈im*(gr-gr') rtol=2e-13
        @test number_weights_A(covariance_A(gf,kf),(wᴱ=[1.0],wᵏ=[1.0]),2)[1]>0
        # Weak, complex noncommuting source: independent Big192 resolvent + loops.
        weak=1e-27 .* K;@test norm(H*weak-weak*H)>0
        oracle=setprecision(192) do
            hb=Complex{BigFloat}.(H);kb=Complex{BigFloat}.(weak)
            gb=inv(3I-hb+BigFloat(0.5)*im*kb);ab=zeros(Complex{BigFloat},2,2)
            for a in 1:2,b in 1:2,i in 1:2,j in 1:2;ab[a,b]+=gb[a,i]*kb[i,j]*conj(gb[b,j]);end
            ComplexF64.(ab)
        end
        gf[1,1,:,:]=inv(3I-H+0.5im*weak);kf[1,1,:,:]=weak
        ac=covariance_A(gf,kf)[1,1,:,:]
        @test norm(ac-oracle)/norm(oracle)<1e-13
        @test real(tr(ac))>0 && eigmin(Hermitian(ac))>0
    end
end
function run_preflight(input,output)
    cfg=TOML.parsefile(input);cfg["runtime_authorized"]===true || error("operator-only implementation freeze")
    base=dirname(abspath(input));path(k)=normpath(joinpath(base,cfg[k]))
    for (p,s) in (("source_input_path","source_input_sha256"),("positive_seed_path","positive_seed_sha256"),("model_path","model_sha256"))
        E.checked_digest(path(p),cfg[s])
    end
    options=TOML.parsefile(path("source_input_path"));options["_input_directory"]=dirname(path("source_input_path"))
    core=dirname(dirname(pathof(QCLNEGF)))
    readchomp(`git -C $core rev-parse HEAD`)==options["expected_core_head"] || error("Core differs")
    readchomp(`git -C $core rev-parse $("HEAD^{tree}")`)==options["expected_core_tree"] || error("Core tree differs")
    string(VERSION)==options["expected_julia_version"] || error("Julia version differs")
    isempty(readchomp(`git -C $core status --porcelain --untracked-files=no`)) || error("Core dirty")
    Threads.nthreads()==4 || error("J4 required");BLAS.set_num_threads(1)
    member=only(filter(x->startswith(x,"0::"),readlines("/proc/self/cgroup")))
    cg=joinpath("/sys/fs/cgroup",lstrip(split(member,':';limit=3)[3],'/'))
    for (field,key) in (("memory.max","memory_max_bytes"),("memory.high","memory_high_bytes"))
        strip(read(joinpath(cg,field),String))==string(cfg[key]) || error("OS cap differs: $field")
    end
    strip(read(joinpath(cg,"memory.swap.max"),String))=="0" || error("Swap0 required")
    cpu=split(strip(read(joinpath(cg,"cpu.max"),String)))
    cpu[1]!="max" && parse(Int,cpu[1])/parse(Int,cpu[2])<=4 || error("CPU cap exceeds400%")
    limit=split(only(filter(line->startswith(line,"Max file size"),readlines("/proc/self/limits"))))
    limit[4]!="unlimited" && parse(Int,limit[4])<=cfg["output_max_bytes"] || error("LimitFSIZE exceeds32MiB")
    ispath(output) && error("Refusing overwrite");mkpath(output);cp(input,joinpath(output,"input.toml"));cp(@__FILE__,joinpath(output,"preflight.jl"))
    started=time();checkwall()=time()-started<cfg["wall_seconds"] || error("preflight wall cap")
    report=Dict{String,Any}("scope"=>"operator-only; new arithmetic identity, no AA/projector/solver/native acceptance",
        "source_sha256"=>E.file_digest(@__FILE__),"input_sha256"=>E.file_digest(input),"maps_started"=>0,"maps_completed"=>0)
    payload=Dict{String,Any}();chosen=Int[];c=nothing;gamma=nothing
    try
        tiny_checks();report["tiny_checks"]="pass"
        c0=E.build_context(options)
        HDF5.h5open(path("model_path"),"r") do f;E.verify_problem_model(c0,f);end
        gamma=HDF5.h5open(path("positive_seed_path"),"r") do f
            read(f["schema"])=="campaign133-warm-GammaK-ghost0-v1" || error("seed schema differs")
            read(f["energy_scale_eV"])==c0.problem.scales.E₀_eV && read(f["length_scale_m"])==c0.problem.scales.L₀_m && read(f["spin_degeneracy"])==c0.problem.physical.g_s || error("seed scales differ")
            for name in fieldnames(typeof(c0.problem.grids))
                actual=getfield(c0.problem.grids,name);actual isa BitArray && (actual=Int8.(actual))
                isequal(actual,read(f["axes/grid_"*string(name)])) || error("seed axes differ: $name")
            end
            E.read_complex(f,"Gamma_seed")::Field
        end
        seedcone=E.cone(gamma);seedcone.minimum_eigenvalue>0 || error("operator qualification requires the original strictly positive seed")
        report["input_Gamma_cone"]=Dict(string(k)=>v for (k,v) in pairs(seedcone))
        hH=copy(c0.cavity.hlocal);dh=norm(hH-c0.h);hscale=norm(c0.h)
        E.relative(dh,hscale)<=1e-10 || error("h symmetry restoration exceeds existing arithmetic target")
        cavity=N._cavity_production_plan(c0.problem,hH,c0.options)
        c=E.Context(c0.problem,c0.cache,c0.options,hH,cavity,c0.kBT,c0.target)
        report["h_correction"]=Dict("absolute"=>dh,"scale"=>hscale,"relative"=>E.relative(dh,hscale),"side_repeat_defect"=>norm(cavity.hlocal-hH))
        payload["h_raw"]=c0.h;payload["h_H"]=hH;payload["Tplus"]=c.problem.basis.T₊
        checkwall();data=spectrum(c,gamma);GC.gc(true)
        for k in (:GR,:A,:Asub,:physical,:effective,:localR,:totalR,:Lambda_raw,:Lambda_H_preassembly);payload[string(k)]=getproperty(data,k);end
        g=c.problem.grids;q=number_weights_A(data.A,g,c.problem.physical.g_s);qsub=number_weights_A(data.Asub,g,c.problem.physical.g_s)
        payload["q_cov"]=q;payload["q_sub"]=qsub
        chosen=sort(unique([argmin(q),argmin(qsub),length(q)]))
        report["q_cov"]=qinfo(q,c.target);report["q_sub"]=qinfo(qsub,c.target)
        for (name,A) in (("A_cov",data.A),("A_sub",data.Asub))
            witness=E.cone(A;require_nonzero=false);report[name*"_cone"]=Dict(string(k)=>v for (k,v) in pairs(witness));push!(chosen,witness.energy_index)
        end
        report["same_GR_spectral_difference"]=Dict(string(k)=>v for (k,v) in pairs(E.norms(data.A,data.Asub,g.wᴱ,g.wᵏ)))
        report["Lambda_projection"]=Dict(string(k)=>v for (k,v) in pairs(data.correction));report["PI_roundoff"]=data.roundoff
        report["maximum_GR_magnitude"]=maximum(abs,data.GR);report["maximum_Dyson_condition"]=maximum(data.condition)
        fitted=E.fit_mu(g.ε,q,c.target,c.kBT);push!(chosen,argmin(abs.(g.ε.-fitted.mu)))
        green,_=E.fd_green(data.GR,data.A,data.condition,data.scale,g.ε,fitted.mu,c.kBT)
        report["mu_cov"]=Dict(string(k)=>v for (k,v) in pairs(fitted))
        report["same_A_number_check"]=N._number_functional_production(green.Gˡ,g,c.problem.physical.g_s,c.options)/c.target-1
        checkwall();report["maps_started"]+=1;raw,channels,_=E.phi(c,green;keep_channels=true);report["maps_completed"]+=1;payload["Phi_cov_raw"]=raw
        phi_cone=E.cone(raw);report["Phi_cov_cone"]=Dict(string(k)=>v for (k,v) in pairs(phi_cone));push!(chosen,phi_cone.energy_index)
        report["raw_Gamma_residual"]=Dict(string(k)=>v for (k,v) in pairs(E.norms(raw,gamma,g.wᴱ,g.wᵏ)))
        try
            subfit=E.fit_mu(g.ε,qsub,c.target,c.kBT)
            control,_=E.fd_green(data.GR,data.Asub,data.condition,data.scale,g.ε,subfit.mu,c.kBT)
            checkwall();report["maps_started"]+=1;ctrl,_,_=E.phi(c,control);report["maps_completed"]+=1;payload["Phi_sub_raw"]=ctrl
            report["mu_sub"]=Dict(string(k)=>v for (k,v) in pairs(subfit));report["matched_Phi_difference"]=Dict(string(k)=>v for (k,v) in pairs(E.norms(raw,ctrl,g.wᴱ,g.wᵏ)))
            report["Phi_sub_cone"]=Dict(string(k)=>v for (k,v) in pairs(E.cone(ctrl)))
        catch err
            err isa E.Infeasible || rethrow();report["matched_sub_control"]=Dict("not_measured"=>true,"reason"=>sprint(showerror,err))
        end
        checkwall();report["maps_started"]+=1;report["native_comparability"]=native_compare(c,green,raw,channels,fitted.mu);report["maps_completed"]+=1
        report["execution"]="complete"
    catch err
        report["execution"]="stopped";report["reason"]=sprint(showerror,err,catch_backtrace());rethrow()
    finally
        if c!==nothing && gamma!==nothing && !isempty(chosen)
            chosen=sort(unique(chosen));g=c.problem.grids
            HDF5.h5open(joinpath(output,"selected.h5"),"w") do f
                f["schema"]="campaign133-positive-source-operator-selected-v1";f["energy"]=g.ε;f["weight_energy"]=g.wᴱ;f["weight_k"]=g.wᵏ
                f["selected_energy_indices"]=chosen;f["energy_scale_eV"]=c.problem.scales.E₀_eV;f["target"]=c.target;f["spin_degeneracy"]=c.problem.physical.g_s
                E.write_complex(f,"Gamma",gamma[chosen,:,:,:])
                for (key,X) in payload
                    if X isa Field;E.write_complex(f,key,X[chosen,:,:,:])
                    elseif X isa Array{ComplexF64};E.write_complex(f,key,X)
                    else;f[key]=X;end
                end
            end
        end
        report["elapsed_seconds"]=time()-started
        open(joinpath(output,"summary.toml"),"w") do io;TOML.print(io,report);end
        E.data_bytes(output)<=cfg["output_max_bytes"] || error("small-output cap exceeded; artifacts preserved")
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    ARGS==["--tiny"] ? tiny_checks() : (length(ARGS)==2 ? run_preflight(abspath(ARGS[1]),abspath(ARGS[2])) : error("preflight.jl --tiny OR INPUT NEWOUTPUT"))
end
