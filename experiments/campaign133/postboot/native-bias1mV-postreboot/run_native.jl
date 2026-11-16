# Fixed-U full-Keldysh native pilot. No Runner, Matsubara, FD projection, warm import or Pkg.
const PROCESS_STARTED=time()
using QCLNEGF, Unitful, LinearAlgebra, HDF5, TOML, SHA, Dates, Test
const N=QCLNEGF.QCLNumerics
const BASE=normpath(joinpath(@__DIR__, "..", "..", "equilibrium-gamma-stableA-cont64", "attempt01"))
sha(p)=bytes2hex(open(sha256,p))
for (f,d) in (("gamma_functions.jl","50bcfaaeb3c399baaaa2a4b9b9e023167e97b4d87c6af95d83353cb453072211"),
 ("warm_seed.jl","68fd8582334a1ed3f257efcc0b1096a429361dd73aa5be3be415dd79ddcbc157"),
 ("stable_spectral.jl","eeca1f29710d9bb60bf6ce70ce0ebd1620fb86d4c77314b4bd8ffa9a5bcd5aa5"),
 ("gamma_restart.jl","2a805be29aeef33d7a31356850b7fdaa80da7e68506d8f96b980a8108bd56a00"))
    sha(joinpath(BASE,f))==d || error("serialization helper differs")
end
include(joinpath(BASE,"gamma_functions.jl"));const GP=EquilibriumGammaPilot
plain(x::Union{Symbol,AbstractString,Number})=x isa Symbol ? string(x) : x
plain(x::Complex)=Dict("real"=>real(x),"imag"=>imag(x))
plain(x::Nothing)="not_measured"
plain(x::AbstractVector)=plain.(x)
plain(x::AbstractDict)=Dict{String,Any}(string(k)=>plain(v) for (k,v) in x)
plain(x::AbstractMatrix)=Dict("shape"=>collect(size(x)),"rows"=>[[plain(x[i,j]) for j in axes(x,2)] for i in axes(x,1)])
plain(x)=Dict{String,Any}(string(f)=>plain(getfield(x,f)) for f in fieldnames(typeof(x)))
write_toml(p,x)=open(io->TOML.print(io,x),p,"w")
stop_due(c)=time()-PROCESS_STARTED>=c["overall_soft_seconds"]-c["save_reserve_seconds"] || Dates.now(Dates.UTC)>=DateTime(2026,11,16,13)-Second(c["save_reserve_seconds"])
function policy(c;tiny=false)
    a=AlgorithmOptions(kernel_build=:direct,embedding=:finite_chain,embedding_periods=4,mixing=:linear,
        energy_shift=:sparse_plan,hilbert=:product_integration,elastic_retarded=:direct_green)
    po=ProductionOptions(algorithms=a,worker_count=4,hilbert_columns=1,parallel_backend=:threads,
        memory_budget_bytes=c["memory_max_bytes"]-c["reserve_bytes"],verify_fft_roundoff=true,
        checkpoint_every_scba=0,checkpoint_every_outer=0,progress_every_scba=1,
        physics_markers=SCBAPhysicsMarkerPolicy(cadence=1,max_spectral_blocks=64))
    n=NumericalParameters(N_z=tiny ? 32 : 96,N_b=tiny ? 2 : 5,P_basis=3,E_min=-0.3u"eV",E_max=0.9u"eV",N_E=tiny ? 65 : 24001,
        M_E=80u"meV",k_max=0.6u"nm^-1",N_k=tiny ? 3 : 11,N_φ=tiny ? 8 : 32,qz_max=10u"nm^-1",N_qz=tiny ? 9 : 81,η_seed=1.25u"meV")
    so=SolverOptions(α_Σ=0.2,max_scba=tiny ? 2 : c["max_scba"],max_poisson=1,
        convergence=ConvergencePolicy(stagnation_window=0,stagnation_relative_improvement=0.0,
            minimum_scba_iterations=2,required_consecutive_scba_passes=3))
    return a,po,n,so
end
function gate(c;actual=true)
    core=dirname(dirname(pathof(QCLNEGF)))
    readchomp(`git -C $core rev-parse HEAD`)==c["core_head"] || error("Core HEAD")
    readchomp(`git -C $core rev-parse $("HEAD^{tree}")`)==c["core_tree"] || error("Core tree")
    isempty(readchomp(`git -C $core status --porcelain --untracked-files=no`)) || error("Core dirty")
    string(VERSION)=="1.13.4" && Threads.nthreads()==4 || error("Julia/J4 identity")
    c["voltage_mV"] in (0.0,1.0) && c["max_scba"]==32 || error("bounded matched case")
    c["memory_max_bytes"]==47*1024^3÷4 && c["reserve_bytes"]==1400*1024^2 || error("resource envelope")
    BLAS.set_num_threads(1);BLAS.get_num_threads()==1 || error("BLAS1")
    if actual
        c["runtime_authorized"]===true || error("implementation-only freeze")
        member=only(filter(x->startswith(x,"0::"),readlines("/proc/self/cgroup")))
        cg=joinpath("/sys/fs/cgroup",lstrip(split(member,':';limit=3)[3],'/'))
        for (f,k) in (("memory.max","memory_max_bytes"),("memory.high","memory_high_bytes"),("memory.swap.max","swap_max_bytes"))
            strip(read(joinpath(cg,f),String))==string(c[k]) || error("OS $f")
        end
        cpu=split(strip(read(joinpath(cg,"cpu.max"),String)))
        cpu[1]!="max" && parse(Int,cpu[1])/parse(Int,cpu[2])<=4 || error("CPU400%")
        lim=split(only(filter(x->startswith(x,"Max file size"),readlines("/proc/self/limits"))))
        lim[4]!="unlimited" && parse(Int,lim[4])<=c["single_file_max_bytes"] || error("OS file cap")
    end
    return core
end
function plan(c;tiny=false)
    a,po,n,so=policy(c;tiny);nodes=collect(range(-3.0,9.0;length=n.N_E))
    e=estimate_production_memory(n,5;dense_mechanism_count=3,options=po,worker_capacity=4,energy_nodes=nodes)
    e.peak_bytes<=po.memory_budget_bytes || error("native peak+reserve exceeds hard cap")
    return e
end
struct AlignedStop <: Exception
    state::SCBAResult
end
function save_terminal(p,r,out,c)
    shape=size(r.green.A);expected=29*GP.expected_padded_complex_bytes(shape)+2*8*shape[1]*shape[2]+32*1024^2
    GP.atomic_h5(joinpath(out,"state.h5"),c;expected_bytes=expected) do f
        f["schema"]="campaign133-native-bias-fullstate-v1";f["voltage_per_period_mV"]=c["voltage_mV"]
        f["status"]=string(r.status);f["converged"]=r.converged;f["quality"]=string(r.quality)
        f["energy"]=p.grids.ε;f["weight_energy"]=p.grids.wᴱ;f["weight_k"]=p.grids.wᵏ;f["hartree"]=zeros(p.numerical.N_z)
        f["energy_scale_eV"]=p.scales.E₀_eV;f["length_scale_m"]=p.scales.L₀_m;f["spin_degeneracy"]=p.physical.g_s
        GP.write_green(f,r.green);families=copy(r.scattering)
        merge!(families,Dict(:embedding=>r.embedding,:embedding_plus=>r.embedding_plus,:embedding_minus=>r.embedding_minus))
        group=HDF5.create_group(f,"selfenergy");gamma=zeros(ComplexF64,shape)
        for (id,family) in families
            sub=HDF5.create_group(group,string(id))
            for (name,X) in (("retarded",family.Σᴿ),("lesser",family.Σˡ),("greater",family.Σᵍ));GP.write_complex(sub,name,X);end
            if id in p.kernels.enabled
                @inbounds @simd for i in eachindex(gamma);gamma[i]+=im*(family.Σᵍ[i]-family.Σˡ[i]);end
            end
        end
        GP.write_complex(f,"Gamma_local_kinetic",gamma)
        f["restart_contract_toml"]=sprint(io->TOML.print(io,r.restart_contract))
    end
end
function run_case(input,out;tiny=false)
    c=TOML.parsefile(input);core=gate(c;actual=!tiny);estimate=plan(c;tiny)
    ispath(out) && error("new output required");mkpath(out);cp(input,joinpath(out,"input.toml"));cp(@__FILE__,joinpath(out,"run_native.jl"))
    stop_due(c) && error("no new case: overall/UTC13 save reserve")
    started=time();a,po,n,so=policy(c;tiny);summary=Dict{String,Any}("execution"=>"started","scientific_acceptance"=>false,
        "initialization"=>"native cold seed; FD only initialization, full Keldysh updates thereafter","EE"=>"none","maps_completed"=>0)
    write_toml(joinpath(out,"identity.toml"),Dict("core"=>core,"core_head"=>c["core_head"],"input_sha256"=>sha(input),
        "script_sha256"=>sha(@__FILE__),"project"=>Base.active_project(),"manifest_sha256"=>sha(joinpath(dirname(Base.active_project()),"Manifest.toml")),"estimate_peak_bytes"=>estimate.peak_bytes,"native_budget_bytes"=>po.memory_budget_bytes))
    try
        physical=reference_parameters(V_period=0u"mV",Tᴸ=70u"K",Tᴸᴼ=70u"K",ΔV_alloy=.60u"eV",Ω₀=4.5168e-29u"m^3")
        models=PhysicalModelOptions(electron_electron=ElectronElectronOptions(mode=:none))
        built=build_configured_scattering_problem(physical=physical,numerical=n,scales=ScaleSystem(),scattering=ScatteringOptions(alloy=true),
            algorithms=a,kernel_options=ProductionKernelOptions(),physical_models=models)
        p0=built.problem;p=retarget_problem(p0;V_period=c["voltage_mV"]*u"mV",energy_shift=:sparse_plan)
        Set(p.kernels.enabled)==Set((:LO,:acoustic,:impurity,:IFR,:alloy)) && p.models.electron_electron.mode===:none || error("all5/noEE")
        sp=N._scaled_physics(p.physical,p.scales);isapprox(sp.Eᵖ,c["voltage_mV"]*1e-3/p.scales.E₀_eV;rtol=1e-14,atol=0) || error("period bias units")
        p.basis===p0.basis && p.kernels.K===p0.kernels.K || error("fixed zero-field basis/kernel reuse")
        cache=build_production_cache(p;options=po);h=project_hamiltonians(p,zeros(n.N_z));ctx=GP.Context(p,cache,po,h,nothing,Float64(ustrip(u"eV",CODATA.kᴮₑᵥ*p.physical.Tᴸ))/p.scales.E₀_eV,sp.Nᴰ²ᴰ)
        GP.save_model(ctx,out,c);summary["model_sha256"]=sha(joinpath(out,"model.h5"))
        po=N.with_production_options(po;checkpoint_request=req->req.stage===:scba) # borrowed aligned state, no retained copy
        fields=[f for f in fieldnames(SCBAIteration) if fieldtype(SCBAIteration,f)<:Real]
        result=open(joinpath(out,"history.csv"),"w") do io
            println(io,join(string.(fields),","))
            observer=row->begin
                println(io,join([string(getfield(row,f)) for f in fields],","));flush(io)
                summary["maps_completed"]+=1
                row.physical_markers!==nothing && write_toml(joinpath(out,"markers-$(lpad(row.ν,3,'0')).toml"),plain(row.physical_markers))
                row.witness!==nothing && write_toml(joinpath(out,"witness-$(lpad(row.ν,3,'0')).toml"),plain(row.witness))
            end
            try
                solve_scba_production(p,zeros(n.N_z);options=so,production_options=po,cache=cache,
                    history_observer=observer,iteration_callback=state->begin
                        row=last(state.history);plus=N._channel_collision_marker(p,state.green,:embedding_plus,:mixed,state.embedding_plus)
                        minus=N._channel_collision_marker(p,state.green,:embedding_minus,:mixed,state.embedding_minus)
                        internal=filter(x->x.state_kind===:mixed,row.physical_markers.collisions)
                        signed=sum(x->x.particle_signed,internal)+plus.particle_signed+minus.particle_signed
                        scale=sum(x->x.particle_absolute,internal)+plus.particle_absolute+minus.particle_absolute
                        powerunit=p.scales.E₀_J*p.scales.E₀_eV/(Float64(ustrip(u"eV*s",CODATA.ħₑᵥ))*p.scales.L₀_m^2)
                        Psc=sum(x->x.energy_signed,internal)*powerunit;Pf=row.J*c["voltage_mV"]*1e-3
                        Ps=sum(x->abs(x.energy_signed),internal)*powerunit+abs(Pf)
                        write_toml(joinpath(out,"balance-$(lpad(row.ν,3,'0')).toml"),Dict("plus"=>plain(plus),"minus"=>plain(minus),
                            "particle_signed"=>signed,"particle_scale"=>scale,"particle_relative"=>GP.relative(abs(signed),scale),
                            "measure"=>"gs*wE*wk/(2pi), mixed internal+embedding signed sum; scaled rate; no hbar conversion, diagnostic not acceptance",
                            "J_unit"=>"SCBAIteration.J A/m^2","Psc_W_m2"=>Psc,"Pfield_W_m2"=>Pf,
                            "power_balance_abs_W_m2"=>abs(Psc+Pf),"power_scale_W_m2"=>Ps,"power_relative_diagnostic"=>GP.relative(abs(Psc+Pf),Ps)))
                        (stop_due(c) || time()-started>=c["solver_seconds"]) && throw(AlignedStop(state))
                    end)
            catch e
                e isa AlignedStop || rethrow();e.state
            end
        end
        GC.gc(true);save_terminal(p,result,out,c)
        summary["execution"]="complete";summary["status"]=string(result.status);summary["converged"]=result.converged
        summary["quality"]=string(result.quality);summary["state_sha256"]=sha(joinpath(out,"state.h5"))
        summary["last_row"]=Dict(string(f)=>getfield(last(result.history),f) for f in fields)
    catch e
        summary["execution"]="stopped";summary["reason"]=sprint(showerror,e,catch_backtrace());rethrow()
    finally
        summary["elapsed_seconds"]=time()-started;summary["data_bytes"]=GP.data_bytes(out)
        write_toml(joinpath(out,"summary.toml"),summary)
    end
    return summary # only scalars/strings; every full native object dies at function barrier
end
function run_matched(input,out)
    c=TOML.parsefile(input);ispath(out) && error("new matched output required");mkpath(out);cp(input,joinpath(out,"input.toml"))
    cases=Dict{String,Any}[];summary=Dict{String,Any}("cases"=>cases,"scientific_acceptance"=>false)
    try
        c["case_voltages_mV"]==[0.0,1.0] || error("matched V0 then V1 required")
        for voltage in c["case_voltages_mV"]
            stop_due(c) && error("no new case before UTC13/overall save reserve")
            cfg=copy(c);cfg["voltage_mV"]=voltage;cfg["identity"]=c["identity"]*"-V$(Int(voltage))"
            cfg["data_max_bytes"]=min(c["data_max_bytes"],c["overall_data_max_bytes"]-GP.data_bytes(out)-1024^2)
            path=joinpath(out,"case$(Int(voltage))-input.toml");write_toml(path,cfg)
            metadata=run_case(path,joinpath(out,"case$(Int(voltage))"));push!(cases,metadata);GC.gc(true)
        end
        GP.data_bytes(out)<=c["overall_data_max_bytes"] || error("combined8GiB cap; raw preserved")
        summary["execution"]="complete"
    catch e
        summary["execution"]="stopped";summary["reason"]=sprint(showerror,e,catch_backtrace());rethrow()
    finally
        summary["elapsed_process_seconds"]=time()-PROCESS_STARTED;summary["bytes"]=GP.data_bytes(out)
        write_toml(joinpath(out,"summary.toml"),summary)
    end
end
function tiny_checks()
    c=TOML.parsefile(joinpath(@__DIR__,"matched.toml"))
    @testset "native all5 full-Keldysh small V0/V1 IO" begin
        mktempdir() do dir
            for v in (0.0,1.0)
                cfg=copy(c);cfg["voltage_mV"]=v;input=joinpath(dir,"i$(v).toml");write_toml(input,cfg)
                out=joinpath(dir,"v$(v)");meta=run_case(input,out;tiny=true)
                @test meta["maps_completed"]==2
                h5open(joinpath(out,"state.h5"),"r") do f
                    @test read(f["schema"])=="campaign133-native-bias-fullstate-v1"
                    @test length(keys(f["selfenergy"]))==8
                    A=GP.read_complex(f["green"],"A");L=GP.read_complex(f["green"],"GL");G=GP.read_complex(f["green"],"GG")
                    @test size(A)==(65,3,2,2) && norm(L)>0 && norm(G)>0
                    @test norm(im.*(G.-L).-A)/norm(A)<1e-10
                    @test all(isfinite,GP.read_complex(f,"Gamma_local_kinetic"))
                end
                @test isfile(joinpath(out,"model.h5")) && isfile(joinpath(out,"balance-002.toml"))
            end
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS==["--tiny"];tiny_checks()
    elseif length(ARGS)==2 && ARGS[1]=="--plan"
        c=TOML.parsefile(ARGS[2]);c["voltage_mV"]=0.0;gate(c;actual=false);show(stdout,plan(c));println()
    else
        length(ARGS)==2 || error("--plan INPUT OR INPUT NEWOUT");run_matched(abspath(ARGS[1]),abspath(ARGS[2]))
    end
end
