# Included by run_control.jl; frozen Core primitives, private map-only admission.
record(x)=Dict{String,Any}(string(k)=>v for (k,v) in pairs(x))
fieldbytes(p)=16*p.numerical.N_E*p.numerical.N_k*p.numerical.N_b^2
function gamma_cache(p,o)
    full=estimate_production_memory(p;options=o) # native scope kept honestly
    envelope=24fieldbytes(p)+3*1024^3
    envelope<=o.memory_budget_bytes || error("Gamma-only envelope exceeds real cap")
    o.algorithms.hilbert===:product_integration && o.algorithms.energy_shift===:sparse_plan || error("cache policy")
    N._uniform_product_integration_geometry(p.grids.ε)!==nothing || error("uniform FFT required")
    ops=Dict{Symbol,N.AbstractProductionKernel}(s=>N._production_kernel(p.kernels.K[s],o.algorithms) for s in p.kernels.enabled)
    sp=N._scaled_physics(p.physical,p.scales)
    shifts=map(d->build_shift_plan(p.grids.ε,d),(sp.Eᵖ,-sp.Eᵖ,sp.ħωᴸᴼ,-sp.ħωᴸᴼ))
    plan=N._ProductionPIFFTHilbertPlan(p.grids.ε,p.numerical.N_k*p.numerical.N_b^2;column_chunk=1,worker_count=1)
    cache=ProductionCache(shifts...,ops,nothing,plan,full,p.grids,p.kernels.K,o.algorithms)
    N._cache_source_contract(cache,p); N._cache_shift_contract(cache,p)
    payload=sum(N._kernel_operator_bytes(v) for v in values(ops))+
        production_fft_workspace_bytes(plan)+production_fft_workspace_bytes(plan;verification=true)+
        sum(Base.summarysize(s) for s in shifts)
    payload<256*1024^2 || error("unexpected cache payload")
    ledger=Dict("scope"=>"Gamma-only two-map; no native exit", "native_full_peak_bytes"=>full.peak_bytes,
        "Gamma_envelope_bytes"=>envelope,"cap_bytes"=>o.memory_budget_bytes,"cache_payload_bytes"=>payload,
        "full_field_bytes"=>fieldbytes(p),"history_fields"=>0,"native_guard_unchanged"=>true)
    return cache,ledger
end
function cache_equal(a,b,p,G)
    for name in (:W̃₊ᴱᵖ,:W̃₋ᴱᵖ,:W̃₊ᴸᴼ,:W̃₋ᴸᴼ)
        x,y=getfield(a,name),getfield(b,name)
        all(getfield(x,s)==getfield(y,s) for s in fieldnames(typeof(x))) || error("shift cache differs")
    end
    for s in p.kernels.enabled
        x,y=a.kernels[s],b.kernels[s]
        typeof(x)===typeof(y) && all(getfield(x,k)==getfield(y,k) for k in fieldnames(typeof(x))) || error("kernel cache differs")
    end
    L,ra=N._production_pi_fft_hilbert_transform(G,p.grids.ε,p.grids.wᴱ,a.product_hilbert_operator)
    M,rb=N._production_pi_fft_hilbert_transform(G,p.grids.ε,p.grids.wᴱ,b.product_hilbert_operator)
    L==M || error("coarse public/private FFT values differ")
    return L,Dict("exact_value_equality"=>true,"native_roundoff"=>ra,"private_roundoff"=>rb)
end
function refined_context(c)
    p=c.problem; a=p.numerical
    n=NumericalParameters(a.N_z,a.N_b,a.P_basis,a.E_min,a.E_max,2*a.N_E-1,a.M_E,a.k_max,a.N_k,a.N_φ,a.qz_max,a.N_qz,a.η_seed)
    g=build_grids(p.physical,n,p.scales)
    g.ε[1:2:end]==p.grids.ε || error("canonical refinement not bitwise nested")
    for s in (:x,:wˣ,:κ,:wᵏ,:qᶻ,:wᑫᶻ,:φ); getfield(g,s)==getfield(p.grids,s) || error("non-energy axis changed"); end
    sp=N._scaled_physics(p.physical,p.scales)
    shifts=map(d->build_shift_plan(g.ε,d),(sp.Eᵖ,-sp.Eᵖ,sp.ħωᴸᴼ,-sp.ħωᴸᴼ))
    q=NEGFProblem(p.physical,n,p.scattering,p.scales,g,p.profiles,p.basis,p.kernels,shifts...,:nodal_linear,p.models)
    cache,ledger=gamma_cache(q,c.options)
    refused=false
    try; N._validated_production_cache(cache,q,c.options)
    catch e
        e isa ArgumentError && occursin("estimated peak memory",sprint(showerror,e)) || rethrow()
        refused=true
    end
    full_exceeds=cache.estimate.peak_bytes>c.options.memory_budget_bytes
    refused==full_exceeds || error("native admission behavior differs")
    ledger["native_validation_refused_for_full_budget"]=refused
    cavity=N._cavity_production_plan(q,c.h,c.options)
    ledger["cavity_plan_summary_bytes"]=Base.summarysize(cavity)
    ledger["cache_payload_bytes"]+ledger["cavity_plan_summary_bytes"]<=512*1024^2 || error("unexpected cache/cavity payload")
    return GP.Context(q,cache,c.options,c.h,cavity,c.kBT,c.target),ledger
end
function prolong_rows(f,name,rows,source_energy,target_energy)
    nsource=length(source_energy)
    lo=(first(rows)+1)÷2; hi=cld(last(rows)+1,2)
    1<=lo<=hi<=nsource || error("prolong row bounds")
    X=GP.read_complex_rows(f,name,lo:hi); Y=zeros(ComplexF64,length(rows),size(X)[2:end]...)
    for (i,e) in enumerate(rows)
        j=(e+1)÷2-lo+1
        if isodd(e); Y[i,:,:,:].=X[j,:,:,:]
        else
            global_j=e÷2
            theta=(target_energy[e]-source_energy[global_j])/(source_energy[global_j+1]-source_energy[global_j])
            isfinite(theta) && 0<=theta<=1 || error("nonconvex actual-coordinate interpolation")
            Y[i,:,:,:].=(1-theta).*X[j,:,:,:].+theta.*X[j+1,:,:,:]
        end
    end
    return Y
end
function prolong_field(f,name,source_energy,target_energy)
    Y=zeros(ComplexF64,length(target_energy),size(f[name*"_real"])[2:end]...)
    for s in 1:1024:size(Y,1)
        rows=s:min(s+1023,size(Y,1)); Y[rows,:,:,:].=prolong_rows(f,name,rows,source_energy,target_energy)
    end
    return Y
end
function one_map(c,G,mu0,started,checkpoint)
    started[]<2 || error("two-map start cap"); started[]+=1; checkpoint("Pi_started")
    t=time(); input=GP.cone(G)
    GR,condition,scale,roundoff=GP.retarded_state(c,G); checkpoint("retarded_complete")
    A=spectral_function(GR); acone=GP.cone(A;require_nonzero=false)
    q=N._seed_number_weights_production(GR,c.problem.grids,c.problem.physical.g_s;options=c.options)
    fixedN,fixedDerivative=GP.number_mu(c.problem.grids.ε,q,mu0,c.kBT) # BEFORE refit
    fit=GP.fit_mu(c.problem.grids.ε,q,c.target,c.kBT)
    green,pair=GP.fd_green(GR,A,condition,scale,c.problem.grids.ε,fit.mu,c.kBT)
    Phi,channels,sizes=GP.phi(c,green;keep_channels=false); isempty(channels) || error("no retained channels")
    raw=GP.norms(Phi,G,c.problem.grids.wᴱ,c.problem.grids.wᵏ)
    stats=Dict("seconds"=>time()-t,"raw_residual"=>record(raw),"mu"=>record(fit),
        "fixed_mu_diagnostic"=>Dict("mu"=>mu0,"N"=>fixedN,"dN_dmu"=>fixedDerivative,"relative_error"=>(fixedN-c.target)/c.target),
        "input_cone"=>record(input),"Phi_cone"=>record(GP.cone(Phi)),"A_cone"=>record(acone),
        "PI_roundoff"=>roundoff,"FD_pair_error"=>pair,"pair_A_error"=>GP.pair_A_error(green),
        "max_condition"=>maximum(condition),"min_Dyson_scale"=>minimum(scale),
        "channel_Gamma_norms"=>sizes,"stencil"=>GP.stencil_witness(c,fit.mu))
    checkpoint("Pi_complete"); return Phi,stats # GR/FD/cavity arrays cannot escape
end
function compare_file(X,f,name,we,wk;restricted=false)
    d=GP.WarmNorm(); n=restricted ? (size(X,1)+1)÷2 : size(X,1)
    for s in 1:1024:n
        rows=s:min(s+1023,n); ref=GP.read_complex_rows(f,name,rows)
        Y=restricted ? X[2 .* rows .- 1,:,:,:] : X[rows,:,:,:]
        GP.accumulate_norm!(d,ref,Y,we,wk,s-1)
    end
    return GP.warm_norm_dict(d)
end
function save_fields(path,cfg,pairs;energy,weights)
    bytes=sum(GP.expected_padded_complex_bytes(size(A)) for (_,A) in pairs)+16*1024^2
    GP.atomic_h5(path,cfg;expected_bytes=bytes) do f
        f["energy"]=energy; f["weight_energy"]=weights
        for (name,A) in pairs; GP.write_complex(f,name,A); end
    end
end
