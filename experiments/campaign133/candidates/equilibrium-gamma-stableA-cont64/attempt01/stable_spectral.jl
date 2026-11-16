# Experimental finite-campaign spectral functional. Native Context is immutable.
# covariance_A and number_weights_A copied byte-for-byte from qualified preflight.
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

function stable_context(c::Context)
    hH=copy(c.cavity.hlocal) # exact native F64 sidecell policy, now central too
    cavity=N._cavity_production_plan(c.problem,hH,c.options)
    isequal(cavity.hlocal,hH) || error("hH all9 side repeat differs")
    return Context(c.problem,c.cache,c.options,hH,cavity,c.kBT,c.target)
end
function diagnostic_cone(X::Field)
    try
        witness=cone(X;require_nonzero=false)
        return Dict{String,Any}("measured"=>true,"admissible"=>true,
            (string(k)=>v for (k,v) in pairs(witness))...)
    catch err
        err isa Infeasible || rethrow()
        return Dict{String,Any}("measured"=>true,"admissible"=>false,"reason"=>sprint(showerror,err))
    end
end
function stable_q_info(q,target)
    negative=count(<(0.0),q);mass=-sum(x->min(x,0.0),q)
    return Dict{String,Any}("minimum"=>minimum(q),"minimum_index"=>argmin(q),
        "negative_count"=>negative,"negative_absolute_mass"=>mass,"negative_over_target"=>mass/target,
        "capacity"=>sum(q),"all_finite"=>all(isfinite,q))
end
function _stable_spectral_state(c::Context,gamma::Field;output=nothing,cfg=nothing,call=0)
    s=stable_context(c);g=c.problem.grids
    Lraw,roundoff=N._production_pi_fft_hilbert_transform(gamma,g.ε,g.wᴱ,c.cache.product_hilbert_operator)
    LH=hermitian_copy(Lraw);correction=norms(Lraw,LH,g.wᴱ,g.wᵏ)
    correction.relative<=1e-10 && correction.weighted_relative<=1e-10 ||
        throw(Infeasible("PI real-part symmetry correction exceeds original arithmetic target"))
    dh=norm(s.h-c.h);hscale=norm(c.h)
    relative(dh,hscale)<=1e-10 || throw(Infeasible("hH correction exceeds original arithmetic target"))
    R=copy(LH)
    @inbounds @simd for i in eachindex(R);R[i]-=im*gamma[i]/2;end
    probe=SelfEnergyFamily(R,im.*gamma,zeros(ComplexF64,size(gamma)))
    emb,plus,minus=N._cavity_embedding_production(c.problem,probe,s.h,c.options;plan=s.cavity)
    totalR=R+emb.Σᴿ;physical=gamma.-im.*emb.Σˡ
    GR,condition,scale=N._retarded_green_production(g.ε,s.h,totalR;η=0.0,options=c.options)
    all(isfinite,GR) && all(isfinite,condition) && all(isfinite,scale) || throw(Infeasible("nonfinite Dyson inverse/condition"))
    A=covariance_A(GR,physical);Asub=spectral_function(GR)
    q=number_weights_A(A,g,c.problem.physical.g_s);qsub=number_weights_A(Asub,g,c.problem.physical.g_s)
    effective=similar(physical)
    @inbounds for b in axes(GR,4),a in axes(GR,3),m in axes(GR,2),e in axes(GR,1)
        effective[e,m,a,b]=im*(s.h[m,a,b]+totalR[e,m,a,b]-conj(s.h[m,b,a]+totalR[e,m,b,a]))
    end
    ledger=Dict{String,Any}("operator"=>"native_covariance_exact_plan_hlocal_all9",
        "native_exit_h_scope"=>"original raw-h Context; reduced hH policy is not a native gate override",
        "covariance_arithmetic"=>"generic F64 congruence, not certified Gram factor arithmetic",
        "h_correction"=>Dict("absolute"=>dh,"scale"=>hscale,"relative"=>relative(dh,hscale)),
        "Lambda_projection"=>Dict(string(k)=>v for (k,v) in pairs(correction)),
        "A_cov_cone"=>diagnostic_cone(A),"A_sub_cone"=>diagnostic_cone(Asub),
        "q_cov"=>stable_q_info(q,c.target),"q_sub"=>stable_q_info(qsub,c.target),
        "same_GR_spectral_difference"=>Dict(string(k)=>v for (k,v) in pairs(norms(A,Asub,g.wᴱ,g.wᵏ))),
        "effective_vs_physical_source"=>Dict(string(k)=>v for (k,v) in pairs(norms(effective,physical,g.wᴱ,g.wᵏ))),
        "maximum_GR_magnitude"=>maximum(abs,GR),"roundoff"=>roundoff)
    if output!==nothing
        chosen=sort(unique([argmin(q),argmin(qsub),length(q)]))
        for key in ("A_cov_cone","A_sub_cone")
            get(ledger[key],"admissible",false) && push!(chosen,ledger[key]["energy_index"])
        end
        chosen=sort(unique(chosen))
        selected_shape=(length(chosen),size(gamma)[2:end]...)
        expected=10*expected_padded_complex_bytes(selected_shape)+4*1024^2
        atomic_h5(joinpath(output,"spectral-selected-$(lpad(call,3,'0')).h5"),cfg;expected_bytes=expected) do f
            f["schema"]="campaign133-stableA-spectral-selected-v1";f["selected_energy_indices"]=chosen
            f["energy"]=g.ε;f["weight_energy"]=g.wᴱ;f["weight_k"]=g.wᵏ
            f["q_cov"]=q;f["q_sub"]=qsub;f["energy_scale_eV"]=c.problem.scales.E₀_eV
            write_complex(f,"h_raw",c.h);write_complex(f,"h_H",s.h);write_complex(f,"Tplus",c.problem.basis.T₊)
            for (name,X) in (("Gamma",gamma),("Lambda_raw",Lraw),("Lambda_H_preassembly",LH),
                            ("localR",R),("totalR",totalR),("GR",GR),("A_cov",A),("A_sub",Asub),
                            ("Gamma_physical_total",physical),("Gamma_eff",effective))
                write_complex(f,name,X[chosen,:,:,:])
            end
            f["ledger_toml"]=sprint(io->TOML.print(io,ledger))
        end
        open(joinpath(output,"spectral-ledger-$(lpad(call,3,'0')).toml"),"w") do io;TOML.print(io,ledger);end
    end
    return GR,A,condition,scale,roundoff,q,ledger
end
stable_spectral_state(c::Context,gamma::Field;kwargs...) = inverse_domain(()->_stable_spectral_state(c,gamma;kwargs...))
