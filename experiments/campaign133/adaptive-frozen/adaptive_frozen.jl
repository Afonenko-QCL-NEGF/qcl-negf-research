using QCLNEGF, QuadGK, HDF5, LinearAlgebra, TOML, SHA, Dates, Unitful

function main()
    BLAS.set_num_threads(1)
    out=@__DIR__; src=joinpath(out,"../raw/S01-fixedU32/state.h5")
    base=joinpath(out,"../frozen-spectrum/integrals.h5")
    started=time(); startedutc=string(Dates.now(Dates.UTC))
    source_sha=bytes2hex(open(sha256,src))
    @assert source_sha=="155f8d296715f8e8a4ca951a0ed07aa392742112bd8b3b934beea54304374a2a"
    ne,nk,nb=30721,11,5; block=nb^2; half=block*nk; E0=.1
    mu=.008519932576353691/E0; kBT=QCLNEGF._electronvolts(CODATA.kᴮₑᵥ*70u"K")/E0
    target=.045000000000000005; gs=2.
    h,energy,wk,originalM,originalF,fineM,fineF=h5open(base,"r") do f
        z(g,name)=complex.(read(g[name*"_real"]),read(g[name*"_imag"]))
        (z(f,"h"),read(f["original_trapezoid_30721/energy_dimensionless"]),read(f["weight_k"]),
         z(f["original_trapezoid_30721"],"integral_A_over_2pi"),z(f["original_trapezoid_30721"],"integral_fA_over_2pi"),
         z(f["refined_trapezoid_61441"],"integral_A_over_2pi"),z(f["refined_trapezoid_61441"],"integral_fA_over_2pi"))
    end
    @assert bytes2hex(sha256(reinterpret(UInt8,vec(h))))=="8cde8ea024d9e690bf0396428574496734d7bad5962fc00617b9b380ae1194a7"
    delta=energy[2]-energy[1]
    seeds=collect(range(energy[1],energy[end];length=65))
    peaks=Vector{NamedTuple}()
    sigma=zeros(ComplexF64,ne,nk,nb,nb)
    families=["IFR","LO","acoustic","alloy","embedding","impurity"]
    h5open(src,"r") do f
        @assert read(f["energy"])==energy && read(f["weight_k"])==wk
        @assert all(iszero,read(f["hartree"]))
        for m=1:nk
            for name in families
                group=f["selfenergy/"*name]
                re=group["retarded_real"][:,m,:,:];iv=group["retarded_imag"][:,m,:,:]
                @inbounds for b=1:nb,a=1:nb,e=1:ne
                    sigma[e,m,a,b]+=complex(re[e,a,b],iv[e,a,b])
                end
            end
            Ar=f["green/A_real"][:,m,:,:]
            GLi=f["green/GL_imag"][:,m,:,:]
            traceA=[sum(Ar[e,a,a] for a=1:nb) for e=1:ne]
            traceFA=[QCLNEGF._fermi(energy[e],mu,kBT)*traceA[e] for e=1:ne]
            traceGn=[sum(GLi[e,a,a] for a=1:nb) for e=1:ne]
            for (label,v) in zip(("A","fA","actualGn"),(traceA,traceFA,traceGn))
                threshold=.001*maximum(v)
                indices=[e for e=2:ne-1 if v[e]>=v[e-1] && v[e]>v[e+1] && v[e]>=threshold]
                sort!(indices;by=e->v[e],rev=true)
                for e in Iterators.take(indices,3)
                    left=e;right=e
                    while left>1 && v[left]>.5v[e];left-=1;end
                    while right<ne && v[right]>.5v[e];right+=1;end
                    radius=max((energy[right]-energy[left])/2,2delta)
                    append!(seeds,clamp.([energy[e]-radius,energy[e],energy[e]+radius],energy[1],energy[end]))
                    push!(peaks,(k=m,label=label,index=e,energy_eV=energy[e]*E0,value=v[e],halfheight_left_eV=energy[left]*E0,halfheight_right_eV=energy[right]*E0,seed_radius_eV=radius*E0))
                end
            end
            println("source_k=",m," elapsed_s=",time()-started);flush(stdout)
        end
    end
    sort!(unique!(seeds))
    open(joinpath(out,"peaks.csv"),"w") do io
        println(io,"k,label,source_index,energy_eV,trace_height,halfheight_left_eV,halfheight_right_eV,seed_radius_eV")
        for row in peaks;println(io,join(values(row),","));end
    end
    scaleF=gs.*wk.*sqrt(nb)./target
    function science_norm(v)
        ma=0.;sf=0.
        for m=1:nk
            a=0.;b=0.
            @inbounds for j=(m-1)*block+1:m*block
                a+=abs2(v[j]);b+=abs2(v[half+j])
            end
            ma=max(ma,sqrt(a)/sqrt(nb));sf+=sqrt(b)
        end
        return max(ma,sf)
    end
    cache=Dict{Float64,Vector{ComplexF64}}()
    D=zeros(ComplexF64,nb,nb);S=similar(D)
    workspace=QCLNEGF.QCLNumerics._ProductionInverseWorkspace(nb)
    minimum_sv=fill(Inf,nk);max_condition=zeros(nk)
    calls=Ref(0)
    function integrand(en)
        calls[]+=1
        haskey(cache,en) && return cache[en]
        length(cache)<16384 || error("hard unique point cap reached before quadrature return")
        e=clamp(searchsortedlast(energy,en),1,ne-1)
        t=(en-energy[e])/(energy[e+1]-energy[e])
        @assert 0<=t<=1
        fd=QCLNEGF._fermi(en,mu,kBT)
        value=zeros(ComplexF64,2half)
        for m=1:nk
            @inbounds for b=1:nb,a=1:nb
                sig=(1-t)*sigma[e,m,a,b]+t*sigma[e+1,m,a,b]
                D[a,b]=(a==b ? en : 0.)-h[m,a,b]-sig
                S[a,b]=D[a,b]
            end
            singular=svdvals!(S)
            minimum_sv[m]=min(minimum_sv[m],singular[end])
            max_condition[m]=max(max_condition[m],singular[1]/singular[end])
            QCLNEGF.QCLNumerics._production_inverse!(D,workspace)
            @inbounds for b=1:nb,a=1:nb
                j=(m-1)*block+(b-1)*nb+a
                value[j]=im*(D[a,b]-conj(D[b,a]))/(2π)
                value[half+j]=scaleF[m]*fd*value[j]
            end
        end
        cache[en]=value
        return value
    end
    unpack(v)=begin
        M=zeros(ComplexF64,nk,nb,nb);F=similar(M)
        for m=1:nk,b=1:nb,a=1:nb
            j=(m-1)*block+(b-1)*nb+a
            M[m,a,b]=v[j];F[m,a,b]=v[half+j]/scaleF[m]
        end
        M,F
    end
    charge(F)=gs*sum(wk[m]*sum(real(F[m,a,a]) for a=1:nb) for m=1:nk)
    rsum(M)=maximum(norm(Matrix(view(M,m,:,:))-I)/sqrt(nb) for m=1:nk)
    compare(M,F,refM,refF)=Dict("max_k_M_absolute_Frobenius_over_sqrtNb"=>maximum(norm(view(M,m,:,:)-view(refM,m,:,:))/sqrt(nb) for m=1:nk),"max_k_M_relative_Frobenius"=>maximum(norm(view(M,m,:,:)-view(refM,m,:,:))/norm(view(refM,m,:,:)) for m=1:nk),"max_k_F_absolute_Frobenius_over_sqrtNb"=>maximum(norm(view(F,m,:,:)-view(refF,m,:,:))/sqrt(nb) for m=1:nk),"max_k_F_relative_Frobenius"=>maximum(norm(view(F,m,:,:)-view(refF,m,:,:))/norm(view(refF,m,:,:)) for m=1:nk),"global_M_relative_Frobenius"=>norm(M-refM)/norm(refM),"N_difference_over_target"=>(charge(F)-charge(refF))/target)
    results=Dict{String,Any}()
    summary=Dict{String,Any}("source_sha256"=>source_sha,"h_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(h)))),"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"core_head"=>readchomp(`git -C $(dirname(dirname(pathof(QCLNEGF)))) rev-parse HEAD`),"quadgk_version"=>string(pkgversion(QuadGK)),"quadgk_source"=>pathof(QuadGK),"quadgk_source_sha256"=>bytes2hex(sha256(read(pathof(QuadGK)))),"quadgk_adapt_sha256"=>bytes2hex(sha256(read(joinpath(dirname(pathof(QuadGK)),"adapt.jl")))),"project_sha256"=>bytes2hex(sha256(read(Base.active_project()))),"manifest_sha256"=>bytes2hex(sha256(read(joinpath(dirname(Base.active_project()),"Manifest.toml")))),"start_utc"=>startedutc,"mu_eV"=>mu*E0,"temperature_K"=>70.,"target_dimensionless"=>target,"E0_eV"=>E0,"kBT_dimensionless"=>kBT,"families_summed_once"=>families,"initial_boundaries"=>length(seeds),"retained_peaks"=>length(peaks),"rule"=>"QuadGK order7 Gauss/Kronrod7/15, vector integrand all11k, absolute diagnostic target, rtol0","science_norm"=>"max(max_k ||M_k||F/sqrt5, sum_k ||scaledF_k||F); scaledF_k=gs*wk*sqrt5/target*F_k. Triangle and |trF|<=sqrt5||F||F bound absolute N error/target by second norm. Both full A and fA matrices are integrated.","limitations"=>"embedded error heuristic; seeds selected from sampled source maxima cannot exclude hidden peaks; frozen piecewise linear SigmaR interpolant, no SCBA/LO/KK grid changes, not continuum/window/solver acceptance","results"=>results)
    h5open(joinpath(out,"integrals.h5"),"w") do f
        f["schema"]="campaign133-adaptive-frozen-v1"
        f["seed_boundaries_dimensionless"]=seeds;f["weight_k"]=wk;f["scaleF"]=scaleF
        f["domain_dimensionless"]=[energy[1],energy[end]]
    end
    seg=nothing
    for (stage,tol) in enumerate((1e-5,1e-6))
        before=length(cache);startcalls=calls[];passstart=time()
        initial_evals=stage==1 ? 15*(length(seeds)-1) : 15*length(seg)
        # QuadGK maxevals is approximate and may overshoot by at most one 30-call split.
        maxevals=stage==1 ? 10000 : initial_evals+(16384-length(cache))-30
        println("stage=",stage," atol=",tol," initial_evals=",initial_evals," maxevals=",maxevals," unique_before=",before);flush(stdout)
        Iq,err,nextseg = stage==1 ? quadgk_segbuf(integrand,seeds;atol=tol,rtol=0.,order=7,maxevals=maxevals,norm=science_norm) : quadgk_segbuf(integrand,energy[1],energy[end];atol=tol,rtol=0.,order=7,maxevals=maxevals,norm=science_norm,eval_segbuf=seg)
        seg=nextseg
        M,F=unpack(Iq)
        bounds=sort([(s.a,s.b,s.E) for s in seg];by=first)
        label=stage==1 ? "target_1e-5" : "target_1e-6"
        item=Dict{String,Any}("diagnostic_absolute_target"=>tol,"embedded_error_estimate"=>err,"diagnostic_target_reached"=>err<=tol,"status"=>err<=tol ? "estimated_target_reached" : "evaluation_budget_reached","unique_common_points_total"=>length(cache),"unique_common_points_added"=>length(cache)-before,"total_resolvent_inverses"=>nk*length(cache),"pass_function_calls_including_cache"=>calls[]-startcalls,"pass_seconds"=>time()-passstart,"panels"=>length(seg),"minimum_panel_width_eV"=>minimum(s.b-s.a for s in seg)*E0,"maximum_panel_width_eV"=>maximum(s.b-s.a for s in seg)*E0,"N_fA_dimensionless"=>charge(F),"N_over_target"=>charge(F)/target,"r_sum"=>rsum(M),"comparison_original30721"=>compare(M,F,originalM,originalF),"comparison_refined61441"=>compare(M,F,fineM,fineF))
        if stage==2
            oldM,oldF=h5open(joinpath(out,"integrals.h5"),"r") do f
                g=f["target_1e-5"]
                (complex.(read(g["M_real"]),read(g["M_imag"])),complex.(read(g["F_real"]),read(g["F_imag"])))
            end
            item["comparison_previous_target"]=compare(M,F,oldM,oldF)
        end
        results[label]=item
        # Record exact evaluation arithmetic of QuadGK evalrule, with positive Kronrod weights.
        xrule,wrule,_=kronrod(7)
        finalnodes=Float64[];finalweights=Float64[]
        for panel in bounds
            a,b,_=panel;st=.5*(b-a)
            push!(finalnodes,a+st);push!(finalweights,st*wrule[end])
            for j=1:length(xrule)-1
                push!(finalnodes,a+(1+xrule[j])*st);push!(finalweights,st*wrule[j])
                push!(finalnodes,a+(1-xrule[j])*st);push!(finalweights,st*wrule[j])
            end
        end
        check=zeros(ComplexF64,2half)
        @assert all(haskey(cache,x) for x in finalnodes)
        for (x,w) in zip(finalnodes,finalweights);check.+=w.*cache[x];end
        item["saved_rule_reduction_defect_science_norm"]=science_norm(check-Iq)
        item["final_quadrature_nodes"]=length(finalnodes)
        item["final_weight_sum_dimensionless"]=sum(finalweights)
        h5open(joinpath(out,"integrals.h5"),"r+") do f
            group=create_group(f,label)
            group["M_real"]=real.(M);group["M_imag"]=imag.(M);group["F_real"]=real.(F);group["F_imag"]=imag.(F)
            group["panel_bounds_and_error"]=reduce(hcat,collect.(bounds))
            group["all_unique_evaluated_points_dimensionless"]=sort(collect(keys(cache)))
            group["quadrature_nodes_dimensionless"]=finalnodes;group["quadrature_weights_dimensionless"]=finalweights
            group["minimum_sampled_D_singular_value_by_k"]=copy(minimum_sv)
            group["maximum_sampled_D_condition_number_by_k"]=copy(max_condition)
        end
        summary["unique_points_total"]=length(cache);summary["elapsed_s"]=time()-started
        summary["minimum_sampled_D_singular_value"]=minimum(minimum_sv);summary["maximum_sampled_D_condition_number"]=maximum(max_condition)
        open(joinpath(out,"summary.toml"),"w") do io;TOML.print(io,summary);end
        println("stage_result=",item);flush(stdout)
    end
    println("done_utc=",Dates.now(Dates.UTC)," elapsed_s=",time()-started," unique_points=",length(cache));flush(stdout)
end
main()
