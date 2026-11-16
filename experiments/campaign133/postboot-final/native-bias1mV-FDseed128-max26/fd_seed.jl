# Experimental ONE-TIME global FD INITIAL GUESS, including at V1. No map projection.
array_sha(X::Array{ComplexF64,4})=bytes2hex(sha256(reinterpret(UInt8,vec(X))))
function seed_cone(X::Array{ComplexF64,4};retarded=false)
    ne,nk,nb,_=size(X);nw=min(Threads.nthreads(),ne*nk)
    reports=Vector{Dict{String,Any}}(undef,nw)
    Threads.@threads :static for worker in 1:nw
        B=zeros(ComplexF64,nb,nb);H=similar(B);D=similar(B)
        low=Inf;worst=0.0;maxherm=0.0;bad=0;where_min=[0,0];where_bad=[0,0];minscale=0.0;minbudget=0.0
        for index in worker:nw:ne*nk
            e=(index-1)%ne+1;m=(index-1)÷ne+1
            @inbounds for b in 1:nb,a in 1:nb
                B[a,b]=retarded ? im*(X[e,m,a,b]-conj(X[e,m,b,a])) : X[e,m,a,b]
            end
            if !all(isfinite,B)
                bad+=1;where_bad=[e,m];continue
            end
            scale=opnorm(B);budget=N.psd_error_budget(scale,nb)
            @inbounds for b in 1:nb,a in 1:nb
                H[a,b]=(B[a,b]+conj(B[b,a]))/2;D[a,b]=B[a,b]-conj(B[b,a])
            end
            herm=norm(D)/2;lowhere=eigmin(Hermitian(H))
            if lowhere<low;low=lowhere;where_min=[e,m];minscale=scale;minbudget=budget;end
            maxherm=max(maxherm,herm)
            excess=max(herm-budget,-lowhere-budget,0.0)
            if excess>0;bad+=1;if excess>worst;worst=excess;where_bad=[e,m];end;end
        end
        reports[worker]=Dict{String,Any}("minimum_eigenvalue"=>low,"minimum_location"=>where_min,
            "minimum_block_norm"=>minscale,"minimum_budget"=>minbudget,"maximum_hermiticity_defect"=>maxherm,
            "failed_blocks"=>bad,"maximum_budget_excess"=>worst,"failure_location"=>where_bad)
    end
    j=argmin([r["minimum_eigenvalue"] for r in reports]);answer=copy(reports[j])
    k=argmax([r["maximum_budget_excess"] for r in reports])
    answer["maximum_budget_excess"]=reports[k]["maximum_budget_excess"]
    answer["failure_location"]=reports[k]["failure_location"]
    answer["failed_blocks"]=sum(r["failed_blocks"] for r in reports)
    answer["maximum_hermiticity_defect"]=maximum(r["maximum_hermiticity_defect"] for r in reports)
    answer["workers"]=reports;answer["checked_blocks"]=ne*nk
    answer["passed"]=answer["failed_blocks"]==0
    answer["budget"]="Coree171 psd_error_budget(opnorm(block),Nb), default construction_scale=1; no clip/floor"
    return answer
end
function family_list(s)
    entries=sort!(collect(s.scattering);by=x->string(first(x)))
    append!(entries,[:embedding=>s.embedding,:embedding_plus=>s.embedding_plus,:embedding_minus=>s.embedding_minus])
    return entries
end
function seed_pair!(L::Array{ComplexF64,4},G::Array{ComplexF64,4},X::Array{ComplexF64,4},f,h;retarded=false,weights=nothing)
    delta=0.0;scale=0.0;weighted_delta=0.0;weighted_scale=0.0;identity=0.0
    @inbounds for b in axes(X,4),a in axes(X,3),m in axes(X,2),e in axes(X,1)
        gamma=retarded ? im*(X[e,m,a,b]-conj(X[e,m,b,a])) : X[e,m,a,b]
        newL=im*f[e]*gamma;newG=-im*h[e]*gamma
        d=abs2(newL-L[e,m,a,b])+abs2(newG-G[e,m,a,b]);s=abs2(L[e,m,a,b])+abs2(G[e,m,a,b])
        delta+=d;scale+=s
        if weights!==nothing;w=weights[1][e]*weights[2][m];weighted_delta+=w*d;weighted_scale+=w*s;end
        identity=max(identity,abs(im*(newG-newL)-gamma))
        L[e,m,a,b]=newL;G[e,m,a,b]=newG
    end
    return Dict{String,Any}("change_absolute"=>sqrt(delta),"old_scale"=>sqrt(scale),"change_relative"=>GP.relative(sqrt(delta),sqrt(scale)),
        "weighted_change_absolute"=>sqrt(weighted_delta),"weighted_old_scale"=>sqrt(weighted_scale),"pair_identity_matrix_max"=>identity,
        "new_lesser_sha256"=>array_sha(L),"new_greater_sha256"=>array_sha(G))
end
function embedding_additivity(s)
    ans=Dict{String,Any}()
    for component in (:Σᴿ,:Σˡ,:Σᵍ)
        T=getfield(s.embedding,component);P=getfield(s.embedding_plus,component);M=getfield(s.embedding_minus,component)
        d=0.0;scale=0.0
        @inbounds for i in eachindex(T);d+=abs2(T[i]-P[i]-M[i]);scale+=abs2(T[i]);end
        ans[string(component)]=Dict("absolute"=>sqrt(d),"scale"=>sqrt(scale),"relative"=>GP.relative(sqrt(d),sqrt(scale)))
    end
    return ans
end
function one_time_fd!(ctx,state,out,c)
    p=ctx.problem;g=p.grids;A=state.green.A;entries=family_list(state)
    ledger=Dict{String,Any}("mode"=>"one-time global FD INITIAL GUESS at finite bias; never projection of later maps",
        "PSD_checks"=>Dict{String,Any}(),"retarded_unchanged"=>Dict{String,Any}(),"changes"=>Dict{String,Any}(),
        "constructed_first_lambda_current"=>"lambda near1/current near0 are initialization identities, NOT physical acceptance",
        "embedding_additivity_before"=>embedding_additivity(state))
    function persist();write_toml(joinpath(out,"FD-initialization.toml"),ledger);end
    # Complete all source checks BEFORE changing any L/G field.
    checks=ledger["PSD_checks"];checks["A"]=seed_cone(A);persist()
    checks["A"]["passed"] || error("FDseed retained A failed native PSD/Hermitian budget; no repair")
    for (name,family) in entries
        stop_due(c) && error("FD precheck stopped by overall/UTC save reserve before mutation")
        checks[string(name)]=seed_cone(family.Σᴿ;retarded=true);persist()
        checks[string(name)]["passed"] || error("FDseed retained GammaR failed native PSD/Hermitian budget: $name; no repair")
    end
    q=zeros(Float64,length(g.ε));before=0.0
    @inbounds for m in axes(A,2),e in axes(A,1),a in axes(A,3)
        w=p.physical.g_s*g.wᴱ[e]*g.wᵏ[m]/(2π)
        q[e]+=w*real(A[e,m,a,a]);before+=w*real(-im*state.green.Gˡ[e,m,a,a])
    end
    stop_due(c) && error("FD precheck stopped before number fit/mutation")
    fit=GP.fit_mu(g.ε,q,ctx.target,ctx.kBT)
    f=zeros(length(q));hole=similar(f)
    for e in eachindex(q);f[e],hole[e]=GP.fermi_pair((g.ε[e]-fit.mu)/ctx.kBT);end
    ledger["number_fit"]=plain(fit);ledger["mu_eV"]=fit.mu*p.scales.E₀_eV;ledger["kBT_scaled"]=ctx.kBT
    ledger["charge_before"]=before;ledger["q"]=q;ledger["f"]=f;ledger["hole"]=hole
    ledger["A_sha256"]=array_sha(A);ledger["GR_sha256"]=array_sha(state.green.Gᴿ)
    for (name,family) in entries
        oldR=array_sha(family.Σᴿ);oldL=array_sha(family.Σˡ);oldG=array_sha(family.Σᵍ)
        change=seed_pair!(family.Σˡ,family.Σᵍ,family.Σᴿ,f,hole;retarded=true,weights=(g.wᴱ,g.wᵏ))
        change["old_lesser_sha256"]=oldL;change["old_greater_sha256"]=oldG
        ledger["changes"][string(name)]=change
        ledger["retarded_unchanged"][string(name)]=array_sha(family.Σᴿ)==oldR
    end
    ledger["changes"]["green"]=seed_pair!(state.green.Gˡ,state.green.Gᵍ,A,f,hole;weights=(g.wᴱ,g.wᵏ))
    ledger["embedding_additivity_after"]=embedding_additivity(state)
    ledger["A_unchanged"]=array_sha(A)==ledger["A_sha256"];ledger["GR_unchanged"]=array_sha(state.green.Gᴿ)==ledger["GR_sha256"]
    all(values(ledger["retarded_unchanged"])) && ledger["A_unchanged"] && ledger["GR_unchanged"] || error("FDseed changed retained retarded arrays")
    ledger["completed"]=true;persist();return state
end
function fd_seed_tiny_checks()
    @testset "one-time FD source cone and exact complex pairing" begin
        X=reshape(ComplexF64[2,im,-im,3],1,1,2,2);L=zeros(ComplexF64,size(X));G=zeros(ComplexF64,size(X))
        @test seed_cone(X)["passed"]
        R=-0.5im.*X;R0=copy(R)
        for voltage in (0.0,1.0)
            f,h=GP.fermi_pair(-2.0);seed_pair!(L,G,R,[f],[h];retarded=true)
            @test R==R0
            @test -im.*L ≈ f.*X && im.*G ≈ h.*X
            @test im.*(G-L) ≈ X
        end
        Y=copy(X);Y[1,1,1,1]=-1.0;@test !seed_cone(Y)["passed"]
        Y=copy(X);Y[1,1,1,2]+=0.01;@test !seed_cone(Y)["passed"]
        Y=copy(X);Y[1,1,1,1]=NaN;@test !seed_cone(Y)["passed"]
        @test !seed_cone(-R;retarded=true)["passed"]
        # Real two-block/worker scan: min eigenvalue and largest Hermiticity excess differ.
        Z=zeros(ComplexF64,2,1,2,2)
        Z[1,1,1,1]=-2;Z[1,1,2,2]=-2
        Z[2,1,1,1]=1;Z[2,1,2,2]=1;Z[2,1,1,2]=10;Z[2,1,2,1]=-10
        witness=seed_cone(Z)
        @test length(witness["workers"])==2 && witness["failed_blocks"]==2
        @test witness["minimum_eigenvalue"]==-2 && witness["minimum_location"]==[1,1]
        @test witness["maximum_budget_excess"]>10 && witness["failure_location"]==[2,1]
        @test_throws GP.Infeasible GP.fit_mu([0.0,1.0],[0.0,0.0],1.0,0.1)
        f,h=GP.fermi_pair(-100.0);@test h>0 && f==1.0
    end
end
