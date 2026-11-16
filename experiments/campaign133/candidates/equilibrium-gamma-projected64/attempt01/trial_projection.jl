# Algorithmic AA proposal only. No residual/map/channel ever calls this operator.
"Worker-owned LAPACK zheev eigenvectors; same ABI as accepted Core eigenvalues helper."
function projection_eigenvectors!(H,values,work,rwork,info)
    n=LinearAlgebra.BlasInt(size(H,1));lwork=LinearAlgebra.BlasInt(length(work));info[1]=0
    ccall((LinearAlgebra.BLAS.@blasfunc(zheev_),LinearAlgebra.LAPACK.libblastrampoline),Cvoid,
        (Ref{UInt8},Ref{UInt8},Ref{LinearAlgebra.BlasInt},Ptr{ComplexF64},Ref{LinearAlgebra.BlasInt},
         Ptr{Float64},Ptr{ComplexF64},Ref{LinearAlgebra.BlasInt},Ptr{Float64},Ptr{LinearAlgebra.BlasInt},Clong,Clong),
        'V','U',n,H,n,values,work,lwork,rwork,info,1,1)
    info[1]==0 || throw(Infeasible("AA trial projection LAPACK info=$(info[1])"))
end
function projection_worker!(Y::Field,X::Field,we::Vector{Float64},wk::Vector{Float64},first,last)
    ne,nk,nb,_=size(X)
    raw=zeros(ComplexF64,nb,nb);U=similar(raw);P=similar(raw);H=similar(raw)
    values=zeros(Float64,nb);work=zeros(ComplexF64,max(1,3nb));rwork=zeros(Float64,max(1,5nb));info=zeros(LinearAlgebra.BlasInt,1)
    d2=0.0;dw2=0.0;s2=0.0;sw2=0.0;expected2=0.0;herm2=0.0;tracechange=0.0;weightedtrace=0.0;tracebefore=0.0;traceafter=0.0
    changed=0;unchanged=0;negative=0;zero=0;ranks=zeros(Int,nb+1)
    minlambda=Inf;minwhere=(0,0);maxdelta=0.0;deltawhere=(0,0)
    maxeigres=0.0;eigwhere=(0,0);maxorth=0.0;orthwhere=(0,0)
    maxeigenabs=0.0;maxeigscale=0.0
    @inbounds for linear in first:last
        e=(linear-1)%ne+1;k=(linear-1)÷ne+1;exact=true;scale2=0.0;anti2=0.0
        for b in 1:nb,a in 1:nb
            z=X[e,k,a,b];other=conj(X[e,k,b,a]);raw[a,b]=z
            H[a,b]=(z+other)/2;U[a,b]=H[a,b];exact &= z==other
            scale2+=abs2(z);anti2+=abs2(z-other)
        end
        projection_eigenvectors!(U,values,work,rwork,info)
        all(isfinite,values) && all(isfinite,U) || throw(Infeasible("nonfinite AA projection eigensystem"))
        low=values[1]
        if low<minlambda;minlambda=low;minwhere=(e,k);end
        rank=count(>(0.0),values);ranks[rank+1]+=1
        negative+=count(<(0.0),values);zero+=count(==(0.0),values)
        neg2=0.0
        for j in 1:nb;values[j]<0 && (neg2+=abs2(values[j]));end
        expected2+=neg2+anti2/4;herm2+=anti2
        eig2=0.0;orth2=0.0
        for a in 1:nb,b in 1:nb
            eigen=0.0+0.0im;orth=0.0+0.0im
            for j in 1:nb;eigen+=H[a,j]*U[j,b];orth+=conj(U[j,a])*U[j,b];end
            eig2+=abs2(eigen-U[a,b]*values[b]);orth2+=abs2(orth-(a==b ? 1.0 : 0.0))
        end
        eigabs=sqrt(eig2);eigscale=sqrt(scale2);eigrel=relative(eigabs,eigscale);orth=sqrt(orth2)
        if eigrel>maxeigres;maxeigres=eigrel;eigwhere=(e,k);maxeigenabs=eigabs;maxeigscale=eigscale;end
        if orth>maxorth;maxorth=orth;orthwhere=(e,k);end
        # Exact-Hermitian feasible block is retained bit for bit, without reconstruction.
        if exact && low>=0
            copyto!(P,raw);unchanged+=1
        else
            for b in 1:nb,a in 1:nb
                z=0.0+0.0im
                for j in 1:nb;z+=U[a,j]*max(values[j],0.0)*conj(U[b,j]);end
                P[a,b]=z
            end
            # Explicit Hermitian representation of the reconstruction; no spectral floor.
            for b in 1:nb,a in 1:b
                z=(P[a,b]+conj(P[b,a]))/2;P[a,b]=z;P[b,a]=conj(z)
            end
            changed+=1
        end
        blockdelta=0.0;trace=0.0
        for b in 1:nb,a in 1:nb
            Y[e,k,a,b]=P[a,b];blockdelta+=abs2(P[a,b]-raw[a,b])
            if a==b
                trace+=real(P[a,b]-raw[a,b]);tracebefore+=real(raw[a,b]);traceafter+=real(P[a,b])
            end
        end
        weight=we[e]*wk[k];d2+=blockdelta;dw2+=weight*blockdelta;s2+=scale2;sw2+=weight*scale2
        tracechange+=trace;weightedtrace+=weight*trace
        if sqrt(blockdelta)>maxdelta;maxdelta=sqrt(blockdelta);deltawhere=(e,k);end
    end
    return (;d2,dw2,s2,sw2,expected2,herm2,tracechange,weightedtrace,tracebefore,traceafter,changed,unchanged,negative,zero,ranks,
        minlambda,minwhere,maxdelta,deltawhere,maxeigres,eigwhere,maxorth,orthwhere,maxeigenabs,maxeigscale)
end
"Blockwise nearest PSD(Herm X); X is never modified; no positive eigenvalue floor."
function project_psd_trial(X::Field,we::Vector{Float64},wk::Vector{Float64})
    all(isfinite,X) || throw(Infeasible("nonfinite AA trial before projection"))
    Y=similar(X);workers=min(4,Threads.nthreads(),size(X,1)*size(X,2));records=Vector{Any}(undef,workers)
    N._production_worker_ranges!(size(X,1)*size(X,2),workers) do worker,first,last
        records[worker]=projection_worker!(Y,X,we,wk,first,last)
    end
    sumfield(key)=sum(r->getproperty(r,key),records)
    low=argmin([r.minlambda for r in records]);delta=argmax([r.maxdelta for r in records])
    eig=argmax([r.maxeigres for r in records]);orth=argmax([r.maxorth for r in records])
    changed=sumfield(:changed);absolute=sqrt(sumfield(:d2));scale=sqrt(sumfield(:s2))
    weighted=sqrt(sumfield(:dw2));weightedscale=sqrt(sumfield(:sw2))
    ledger=Dict{String,Any}("operator"=>"P_PSD(Herm AAtrial), eigenvalue floor exactly zero",
        "changed_blocks"=>changed,"unchanged_exact_feasible_blocks"=>sumfield(:unchanged),
        "negative_input_eigenvalues"=>sumfield(:negative),"zero_input_eigenvalues"=>sumfield(:zero),
        "output_positive_rank_counts"=>sum((r.ranks for r in records)),
        "absolute"=>absolute,"scale"=>scale,"relative"=>relative(absolute,scale),
        "weighted_absolute"=>weighted,"weighted_scale"=>weightedscale,"weighted_relative"=>relative(weighted,weightedscale),
        "expected_orthogonal_distance"=>sqrt(sumfield(:expected2)),"input_antihermitian_norm"=>sqrt(sumfield(:herm2)),
        "trace_before"=>sumfield(:tracebefore),"trace_after"=>sumfield(:traceafter),
        "trace_change"=>sumfield(:tracechange),"weighted_trace_change"=>sumfield(:weightedtrace),
        "distance_formula_absolute_defect"=>abs(absolute-sqrt(sumfield(:expected2))),
        "distance_formula_relative_defect"=>relative(abs(absolute-sqrt(sumfield(:expected2))),sqrt(sumfield(:expected2))),
        "minimum_input_eigenvalue"=>records[low].minlambda,"minimum_input_location"=>collect(records[low].minwhere),
        "maximum_block_change"=>records[delta].maxdelta,"maximum_change_location"=>collect(records[delta].deltawhere),
        "maximum_eigen_backward_relative"=>records[eig].maxeigres,"eigen_backward_location"=>collect(records[eig].eigwhere),
        "eigen_backward_absolute_at_maximum"=>records[eig].maxeigenabs,"eigen_backward_scale_at_maximum"=>records[eig].maxeigscale,
        "maximum_eigenvector_orthogonality_defect"=>records[orth].maxorth,"orthogonality_location"=>collect(records[orth].orthwhere),
        "arithmetic_scope"=>"eigen/reconstruction rounding is measured, not a strict computed-PSD certificate; actual trial is guarded by original cone and fresh Pi domain",
        "unchanged_field_identity"=>changed==0)
    points=unique(filter(p->p!=(0,0),[records[low].minwhere,records[delta].deltawhere,records[eig].eigwhere,records[orth].orthwhere]))
    return (changed==0 ? X : Y),ledger,points
end
function project_aa_trial(X::Field,we,wk,output,cfg,nextcall)
    started=time();Y,ledger,points=project_psd_trial(X,we,wk)
    index=count(name->startswith(name,"projection-") && endswith(name,".toml"),readdir(output))+1
    ledger["projection_attempt"]=index;ledger["next_pi_call"]=nextcall;ledger["seconds"]=time()-started
    ledger["full_pair_preserved"]=index<=4 # declared before the experiment
    path=joinpath(output,"projection-$(lpad(index,3,'0')).h5")
    expected=index<=4 ? 2*expected_padded_complex_bytes(size(X))+8*1024^2 : 2*length(points)*size(X,3)^2*16+1024^2
    atomic_h5(path,cfg;expected_bytes=expected) do f
        f["schema"]="campaign133-AAtrial-PSD-projection-v1";f["projection_attempt"]=index;f["next_pi_call"]=nextcall
        f["ledger_toml"]=sprint(io->TOML.print(io,ledger))
        if index<=4
            write_complex(f,"raw_trial",X);write_complex(f,"projected_trial",Y)
        else
            raw=zeros(ComplexF64,length(points),size(X,3),size(X,4));projected=similar(raw)
            for (i,(e,k)) in enumerate(points);raw[i,:,:]=X[e,k,:,:];projected[i,:,:]=Y[e,k,:,:];end
            f["selected_energy_indices"]=[p[1] for p in points];f["selected_k_indices"]=[p[2] for p in points]
            write_complex(f,"raw_trial_selected",raw);write_complex(f,"projected_trial_selected",projected)
        end
    end
    open(joinpath(output,"projection-$(lpad(index,3,'0')).toml"),"w") do io;TOML.print(io,ledger);end
    return Y,ledger
end
