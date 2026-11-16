# Selected intended-model qualifier only. No Core/Runner/Pkg/full-grid calls.
using LinearAlgebra, HDF5, TOML, SHA, Test, Dates
const RQ=Rational{BigInt}
const CQ=Complex{RQ}
digest(path)=bytes2hex(open(sha256,path))
cx(f,name)=complex.(read(f[name*"_real"]),read(f[name*"_imag"]))
function putcx(f,name,X)
    f[name*"_real"]=real.(ComplexF64.(X));f[name*"_imag"]=imag.(ComplexF64.(X))
end
function checked(path,expected)
    occursin(r"^[0-9a-f]{64}$",expected) || error("source SHA must be frozen")
    actual=digest(path);actual==expected || error("source SHA mismatch: $path")
    return actual
end
rel(a,b)=b==0 ? (a==0 ? zero(a) : oftype(a,Inf)) : a/b
qcomplex(z::ComplexF64)=CQ(RQ(real(z)),RQ(imag(z))) # Base exact, tol=0 conversion
function rational_record(x::RQ)
    decimal=setprecision(BigFloat,256) do
        string(BigFloat(numerator(x))/BigFloat(denominator(x)))
    end
    Dict{String,Any}("numerator"=>string(numerator(x)),"denominator"=>string(denominator(x)),
        "sign"=>(x<0 ? -1 : x>0 ? 1 : 0),"decimal256"=>decimal)
end
"Exact no-pivot Hermitian LDL; only positive pivots plus exact reconstruction certify PD."
function exact_ldl(X::Matrix{ComplexF64})
    all(isfinite,X) || error("nonfinite source Gamma")
    Q=qcomplex.(X);n=size(Q,1);size(Q,2)==n || error("nonsquare Gamma")
    rec=Dict{String,Any}("matrix_dimension"=>n,"exact_Hermitian"=>(Q==Q'),"strict_PD_certified"=>false,
        "reconstruction_exact"=>false,"pivots"=>Dict{String,Any}[],
        "scope"=>"exact binary Float64 selected Gamma only; no full-grid or input-uncertainty certificate")
    if Q!=Q'
        rec["reason"]="Gamma is not exactly Hermitian; no repair or PD claim"
        return rec
    end
    L=zeros(CQ,n,n);d=zeros(RQ,n)
    for j in 1:n
        L[j,j]=one(RQ);p=Q[j,j]
        for k in 1:j-1;p-=L[j,k]*d[k]*conj(L[j,k]);end
        iszero(imag(p)) || error("internal LDL imaginary pivot")
        d[j]=real(p);push!(rec["pivots"],rational_record(d[j]))
        if d[j]<=0
            rec["reason"]="nonpositive exact leading LDL pivot; strict PD not certified"
            return rec
        end
        for i in j+1:n
            z=Q[i,j]
            for k in 1:j-1;z-=L[i,k]*d[k]*conj(L[j,k]);end
            L[i,j]=z/d[j]
        end
    end
    reconstructed=zeros(CQ,n,n)
    for b in 1:n,a in 1:n,k in 1:n
        reconstructed[a,b]+=L[a,k]*d[k]*conj(L[b,k])
    end
    rec["reconstruction_exact"]=reconstructed==Q
    rec["reconstruction_exact"] || error("internal exact LDL reconstruction failed")
    rec["strict_PD_certified"]=true;rec["reason"]="positive rational pivots and exact LDLdagger reconstruction"
    return rec
end
function distance(X,Y)
    numerator=norm(X-Y);denominator=norm(Y);ratio=rel(numerator,denominator)
    Dict{String,Any}("absolute"=>Float64(numerator),"scale"=>Float64(denominator),"relative"=>Float64(ratio),
        "absolute_decimal"=>string(numerator),"scale_decimal"=>string(denominator),"relative_decimal"=>string(ratio))
end
"All nine cells use the same h_H; split primitives promoted BEFORE any assembly."
function intended_oracle(energy,hH,LH,Gamma,Tplus;bits=192,periods=4,native=nothing)
    setprecision(BigFloat,bits) do
        CT=Complex{BigFloat};H,L,K,T=CT.(hH),CT.(LH),CT.(Gamma),CT.(Tplus)
        n=size(H,1);cells=2periods+1;N=cells*n;D=zeros(CT,N,N)
        block=-H-L+(im/BigFloat(2))*K
        for a in 1:n;block[a,a]+=BigFloat(energy);end
        for cell in 1:cells
            r=(cell-1)*n+1:cell*n;D[r,r]=block
            if cell<cells
                s=cell*n+1:(cell+1)*n;D[r,s]=-T;D[s,r]=-T'
            end
        end
        center=periods*n+1:(periods+1)*n;rhs=zeros(CT,N,n)
        for a in 1:n;rhs[center[a],a]=1;end
        # Adjoint solve gives CENTRAL ROW, needed for repeated-source covariance.
        Y=lu(Matrix(D'))\rhs;row=Matrix(Y');G=Matrix(row[:,center]);A=zeros(CT,n,n)
        for cell in 1:cells
            B=Matrix(row[:,(cell-1)*n+1:cell*n]);A+=B*K*B'
        end
        Asub=im*(G-G');residual=norm(D'*Y-rhs);rescale=norm(D)*norm(Y)+norm(rhs)
        traces=Dict{String,Any}("A_product_trace"=>string(real(tr(A))),"A_sub_trace"=>string(real(tr(Asub))),
            "G_real_trace"=>string(real(tr(G))),"G_imag_trace"=>string(imag(tr(G))),
            "A_antiHermitian_F"=>string(norm(A-A')),"LU_backward_absolute"=>string(residual),
            "LU_backward_scale"=>string(rescale),"LU_backward_relative"=>string(rel(residual,rescale)))
        metrics=Dict{String,Any}("Big_sub_vs_product"=>distance(Asub,A),
            "condition2_F64_downcast_D"=>cond(ComplexF64.(D)),
            "condition_scope"=>"Float64 SVD condition estimate of downcast intended D; not a certified Big condition bound")
        if native!==nothing
            metrics["Double_A_cov_vs_Big_product"]=distance(CT.(native.A),A)
            metrics["Double_Asub_vs_Big_product"]=distance(CT.(native.Asub),A)
            metrics["Double_GR_vs_Big_GR"]=distance(CT.(native.GR),G)
            metrics["comparison_scope"]="includes split-vs-F64-assembled Sigma/cavity/GR/A arithmetic; does not isolate LU alone"
        end
        return (G=G,A=A,Asub=Asub,row=row,traces=traces,metrics=metrics)
    end
end
function tiny_checks()
    @testset "intended split positive-source oracle" begin
        @test RQ(0.1)==BigInt(3602879701896397)//BigInt(36028797018963968)
        K=ComplexF64[2 0.5+im;0.5-im 3];proof=exact_ldl(K)
        @test proof["strict_PD_certified"] && proof["reconstruction_exact"]
        @test proof["pivots"][2]["numerator"]=="19" && proof["pivots"][2]["denominator"]=="8"
        @test !exact_ldl(ComplexF64[1 0;0 -1])["strict_PD_certified"]
        bad=copy(K);bad[1,2]+=1e-12
        @test !exact_ldl(bad)["exact_Hermitian"]
        h=ComplexF64[0.1 0.03+0.02im;0.03-0.02im 0.4]
        L=ComplexF64[0.02 0.01im;-0.01im -0.03];T=ComplexF64[0.01 0.03im;0.02+0.01im 0.015]
        Gamma=0.01 .* K;energy=1.3;D=energy*I-h-L+0.5im*Gamma
        side=inv(D);Gc=inv(D-T*side*T'-T'*side*T)
        source=Gamma+T*side*Gamma*side'*T'+T'*side*Gamma*side'*T
        analyticA=Gc*source*Gc';o=intended_oracle(energy,h,L,Gamma,T;periods=1)
        @test norm(ComplexF64.(o.G)-Gc)/norm(Gc)<2e-13
        @test norm(ComplexF64.(o.A)-analyticA)/norm(analyticA)<2e-13
        @test o.metrics["Big_sub_vs_product"]["relative"]<1e-40
        @test maximum(abs,imag.(analyticA))>1e-6 # nontrivial complex offdiagonal
        weak=1e-25 .* Gamma;w=intended_oracle(energy,h,L,weak,T;periods=1)
        @test exact_ldl(weak)["strict_PD_certified"]
        @test w.metrics["Big_sub_vs_product"]["relative"]<1e-20
        @test real(tr(w.A))>0
        @test TOML.parse(sprint(io->TOML.print(io,proof)))["pivots"][2]["numerator"]=="19"
    end
end
function run_oracle(input,output)
    cfg=TOML.parsefile(input);cfg["runtime_authorized"]===true || error("write-only freeze; runtime disabled")
    string(VERSION)==cfg["expected_julia_version"] && Threads.nthreads()==1 || error("Julia/thread identity differs")
    cfg["precision_bits"]==192 && cfg["repeat_bits"]==384 || error("precision identity differs")
    base=dirname(abspath(input));path(key)=normpath(joinpath(base,cfg[key]))
    for (p,s) in (("selected_path","selected_sha256"),("source_input_path","source_input_sha256"),
                  ("source_summary_path","source_summary_sha256"),("source_script_path","source_script_sha256"))
        checked(path(p),cfg[s])
    end
    prior=TOML.parsefile(path("source_summary_path"))
    prior["execution"]=="complete" && prior["maps_started"]==prior["maps_completed"]==3 || error("source preflight incomplete")
    prior["source_sha256"]==cfg["source_script_sha256"] && prior["input_sha256"]==cfg["source_input_sha256"] || error("source summary bindings differ")
    ispath(output) && error("refusing overwrite");mkpath(output);cp(input,joinpath(output,"input.toml"));cp(@__FILE__,joinpath(output,"oracle.jl"))
    BLAS.set_num_threads(1);started=time();checkwall()=time()-started<110 || error("internal110s wall stop")
    report=Dict{String,Any}("execution"=>"running","script_sha256"=>digest(@__FILE__),"input_sha256"=>digest(input),
        "selected_sha256"=>cfg["selected_sha256"],"scope"=>"selected intended split model, all9 same h_H; no full-grid or SCBA acceptance",
        "scientific_acceptance"=>false,
        "nodes"=>Dict{String,Any}(),"selected_q"=>Dict{String,Any}(),"precision_bits"=>192,
        "precision_scope"=>"192 arithmetic, ONE384 repeat at worst selected DoubleAcov relative-error node; no interval or all-node384 certificate")
    try
        tiny_checks();report["tiny_checks"]="pass"
        HDF5.h5open(path("selected_path"),"r") do f
            read(f["schema"])=="campaign133-positive-source-operator-selected-v1" || error("selected schema differs")
            indices=Int.(read(f["selected_energy_indices"]));indices==cfg["expected_selected_energy_indices"] || error("selected E scope differs")
            energy,we,wk=read(f["energy"]),read(f["weight_energy"]),read(f["weight_k"])
            gs=read(f["spin_degeneracy"]);E0=read(f["energy_scale_eV"]);target=read(f["target"])
            length(energy)==length(we)==30704 && length(wk)==11 || error("grid metadata shape differs")
            all(isfinite,energy) && all(w->isfinite(w)&&w>0,we) && all(w->isfinite(w)&&w>0,wk) || error("nonfinite/nonpositive quadrature metadata")
            K,L,G,A,S=(cx(f,key) for key in ("Gamma","Lambda_H_preassembly","GR","A","Asub"))
            H,T=cx(f,"h_H"),cx(f,"Tplus");ns,nk,nb,nb2=size(K)
            (ns,nk,nb,nb2)==(4,11,5,5) && ns*nk<=77 || error("expected selected44 logical(E,k,a,b) blocks")
            size(L)==size(G)==size(A)==size(S)==size(K) || error("selected field axes differ")
            size(H)==(11,5,5) && size(T)==(5,5) || error("h/T axes differ")
            all(x->all(isfinite,x),(K,L,G,A,S,H,T)) || error("nonfinite selected inputs")
            qcov,qsub=read(f["q_cov"]),read(f["q_sub"]);acc=Dict{Int,Vector{BigFloat}}();counts=Dict{Int,Int}()
            worst=nothing;worst_error=-Inf;worst_sub=nothing;worst_suberror=-Inf
            traces=String["energy_index,k_index,A_product_trace_192,A_sub_trace_192,Double_A_cov_trace,Double_A_sub_trace"]
            for j in 1:ns,m in 1:nk
                checkwall();e=indices[j];name="e$(e)_k$(m)"
                h=Matrix(H[m,:,:]);lambda=Matrix(L[j,m,:,:]);gamma=Matrix(K[j,m,:,:])
                native=(GR=Matrix(G[j,m,:,:]),A=Matrix(A[j,m,:,:]),Asub=Matrix(S[j,m,:,:]))
                proof=exact_ldl(gamma)
                node=Dict{String,Any}("energy_index"=>e,"k_index"=>m,"energy_eV"=>energy[e]*E0,"Gamma_exact_LDL"=>proof,
                    "h_H_exact_Hermitian"=>(qcomplex.(h)==qcomplex.(h)'),"Lambda_pre_exact_Hermitian"=>(qcomplex.(lambda)==qcomplex.(lambda)'))
                report["nodes"][name]=node
                proof["exact_Hermitian"] && node["h_H_exact_Hermitian"] && node["Lambda_pre_exact_Hermitian"] || error("intended Hermitian-input contract failed; no repair")
                o=intended_oracle(energy[e],h,lambda,gamma,T;native)
                node["traces192"]=o.traces;node["metrics192"]=o.metrics
                eigenvalues=eigvals(Hermitian(gamma));node["Gamma_Float_eigenvalues"]=eigenvalues
                node["Gamma_relative_eigensolve_scale"]=64*nb*eps(Float64)*norm(gamma)
                node["Gamma_eigenvalue_scope"]="Float comparison only; exact positive LDL is the selected PD proof"
                setprecision(BigFloat,192) do
                    sums=get!(acc,e,BigFloat[0,0]);weight=BigFloat(Float64(gs)/(2pi))*BigFloat(we[e])*BigFloat(wk[m])
                    sums[1]+=weight*parse(BigFloat,o.traces["A_product_trace"]);sums[2]+=weight*parse(BigFloat,o.traces["A_sub_trace"])
                end
                counts[e]=get(counts,e,0)+1
                push!(traces,join((e,m,o.traces["A_product_trace"],o.traces["A_sub_trace"],real(tr(native.A)),real(tr(native.Asub))),","))
                cov_error=o.metrics["Double_A_cov_vs_Big_product"]["relative"];suberror=o.metrics["Double_Asub_vs_Big_product"]["relative"]
                saved=(name=name,e=e,m=m,h=h,L=lambda,K=gamma,T=T,native=native,result=o)
                if worst===nothing || cov_error>worst_error;worst=saved;worst_error=cov_error;end
                if worst_sub===nothing || suberror>worst_suberror;worst_sub=saved;worst_suberror=suberror;end
            end
            checkwall();high=intended_oracle(energy[worst.e],worst.h,worst.L,worst.K,worst.T;bits=384,native=worst.native)
            repeat=setprecision(BigFloat,384) do
                Dict{String,Any}("A_product_trace"=>distance(parse(BigFloat,worst.result.traces["A_product_trace"]),parse(BigFloat,high.traces["A_product_trace"])),
                    "A_sub_trace"=>distance(parse(BigFloat,worst.result.traces["A_sub_trace"]),parse(BigFloat,high.traces["A_sub_trace"])),
                    "traces384"=>high.traces,"metrics384"=>high.metrics)
            end
            report["repeat384"]=Dict{String,Any}("node"=>worst.name,"selection_rule"=>"maximum selected DoubleAcov/Big intended-A relative error","result"=>repeat)
            for e in indices
                report["selected_q"]["e$(e)"]=Dict{String,Any}("k_count"=>counts[e],"q192_product"=>string(acc[e][1]),"q192_subtraction"=>string(acc[e][2]),
                    "stored_Double_q_cov"=>qcov[e],"stored_Double_q_sub"=>qsub[e],"q192_product_over_target"=>Float64(acc[e][1])/target,
                    "measure"=>"all11k;192 sums, exact promoted Double gs/(2pi),wE,wk; not a whole energy-grid density integral")
            end
            HDF5.h5open(joinpath(output,"critical.h5"),"w") do out
                out["schema"]="campaign133-positive-source-selected-big-critical-v1"
                seen=Set{String}()
                for block in (worst,worst_sub)
                    block.name in seen && continue;push!(seen,block.name);group=HDF5.create_group(out,block.name)
                    group["energy_index"]=block.e;group["k_index"]=block.m
                    for (key,X) in (("h_H",block.h),("Lambda_H_preassembly",block.L),("Gamma",block.K),("Tplus",block.T),
                                    ("Double_GR",block.native.GR),("Double_A_cov",block.native.A),("Double_A_sub",block.native.Asub),
                                    ("Big192_GR_downcast",block.result.G),("Big192_A_product_downcast",block.result.A),
                                    ("Big192_A_sub_downcast",block.result.Asub),("Big192_central_row_downcast",block.result.row))
                        putcx(group,key,X)
                    end
                end
            end
            write(joinpath(output,"traces.csv"),join(traces,"\n")*"\n")
            report["exact_PD_certified_count"]=count(n->n["Gamma_exact_LDL"]["strict_PD_certified"],Base.values(report["nodes"]))
            report["selected_block_count"]=ns*nk;report["all_selected_exact_PD"]=report["exact_PD_certified_count"]==ns*nk
            report["worst_cov_relative"]=worst_error;report["worst_sub_relative"]=worst_suberror
            report["critical_nodes"]=Dict("worst_cov"=>worst.name,"worst_sub"=>worst_sub.name)
            report["matrix_storage_scope"]="critical HDF matrices are Float64 downcasts;192/384 traces and residuals remain decimal strings"
            report["execution"]="complete"
        end
    catch error
        report["execution"]="stopped";report["reason"]=sprint(showerror,error,catch_backtrace());rethrow()
    finally
        report["elapsed_seconds"]=time()-started
        open(joinpath(output,"summary.toml"),"w") do io;TOML.print(io,report);end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    ARGS==["--tiny"] ? tiny_checks() : (length(ARGS)==2 ? run_oracle(abspath(ARGS[1]),abspath(ARGS[2])) : error("--tiny OR INPUT NEW_OUTPUT"))
end
