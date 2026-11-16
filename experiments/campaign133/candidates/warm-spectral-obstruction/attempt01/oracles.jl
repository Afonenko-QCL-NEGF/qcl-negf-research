using LinearAlgebra, HDF5, TOML, SHA, Test
digest(path)=bytes2hex(open(sha256,path))
cx(f,name)=complex.(read(f[name*"_real"]),read(f[name*"_imag"]))
function putcx(f,name,X)
    f[name*"_real"]=real.(ComplexF64.(X));f[name*"_imag"]=imag.(ComplexF64.(X))
end
rel(n,d)=d==0 ? (n==0 ? 0.0 : Inf) : n/d
function distance(X,Y)
    Dict("absolute"=>Float64(norm(X-Y)),"reference_norm"=>Float64(norm(Y)),
        "relative"=>Float64(rel(norm(X-Y),norm(Y))))
end
function psd_info(X)
    H=(X+X')/2;scale=norm(X);values=eigvals(Hermitian(ComplexF64.(H)))
    Dict("norm_F"=>Float64(scale),"hermiticity_F"=>Float64(norm(X-X')),
        "minimum_eigenvalue_Herm"=>first(values),"eigenvalues_Herm"=>values,
        "relative_eigensolve_scale"=>64*size(X,1)*eps(Float64)*Float64(scale))
end

"Exact-double primitives are converted BEFORE operations at requested precision."
function chain_data(energy,h,Tplus,localR,gamma,CT; periods=4,herm_h=false,herm_lambda=false)
    # Actual finite_chain side cells use the F64 plan.hlocal; central uses raw h.
    # Freeze that F64 construction BEFORE promoting exact doubles to BigFloat.
    side64=(ComplexF64.(h)+ComplexF64.(h)')/2
    H,T,S,K=CT.(h),CT.(Tplus),CT.(localR),CT.(gamma)
    side=CT.(side64);herm_h && (H=side)
    if herm_lambda
        L=S+convert(CT,im)*K/2
        S=(L+L')/2-convert(CT,im)*K/2
    end
    nb=size(H,1);cells=2periods+1;n=cells*nb;D=zeros(CT,n,n)
    effective=Matrix{CT}[]
    for j in 1:cells
        onsite=j==periods+1 ? H : side
        localD=-onsite-S
        for a in 1:nb;localD[a,a]+=convert(CT,energy);end # native assembly order
        rows=(j-1)*nb+1:j*nb;D[rows,rows]=localD
        push!(effective,convert(CT,im)*((onsite+S)-(onsite+S)'))
        if j<cells
            next=j*nb+1:(j+1)*nb
            D[rows,next]=-T;D[next,rows]=-T'
        end
    end
    center=periods*nb+1:(periods+1)*nb
    return D,center,K,effective
end
function central_products(row,K,effective)
    nb=size(K,1);physical=zeros(eltype(row),nb,nb);eff=similar(physical);fill!(eff,0)
    for (cell,firstcol) in enumerate(1:nb:size(row,2))
        block=row[:,firstcol:firstcol+nb-1]
        physical+=block*K*block';eff+=block*effective[cell]*block'
    end
    return physical,eff
end
function float_oracle(energy,h,T,S,K;periods=4,herm_h=false,herm_lambda=false)
    D,c,g,eff=chain_data(energy,h,T,S,K,ComplexF64;periods,herm_h,herm_lambda)
    full=inv(D);row=full[c,:];R=row[:,c];A=im*(R-R')
    physical,effective=central_products(row,g,eff)
    return (R=R,A=A,physical=physical,effective=effective,row=row,
        residual=norm(D*full-I)/sqrt(size(D,1)))
end
function big_oracle(energy,h,T,S,K;periods=4,herm_h=false,herm_lambda=false,bits=192)
    setprecision(BigFloat,bits) do
        D,c,g,eff=chain_data(energy,h,T,S,K,Complex{BigFloat};periods,herm_h,herm_lambda)
        nb=length(c);rhs=zeros(eltype(D),size(D,1),nb)
        for a in 1:nb;rhs[c[a],a]=1;end
        # Five RHS of D† give the CENTRAL ROW of GR. D\P would give a column.
        Y=lu(Matrix(D'))\rhs;row=Matrix(Y');R=row[:,c]
        A=im*(R-R');physical,effective=central_products(row,g,eff)
        residual=norm(D'*Y-rhs)/norm(rhs)
        traces=Dict("A_trace"=>string(real(tr(A))),"physical_trace"=>string(real(tr(physical))),
            "effective_trace"=>string(real(tr(effective))),"LU_residual"=>string(residual),
            "A_minus_effective_F"=>string(norm(A-effective)),"A_minus_physical_F"=>string(norm(A-physical)))
        return (R=ComplexF64.(R),A=ComplexF64.(A),physical=ComplexF64.(physical),
            effective=ComplexF64.(effective),row=ComplexF64.(row),traces=traces)
    end
end
function tiny_checks()
    @testset "45-chain oracle small controls" begin
        h=fill(0.2+0im,1,1);T=fill(0.1+0im,1,1);K=fill(1e-25+0im,1,1);S=-im*K/2
        f=float_oracle(2.0,h,T,S,K;periods=1)
        b=big_oracle(2.0,h,T,S,K;periods=1)
        z=2.0-h[1]-S[1];analytic=inv(z-2*T[1]^2/z)
        @test f.R[1]≈analytic rtol=1e-14
        @test b.A≈b.physical rtol=1e-14
        bad=S+fill(1e-18im,1,1)
        raw=big_oracle(2.0,h,T,bad,K;periods=1)
        projected=big_oracle(2.0,h,T,bad,K;periods=1,herm_lambda=true)
        @test real(tr(raw.A))<0 && real(tr(raw.physical))>0
        @test raw.A≈raw.effective rtol=1e-14
        @test real(tr(projected.A))>0
        hc=ComplexF64[0 1+im;1-im 2];tc=ComplexF64[0.1 0.02im;0.03 0.2]
        gc=0.02*Matrix{ComplexF64}(I,2,2);sc=-im*gc/2
        fc=float_oracle(3.0,hc,tc,sc,gc;periods=1);bc=big_oracle(3.0,hc,tc,sc,gc;periods=1)
        @test fc.R≈bc.R rtol=2e-13
        @test bc.A≈bc.physical rtol=2e-14
        # Protect mixed native contract: ONLY central hraw, two side hHerm.
        hn=copy(hc);hn[1,2]+=1e-5;hs=(hn+hn')/2
        mixed=float_oracle(5.0,hn,tc,sc,gc;periods=1)
        side=inv(5.0I-hs-sc)
        schur=inv(5.0I-hn-sc-tc*side*tc'-tc'*side*tc)
        @test mixed.R≈schur rtol=2e-14
        bm=big_oracle(5.0,hn,tc,sc,gc;periods=1)
        @test bm.A≈bm.effective rtol=2e-14
    end
end

function run_oracles(input,output)
    cfg=TOML.parsefile(input);cfg["runtime_authorized"]===true || error("Diagnostic freeze; runtime disabled")
    base=dirname(abspath(input));source=normpath(joinpath(base,cfg["selected_path"]))
    model=normpath(joinpath(base,cfg["cold_model_path"]))
    digest(source)==cfg["selected_sha256"] || error("selected input SHA differs")
    digest(model)==cfg["cold_model_sha256"] || error("cold model SHA differs")
    ispath(output) && error("Refusing to overwrite diagnostic output")
    mkpath(output);cp(input,joinpath(output,"input.toml"));cp(@__FILE__,joinpath(output,"oracles.jl"))
    started=time();BLAS.set_num_threads(1)
    checkwall()=time()-started<300 || error("oracle wall cap exhausted")
    report=Dict{String,Any}("scope"=>"selected exact-double 45-chain oracles only; no retarded full pass/mu/Phi/SCBA",
        "selected_sha256"=>cfg["selected_sha256"],"cold_model_sha256"=>cfg["cold_model_sha256"],
        "script_sha256"=>digest(@__FILE__),"precision_bits"=>192,"nodes"=>Dict{String,Any}())
    try
        tiny_checks()
        Tminus=HDF5.h5open(model,"r") do f;cx(f,"basis_T₋");end
        HDF5.h5open(source,"r") do f
            read(f["schema"])=="campaign133-warm-q-selected-v1" || error("selected schema differs")
            energy,we,wk=read(f["energy"]),read(f["weight_energy"]),read(f["weight_k"])
            selected=Int.(read(f["selected_energy_indices"]));E0=read(f["energy_scale_eV"]);gs=read(f["spin_degeneracy"])
            q,qe,qp=read(f["q_raw"]),read(f["q_eff"]),read(f["q_physical"])
            h,T=cx(f,"h"),cx(f,"Tplus");Gamma,S,Rtotal,GR=(cx(f,key) for key in ("Gamma","localR","totalR","GR"))
            Asub,Aeff,Aphysical=(cx(f,key) for key in ("A_sub","A_eff","A_physical"))
            Beff,Bphysical=(cx(f,key) for key in ("Gamma_eff","Gamma_physical_total"))
            nb=size(Gamma,3);nk=length(wk);nb==5 && nk==11 || error("expected S01 5x11 frame")
            report["Tminus_minus_Tplus_adjoint_F"]=norm(Tminus-T')
            report["native_coupling_contract"]="accepted Core production_operators.jl:370,384 finite_chain uses left adjoint(Tplus), right Tplus; basis Tminus is recorded only as diagnostic"
            bad=findall(<(0.0),q);isempty(bad) && error("no negative q in supplied obstruction witness")
            allk_energies=unique([first(bad),argmin(q)])
            all(e->e in selected,allk_energies) || error("first/worst bad energy missing")
            nodes=[(findfirst(==(e),selected),m) for e in allk_energies for m in 1:nk]
            push!(nodes,(length(selected),nk));push!(nodes,(argmin(abs.(energy[selected].-0.085)),1));unique!(nodes)
            variants=(("raw",false,false),("HermLambda",false,true),("HermH",true,false),("HermBoth",true,true))
            qrows=Dict{String,Any}();qhigh=Dict{String,Vector{BigFloat}}()
            worstj=findfirst(==(argmin(q)),selected)
            worstm=argmin([wk[m]*real(tr(Asub[worstj,m,:,:])) for m in 1:nk])
            rowscsv=String["energy_index,k_index,energy_eV,variant,float_A_trace,big_A_trace,big_physical_trace,big_effective_trace"]
            HDF5.h5open(joinpath(output,"witnesses.h5"),"w") do out
                out["schema"]="campaign133-warm-selected-45-oracles-v1"
                for (j,m) in nodes
                    checkwall();e=selected[j];name="e$(e)_k$(m)";group=HDF5.create_group(out,name)
                    k=Matrix(Gamma[j,m,:,:]);s=Matrix(S[j,m,:,:]);hm=Matrix(h[m,:,:]);native=Matrix(GR[j,m,:,:])
                    datum=Dict{String,Any}("energy_index"=>e,"k_index"=>m,"energy_eV"=>energy[e]*E0,
                        "input_Gamma"=>psd_info(k),"h_antiHerm_F"=>norm(hm-hm'),
                        "reconstructed_Lambda_antiHerm_F"=>norm((s+im*k/2)-(s+im*k/2)'),
                        "local_broadening_minus_Gamma_F"=>norm(im*(s-s')-k),
                        "same_G_sub_vs_effective"=>distance(Asub[j,m,:,:],Aeff[j,m,:,:]),
                        "same_G_sub_vs_physical"=>distance(Asub[j,m,:,:],Aphysical[j,m,:,:]),"variants"=>Dict{String,Any}())
                    for (label,X) in (("Gamma",k),("localR",s),("h",hm),("native_side_h_F64",(hm+hm')/2),("Tplus",T),("totalR",Rtotal[j,m,:,:]),
                                      ("native_GR",native),("A_sub",Asub[j,m,:,:]),("A_eff",Aeff[j,m,:,:]),("A_physical",Aphysical[j,m,:,:]),
                                      ("Gamma_eff",Beff[j,m,:,:]),("Gamma_physical_total",Bphysical[j,m,:,:]))
                        putcx(group,label,X)
                    end
                    for (label,hh,ll) in variants
                        checkwall();fd=float_oracle(energy[e],hm,T,s,k;herm_h=hh,herm_lambda=ll)
                        bd=big_oracle(energy[e],hm,T,s,k;herm_h=hh,herm_lambda=ll)
                        sub=HDF5.create_group(group,label)
                        for (key,X) in (("float_R",fd.R),("float_A",fd.A),("float_physical",fd.physical),("float_effective",fd.effective),
                                        ("big_R",bd.R),("big_A",bd.A),("big_physical",bd.physical),("big_effective",bd.effective),("big_central_row",bd.row))
                            putcx(sub,key,X)
                        end
                        metric=Dict{String,Any}("big_exact_traces"=>bd.traces,"float_inverse_residual"=>fd.residual,
                            "float_R_vs_native"=>distance(fd.R,native),"big_R_vs_native"=>distance(bd.R,native),
                            "float_vs_downcast_big_A"=>distance(fd.A,bd.A),"big_A_vs_physical"=>distance(bd.A,bd.physical),
                            "big_A_vs_effective"=>distance(bd.A,bd.effective),"big_A"=>psd_info(bd.A),"big_physical"=>psd_info(bd.physical))
                        datum["variants"][label]=metric
                        if j==worstj && m==worstm
                            high=big_oracle(energy[e],hm,T,s,k;herm_h=hh,herm_lambda=ll,bits=384)
                            metric["repeat384_traces"]=high.traces
                            metric["trace192_vs384"]=setprecision(BigFloat,384) do
                                a=parse(BigFloat,bd.traces["A_trace"]);b=parse(BigFloat,high.traces["A_trace"])
                                Dict("absolute"=>string(abs(a-b)),"relative"=>string(rel(abs(a-b),abs(b))))
                            end
                        end
                        key="e$(e)_"*label
                        row=get!(qrows,key,Dict{String,Any}("energy_index"=>e,"variant"=>label,"k_count"=>0,"q_downcast_matrix_A"=>0.0,"q_downcast_matrix_physical"=>0.0,"q_downcast_matrix_effective"=>0.0))
                        w=(gs/(2pi))*we[e]*wk[m];row["k_count"]+=1
                        row["q_downcast_matrix_A"]+=w*real(tr(bd.A));row["q_downcast_matrix_physical"]+=w*real(tr(bd.physical));row["q_downcast_matrix_effective"]+=w*real(tr(bd.effective))
                        setprecision(BigFloat,192) do
                            accumulator=get!(qhigh,key,BigFloat[0,0,0])
                            weight=BigFloat(gs/(2pi))*BigFloat(we[e])*BigFloat(wk[m])
                            for (i,k) in enumerate(("A_trace","physical_trace","effective_trace"))
                                accumulator[i]+=weight*parse(BigFloat,bd.traces[k])
                            end
                        end
                        push!(rowscsv,join((e,m,energy[e]*E0,label,real(tr(fd.A)),bd.traces["A_trace"],bd.traces["physical_trace"],bd.traces["effective_trace"]),","))
                    end
                    report["nodes"][name]=datum
                end
            end
            for (key,row) in qrows
                e=row["energy_index"];row["scope"]=row["k_count"]==nk ? "all k: whole q at selected E" : "selected k only: partial mass, not whole q"
                row["q_192_A"],row["q_192_physical"],row["q_192_effective"]=string.(qhigh[key])
                row["native_whole_q_raw"]=q[e];row["native_whole_q_eff"]=qe[e];row["native_whole_q_physical"]=qp[e]
            end
            report["selected_q_aggregates"]=qrows
            report["precision_scope"]="signed q uses192bit trace strings and192bit weighted sums; matrices in HDF areFloat64 downcasts;384repeat only strongest negative k at worst q energy, other nodes have192bit LU/identity residual evidence, not per-node precision convergence"
            report["global_q_negative_count"]=length(bad);report["global_q_min"]=minimum(q)
            report["global_negative_absolute_mass"]=-sum(q[bad]);report["global_signed_mass"]=sum(q);report["global_absolute_mass"]=sum(abs,q)
            report["execution"]="complete";write(joinpath(output,"traces.csv"),join(rowscsv,"\n")*"\n")
        end
    catch error
        report["execution"]="stopped";report["reason"]=sprint(showerror,error,catch_backtrace());rethrow()
    finally
        report["elapsed_seconds"]=time()-started
        open(joinpath(output,"summary.toml"),"w") do io;TOML.print(io,report);end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    ARGS==["--tiny"] ? tiny_checks() : (length(ARGS)==2 ? run_oracles(abspath(ARGS[1]),abspath(ARGS[2])) : error("oracles.jl --tiny OR INPUT NEW_OUTPUT"))
end
