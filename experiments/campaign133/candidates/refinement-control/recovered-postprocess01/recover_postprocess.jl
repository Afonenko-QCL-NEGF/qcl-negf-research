# Recovery ONLY: stdlib/HDF5, no Core/Γ module/context/Π/inverse imports.
using HDF5, LinearAlgebra, TOML, SHA, Test
sha_file(p)=bytes2hex(open(sha256,p))
sha_array(A)=bytes2hex(sha256(reinterpret(UInt8,vec(A))))
readrows(f,name,rows)=complex.(f[name*"_real"][rows,:,:,:],f[name*"_imag"][rows,:,:,:])
ratio(a,b)=b==0 ? (a==0 ? 0.0 : Inf) : a/b
function source_precondition(s)
    s["Pi_started"]==2 && s["Pi_completed"]==2 || error("both Π must already be complete")
    s["execution"]=="stopped" && s["stage"]=="Pi_complete" || error("source failure stage not admitted")
    s["reason"]=="InexactError: Int64(1.2)" || error("only known postprocessing failure admitted")
    haskey(s,"coarse_map") && haskey(s,"fine_map") || error("both map records required")
    return true
end
fine_index(e::Integer)=2 * e - 1 # explicit *, never a numeric literal adjacent to e
function selected_indices(n,anchors,step,halo)
    n>=2 && step>=1 && halo>=0 || error("index controls")
    out=Int[]
    for a in anchors,shift in -cld(n,step):cld(n,step),s in -halo:halo,delta in -1:1
        e=a+step*shift+s; q=fine_index(e)+delta
        1<=e<=n && 1<=q<=2 * n-1 && push!(out,q)
    end
    return sort!(unique!(out))
end
source_rows(rows)=((first(rows)+1)÷2):cld(last(rows)+1,2)
function prolong_chunk(X,offset,rows,Ec,Ef)
    Y=zeros(ComplexF64,length(rows),size(X)[2:end]...)
    for (i,e) in enumerate(rows)
        j=(e+1)÷2; localj=j-offset+1
        if isodd(e); Y[i,:,:,:].=X[localj,:,:,:]
        else
            theta=(Ef[e]-Ec[j])/(Ec[j+1]-Ec[j])
            isfinite(theta) && 0<=theta<=1 || error("nonconvex actual-node weight")
            Y[i,:,:,:].=(1-theta).*X[localj,:,:,:].+theta.*X[localj+1,:,:,:]
        end
    end
    return Y
end
mutable struct NormRecord
    n2::Float64; s2::Float64; nw2::Float64; sw2::Float64
    maxblock2::Float64; energy::Int; k::Int
end
NormRecord()=NormRecord(0,0,0,0,0,0,0)
function accumulate!(d,ref,value,we,wk,rows;mode=:all)
    for (i,e) in enumerate(rows)
        mode===:even && isodd(e) && continue
        mode===:coarse && iseven(e) && continue
        index=mode===:coarse ? (e+1)÷2 : e
        for k in axes(ref,2)
            bn=0.0; bs=0.0
            for b in axes(ref,4),a in axes(ref,3)
                bn+=abs2(value[i,k,a,b]-ref[i,k,a,b]); bs+=abs2(ref[i,k,a,b])
            end
            w=we[index]*wk[k]; d.n2+=bn; d.s2+=bs; d.nw2+=w*bn; d.sw2+=w*bs
            if bn>d.maxblock2; d.maxblock2=bn; d.energy=index; d.k=k; end
        end
    end; return d
end
function metrics(d)
    Dict{String,Any}("absolute"=>sqrt(d.n2),"scale"=>sqrt(d.s2),"relative"=>ratio(sqrt(d.n2),sqrt(d.s2)),
        "weighted_absolute"=>sqrt(d.nw2),"weighted_scale"=>sqrt(d.sw2),"weighted_relative"=>ratio(sqrt(d.nw2),sqrt(d.sw2)),
        "maximum_block_absolute"=>sqrt(d.maxblock2),"energy_index"=>d.energy,"k_index"=>d.k)
end
function replay(current,reference,tol)
    isfinite(current) && isfinite(reference) && isfinite(tol) && tol>=0 || return false
    current==reference && return true
    return reference!=0 && abs(current-reference)<=tol*abs(reference)
end
function trapezoid(E)
    w=fill((last(E)-first(E))/(length(E)-1),length(E)); w[1]/=2; w[end]/=2; w
end
function tiny_checks()
    @testset "recovery index/stream controls" begin
        @test fine_index(3)==5 && fine_index(9785)==19569
        @test selected_indices(9,[3],99,0)==[4,5,6]
        @test 15813 in selected_indices(30704,[9785],939,0)
        @test_throws InexactError Int(2e-1+1) # regression: bad spelling evaluates to1.2
        s=Dict("Pi_started"=>2,"Pi_completed"=>2,"execution"=>"stopped","stage"=>"Pi_complete",
            "reason"=>"InexactError: Int64(1.2)","coarse_map"=>Dict(),"fine_map"=>Dict())
        @test source_precondition(s)
        t=copy(s); t["Pi_completed"]=1; @test_throws Exception source_precondition(t)
        t=copy(s); t["reason"]="other"; @test_throws Exception source_precondition(t)
        Ec=collect(range(-.3,.9;length=5)); Ef=collect(range(-.3,.9;length=9))
        X=reshape(ComplexF64.(1:40),5,2,2,2).*(1+.2im); P=zeros(ComplexF64,9,2,2,2)
        for firste in 1:2:9
            rows=firste:min(firste+1,9); sr=source_rows(rows)
            P[rows,:,:,:].=prolong_chunk(X[sr,:,:,:],first(sr),rows,Ec,Ef)
        end
        @test P[1:2:end,:,:,:]==X
        V=P.+.03im; wf=trapezoid(Ef); wk=[.2,.8]; d=NormRecord(); mid=NormRecord(); c=NormRecord()
        for firste in 1:2:9
            rows=firste:min(firste+1,9)
            accumulate!(d,P[rows,:,:,:],V[rows,:,:,:],wf,wk,rows)
            accumulate!(mid,P[rows,:,:,:],V[rows,:,:,:],wf,wk,rows;mode=:even)
            accumulate!(c,P[rows,:,:,:],V[rows,:,:,:],trapezoid(Ec),wk,rows;mode=:coarse)
        end
        @test metrics(d)["absolute"]≈norm(V-P) atol=1e-14
        expected=sum(wf[e]*wk[k]*sum(abs2,V[e,k,:,:]-P[e,k,:,:]) for e in 1:9,k in 1:2)
        @test metrics(d)["weighted_absolute"]^2≈expected atol=1e-14
        @test metrics(mid)["absolute"]≈norm(V[2:2:end,:,:,:]-P[2:2:end,:,:,:]) atol=1e-14
        @test metrics(c)["absolute"]≈norm(V[1:2:end,:,:,:]-X) atol=1e-14
        @test metrics(d)["energy_index"] isa Int && eltype(selected_indices(9,[3],99,0))===Int
        @test replay(0.0,0.0,1e-10) && !replay(1e-12,0.0,1e-10)
        @test !replay(NaN,1.0,1e-10) && !replay(1.0,1.0,NaN)
    end
end
function run_recovery(input,out)
    cfg=TOML.parsefile(input); cfg["runtime_authorized"]===true || error("ROOT runtime authorization required")
    Threads.nthreads()==1 || error("Julia1 required"); BLAS.set_num_threads(1)
    0<cfg["soft_seconds"]<=170 && 0<cfg["output_cap_bytes"]<=64*1024^2 || error("recovery resource controls")
    cfg["memory_max_bytes"]==3*1024^3 && cfg["external_seconds"]==180 || error("recovery envelope")
    cfg["chunk_rows"] isa Int && 1<=cfg["chunk_rows"]<=4096 || error("stream chunk control")
    cfg["anchor_indices"]==[9785] && cfg["E0_eV"]==0.1 || error("recovery identity controls")
    entry=only(filter(s->startswith(s,"0::"),readlines("/proc/self/cgroup")))
    cg=joinpath("/sys/fs/cgroup",lstrip(split(entry,':';limit=3)[3],'/'))
    strip(read(joinpath(cg,"memory.max"),String))==string(cfg["memory_max_bytes"]) || error("external memory cap")
    strip(read(joinpath(cg,"memory.swap.max"),String))=="0" || error("swap cap")
    cpu=split(read(joinpath(cg,"cpu.max"),String)); cpu[1]!="max" && parse(Int,cpu[1])<=parse(Int,cpu[2]) || error("CPU cap")
    ispath(out) && error("new recovery directory required"); mkpath(out)
    source=normpath(joinpath(@__DIR__,cfg["source_attempt"])); summarypath=joinpath(source,"summary.toml")
    old=TOML.parsefile(summarypath); source_precondition(old)
    sha_file(summarypath)==cfg["source_summary_sha256"] || error("immutable source summary SHA differs")
    begun=time(); result=Dict{String,Any}("status"=>"running_postprocess_only","new_Pi_calls"=>0,
        "original_Pi_started"=>2,"original_Pi_completed"=>2,"original_failure"=>old["reason"],
        "source_summary_sha256"=>sha_file(summarypath),"postprocessor_sha256"=>sha_file(@__FILE__),"input_sha256"=>sha_file(input))
    function checkpoint(stage)
        result["stage"]=stage; result["seconds"]=time()-begun; result["maxrss_bytes"]=Sys.maxrss()
        open(joinpath(out,"summary.toml"),"w") do io; TOML.print(io,result); end
        sum(filesize(joinpath(d,n)) for (d,_,fs) in walkdir(out) for n in fs;init=0)<=cfg["output_cap_bytes"] || error("output cap")
        time()-begun<=cfg["soft_seconds"] || error("soft time cap")
        println("stage=",stage," elapsed=",time()-begun); flush(stdout)
    end
    try
        hashes=Dict{String,Any}()
        cp(@__FILE__,joinpath(out,"recover_postprocess.jl")); cp(input,joinpath(out,"input.toml"))
        paths=Dict("coarse.h5"=>joinpath(source,"coarse.h5"),"fine-input.h5"=>joinpath(source,"fine-input.h5"),
            "fine.h5"=>joinpath(source,"fine.h5"),"model.h5"=>normpath(joinpath(@__DIR__,cfg["model"])))
        for (name,path) in paths
            digest=sha_file(path); expected=cfg["expected_sha256"][name]
            expected=="COMPUTE_AT_RUNTIME" || digest==expected || error("source SHA differs: "*name)
            hashes[name]=Dict("sha256"=>digest,"scope"=>expected=="COMPUTE_AT_RUNTIME" ? "ROOT-authorized runtime whole-file hash" : "expected SHA matched")
        end
        for (name,bytes) in (("coarse.h5",3*30704*11*25*16),("fine-input.h5",61407*11*25*16),("fine.h5",61407*11*25*16))
            filesize(paths[name])>=bytes || error("source field file smaller than raw dataset payload")
        end
        for name in ("run_control.jl","operators.jl","tiny_checks.jl")
            sha_file(joinpath(source,name))==old["source_code_sha256"][name] || error("immutable attempt source differs")
        end
        sha_file(joinpath(source,"input.toml"))==old["input_sha256"] || error("attempt input differs")
        old["identity"]["core_head"]==cfg["core_head"] || error("original Core differs")
        result["source_hashes"]=hashes; checkpoint("whole_file_identity")
        wk=h5open(paths["model.h5"],"r") do m; read(m["grid_wᵏ"]); end
        sha_array(wk)==old["coarse_model_binding"]["grid_wᵏ"] || error("k weights digest differs")
        all(isfinite,wk) && all(>(0),wk) && length(wk)==11 || error("k weights")
        h5open(paths["coarse.h5"],"r") do coarse
          h5open(paths["fine-input.h5"],"r") do gamma
           h5open(paths["fine.h5"],"r") do fine
            Ec=read(coarse["energy"]); Ef=read(fine["energy"]); wc=read(coarse["weight_energy"]); wf=read(fine["weight_energy"])
            length(Ec)==30704 && length(Ef)==61407 && Ef[1:2:end]==Ec || error("nested source axes")
            all(isfinite,Ec) && all(>(0),diff(Ec)) && all(isfinite,Ef) && all(>(0),diff(Ef)) || error("finite sorted axes")
            read(gamma["energy"])==Ef && read(gamma["weight_energy"])==wf || error("fine input/output axes differ")
            wc==trapezoid(Ec) && wf==trapezoid(Ef) || error("quadrature weights differ")
            for (f,names,ne) in ((coarse,("Gamma","Phi","Lambda"),30704),(gamma,("Gamma",),61407),(fine,("Phi",),61407))
                for name in names,suffix in ("_real","_imag"); size(f[name*suffix])==(ne,11,5,5) || error("incomplete field: "*name); end
            end
            allD=NormRecord(); midD=NormRecord(); coarseD=NormRecord(); prolongG=NormRecord()
            rawC=NormRecord(); rawF=NormRecord(); identity=NormRecord()
            selected=selected_indices(30704,cfg["anchor_indices"],939,3); selectedset=Set(selected)
            open(joinpath(out,"selected-defect.csv"),"w") do io
                println(io,"e,k,E_eV,raw_r2,collocation_defect_r2,Phi_reference_r2,weighted_raw_r2,weighted_defect_r2,Gamma_trace_perNb_eV")
                for start in 1:cfg["chunk_rows"]:61407
                    rows=start:min(start+cfg["chunk_rows"]-1,61407); sr=source_rows(rows)
                    Cg=readrows(coarse,"Gamma",sr); Cp=readrows(coarse,"Phi",sr)
                    Pg=prolong_chunk(Cg,first(sr),rows,Ec,Ef); Pp=prolong_chunk(Cp,first(sr),rows,Ec,Ef)
                    G=readrows(gamma,"Gamma",rows); P=readrows(fine,"Phi",rows)
                    all(A->all(isfinite,A),(Cg,Cp,Pg,Pp,G,P)) || error("nonfinite persisted field")
                    accumulate!(allD,Pp,P,wf,wk,rows); accumulate!(midD,Pp,P,wf,wk,rows;mode=:even)
                    accumulate!(coarseD,Pp,P,wc,wk,rows;mode=:coarse); accumulate!(prolongG,Pg,G,wf,wk,rows)
                    accumulate!(rawC,Pp,Pg,wc,wk,rows;mode=:coarse); accumulate!(rawF,P,G,wf,wk,rows)
                    accumulate!(identity,P-G,(Pp-Pg)+(P-Pp),wf,wk,rows)
                    for (i,e) in enumerate(rows)
                        e in selectedset || continue
                        for k in 1:11
                            r2=sum(abs2,P[i,k,:,:]-G[i,k,:,:]); d2=sum(abs2,P[i,k,:,:]-Pp[i,k,:,:]); w=wf[e]*wk[k]
                            println(io,join((e,k,Ef[e]*cfg["E0_eV"],r2,d2,sum(abs2,Pp[i,k,:,:]),w*r2,w*d2,real(tr(G[i,k,:,:]))/5*cfg["E0_eV"]),','))
                        end
                    end
                    checkpoint("stream_rows_"*string(last(rows)))
                end
            end
            for (label,m,key) in (("coarse_raw_replay",rawC,"coarse_map"),("fine_raw_replay",rawF,"fine_map"))
                current=metrics(m); reference=old[key]["raw_residual"]
                for field in ("absolute","scale","relative","weighted_absolute","weighted_scale","weighted_relative")
                    replay(current[field],reference[field],cfg["replay_relative"]) || error("raw norm replay differs: "*label*field)
                end; result[label]=current
            end
            result["Gamma_prolongation_replay"]=metrics(prolongG)
            metrics(prolongG)["relative"]<=cfg["replay_relative"] || error("persisted fine Gamma differs from P source Gamma")
            metrics(prolongG)["weighted_relative"]<=cfg["replay_relative"] || error("weighted Gamma prolongation replay")
            metrics(identity)["relative"]<=cfg["replay_relative"] || error("residual collocation identity")
            result["full_fine_defect"]=metrics(allD); result["new_midpoint_defect"]=metrics(midD)
            result["coarse_restriction_defect"]=metrics(coarseD); result["residual_identity"]=metrics(identity)
            result["selected_fine_indices"]=selected; result["E0_eV"]=cfg["E0_eV"]
           end
          end
        end
        result["original_number_diagnostics_not_recomputed"]=Dict("coarse"=>old["coarse_map"]["fixed_mu_diagnostic"],"fine"=>old["fine_map"]["fixed_mu_diagnostic"])
        result["status"]="completed_streaming_recovery_no_new_maps"; checkpoint("complete")
    catch e
        result["status"]="recovery_failed"; result["error"]=sprint(showerror,e); rethrow()
    finally
        result["seconds"]=time()-begun; result["maxrss_bytes"]=Sys.maxrss()
        open(joinpath(out,"summary.toml"),"w") do io; TOML.print(io,result); end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS==["--tiny"]; tiny_checks()
    else
        length(ARGS)==2 || error("usage: recover_postprocess.jl --tiny | RECOVERY_INPUT NEW_OUTPUT")
        run_recovery(ARGS[1],abspath(ARGS[2]))
    end
end
