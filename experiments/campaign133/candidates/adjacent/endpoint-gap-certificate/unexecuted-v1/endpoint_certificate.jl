# Read-only frozen-hat endpoint/strip enclosure. No Green/Phi/SCBA calculation.
const START=time()
using LinearAlgebra,HDF5,TOML,SHA,Dates
const V4PATH=normpath(joinpath(@__DIR__,"../gamma-density-hermitian-v4/matsubara_oracle.jl"))
bytes2hex(open(sha256,V4PATH))=="814dff509cef6cdb1e0e1d1a762cabbd885a761c2a5b4964f80373df2b01b21a" || error("V4 source")
module V4
include("../gamma-density-hermitian-v4/matsubara_oracle.jl")
end
const Iv=Tuple{Float64,Float64}
Q(x::Float64)=rationalize(BigInt,x;tol=0)
point(x)=(Float64(x),Float64(x))
function outward(l,u)
    all(isfinite,(l,u)) || error("nonfinite interval")
    p,q=prevfloat(l),nextfloat(u);all(isfinite,(p,q)) || error("interval overflow");return p,q
end
iadd(a::Iv,b::Iv)=outward(a[1]+b[1],a[2]+b[2])
ineg(a::Iv)=(-a[2],-a[1])
isub(a::Iv,b::Iv)=iadd(a,ineg(b))
function imul(a::Iv,b::Iv)
    p=(a[1]*b[1],a[1]*b[2],a[2]*b[1],a[2]*b[2]);return outward(minimum(p),maximum(p))
end
function ipdiv(a::Iv,b::Iv)
    a[1]>=0 && b[1]>0 || error("nonpositive division domain")
    l,u=outward(a[1]/b[2],a[2]/b[1]);return max(0.0,l),u
end
function qinterval(l,u=l)
    x,y=Float64(l),Float64(u)
    while Q(x)>l;x=prevfloat(x);end
    while Q(y)<u;y=nextfloat(y);end
    return x,y
end
function constants()
    function atanbox(q)
        s=sum((-1)^m//(big(q)^(2m+1)*(2m+1)) for m in 0:31)
        next=1//(big(q)^65*65);return s,s+next # 32 terms, next term positive
    end
    a,b=atanbox(5);c,d=atanbox(239);pi=qinterval(16a-4d,16b-4c)
    s=2sum(1//(big(3)^(2m+1)*(2m+1)) for m in 0:39)
    tail=(2//(big(3)^81*81))/(1-1//9);log2=qinterval(s,s+tail)
    return pi,ipdiv(log2,pi)
end
function kernel(t::Iv,pi::Iv,k1::Iv;one=false)
    one && return k1
    t[1]>1 || error("kernel requires t>1 or exact endpoint1")
    u=ipdiv(point(1),t);u2=imul(u,u);p=u;s=point(0)
    for m in 0:23
        s=iadd(s,ipdiv(p,point((2m+1)*(2m+2))));p=imul(p,u2)
    end
    tail=ipdiv(p,imul(point(49*50),isub(point(1),u2)))
    return ipdiv((s[1],iadd(s,tail)[2]),pi)
end
function scale(c::Iv,x::Float64)
    return x>=0 ? outward(c[1]*x,c[2]*x) : outward(c[2]*x,c[1]*x)
end
function gap(H,L,U,edge::Iv,sgn)
    n=size(H,1);nb=size(L,1);result=Inf
    for i in 1:n
        cell=(i-1)÷nb;a=(i-1)%nb+1;radius=0.0;diagonal=0.0
        for j in 1:n
            b=(j-1)%nb+1;localblock=(j-1)÷nb==cell
            re=point(sgn*real(H[i,j]));im=point(sgn*imag(H[i,j]))
            if localblock
                re=isub(re,(real(L[a,b]),real(U[a,b])));im=isub(im,(imag(L[a,b]),imag(U[a,b])))
            end
            if i==j
                diagonal=iadd(re,sgn==1 ? ineg(edge) : edge)[1]
            else
                term=nextfloat(max(abs(re[1]),abs(re[2]))+max(abs(im[1]),abs(im[2])))
                radius=nextfloat(radius+term)
            end
        end
        result=min(result,prevfloat(diagonal-radius))
    end
    return result
end
function area(c,h,x)
    x<=c-h && return zero(h)
    x>=c+h && return h
    return x<c ? (x-c+h)^2/(2h) : h-(c+h-x)^2/(2h)
end
function ghostbound(left,right,dL,dR,pi)
    dL>0 && dR>0 || error("nonpositive strip gap; never abs/floor")
    dl=prevfloat(dL^2);dr=prevfloat(dR^2);dl>0 && dr>0 || error("squared gap underflow")
    beta=nextfloat(nextfloat(left/dl)+nextfloat(right/dr));beta=nextfloat(beta/prevfloat(2pi[1]))
    isfinite(beta) || error("mass bound overflow");return beta
end
function main()
    length(ARGS)==2 || error("usage: AUTHORIZED_INPUT NEW_OUTPUT")
    config,output=abspath.(ARGS);o=TOML.parsefile(config);base=@__DIR__
    frozen=TOML.parsefile(joinpath(base,"input.toml"))
    all(k=="execution_authorized" || (haskey(o,k) && isequal(o[k],v)) for (k,v) in frozen) || error("frozen option changed")
    o["execution_authorized"]===true || error("write-only/HOLD")
    Threads.nthreads()==1 && !get_zero_subnormals() && rounding(Float64)==RoundNearest || error("IEEE/J1 domain")
    BLAS.set_num_threads(1);ispath(output) && error("output exists")
    member=only(filter(s->startswith(s,"0::"),readlines("/proc/self/cgroup")))
    cg=joinpath("/sys/fs/cgroup",lstrip(split(member,':';limit=3)[3],'/'))
    strip(read(joinpath(cg,"memory.max"),String))=="2147483648" && strip(read(joinpath(cg,"memory.swap.max"),String))=="0" || error("2GiB/Swap0")
    cpu=split(read(joinpath(cg,"cpu.max"),String));cpu[1]!="max" && parse(Int,cpu[1])<=parse(Int,cpu[2]) || error("CPU100")
    cfg=V4.guarded_file(base,o,"V4_input");c=TOML.parsefile(cfg);cb=dirname(cfg)
    summarypath=V4.guarded_file(base,o,"V4_summary");summary=TOML.parsefile(summarypath)
    summary["status"]=="completed" && summary["script_sha256"]==c["oracle_script_sha256"] && summary["input_sha256"]==o["V4_input_sha256"] || error("actual V4 binding")
    V4.sha_file(Base.active_project())==c["project_sha256"] && V4.sha_file(joinpath(dirname(Base.active_project()),"Manifest.toml"))==c["manifest_sha256"] || error("environment")
    keys=("terminal_h5","model_h5","source_input","source_identity","source_writer","map_summary","numerical_parent","checkpoint","spectral_selected","spectral_writer","source_ledger","source_summary")
    paths=Dict(k=>V4.guarded_file(cb,c,k) for k in keys);pd=V4.require_all_PD_certificate(c,cb,paths)
    mkpath(output);cp(config,joinpath(output,"input.toml"));cp(@__FILE__,joinpath(output,"endpoint_certificate.jl"))
    cp(V4PATH,joinpath(output,"V4-matsubara_oracle.jl"));cp(cfg,joinpath(output,"V4-input.toml"));cp(summarypath,joinpath(output,"V4-summary.toml"))
    report=Dict{String,Any}("status"=>"running","scientific_acceptance"=>false,"scope"=>"same V4 dyadic-centre hats; IEEE outward enclosures, not hardware formal verification",
        "input_sha256"=>V4.sha_file(config),"script_sha256"=>V4.sha_file(@__FILE__),"V4_input_sha256"=>o["V4_input_sha256"],"V4_summary_sha256"=>o["V4_summary_sha256"])
    checkstage(s)=(time()-START<120 || error("120s cap"))
    try
        f=V4.source_frame(paths,c,TOML.parsefile(paths["source_input"]),TOML.parsefile(paths["map_summary"]),checkstage)
        gamma=f.gamma::Array{ComplexF64,4};centres=f.canonical::Vector{Float64};step=f.step::Float64
        hcell=f.h::Array{ComplexF64,3};T=f.basis["Tplus"]::Matrix{ComplexF64}
        ne,nk,nb=size(gamma)[1:3];pi,k1=constants();h=point(step)
        all(0<isub(point(centres[j+1]),point(centres[j]))[1] && isub(point(centres[j+1]),point(centres[j]))[2]<prevfloat(2f.step) for j in 1:ne-1) || error("hat support gap")
        ell=qinterval(Q(first(centres))-Q(f.step));rr=qinterval(Q(last(centres))+Q(f.step));A=first(f.e);B=last(f.e)
        ell[2]<A<B<rr[1] || error("native window/support order")
        M=ipdiv(point(1),pi)[2];lo=zeros(ComplexF64,4,nk,nb,nb);hi=similar(lo);fill!(hi,0)
        left=zeros(nk);right=zeros(nk);qh=Q(f.step);qa=Q(A);qb=Q(B)
        for j in 1:ne
            j%1024==0 && checkstage("endpoint_pass")
            tl=ipdiv(isub(point(centres[j]),ell),h);tr=ipdiv(isub(rr,point(centres[j])),h)
            dl=ipdiv((max(0.0,isub(point(centres[j]),point(A))[1]),max(0.0,isub(point(centres[j]),point(A))[2])),h)
            dr=ipdiv((max(0.0,isub(point(B),point(centres[j]))[1]),max(0.0,isub(point(B),point(centres[j]))[2])),h)
            coeff=(kernel(tl,pi,k1;one=j==1),kernel(tr,pi,k1;one=j==ne),
                point(dl[1]>1 ? kernel(dl,pi,k1)[2] : M),point(dr[1]>1 ? kernel(dr,pi,k1)[2] : M))
            al=isub(point(centres[j]),h)[1]<A ? qinterval(area(Q(centres[j]),qh,qa))[2] : 0.0
            ar=iadd(point(centres[j]),h)[2]>B ? qinterval(qh-area(Q(centres[j]),qh,qb))[2] : 0.0
            for k in 1:nk
                trace=0.0
                for b in 1:nb,a in 1:nb
                    z=gamma[j,k,a,b];a==b && (trace=nextfloat(trace+real(z)))
                    for q in 1:4
                        re=scale(coeff[q],real(z));im=scale(coeff[q],imag(z))
                        lo[q,k,a,b]=complex(prevfloat(real(lo[q,k,a,b])+re[1]),prevfloat(imag(lo[q,k,a,b])+im[1]))
                        hi[q,k,a,b]=complex(nextfloat(real(hi[q,k,a,b])+re[2]),nextfloat(imag(hi[q,k,a,b])+im[2]))
                    end
                end
                al>0 && (left[k]=nextfloat(left[k]+nextfloat(trace*al)))
                ar>0 && (right[k]=nextfloat(right[k]+nextfloat(trace*ar)))
            end
        end
        records=Dict{String,Any}();sheet=0.0;allgaps=true
        for k in 1:nk
            H=V4.chain_hamiltonian(Matrix(view(hcell,k,:,:)),T,4)
            gaps=[gap(H,Matrix(lo[q,k,:,:]),Matrix(hi[q,k,:,:]),(ell,rr,point(A),point(B))[q],q%2==1 ? 1 : -1) for q in 1:4]
            all(isfinite,gaps) || error("nonfinite Gershgorin gap")
            r=Dict{String,Any}("lower_gaps"=>gaps,"status"=>"unavailable_nonpositive_Gershgorin_bound")
            if all(>(0),gaps)
                beta=ghostbound(left[k],right[k],gaps[3],gaps[4],pi)
                r["status"]="positive_interval_Gershgorin_gaps";r["ghost_mass_operator_upper"]=beta
                sheet=nextfloat(sheet+nextfloat(f.wk[k]*beta))
            else;allgaps=false;end
            records[string(k)]=r
        end
        isfinite(sheet) || error("sheet-bound overflow")
        report["all_gaps_positive"]=allgaps;report["per_k"]=records;report["PD_binding"]=pd;report["array_sha256"]=c["array_sha256"]
        report["gap_order"]=["exterior_left","exterior_right","strip_left","strip_right"];report["units"]=c["units"]
        report["constant_intervals"]=Dict("pi"=>collect(pi),"k1_log2_over_pi"=>collect(k1),"M_upper"=>M)
        if allgaps
            report["ghost_sheet_number_upper"]=nextfloat(nextfloat(f.gs*nb)*sheet)
            report["ghost_sheet_over_target_upper"]=nextfloat(report["ghost_sheet_number_upper"]/c["units"]["target_dimensionless"])
        end
        h5open(joinpath(output,"endpoint_enclosures.h5"),"w") do g
            V4.write_complex(g,"positive_weighted_Gamma_lower",lo);V4.write_complex(g,"positive_weighted_Gamma_upper",hi)
            g["support_left_interval"]=collect(ell);g["support_right_interval"]=collect(rr);g["native_window"]=[A,B]
            g["q_order"]=["outer_left_positive","outer_right_positive","strip_left_upper","strip_right_upper"]
        end
        report["status"]="completed"
    catch e;report["status"]="failed";report["reason"]=sprint(showerror,e);rethrow()
    finally
        report["elapsed_seconds"]=time()-START
        open(joinpath(output,"summary.toml"),"w") do io;TOML.print(io,report);end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;main();end
