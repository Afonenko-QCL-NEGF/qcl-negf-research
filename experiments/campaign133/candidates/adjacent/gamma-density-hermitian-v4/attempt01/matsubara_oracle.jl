# Write-only prototype. Execute only after root review and resource-green.
using LinearAlgebra, HDF5, TOML, SHA, Dates

sha_file(p) = bytes2hex(open(sha256, p))
sha_array(x) = bytes2hex(sha256(reinterpret(UInt8, vec(x))))
read_complex(g, n) = complex.(read(g[n*"_real"]), read(g[n*"_imag"]))
read_k(g, n, k) = complex.(g[n*"_real"][:,k,:,:], g[n*"_imag"][:,k,:,:])
herm(X) = (X + X') / 2
function write_complex(g, n, x)
    g[n*"_real"] = real.(x); g[n*"_imag"] = imag.(x)
end
function guarded_file(base, c, key)
    p = normpath(joinpath(base,c[key]))
    sha_file(p) == c[key*"_sha256"] || error("file guard failed: $key")
    return p
end

# Integral of a full unit triangular hat, with the 1/(2π) already included.
function hat_series(t)
    u = inv(t); u2 = u*u
    return u * (1 + u2/6 + u2*u2/15) / (2*oftype(real(t),π))
end
function hat_remainder(t)
    r = abs(t); r > 1 || error("hat series outside convergence domain")
    return inv(56*oftype(r,π)*r^7*(1-inv(r*r)))
end
function hat_exact_big(t)
    setprecision(256) do
        z = Complex{BigFloat}(t)
        F(x) = x*log(x)
        return (F(z+1)-2F(z)+F(z-1))/(2big(π))
    end
end
function kernel_controls(min_imag, real_extent)
    records = Dict{String,Any}[]
    for t in (complex(0.0,min_imag),complex(real_extent,min_imag),
              complex(-real_extent,min_imag),complex(0.0,1.0e6),
              complex(real_extent/2,3min_imag))
        exact = hat_exact_big(t)
        math_error, bound = setprecision(256) do
            tb = Complex{BigFloat}(t)
            Float64(abs(hat_series(tb)-exact)), Float64(hat_remainder(tb))
        end
        # Separate the finite-series truncation from Float64 evaluation roundoff.
        float_error = Float64(abs(Complex{BigFloat}(hat_series(t))-exact))
        math_error <= bound*(1+1e-12) || error("BigFloat hat remainder control failed")
        float_error <= bound + 64eps(Float64)*abs(hat_series(t)) ||
            error("Float64 hat coefficient control failed")
        push!(records,Dict("t_real"=>real(t),"t_imag"=>imag(t),
            "coefficient_real"=>real(hat_series(t)),"coefficient_imag"=>imag(hat_series(t)),
            "BigFloat_exact_real_decimal"=>string(real(exact)),
            "BigFloat_exact_imag_decimal"=>string(imag(exact)),
            "series_truncation_absolute"=>math_error,"analytic_remainder_bound"=>bound,
            "Float64_absolute_error"=>float_error))
    end
    return records
end

# Real boundary coefficient of the same hat. Small lags use independent
# BigFloat exact arithmetic; large lags use a rapidly convergent moment sum.
function real_hat_lag(lag::Int)
    lag==0 && return 0.0
    lag>0 || return -real_hat_lag(-lag)
    if lag < 8
        return setprecision(256) do
            F(x)=iszero(x) ? big(0.0) : x*log(x)
            t=BigFloat(lag)
            Float64((F(t+1)-2F(t)+F(t-1))/(2big(π)))
        end
    end
    u=inv(Float64(lag)); term=u; value=0.0
    for m in 0:12
        value += 2term/((2m+1)*(2m+2))
        term *= u*u
    end
    return value/(2π)
end

function ghost_strip_witness(H,La,Lb,Gfirst,Glast,a,b,h,nb,periods,mu,kBT)
    count=2periods+1; M=asinh(1.0)/π
    BL=-La+M*Gfirst; BR=Lb+M*Glast
    DL=copy(H); DR=-H
    for cell in 1:count
        r=(cell-1)*nb+1:cell*nb
        DL[r,r] .-= a*I+BL
        DR[r,r] .+= b*I-BR
    end
    dL=eigmin(Hermitian(herm(DL))); dR=eigmin(Hermitian(herm(DR)))
    record=Dict{String,Any}("dL_full_dimensionless"=>dL,"dR_full_dimensionless"=>dR,
        "Gamma_first_Frobenius"=>norm(Gfirst),"Gamma_last_Frobenius"=>norm(Glast),
        "status"=>"unavailable_nonpositive_margin",
        "qualification"=>"floating witness conditional on exact PSD/input/error assumptions; not interval certificate")
    if dL>0 && dR>0
        left=h*norm(Gfirst)/(4π*dL^2); right=h*norm(Glast)/(4π*dR^2)
        record["status"]="conditional_floating_witness"
        record["left_mass_operator_upper_witness"]=left
        record["right_mass_operator_upper_witness"]=right
        record["mass_normalized_Frobenius_upper_witness"]=left+right
        record["FD_density_operator_upper_witness"]=fermi(a-h-mu,kBT)*left+fermi(b-mu,kBT)*right
    end
    return (;record,La,Lb,BL,BR)
end

function chain_hamiltonian(h, T, periods)
    nb = size(h,1); count = 2periods+1
    H = zeros(ComplexF64,count*nb,count*nb)
    for cell in 1:count
        r = (cell-1)*nb+1:cell*nb
        H[r,r] .= h
        if cell < count
            q = cell*nb+1:(cell+1)*nb
            H[r,q] .= T; H[q,r] .= T'
        end
    end
    return H
end
function central_selector(nb, periods)
    r = periods*nb+1:(periods+1)*nb
    B = zeros(ComplexF64,(2periods+1)*nb,nb)
    B[r,:] .= Matrix{ComplexF64}(I,nb,nb)
    return r,B
end

# No fitting or number normalization: source μ is a fixed input.
fermi(x,kBT) = x >= 0 ? exp(-x/kBT)/(1+exp(-x/kBT)) : inv(1+exp(x/kBT))
function odd_tails(checkpoints, beta)
    setprecision(256) do
        b = BigFloat(beta); s2 = big(0.0); s4 = big(0.0)
        result = Dict{Int,Tuple{Float64,Float64}}()
        for n in 0:maximum(checkpoints)-1
            q = BigFloat(2n+1); s2 += inv(q^2); s4 += inv(q^4)
            if n+1 in checkpoints
                result[n+1] = (Float64((b/big(π))^2*(big(π)^2/8-s2)),
                               Float64((b/big(π))^4*(big(π)^4/96-s4)))
            end
        end
        return result
    end
end
function moment_blocks(H,h,S0,S1,central,mu)
    nb = size(h,1); eye = Matrix{ComplexF64}(I,nb,nb)
    H2 = H*H; H3 = H2*H
    M1 = Matrix(H[central,central])
    M2 = Matrix(H2[central,central])+S0
    M3 = Matrix(H3[central,central])+h*S0+S0*h+S1
    return (;M1,M2,M3,a1=M1-mu*eye,a3=M3-3mu*M2+3mu^2*M1-mu^3*eye)
end
function density_estimates(acc,mom,beta,tails)
    nb = size(acc,1); raw = Matrix{ComplexF64}(I,nb,nb)/2 + acc
    m1 = raw-(2/beta)*tails[1]*mom.a1
    m3 = m1+(2/beta)*tails[2]*mom.a3
    return (;raw,m1,m3)
end

function bare_control(H,central,B,mu,kBT,checkpoints,tails,checkstage)
    beta = inv(kBT); nb = length(central)
    h = Matrix(H[central,central]); zero = zeros(ComplexF64,nb,nb)
    mom = moment_blocks(H,h,zero,zero,central,mu)
    ev = eigen(Hermitian(H))
    Vc = ev.vectors[central,:]
    exact = Vc*Diagonal(fermi.(ev.values.-mu,kBT))*Vc'
    acc = zeros(ComplexF64,nb,nb); records = Dict{String,Any}()
    for n in 0:maximum(checkpoints)-1
        z = mu+im*(2n+1)*π/beta
        Gc = (lu(z*I-H)\B)[central,:]
        acc .+= (2/beta).*herm(Gc)
        if n+1 in checkpoints
            checkstage("bare_control_N$(n+1)")
            rho = density_estimates(acc,mom,beta,tails[n+1])
            records[string(n+1)] = Dict("M1_absolute_Frobenius"=>norm(rho.m1-exact),
                "M3_absolute_Frobenius"=>norm(rho.m3-exact),
                "exact_scale_Frobenius"=>norm(exact))
        end
    end
    return (;records,exact,maximum_imaginary_entry=maximum(abs,imag.(exact)),
        final_error=records[string(maximum(checkpoints))]["M3_absolute_Frobenius"])
end

# Terminal-source adapter only; the numerical backend above/below is unchanged.
# Hashes describe Julia logical arrays in (E,k,a,b) order, not HDF5 physical axes.
function checked_array(x,c,name)
    all(isfinite,x) || error("nonfinite array: $name")
    sha_array(x)==c["array_sha256"][name] || error("logical array hash differs: $name")
    return x
end
require_finite_metadata(values...) = all(isfinite,values) || error("nonfinite physical metadata or binding intermediate")
function require_all_PD_certificate(c,base,paths)
    c["execution_authorized"]===true || error("actual Matsubara execution is not authorized")
    certpath=guarded_file(base,c,"pd_certificate"); proofpath=guarded_file(base,c,"pd_certificate_input")
    cert=TOML.parsefile(certpath); proof=TOML.parsefile(proofpath)
    shape=[c["expected_NE"],c["expected_Nk"],c["expected_Nb"],c["expected_Nb"]]; count=shape[1]*shape[2]
    count==337744 && shape==[30704,11,5,5] || error("wrong certificate domain")
    cert["execution"]=="complete" && cert["mode"]=="all" && cert["all_blocks_strict_PD_certified"]===true || error("no completed ALL PD certificate")
    !haskey(cert,"failure") && cert["expected_all_blocks"]==cert["requested_blocks"]==count || error("partial/failed certificate")
    r=cert["results"]; r["checked_blocks"]==r["counts"]["pd"]==count || error("PD count")
    all(r["counts"][k]==0 for k in ("nonpositive","nonhermitian","nonfinite")) && isempty(r["failure_witnesses_first8"]) || error("PD failures")
    r["checked_exact_divrem_components"]==60count || error("inconsistent exact divrem count")
    cert["input_sha256"]==c["pd_certificate_input_sha256"] && proof["mode"]=="all" && proof["shape"]==shape && proof["source_call"]==64 || error("proof input domain")
    proof["source_identity"]==TOML.parsefile(paths["source_input"])["identity"] || error("proof source identity")
    cert["script_sha256"]==c["pd_runner_sha256"] && cert["math_sha256"]==c["pd_math_sha256"] || error("PD implementation identity")
    for (key,source) in (("gamma","checkpoint"),("model","model_h5"),("input","source_input"),("identity","source_identity"),
                         ("operator","source_writer"),("spectral","spectral_writer"),("ledger","source_ledger"),("summary","source_summary"),("map","map_summary"))
        cert["source_sha256"][key]==c[source*"_sha256"]==proof["bindings"][key]["sha256"] || error("certificate source binding: $key")
    end
    return Dict("certificate_sha256"=>c["pd_certificate_sha256"],"checkpoint_sha256"=>c["checkpoint_sha256"],
        "checked_blocks"=>count,"scope"=>"exact stored binary Gamma nodes only; not SCBA convergence or total forward-error certificate")
end
function checkpoint_terminal_equal(paths,shape,checkstage)
    h5open(paths["checkpoint"],"r") do checkpoint
        read(checkpoint["schema"])=="campaign133-gamma-checkpoint-v1" && read(checkpoint["pi_calls_total"])==64 &&
            read(checkpoint["identity"])==TOML.parsefile(paths["source_input"])["identity"] || error("checkpoint identity")
        h5open(paths["terminal_h5"],"r") do terminal
            for suffix in ("_real","_imag")
                name="Gamma"*suffix; size(checkpoint[name])==size(terminal[name])==shape &&
                    eltype(checkpoint[name])==eltype(terminal[name])==Float64 || error("checkpoint axes/type")
                for first in 1:1024:shape[1]
                    checkstage("checkpoint_Gamma_exact_equality"); rows=first:min(first+1023,shape[1])
                    isequal(checkpoint[name][rows,:,:,:],terminal[name][rows,:,:,:]) || error("terminal Gamma differs from certified checkpoint")
                end
            end
        end
    end
    return true
end
function source_frame(paths,c,scfg,mapinfo,checkstage)
    # Owned paths of the CURRENT campaign writers, not a configurable schema.
    ds=Dict("algorithms_toml"=>"algorithms_toml","energy"=>"grid_ε",
        "weight_energy"=>"grid_wᴱ","momentum"=>"grid_κ","weight_momentum"=>"grid_wᵏ",
        "hcell"=>"h_fixedU0")
    basis_paths=Dict("Phi"=>"basis_Φ","chi"=>"basis_χ","H0"=>"basis_H₀",
        "Z"=>"basis_Z","Minv"=>"basis_M⁻¹","Tplus"=>"basis_T₊","Tminus"=>"basis_T₋")
    units=c["units"]
    require_finite_metadata(Float64(units["E0_eV"]),Float64(units["L0_m"]),
        Float64(units["spin_degeneracy"]),Float64(units["target_dimensionless"]),
        Float64(scfg["kmax_per_nm"]),Float64(c["temperature_K"]),Float64(c["voltage_mV"]),
        Float64(c["PI_origin_dimensionless"]),Float64(c["PI_spacing_dimensionless"]),
        Float64(c["mu_dimensionless"]))
    c["logical_axes"]==["energy","momentum","orbital_row","orbital_column"] || error("unsupported logical axes")
    ne,nk,nb=c["expected_NE"],c["expected_Nk"],c["expected_Nb"]
    frame=h5open(paths["model_h5"],"r") do file
        read(file["schema"])=="campaign133-gamma-model-v2" || error("model schema differs")
        algorithms=TOML.parse(read(file[ds["algorithms_toml"]]))
        algorithms["hilbert"]=="product_integration" || error("model is not PI")
        algorithms["embedding"]=="finite_chain" && algorithms["embedding_periods"]==4 || error("model boundary differs")
        e=checked_array(read(file[ds["energy"]]),c,"energy")
        we=checked_array(read(file[ds["weight_energy"]]),c,"weight_energy")
        k=checked_array(read(file[ds["momentum"]]),c,"momentum")
        wk=checked_array(read(file[ds["weight_momentum"]]),c,"weight_momentum")
        h_raw=read_complex(file,ds["hcell"])
        all(isfinite,h_raw) && sha_array(h_raw)==c["h_raw_sha256"] || error("model h_raw binding")
        basis=Dict(name=>checked_array(read_complex(file,path),c,name) for (name,path) in basis_paths)
        return (;e,we,k,wk,h_raw,basis,E0=Float64(units["E0_eV"]),L0=Float64(units["L0_m"]),gs=Float64(units["spin_degeneracy"]))
    end
    h=h5open(paths["spectral_selected"],"r") do file
        read(file["schema"])=="campaign133-stableA-spectral-selected-v1" || error("saved h_H schema")
        isequal(read_complex(file,"h_raw"),frame.h_raw) || error("spectral h_raw differs from model")
        isequal(read(file["energy"]),frame.e) && isequal(read(file["weight_energy"]),frame.we) && isequal(read(file["weight_k"]),frame.wk) || error("spectral axes")
        read(file["energy_scale_eV"])==frame.E0 && isequal(read_complex(file,"Tplus"),frame.basis["Tplus"]) || error("spectral units/coupling")
        saved=checked_array(read_complex(file,"h_H"),c,"hcell")
        size(saved)==(nk,nb,nb) && eltype(saved)==ComplexF64 || error("h_H shape/type")
        all(saved[k,a,b]==conj(saved[k,b,a]) for k in 1:nk,a in 1:nb,b in 1:nb) || error("saved h_H is not EXACT Hermitian")
        TOML.parse(read(file["ledger_toml"]))==TOML.parsefile(paths["source_ledger"]) || error("saved spectral ledger")
        saved
    end
    frame=(;frame...,h)
    length(frame.e)==length(frame.we)==ne && length(frame.k)==length(frame.wk)==nk || error("model axes differ")
    size(frame.h)==(nk,nb,nb) || error("model hcell shape differs")
    all(size(frame.basis[name])==(nb,nb) for name in ("H0","Z","Minv","Tplus","Tminus")) || error("basis operator shape differs")
    frame.E0>0 && frame.L0>0 && frame.gs>0 || error("invalid units")
    Float64(scfg["kmax_per_nm"])>0 || error("source kmax must be positive for L0 binding")
    expected_kmax=Float64(scfg["kmax_per_nm"])*1e9*frame.L0
    kmax_defect=abs(maximum(frame.k)-expected_kmax)
    kmax_tolerance=128eps(Float64)*max(1.0,abs(expected_kmax))
    require_finite_metadata(expected_kmax,kmax_defect,kmax_tolerance)
    kmax_defect<=kmax_tolerance && minimum(frame.k)>=0 || error("L0/source kmax/model kappa mismatch")
    all(>(0),frame.we) && all(>=(0),frame.wk) || error("invalid integration measures")
    all(frame.basis["Tminus"][a,b]==conj(frame.basis["Tplus"][b,a]) for a in 1:nb,b in 1:nb) || error("Tminus not EXACT adjoint")
    norm(frame.basis["H0"]-frame.basis["H0"]')<1e-12 || error("H0 not Hermitian")
    a=Float64(c["PI_origin_dimensionless"]); step=Float64(c["PI_spacing_dimensionless"])
    step>0 || error("nonpositive PI spacing")
    canonical=[a+(j-1)*step for j in 1:ne]
    geometry_error=maximum(abs,frame.e-canonical)
    require_finite_metadata(geometry_error,maximum(abs,canonical))
    geometry_error<=64eps(Float64)*max(1.0,maximum(abs,frame.e)) || error("nonuniform model grid")
    mu=h5open(paths["terminal_h5"],"r") do file
        read(file["schema"])=="campaign133-gamma-reduced-terminal-v1" || error("terminal schema differs")
        Float64(read(file["mu"]))
    end
    isfinite(mu) && mu==c["mu_dimensionless"] || error("stored mu freeze differs")
    checkstage("load_frozen_aggregate_local_Gamma")
    checkpoint_terminal_equal(paths,(ne,nk,nb,nb),checkstage)
    gamma=h5open(paths["terminal_h5"],"r") do file
        read(file["schema"])=="campaign133-gamma-reduced-terminal-v1" || error("terminal schema differs")
        expected=(ne,nk,nb,nb)
        size(file["Gamma_real"])==size(file["Gamma_imag"])==expected || error("Gamma logical shape differs")
        checked_array(read_complex(file,"Gamma"),c,"Gamma")
    end
    # Loader comparator only: stream one k slice of saved A through the SAME
    # original positive real-axis measure, at stored mu; no fit or new Green.
    realrho=zeros(ComplexF64,nk,nb,nb)
    kBT=(1.380649e-23/1.602176634e-19)*c["temperature_K"]/frame.E0
    require_finite_metadata(kBT)
    h5open(paths["terminal_h5"],"r") do file
        read(file["schema"])=="campaign133-gamma-reduced-terminal-v1" || error("terminal schema differs")
        group=file["green"]; name="A"
        size(group[name*"_real"])==size(group[name*"_imag"])==size(gamma) || error("reference A shape differs")
        for k in 1:nk
            A=read_k(group,name,k)
            sha_array(A)==c["reference_A_k_sha256"][string(k)] || error("reference A slice hash differs")
            for e in 1:ne,b in 1:nb,a in 1:nb
                realrho[k,a,b]+=(frame.we[e]/(2π))*fermi(frame.e[e]-mu,kBT)*A[e,a,b]
            end
        end
    end
    saved_number=frame.gs*sum(frame.wk[k]*real(tr(Matrix(view(realrho,k,:,:)))) for k in 1:nk)
    map_number=Float64(mapinfo["mu"]["number"])
    require_finite_metadata(saved_number,map_number)
    isfinite(map_number) && map_number>0 && Float64(units["target_dimensionless"])>0 || error("invalid source number/target")
    number_scale=max(abs(map_number),abs(Float64(units["target_dimensionless"])),eps(Float64))
    number_tolerance=2e-11*number_scale
    require_finite_metadata(number_scale,number_tolerance,saved_number-map_number)
    abs(saved_number-map_number)<=number_tolerance || error("gs/saved A/map number mismatch")
    binding_witness=Dict("source_kmax_per_nm"=>scfg["kmax_per_nm"],"model_kappa_max"=>maximum(frame.k),
        "L0_kmax_absolute_defect"=>kmax_defect,"L0_kmax_absolute_tolerance"=>kmax_tolerance,
        "saved_A_fA_number"=>saved_number,"map_mu_number"=>map_number,
        "number_absolute_defect"=>abs(saved_number-map_number),"number_comparison_scale"=>number_scale,
        "number_absolute_tolerance"=>number_tolerance,"number_scope"=>"unit/input binding, no fit and no quadrature acceptance")
    refs=Dict("stored_same_Gamma_real_PI_fA"=>realrho)
    return (;frame...,refs,gamma,a,step,canonical,geometry_error,mu,binding_witness)
end

function primary_norm(delta,reference)
    den = norm(reference)
    d = Dict{String,Any}("absolute_Frobenius"=>norm(delta),"reference_Frobenius"=>den)
    if den > 0; d["relative_Frobenius"]=norm(delta)/den
    else; d["relative_status"]="unavailable_zero_denominator"; end
    return d
end
sheet(rho,frame) = frame.gs*sum(frame.wk[k]*real(tr(Matrix(view(rho,k,:,:)))) for k in eachindex(frame.k))

function main()
    length(ARGS)==2 || error("usage: matsubara_oracle.jl INPUT.toml NEW_OUTPUT_DIRECTORY")
    config,output = abspath.(ARGS); c=TOML.parsefile(config); base=dirname(config)
    ispath(output) && error("output already exists")
    mkpath(output); cp(config,joinpath(output,"input.toml")); cp(@__FILE__,joinpath(output,"matsubara_oracle.jl"))
    started=time(); stage=Ref("source_guards")
    summary=Dict{String,Any}("status"=>"running","scope"=>"fixed frozen-PI off-axis density, no SCBA/fit",
        "started_utc"=>string(Dates.now(Dates.UTC)),"script_sha256"=>sha_file(@__FILE__),
        "input_sha256"=>sha_file(config),"completed_checkpoints"=>Int[])
    checkstage(name) = begin
        stage[]=name
        time()-started <= c["time_limit_seconds"] || error("phase/chunk wall-time stop")
        bytes=sum(filesize(joinpath(r,f)) for (r,_,fs) in walkdir(output) for f in fs;init=0)
        bytes <= c["maximum_output_bytes"] || error("output cap exceeded")
        println("stage=",name," elapsed_s=",time()-started); flush(stdout)
    end
    try
        c["execution_authorized"]===true || error("V4 actual execution remains on HOLD")
        Threads.nthreads()==c["Julia_threads"] || error("Julia thread count differs")
        BLAS.set_num_threads(c["BLAS_threads"])
        project=Base.active_project(); manifest=joinpath(dirname(project),"Manifest.toml")
        sha_file(project)==c["project_sha256"] || error("active project differs")
        sha_file(manifest)==c["manifest_sha256"] || error("active Manifest differs")
        paths=Dict(key=>guarded_file(base,c,key) for key in
            ("terminal_h5","model_h5","source_input","source_identity","source_writer","map_summary","numerical_parent",
             "checkpoint","spectral_selected","spectral_writer","source_ledger","source_summary"))
        ident=TOML.parsefile(paths["source_identity"]); scfg=TOML.parsefile(paths["source_input"])
        mapinfo=TOML.parsefile(paths["map_summary"])
        # These exact writers own model-v2 and terminal-v1 in one attempt dir.
        owner_dir=dirname(realpath(paths["terminal_h5"]))
        all(dirname(realpath(paths[key]))==owner_dir for key in
            ("model_h5","source_input","source_identity","source_writer","map_summary","checkpoint","spectral_selected","spectral_writer","source_ledger","source_summary")) || error("source artifacts do not share the writer attempt directory")
        ident["source_sha256"]["gamma_functions.jl"]==c["source_writer_sha256"] || error("source writer/identity mismatch")
        ident["source_sha256"]["stable_spectral.jl"]==c["spectral_writer_sha256"] || error("spectral writer/identity mismatch")
        basename(paths["source_writer"])=="gamma_functions.jl" || error("unsupported source writer")
        ident["core_head"]==c["source_core_head"] && ident["core_tree"]==c["source_core_tree"] || error("source Core identity differs")
        ident["input_sha256"]==c["source_input_sha256"] || error("source identity/input mismatch")
        isempty(ident["core_dirty"]) || error("source Core was dirty")
        scfg["temperature_K"]==c["temperature_K"]==70 && c["voltage_mV"]==0 || error("bath/bias differs")
        all(scfg[key]==c[dest] for (key,dest) in (("NE","expected_NE"),("Nk","expected_Nk"),("Nb","expected_Nb"))) || error("source input shape differs")
        c["periods_each_side"]==4 || error("boundary differs")
        c["h_policy"]=="saved_binary_Hermitian_h_H_all9" && scfg["project_AA_trial"]===false || error("wrong h/Gamma policy")
        scfg["spectral_operator"]=="native_covariance_exact_plan_hlocal_all9" && scfg["pi_calls_total"]==mapinfo["call"]==64 || error("not stableA terminal64")
        ledger=TOML.parsefile(paths["source_ledger"]); source_summary=TOML.parsefile(paths["source_summary"])
        ledger==mapinfo["spectral_operator_ledger"] && ledger["operator"]==scfg["spectral_operator"] || error("spectral ledger mismatch")
        source_summary["execution"]=="complete" && source_summary["pi_calls"]==64 && source_summary["reduced"]["last_map"]==mapinfo || error("source terminal map ownership")
        pd_witness=require_all_PD_certificate(c,base,paths)
        mapinfo["absolute_E0_eV"]==c["units"]["E0_eV"] || error("source map energy unit differs")
        mapinfo["mu"]["target"]==c["units"]["target_dimensionless"] || error("source target differs")
        frame=source_frame(paths,c,scfg,mapinfo,checkstage)
        frame.mu==mapinfo["mu"]["mu"] || error("mu and source map are not aligned")
        prior=Dict("target_dimensionless"=>Float64(c["units"]["target_dimensionless"]))
        c["mu_eV"]=frame.mu*frame.E0 # metadata conversion only; solve uses stored mu below.
        ne,nk,nb=size(frame.gamma)[1:3]; periods=c["periods_each_side"]
        central,B=central_selector(nb,periods); checkpoints=Int.(c["checkpoints"])
        issorted(checkpoints) && all(>(0),checkpoints) || error("bad checkpoints")
        # Exact SI defining constants, converted to the source dimensionless energy.
        kB_eV=1.380649e-23/1.602176634e-19
        kBT=kB_eV*c["temperature_K"]/frame.E0; beta=inv(kBT); mu=frame.mu
        require_finite_metadata(kBT,beta,mu,c["mu_eV"])
        tails=odd_tails(checkpoints,beta)
        Hs=[chain_hamiltonian(Matrix(view(frame.h,k,:,:)),frame.basis["Tplus"],periods) for k in 1:nk]
        maximum(norm(H-H') for H in Hs) <= 1e-11 || error("chain H is non-Hermitian")
        summary["identity"]=Dict("source_sha256"=>c["terminal_h5_sha256"],"gamma_sha256"=>sha_array(frame.gamma),
            "hcell_sha256"=>sha_array(frame.h),"model_h5_sha256"=>c["model_h5_sha256"],
            "terminal_h5_sha256"=>c["terminal_h5_sha256"],"source_core_head"=>c["source_core_head"],
            "source_core_tree"=>c["source_core_tree"],"source_identity_sha256"=>c["source_identity_sha256"],
            "map_summary_sha256"=>c["map_summary_sha256"],"source_writer_sha256"=>c["source_writer_sha256"],
            "owned_terminal_paths"=>["Gamma_real","Gamma_imag","mu","green/A_real","green/A_imag"],
            "logical_axes"=>c["logical_axes"],"basis_and_grid_array_sha256"=>c["array_sha256"],
            "project_sha256"=>sha_file(project),"manifest_sha256"=>sha_file(manifest),
            "julia_version"=>string(VERSION),"Julia_threads"=>Threads.nthreads(),"BLAS_threads"=>BLAS.get_num_threads(),
            "E0_eV"=>frame.E0,"L0_m"=>frame.L0,"spin_degeneracy"=>frame.gs,
            "mu_eV"=>c["mu_eV"],"mu_dimensionless"=>mu,"kBT_dimensionless"=>kBT,
            "origin"=>frame.a,"spacing"=>frame.step,"canonical_geometry_defect"=>frame.geometry_error,
            "source_families"=>c["local_names"],"excluded"=>["embedding","embedding_plus","embedding_minus"])
        summary["source_binding_witness"]=frame.binding_witness
        summary["global_exact_PD_binding"]=pd_witness
        summary["h_policy"]="actual saved EXACT Hermitian h_H, all nine cells; raw h crosschecked only"
        checkstage("independent_kernel_and_bare_controls")
        min_imag=π/beta/frame.step
        require_finite_metadata(min_imag)
        summary["kernel_controls"]=kernel_controls(min_imag,maximum(abs,frame.canonical.-mu)/frame.step)
        summary["minimum_abs_t_lower_bound"]=min_imag
        # Complex control makes using real.(G) instead of Herm(G) observable.
        Hcomplex=Matrix(Diagonal(complex.(range(-0.4,0.4,length=size(B,1)))))
        for j in 1:size(B,1)-1; Hcomplex[j,j+1]=0.03+0.02im; Hcomplex[j+1,j]=0.03-0.02im; end
        ctrl=bare_control(Hcomplex,central,B,mu,kBT,checkpoints,tails,checkstage)
        scalar_records=Dict{String,Any}()
        for energy in (-0.4,mu,0.4)
            scalar=bare_control(reshape(ComplexF64[energy],1,1),1:1,ones(ComplexF64,1,1),
                mu,kBT,checkpoints,tails,checkstage)
            scalar.final_error<=c["control_absolute_Frobenius_tolerance"] || error("scalar Fermi control failed")
            scalar_records[string(energy)]=scalar.records
        end
        ctrl.final_error<=c["control_absolute_Frobenius_tolerance"] || error("complex bare-chain control failed")
        ctrl.maximum_imaginary_entry>1e-5 || error("complex control lacks imaginary offdiagonal signal")
        bare_source=bare_control(Hs[1],central,B,mu,kBT,checkpoints,tails,checkstage)
        bare_source.final_error<=c["control_absolute_Frobenius_tolerance"] || error("source bare-chain control failed")
        summary["bare_complex_controls"]=ctrl.records; summary["bare_source_k1_controls"]=bare_source.records
        summary["scalar_fermi_controls"]=scalar_records
        summary["bare_complex_exact_imaginary_entry_maximum"]=ctrl.maximum_imaginary_entry
        checkstage("selfenergy_moments")
        moments=Vector{Any}(undef,nk); gamma_norm_sums=zeros(nk)
        ghost=Vector{Any}(undef,nk); endpoint_coeff=[real_hat_lag(j) for j in 0:ne-1]
        for k in 1:nk
            S0=zeros(ComplexF64,nb,nb); S1=similar(S0); fill!(S1,0)
            La=zeros(ComplexF64,nb,nb); Lb=similar(La); fill!(Lb,0)
            for j in 1:ne
                block=Matrix(view(frame.gamma,j,k,:,:)); gamma_norm_sums[k]+=norm(block)
                S0 .+= (frame.step/(2π)).*block
                S1 .+= (frame.step*frame.canonical[j]/(2π)).*block
                La .-= endpoint_coeff[j].*block
                Lb .+= endpoint_coeff[ne+1-j].*block
            end
            moments[k]=merge(moment_blocks(Hs[k],Matrix(view(frame.h,k,:,:)),S0,S1,central,mu),(;S0,S1))
            ghost[k]=ghost_strip_witness(Hs[k],La,Lb,Matrix(view(frame.gamma,1,k,:,:)),
                Matrix(view(frame.gamma,ne,k,:,:)),frame.a,last(frame.canonical),frame.step,nb,periods,mu,kBT)
        end
        summary["ghost_strip_witnesses"]=Dict("k$(k)"=>ghost[k].record for k in 1:nk)
        small=joinpath(output,"density_witnesses.h5")
        h5open(small,"w") do file
            file["schema"]="campaign133-frozen-PI-Matsubara-density-v1"
            file["energy_dimensionless"]=frame.e; file["k_dimensionless"]=frame.k
            file["weight_k_dimensionless"]=frame.wk
            write_complex(file,"hcell",frame.h); write_complex(file,"Tplus",frame.basis["Tplus"])
            write_complex(file,"bare_complex_exact_density",ctrl.exact)
            write_complex(file,"bare_source_k1_exact_density",bare_source.exact)
            mg=create_group(file,"moments")
            for k in 1:nk
                kg=create_group(mg,"k$(k)")
                for name in (:S0,:S1,:M1,:M2,:M3,:a1,:a3)
                    write_complex(kg,string(name),getproperty(moments[k],name))
                end
                for name in (:La,:Lb,:BL,:BR)
                    write_complex(kg,string(name),getproperty(ghost[k],name))
                end
            end
            for label in keys(frame.refs); write_complex(create_group(file,label),"rho",frame.refs[label]); end
            create_group(file,"checkpoints"); create_group(file,"resolvent_witnesses")
        end
        chunk=c["frequency_chunk"]; W=zeros(ComplexF64,ne,chunk)
        sigma=zeros(ComplexF64,chunk,nk*nb^2); Gammaflat=reshape(frame.gamma,ne,:)
        acc=zeros(ComplexF64,nk,nb,nb); omitted_sigma_bounds=zeros(maximum(checkpoints),nk)
        max_backward=0.0; max_herm_defect=0.0; max_sigma_noncausal=0.0
        previous=nothing
        for firstn in 0:chunk:maximum(checkpoints)-1
            checkstage("Matsubara_chunk_n$(firstn)")
            nfreq=min(chunk,maximum(checkpoints)-firstn)
            for q in 1:nfreq,j in 1:ne
                z=mu+im*(2(firstn+q-1)+1)*π/beta
                W[j,q]=hat_series((z-frame.canonical[j])/frame.step)
            end
            mul!(view(sigma,1:nfreq,:),transpose(view(W,:,1:nfreq)),Gammaflat)
            for q in 1:nfreq
                n=firstn+q-1; omega=(2n+1)*π/beta; z=mu+im*omega
                for k in 1:nk
                    S=zeros(ComplexF64,nb,nb)
                    for b in 1:nb,a in 1:nb; S[a,b]=sigma[q,k+(a-1)*nk+(b-1)*nk*nb]; end
                    D=z*I-Hs[k]
                    for cell in 1:2periods+1
                        r=(cell-1)*nb+1:cell*nb; D[r,r] .-= S
                    end
                    X=lu(D)\B; G=Matrix(X[central,:]); Hpart=herm(G)
                    acc[k,:,:] .+= (2/beta).*Hpart
                    backward=norm(D*X-B)/(norm(D)*norm(X)+norm(B))
                    max_backward=max(max_backward,backward)
                    max_herm_defect=max(max_herm_defect,norm(Hpart-Hpart'))
                    omitted_sigma_bounds[n+1,k]=hat_remainder(complex(0.0,omega/frame.step))*gamma_norm_sums[k]
                    if n in (0,127,255,511,1023)
                        imS=(S-S')/(2im)
                        max_sigma_noncausal=max(max_sigma_noncausal,eigmax(Hermitian(herm(imS))))
                        h5open(small,"r+") do file
                            w=create_group(file["resolvent_witnesses"],"n$(n)_k$(k)")
                            w["omega_dimensionless"]=omega; w["backward_error"]=backward
                            write_complex(w,"Sigma",S); write_complex(w,"Gc",G)
                            w["omitted_Sigma_Frobenius_bound"]=omitted_sigma_bounds[n+1,k]
                            w["zSigma_minus_S0_Frobenius"]=norm(z*S-moments[k].S0)
                            w["z2Sigma_minus_zS0_minus_S1_Frobenius"]=norm(z^2*S-z*moments[k].S0-moments[k].S1)
                            w["HermG_minus_M3_asymptotic_Frobenius"]=norm(Hpart+moments[k].a1/omega^2-moments[k].a3/omega^4)
                        end
                    end
                end
                if n+1 in checkpoints
                    nn=n+1; checkstage("save_density_N$(nn)")
                    raw=similar(acc); r1=similar(acc); r3=similar(acc)
                    for k in 1:nk
                        r=density_estimates(Matrix(view(acc,k,:,:)),moments[k],beta,tails[nn])
                        raw[k,:,:].=r.raw; r1[k,:,:].=r.m1; r3[k,:,:].=r.m3
                    end
                    record=Dict{String,Any}("N"=>nn,"tail_T2"=>tails[nn][1],"tail_T4"=>tails[nn][2],
                        "M1_minus_M3_absolute_Frobenius"=>norm(r1-r3),
                        "M1_N_over_target"=>sheet(r1,frame)/prior["target_dimensionless"],
                        "M3_N_over_target"=>sheet(r3,frame)/prior["target_dimensionless"])
                    for (label,ref) in frame.refs
                        d=r3-ref
                        record["M3_minus_"*label]=merge(primary_norm(d,ref),Dict(
                            "per_k_absolute_Frobenius"=>[norm(Matrix(view(d,k,:,:))) for k in 1:nk],
                            "signed_delta_N_over_target"=>(sheet(r3,frame)-sheet(ref,frame))/prior["target_dimensionless"]))
                    end
                    if previous!==nothing
                        record["M3_minus_previous"]=merge(primary_norm(r3-previous,previous),Dict(
                            "per_k_absolute_Frobenius"=>[norm(Matrix(view(r3-previous,k,:,:))) for k in 1:nk],
                            "signed_delta_N_over_target"=>(sheet(r3,frame)-sheet(previous,frame))/prior["target_dimensionless"]))
                    end
                    h5open(small,"r+") do file
                        group=create_group(file["checkpoints"],"N$(nn)")
                        write_complex(group,"raw_partial_density",raw); write_complex(group,"rho_M1",r1)
                        write_complex(group,"rho_M3",r3)
                    end
                    summary["N$(nn)"]=record; push!(summary["completed_checkpoints"],nn)
                    previous=copy(r3)
                    summary["elapsed_seconds"]=time()-started
                    open(joinpath(output,"summary.toml"),"w") do io; TOML.print(io,summary); end
                end
            end
        end
        h5open(small,"r+") do file; file["kernel_series_omitted_Sigma_Frobenius_bounds"]=omitted_sigma_bounds; end
        summary["maximum_selected_solve_backward_error"]=max_backward
        summary["maximum_accumulated_Herm_part_defect"]=max_herm_defect
        summary["maximum_selected_Sigma_positive_imag_eigenvalue"]=max_sigma_noncausal
        summary["tail_status"]="M1/M3 and N refinement witnesses; no formal total tail error bound"
        summary["series_status"]="analytic omitted-series Sigma norm bound excludes Float64 summation/evaluation errors"
        summary["SCBA_convergence"]="not_measured"
        summary["atom_and_ghost_attribution"]="not established by density difference alone"
        summary["status"]="completed"
    catch err
        summary["status"]="failed"; summary["error"]=sprint(showerror,err)
        rethrow()
    finally
        summary["current_stage"]=stage[]; summary["elapsed_seconds"]=time()-started
        summary["finish_utc"]=string(Dates.now(Dates.UTC))
        open(joinpath(output,"summary.toml"),"w") do io; TOML.print(io,summary); end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__; main(); end
