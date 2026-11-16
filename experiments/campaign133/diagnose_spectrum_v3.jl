using QCLNEGF, Unitful, LinearAlgebra, HDF5, TOML, SHA
function main()
BLAS.set_num_threads(1)
length(ARGS)==3 || error("usage: diagnose_spectrum.jl INPUT STATE NEW_OUTPUT")
input,state,output=abspath.(ARGS)
ispath(output) && error("output exists")
mkpath(output)
c=TOML.parsefile(input)
p=reference_parameters(V_period=c["voltage_mV"]*u"mV",Tᴸ=c["temperature_K"]*u"K",Tᴸᴼ=c["temperature_K"]*u"K",ΔV_alloy=get(c,"alloy",false) ? .60u"eV" : nothing,Ω₀=get(c,"alloy",false) ? 4.5168e-29u"m^3" : nothing)
n=NumericalParameters(N_z=c["Nz"],N_b=c["Nb"],P_basis=c["basis_periods"],E_min=c["energy_min_eV"]*u"eV",E_max=c["energy_max_eV"]*u"eV",N_E=c["NE"],M_E=c["trust_margin_meV"]*u"meV",k_max=c["kmax_per_nm"]*u"nm^-1",N_k=c["Nk"],N_φ=c["Nphi"],qz_max=c["qzmax_per_nm"]*u"nm^-1",N_qz=c["Nqz"],η_seed=c["seed_meV"]*u"meV")
s=ScaleSystem()
g=build_grids(p,n,s)
readcomplex(group,name)=complex.(read(group[name*"_real"]),read(group[name*"_imag"]))
A,Gn,Gp,R=h5open(state,"r") do file
    @assert read(file["energy"])==g.ε && read(file["weight_energy"])==g.wᴱ
    @assert read(file["k"])==g.κ && read(file["weight_k"])==g.wᵏ
    A=readcomplex(file["green"],"A")
    Gn=readcomplex(file["green"],"GL"); Gn .*= -im
    Gp=readcomplex(file["green"],"GG"); Gp .*= im
    R=zeros(ComplexF64,size(A))
    for name in keys(file["selfenergy"])
        name in ("embedding_plus","embedding_minus") && continue
        piece=readcomplex(file["selfenergy/"*name],"retarded")
        R .+= piece
        piece=nothing
        GC.gc()
    end
    (A,Gn,Gp,R)
end
q=[p.g_s/(2π)*g.wᴱ[e]*sum(g.wᵏ[m]*sum(real(A[e,m,a,a]) for a=1:n.N_b) for m=1:n.N_k) for e=1:n.N_E]
target=QCLNEGF._scaled_physics(p,s).Nᴰ²ᴰ
kBT=QCLNEGF._electronvolts(CODATA.kᴮₑᵥ*p.Tᴸ)/s.E₀_eV
fit=QCLNEGF.QCLNumerics._marker_spectral_chemical_potential(g.ε,q,target,kBT)
fit.status===:available || error("mu fit unavailable: $(fit.status)")
f=[QCLNEGF._fermi(e,fit.mu,kBT) for e in g.ε]
less=zeros(ComplexF64,size(R))
greater=similar(R)
for m=1:n.N_k,e=1:n.N_E,b=1:n.N_b,a=1:n.N_b
    gamma=im*(R[e,m,a,b]-conj(R[e,m,b,a]))
    greater[e,m,a,b]=-im*gamma
end
family=SelfEnergyFamily(R,less,greater)
policy=SCBAPhysicsMarkerPolicy(max_spectral_blocks=2048)
open(joinpath(output,"linewidth.csv"),"w") do io
    println(io,"weight,q10_gamma_over_dE,q50_gamma_over_dE,q90_gamma_over_dE,underresolved_weight,blocks,status")
    for label in ("spectral_A","thermal_fA","actual_Gn")
        spectrum = label=="actual_Gn" ? Gn : label=="spectral_A" ? A : copy(A)
        if label=="thermal_fA"
            for m=1:n.N_k,e=1:n.N_E,b=1:n.N_b,a=1:n.N_b
                spectrum[e,m,a,b]*=f[e]
            end
        end
        r=QCLNEGF.QCLNumerics._sampled_linewidth_markers(spectrum,family,g,policy)
        println(io,join((label,r.q10,r.q50,r.q90,r.underresolved,r.blocks,r.status),","))
        spectrum=nothing
        GC.gc()
    end
end
function integrals(stride)
    indices=1:stride:n.N_E
    (n.N_E-1)%stride==0 || error("incommensurate stride")
    step=(g.ε[2]-g.ε[1])*stride
    matrices=zeros(ComplexF64,n.N_k,n.N_b,n.N_b)
    charge,fdcharge=0.0,0.0
    for m=1:n.N_k,(j,e) in enumerate(indices)
        weight=step*(j==1 || j==length(indices) ? .5 : 1.0)
        for b=1:n.N_b,a=1:n.N_b
            matrices[m,a,b]+=weight/(2π)*A[e,m,a,b]
        end
        charge+=p.g_s/(2π)*weight*g.wᵏ[m]*sum(real(Gn[e,m,a,a]) for a=1:n.N_b)
        fdcharge+=p.g_s/(2π)*weight*g.wᵏ[m]*f[e]*sum(real(A[e,m,a,a]) for a=1:n.N_b)
    end
    return matrices,charge,fdcharge
end
fine,chargefine,fdfine=integrals(1)
open(joinpath(output,"quadrature.csv"),"w") do io
    println(io,"stride,nodes,step_eV,max_r_sum,relative_integral_difference_from_fine,actual_charge,actual_charge_relative_to_target,fd_charge,fd_charge_relative_to_target")
    for stride in (1,2,4,8)
        M,N,Nfd=integrals(stride)
        residual=maximum(norm(Matrix(view(M,m,:,:))-I)/sqrt(n.N_b) for m=1:n.N_k)
        println(io,join((stride,length(1:stride:n.N_E),(g.ε[2]-g.ε[1])*stride*s.E₀_eV,residual,norm(M-fine)/norm(fine),N,N/target-1,Nfd,Nfd/target-1),","))
    end
end
delta_n,delta_p,anorm,fnorm,identity=0.0,0.0,0.0,0.0,0.0
for m=1:n.N_k,e=1:n.N_E,b=1:n.N_b,a=1:n.N_b
    w=p.g_s/(2π)*g.wᴱ[e]*g.wᵏ[m]
    z=A[e,m,a,b]
    delta_n+=w*abs2(Gn[e,m,a,b]-f[e]*z)
    delta_p+=w*abs2(Gp[e,m,a,b]-(1-f[e])*z)
    anorm+=w*abs2(z)
    fnorm+=w*abs2(f[e]*z)
    identity+=w*abs2(Gn[e,m,a,b]+Gp[e,m,a,b]-z)
end
summary=Dict("mu_eV"=>fit.mu*s.E₀_eV,"actual_charge"=>chargefine,"fd_charge"=>fdfine,"target_charge"=>target,"fdt_pair_over_A"=>sqrt((delta_n+delta_p)/anorm),"occupied_defect_over_thermal_occupied"=>sqrt(delta_n/fnorm),"paired_spectral_identity_relative"=>sqrt(identity/anorm),"state_sha256"=>bytes2hex(open(sha256,state)),"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"note"=>"Subsampling is quadrature on one fixed saved fine-grid spectrum, not independently solved coarse grids. Linewidths are sampled modal diagnostics, not FWHM or acceptance gates.")
open(joinpath(output,"summary.toml"),"w") do io
    TOML.print(io,summary)
end
cp(@__FILE__,joinpath(output,"diagnose_spectrum.jl"))
println(summary)
end
main()
