using QCLNEGF, Unitful, LinearAlgebra, HDF5, TOML, SHA, Dates

function main()
    BLAS.set_num_threads(1)
    out = @__DIR__
    state = joinpath(out,"../raw/S01-fixedU32/state.h5")
    input = joinpath(out,"source-input.toml")
    started = time()
    cfg = TOML.parsefile(input)
    srcsha = bytes2hex(open(sha256,state))
    srcsha == "155f8d296715f8e8a4ca951a0ed07aa392742112bd8b3b934beea54304374a2a" || error("source SHA mismatch")
    p=reference_parameters(V_period=cfg["voltage_mV"]*u"mV",Tᴸ=cfg["temperature_K"]*u"K",Tᴸᴼ=cfg["temperature_K"]*u"K",ΔV_alloy=.60u"eV",Ω₀=4.5168e-29u"m^3")
    n=NumericalParameters(N_z=cfg["Nz"],N_b=cfg["Nb"],P_basis=cfg["basis_periods"],E_min=cfg["energy_min_eV"]*u"eV",E_max=cfg["energy_max_eV"]*u"eV",N_E=cfg["NE"],M_E=cfg["trust_margin_meV"]*u"meV",k_max=cfg["kmax_per_nm"]*u"nm^-1",N_k=cfg["Nk"],N_φ=cfg["Nphi"],qz_max=cfg["qzmax_per_nm"]*u"nm^-1",N_qz=cfg["Nqz"],η_seed=cfg["seed_meV"]*u"meV")
    s=ScaleSystem(); g=build_grids(p,n,s)
    profiles=build_profiles(p,n,s,g)
    basis=build_basis(p,n,s,g,profiles;localization=:pzp)
    h=project_hamiltonians(basis,g,profiles,p,s,zeros(n.N_z))
    target=QCLNEGF._scaled_physics(p,s).Nᴰ²ᴰ
    mu_eV=8.519932576353691e-3
    mu=mu_eV/s.E₀_eV
    kBT=QCLNEGF._electronvolts(CODATA.kᴮₑᵥ*p.Tᴸ)/s.E₀_eV
    ne=n.N_E; nk=n.N_k; nb=n.N_b
    step=(g.ε[end]-g.ε[1])/(ne-1)
    grids=[g.ε,collect(range(g.ε[1],g.ε[end];length=2ne-1)),[g.ε[1]+(j-.5)*step for j=1:ne-1]]
    weights=[g.wᴱ,fill(step/2,2ne-1),fill(step,ne-1)]
    weights[2][[1,end]]./=2
    labels=["original_trapezoid_30721","refined_trapezoid_61441","midpoint_phase_30720"]
    mats=[zeros(ComplexF64,nk,nb,nb) for _=1:3]
    fmats=[zeros(ComplexF64,nk,nb,nb) for _=1:3]
    storedmat=zeros(ComplexF64,nk,nb,nb)
    storedfmat=similar(storedmat); fill!(storedfmat,0)
    alignment=zeros(nk,4)
    families=String[]
    D=zeros(ComplexF64,nb,nb); A=similar(D)
    work=QCLNEGF.QCLNumerics._ProductionInverseWorkspace(nb)
    println("start_utc=",Dates.now(Dates.UTC)," source_sha=",srcsha," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads());flush(stdout)
    h5open(state,"r") do file
        @assert read(file["energy"])==g.ε && read(file["weight_energy"])==g.wᴱ
        @assert read(file["k"])==g.κ && read(file["weight_k"])==g.wᵏ
        @assert all(iszero,read(file["hartree"]))
        append!(families,sort(filter(x->!(x in ["embedding_plus","embedding_minus"]),collect(keys(file["selfenergy"])))))
        @assert length(families)==6 && "embedding" in families
        for m=1:nk
            sigma=zeros(ComplexF64,ne,nb,nb)
            for name in families
                group=file["selfenergy/"*name]
                re=group["retarded_real"][:,m,:,:]
                imv=group["retarded_imag"][:,m,:,:]
                @inbounds for i in eachindex(sigma)
                    sigma[i]+=complex(re[i],imv[i])
                end
            end
            Rstored=complex.(file["green/GR_real"][:,m,:,:],file["green/GR_imag"][:,m,:,:])
            Astored=complex.(file["green/A_real"][:,m,:,:],file["green/A_imag"][:,m,:,:])
            for test=1:3
                energies=grids[test]; w=weights[test]
                for j in eachindex(energies)
                    e,t = test==1 ? (j,0.) : test==2 ? ((j+1)÷2,iseven(j) ? .5 : 0.) : (j,.5)
                    en=energies[j]
                    @inbounds for b=1:nb,a=1:nb
                        sig = iszero(t) ? sigma[e,a,b] : (1-t)*sigma[e,a,b]+t*sigma[e+1,a,b]
                        D[a,b]=(a==b ? en : 0.)-h[m,a,b]-sig
                    end
                    QCLNEGF.QCLNumerics._production_inverse!(D,work)
                    @inbounds for b=1:nb,a=1:nb
                        A[a,b]=im*(D[a,b]-conj(D[b,a]))
                    end
                    fd=QCLNEGF._fermi(en,mu,kBT)
                    @inbounds for b=1:nb,a=1:nb
                        mats[test][m,a,b]+=w[j]/(2π)*A[a,b]
                        fmats[test][m,a,b]+=w[j]/(2π)*fd*A[a,b]
                    end
                    if test==1
                        dr,da,nr,na=0.,0.,0.,0.
                        @inbounds for b=1:nb,a=1:nb
                            dr+=abs2(D[a,b]-Rstored[j,a,b]);nr+=abs2(Rstored[j,a,b])
                            da+=abs2(A[a,b]-Astored[j,a,b]);na+=abs2(Astored[j,a,b])
                            storedmat[m,a,b]+=w[j]/(2π)*Astored[j,a,b]
                            storedfmat[m,a,b]+=w[j]/(2π)*fd*Astored[j,a,b]
                        end
                        alignment[m,1]=max(alignment[m,1],sqrt(dr)/max(sqrt(nr),1e-300))
                        alignment[m,2]=max(alignment[m,2],sqrt(da)/max(sqrt(na),1e-300))
                        alignment[m,3]=max(alignment[m,3],sqrt(dr))
                        alignment[m,4]=max(alignment[m,4],sqrt(da))
                    end
                end
            end
            println("k=",m," alignmentGR=",alignment[m,1]," alignmentA=",alignment[m,2]," elapsed_s=",time()-started);flush(stdout)
        end
    end
    charge(M)=p.g_s*sum(g.wᵏ[m]*sum(real(M[m,a,a]) for a=1:nb) for m=1:nk)
    rsum(M)=maximum(norm(Matrix(view(M,m,:,:))-I)/sqrt(nb) for m=1:nk)
    result=Dict{String,Any}()
    for j=1:3
        result[labels[j]]=Dict("nodes"=>length(grids[j]),"domain_min_eV"=>g.ε[1]*s.E₀_eV,"domain_max_eV"=>g.ε[end]*s.E₀_eV,"node_min_eV"=>grids[j][1]*s.E₀_eV,"node_max_eV"=>grids[j][end]*s.E₀_eV,"weight_sum_eV"=>sum(weights[j])*s.E₀_eV,"r_sum"=>rsum(mats[j]),"N_fA_dimensionless"=>charge(fmats[j]),"N_over_target"=>charge(fmats[j])/target,"integral_relative_difference_from_original"=>norm(mats[j]-mats[1])/norm(mats[1]),"N_relative_difference_from_original"=>charge(fmats[j])/charge(fmats[1])-1)
    end
    summary=Dict{String,Any}("source_state_sha256"=>srcsha,"source_input_sha256"=>bytes2hex(sha256(read(input))),"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"h_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(h)))),"core_head"=>readchomp(`git -C $(dirname(dirname(pathof(QCLNEGF)))) rev-parse HEAD`),"core_source"=>pathof(QCLNEGF),"julia_version"=>string(VERSION),"families_summed_once"=>families,"embedding_excluded"=>["embedding_plus","embedding_minus"],"mu_eV"=>mu_eV,"temperature_K"=>70.,"target_dimensionless"=>target,"E0_eV"=>s.E₀_eV,"max_alignment_GR_relative_Frobenius"=>maximum(alignment[:,1]),"max_alignment_A_relative_Frobenius"=>maximum(alignment[:,2]),"stored_A_r_sum"=>rsum(storedmat),"stored_A_N_over_target"=>charge(storedfmat)/target,"elapsed_s"=>time()-started,"interpolation"=>"piecewise linear complex entries; no extrapolation; fixed total selfenergy, fixed U=0 Hamiltonian; no seed, SCBA update, normalization, pins or guards changed","r_sum_definition"=>"max_k || integral A(E,k) dE/(2pi) - I_Nb ||_F / sqrt(Nb); dimensionless; exact core validation.jl:_spectral_sum_residual","N_definition"=>"g_s sum_k weight_k tr[integral f(E;fixed_mu,70K) A(E,k)dE/(2pi)]; dimensionless sheet density; no normalization","quadrature_note"=>"original/refined trapezoid includes both domain endpoints; midpoint phase uses 30720 full-cell weights over same source domain, nodes half step inward; phase control is not a baseline NEGF grid refinement","scope"=>"frozen-resolvent diagnostic on nonlinear iterate32, unconverged source; not fine SCBA solution, not energy-window accuracy proof","results"=>result)
    h5open(joinpath(out,"integrals.h5"),"w") do f
        f["schema"]="campaign133-frozen-resolvent-v1"
        f["h_real"]=real.(h);f["h_imag"]=imag.(h)
        f["k"]=g.κ;f["weight_k"]=g.wᵏ
        f["alignment_by_k"]=alignment
        for j=1:3
            group=create_group(f,labels[j]);group["energy_dimensionless"]=grids[j];group["weight_energy_dimensionless"]=weights[j]
            group["integral_A_over_2pi_real"]=real.(mats[j]);group["integral_A_over_2pi_imag"]=imag.(mats[j])
            group["integral_fA_over_2pi_real"]=real.(fmats[j]);group["integral_fA_over_2pi_imag"]=imag.(fmats[j])
        end
        f["stored_integral_A_over_2pi_real"]=real.(storedmat);f["stored_integral_A_over_2pi_imag"]=imag.(storedmat)
        f["stored_integral_fA_over_2pi_real"]=real.(storedfmat);f["stored_integral_fA_over_2pi_imag"]=imag.(storedfmat)
    end
    open(joinpath(out,"summary.toml"),"w") do io;TOML.print(io,summary);end
    println(summary)
end
main()
