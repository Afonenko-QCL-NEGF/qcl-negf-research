# Source-frame reconstruction retained from the reviewed frozen PI control.
# No scientific arrays are read until run_probe.main() passes identity guards.
sha_file(path) = bytes2hex(open(sha256, path))
sha_array(x) = bytes2hex(sha256(reinterpret(UInt8, vec(x))))
read_complex(group,name) = complex.(read(group[name*"_real"]),read(group[name*"_imag"]))
read_k(group,name,k) = complex.(group[name*"_real"][:,k,:,:],group[name*"_imag"][:,k,:,:])
function write_complex(group,name,x)
    group[name*"_real"]=real.(x);group[name*"_imag"]=imag.(x)
end
const LOCAL_NAMES = ["IFR", "LO", "acoustic", "alloy", "impurity"]
function parameters(c)
    p = reference_parameters(V_period=c["voltage_mV"]*u"mV", Tᴸ=c["temperature_K"]*u"K",
        Tᴸᴼ=c["temperature_K"]*u"K", ΔV_alloy=.60u"eV", Ω₀=4.5168e-29u"m^3")
    n = NumericalParameters(N_z=c["Nz"], N_b=c["Nb"], P_basis=c["basis_periods"],
        E_min=c["energy_min_eV"]*u"eV", E_max=c["energy_max_eV"]*u"eV", N_E=c["NE"],
        M_E=c["trust_margin_meV"]*u"meV", k_max=c["kmax_per_nm"]*u"nm^-1", N_k=c["Nk"],
        N_φ=c["Nphi"], qz_max=c["qzmax_per_nm"]*u"nm^-1", N_qz=c["Nqz"], η_seed=c["seed_meV"]*u"meV")
    return p, n, ScaleSystem()
end

function source_frame(file, c)
    p, n, s = parameters(c)
    g = build_grids(p,n,s)
    for (name, actual) in (("energy",g.ε),("weight_energy",g.wᴱ),("k",g.κ),
                           ("weight_k",g.wᵏ),("x",g.x),("weight_x",g.wˣ))
        read(file[name]) == actual || error("source axis differs: $name")
    end
    read(file["schema"]) == "campaign133-inner-raw-v3" || error("source is not v3")
    read(file["energy_scale_eV"]) == s.E₀_eV || error("source E0 differs")
    read(file["length_scale_m"]) == s.L₀_m || error("source L0 differs")
    read(file["spin_degeneracy"]) == p.g_s || error("source spin differs")
    all(iszero,read(file["hartree"])) || error("source Hartree differs")
    Set(keys(file["selfenergy"])) == Set([LOCAL_NAMES;"embedding";"embedding_plus";"embedding_minus"]) || error("source families differ")
    profiles = build_profiles(p,n,s,g)
    rebuilt = build_basis(p,n,s,g,profiles; localization=:pzp)
    mapping = Dict(:Φ=>"Phi",:χ=>"chi",:H₀=>"H0",:Z=>"Z",:M⁻¹=>"Minv",:T₊=>"Tplus",:T₋=>"Tminus")
    arrays = Dict(name=>read_complex(file["basis"], stored) for (name,stored) in mapping)
    comparison = Dict(string(name)=>Dict("source_sha256"=>sha_array(array),
        "rebuild_absolute_Frobenius"=>norm(array-getfield(rebuilt,name)),
        "rebuild_scale_Frobenius"=>norm(array)) for (name,array) in arrays)
    # Numerical source operators are bound exactly to v3. Only unpersisted
    # diagnostic metadata comes from the matching rebuilt basis container.
    basis = BasisData((get(arrays,name,getfield(rebuilt,name)) for name in fieldnames(BasisData))...)
    empty_kernels = KernelSet(Dict{Symbol,Array{ComplexF64,6}}(),Dict{Symbol,Float64}(),zeros(ComplexF64,0,0,0),Symbol[])
    sp = QCLNEGF._scaled_physics(p,s)
    shifts = map(δ->build_shift_plan(g.ε,δ),(sp.Eᵖ,-sp.Eᵖ,sp.ħωᴸᴼ,-sp.ħωᴸᴼ))
    container = NEGFProblem(p,n,ScatteringOptions(LO=false,acoustic=false,impurity=false,IFR=false,alloy=false),
        s,g,profiles,basis,empty_kernels,shifts...,:nodal_linear)
    h = project_hamiltonians(container,zeros(n.N_z))
    return container, h, comparison
end

