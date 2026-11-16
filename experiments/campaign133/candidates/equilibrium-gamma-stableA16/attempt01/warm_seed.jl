# Included INSIDE EquilibriumGammaPilot; no execution on include.
const WARM_SOURCE_SHA = "a690ba56918d4c9803d70d759f2c228f4e9133fdc7f4bcede132f19229c0042f"
const WARM_MODEL_SHA = "6b41f93e29f2482a814f82a9db69ad8c8ff3bf7d9eb99ccfed4a510244e364b7"
const WARM_IDS = (:IFR,:LO,:acoustic,:alloy,:impurity)
const BASIS_KEYS = (("Phi","Φ"),("chi","χ"),("H0","H₀"),("Z","Z"),
                    ("Minv","M⁻¹"),("Tplus","T₊"),("Tminus","T₋"))

warm_path(cfg,key) = normpath(joinpath(cfg["_input_directory"],cfg[key]))
file_digest(path) = bytes2hex(open(sha256,path)) # SHA IO streaming, not read(path)
function checked_digest(path,expected)
    actual = file_digest(path)
    actual == expected || error("SHA mismatch: $path, actual=$actual")
    return actual
end
function array_digest(A::Array)
    bytes2hex(sha256(reinterpret(UInt8,vec(A))))
end
read_complex(group,name) = complex.(read(group[name*"_real"]),read(group[name*"_imag"]))
function read_complex_rows(group,name,rows::UnitRange{Int})::Field
    complex.(group[name*"_real"][rows,:,:,:],group[name*"_imag"][rows,:,:,:])
end
function create_complex_chunks(group,name,shape,rows)
    chunk = (min(rows,shape[1]),shape[2:end]...)
    return (HDF5.create_dataset(group,name*"_real",Float64,shape;chunk),
            HDF5.create_dataset(group,name*"_imag",Float64,shape;chunk))
end
function write_complex_rows(datasets,X::Field,rows::UnitRange{Int})
    datasets[1][rows,:,:,:] = real.(X)
    datasets[2][rows,:,:,:] = imag.(X)
end
close_complex(datasets) = foreach(close,datasets)

"Common scalar convex weights for ALL k and complex matrix entries."
function warm_interpolation_plan(source::Vector{Float64},target::Vector{Float64})
    length(source)>=2 && all(isfinite,source) && all(diff(source).>0) || error("bad source nodes")
    all(isfinite,target) && all(diff(target).>0) || error("bad target nodes")
    n = length(source)
    ghostleft,ghostright = source[1]-(source[2]-source[1]),source[n]+(source[n]-source[n-1])
    first(target)>=ghostleft && last(target)<=ghostright || error("target exceeds explicit adjacent zero ghosts")
    left,right = zeros(Int,length(target)),zeros(Int,length(target))
    weight = zeros(Float64,length(target)); ghost = zeros(Int8,length(target))
    for e in eachindex(target)
        j = searchsortedlast(source,target[e])
        l = min(j,n); r=l+1
        a = l==0 ? ghostleft : source[l]
        b = r==n+1 ? ghostright : source[r]
        theta = (target[e]-a)/(b-a)
        0<=theta<=1 || error("nonconvex interpolation row $e")
        left[e],right[e],weight[e] = l,r,theta
        ghost[e] = l==0 || r==n+1
    end
    return (left=left,right=right,weight_right=weight,ghost=ghost,
            ghost_left=ghostleft,ghost_right=ghostright)
end
function interpolate_gamma(X::Field,first_source::Int,plan,rows::UnitRange{Int},nsource::Int)
    Y = zeros(ComplexF64,length(rows),size(X,2),size(X,3),size(X,4))
    @inbounds for (local_e,e) in enumerate(rows)
        l,r,w = plan.left[e],plan.right[e],plan.weight_right[e]
        for b in axes(X,4),a in axes(X,3),m in axes(X,2)
            zl = l==0 ? 0.0+0im : X[l-first_source+1,m,a,b]
            zr = r==nsource+1 ? 0.0+0im : X[r-first_source+1,m,a,b]
            Y[local_e,m,a,b] = (1-w)*zl+w*zr
        end
    end
    return Y
end
function warm_broadenings(R::Field,L::Field,G::Field)
    K,Q = similar(R),similar(R)
    @inbounds for b in axes(R,4),a in axes(R,3),m in axes(R,2),e in axes(R,1)
        K[e,m,a,b] = im*(G[e,m,a,b]-L[e,m,a,b])
        Q[e,m,a,b] = im*(R[e,m,a,b]-conj(R[e,m,b,a]))
    end
    return K,Q
end

mutable struct WarmNorm
    n2::Float64; s2::Float64; nw2::Float64; sw2::Float64
    block_max::Float64; energy_index::Int; k_index::Int
end
WarmNorm() = WarmNorm(0,0,0,0,0,0,0)
function accumulate_norm!(d::WarmNorm,X::Field,Y::Field,we,wk,offset::Int)
    @inbounds for m in axes(X,2),e in axes(X,1)
        bn,bs=0.0,0.0
        for b in axes(X,4),a in axes(X,3)
            bn += abs2(X[e,m,a,b]-Y[e,m,a,b]); bs += abs2(X[e,m,a,b])
        end
        w=we[e+offset]*wk[m]
        d.n2+=bn; d.s2+=bs; d.nw2+=w*bn; d.sw2+=w*bs
        if sqrt(bn)>d.block_max
            d.block_max=sqrt(bn); d.energy_index=e+offset; d.k_index=m
        end
    end
end
function warm_norm_dict(d::WarmNorm)
    Dict{String,Any}("absolute"=>sqrt(d.n2),"scale"=>sqrt(d.s2),"relative"=>relative(sqrt(d.n2),sqrt(d.s2)),
        "weighted_absolute"=>sqrt(d.nw2),"weighted_scale"=>sqrt(d.sw2),
        "weighted_relative"=>relative(sqrt(d.nw2),sqrt(d.sw2)),
        "maximum_block_absolute"=>d.block_max,"energy_index"=>d.energy_index,"k_index"=>d.k_index)
end
function keep_cone!(records::Dict{String,Any},name::String,X::Field,offset::Int)
    v=Dict{String,Any}(string(k)=>value for (k,value) in pairs(cone(X;require_nonzero=false)))
    v["energy_index"]+=offset
    if haskey(records,name)
        old=records[name]
        margin=min(old["minimum_margin"],v["minimum_margin"])
        defect=max(old["hermiticity_absolute"],v["hermiticity_absolute"])
        old["minimum_eigenvalue"]<v["minimum_eigenvalue"] && (v=old)
        v["minimum_margin"],v["hermiticity_absolute"]=margin,defect
    end
    records[name]=v
end

function require_equal(name,A::Array,B::Array,digests::Dict{String,String})
    size(A)==size(B) && eltype(A)==eltype(B) && array_digest(A)==array_digest(B) ||
        error("raw model mismatch: $name")
    digests[name]=array_digest(B)
end
function verify_source_model(source,target)
    digests=Dict{String,String}()
    for key in ("energy","weight_energy")
        digests["source_"*key]=array_digest(read(source[key]))
    end
    for key in ("grid_ε","grid_wᴱ")
        digests["target_"*key]=array_digest(read(target[key]))
    end
    for (old,new) in (("k","κ"),("weight_k","wᵏ"),("x","x"),("weight_x","wˣ"))
        require_equal("grid_"*new,read(source[old]),read(target["grid_"*new]),digests)
    end
    for (old,new) in BASIS_KEYS
        require_equal("basis_"*new,read_complex(source["basis"],old),read_complex(target,"basis_"*new),digests)
    end
    Set(keys(source["kernels"]))==Set(string.(WARM_IDS)) || error("source local channels differ")
    Set(keys(target["kernels"]))==Set(string.(WARM_IDS)) || error("target local channels differ")
    for id in WARM_IDS
        old,new=source["kernels/"*string(id)],target["kernels/"*string(id)]
        read(old["qK"])==read(new["qK"]) || error("kernel scale differs: $id")
        digests["qK_"*string(id)]=array_digest([read(new["qK"])])
        require_equal("kernel_"*string(id),read_complex(old,"normalized_K"),read_complex(new,"K"),digests)
    end
    all(iszero,read(source["hartree"])) || error("source Hartree is not zero")
    read(source["energy_scale_eV"])==0.1 && read(source["length_scale_m"])==1e-8 &&
        read(source["spin_degeneracy"])==2 || error("source scale/spin contract differs")
    return digests
end
function verify_problem_model(c::Context,target)
    digests=Dict{String,String}();p=c.problem
    for name in fieldnames(typeof(p.grids))
        actual=getfield(p.grids,name)
        actual isa BitArray && (actual=Int8.(actual))
        require_equal("grid_"*string(name),actual,read(target["grid_"*string(name)]),digests)
    end
    for name in fieldnames(typeof(p.profiles))
        require_equal("profile_"*string(name),getfield(p.profiles,name),read(target["profile_"*string(name)]),digests)
    end
    for (_,name) in BASIS_KEYS
        require_equal("basis_"*name,getfield(p.basis,Symbol(name)),read_complex(target,"basis_"*name),digests)
    end
    require_equal("h_fixedU0",c.h,read_complex(target,"h_fixedU0"),digests)
    Set(p.kernels.enabled)==Set(WARM_IDS) || error("warm local channels differ")
    for id in WARM_IDS
        group=target["kernels/"*string(id)]
        p.kernels.qᴷ[id]==read(group["qK"]) || error("target kernel scale differs: $id")
        digests["qK_"*string(id)]=array_digest([p.kernels.qᴷ[id]])
        require_equal("kernel_"*string(id),p.kernels.K[id],read_complex(group,"K"),digests)
    end
    return digests
end
function verify_model_inputs(cfg)
    source=TOML.parsefile(warm_path(cfg,"source_input_path"))
    target=TOML.parsefile(warm_path(cfg,"target_input_path"))
    pairs=(("Nz","Nz"),("Nb","Nb"),("basis_periods","basis_periods"),
           ("Nk","Nk"),("kmax_per_nm","kmax_per_nm"),("Nphi","Nphi"),
           ("Nqz","Nqz"),("qzmax_per_nm","qzmax_per_nm"),
           ("temperature_K","temperature_K"),("trust_margin_meV","trust_margin_meV"))
    for (a,b) in pairs
        source[a]==target[b] || error("input physical/numerical mismatch: $a")
    end
    source["voltage_mV"]==0 && source["alloy"]===true && source["embedding"]=="finite_chain" &&
        source["embedding_periods"]==4 && source["NE"]==30721 && target["NE"]==30704 || error("input scheme differs")
    source["energy_min_eV"]==target["energy_min_eV"] && source["energy_max_eV"]==0.9 &&
        target["energy_max_eV"]==0.9000001064962726 || error("window identity differs")
    old=TOML.parsefile(warm_path(cfg,"source_identity_path"))
    old["core_head"]=="f50fa9bdfc04e299f82b6260f70fd8f829d838f3" &&
        old["runner_head"]=="18f45edb3bd086ebc7d5c984d68266a937358d4e" || error("source identity differs")
    old["input_sha256"]==cfg["source_input_sha256"] || error("source input/identity differs")
    current=TOML.parsefile(warm_path(cfg,"target_identity_path"))
    current["core_head"]==cfg["expected_core_head"] && current["core_tree"]==cfg["expected_core_tree"] || error("target source differs")
    return nothing
end
function same_pilot_options(cfg,cold)
    for (key,value) in cold
        key in ("identity","runtime_authorized") && continue
        if key=="affine_l1_max"
            value==10.0 && cfg[key]==4096.0 || error("approved affine control delta differs")
            continue
        elseif key=="data_max_bytes"
            value==6*1024^3 && cfg[key]==8*1024^3 || error("approved storage-only delta differs")
            continue
        end
        haskey(cfg,key) && isequal(cfg[key],value) || error("warm changed a frozen cold option: $key")
    end
    return nothing
end

function build_source_gamma(sourcepath,output,cfg,checkwall,report)
    HDF5.h5open(sourcepath,"r") do source
        energy,weights,wk=read(source["energy"]),read(source["weight_energy"]),read(source["weight_k"])
        shape=size(source["selfenergy/LO/lesser_real"])
        shape==(30721,11,5,5) || error("Julia HDF source axis order/shape differs")
        records=Dict{String,Any}();stats=Dict(id=>WarmNorm() for id in (WARM_IDS...,:total))
        atomic_h5(joinpath(output,"source-gamma.h5"),cfg;expected_bytes=16prod(shape)+8*1024^2) do file
            file["schema"]="campaign133-source-GammaK-v1"
            file["source_state_sha256"]=cfg["source_state_sha256"]
            file["energy"]=energy;file["weight_energy"]=weights;file["weight_k"]=wk
            datasets=create_complex_chunks(file,"GammaK_raw",shape,cfg["chunk_energies"])
            try
                for start in 1:cfg["chunk_energies"]:shape[1]
                    checkwall(); rows=start:min(start+cfg["chunk_energies"]-1,shape[1])
                    report["source_chunk_started"]=start
                    sumK=zeros(ComplexF64,length(rows),shape[2:end]...);sumR=similar(sumK);fill!(sumR,0)
                    for id in WARM_IDS
                        group=source["selfenergy/"*string(id)]
                        R,L,G=(read_complex_rows(group,name,rows) for name in ("retarded","lesser","greater"))
                        K,Q=warm_broadenings(R,L,G)
                        accumulate_norm!(stats[id],K,Q,weights,wk,start-1)
                        keep_cone!(records,string(id)*"_GammaK",K,start-1)
                        sumK .+= K;sumR .+= Q
                    end
                    write_complex_rows(datasets,sumK,rows)
                    accumulate_norm!(stats[:total],sumK,sumR,weights,wk,start-1)
                    keep_cone!(records,"total_GammaK",sumK,start-1)
                end
            finally
                close_complex(datasets)
            end
        end
        report["source_cones"]=records
        report["source_GammaK_vs_GammaR"]=Dict(string(k)=>warm_norm_dict(v) for (k,v) in stats)
        for (id,d) in stats
            metric=warm_norm_dict(d)
            metric["relative"]<=cfg["arithmetic_relative"] && metric["weighted_relative"]<=cfg["arithmetic_relative"] ||
                error("source retarded/kinetic identity exceeds arithmetic target: $id")
        end
    end
end

function build_target_seed(targetpath,output,cfg,checkwall,report)
    HDF5.h5open(targetpath,"r") do target
        HDF5.h5open(joinpath(output,"source-gamma.h5"),"r") do source
            old,new=read(source["energy"]),read(target["grid_ε"])
            plan=warm_interpolation_plan(old,new);shape=(length(new),11,5,5)
            records=Dict{String,Any}();correction=WarmNorm()
            we,wk=read(target["grid_wᴱ"]),read(target["grid_wᵏ"])
            atomic_h5(joinpath(output,"gamma-seed.h5"),cfg;expected_bytes=32prod(shape)+12*1024^2) do file
                file["schema"]="campaign133-warm-GammaK-ghost0-v1"
                file["source_state_sha256"]=cfg["source_state_sha256"]
                file["target_model_sha256"]=cfg["target_model_sha256"]
                file["energy_scale_eV"]=0.1;file["length_scale_m"]=1e-8;file["spin_degeneracy"]=2
                axes=HDF5.create_group(file,"axes")
                for key in keys(target)
                    startswith(key,"grid_") && write_metadata(axes,key,read(target[key]))
                end
                mapping=HDF5.create_group(file,"interpolation")
                mapping["source_energy"]=old
                for key in propertynames(plan);write_metadata(mapping,string(key),getproperty(plan,key));end
                raw=create_complex_chunks(file,"Gamma_import_raw",shape,cfg["chunk_energies"])
                work=create_complex_chunks(file,"Gamma_seed",shape,cfg["chunk_energies"])
                try
                    for start in 1:cfg["chunk_energies"]:length(new)
                        checkwall();rows=start:min(start+cfg["chunk_energies"]-1,length(new))
                        report["target_chunk_started"]=start
                        firstold=max(1,minimum(plan.left[rows]));lastold=min(length(old),maximum(plan.right[rows]))
                        X=read_complex_rows(source,"GammaK_raw",firstold:lastold)
                        Y=interpolate_gamma(X,firstold,plan,rows,length(old))
                        Z=hermitian_copy(Y) # explicit representation only, never eigen-clipping
                        write_complex_rows(raw,Y,rows);write_complex_rows(work,Z,rows)
                        keep_cone!(records,"raw_import",Y,start-1);keep_cone!(records,"working_seed",Z,start-1)
                        accumulate_norm!(correction,Y,Z,we,wk,start-1)
                    end
                finally
                    close_complex(raw);close_complex(work)
                end
            end
            report["target_cones"]=records;report["representation_conversion"]=warm_norm_dict(correction)
            report["last_source_weight"]=1-last(plan.weight_right)
            report["last_ghost_weight"]=last(plan.weight_right)
            report["ghost_right_scaled"]=plan.ghost_right
        end
    end
end

function load_warm_seed(c::Context,cfg,output)
    cfg["expected_source_state_sha256"]==WARM_SOURCE_SHA || error("warm source contract differs")
    cfg["target_model_sha256"]==WARM_MODEL_SHA || error("warm target contract differs")
    manifestpath=warm_path(cfg,"seed_manifest_path")
    checked_digest(manifestpath,cfg["expected_seed_manifest_sha256"])
    manifest=TOML.parsefile(manifestpath)
    manifest["schema"]=="campaign133-warm-seed-manifest-v1" || error("manifest schema differs")
    manifest["source_state_sha256"]==WARM_SOURCE_SHA && manifest["target_model_sha256"]==WARM_MODEL_SHA || error("seed provenance differs")
    manifest["execution"]=="complete" || error("builder did not complete")
    manifest["core_head"]==cfg["expected_core_head"] && manifest["core_tree"]==cfg["expected_core_tree"] || error("builder Core differs")
    coldinput=joinpath(dirname(manifestpath),"target-input.toml")
    checked_digest(coldinput,manifest["target_input_sha256"])
    same_pilot_options(cfg,TOML.parsefile(coldinput))
    seedpath=joinpath(dirname(manifestpath),manifest["seed_file"])
    checked_digest(seedpath,manifest["seed_sha256"])
    modelpath=warm_path(cfg,"target_model_path")
    checked_digest(modelpath,cfg["target_model_sha256"])
    fingerprints=HDF5.h5open(modelpath,"r") do target;verify_problem_model(c,target);end
    gamma=HDF5.h5open(seedpath,"r") do seed
        read(seed["schema"])=="campaign133-warm-GammaK-ghost0-v1" || error("seed schema differs")
        read(seed["source_state_sha256"])==WARM_SOURCE_SHA && read(seed["target_model_sha256"])==WARM_MODEL_SHA || error("seed source differs")
        read(seed["energy_scale_eV"])==c.problem.scales.E₀_eV && read(seed["length_scale_m"])==c.problem.scales.L₀_m &&
            read(seed["spin_degeneracy"])==c.problem.physical.g_s || error("seed scale/spin differs")
        for name in fieldnames(typeof(c.problem.grids))
            actual=getfield(c.problem.grids,name);actual isa BitArray && (actual=Int8.(actual))
            isequal(actual,read(seed["axes/grid_"*string(name)])) || error("seed axes differ: $name")
        end
        expected=(c.problem.numerical.N_E,c.problem.numerical.N_k,c.problem.numerical.N_b,c.problem.numerical.N_b)
        size(seed["Gamma_seed_real"])==expected || error("seed shape/Julia axis order differs")
        read_complex(seed,"Gamma_seed")::Field
    end
    witness=cone(gamma)
    record=Dict{String,Any}("initialization"=>cfg["initialization"],"manifest_sha256"=>cfg["expected_seed_manifest_sha256"],
        "seed_sha256"=>manifest["seed_sha256"],"source_state_sha256"=>WARM_SOURCE_SHA,
        "seed_cone"=>Dict(string(k)=>v for (k,v) in pairs(witness)),"actual_problem_fingerprints"=>fingerprints,
        "scope"=>"initial Gamma only; fresh current mu and empty AA history; not exact native restart")
    open(joinpath(output,"warm-initialization.toml"),"w") do io;TOML.print(io,record);end
    cp(manifestpath,joinpath(output,"seed-manifest.toml"))
    return gamma
end
