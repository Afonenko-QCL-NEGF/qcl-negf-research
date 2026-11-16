# Source-free tests; definitions only until ROOT requests --tiny.
function tiny_checks()
    cfg=Dict{String,Any}("temperature_K"=>70.0,"Nz"=>24,"Nb"=>2,"basis_periods"=>1,
        "NE"=>17,"energy_min_eV"=>-.3,"energy_max_eV"=>.9000001064962726,
        "trust_margin_meV"=>80.0,"kmax_per_nm"=>.6,"Nk"=>3,"Nphi"=>4,"Nqz"=>5,
        "qzmax_per_nm"=>10.0,"seed_halfwidth_meV"=>1.25,"workers"=>1,"memory_max_bytes"=>13019119616)
    @testset "private cache equality and native guard isolation" begin
        c=GP.build_context(cfg); cache,ledger=gamma_cache(c.problem,c.options)
        @test ledger["Gamma_envelope_bytes"]<cfg["memory_max_bytes"]
        @test cache.estimate.peak_bytes==estimate_production_memory(c.problem;options=c.options).peak_bytes
        @test cache.fft_hilbert_plan===nothing
        @test cache.product_hilbert_operator isa N._ProductionPIFFTHilbertPlan
        @test all(s->s isa EnergyShiftPlan,(cache.W̃₊ᴱᵖ,cache.W̃₋ᴱᵖ,cache.W̃₊ᴸᴼ,cache.W̃₋ᴸᴼ))
        G=zeros(ComplexF64,17,3,2,2); A=ComplexF64[1 .2im;-.2im 2]
        for e in 1:17,k in 1:3; G[e,k,:,:].=sinpi((e-1)/16)^2 .* A; end
        G[1,:,:,:].=0; G[end,:,:,:].=0
        L,equal=cache_equal(c.cache,cache,c.problem,G)
        @test equal["exact_value_equality"]
        @test N._validated_production_cache(cache,c.problem,c.options) isa ProductionCache
        low=N.with_production_options(c.options;memory_budget_bytes=1)
        @test_throws ArgumentError N._validated_production_cache(cache,c.problem,low)
        f,fl=refined_context(c)
        @test f.problem.numerical.N_E==33
        @test f.problem.grids.ε[1:2:end]==c.problem.grids.ε
        @test f.problem.kernels===c.problem.kernels && f.problem.basis===c.problem.basis
        @test f.h===c.h && f.options.memory_budget_bytes==c.options.memory_budget_bytes
        @test fl["native_validation_refused_for_full_budget"]==false
        P=mktempdir() do dir
            path=joinpath(dir,"fixture.h5")
            h5open(path,"w") do file; GP.write_complex(file,"Gamma",G); end
            h5open(path,"r") do file
                prolong_field(file,"Gamma",c.problem.grids.ε,f.problem.grids.ε)
            end
        end
        LF,_=N._production_pi_fft_hilbert_transform(P,f.problem.grids.ε,f.problem.grids.wᴱ,f.cache.product_hilbert_operator)
        @test LF[1:2:end,:,:,:]≈L atol=2e-13
        # Convex matrix interpolation preserves source values; no repair.
        @test P[1:2:end,:,:,:]==G
        @test GP.cone(P).minimum_margin>=0
        @test record(GP.norms(G,G,c.problem.grids.wᴱ,c.problem.grids.wᵏ))["absolute"]==0
    end
end
