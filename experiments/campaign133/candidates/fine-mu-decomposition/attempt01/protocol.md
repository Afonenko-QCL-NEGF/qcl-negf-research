# Fine μ-attribution: frozen spectrum, two thermal LO evaluations

WRITE-ONLY candidate; no Julia/Pkg/H5 execution by the author. ROOT/reviewer authorizes and owns any bounded run. Source `.build/campaign133/refinement-control/attempt01/{coarse.h5,fine-input.h5,fine.h5,summary.toml}` remains immutable; separate recovery postprocessing has a different owner and files. No new SCBA iteration, full Φ/Pi, μ-fit, grid choice, mixing, clipping or density normalization.

Question: which part of the saved collocation defect `D=Φf(PΓc,μf)−PΦc(Γc,μc)` is solely the scalar chemical-potential change on the SAME fine spectrum? Elastic Γ contractions depend on A, not μ. For thermal LO at fixed A/cache,

`M = ΓLO,f(μf) − ΓLO,f(μc)`;
`D_fixedμ = D − M = Φf(PΓc,μc) − PΦc(Γc,μc)`.

The equality is exact for this finite operator; evaluated fields have Float64 roundoff. This controls μ attribution, not stationarity, convergence or discretization acceptance. Retarded physics/cavity, fine Γ, PI, weights, spin, LO population and all kernels stay identical. We retain saved fitted μf/μc; μ is not re-fit.

High-level Julia algorithm: recreate the frozen coarse model once, use the unchanged frozen `refined_context`/`gamma_cache`; load fine Γ; `GR,…=GP.retarded_state(c,Γ)` once; `A=spectral_function(GR)` once; then for each stored μ call `fd_green(GR,A,…,μ)` and native `_production_lo_contraction` once. Form `im*(Σg−Σl)` and subtract in place. Finally stream `fine Φ`, actual-coordinate `P coarse Φ`, and fine Γ in 512-energy blocks, computing D/M/F and raw R/R−M over **all complex entries**. No static channel needs reevaluation.

## Identity/admission

Immutable Core e171caa60cce676c6c4238dd671e7c527f2bd977, clean tracked source, same Julia version as parent. `run_mu.jl` verifies frozen refinement source before include: run_control SHA5084c1f1…, operators SHA7d71dadc…, tiny_checks SHA565a8d18…. Included Gamma/warm/restart definitions are independently checked against fixed literal digests before loading; source functions are defined but neither original run_control nor one_map is called. Runtime summary records included and new code hashes.

ROOT fills all `UNSET` WHOLE hashes in a **new** runtime input: three H5, immutable summary, original parent model/input/identity. All are streaming-SHA checked BEFORE opening ANY source HDF5/reading fields. Model arrays are then bitwise verified through original `verify_problem_model`; canonical fine axes/nesting and fine Γ array digest are checked against saved summary. WHOLE hashes are streamed again at the end. Output never writes any source directory.

Keep the original numerical and ACTUAL OS limits equal: memory.max=13019119616B (12.125GiB), memory.high=12348030976B (11.5GiB), swap=0, reserve1400M per ROOT, Julia1/BLAS1, CPU≤100%. This matches the frozen constructor admission; **not** an 8GiB scope paired with a 12GiB numerical budget. The original fine Γ-only cache ledger is preserved; native full-solver validation still refuses its original over-budget estimate. This diagnostic does not call native solver or weaken guard logic. Previous actual ~5.4GiB/57s two-map result is a cost witness, not admission authority. At most one fresh retarded evaluation and two LO evaluations, zero full Φ/Pi; counters increment before expensive calls. Soft wall290s/external hard300s, data64MiB/file32MiB. External scope must enforce these values, heap hints alone do not.

## Required replay controls

1. Fresh q is the SAME frozen `_seed_number_weights_production(GR,…)` functional. At stored μf reproduce saved fine-map fitted N; at stored μc reproduce saved fine-map fixed-μ N. Both errors over target must be≤original arithmetic guard1e−10. Record dN/dμ, capacity and PI roundoff; no refit.
2. First LO uses stored μf. Its unweighted full-field Frobenius norm must reproduce saved `fine_map.channel_Gamma_norms.LO` within1e−10. This is the same `im*(greater−lesser)` functional as GP.phi, with qK/wk included once by native LO. Second LO changes μ only.
3. All output norms include unweighted and `wE_f*wk` weighted squared norms and real complex inner products. Report `||D||`, `||M||`, `||F||`, `||R||`, `||R−M||`, denominators `||PΦc||`, D·M,D·F,M·F and the squared decomposition identity. No floor: 0/0→0 by declared GP.relative convention, nonzero/0→Inf. Old-node and midpoint partitions use fine-grid weights; they are not coarse quadrature norms. Inner products are required because norm ratios cannot measure cancellation/attribution.
4. `selected-LO-mu.csv` samples old anchor9785 plus coarse LO sidebands939, local±3 and fine-neighbor±1, all k. Integer `2*e-1` indexing avoids the previous Julia scientific-literal trap `2e-1`. It reports block D/M/F squared norms, real inner products, and complex trace of M; no full field outputs.

Tiny controls (`--tiny-mu`, ROOT only): native thermal LO on a small independent complex PSD fixture against the independent Bose-cancelled analytic delta `K[shift+(Δf A)−shift−(Δf A)]`, affine f mixture/zero delta; actual nonuniform-coordinate prolongation and all-entry complex weighted norm/inner-product identities. Tiny mode does not read preserved source H5 or call the large control; its own temporary H5 only exercises prolongation. No packages/install are needed. These controls remain unexecuted until ROOT requests them.

Prospective commands (resource scope supplied by ROOT, separate tiny budget):

```sh
julia --startup-file=no --check-bounds=yes --threads=1 --project=.build/campaign133/environment/accepted .build/campaign133/fine-mu-decomposition/run_mu.jl --tiny-mu
julia --startup-file=no --threads=1 --project=.build/campaign133/environment/accepted .build/campaign133/fine-mu-decomposition/run_mu.jl NEW_BOUND_RUNTIME_INPUT NEW_OUTPUT_DIRECTORY
```

Outputs: copied small code/protocol/input, TOML summary/identity/counters, selected CSV only. There is no end-to-end solver or physical acceptance claim. An unexplained q/N/LO replay failure stops before attribution is interpreted; raw evidence and error reason are preserved.

V2 minimal performance/correctness revision (unexecuted): exact V1 source/input/protocol/checksums are archived in `v1-unexecuted/` (main f3bfd445…). The hot stream now updates a concrete `Matrix{Float64}(18,3)` with fixed scalar slots for full/old/midpoint partitions; all original formulas, per-slot summation order, output labels, denominators and rates remain the same. Named Dicts are constructed only once after streaming. This removes per-entry dynamic string concatenation/hash updates; no timing/allocation result is claimed before runtime. Every loaded fine/Pcoarse/Gamma chunk now rejects nonfinite entries before accumulation. Tiny controls compare all18 raw slots bitwise with the old dictionary update on the identical complex fixture, retain independent analytic norm/inner-product checks, and reject NaN in each of the three input sources. Source/kernel physics, stored μ and all call/resource limits are unchanged.
