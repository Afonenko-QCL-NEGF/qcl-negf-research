# Frozen collocation control: write-only candidate, static review before ROOT runtime

Новый отдельный harness; исходный Γ128/Core/исторические отчёты не менялись.
Автор не запускал Julia, Pkg, новые карты или HDF5 arrays. `input.toml`
оставлен `runtime_authorized=false`. Source-free `--tiny` — отдельный ROOT
этап до scientific data. После static/tiny/source freeze ROOT выбирает
resource gap и новую runtime-input копию с явно указанной авторизацией.

Научный вопрос: **какова grid collocation sensitivity Φ на одном и том же
Γ128 input?** Не проверяется nonlinear convergence. Read-only обоснование,
allocated/live-set и future LO-closed nonuniform alternative:
`../adjacent/residual-support-refinement-plan.md`.

## Численный объект

Frozen Γ module `equilibrium-gamma-hermitian-cont128/attempt01/gamma_functions.jl`,
SHA57f2c277…; include загружает definitions, без исполнения его128-loop.
Переиспользуются `retarded_state`, `phi`, `cone`, `fit_mu`, `fd_green`,
`stencil_witness`, streaming norms/writer. Harness не вызывает `evaluate`,
native compare/handoff/public solver, AA, mixing или Poisson. Included parent
module содержит их definitions, но CLI не имеет такого execution path.

Source terminal `fb7cea34…`, model47682acc…, map1288a34e3ac… связываются
через SHA и source identity; included warm_seed/gamma_restart также bound
source identity digests. Coarse model/basis/h/kernels/axes сравниваются с
actual rebuilt admitted native context raw array hashes. Поле Γ читается
из terminal без repair; tiny admitted negative eigenvalues сохраняются
как в parent `cone`/native backward policy. Никакого eigenvalue clipping,
zeroing q или η floor. Domain failure — bounded negative outcome.

Две Π starts TOTAL: coarse fresh replay при worker1 и fine fresh map.
`started` увеличивается до domain/retarded работы; failed starts входят в2.
После каждого completion возвращается только raw Φ и small scalars;
Green/cavity-side/LO buffers не переживают barrier. Ни итераций, ни retries.

Coarse NE30704 → fine61407, те же endpoints −.3/.9000001064962726eV,
Ω=.0367eV, m939→1878. Canonical fine odd nodes обязаны bitwise совпасть
с source. Scalar convex interpolation Γ использует **actual coordinates**
для каждого even node; weight не форсируется к0.5. Одни веса для всех
complex entries/k, source odd values копируются. Новые grids/wE/compact
shift/cavity/cache consistent с NumericalParameters. Basis/profiles/K/h
переиспользуются из bound source frame, energy-dependent arrays — fresh.

В каждой Π сначала на свежем GR вычисляется N(source μ) как diagnostic;
только затем fitting собственного current μ и raw Φ. `μfit→Ntarget` не
называется convergence. Π_f включает изменение quadrature/current μ,
поэтому defect не выделяет уникально LO/KK/μ/Dyson причину.

## Cache construction admission отдельно от native workload admission

Coarse использует public `GP.build_context` → `build_production_cache`
с реальным cap12.125GiB. Дополнительный private cache должен иметь **точно
равные** shift fields/kernel fields и все full FFTΛ values на source Γ;
это control до fine, без дополнительной Φ.

Private `gamma_cache` строит owning ProductionCache из frozen native
kernel/sparse-shift/PIFFT constructors. Plain nodal FFT plan=nothing;
один selected PIFFT plan с primary/doubled-padding verification, workers1,
columns1. Не создаются dense H/dense shift matrices/full-family histories.
Source grid/kernel ownership contracts сохранены. В cache.estimate остаётся
**настоящий native full estimate**, хотя он около24GiB для fine. Options
memory_budget=actual12.125GiB не увеличивается. Отдельный Γ stage envelope
`24F+3GiB≈9.04GiB` admission проверяется до allocations; payload/cavity
измеряются отдельно и должны вместе быть≤.5GiB. F=257.674MiB, hist0;
native primitive live-set peak≈13F плюс cache/scalars. Это practical
envelope/reserve, не proof RSS. External MemoryMax — enforcement.

Fine `_validated_production_cache` отдельно проверяется: если native full
estimate превышает cap, он обязан вернуть свой прежний memory refusal.
Никакого вызова full solver при этом нет. Default Core builder/validator
не менялись. Tiny fixture проверяет равенство private/public cache,
fullguard refusal при малом budget, owning source identities, nested
geometry, actual-coordinate HDF prolongation и HPI invariance zero-edge hats.

## HPI/ghost control

Coarse native/private FFTΛ equality проверяется по всем entries; output
coarse Lambda сохраняется. Fine HPI(PΓ), restricted на old nodes, сравнивается
с coarse Lambda в исходной coarse measure. Relative arithmetic diagnostic
≤1e−10, без ослабления original map tolerances. Всего3 FFT controls outside
Π (native/private coarse +fine), плюс один internal HPI в каждой из2Π;
они не считаются бесплатными по времени, но не переименовываются в Φ maps.

При refinement default FFT ghost width h→h/2: **explicit boundary variant**.
Continuous physical-window linear Γ сохраняется вплоть до floating
interpolation errors, а старые outer ghost tapers — нет. Сохраняется
actual source endpoint norm bound
`log2/(2π)(||Γfirst||F+||Γlast||F)` на δΛ внутри physical window,
отдельно от measured HPI discrepancy/roundoff. Это kernel exact-arithmetic
bound на taper delta, не whole-program/resolvent certificate. Нельзя выдать
маленький HPI difference за unchanged full model.

## Выходы и сравнение

- `coarse.h5`: Gamma, Lambda, raw Phi, energy/wE; ни Green, ни families.
- `fine-input.h5`: Gamma/axes до второго Π, immutable checkpoint.
- `fine.h5`: raw Phi/axes после второй Π.
- `summary.toml`: starts/completed, fixed-μ N и fitted μ, raw norms, cones,
  FD pair/A identity, channels, Dyson condition/scale, actual LO/KMS/adjoint
  stencil witness, cache/native admission, source hashes, timings/RSS.
- `selected-defect.csv`: anchor9785 ±3 coarse nodes, replicated по всем
  ±n939 ladders внутри window; fine old/midpoint neighbours. Все11 k,
  raw r², collocation D², reference Φ² и Γtrace. Полный D_f norm вычисляется
  streamed по **всем fine rows**, отдельно new-even-node norm и coarse
  restriction norm. Weights — actual owning grid, source criterion не меняется.

`D_f=Φ_f(PΓ)−PΦ_c(Γ)`. Coarse restriction `D_c=RΦ_f−Φ_c`, denominator и
measure source Φ_c/coarse wE. Сохраняются absolute/weighted/unweighted
norms без additive floor. На new midpoints отдельно вычисляется defect,
который old-node comparison мог бы пропустить. CSV не является заменой
полных raw norms. Согласие coarse Φ с saved source≤1e−10 должно быть до
fine; HPI mismatch также останавливает control. Никакого native acceptance.

## Ресурсы и команды ROOT

Один процесс J1/BLAS1/CPU100%, hard memory13019119616B=12.125GiB,
High12348030976B=11.5GiB, Swap0. Root admission до запуска:
MemAvailable≥hard+1400MiB, свободный диск≥3GiB+margin, других heavy нет.
External900s, soft840s, ≤2 map starts включая failures; file2GiB/data3GiB.
Readback memory/high/swap/CPU и LimitFSIZE проверяется source_gate.
Writer учитывает padded HDF chunks. Original source/output не удаляются.
Checkpoint `summary.toml` сохраняется между stages; external termination
может оставить последний завершённый checkpoint и .tmp raw writer.

До actual ROOT запуска нужно static review и `--tiny`, затем SHA freeze.
Пример source-free preflight, существующая среда accepted, **без Pkg**:

```bash
timeout 180s env JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --startup-file=no --project=.build/campaign133/environment/accepted \
  .build/campaign133/refinement-control/run_control.jl --tiny
```

ROOT задаёт для tiny отдельный конечный memory scope. Actual после review
использует runtime copy input с `runtime_authorized=true`, original copy
сохраняется; SHA runtime copy попадает в output. Пример внешней оболочки:

```bash
systemd-run --user --scope -p MemoryMax=13019119616 -p MemoryHigh=12348030976 \
  -p MemorySwapMax=0 -p CPUQuota=100% \
  timeout --signal=TERM 900s bash -c 'ulimit -f 2097152; exec env JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=.build/campaign133/environment/accepted .build/campaign133/refinement-control/run_control.jl .build/campaign133/refinement-control/runtime-input.toml .build/campaign133/refinement-control/attempt01'
```

Здесь bash `ulimit -f` использует1024-byte blocks; фактический readback
/proc/self/limits авторитетен. Existing output отклоняется. Outcome любого
запуска — frozen collocation diagnostic, не mesh adequacy certificate,
не доказательство причины stagnation и не SCBA pass.
