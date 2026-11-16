# FD64 transform-only harness: pending ROOT runtime decision

Статус `written_not_executed`: автор не запускал Julia, Pkg, HDF5 reads или
actual numerical probe. Исходный module a74ba641… и tests3be8f0d8… не менялись.
`run_frozen.jl` — отдельный harness,156LOC. Это дополнение к `protocol.md`
с более узкой заранее объявленной sampling policy, не изменение frozen design.

Порядок ROOT:

1. Source-free `--preflight`: include tiny tests, load QCLNEGF, bind typed
   aliases `CT.PanelTree`, `PIPlan`, `PIFunction`. Никаких source arrays.
2. По resource gap после warm-Γ/IFR приоритета — actual new-output run.
   Before source HDF5: frozen module/tests SHA, accepted imported Core
   `e171caa…` clean, prior summary SHA, source/writer/input/markers identities.
3. Source SHA проверяется **один раз** до H5open: file3,797,438,207bytes.
   Original33-line loader не изменён. Harness воспроизводит его локальный
   local-R extraction recipe без трёх повторных hashes3.8GB. Sequential
   k1/6/11 hyperslabs; full all-k Γ не читается и не выделяется.
4. Basis H0/Z/Minv/T±5×5 и Phi/chi96×5 сверяются по dims и persisted array
   digests existing PI control. Raw Γ определяется `i(sumΣR−sumΣR†)`;
   embedding excluded. Никаких repair/Hermitization/clipping source.

Aggregate Γ reference SHA `edb8530e…` связывается с existing prior summary,
но **не заявляется recomputed selected-k hash**. Each gamma_k_sha256 и
energy_sha256 сохраняются отдельно. Полный source-file hash + reviewed
writer/extraction definition связывают selected Γ с исходным артефактом.

Input Hermitian/PSD numerical allowance явно равен
`64eps(Float64)*Nb*max(max_e||Γe||F,eps(Float64))`, Nb5, в source native
energy units. Это размерная backward-error-scale allowance eigensolve,
не global/interval PSD proof. Сохраняются actual Hermitian defect,
minimum computed eigenvalue и численное значение allowance. Непрохождение
останавливает run; source не изменяется.

Minimal pseudocode:

```text
tiny include + type aliases → frozen code/provenance guards → source hash
for k=1,6,11: load original Γ(E,k); build p24 tree; full accepted FFT PI
  cache independent BigFloat192 F(z−x) hat weights for four declared targets
  for τ=1e−9,1e−11: evaluate tree on16 source nodes +4 offgrid targets
    save full5×5 values/reference + ledger + runtime/RSS → diagnostic gate
close source → completed_diagnostic_only, or checkpointed failed/budget stop
```

Source nodes ровно16, объявлены в `frozen-input.toml`: endpoints/neighbours,
широкое равномерное покрытие и occupied-near-μ node7898. Offgrid4:
midpoints after source indices7681,15361, left ghost-strip midpoint и
`μ/E0+0.2i`. BigFloat192 reference использует **независимый positive-hat
divided difference F(t)=t log t**, а для real PV — F(t)=t log|t|.
F cache строится однажды по всем original source+ghost knots для target;
готовые scalar weights reusable across k. Не используются tree moments,
`panel_integral`, Core Hilbert или generic adaptive integration. Scalar
weights округляются в ComplexF64; обычный transpose, без conjugation;
forward/reverse contraction sensitivity сохраняется отдельно.

Real offgrid reference добавляет boundary `−iΓ(E)/2` ровно один раз через
независимую interpolation original knots. Complex target — ordinary Cauchy
без boundary term. Полная source Γ30721 участвует в каждой reference.
FFT reference считается на полном original grid на одном k/25 columns,
но tree sampled только в declared16 indices. Dense H/full solver cache
не строятся.

Оговорка геометрии: tree использует actual stored Float64 knots и adjacent
original ghost gaps. Core FFT PI использует canonical origin/spacing и
Toeplitz ideal uniform hats; их исходная floating coordinate discrepancy
сохраняется в FFT metadata, actual/ideal ghost coordinates и max difference.
Согласие в численном allowance — diagnostic witness, не доказательство
математического тождества actual rounded и ideal uniform geometries.

Gate каждой full5×5 matrix:
`norm(value−reference) <= ledger.far_truncation_estimate +
max(1e−11,1e−10*norm(reference))` в native source units E0=.1eV.
Roundoff sensitivities, moment reverse witness и far truncation идут
отдельными полями, без их переименования в verified forward bound.
`formal_whole_certificate=false`. Gate failure сохраняет матрицы/ledger
и останавливает run; allowance не изменяется. При первых failures второй
τ или следующие k могут остаться `not_measured`.

Предложенные внешние cap:120s wall,1536MiB RSS,Swap0,CPU100%,J1/BLAS1,
file32MiB; internal soft110s. ROOT выбирает конечный MemoryMax после
preflight; если JIT требует2GiB, это явно новый enforced cap и записывается
в admission/tool log. Script `Sys.maxrss()` и elapsed включают прошедший
process lifetime/RSS, soft timer начинается после tiny include; внешний
timeout120s покрывает также imports/JIT/tiny. Отдельные source-free preflight
лимиты60s/1536MiB. Никакого background запуска автором.

Команды ниже — для ROOT при свободном ресурсе; переменную executable ROOT
может заменить существующим Julia path. Никакого Pkg resolve/instantiate.

```bash
env JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --startup-file=no --project=.build/campaign133/environment/accepted \
  .build/campaign133/adjacent/cauchy-tree-prototype/run_frozen.jl --preflight

systemd-run --user --scope -p MemoryMax=1536M -p MemorySwapMax=0 -p CPUQuota=100% \
  timeout --signal=TERM 120s env JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --startup-file=no --project=.build/campaign133/environment/accepted \
  .build/campaign133/adjacent/cauchy-tree-prototype/run_frozen.jl \
  .build/campaign133/adjacent/cauchy-tree-prototype/frozen-attempt01
```

ROOT обязан также ограничить preflight штатным scope/timeout, если запускает
его новым процессом; сокращённая первая команда показывает только Julia
invocation. Existing output directory script отвергает. Новая попытка имеет
новый directory/admission, исходный output не перезаписывается.

Любой исход сообщает только transform implementation comparison.
No Green inversions/density/Matsubara match/SCBA/KMS solve. Никакой adaptive
economy claim; memory/time до actual root run остаются `not_measured`.
