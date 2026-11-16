# Bounded fine2x Γ stage: write-only, ROOT static/tiny/source freeze first

Новый standalone harness, `<250LOC`, reuse frozen definitions. Автор Julia,
Pkg, HDF5 numerical arrays и новые Π не запускал. Input runtime=false.
Source — **unprojected stableA-cont64 Γ-last64**, не projected campaign и
не прежний subtractive-A Γ128. Source whole hashes заданы ROOT в
`../stableA-cont64-source-digests.json`, зафиксированы в input.

## Уравнения и execution identity

Frozen Core e171, source Γ module50bcfaae…, stable spectral definition
eeca1f29…; streaming8df7885c… реализует ту же native covariance A по
plan.hlocal all9. Preinclude hashes проверяются **до include** всех4
definitions, stream и owning cache helper7d71dadc… из immutable предыдущего
attempt. Qualifier c48d6e4f… имеет original arithmetic equivalence pass;
runtime gate связывает его source hashes/Γ payload с этим же Γ-last64.
ROOT независимый result review qualifier остаётся prerequisite actual;
если qualifier признан невалидным, actual stage не запускается.

Во всех maps:

`Γ → HPIΓ−iΓ/2 → fresh finite4 hH resolvent/covariance A → q(A) → current μ
fit → FD partners → native five-channel raw Φ`.

Не меняются h/K/basis/k quadrature, window, LO bath, channels, finite4,
current-μ equation и raw criterion. Статический Γ physical/source operator
не превращается в direct_green. Νовая fine canonical grid30704→61407
с теми же endpoints, nominal Ω/h939→1878. Actual Float64 old nodes должны
совпасть с fine odd nodes; otherwise stop. PΓ — common convex weights из
actual source/fine coordinates, без Γ repair/rescale/eigenclip. Fine wE,
compact shifts/cache/cavity пересобраны; старые weights не наследуются.

**History3→1** — отдельная execution identity; не делать isolated
acceleration/convergence inference из сравнения source m3 и fine m1.
`project_AA_trial=false`. Рабочая representation сохранена от parent:
`GP.gamma_residuals(rawΦ,Γ)` даёт Herm(rawΦ) для AA и отдельный conversion
ledger; `F_AA=Herm(rawΦ)−Γ`. **RAW** Φ не PSD-project/clip, сохраняется и
задаёт raw residual/trigger. Hermitian representation restore не выдаётся
за PSD projection. Исходный raw norm не заменяется working norm.

AAII history1, real weighted QR coefficients, β1, fallback α.2,
condition≤1e8, affine l1≤4096. Native source `aa_trial` reuse.
Каждый trial проходит исходный Γ cone/backward domain, без projector.
AAII domain failure допускает только один convex fallback от последнего
valid Γ/F и reset history; fallback считается в budget. Если fails fallback
или non-AA trial, stop; последнее aligned Γ/rawΦ сохраняется.

## Call accounting и этапы

Жёстко максимум **1 coarse initial +12 fine candidate evaluations +1
terminal coarse restriction =14 attempts**. Conservative attempts считаются
до cone, поэтому precone rejects тоже расходуют12; отдельно для каждого
group записываются `entered_Pi` (после Γ cone), `completed_Phi` (возврат
native Φ), `successful` (все domain gates и запись). Это позволяет отделить
реальную physical map работу от rejected attempts и честно сравнивать counts.

1. Source model/arrays/manifest/Γ payload binding. Public admitted coarse
   cache vs private exact shift/kernel/full FFTΛ equality, без Φ controls.
   Один fresh streamed coarse Φ, replay saved source rawΦ64≤1e−10.
2. Same bound frame → fine objects/private cache → PΓ. Один seed fine Π,
   затем до12 attempts total с history1. Strict raw+weighted1e−9 ×3
   streak считается с seed и сбрасывается после failure. Это reduced
   progress only, **никакого native exit/SCBA acceptance**.
3. Перед terminal stage history/QR/F/trial lifetimes освобождаются.
   Γterminal restricted на old nodes → один fresh coarse Π. Сравнивается
   coarse Φ(RΓterminal) с RΦfine(Γterminal), coarse wE/wk и current-μ каждой
   map. Это diagnostic cross-grid response, не дополнительная итерация.

Fine loop прекращается при750s, резервируя время terminal/checkpoint
в общем soft840/external900s. Если domain/API/wall stop не позволяет
terminal map, он остаётся not_measured; attempts никогда не расширяются.
Сохранение `fine-terminal` идёт в finally для последнего aligned результата.

## Cache/workload memory contract

Reuse private owning cache7d71: compact shift plans/native kernel maps,
selected PIFFT only, no dense H/unused nodal FFT. PI plan worker1/columns1
явно отличается от streamed/collision options/J4; coarse exact equality
контролирует значения. Fine public `build_production_cache` обязан вернуть
**прежний native full-peak budget refusal до allocations** при real cap.
Ни public guard, ни numeric cap не повышаются; честный native estimate
~25.84GB сохраняется в private cache. Native solver не вызывается.

Ffine=270190800B=257.674MiB. Новая stage admission:
`31F+3GiB≈10.80GiB ≤12.125GiB`. Это conservative history1/temporary/runtime
envelope, не blindly вычтенная native estimate. Constructor ledger24F+3GiB
сохраняется отдельно; новый31F envelope покрывает4 retained Γ/F fields,
aligned rawΦ, trial/QR/conversion temporaries и initial/terminal snapshots.
Stream4F относится только к Γ/S/GR/A sweep, **не всему процессу**;
PI/GreenFD/Φ/LO/history/cache/JIT/HDF живут дополнительно. Original family
cavity9F arrays не materialize в streaming primitive. Actual RSS измеряется.

Default Core builder/validator остается e171 unchanged. Private stage
guard не разрешает full-native workload с тем же cap; public fine refusal
явно проверяется и записывается. Контексты coarse/fine retained только как
малые frame/cache/plans; old coarse Green/Phi не удерживаются.

## Артефакты и scientific limits

- `fine-seed-Gamma.h5`: regridded Γ до первого fine map, включая failure case;
  `fine-initial.h5`, `fine-terminal.h5`: Γ и **rawΦ**, energy/wE.
- `coarse-restriction-terminal.h5`: restricted Γ/fresh coarse rawΦ.
- Per attempt `response_<group>_<n>.h5`: полный q, axes, predeclared≤128
  selected spectral blocks. Записывается **до μ fit**, чтобы q/domain failure
  не скрывал spectral evidence. μ — в linked successful map TOML; иначе
  явно not_measured, не ноль.
- Per completed Φ `map_...toml`: raw/working/conversion norms, μ/number,
  cones/backward allowances, pair/A, conditioning, channel norms, stream
  source/Λ/h ledger. First post-Φ domain failure сохраняет Γ/rawΦ snapshot.
- Events: AA coefficients, suffix reasons, attempted/fallback/rejection
  reasons; summary содержит counts и последний strict reduced streak.

Source replicas old spectral indices и chosen fine old/midpoint neighbours
ограничены pick(): source anomaly382/525/8040, μ-region7894/7907,
residual9785 и LO descendants/ancestors8846/10724/11663 и широкие endpoints.
Новые source rows выбираются с явным **`2 * e`**, без Julia exponent trap.
Смена default ghost taper h→h/2 — explicit boundary variant. Actual
source endpoint bound `log2/(2π)(||Γfirst||F+||Γlast||F)` записывается как
uniform-in-E combined-k bound на Λ delta, **не** global energy norm или
resolvent/density bound. Нельзя назвать всю causal continuation unchanged.

Fine finite-iteration progress не доказывает достаточность h/2, continuum
solution или native equations. Current-μ fitting не acceptance. Native gates
и source1e−9 tolerance не ослабляются. Новая map trajectory не принимает
старую задачу задним числом; source/all earlier failures сохраняются.

## Tiny/static/ROOT execution

Source-free `--tiny` проверяет private/public PI equality J4/plan1,
public memory refusal, old stable covariance versus stream GR/A/q,
raw-vs-Hermitian representation semantics, known fine indices19569,
≤128 sample bound,31F envelope и analytic affine fixed point для AA1.
Автор tests не запускал. Static/source freeze необходимы до ROOT tiny;
runtime input сохраняется отдельно с authorization=true только после
qualifier result review и tiny pass.

External hard13019119616B/High12348030976B/Swap0, CPU400%, J4/BLAS1,
heap hint8G,900s/soft840, file2GiB/data4GiB. Source gate читает реальные
cgroup controls и finite LimitFSIZE. ROOT до запуска проверяет MemAvailable
≥hard+1400MiB и disk≥4GiB+margin, исключает overlapping heavy computation.
Old artifacts не удаляются. Data/file caps включают HDF padding.

```bash
timeout 180s env JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  julia --startup-file=no --project=.build/campaign133/environment/accepted \
  .build/campaign133/equilibrium-gamma-fine2x-m1/run_fine.jl --tiny
```

ROOT задаёт tiny отдельный finite memory scope. Actual пример, только
после review и отдельной runtime-input copy:

```bash
systemd-run --user --scope -p MemoryMax=13019119616 -p MemoryHigh=12348030976 \
  -p MemorySwapMax=0 -p CPUQuota=400% \
  timeout --signal=TERM 900s bash -c 'ulimit -f 2097152; exec env JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --heap-size-hint=8G --project=.build/campaign133/environment/accepted .build/campaign133/equilibrium-gamma-fine2x-m1/run_fine.jl .build/campaign133/equilibrium-gamma-fine2x-m1/runtime-input.toml .build/campaign133/equilibrium-gamma-fine2x-m1/attempt01'
```

Не запускать повторные maps для восстановления только presentation.
Если postprocess fails после successful maps, сохраняются raw snapshots
и counts; восстановление выполняется отдельным read-only postprocessor.
