# Γ-pilot16 V2: bounded math delta review

**Решение: V1 F1–F5 устранены на уровне проверенного source/delta;
нового math blocker для первого bounded pilot не найдено.** Это review
формул и evidence, не numerical acceptance S01 и не самостоятельный
resource-green. Actual S01 maps рецензент не запускал; Julia/Pkg/Core и
harness не менял. Overall guards отдельно проверяет b0.

## Frozen identity и имеющееся свидетельство

`sha256sum -c equilibrium-gamma-solver/SHA256SUMS`: все восемь OK.

| Файл | SHA256 |
| --- | --- |
| run_gamma.jl | `abda3fb908137eba303ac92aac35dcbf6e312a0dc6edb0081c39be72d7c0b46e` |
| gamma_functions.jl | `4475101718dab32cb5f4c4446390a496caaf2731366ca6e8e503c76232970f46` |
| checks.jl | `0942e463bc5ff5e26e721d17fb6dc9e72abec6f416d3b3187b029857f7c0f494` |
| pilot16.toml | `5dadb59f19fddbf0aa483687b7addb8a792d0e97681456022585b47796be3498` |
| protocol.md | `e38b511fd106b382a12771854218b5d78dcbae30ab8624fb0839401143fa733b` |
| handoff.md | `03cbe7c7752628fd1221bff5d9f7b2bcea9a5dc7fcd94c917d03e373140e10b7` |
| algorithm.md | `1382608befb3a70661d3c5cfc6593fdd9657138778eeaa9dd6eefd4ca7d428dd` |
| v1-to-v2.patch | `01d122e50c02bde91029223c6c807a86b1dd02e13cf1797174048157aae843a3` |

Native source HEAD `e171caa60cce676c6c4238dd671e7c527f2bd977`, tree
`742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8`, tracked clean — прежний.
Archive `unexecuted-v1` сохраняет исходные main/functions/checks SHA из
`gamma-math-code.md`. Diff input меняет только identity v1→v2; physical
parameters, grid, memory/time/call caps не меняет.

Прочитан уже созданный Root log
`environment/accepted/gamma-small-checks-v2.log`, SHA256
`aba6488b39dde5174573153ba0d2dafbd8b6222d6b678529912ad202213280a1`:
75/75 small controls, marker `GammaSmallV2ChecksOK_NoActualS01Maps`, exit0.
Это external evidence Root, не мой запуск; он устраняет parse/runtime
blocker small suite, но не исполняет actual Π/LO/KMS/cavity S01 precheck.

## F1–F5: проверка исправлений

**F1 parse.** `checks.jl:143` содержит однозначный `0.2 .* x .* x`.
Root log подтверждает parse/выполнение small suite на V2. Старый исход
не затёрт и actual S01 map за него не засчитан.

**F2 raw residual.** `gamma_residuals:170` возвращает четыре раздельных
объекта: Hermitian working Φ, raw `norms(rawphi,gamma)`, working residual
и representation error `norms(rawphi,mapped)`. В `evaluate:583–600`
raw absolute/denominator/weighted norms управляют streak/trigger;
`reduced_stage` best score тоже читает raw stats. Return F для AA остаётся
Herm(Φraw)−Γ. Таким образом AA работает в real Hermitian coordinates,
а anti-Hermitian raw defect не исчезает из scientific trigger.
Новый tiny-asymmetric control проверяет именно прежний failure mode:
cone-admissible imaginary defect может доминировать raw relative norm,
хотя working residual равен0. Это содержательный regression control.

**F3 known inverse failure.** `known_inverse_failure:18–31` принимает
SingularException(info≥0), LAPACKException(info>0) и только непустые
TaskFailed/Composite/Captured wrappers, **все** вложенные причины которых
также numeric. Negative LAPACK info, MethodError/BoundsError, empty/mixed
wrapper и arbitrary ErrorException не принимаются. `inverse_domain`
оборачивает только retarded-state barrier; broad map/API catch не появился.
Обычный Infeasible от μ/capacity и ранее имел свой путь. Existing outer
catch даёт один counted convex fallback только rejected AA trial,
с очисткой history и streak; linear/fallback failure честно останавливает
attempt. Новые small controls покрывают numeric и mixed/programming
wrappers. Это restrictive recovery policy, не попытка repair matrix.

**F4 KMS primaries.** `kms_metrics:479–506` сохраняет full-domain joint
L/G numerator/scale/relative, weighted wE*wk counterparts, maximum-block
absolute/scale/locator и maximum-block relative. Знаки
Σ<−ifΓ и Σ>+i hole Γ правильны, denominator joint native Σ</> scale.
Zero policy прежняя, additive floor нет. Named scope честно ограничен
initial seed/current μ; independent stencil witness не выдан за uniform
theorem будущих Γ/μ. Все scalars считаются на уже имеющихся arrays.

**F5 per-channel/component audit.** `component_metrics:839` использует
norms(fresh,stored) для R/L/G: numerator fresh−stored, denominator fresh.
`audit_native:908–921` включает все пять local channels и embedding
total/plus/minus. Absolute/weighted primaries сохранены, diagnostic
explicit zero policy не подменяет owning native max rΣ с прежним floor.
Retarded/kinetic identity и signed collision metrics по этим families
сохранились. Full fresh candidate map повторно ради diagnostics не строится.

## Остальная delta и invariants

AA QR теперь сохраняет firstpair, R diagonal и initial weighted column
scales; формула real weighted type-II не менялась. Doc согласован с code:
это **nonmonotone feasibility-safeguarded** AA, growth raw residual сам
по себе не вызывает rejection. Не введён невыполняемый merit threshold.
Negative coefficients по-прежнему требуют cone/domain admission и fresh Π.

FDR output отдельно называет raw stored-total GRΣ</>GA и published paired
native Green. Full/weighted joint L/G, occupied-only lesser и trusted
weighted norms имеют явно записанные denominators; occupied norm не
размывается large empty spectrum. Trusted mask влияет только на diagnostic
accumulation. Diagnostic μ fit не меняет native arrays, raw λ или rK.
Additional A-normalized ratio явно именован отдельно от FD-normalized ratio.

BitArray metadata материализуется в dense Int8, как записано в model
schema; copied algorithms TOML добавляет provenance actual defaults.
Это I/O representation fix, не изменение grid mask/physics.

Native accounting теперь различает admitted/started materialization,
completed materialization, entered solver, observed completed iterations,
started/completed/report-written audit. На interrupted solver observed
count не объявляется exact started count: сохраняются lower/upper bounds
с reserved solver cap30, общий structural cap32. Source entrypoints явно
фиксируют 1+30+1; callback по-прежнему вызывается native после measured
candidate, до следующего mixing. Эти counters относятся к whole native
map attempts, не к числу отдельных BLAS/cavity/Dyson kernels.

Π equation, Γseed=2ηseedI, stable holes/current μ, local-only Γ,
static4 K[A], native LO Γ, same uniform PI, η=0 и same-input cavity
comparison не изменены. Native handoff по-прежнему заново materializes
five physical families/embedding, вызывает native raw Keldysh+paired
constraint и returns raw λ; history empty. Затем native linear α=.2,
strict max30/required3passes, без FD forcing, и одна fresh final audit.
Никакие reduced/μ-fit passes в native history не сконструированы.

Остаются ограничения V1 math review: declared PSD backward allowance
не является exact PSD/causality certificate; actual Float64 LO KMS
и consumer equivalence проверяются только первым actual precheck;
малый total Γ residual не доказывает weak-channel/native kinetic solution.
Mandatory final native raw rΣ/rK/λ/candidateλ/PSD/causality audit не
заменяется trigger. Small controls или bounded runtime success не
принимают Poisson, discretization adequacy и experiment validation.

Готовность к первому pilot относится только к проверенной mathematical
delta и имеющемуся small-suite evidence. Resource admission/overall
review/старт процесса остаются за ведущим; actual scientific outcome
ещё **not_measured**.
