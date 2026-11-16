# Независимый static review: frozen Γ128 uniform refinement control

**Scientific path статически пригоден для bounded two-map diagnostic.**
Новых блокеров не найдено; root должен завершить tiny v2 и собственный
OS resource admission до actual. Это не SCBA solve/acceptance, не
continuous-number oracle и не доказательство причины stagnation.
Reviewer Julia/Pkg/source HDF5 arrays не запускал/не читал и harness
не редактировал. Читались frozen sources, малые TOML и ранее сохранённые
endpoint CSV scalars. Применён qcl-physics-audit.

## Freeze/source bindings

Все пять SHA256SUMS entries независимо пересчитаны. Current freeze:

| Artifact | SHA256 |
|---|---|
| run_control.jl |5084c1f1f97e6b90bb8a369e3a212d6e68f94f9930a2aaa467406728d7d484cb|
| operators.jl |7d71dadce03edc17bcd7c68bae2fd5d1057e17cf95aae66432859ce4caa4a4dd|
| tiny_checks.jl v2 |565a8d182372c0d5cfab45f7598df10a0f593c14cfd39445140aa69438e33120|
| write-only input |ef917dcbb6cff39ceaf84aa5fdc1b53413dd733bdd58e01d32da984fd949d6e5|
| protocol.md |d0c3f35ce7d4926c0ff2d56e9fe96d248f15fa3aeb55682d84da68a9ef4287fa|

Source map/input/identity/gamma_functions hashes independently match
input; included warm_seed/gamma_restart hashes match source identity.
Large model/terminal SHA совпадают с previously audited cont128 bindings;
reviewer повторно их не читал. Actual code проверит все declared SHA,
source Julia version, imported Core e171/clean, coarse rebuilt model
grids/profiles/basis/h/K/qK exact equality, terminal shape и μ до карты.
Actual environment path/Manifest SHA сохраняются в output; они не
принимаются за scientific validation. Runtime input остаётся отдельной
root-owned копией; исходный `runtime_authorized=false` сохранён.

Include parent грузит definitions; actual CLI не вызывает parent128-loop,
evaluate/native compare, native full solver, mixing/AA или Poisson.
Только два syntactic `one_map` call sites. started увеличивается до
cone/retarded work и жёстко ограничен2; completion только после return.
Неудачный start также расходует cap, retry path отсутствует.

## API/cache и native guard

Проверены actual accepted Core constructors/fields и QCLNumerics namespace:
NumericalParameters positional constructor, ModelGrids/NEGFProblem,
ProductionCache, private kernel/PIFFT constructors, source/shift contracts,
validated cache, seed number weights и static/LO/retarded primitives.
Private cache точно повторяет owning public fields и operator builders;
full estimate берётся из `estimate_production_memory(problem)`, который
использует те же detected mechanism counts/compact-shift/PIFFT flags,
что post-detection public builder.

Coarse public/private equality проверяет все shift/kernel fields и
**все entries** двух actual FFTΛ evaluations на source Γ до fine.
Fine problem/cache/cavity принадлежат новому grid object, используют
новые sparse shifts и uniform PIFFT, прежние microscopic K/basis/h.
В cache.estimate остаётся настоящий full native estimate; options budget
не повышается. Fine `_validated_production_cache` обязан сохранить memory
refusal, если full estimate превышает cap. Его отказ проверяется отдельно,
а actual native primitives сами не вызывают full-solver validator.

Это отдельное admission map-only workload, не заявление, что24GiB full
solver помещается в12.125GiB. Core/default builder/guard не редактированы.
Tiny v1 root обнаружил namespace error после7 tests; v2 patch только
меняет `with_production_options` на `N.with_production_options` в tiny.
Это корректное фактическое место helper. Run/operators unchanged;
V1 archive и patch сохранены. Reviewer v2 runtime не выполнял. Его17
assertions проверяют cache/guard/nesting/prolongation/zero-edge HPI,
не source Φ или convergence; synthetic HDF fixture не читает source HDF.

## Grid, frozen Γ и μ/map stages

Coarse30704→fine61407=2N−1 строится native `build_grids`, те же endpoints.
Fine odd ε обязаны bitwise равняться source; остальные axes/weightsκ
сверяются exact. K,profiles,basis,h сохраняются, energy weights/trusted
mask/shift/cavity пересоздаются согласованно. Declared LO939→1878
следует из half step и той же physical36.7meV; actual minor interpolation,
KMS и weighted-adjoint witnesses сохраняются, integer snapping нет.

`prolong_rows` bounds корректны и для even последнего row chunk:
считывается дополнительный coarse right neighbour. Odd values копируются;
even используют один actual-coordinate scalar θ∈[0,1] для всех complex
matrix entries/k. Γ не Hermitize/clip/rescale. Convex prolongation
сохраняет piecewise-linear interior в exact arithmetic; floating source
и endpoint model оговорены отдельно. Input/fine A/Φ cone проверки
наследуют source backward-error policy, не превращают её в strict PSD
certificate и не снимают прежнюю Mats/causality qualification.

Каждый `one_map`: fresh native retarded state → spectral A → q_e из
GR с **source g_s·wE·wk/(2π)** → N(source μ) diagnostic → own current μ
fit → exact FD pair → один raw five-channel Φ. Fixed-μ N вычислен до
refit. Fit→target не называется convergence. Raw map residual denominator
есть raw Φ; Φ не заменяется Hermitian working map. Fine defect включает
quadrature/Dyson/embedding/current-μ эффекты, не изолирует один LO канал.

Coarse fresh raw Φ должен replay saved raw_Phi ≤1e−10 до fine.
Control D_f=Φ_f(PΓ)−PΦ_c считается на **всех fine rows** с reference
PΦ_c и fine weights; отдельно все even rows и old-node restriction
с referenceΦ_c/coarse weights. Selected CSV не подменяет full norms.
Anchor ±LO ladders/nearby rows покрывают интересующий support, но не
обещают обнаружить все новые peaks. Fixed .1E0 и Nb5 в CSV соответствуют
exact source model binding; это не общий arbitrary-model harness.

## HPI/ghost meaning

Одинаковые public/private coarse FFT values — arithmetic control.
Fine HPI(PΓ) old-node comparison использует coarse measure; ≤1e−10
gate проверяет выбранный source, а не универсальную invariance theorem.
ВнеΠ выполняются3 HPI operator calls; каждая с native doubled-padding
verification, так что это не утверждение «только3 низкоуровневых FFT».

Ghost taper ширина h→h/2 — **явный boundary variant**. Для E внутри
physical window bound `log(2)/(2π)(‖Γfirst‖F+‖Γlast‖F)` корректен:
разность left tapers даёт maximum integral log2 при E=a, right аналогично.
Runtime выражение с views на все k даёт uniform-in-E **combined-k
Frobenius** bound. Это не global-energy L2 norm, relative HPI bound
или bound на resolvent/Φ amplification. Такая размерность должна быть
сохранена при публикации summary scalar `ghost_model_change_bound_native`.

Независимо из saved perE endpoint scalars и triangle inequality
‖Γedge‖≤√raw_ref²+√raw_residual² получается тот же upper bound
5.159247206455558e−20 source-energy units. Runtime запишет actual
endpoint-norm bound; маленькая величина не разрешает объявить весь
continuous Γ/ghost model неизменным или объяснить fine-map outcome
исключительно interior mesh.

## Lifetimes/ресурсная область

Fine F=270190800B≈257.674MiB. Separate admission24F+3GiB=
9.039235GiB, hard12.125GiB. Наибольший явно прослеженный primitive
live set: Γinput+R+shared zero dummy + cavity plus/minus/total9F +GR
до≈13F; thermal LO существенно ниже24F. Histories=0, channels не
retained. Context/retarded/one_map function barriers не возвращают
Green/LO/cavity fields; перед fine большие coarse Γ/Φ/Λ references
сняты и GC выполнен. Fine и coarse cache/context кратко пересекаются
при rebuild, но без coarse full fields; общий K/basis ownership shared.
Cache/cavity actual payload ограничен .5GiB. Dead-but-not-collected
temporaries/JIT/HDF writer дополнительно учтены reserve; это static
live-set estimate, **не доказанный RSS bound**.

Resource gate читает cgroup memory.max/high/swap/cpu и actual file soft
limit; Julia1/BLAS1,CPU≤100%,hard13019119616B,Swap0 проверяются.
Soft840s проверяется между stages; external900s обязан ограничивать
также include/import/блокирующий primitive. Data3GiB/file2GiB и padded
atomic HDF writer контролируют output; существующий out отвергается.
Summary/failed outcome сохраняются в finally; external kill может
оставить последний checkpoint/.tmp и не обязан пройти finally.

Scope/admission MemAvailable/disk/no-other-heavy остаются обязанностью
root. No extra map/native exit допускается даже при экономии ресурсов.
Measured tiny/RSS и actual two-map outcomes пока не приняты этим
статическим review; после них вывод только о frozen discrete collocation
sensitivity, без SCBA/mesh adequacy/causal-stall acceptance.

## Root tiny v2 evidence и current decision

Прочитан root `tiny-v2-root.log`:17/17 PASS,exit0,14.36s,
peak818528KiB,Swap0; run5084/operators7d71 неизменны, tiny565a.
Current SHA256SUMS повторно verified all match после one-line fix.
Это закрывает actual tiny API/cache/nesting fixture, не replay source
HDF5/Φ, fine source-grid admission или two-map scientific outcome.
Reviewer runtime не выполнял. Для последующего отдельного actual
two-map window статических блокеров нет; actual пока не запускался,
согласованный root план не допускает overlap с другими heavy/oracle
qualification. Успешный tiny не расширяет2-map/resource cap.
