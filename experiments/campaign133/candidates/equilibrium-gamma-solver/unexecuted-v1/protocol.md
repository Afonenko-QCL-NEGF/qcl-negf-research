# Equilibrium-reduced Γ SCBA: pilot16 implementation protocol

**Статус: IMPLEMENTED, NOT EXECUTED.** Root разрешил только implementation
harness в этом каталоге. Julia/Pkg/compute не запускались. Core/Runner и их
frozen worktrees не менялись. `pilot16.toml` имеет `runtime_authorized=false`;
после static review Root отдельно разрешает runtime и передаёт accepted
Core environment. Запуск не следует из наличия этих файлов.

Accepted Core: `e171caa60cce676c6c4238dd671e7c527f2bd977`, tree
`742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8` (PR30, Bv2+C). Фактически
импортированный path/HEAD/tree/dirty и hashes файлов/входа проверяются при
запуске и сохраняются; source path задаёт штатная Julia environment Root.

## 1. Какой вопрос и какая identity

Проверить, помогает ли устранение thermal occupation/отдельных channel-history
degrees of freedom найти **равновесный fixed point** исходного пятиканального
finite-window PI map. Решается один LOCAL irreducible scattering broadening
field Γ(E,k) размером NE×Nk×Nb×Nb, не total broadening вместе с embedding.
Boundary embedding вычисляется отдельно из local self-energy.

All-five S01: V=0, TL=TLO70K, LO/acoustic/impurity/IFR/alloy; e–e none;
fixed screening, parabolic dispersion, unchanged physical kernels/donors;
Nz96/Nb5/Pbasis3, Nk11/kmax.60nm^-1, Nφ32/Nqz81/qzmax10nm^-1;
finite_chain4 с каждой стороны, U=0. Это fixed-Hartree этап, не Poisson solve.

Grid **m939**: Ω=.0367eV, Emin=−.3eV, NE30704,
Emax=.9000001064962726eV, h≈39.0841320554μeV. Actual nodes и native stencils
сохраняются; никакого snapping, wrap или округления Ω. Это отдельная
grid identity, не в точности прежнее окно [−.3,.9]eV/NE30721.

Retarded reconstruction: **uniform product integration с existing exterior
zero ghost hats**, `retarded_real_part=:kramers_kronig`,
`elastic_retarded=:hilbert`, full matrices/resolved k/direct kernels/dense BLAS.
Используется accepted Core PR30 с uniform PI FFT Bv2 и C selector,
который явно **hilbert**. C `direct_green` не эквивалентен этой
Γ-only редукции; никакого silent выбора C по наличию API.

Метод обновления Γ — **custom Gamma Anderson-II**, отдельная execution identity.
`ProductionOptions.algorithms.mixing=:linear` нужен только существующим
operators/cache; он не выдаётся за реально выполненный nonlinear algorithm.
Native full-state Anderson не вызывается и его history не создаётся.

## 2. Π(Γ): retarded state, current μ и reduced map

Предпосылки, подтверждённые equilibrium_oracle: Γ Hermitian/PSD, PI matrix
одна и та же для всех columns/channels, V=0, один common μ/T, thermal LO,
actual static K комплексно-линейны/adjoint-preserving и exact LO shift
согласован с Ω. На actual Float64 stencils последнее проверяется измерением,
а не объявляется exact из имени m939.
Полный conditional proof/domain/causal-span разбор:
`../pi-causality-window-certificate.md`, AppendixA1–A6. Нижеследующий протокол
сверен с ним после получения; proof не заменяет actual-stencil/native checks.

Один evaluation (без смешивания families внутри):

1. Проверить Γ: finite, Hermiticity и full-matrix PSD по каждому E,k.
   Нет eigenvalue clipping/η floor. Если нарушение больше native backward
   arithmetic budget, trial отклоняется. Записать minimum eigenvalue,
   norm/budget и locator худшего block.
2. Local ΣRΓ=H_PI[Γ]−iΓ/2. H_PI задаёт causal piecewise-linear positive
   matrix spectral density с тем же finite support/ghost taper; ни density
   normalization, ни дополнительный seed imaginary term не добавляются.
3. Native cavity4 из **этой local ΣRΓ**, затем native Dyson центрального
   full GR с η=0 и A=i(GR−GA). Для R-only вычисления native cavity API
   может получить временные local lesser=greater=0: cavityR от них не
   зависит. Эти нулевые fields — computational placeholders, не physical
   Keldysh family, не raw state и не scientific zero; они сразу освобождаются.
4. По current A найти **новый** общий μ, N[fA]=NDsheet, при T70K.
   Проверить positive spectral number weights, bracket/capacity
   0<ND<C[A] и finite root. Записать μ, Nfit/target−1, capacity/target,
   dN/dμ=Σ q_e f_e(1−f_e)/kBT и root bracket/termination.
   Γ PSD не гарантирует capacity/root feasibility. Использовать μ прошлого
   trial или подменять failed root density scaling нельзя.
5. Построить GFD<=ifA, GFD>=−i(1−f)A; малое hole filling считать stable
   complementary logistic, не вычитанием почти единичного f, и сохранять
   measured f+hole−1/Green-identity error. Это явная thermal restriction
   custom map, не native paired-normalized state. Она действует только в
   reduced stage и обязательно записывается в identity.
6. Static4: Γ_s=K_s[A] с прежними qK,wk. LO: existing native
   `_production_lo_contraction(GFD)` → ΣLO</>; ΓLO=i(ΣLO>−ΣLO<).
   Φ(Γ)=Σ_s Γ_s. Не добавлять embedding в Γmap и не применять Hilbert
   к kinetic input повторно. LO и e–e rules не переписываются.

Выход Π: GRA/GFD/currentμ, raw Γmap/per-channel broadening, feasibility
witnesses и raw **unmixed** F=Φ(Γ)−Γ. Cavity from fresh Γmap не обязан
совпадать с cavity from inputΓ до fixed point; не требовать это как control.

В точном equilibrium-invariant случае Φ(Γ*)=Γ* позволяет materialize
channel Γ_s, ΣRs=H_PIΓ_s−iΓ_s/2 и Σs</>=ifΓ_s/−i(1−f)Γ_s.
Линейность PI даёт сумму retarded local families; cavity FDR и number fit
дают native Graw=GFD, λ=1. Это доказательство применимости редукции, **не**
доказательство сходимости итераций или uniqueness thermal ensemble.
Actual Float64 KMS/edge/PSD и strict native audits остаются решающими.

## 3. Anderson-II на одном Hermitian field

F_j=Φ(Γ_j)−Γ_j, ΔΓ_i=Γ_{i+1}−Γ_i, ΔF_i=F_{i+1}−F_i.
History **m=3** (максимум4 past states/residuals); m5 не включается
автоматически. Все scalar coefficients вещественные. Working inner product:
Re Σ_E,k wE wk Tr(X†Y), положительные weights; сохранять также native-style
unweighted residual, поскольку последний нужен для сопоставления с Core.

Решить min_γ||F_j−[ΔF]γ||_W; candidate

```
ΓAA = Γj + βFj − ([ΔΓ]+β[ΔF])γ,  β=1.
```

Это difference form Anderson-II. Не использовать normal equations с молчаливым
изменением кондиционирования: two-pass modified Gram–Schmidt/QR на максимум3
разностных columns, маленький R solve; dropping oldest history при rank loss
или cond(R)>1e8. Сохранять singular/diagonal scales, retained columns и
coefficients. Предел sum|affine coefficients|≤10 — заранее объявленный
algorithm guard, не scientific tolerance. Нет пригодной history → linear step.

1. AA trial сначала проходит full Γ Hermitian/PSD/nonfinite check.
2. Только после допуска строится **новая Π(trial)**; старые GR/μ не
   смешиваются с новым Γ. Если Π trial не имеет finite inverse/μ/capacity,
   отклонить trial, очистить history, попробовать **один** convex fallback
   Γlin=.8Γj+.2Φ(Γj). Этот fresh Π тоже потребляет evaluation/wall budget.
3. Если fallback inadmissible, остановиться с raw failure; convexity PSD
   не гарантирует μ feasibility. Никаких hidden diagonal regulator, matrix
   repair или бесконечного backtracking.
4. Сильный raw regression/лимит evaluations завершают bounded attempt с
   best/last data, не расширяют окно. Global convergence AA не обещается.

Hermitian coordinates: diagonal real, upper triangle complex, lower определяется
adjoint. Raw output Φ сохраняется до representation conversion; если его
anti-Hermitian defect превышает native arithmetic budget, stop. Допустимое
округление при packing/unpacking измеряется и записывается, eigenvalues не
исправляются. All-native terminal audit использует исходные raw operators.

Reduced stopping **только handoff trigger**: rΓ=||Φ−Γ||F/||Φ||F<1e−9
и weighted relative norm<1e−9 три последовательных valid evaluations;
обе абсолютные нормы, оба denominators и E0_eV сохраняются. Нет additive
floor: 0/0=0, positive/0=Inf. Сам Γ=0 отклоняется как inadmissible input,
а capacity/root tests обязательны независимо от relative norm;
τΓ строже native self-energy target1e−8, но не подменяет channel-relative
rΣ или rK. μ/N/PSD feasibility обязательны. Малый mixed step не используется.

[Walker & Ni, SIAM J. Numer. Anal.49 (2011),1715–1735](https://doi.org/10.1137/10078356X),
[проверенный primary PDF](https://users.wpi.edu/~walker/Papers/Walker-Ni,SINUM,V49,1715-1735.pdf):
AlgorithmAA/Eq1.1–1.2 задают constrained residual least squares и relaxation;
Sec3 TypeII и Sec4 — difference formulation, QR/conditioning/history removal.
Их линейные GMRES equivalence results не распространяются автоматически на
наш нелинейный PSD/μ-constrained map. Указанные safeguards и limits — наш
explicit variant, не цитируемая из статьи theorem для QCL.

## 4. API bindings (используются, а не переписываются)

Символы доступны через existing `QCLNEGF.QCLNumerics`/root imports;
private-port harness остаётся bounded campaign artifact, не новым public API.

| Действие | Реальная binding |
| --- | --- |
| Исходная модель/actual grid/kernels | `reference_parameters`, `NumericalParameters`, `build_configured_scattering_problem`, `project_hamiltonians(problem,zeros(Nz))` |
| Cache/PI plan | `build_production_cache(problem;options=po)`, где po hilbert=:product_integration, mixing=:linear, workers4/hilbert_columns1; проверить cache.product_hilbert_operator именно `_ProductionPIFFTHilbertPlan`, не dense H |
| PI Γ→Λ | `_production_pi_fft_hilbert_transform(Γ,g.ε,g.wᴱ,cache.product_hilbert_operator)`; returned roundoff отдельным scalar. Не подавать Γ какΣ<, если можно вызвать этот owning port напрямую |
| Cavity R | `_cavity_production_plan(problem,h,po)` once; `_cavity_embedding_production(problem,local_family,h,po;plan)`; extract R и освободить placeholder L/G arrays |
| Dyson/A | `_retarded_green_production(g.ε,h,localR+embeddingR;η=0,options=po)`, `spectral_function(GR)` |
| Common μ | `_seed_number_weights_production(GR,g,gs;options=po)` + `_marker_spectral_chemical_potential` либо existing bracketed root on `_seed_number_from_weights`; check returned status/actual root, no repeated array allocations per μ |
| GFD container | explicit ifA/−i(1−f)A, `GreenState` or `_green_state_production(...;greater=GG)`; no `_normalize_keldysh_pair!` in reduced stage |
| Static Γ_s | `_selected_static_contraction(problem,cache,id,A,qK_s,po)` or preallocated `_production_static_contraction!` on actual cached operator; same axes/qK/wk |
| LO Γ | `_production_lo_contraction(cache.kernels[:LO],GFD,problem,cache,qKLO,po)`; actual sparse ±Ω plans |
| Independent native compare | `_scattering_candidate_production(problem,GFD,cache,po;return_roundoff=true)` + `_embedding_self_energy_selected(...;scattering=native,h)` |
| Native exit | `_current_green_production`, `_solver_restart_contract`, `SCBAResult`, `solve_scba_production`, `fresh_fixed_point_metrics` with actual native candidate |

Перед runtime assert actual `pathof`, HEAD/tree/dirty patch, new approved PI
source/ref, exact algorithms/identity, model IDs, all five kernels and no e–e.
Current C/Runner refs не перепиниваются автоматически этим design.

## 5. До nonlinear loop: cheap crosschecks и области их независимости

**A. Новые small checks** (не меняют S01): CP complex2×2 rational GR и
independent explicit VXV† loops against actual literal + production static
contraction; ΓR=K[A]. Для custom engine: weighted QR/type-II vs independent
dense real least squares; nonlinear positive scalar fixed point; negative,
non-Hermitian и zero Γ rejection; analytic one-level μ, infeasible capacity,
stable occupied/hole factors; tiny HDF hyperslab roundtrip для2/3/4/6D и
16-call budget rejection. Это выполняется `small_checks(cfg)` до model
build. Новые small native LO/full-chain/PI fixture в harness не дублируются:
они опираются на existing accepted Core focused tests и frozen PI control.
Actual five-channel/cavity consumer проверяет следующий precheck.

**B. Actual S01 m939 before SCBA:** одна Π(Γseed), один full-native
five-channel fresh candidate на его GFD. Compare Γs_reduced against
i(Σs>−Σs<), native sumΣRs against PI[Γmap]−iΓmap/2, full raw KMS and Γ
Hermitian/PSD. Это actual consumer equivalence, не независимый oracle самих
общих native contractions. Проверить также cavityR equality для **одной и той
же localΣR** при dummyL/G0 и physicalFD localL/G currentμ, а не междуΓinput
иΓmap. Полную native fresh embedding построить из native fresh families,
сохранить отдельно от input cavity.

Arithmetic target comparisons≤1e−10 relative с явно сохранёнными numerator/
scale и explicit zero policy и matrix maxima; physical KMS expectation сопоставить с independent
actual-stencil coefficient prediction из `commensurate-grid-protocol.md`.
Если capacity/inverse/KMS equivalence не проходит, nonlinear loop не начинается.
Новый tolerance по observed failure не подбирается без новой identity.

Seed: local constant Γseed=2ηseed I =2.5meV I (=.025I scaled E0=.1eV).
Это **начальное поле неизвестной**, не persistent η в Dyson. Никакой
неявной загрузки непринятого S01 state вместо seed. Γ=0 не используется,
поскольку represented capacity может быть нулевой при промахе poles.

## 6. После reduced trigger: mandatory unprojected native validation

1. На последнем GFD получить все **пять native** families через existing
   candidate routine и fresh finite-chain embedding±/total. ChannelR uses PI,
   directGreen выключен. Не объявлять Γ sum decomposition найденным solution
   без этого шага: individual relative residuals могут усиливать малые channels.
2. `_current_green_production(problem,h,scattering,embedding,po)` строит
   **native** GR и independently evaluated raw Keldysh pair + existing paired
   constraint, возвращает raw λ/total. После этого μ **не подбирается**, FD
   не навязывается и fields не заменяются fA. Так materialized Green aligned
   со stored native families, даже если reduced trigger был лишь приближённым.
3. Сохранить native warm state `converged=false, quality=:unresolved`, history
   empty, actual `_solver_restart_contract(strict_options,po.algorithms)` и
   empty **linear** mixer state. Не записывать custom Γ history как native
   exact Anderson continuation. Дальше public `solve_scba_production(initial=...)`
   с max30, strict mode/disabled approximate quality/required3 passes,
   αΣ=.2. Общий native cap32 включает 1 materialization map +≤30 solver maps
   +1 fresh final audit map; это не32 iterations сверх materialization/audit.
4. Обязательны fresh rK/rΣ/λ/candidateλ, rD/rA, PSD/causality, per-channel
   kinetic/retarded identity, weighted signed collision/power и full
   raw/normalized FDR с independent commonμ diagnostic. FDR diagnostic теперь
   только наблюдает native state и ничего не переписывает.
5. Только native SCBA strict result можно принять как fixed-U inner
   convergence. Poisson/history отсутствуют: полного NEGFSolution.converged,
   discretization adequacy/experiment и unique thermalization не заявлять.
   r_sum/window/basis/k/embedding sensitivity остаются отдельными gates.

Если Γ trigger не достигнут или native strict stage fails, сохранить оба
исхода отдельно. Reduced convergence никогда не получает старое scientific
baseline имя и не превращается в native pass задним числом.

## 7. Замороженный budget и память

VM20.6GiB, available около14GiB (не прежняя оценка25GB). Один process,
Julia4/BLAS1, heap hint8G, CPUQuota400%, `MemoryMax=13287555072`B
(**12.375GiB**), `MemoryHigh=12348030976`B (11.5GiB), SwapMax0;
Root ставит штатные systemd limits, RuntimeMaxSec1800, LimitFSIZE5GiB,
TERM→KILL10s. Source gate читает actual cgroup memory.max/high/swap,
cpu.max и /proc/self/limits, отказывается при несогласованном admission.
Новых heavy после05UTC нет. Параллельного heavy physics job нет.

- Prechecks≤180s, ≤6 actual calls. Actual reduced seed Π + full native
  candidate + physical same-input cavity считаются как3 actual calls;
  внутренние primitive contractions не переименовываются в бесплатные maps.
- **Всего≤16 Π calls**, включая precheck seed и все rejected/infeasible
  Π trials/fallbacks. После seed остаётся≤15. Reduced stage≤600s; нет
  автоматического расширения до128/192. PSD rejection до Π сохраняется как
  rejected trial, но вычислений Π не выполняет.
- Native stage≤900s и≤32 native maps (1+30+1), только по trigger. Никакой
  μ/FD projection в native stage; posthoc FDR fit только диагностический.
- Общий wall≤1800s с момента harness budget; внешний systemd cap охватывает
  также import/JIT. Числовой call budget проверяется ДО начала primitive/map.
  Внутри native public loop wall check запрашивает штатный aligned checkpoint;
  callback сохраняет/возвращает borrowed state посредством специального
  exception, не возвращаясь к in-place mixing. Stage wall может быть превышен
  на один уже начатый native primitive/map; внешний OS cap остаётся конечным.
- Один файл≤5GiB, все данные pilot≤6GiB, stdout finite у Root. Writer проверяет
  actual bytes; LimitFSIZE ограничивает один файл независимо от кода. Root
  до запуска проверяет disk headroom≥6GiB+margin. Старые raw не удаляются.

Один ComplexF64 field =135097600B≈128.84MiB. m3 history — максимум4 Γ/F
pairs, 8 fields≈1.01GiB; QR≤3 columns≈386.5MiB. Дополнительно живут current
Γ/Φ/residual/trial/best, Green4, localR/dummy/cavity temporaries, LO buffers,
cache, JIT/GC и HDF writer. Нет утверждения, что весь процесс занимает один
field или менее3GiB. Примерный Γ-map live-array ledger — порядка25–35 fields
(3.1–4.4GiB) плюс Julia/native/cache; actual RSS до исполнения not_measured.
Precheck временно хранит5 channelΓ +15 native Σ fields, зато history/QR ещё
отсутствуют. QR и map stage не держат одновременно две full-state histories.

Native full-state Anderson с4state/residual pairs потребовал бы168 fields
только истории (~21.1GiB), поэтому не используется. `mixing=:linear` cache
policy явно отличается от custom Γ-AA identity.

Core `ProductionOptions.memory_budget_bytes` равен реальному hard cap, guard
**не обходится**. Его native full-solver estimate не складывается слепо с
Γ-history: эти stages имеют разный live set. Root сравнивает собственный
Γ-stage ledger и native estimate с OS cap/reserve. Source guard refusal —
honest admission stop, не повод выдать Core фиктивный budget16GiB.
`reduced_stage` возвращает только Γ и small scalar summary; все histories,
QR scratch и old Green находятся за function barrier. Затем `GC.gc(true)`
до `native_stage`; warm/native scratch тоже освобождается перед final audit.

## 8. Артефакты и ограничения реализации

`gamma_functions.jl`: typed `Context{P,C,O,H,CP}`, `Field=Array{ComplexF64,4}`,
Π/μ/PSD/AA, raw diagnostics и native handoff. Большие loops находятся в
специализированных функциях, не в top-level main с abstract Dict values.
`checks.jl`: small behavioral/independent controls, не выполняются при include.
`run_gamma.jl`: gate → checks → build → precheck/loop → function barrier →
conditional native stage. Никакого Pkg/Runner/new public framework.

Каждый map имеет immutable map-NNN.toml (μ/capacity/dN_dμ, absolute/relative
weighted/unweighted residuals, cone/budget locators, channel norms, cost),
каждый accepted/rejected attempt отдельный TOML с AA coefficients/condition.
Γ-first/intermediate(call8)/best/last — примерно129MiB каждый, только Γ с
metadata и recipe fresh-history continuation, **не native exact restart**.
Actual-one-map raw содержит rawΦ/nativeΓsum/nativeΣRsum (~387MiB); scalar
control хранит per-channel KMS/norms, actual shift coefficient witness и
weighted-adjoint defect. Original stencils/weights/basis/h/kernels сохраняются
в model.h5. Raw failed PSD Φ/trials сохраняются без clipping.

При cap/trigger сохраняется reduced-terminal.h5: Γ +rawΦ +GR/GL/GG/A (~773MiB).
При domain failure перед завершённым last-map может отсутствовать last G;
raw failed Γ, предыдущие Γ checkpoints и failure reason остаются. Повторная Π
ради output не выполняется сверх cap. Warm native private state не дублируется
на диске: конечный native-terminal.h5 содержит все Green4+5families15+
embedding/±9 fields (~3.52GiB), restart contract, status и native history CSV.
Native fresh candidate массивы не дублируются вторым3.52GiB файлом: их raw
residuals/collision/FDR/identity diagnostics сохраняются, а точное повторение
возможно из frozen source+stored G+model. Это ограничение data cap явно.
Chunked HDF hyperslabs имеют writer buffer≤~8MiB. Admission до записи
проверяет raw array byte estimate+metadata reserve против file/data caps,
после записи — actual bytes. Если rejected-trial evidence заняло слишком
много места для полного native output, handoff получает честный admission
stop до native allocation. Имеющиеся artifacts не перезаписываются; tmp
при failure сохраняется. Signed collision означает incoming−outgoing,
Tr(Σ< G>−Σ>G<), measure wE*wk/(2π); spin/ħ/physical-rate conversion здесь
не применяются, conservation сравнивается с absolute incoming/outgoing scale.

`scientific_acceptance=false` остаётся default summary: возможный native
strict inner pass требует независимого разбора fresh audit и не принимает
Poisson/discretization/experiment. Если native audit не выполнен из-за wall,
он получает `not_measured`, не pass. Теорема Γ reduction/Walker–Ni не обещает
сходимости этого bounded nonlinear attempt.

До runtime: Root independent static review frozen code+pilot16, resourcegreen,
accepted environment с Core e171/tree742; затем отдельная разрешённая копия
input меняет только `runtime_authorized=true` и получает собственный SHA.
Никаких Core/Runner implementation/source-pin edits этим harness.
