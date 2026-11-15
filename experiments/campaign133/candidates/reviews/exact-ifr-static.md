# Независимый review exact IFR structured contraction

Текущий статус: **frozen implementation/test-v2 и same-field benchmark приняты
в узкой области Core IFR optimization**. 457 current distinct assertions pass;
measured IFR apply median speedup2.196× при all-entry error порядка roundoff.
Итоговый result/границы — `reviews/exact-ifr-benchmark-result.md`. Exact freeze
и evidence ниже; прежние pending статусы сохранены как история review.
Accepted source прочитан на Core
`e171caa60cce676c6c4238dd671e7c527f2bd977`, worktree
`.build/campaign133/core-accepted-science` clean. Новый isolated worktree при
начале review ещё не существовал; автор затем сообщил создание
`.build/campaign133/core-ifr-separable` и продолжающуюся реализацию.
Окончательный code verdict будет дан только на frozen ref. Julia/Pkg/runtime
не запускались, source не редактировался; **runtime_not_measured**.

Прочитаны `adjacent/exact-elastic-contractions.md`, `adjacent/exact-ifr-design.md`
и фактические canonical/reference builders, normalization, flattening,
generic cache dispatch, selected consumers и memory estimator. Физическое
утверждение проверено по actual implementation; публикации в author memo
не использованы как замена этой привязки.

## Алгебра на accepted коде

Reference `kernels.jl:162–225` задаёт interface interpolation, seam z=0,
signed band jump и
`D_j[a,c] = d_j*conj(χ_j[a])*χ_j[c]`, `d_j=ΔEc_j/L0` real.
В optimized `kernels.jl:458–494` тот же covariance суммируется один раз,
а Gaussian finite-φ mean одинаков для всех interfaces и отдельно зависит
от (m,m′). Общие Δ/Λ в PhysicalParameters здесь действительно общие;
никакой continuum-angle/Bessel замены для тождества не требуется.

Положив `u_j=conj(χ_j)`, имеем `D_j=d_j u_j u_j†` и

```text
Φ(X)_m = Σ_(m′,j) w_m′ S_mm′ D_j X_m′ D_j†
       = Σ_j d_j² u_j [Σ_m′ S_mm′ w_m′ (u_j† X_m′ u_j)] u_j†.
```

Равенство комплексно-линейно для **любого** complex X, не только Hermitian
или PSD. В basis-fast pack
`R[j,(c,d)]=conj(u_j[c])*u_j[d]=χ_j[c]*conj(χ_j[d])`,
`L[(a,b),j]=d_j²*u_j[a]*conj(u_j[b])` с согласованным размещением общего
spectrum/normalization scalar. Ошибка transpose вместо adjoint здесь
изменила бы карту на complex off-diagonal inputs.

Accepted `_production_kernel`, `production_cache.jl:369–407`, имеет rows
`(m,a,b)`, columns `(m′,c,d)`, fastest m. Поэтому явная полная matrix была бы
`B_basis ⊗ S_radial`; при basis-fast packing меняется порядок Kronecker
factors. Трёхстадийное применение не должно молча reshape-ить один порядок
в другой без packing/permutation.

Это сумма rank1 **interface vertices**, а не утверждение global rank1 всего
IFR superoperator. При физических `w_m′>=0`, `S_mm′>=0` карта имеет Kraus
форму с `sqrt(w_m′ S_mm′)D_j`, поэтому CP между direct sums momentum blocks.
PSD radial matrix S как matrix не требуется — достаточно nonnegative entries.
Trace preservation или convergence SCBA из этого не следуют.

Normalized source контракт: `Khat=Kraw/qK`, selected contraction умножает на
external qK ровно один раз. qK является положительным scalar, а не новым
angular/radial weight. Дополнительных wE, spin factors, energy remapping или
FD projection в этой оптимизации быть не должно. Signed d исчезает лишь
в quadratic covariance, не в произвольном физическом coupling.

## Сообщённые автору memory/eligibility challenges

**1. Existing dense estimate не автоматически верхняя граница.**
Предложенный scratch на worker равен
`Nk*chunk*(Nb²+2J)*sizeof(ComplexF64)` против прежних
`2Nk*chunk*Nb²*sizeof(ComplexF64)`. При `2J>Nb²` (например Nb2/J4)
новый scratch больше. Новый payload `2Nb²J+Nk²` complex entries также может
превысить прежний cache payload. Поэтому нельзя одновременно выбирать этот
path для всех J/Nb и оставить прежнюю оценку без eligibility/fallback.

Автор согласовал достаточное ограничение `2J<=Nb²` и payload guard.
Дополнительно отмечено, что baseline payload равен **Nb⁴**, если actual
`_is_momentum_independent(K)` true, даже при Nk>1: при малом Λ finite-angle
exp может округлиться в одинаковые значения. Безопасная простая priority —
оставить existing momentum-independent path раньше IFR specialization;
либо сравнивать payload с фактической reduced/dense baseline веткой.
Оба условия должны быть проверены в frozen code.

**2. Exact recipe admission не доказывает finite новых factors.**
qK должен быть finite positive; balancing/scaling L/R/radial может overflow
или underflow иначе, чем уже хранящийся normalized tensor. Требуется finite
factor guard и fallback к прежнему generic path, без NaN sanitization или
выдачи некорректного structured operator. Автор подтвердил план такого guard.

**3. Canonical admission должно проверить всю карту.**
Требуется exact elementwise comparison каждой normalized entry с прежней
canonical recipe `/source_qK`, а не tolerance, несколько samples или имя
mechanism. Custom/imported K, modified one-ulp entry и reference rounding
должны либо пройти буквально тот же guard, либо сохранить generic fallback.
Raw source K остаётся доступен как oracle/provenance. Exact Float64 recipe
matching допускает roundoff-equivalence применения, не bitwise identity
нового physical-factor arithmetic со старым dense GEMM.

## Достаточные предстоящие behavior oracles

- Полный matrix-unit basis E_cd каждого input momentum на tiny fixture с
  complex χ, unequal positive wk, signed jumps и finite φ: independent
  literal covariance loops против generic dense и selected structured path.
  Одних real/PSD inputs недостаточно для transpose/conjugation диагностики.
- Arbitrary complex GR, complex-linearity, adjoint compatibility,
  PSD/extended-positive либо tiny Choi witness; `i(Φ(GR)−Φ(GR)†)=Φ(A)`.
  Отдельно normalized qK/α/β, complex α/β и β=0 с неинициализированным output.
- Actual canonical production builder против independent reference IFR
  builder на тех же φ nodes; common factoring helper не должен служить
  единственным oracle собственного refactor.
- Basis phase rotation, zero-height/empty-disabled branch, zero/signed jump,
  custom tensor one-ulp mismatch, generic standalone/custom fallback,
  explicit literal/low-rank paths и momentum-independent priority.
- Partial final chunk/reuse, weights, worker count/backend, payload bytes,
  memory guard counterexamples Nb1/Nb2/largeJ и multiworker private scratch.

Error следует хранить как absolute numerator и нормировать на input/operator
scale (weighted sum norm vertices/radial × norm input), особенно около zero
output. Нужны finite-arithmetic допуски, а не clipping PSD или ослабление
scientific acceptance. После tiny qualification отдельно разрешённый warmed
actual-shape timing нужен для speed claim; арифметический MAC ratio и IFR
norm не измеряют wall-time share. Γ/Core scientific trajectories пока не
обновляются и ускорение всей QCL кампании не установлено.

## Frozen code review: patch 234d984d…

Фактический base HEAD worktree `.build/campaign133/core-ifr-separable` —
`e171caa60cce676c6c4238dd671e7c527f2bd977`; source frozen как незакоммиченная
пятифайловая дельта +286/−5. Точный artifact `ifr-separable.patch` SHA256:
`234d984dea45410af8c28ed92c9b114114741295c5380fd4477f2db8b1aa8f5d`.
Hash и `git diff --check` независимо проверены. Прочитан полный source/test
diff и `adjacent/exact-ifr-handoff.md`. Рецензент код/Γ/accepted Core не менял,
Julia/Pkg/runtime не запускал.

**Static verdict frozen patch:** готова к ограниченным runtime checks;
блокирующих ошибок формулы, conjugation, packing/layout, qK normalization,
GEMM aliasing или сохранения memory upper bound не найдено. Это **ещё не
приёмка implementation/PR и не измеренное ускорение**.

### Builder и admission

`optimized/kernels.jl` выделяет `_interface_roughness_kernel_data`, сохраняя
прежние D, order interface covariance sum, `_pair_covariance_to_block`, K0,
finite φ order, max(q²,0), Gaussian exponential и amplitude/φcount arithmetic.
Первоначальный full raw tensor по-прежнему создаётся как radial×block,
нормируется прежней функцией и остаётся в `KernelSet.K`. Запись radial
scalar в отдельную матрицу вместо немедленного заполнения block не меняет
physical angle/boundary model. Добавленные χ vectors — именно `u=conj(χ)`.

`_canonical_ifr_production_kernel` выбирается только для IFR при явном
`contraction=:dense_blas`. Standalone generic `_production_kernel(K)`, literal
и explicitly requested low_rank dispatch не подменяются.

Admission проверяет shape, actual momentum independence, количество interfaces,
scratch/payload eligibility, qK finite positive, positive amplitude и **все**
normalized K entries против повторённой canonical source recipe. Comparison
содержит equality, не tolerance. Расхождение либо невозможность построить
recipe (`ArgumentError`/`DimensionMismatch`) возвращает generic fallback;
неизвестные programming exceptions не проглатываются. Наличие имени IFR само
по себе не даёт fast path. После построения L/R/radial все factors проверяются
на finite values; при отказе используется прежний path.

Canonical exact-entry matching квалифицирует source recipe. Оно не означает
bitwise equality reordered floating arithmetic factors и dense contraction;
для этого требуются предусмотренные roundoff-equivalence oracles.

### Нормировка, complex карта и packing

Constructor задаёт

```text
c_j = (ΔEc_j/L0)*sqrt(amplitude/K0)
L[(a,b),j] = c_j²/qK * u_j[a]*conj(u_j[b])
R[j,(c,d)] = conj(u_j[c])*u_j[d]
S = canonical_radial/amplitude.
```

Следовательно `qK*L*S*R` восстанавливает amplitude/K0, signed-jump square,
independent-interface covariance и finite-angle mean ровно один раз.
wk входит при radial packing ровно один раз; wE/spin/FD/energy shift не
добавлены. Реальный отрицательный jump после квадрата не меняет CP sign.
Формула сохраняет complex off-diagonal X; ни transpose(X), ни conj(X) нет.

Pack `P[(a,b),(e,m)]` явно меняет native m-fast field order на basis-fast
columns. После `T=R*P` weighted scalars пакуются как `X[m′,(e,j)]`, затем
`Y=S*X` и `F[j,(e,m)]` перед `P=L*F`. Scatter возвращает `(e,m,a,b)`.
Перестановки индексов согласованы с literal six-axis order `(m,m′,a,c,d,b)`.

`Y` использует memory scalars только после завершения чтения старого T;
`F` использует memory radial_input только после завершения GEMM Y=S*X.
Каждый GEMM имеет разные input/output storage. Partial final chunk views
охватывают одинаковый contiguous used prefix в обоих reshape headers;
неиспользованный tail не входит в multiplication. P переиспользуется после
чтения старого input. Private arrays создаются внутри каждого run_worker;
parallel chunks пишут непересекающиеся energy ranges.

External α/β semantics совпадают с dense path: β=0 не читает старое Σ,
complex α сохраняется, внешний qK применяется в `factor=α*qK`.

### Memory challenges закрыты по source

Actual `_is_momentum_independent(K)` проверяется **раньше** specialization,
поэтому rounded-flat radial case сохраняет уменьшенный Nb⁴ cache bound.
`2J<=Nb²` даёт
`Nk*chunk*(Nb²+2J)<=2Nk*chunk*Nb²` для трёх реально выделенных scratch arrays.
Reshape/view headers не создают ещё одного full scratch buffer.
Payload guard `2Nb²J+Nk²<=Nk²Nb⁴` сравнивает все три ComplexF64 payload arrays
с dense baseline, применимым после независимого momentum check.
`_kernel_operator_bytes` считает эти реальные payload arrays.

Estimator не уменьшается ни в resident, ни в peak и продолжает считать
старую raw K плюс generic cache/workspace upper bound. Это не экономия
всего process RSS и не изменение external cap. Плохие small-Nb/large-J случаи
переходят обратно к generic implementation.

### Test применимость и оставшиеся qualification шаги

Frozen suite включает independent literal/Kraus tensor с complex phased u,
signed d, unequal positive wk, все matrix units каждого input momentum,
complex scalar multiplication, tiny Choi positivity, PSD/adjoint/retarded
broadening witness и general complex GR. Это сильнее совпадения двух путей,
использующих одну flattening helper. Chunk1/3/16 при NE7, reuse, threads/BLAS,
α complex/β accumulation защищают packing и последний chunk.

Actual canonical production recipe сравнивается с независимым reference
finite-angle IFR builder; selected cache equality с literal проверяется
на random complex input. Actual source не меняется. Custom K mutation,
Nb1 memory fallback и tiny-Λ momentum-independent fallback предусмотрены.
Existing constructor/contraction/low-rank/candidate/memory regressions из
handoff применимы и должны выполняться на этом actual imported source.

Coverage caveats, сообщённые root/author без редактирования frozen source:

- manufactured complex-u test создаёт factors напрямую и не проходит
  canonical constructor. Actual canonical basis может быть real. Tiny
  basis-phase-rotation constructor case из design был бы полезной отдельной
  проверкой complex constructor path; source conjugations здесь проверены
  непосредственно, дефект не найден;
- custom source perturbation равна .01, а не one ulp. Code equality guard
  отвергает любой различающийся entry; one-ulp regression лучше защитила бы
  этот exactness contract от будущего tolerance-based weakening;
- произвольные normalization rescalings/finite-factor refusal не покрыты
  отдельными runtime assertions; проверенный code содержит q/finite guards.

Эти caveats не являются найденной source ошибкой, но должны оставаться
различимыми от выполненного coverage. Runtime пока `not_measured` рецензентом.
После unit qualification для speed claim нужен warmed same-field dense vs
structured measurement с одинаковыми weights/chunks/threads, absolute error
и input/operator scale, actual payload/allocations/time. Макроотношение MAC
и исходный IFR norm не превращать в end-to-end QCL speedup.

## Frozen test delta v2: canonical phase и one-ulp fallback

Независимо проверены incremental test artifact
`ifr-separable-test-v2.patch` SHA256
`35151a838304b0f9439622f64a80d9b73ec45c0b323f42e95022f6e80a36b70c`
и cumulative `ifr-separable.patch`
`8b47d71a2beedee65dfd7cb1620ebc8baa704870db3d9b1fb16d98d6d7cf2d7e`.
Четырёхфайловый implementation diff SHA256
`00e2e7fad176b9411794b8c756f39406b850281b870ec801d83d01658593f751`.
Дельта содержит только 48 дополнительных test lines; source logic не меняется.
`git diff --check` успешен. Прежний implementation verdict остаётся применим.

Phase fixture корректно задаёт χ/Φ→χ/Φ U и H0/Z/Minv/T±→U† A U при
diagonal unitary U. Interface vertex поэтому преобразуется D→U†DU,
а input/output fields сравниваются как U† X U. Это согласованная смена
representation, не изменение IFR physics. Actual normalized canonical
constructor теперь получает complex χ; structured selection подтверждается
меньшим payload, независимый reference finite-angle builder, literal map,
gauge covariance и positive-block witnesses проверяют именно этот path.
Ранее отмеченный complex-constructor coverage gap адресован.

One-ulp fixture восстанавливает исходный normalized K, заменяет реальный
diagonal covariance entry на его `nextfloat`, затем проверяет сохранение
этого entry, generic payload и literal-map equality. Исходный entry конечен
и real nonnegative; mutation действительно отличает source от canonical
recipe и защищает exact equality guard от будущей tolerance-based подмены.
Ранее отмеченный one-ulp coverage gap также адресован.

**Static test-delta verdict: ошибок не найдено.** После сообщения root о
408/408 existing six-suite qualification на неизменном implementation source
достаточно rerun изменённого IFR suite, а не всего набора. Его завершённый
log ещё не рассмотрен здесь; benchmark тоже pending. Code/PR acceptance и
speed claim этим test-only static update не выдаются. Рецензент runtime
не запускал.

## Independent runtime evidence update после test-v2

Прочитаны `.build/campaign133/ifr-native-checks-v2.log`,
`.build/campaign133/ifr-phase-checks.log`, actual `ifr-native-checks.jl` и
`environment/ifr/{Project,Manifest}.toml`. Рецензент Julia/Pkg/compute не запускал.

Первый завершённый native комплект содержит 408/408: исходный IFR suite
301+9 и unchanged broader construction/contraction/literal-low-rank/candidate/
memory suites 54+3+8+22+11. Script assert-ит actual imported Core path,
включает файлы через absolute paths и достигает final marker/exit0.
Wall19.48s, max RSS647496KiB, swap0. Предыдущая root harness include-path
ошибка не использована как zero-test acceptance и не скрыта этим результатом.

После test-only delta выполнен только изменённый IFR suite: **301+58=359/359**,
`IFR_PHASE_DELTA_CHECKS_COMPLETE`, exit0, wall8.59s, max RSS537924KiB, swap0.
Command содержит тот же realpath assert. С учётом неизменных 98 broader
assertions текущий qualified набор — **457 distinct assertions**, а не сумма
двух полностью независимых 408+359 runs. Повторять неизменный broader suite
без новых изменений/сомнений не нужно.

Обе команды явно задают `JULIA_NUM_THREADS=4`, `OPENBLAS_NUM_THREADS=1`,
BLAS.set_num_threads(1), `--startup-file=no`, `--check-bounds=yes`.
Таким образом requested2 worker cases действительно имеют доступный
multiworker pool. Это qualification parallel consumer, не scaling benchmark.

Manifest указывает `../../core-ifr-separable` и Julia1.13.3; modified Core
identity определяется frozen patch/source SHA из предыдущей секции, а не
одним base git HEAD. Проверенные SHA256:

```text
f2846c9d1f56d1a0c0c0fe21b2b5d9c1bfa552e6a23f68e4ab0ba55b61f7996f  environment/ifr/Project.toml
ad9a7e873c9ddd6b31688e9f74a6a39865a8b90f6f2fdc6e591fdfd633b45c13  environment/ifr/Manifest.toml
bbe2fabab2c07e7f5a0eebd74a56935a6a63ae0c3f70d70b0c88fef73606b875  ifr-native-checks-v2.log
ff77d929c4d6460cf8b6bddf0f65ba0b7b8084f9f7a03f785a112ab1115c5ce9  ifr-phase-checks.log
```

Implementation/operator equivalence и указанные guards теперь имеют runtime
evidence. Не измерены отдельно pathological normalization ranges, process
RSS guarantee, IFR share всего Γ map и whole-campaign benefit. Root запросил
same-field warmed benchmark до окончательной приёмки оптимизации; этот этап
ещё pending и не заменяется успешными tiny tests. Γ/scientific closure,
raw convergence criteria и accepted Core source всё ещё не менялись.

## Предварительный static review same-field benchmark harness

Проверены frozen write-only artifacts:

```text
3fff2d7ae9505e9ca35fdc86cf638dec66202c1f42ed1d842bdf9bd3ff319a26  ifr-benchmark.jl
7c398acbd127069b2d4f08dddf0334d89ea5a0256ae6af0911af708dc964adf7  ifr-benchmark-input.toml
a09db8ac90c1239a8f670e2783e79330a5dd746c9771c408ac8df7c6297dd3bd  ifr-benchmark-protocol.md
```

**Harness статически пригоден для заявленного bounded operator-only опыта**
при внешних Root caps4GiB/CPU400%/300s и resource-green admission. Этот verdict
не является результатом benchmark; actual performance пока не измерен в
данном review. Рецензент full HDF5/Julia/Pkg не читал/не запускал.

Script проверяет imported Core realpath, base e171, cumulative/source diff
hashes; затем streaming source-state/input hashes. Модель IFR строится на
NE17 только для energy-independent operator: семь basis operators,
normalized IFR K/qK, wk, x/wx, physical scales и full energy/k axes должны
побайтово совпасть с actual FD64 source. Несовпадение вызывает отказ без
basis alignment, interpolation или изменения raw source. Это позволяет
применить тот же operator к full GR NE30721, не объявляя малую модель
полным solver admission или новой grid qualification.

Обе arms используют один G, qK/wk/chunk64,4Julia/BLAS1 и preallocated full
outputs, α1/β0. Dense arm — прежний generic native dense kernel, structured
arm — фактически выбранный canonical cache. Warmup каждой arm исключён из
timing. Затем5 pairs с чередованием порядка; GC перед pair вне timing.
`@timed` честно включает apply scratch allocation/packing, но не HDF5/hash/
model/cache/warmup. Сохраняются все samples/GC/allocations, median/minimum;
это не statistical confidence interval.

Все complex output entries сравниваются с raw Frobenius/max numerators и
глобальными output scales; limit2e−12 — отдельная arithmetic qualification,
не ослабление scientific tolerance. G hash повторяется после опыта; outputs
полностью не сохраняются. Два returned full outputs остаются потребностью
этого опыта; explicit plan/cache payload не равен process RSS. Внутренний
phase-boundary check285s не заменяет external timeout внутри долгой операции.

Positive result позволил бы утверждать только measured IFR-apply benefit
для этого field/backend/chunk. Без отдельной IFR share measurement нет
оснований переносить его на Γ-map или whole SCBA speedup.

## Final owning documentation delta и publish freeze

После completed benchmark отдельно проверен docs-only patch
`ifr-separable-docs.patch` SHA
`74f69696d5995877789af3da04d36e3701b8c4821573f6d93b6cfaff39df45cf`:
20 новых строк `docs/src/theory/19_production.md`.
Ранг≤1 interface vertex, q/r/Σ formulas, complex adjoints, positive weights,
qK/wk once и generic/memory fallback описание согласованы с reviewed source.
Новых numerical/physical claims сверх exact discrete IFR apply нет.

Actual `git diff -- src test` всё ещё имеет tested/benchmarked SHA
`8b47d71a2beedee65dfd7cb1620ebc8baa704870db3d9b1fb16d98d6d7cf2d7e`;
4-source diff по-прежнему
`00e2e7fad176b9411794b8c756f39406b850281b870ec801d83d01658593f751`.
`git diff --check` успешен. Полный actual publish diff и artifact
`ifr-separable-publish.patch` имеют SHA
`14d3be2d95786d0ad5c5e458e99e24dd8ebcc8ee183a812efdc8f7e92dd72a8a`.

**Final owning Core PR/merge verdict: approve указанный publish freeze.**
Документация не меняет tested code, поэтому runtime retest или научный solve
не нужны. Benchmark freeze8b47 и позднейший docs-added publish14d3 сохраняются
различимыми; docs не объявляются частью выполненного более раннего запуска.
Итоговый result и ограничения см. `exact-ifr-benchmark-result.md`.
