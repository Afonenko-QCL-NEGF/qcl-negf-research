# Точная CP-сохраняющая факторизация упругих contractions

**Вывод:** в Core `e171caa60cce676c6c4238dd671e7c527f2bd977` есть
неиспользованная **при применении** точная IFR separability. На текущем
reference design четыре interface vertices имеют rank1, а angular factor
общий для всех interfaces. Это позволяет заменить dense275×275 contraction
тремя небольшими batched products без SVD, truncation или изменения physics.
Для AC/alloy основная momentum reduction уже реализована; их разложение по
пространственным rank1 vertices при Nb5/Nz96, вероятно, ухудшит стоимость.
Для impurity сопоставимого exact малого ранга из текущей алгебры не следует.

Статус: read-only proof/оценка, не implementation и не benchmark. Julia/Pkg,
сторонние расчёты и pilot не запускались. Frozen Core и harness не менялись.
`core-accepted-science` проверен clean на указанном HEAD; live Core clean f50.

## Фактический путь и что уже сделано

- `_interface_roughness_kernel_separable` уже собирает
  `K_IFR[m,m′,:,:,:,:]=angular_factor[m,m′]*block` для **того же finite φ grid**.
  Но затем возвращает полный six-axis array:
  [optimized/kernels.jl:452–494](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/kernels.jl#L452-L494).
- `_production_kernel` проверяет лишь полную momentum independence. IFR её
  обычно не имеет и получает обычный dense matrix; аналитические factors
  builder не доходят до contraction:
  [production_cache.jl:338–407](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/production_cache.jl#L338-L407).
- AC и alloy копируют один basis block на все `(m,m′)`:
  [reference/kernels.jl:240–310](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/kernels.jl#L240-L310).
  Cache обнаруживает это exact equality и выбирает
  `MomentumIndependentProductionKernel`: сначала `Gbar=Σm′ wm′ Gm′`,
  затем один Nb²×Nb² GEMM, затем broadcast по output m.
  [production_operators.jl:466–525](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/production_operators.jl#L466-L525).
- Dense path уже chunked, с одним `mul!` на chunk и private worker buffers:
  [production_operators.jl:401–463](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/production_operators.jl#L401-L463).
  Предлагать просто «добавить BLAS/batching» — повтор сделанного.

Первичная физическая опора: Lee & Wacker, PRB66,245314,
[DOI](https://doi.org/10.1103/PhysRevB.66.245314),
[полный preprint, Appendix A](https://arxiv.org/html/cond-mat/0212059),
рассматривает covariance vertices и SCBA для QCL scattering. Доказательство
ниже дано именно для **реального кода**, включая его quadratures и common
roughness parameters; публикация не заменяет эту привязку. Математическая
опора CP/Kraus: Choi, Linear Algebra Appl.10,285–290 (1975),
[первичная статья](https://doi.org/10.1016/0024-3795(75)90075-0).

## Короткое доказательство IFR

Для interface j код задаёт
`D_j[a,c]=(ΔEc_j/L0) conj(χ_j[a]) χ_j[c]`.
Положим `u_j=conj(χ_j)`, `d_j=ΔEc_j/L0` (вещественный, может быть отрицателен).
Тогда `D_j=d_j u_j u_j†` — Hermitian rank1. Current design имеет J=4:
[presets.jl:46–57](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/domain/presets.jl#L46-L57).
Interface interpolation, band jumps и seam convention нужно сохранить
буквально, а не переоценивать в другом basis:
[reference/kernels.jl:162–225](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/kernels.jl#L162-L225).

Пусть положительный scalar

`S_mm′ = πΔ²Λ²/K0 * (1/Nφ) Σφ exp[−(Λ/L0)² qbar(mm′,φ)²/4]`.

Это **текущая** finite-angle average, без замены Bessel continuum formula.
Для arbitrary complex matrix field X literal contraction равна

`Φ(X)_m = Σm′,j wm′ S_mm′ D_j X_m′ D_j†`
`         = Σj d_j² u_j [Σm′ S_mm′ wm′ (u_j† X_m′ u_j)] u_j†`.

Следовательно достаточно трёх стадий:

1. `t_jm′=u_j† X_m′ u_j` для всех j,m′,energy;
2. `r_jm=Σm′ S_mm′ wm′ t_jm′`;
3. `Y_m=Σj d_j² r_jm u_j u_j†`.

Это тождество компонентной суммы, не numerical low-rank approximation.
При column-major vec `(a,b)` первая/третья стадии — GEMM с
`R[j,(c,d)]=conj(u_j[c])*u_j[d]` и
`L[(a,b),j]=d_j²*u_j[a]*conj(u_j[b])`.
Conjugation существенна: `u_j†X u_j`, **не** `u_jᵀX u_j`.
Более простой запасной путь без rank1 специализации — exact Kronecker
`B ⊗ S` при native flattening `(m,a,b)` с fastest m; после basis-fast
packing это `S ⊗ B`. Здесь basis superoperator `B=Σj conj(D_j)⊗D_j`.

**CP:** каждый вклад есть `V X V†` с
`V=√(wm′ S_mm′) D_j`; все scalar weights неотрицательны. Поэтому карта из
block-diagonal input direct sum по k′ в output direct sum по k completely
positive. Не нужно предполагать PSD самой radial matrix S: достаточно
неотрицательности её entries. Для PSD Gn/Gp/A сохраняется PSD respective
Σin/Σout/Γ в exact arithmetic. Для arbitrary GR complex linearity и
adjoint compatibility дают `i[Φ(GR)−Φ(GR)†]=Φ(A)`.
Это **не** утверждение trace-preserving quantum channel или самостоятельный
сертификат conservation; equality сохраняет имеющиеся collision identities,
а их применимость всё ещё проверяется на исходной модели.

Generic truncated SVD в cache разлагает другой reshuffled dense map и этого
CP proof не имеет. Здесь factors получены из физических vertices без удаления
ни одного term. Знаки d_j исчезают только как d_j² после правильного
`D X D†`; нельзя заранее заменять знак другого physical coupling по аналогии.

## Нормировка и identity

Формулы выше используют **raw** microscopic kernel: K0 уже находится в S.
В existing API `Khat=Kraw/qK`, затем contraction умножается на `qK>0`.
Нужно выбрать одну эквивалентную convention factors и применить qK ровно
один раз, сохранив external α/β in-place semantics. Normalization source:
[optimized/kernels.jl:512–536](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/kernels.jl#L512-L536).
Нет дополнительного wE, spin factor, energy interpolation или KMS projection.
Exactness означает equality дискретной карты до изменённого floating-point
summation order, не bitwise identity. CP — математическое свойство factors;
малые floating negative eigenvalues требуют прежнего roundoff audit.

Нельзя автоматически применять physical factor recipe к arbitrary imported
или изменённому dense K без provenance/контракта. Для current builder common
Δ/Λ и independent interfaces доказаны source. Для correlated interfaces,
нового basis или interface-dependent spectra нужна новая проверка формулы.

## Затраты: arithmetic estimate, не измеренное ускорение

Обозначения N=Nk, B=Nb, J=interfaces; ниже complex multiply-add equivalents
на **одну энергию и одно поле**, без packing/scalar scaling.

| IFR path | Ведущая стоимость | При N11/B5/J4 |
|---|---:|---:|
| Existing dense `(NB²)²` | N²B⁴ | 75 625 |
| Exact `S⊗B`, два stages | N²B²+NB⁴ | 9 900 |
| Rank1 physical vertices, три stages | J(2NB²+N²) | 2 684 |

Отношения arithmetic counts ≈7.6 и≈28 к dense — **не wall-time speedup**.
Smaller GEMMs, packing и выходные записи могут стать bottleneck. В Γ pilot
число применений карты и доля IFR в total time должны измеряться отдельно.

Existing IFR cache matrix: `275²*16=1,210,000 bytes`.
Fused L/R/S storage: `2*25*4*16+11²*8=4,168 bytes`, плюс маленькая metadata.
Если сохраняется нынешний `KernelSet.K` contract, raw source tensor всё ещё
остаётся в problem; экономится только duplicate cache matrix (~1.2MB).
Нельзя обещать большой выигрыш total RSS: NE×Nk×Nb² state arrays не исчезают.
При chunk64 current dense X/Y scratch≈563200bytes/worker; два rank1 intermediate
`J×Nk×chunk` —≈90112bytes, плюс bounded packing input/output. Нельзя создавать
глобальный дополнительный `NE×Nk×J×components` scratch ради удобного reshape.

Julia план: упаковывать только energy chunks, объединять независимые поля в
columns, выполнять R→radial S→L через обычный `mul!`/reshape/views/private
worker buffers. Никаких новых packages или runtime generator не требуется.
Parallelism оставить по независимым energy chunks с BLAS1 либо явным прежним
`:blas` backend, без nested oversubscription. `mul!` требует неалиасных
input/output, что явно оговорено в
[официальном LinearAlgebra API](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.mul!);
[Julia performance guidance](https://docs.julialang.org/en/v1/manual/performance-tips/#Multithreading-and-linear-algebra)
объясняет риск Julia threads × BLAS threads. Конкретный crossover needs measurement.

## Почему AC/alloy/impurity не дают того же выигрыша

**AC/alloy:** каждый spatial sample тоже имеет rank1 vertex
`Q_g=conj(χ_g)χ_gᵀ`; веса `w_g` и `w_g*xAl(1−xAl)` положительны. Но physical
Kraus sum через Nz96 точек стоит около `2Nz B²=4800` basis operations на
energy после momentum sum, против existing dense reduced `B⁴=625`.
Структурная CP representation полезна для независимой проверки, но это не
хорошее ускорение при Nb5. Current momentum reduction уже снимает N².
AC и alloy не пропорциональны: alloy weight spatially dependent и равен нулю
в well. Можно лишь переиспользовать один `Gbar`/packing для двух различных
basis blocks и batched fields, сохранив отдельные channel outputs/diagnostics.
Это второстепенный выигрыш, не новая low-rank physics.

**Impurity:** optimized direct builder уже делает covariance BLAS:
`V(q)=prefactor*transpose(pair_density)*exp(−q*distance)`,
`C(q)=(V .* donor_weights)*V†`, с исключением только zero-weight donors.
[optimized/kernels.jl:144–186](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/kernels.jl#L144-L186).
V(q,D) — сумма spatial rank1 projectors с различными q-dependent weights и
обычно full-rank Nb×Nb. Она не общий fixed D times scalar radial spectrum.
Kraus apply по всем angle×donor vertices стоит порядка
`Nk²*Nφ*Ndonor*Nb³` вместо cached `Nk²*Nb⁴`; при множестве donor/angle samples
это обычно хуже. Exact semiseparable exponential recurrences могут удешевить
**construction**, но не cached per-iteration map; builder cost уже отдельно
оптимизирован и не должен подменять вопрос о Γ pilot runtime.

При q-lookup convex interpolation **covariances** сохраняет CP. Интерполяция
самих V с последующим outer product не эквивалентна: появляются cross terms
и меняются weights. Нельзя под названием ускорения менять этот numerical
contract. Без дополнительной exact algebraic identity нынешний cached dense
impurity contraction стоит сохранить.

## Малый различающий oracle и решение

Достаточная будущая проверка: NE7/Nk3/Nb3/J4, positive unequal wk, finite φ
как в source; complex interface vectors и signed d_j. Создать literal K
независимой суммой `Dac*conj(Dbd)` и сравнить: literal six-index map, existing
dense BLAS, Kronecker2stage, analytic rank1 3stage. Включить **все** matrix
units E_cd каждого input k′ и arbitrary complex GR; одних real/diagonal PSD
inputs недостаточно для обнаружения transpose/conjugation ошибки.

Отдельно Gn/Gp/A=XX†, phase rotation basis, zero interface height, negative
band jump и normalized qK scale. Проверить Hermiticity/PSD witnesses,
`i(ΣR−ΣR†)=Φ(A)`, α/β output accumulation, partial last chunk/reuse. Для CP
проверить также extended positive block input или Choi matrix tiny map;
operator equality плюс явный Kraus proof сильнее coincidence двух paths.
Нормировать error на `Σm′ wm′ maxm S_mm′ Σj||D_j||² maxm′||X_m′||F`, сохранив
absolute numerator; около нуля это input/operator scale, не output-only ratio.
Не clipping и не exact version/hash assertions.

После tiny equivalence — один разрешённый warmed microcheck той же карты на
S01-shaped frozen fields, с одинаковыми threads/chunks/weights, private scratch,
wall time/allocation/RSS. Предлагаемое окно≤30s/≤512MiB для tiny stage; большой
shape только отдельным ресурсным решением после текущего тяжёлого pilot.
Здесь ни один из этих tests/measurements не выполнен.

**Рекомендация:** ограниченный следующий implementation candidate — только
IFR structured operator, предпочтительно analytic rank1 factors с простым
Kronecker fallback. Preserve raw K baseline и current physics/normalization;
не вводить generic compression. AC/alloy оставить на existing reduced map,
потом при необходимости объединить packing. Impurity оставить dense до нового
различающего свидетельства. Даже подтверждённое ускорение карты не принимает
SCBA convergence, CP всей смешанной nonlinear iteration или реальную QCL model.
