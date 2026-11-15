# Независимый результат review frozen local-Γ nodal/PI control

## Подтверждённый вывод

На одном frozen local broadening несошедшегося FD64 изменение native
retarded reconstruction с nodal FFT на uniform PI существенно меняет
новый spectrum. При той же исходной trapezoid quadrature r_sum уменьшается
с.00492212 до.000406478, примерно в12.11 раза; fixed-μ N[fA] одновременно
увеличивается на.433806%. Γ, h, couplings и μ не подгонялись.

Свежая nodal arm практически воспроизводит archived A: relative weighted
L2 difference1.23e−6, тогда как PI−fresh nodal≈.00808. Поэтому в данном
опыте retarded-operator sensitivity хорошо отделена от archive/fresh-cavity
alignment. Это положительное различающее свидетельство, которого не давало
простое изменение квадратуры frozen source resolvent.

Нельзя принять SCBA, physical Keldysh state или заявить найденную физическую
ошибку baseline. PI задаёт другую finite-window analytic continuation;
уменьшение sum metric одного frozen model не доказывает её физическую
адекватность, nonlinear convergence или достаточность domain/grid.

Ghost-edge D45 и Γ scan дают сильные floating witnesses условий conditional
causality/normalization proof. Формальный numerical certificate, количественная
ghost-strip mass и отсутствие dark states для реальной модели без error
bounds этим не установлены. Runtime или arithmetic blocker результата
review не выявил.

## Primary evidence и независимая проверка

Сначала прочитаны `attempt01/summary.toml`, small HDF5 и archived code;
авторский диагноз для вывода не использован. Source protocol ранее прошёл
отдельный static review. API preflight `PIFrozenDefinitionAndBindingsOK`
выполнен до большого HDF5 I/O.

Пути отсчитываются от `.build/campaign133/`.

- Actual Core B `2a500740a82dfd5b44f6be966ef63dd54010d097`, actual path
  `core-uniform-pi-fft/src/QCLNEGF.jl`.
- Source FD64 SHA256
  `e6fffb997f36397f877aee1e2c98b1dbad7d97417d45458a042b5f074ce57a6d`;
  input SHA256
  `0fa0d0bc2fc27602fa99b8fe0797071eae99c66ab2661dfde05bd00ebe19b86e`.
- Archived `pi-frozen-control/attempt01/frozen_control.jl` SHA256
  `daa23694c74081359010595155d67ee588bd136831cdeed5be645f7d65c26344`
  повторно совпал с pre-reviewed code.
- `integrals_and_witnesses.h5` SHA256
  `ab34e09bb639ec7fb63f2d52df3f723ef40bed787e0f251194d0130a122af72d`;
  summary SHA256
  `b336c2f4dde87b55744e985db3dd8400eb410baf75c941d345bac63bb0f785fc`.

Независимая small-output проверка:

```
timeout 60s python .build/campaign133/reviews/pi_frozen_result_audit.py
```

Script использует только Python stdlib и системный h5dump, читает small
M/F matrices, per-k triangle bounds, h, selected response blocks и22 D45
witnesses; не читает полный source state или новые full arrays. Цена0.95 s,
RSS16,508 KiB. Точные команды и числа сохранены в
`reviews/pi-frozen-result-audit.json`. Julia/Pkg не запускались.

Независимо воспроизведены M/F/r_sum, N, global ΔM, ΔN/target, bounds и
matrix identities в двух selected blocks. Для D45 проверены Hermiticity,
trace/eigenvalue-sum и Frobenius²/eigenvalue²-sum consistency, а signs
дополнительно оценены Gershgorin bounds без повторного LAPACK eigensolve.

## Постановка и frozen identity

Источник — unresolved FD64: V0/U0, TL=TLO70 K, all five local scattering
families, finite_chain4 с каждой стороны, NE30721/Nk11/Nb5,
E∈[−.3,.9] eV, ΔE39.0625 μeV. В experiment не вычислялся ни один K[G].

LocalR source суммирует IFR+LO+acoustic+alloy+impurity ровно один раз,
исключая embedding и±. Frozen Γ=i(localR−localR†); сохранены hashes
localR и Γ. Оба reconstructors читают один Γ array без clipping,
renormalization и дальнейшего обновления. ΓK из sourceΣ</> здесь не
проверяется: scope намеренно retarded-only.

Семь numerical basis operators source v3 точно совпали с matching rebuild;
все reported Frobenius differences0. h SHA совпадает с прежними frozen
controls, Tminus−Tplus†=0. Energy/spatial/k nodes/weights и E0/L0/spin
проверены по source. Новые kernels не строились.

Для каждой fresh arm localΣR=Λ−iΓ/2, затем native cavity4 из этого localΣR
и central Dyson без η. Dummy zero Keldysh inputs нужны только native API;
физическая lesser/greater задача не решалась. В source/new outputs ничего
не ремонтировалось. No SCBA/mixing/Poisson и никакого μ-fit.

## Реконструкция и alignment

Nodal Λ и PI Λ имеют relative anti-Hermitian norms3.90e−16 и4.14e−16.
Reconstructed Γ против frozen Γ совпадает в weighted L2 примерно до
2e−16; occupied-relative discrepancy около7e−15, с сохранёнными small
numerators/scales. Нельзя смешивать этот маленький roundoff с quadrature
или window error.

Nodal rebuilt localΣR отличается от archived localR на relative weighted
L2≈2.42e−16, что дополнительно подтверждает правильную исходную Γ и старую
Hilbert functional. Сохранённая historical mixing арифметика для localR
не создаёт заметного operator discrepancy в этом snapshot.

PI−nodal Λ:

- weighted L2 absolute.000867795 при reference.0787383,
  relative.0110213;
- occupied-f weighted relative.000504537, где в squared norm стоит f².

После свежих cavity/central inversions:

- A weighted relative difference.00807742;
- fA occupied-relative difference.00770206;
- GR weighted relative difference.00776246.

Это измеренная sensitivity response, а не forward-error bound одной arm.
Из разных относительных норм Λ/A нельзя непосредственно вывести condition
number nonlinear map. Native denominator condition maxima1263/1274 не
свидетельствуют о машинной сингулярности и не принимают физическую точность.

В сохранённых max-A/max-fA blocks независимая complex matrix арифметика
восстановила A=i(GR−GR†) точно и Dyson backward residual от7.49e−17 до
1.14e−16. ΓRlocal−frozenΓ relative от1.42e−15 до3.27e−15. Это локальная
проверка units/index/order/subtraction на primary blocks, не проверка
полного спектра между energy nodes.

## Интегралы и стоимость charge изменения

Использованы source positive trapezoid weights:

```
M_k = Σ_e wE[e] A[e,k]/(2π),
F_k = Σ_e wE[e] f(E;fixed source μ,70K) A[e,k]/(2π),
N = gs Σ_k wk[k] Tr F_k,
r_sum = max_k ||M_k−I5||F/sqrt5.
```

Arrays, energy и k scaled; E0=.1 eV, L0=10 nm, target=.045.
Никакого дополнительного E0 в M и никакого wk в per-k r_sum нет.
μ=8.4953934398 meV фиксирована и служит только diagnostic FD weight.

| Arm | r_sum | N[fA]/target |
|---|---:|---:|
|archived source A|.004922121966|1.00000000000001|
|fresh nodal|.004922122187|1.000000047382|
|fresh PI|.000406477752|1.004338110723|

Числа из HDF5 воспроизведены независимо. Global unweighted
`||M_PI−M_nodal||F=.02541294545`; signed charge difference/target
.004338063341, то есть≈.00019521285 solver number. Это **не** raw λ
Keldysh equation: physical G</> и number constraint для fresh arms не
вычислялись. Refitting μ или normalizing N скрыло бы именно измеряемую
fixed-μ sensitivity и должно рассматриваться отдельно.

PI значение r_sum ниже numerical default.001; это проходит только
сравнение одного диагностического числа с прежним gate. Оно не принимает
весь state, другую energy-grid/window accuracy или source SCBA outcome.
Имеется trade-off: более близкая spectral integral к I при изменившемся
fixed-μ occupied number, а не универсальное улучшение всех физических величин.

Для каждого11 k независимо проверены discrete triangle bounds
`||ΔM_k||F≤ΣwE||ΔA||F/(2π)` и аналогичный f-weighted F bound.
Charge bound из trace inequality:

```
|ΔN|/target=.00433806334 ≤ .01016981121.
```

Это верные bounds differences выбранной discrete measure; они не
ограничивают continuum integration error или physical tail mass.

## Γ PD и full45 ghost-edge witnesses

Full nodal Γ scan сообщает minimum eigenvalue2.54639e−15 solver units и
minimum scaled eigenvalue.00205878; anti-Hermitian defect0.
Сохранённый worst-scaled block имеет eigenvalue6.38791e−6 и norm.00310276.
Global absolute minimum относится к другому возможному block и не должен
приписываться этим coordinates. PI/nodal A scans также положительны в
floating arithmetic; их min absolute eigenvalues около4.7e−17, с
min scaled eigenvalue≈2.61e−7. Строгого input/eigensolver error budget
для глобального PD certificate output не содержит.

Для ideal piecewise-linear PSD/PD Γ с adjacent zero ghost nodes Cauchy
model causal и нормирована на полной energy axis. Full9-cell D45(E) с
Hermitian h,T имеет D'(E)≽I снаружи ghost support. Поэтому при точных
предпосылках strict λmaxD(left)<0 и λminD(right)>0 исключают exterior
poles; strict PD Γ внутри support дополнительно исключала бы interior
dark atoms. Это conditional theorem, а не автоматический вывод из double
eigensolver output.

Все22 saved full45 matrices дают правильные floating signs:

- max по k λmax(left)=−3.16959117 solver energy units;
- min по k λmin(right)=+3.47813770.

Независимые Gershgorin bounds ещё проще подтверждают большой numerical
margin на **сохранённых** matrices: max left upper=−1.77562975,
min right lower=+2.08824410. Trace consistency≤5.69e−14,
Frobenius²/eigenvalue² consistency≤1.37e−12. Reported full eigen residuals
≤9.13e−16, anti-Hermitian defects≤1.01e−16. Наличие всех couplings в45×45
принципиально: знак5×5 isolated cell не был бы тем же witness.

Edge energies — ideal ghost coordinates около−.3000390625 и.9000390625 eV.
PI geometry сохраняет max node deviation1.78e−15 solver units и relative
deviation/h4.55e−12. Нет явных interval bounds coefficient/input/geometry
rounding или rigorous global PSD proof. Поэтому результат поддерживает
предпосылки теоремы практически, но не превращается в formal certificate.

Даже формально положительный exterior certificate относился бы к ghost
support[a−h,b+h], а spectral integral считается только на[a,b]. Масса
в двух ghost strips количественно не ограничена. Interior/threshold
singular mass и continuous quadrature error также требуют отдельных
условий/измерений. Нельзя интерпретировать r_sum.000406 как измеренную
missing ghost mass или заявлять отсутствие dark states реальной модели.

## Что теперь доказано и что остаётся

Опыт устанавливает чувствительность данного frozen retarded spectrum к
nodal→PI reconstruction при одной Γ, практически неизменном old-arm
alignment и неизменной интеграционной мере. Он показывает, что одна лишь
frozen-source quadrature sensitivity не описывала весь reconstruction
effect. Это основание исследовать выбранную retarded numerical closure
отдельно от longer SCBA relaxation и от direct-static functional C.

Не установлены: physical exterior continuation, общий causality defect
baseline, nonlinear fixed point PI, LO/static component dominance,
discretization/window convergence, physical Keldysh occupation, μ/Poisson
self-consistency и experiment. Отношение old/new r_sum не доказывает,
что baseline code был неправилен или что PI real QCL результат точнее.
Nodal samples могут допускать другую causal analytic interpolation;
proof compact linear Γ не опровергает её существование.

## Выполнение и ресурс

Primary runtime: script28.68 s, process33.89 s, RSS3,230,448 KiB,
swap0, exit0, Julia4/BLAS1. Это помещается в объявленный OS6 GiB/300 s
scope; output≈1.3 MiB. Full-scientific memory estimate≈12.94e9 bytes
проверялся отдельно против source budget и не выдавался за фактический
retarded-only RSS. API preflight и actual B identity прошли.

Таким образом execution, малая numerical consistency и frozen response
подтверждены; scientific acceptance новых SCBA states остаётся not_measured.
