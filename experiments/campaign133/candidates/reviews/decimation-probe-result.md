# Independent result review: frozen-Γ PI boundary limit

## Вывод

Для данного frozen Γ, исходной energy/k quadrature и фиксированной μ
занятый интеграл N[fA] слабо меняется уже при finite4→8 и практически
стабилен при16→32→∞. **Полный retarded spectrum и его spectral integral
такой сходимости не показывают:** A заметно меняется, а r_sum немонотонен.
Малый r_sum finite4 не выбирает физически правильную boundary length.

Принято узкое свидетельство boundary sensitivity одного frozen retarded
model. Runtime/primary consistency blocker review не выявил. Не приняты:
SCBA/Poisson, infinite interacting fixed point, достаточность energy/k/basis
window, выбор infinite model как физического baseline или объяснение
прежнего raw Keldysh residual.

## Primary identity и исполнение

Пути ниже относительно `.build/campaign133/`. Прочитаны только code,
input/summary/log и small output HDF5; full source state не открывался.

- Run `adjacent/decimation-probe/attempt01`, accepted clean Core
  `e171caa60cce676c6c4238dd671e7c527f2bd977`.
- Actual archived main SHA
  `b713adc8a5682c953a0654b9fc5d26c5de9a262f924c75c8d4ee719703e8e209`;
  algorithms SHA
  `65e6af6ae0e23c0c80e874d5993ec70eb418a4b201932c10241ae225b0b4bb8b`;
  input SHA
  `60cf66758a9ffb4a8565fdca862a24e3f2529c03e341ce0a768c7bc90c42c850`.
- Source FD64 e6fffb997f36397f877aee1e2c98b1dbad7d97417d45458a042b5f074ce57a6d;
  frozen local Γ edb8530e9a59944afe1a447d999bf0c6c42d29d726e16144a5d85bdccdaecc55;
  h 8cde8ea024d9e690bf0396428574496734d7bad5962fc00617b9b380ae1194a7.
- Small HDF5 SHA: boundary integrals
  `6fee16a44e1bf54ef986caf9c489b6b5c01cdb9a95aec1f2a0553223a3d35339`;
  alignment
  `db415f6e51f38eb199364d2db05c1834d936cce4d5d4393728bab694a40b1062`;
  selected oracles
  `e90755c8a1f050eff58c0079a9945308f74a1fdf762dda753d3e888425b7b2f5`;
  summary
  `df7ce9397f8647726aa79e5f9044b027f61932874fcb5894bced433d56475d09`.

Actual main отличается от reviewed v2 только завершающим explicit if-guard.
Мой прежний static review не выявил `@__FILE__ && main()` parse hazard.
Root include failure сохранён в unexecuted-v2parse; guard correction
сохранена отдельным patch и затем root definitions-only include прошёл.
V1/v2 algorithm tiny success не был использован вместо этого include.
Алгоритм/источник/thresholds при исправлении guard не менялись.

Полный diagnostic завершился exit0:49.19s process/41.99s script,
peak2307968KiB, swap0, Julia4/BLAS1, все11 k. Это укладывается в
объявленный6GiB/600s retarded-only scope, не измеряет стоимость полного
solver. Archive/log и прежние unexecuted/parse outcomes сохранены.

## Независимый small-output audit

Скрипт `reviews/decimation_result_audit.py` и первичные пересчёты
`reviews/decimation-result-audit.json`; финальный stdout/time —
`reviews/decimation-result-audit-final.log`. Команда:

```sh
/usr/bin/time -v timeout --signal=TERM --kill-after=5s 90s python .build/campaign133/reviews/decimation_result_audit.py
```

Только stdlib и h5dump, RLIMIT_AS256MiB; audit проверяет, что каждый
допущенный HDF5 output меньше2MiB. Final run1.86s wall, RSS25464KiB,
exit0. Пересчитаны raw per-k M/F, r_sum, sheet N, все семь sensitivity
pairs и77 per-k triangle inequalities. Сохранённые numerator/reference
scales использованы независимо от summary ratios. M/F Hermitian defects
в сохранённой арифметике равны0. Нет Julia/Pkg/new source arrays.

Для selected k1/6/11 повторно собраны45×45 D из сохранённых D0/T,
выполнен независимый stdlib pivoted Gauss–Jordan inverse, а обе Bloch
суммы2048/4096 независимо вычислены простыми5×5 inverses с reuse even
nodes — всего12288 tiny inverses. Это selected witness audit, не новый
energy-grid или SCBA experiment. Все команды HDF5 dataset reads
сохранены в JSON; input/source frame guards выполнял основной harness,
а полный sourceSHA этим reviewer заново не вычислялся.

## Что фиксировалось и какой смысл имеют интегралы

Источник — несошедшийся FD64, zero bias/U=0,70K, Nb5/Nk11/NE30721,
E=[−.3,.9]eV, ΔE39.0625μeV. Sum localR содержит IFR,LO,acoustic,alloy,
impurity ровно по одному разу, исключает embedding/embedding±. Frozen
Γ=i(localR−localR†), затем accepted uniform full-hat PI reconstruction
ΣlocR=ΛPI−iΓ/2. Ни Γ, ни h/T, ни μ не подгонялись.

Все клетки одинаковы и используют этот один local ΣlocR(E,k). Finite
длины означают4/8/16/32 клетки с каждой стороны центральной; infinite
arm использует две independently decimated surfaces. Это разные boundary
variants при одной local functional, а не решение по отдельной Σ каждой
клетки. Lesser/greater state, scattering map и continuity equation тут
не решались. Dummy zero Keldysh buffers в native alignment нужны только
retarded API и не объявлены физическим состоянием.

M_k=Σ_e wE A/(2π), F_k=Σ_e wE f(E;μ,T)A/(2π),
N=gsΣ_k wk TrF_k, r_sum=max_k||M_k−I5||F/√5.
Arrays/weights dimensionless, E0=.1eV, spin2, target=.045.
μ=8.495393439818052meV фиксирована по source marker, нового fit нет.
N[fA] — thermal diagnostic frozen spectrum, **не** raw native Keldysh λ
и не нормированная charge acceptance.

## Actual alignment и independent algebra

В alignment.h5 native finite4 M/F **численно точно** совпали
с сохранёнными prior PI matrices. Независимое чтение actual prior HDF5
подтвердило, что comparison использует именно его fresh_PI, не копию с
другой identity. Native r_sum/N differences0; independent recursion4
против native M/F даёт global Frobenius5.91e−15/1.19e−16.
Independent4−native A weighted relative7.4251e−15.

Selected witnesses:

| k | E, eV | Independent full45 centre vs saved | Independent Bloch4096 vs decimation |
|---|---:|---:|---:|
|1|.0153125|2.20e−14|1.36e−14|
|6|.066171875|5.34e−15|2.25e−15|
|11|−.0191796875|3.33e−16|1.25e−15|

Самостоятельные Bloch matrices совпали с saved4096 до7.26e−16 relative;
2048→4096 independent refinement≤3.16e−15. Direct surface/centre backward
residuals≤8.17e−17, full45 backward≤5.76e−17. Проверены обе ориентации
TgRT† и T†gLT, без принудительного βnew=αnew†.

Это три выбранных occupied peaks, не все worst-conditioned/high-energy
узлы. Global recorded scan прошёл337931 E/k nodes: все11 histograms
суммируются в30721 каждый, max20 iterations<64; stop/equation/causality,
finite-centre и nonfinite failure counts0. Maximum correction/threshold
≤1; max saved surface backward4.66e−14<1e−9. Max surface condition
≈2.24e5: маленькая backward error сама по себе не равна forward bound.
Полный global eigenscan reviewer не повторял; его finite/causality
counters проверены против frozen code/summary, selected equations
пересчитаны независимо.

Combined status корректно требует all completed k, zero failures,
finite integrals/sensitivities и все selected finite/full45/stop/Bloch
refinement/agreement flags. Проверенные данные согласованы с
`surface_finite_and_selected_Bloch_witnesses_pass`; первая непомеренная
θ-refinement NaN имеет refinement_measured=false и не выдана за измерение.
`completed_execution_only` и SCBA/Poisson/discretization/experiment
not_measured остаются отдельными полями.

## Boundary sensitivity: occupied number и весь spectrum различаются

| Boundary | r_sum | N[fA]/target |
|---|---:|---:|
|finite4|.000406477751593|1.004338110722785|
|finite8|.013464469275656|1.004338351765906|
|finite16|.012061196865754|1.004338353403627|
|finite32|.006425472394966|1.004338353406339|
|infinite decimation|.004060506971131|1.004338353406221|

Числа пересчитаны из raw M/F. Occupied sheet number changes/target:
4→8=2.4104312e−7,8→16=1.6377207e−9,16→32=2.71295e−12,
32→∞≈−1.18e−13;4→∞=2.4268344e−7. Последние два difference находятся
на намного меньшем масштабе, чем4→8; не превращаем~1e−13 в certified
absolute error bound continuum observable.

Однако weighted A L2 difference4→∞=0.294749 относительно finite4,
а32→∞=0.174092 относительно finite32. В occupied-fA L2 соответственно
1.4507663e−4 и6.83511e−11; measure wE·wk, occupied squared norm содержит
f². Signed number чувствителен ещё и к trace/energy/k cancellations,
поэтому он не заменяет матричную occupied норму. Global ΔF4→∞=
1.72988584e−7, ΔM4→∞=.0094219270; ΔM32→∞=.0218396806.

Каждый discrete triangle bound проверен независимо. Например
|ΔN4→∞|/target=2.42683e−7 ≤1.3675402e−4;
|ΔN32→∞|/target≈1.18e−13 ≤3.49013e−11. Это bounds выбранных discrete
differences, не доказательство continuum quadrature/tail error.

Отсюда допустимо сообщить: **occupied integral слабо чувствителен к этой
boundary extension при данном frozen Γ и выбранной quadrature; полный
sampled spectrum и матрица его интеграла M заметно чувствительны.** Недопустимо:
«finite32 сошёлся по всему GR/A» или «малый r_sum finite4 доказал правильное
embedding». R_sum меняет худшийk (1→2→11→11→2) и немонотонен.

Изолировать оставшийся spectral-sum deficit этим опытом нельзя: boundary
variant меняет pole/band structure, а fixed energy sampling/window
оставлены прежними. Возможные contributions underresolved weakly damped
features, finite-window mass, ghosts и numerical integration не разделены.
Наблюдения не доказывают конкретное cancellation/tail объяснение и не
разрешают выбрать вариант по наименьшему r_sum задним числом.

## Сопоставление с сохранённым contour oracle

Cross-identity checks подтвердили одинаковые sourceSHA/ΓSHA/h/μ/E0 и
prior PI reference. Contour M3(1024) относится **к finite4 full-hat model**,
а не к infinite arm. Его stored finite4 reference совпадает с нынешним
independent finite4 F до1.19e−16 Frobenius.

На этой matched finite4 паре Ncontour/target=1.004564860555441,
Nreal/target=1.004338110722785, разность2.26749833e−4 target,
matrix difference4.48358739e−5. Измеренная discrete boundary4→∞ number
sensitivity≈934 раза меньше этой contour−real разности. Перенос real arm
на∞ оставляет cross-boundary difference≈2.26507149e−4 target и
matrix norm4.48474074e−5, но это уже **разные boundary models** и не
новый matched quadrature oracle.

В частности, данный опыт не показывает, что прежний finite4 contour−real
gap вызван finite-chain cutoff: обе прежние arms уже имели одну boundary.
Он также не даёт Matsubara density infinite chain и не наследует finite4
ghost/atom bounds автоматически для∞. Formal attribution оставшегося
gap требует прежних causal/input/no-atom/evaluation conditions; actual
continuum/window и SCBA adequacy по-прежнему не установлены.

## Решение для дальнейших утверждений

Публиковать можно qualified frozen-boundary observation, runtime и
primary consistency evidence, сохранив отдельные occupation/spectrum
результаты и source unresolved status. Для исходного научного fixed point
нужен независимый native unconstrained audit; увеличение chain length
само по себе не принимает source SCBA. Для точности spectral sum нужна
отдельная различающая проверка energy quadrature/window/analytic measure,
а не новый выбор boundary по красивому r_sum.
