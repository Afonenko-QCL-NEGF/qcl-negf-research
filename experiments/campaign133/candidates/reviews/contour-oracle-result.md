# Независимый result review frozen-Γ Matsubara oracle

**Подтверждён полезный observable-level результат:** off-axis плотность неизменной frozen full-hat PI-модели практически стабилизируется заметно точнее наблюдаемой разности с исходной real-axis quadrature. Это не новый SCBA result и не формальный certificate причины разности.

## Provenance и исполнение

Attempt `adjacent/contour-oracle/attempt01`, frozen script SHA256 `5249e38a95ebed2b4e487af586251bdda9b4ba0d6060fb5c1bdb96fc10880449`, input `296877b8bf496d048ecf91089121900d6d34e618170f635c63c5560de649ea0d`, summary `8cbf7010df8cb559fca5bf59fa9023a68281318341e2cb532bead6c7d0ea0b5b`, small density H5 `3c56f7d34463d2431cf255721fa9a8c6446ea9f83e97173eef18f7a5eb83c101`. Script совпадает независимо просмотренному static freeze. Actual identity совпадает source FD64 `e6fffb…`, frozen local Γ `edb853…`, hcell `8cde8e…`, prior PI reference `ab34e0…`, μ=.008495393439818052eV, T70K, E0=.1eV, L0=10nm, spin2, Julia1.13.3/1thread/BLAS1. Canonical/source geometry defect1.77636e−15 solver units сохранён. Embedding в Γ не включён; full9-cell/45-state Hamiltonian и fixed μ не изменялись.

Root log: exit0,33.30s process wall,29.47s script, peak981796KiB; все N128/256/512/1024 закрыты. Отдельный no-source preflight complex2×2 даёт error8.12066e−16. Main scalar/complex45/source-H45 Γ=0 controls выполнены; source-k1 M3 error на1024=1.81436e−12 Frobenius, complex45=9.96228e−15, complex density imaginary offdiagonal signal=.06394. Ни fitting, ни SCBA/Poisson вызовов нет. `SCBA_convergence=not_measured` корректен.

## Независимый пересчёт малых артефактов

Сохранены exact воспроизводимый `reviews/contour_result_audit.py` и `contour-result-audit.json` с командами h5dump/dataset paths. Запуск: `timeout 60s python .build/campaign133/reviews/contour_result_audit.py`; только stdlib и small H5,6.44s, peak22020KiB Python (HDF tools отдельными краткими subprocess), без Julia/Core/full source H5. Проверены:

- ВсеρM1/ρM3/raw checkpoints, joint/per-k Frobenius differences и sheet number `2Σk wk Trρ` из HDF, не из summary. Согласие summary до floating accumulation roundoff.
- Tail reconstruction из raw partial density+a1/a3 сохранённых moments — нулевая разность при том же arithmetic expression. Независимый Decimal80 пересчёт odd-frequency T2/T4 совпадает с сохранёнными Float64.
- M1/M2/M3/a1/a3 пересобраны из hcell,Tplus,S0,S1 через явные5×5 matrix products. Максимальное absolute различие3.42496e−14; cross terms h*S0+S0*h и neighbor paths присутствуют.
- Для n0/k1,n127/k1,n511/k6,n1023/k11 сохранённый full45 central Gc независимо воспроизведён двумя4-step5×5 Schur recursions с stdlib pivoted inverse. Relative discrepancies1.06e−16–4.25e−16. Это проверяет central compression/indexing для selected frequencies, не analytic Σ на всех frequencies.
- Ghost β из endpoint norms/margins/h пересчитаны; максимумβL+βR=9.438349685e−20, никаких отрицательных margins или |d| подстановок нет.

## Что говорят density checkpoints

N здесь — число положительных Matsubara frequencies; `number/target` — dimensionless sheet density при **неизменной source μ**, не raw native Keldysh λ.

| N | M1 number/target | M3 number/target | M3 change/target от предыдущего |
|---:|---:|---:|---:|
|128|.992858267237|1.004611344553|—|
|256|1.003097149976|1.004566318265|−4.50262882e−5|
|512|1.004381257690|1.004564904777|−1.41348748e−6|
|1024|1.004541904637|1.004564860555|−4.42219531e−8|

M3 changes при doubling уменьшаются примерно в32 раза, как ожидается для оставшегося O(N−5) tail; M1−M3 matrix correction уменьшается в8 раз. Это сильное numerical asymptotic свидетельство, **не rigorous remaining-error bound**. На1024 разность M1−M3 ещё2.29559e−5 target; её нельзя путать с измеренной M3 refinement error или объявлять error estimate финального M3.

Из reference matrices независимо получено: fresh_PI number/target=1.004338110722785; fresh_nodal=1.000000047381868; archived_source_A=1.000000000000008. На1024:

- M3−fresh_PI: signedΔnumber/target=+2.26749832656e−4 (то есть +.02267498% target), matrix absolute norm4.48358738776e−5, relative norm1.57866755048e−4 относительно ||fresh_PI/F||F=.2840108664.
- M3−fresh_nodal: signedΔnumber/target=+.00456481317357; M3−archive=+.00456486055543. Это разные comparisons: PI vs nodal меняет retarded closure, archive comparison дополнительно несёт исходное history/alignment.
- M3(1024)−M3(512): matrix norm9.70856e−10, signedΔnumber/target−4.42220e−8. Наблюдаемая discrepancy с real PI больше последнего refinement change примерно5128 раз по sheet number и46182 раза по matrix norm. Одна лишь недоделанная M3 Matsubara sum на наблюдаемом asymptotic уровне не выглядит основным объяснением discrepancy.

Matrix discrepancy сосредоточена преимущественно в первых occupied k: per-k norms3.13870e−5,2.08611e−5,2.39014e−5, затем4.30537e−6; последние k дают≈1e−11. Global unweighted matrix norm и weighted sheet number имеют разные меры. Этиρ — equilibrium occupation matrices для frozen model, не stored GL native SCBA trajectory.

## Numerical accuracy и attribution

Минимум |t|≥485.13. BigFloat exact kernel controls показывают analytic first-omitted bound≈8.98749e−22 для наиболее близкого representative t; Float64 coefficient error≈7.12e−21 отдельно и не выдаётся за analytic remainder. На сохранённом n0/k1 full Σ omitted-series Frobenius bound4.25721e−19, затем быстро убывает. Все selected solve backward errors≤4.51104e−17, accumulated Herm-part defect0; selected positive ImΣ eigenvalue violation0. Это полезная согласованность, но не full forward-error bound дляρ: Float64 summation/Γ input uncertainty и propagation в resolvent не сертифицированы. Γ=0 controls не являются oracle energy-dependent Σ; независимыми частями являются hat analytic/BigFloat controls, moments и selected Schur reconstruction.

Matsubara формула включает **всю spectral measure full-hat model**, включая ghost strips и потенциальные atoms. Real fresh_PI/F интегрирует только [a,b] по исходной trapezoid сетке. Поэтому discrepancy нельзя безусловно подписать «чистая quadrature error». Float ghost witnesses дают dL_full≥3.16920,dR_full≥3.47775 solver units и βmax≈9.44e−20: **при exact PSD/input/positive-margin assumptions** ожидаемый continuous ghost contribution на много порядков меньше наблюдаемой density разности. Это не interval certificate этих assumptions. Отсутствие interior dark atoms требует отдельного PD/no-singular-mass основания; floating positive Γ eigenvalues прежнего control сами по себе его не сертифицируют.

Практический вывод сильнее исходного real-vs-real comparison: независимая off-axis integration поддерживает наличие малого, но разрешимого density integration discrepancy в frozen PI-модели. Условно на корректный causal/no-atom input и малую evaluation error она согласуется с real-axis energy quadrature error; bounded evidence не доказывает эту атрибуцию формально. Сам source FD64 остаётся unconverged, Γ/H frozen: это ничего не устанавливает о nonlinear tail SCBA, причинности/PSD будущего fixed point, достаточности basis/k/window/embedding, Poisson или experiment.

## Решение

Принять narrow frozen-observable witness и его стоимость, сохранив условные qualifiers. Не нужен новый большой N только ради уменьшения уже существенно меньшей Matsubara refinement delta: для строгого attribution полезнее propagation/evaluation и PD/no-atom bounds, а для решения исходных уравнений — отдельный unconstrained native audit. Ни spectrum normalization, ни новая μ, ни ретроактивный SCBA pass из этих результатов не следуют.
