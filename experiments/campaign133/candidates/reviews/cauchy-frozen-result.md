**Независимый аудит actual frozen Cauchy transform, 2026-11-16**

Решение: **узкий результат transform-comparison принимается**. Все 120
сохранённых 5×5 сравнений независимо пересчитаны из малых real/imag TOML
матриц и проходят объявленный diagnostic gate. Это результат FD64 frozen
Γ на k1/6/11, 16 source nodes и четырёх off-grid targets при двух τ.
Green functions, integrated density, SCBA, nonuniform-source convergence
и adaptive economy этим опытом не измерялись и не принимаются.

Первичные результаты: `adjacent/cauchy-tree-prototype/frozen-attempt01/summary.toml`,
SHA256 `ad2fab31e9a1e4f96b2cc1200545ab2a18be987aaf43619aefc3ced3c3515672`;
`frozen-attempt01.log`, SHA `fae065aa146738a8696ab215006af5aaaa03c844e22d6f5e55c7021276aa7ceb`.
Frozen script `2e5777af…7765`, input `546a2f57…a63a`, module `a74ba641…257d`,
tests `3be8f0d8…d62f` совпадают с сохранёнными копиями и текущими frozen
файлами; accepted imported Core — e171caa60cce676c6c4238dd671e7c527f2bd977.
SHA всех малых входов и покадровые пересчёты находятся в
`reviews/cauchy-frozen-result-audit.json`; воспроизводимый stdlib-only скрипт
`reviews/cauchy_frozen_result_audit.py`. Он ограничивает address space
128MiB; actual `/proc/self/status` HWM 19648KiB, VmPeak33240KiB.
Γ/HDF5, Julia/Pkg и Core reviewer не открывал/не запускал/не менял.

**Численный результат**

Нормы — Frobenius в исходных безразмерных единицах Σ, E0=.1eV. Для
физической величины Σ нужно умножить их на E0. Они не являются ошибкой
плотности или SCBA residual. В каждом τ имеется 48 source-node сравнений
с native accepted uniform PI FFT и 12 off-grid сравнений с independent
Big192 scalar positive-hat divided-difference reference.

| τfar | Reference | Максимум ||Σtree−Σref||F | Где |
|---|---|---:|---|
| 1e−9 | FFT nodes | 4.6831105e−11 | k11, source index4097 |
| 1e−9 | Independent off-grid | 2.3614255e−11 | k6, μ/E0+0.2i |
| 1e−11 | FFT nodes | 3.2126309e−13 | k11, source index4097 |
| 1e−11 | Independent off-grid | 1.9137409e−13 | k11, μ/E0+0.2i |

Все 60 matched target pairs имеют **тот же reference, z, k и frozen Γ**
между двумя τ. При tighter τ сохранённая разность меньше во всех60:
error_tight/error_loose от0.002047 до0.016678. Это наблюдение выбранных
точек; universal monotonicity, convergence order или bound всего поля
не следуют из него. Максимальная разность outputs двух τ —4.6511889e−11.
Статистика отдельных записей сохранена в audit JSON, а не восстановлена
из итогового scalar pass.

Независимо пересчитаны ||value−ref||F, ||ref||F и gate
`far_estimate+max(1e−11,1e−10||ref||F)` для всех120. Stored norms/allowances
совпали с stdlib `fsum` recomputation в relative2e−13/absolute1e−30
арифметическом допуске; этот допуск относится только к повторению norm,
не заменяет исходный physical/numerical gate. Максимальная доля gate,
занятая actual difference,0.059034 (loose FFT); tight FFT0.017741,
tight off-grid0.011083. Кодовые diagnostic_pass=true подтверждены числами.

**Ledger и reference witnesses**

Для каждой записи finite/nonnegative numerical ledger fields проверены;
`near+far+direct+zero=covered=30722` пересчитано. Это проверка сохранённых
counts; structural partition argument взят из отдельного fixed-source
math review, не заменён суммой counters. Во всех120 far estimate≤τ,
максимум estimate/τ=0.869162. Max far radius/distance≤0.454977<.5,
orders≤24, accepted nodes≤visited nodes. Числа прямых/far panels при
tighter τ показывают ожидаемое уточнение, но не дают самостоятельного
performance/economy вывода.

Reverse-reference witness у off-grid не больше1.697744e−15; у FFT nodes
поле0 — **placeholder** (native FFT doubled-padding witness сохранён
отдельно per-k, около3.5e−16), а не выполненная reverse FFT проверка.
Максимумы ledger sensitivities: moment reverse2.640e−17, polynomial
evaluation2.793e−17, near sum reverse2.013e−19, final sum reverse2.350e−17.
Это order-sensitivity witnesses, не доказанные total forward errors.
`formal_whole_certificate=false`, geometry/log bounds not_certified,
far arithmetic not_outward_rounded проверены и сохранены без переименования.

**Геометрия, источник и единицы**

Все три k имеют один energy digest85b2bfe2…c354. Γ slice digests записаны
отдельно: k1 `385e638a…76c2`, k6 `5d5059f9…24df`, k11 `f2688d0f…7052`.
Аудит сохранил/сверил metadata identities и frozen extraction contract;
**не пересчитывал** source-file/aggregate Γ/slice hashes по большим массивам.
Их actual вычисление — guard/выход выполненного frozen script; completed
status и checkpoint path согласованы со static-reviewed control flow.
Aggregate Γ SHAedb853…ecc55 остаётся prior reference и не выдан за
recomputed selected-k aggregate. Full source e6ff…a6d исходно не менялся.

Tree использовал literal stored knots/adjacent ghosts, native FFT —
canonical uniform hats. Per-k actual ghosts `[-3.0003906249999996,9.000390625]`,
ideal FFT ghosts `[-3.0003906249999996,9.000390625000001]`; max ghost delta
1.7763568394002505e−15. Canonical max node deviation той же величины,
relative к spacing0.000390625 —4.547473508864641e−12. Эти равенства
пересчитаны из малых metadata. Literal field не объявляется тождественным
ideal canonical interpolant в exact arithmetic. Off-grid reference
относится именно к literal hats и обеим нулевым ghost nodes.

Raw input witness сообщает Hermitian defect0, min eigenvalue positive
(на k11 минимум2.546385796281463e−15) и explicit numeric allowance;
representation_changed=false. Witness inequality проверено по записанным
числам; без Γ reviewer не может заново выполнить eigenvalue/scale audit.
Это не formal PSD certificate, clipping/repair не выполнялись по frozen code.

**Исполнение и ресурс**

Сохранены `completed_diagnostic_only`, последний stage complete_k11,
exit0/swap0; script19.286999s, whole process26.84s, peak1026420KiB
(1051054080bytes); script Sys.maxrss1050882048bytes близок, но не обязан
совпасть с поздним whole-process `/usr/bin/time` maximum.
J1/BLAS1. Tiny source-free preflight root:67/67 плюс typed PI aliases,
exit0,7.20s,RSS564820KiB. Reviewer эти процессы не повторял.

Admission JSON явно изменил proposed1536MiB на actual MemoryMax2GiB
ради source import/tree/reference JIT без изменения numerical controls;
фактический launch scope root записал high1.75GiB/Swap0/CPU100%/120s.
Причина и actual cap не скрыты за исходным input proposal. This run
поместился в выделенное окно; сравнительная экономия против FFT,
contour/density implementation или full adaptive SCBA не измерялась.

Следующее действие определяется основным SCBA вопросом и ресурсным окном
ведущего. Повтор этого same-source probe для получения ещё одного pass
не требуется; новые утверждения про off-axis density/Γ evolution/adaptive
grid потребуют собственных фиксированных inputs и различающих проверок.
