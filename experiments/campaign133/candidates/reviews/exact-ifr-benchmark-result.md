# Независимый результат IFR same-field benchmark

**Принят узкий результат:** exact canonical structured IFR apply быстрее
прежнего dense apply примерно **2.20× по медиане** на одном сохранённом
actual FD64 complex GR при прежней карте, весах и normalization. Вместе с
static review и 457 qualified assertions это достаточно для принятия
проверенной Core optimization в отдельный PR/merge. Итоговый committed source
должен совпадать с freeze ниже. Научная модель/критерии не меняются.

Ни actual IFR share всего Γ-map, ни end-to-end SCBA speedup этим опытом
не измерены. Большие HDF5/arrays рецензент не открывал, Julia/Pkg/benchmark
сам не запускал; прочитаны tiny archived source/input/summary и ограниченные
phase/time строки completed log, пересчитана только малая scalar metadata.

## Source и scope действительно совпадают

Base Core `e171caa60cce676c6c4238dd671e7c527f2bd977`, isolated modified
`core-ifr-separable`; cumulative patch
`8b47d71a2beedee65dfd7cb1620ebc8baa704870db3d9b1fb16d98d6d7cf2d7e`,
4-source implementation diff
`00e2e7fad176b9411794b8c756f39406b850281b870ec801d83d01658593f751`.
Оба independently recomputed текущих git diffs совпадают с benchmark gate.
Static source/test qualification — `reviews/exact-ifr-static.md`.

Actual archived benchmark script SHA
`3fff2d7ae9505e9ca35fdc86cf638dec66202c1f42ed1d842bdf9bd3ff319a26`
и input SHA
`7c398acbd127069b2d4f08dddf0334d89ea5a0256ae6af0911af708dc964adf7`
побайтово совпадают с предварительно reviewed freeze. Script содержит
executable realpath/base/patch/implementation gates; completed outcome
подтверждает прохождение этих gates, не только печать ожидаемого имени.

Source HDF5 SHA проверен самим frozen harness streaming pass:
`e6fffb997f36397f877aee1e2c98b1dbad7d97417d45458a042b5f074ce57a6d`.
Source dataset — actual `green/GR`, field shape **30721×11×5×5**, payload
**135172400 B**. Это complex retarded field сохранённого несошедшегося FD64
state, не synthetic random field и не новое accepted решение.

Model NE17 служит только construction energy-independent IFR operator.
Перед apply script требует bitwise equality full source axes/measures/scales,
семи basis operators, normalized IFR K/qK и wk с canonical construction.
При несовпадении нет basis alignment или geometry replacement. Нормировка
qK=.0019883546138492096 и wk применена одинаково обоими arms.
Этот small model не объявляется полной solver grid или обходом full solver
RAM admission. Проверить large source bytes независимо данным read-only
review не пытались; оценён реальный consumer, assertions и completed outcome.

Сохранённые Project/Manifest SHA совпали с actual environment/ifr files:
`f2846c9d1f56d1a0c0c0fe21b2b5d9c1bfa552e6a23f68e4ab0ba55b61f7996f`
и `ad9a7e873c9ddd6b31688e9f74a6a39865a8b90f6f2fdc6e591fdfd633b45c13`.
Julia1.13.3, Julia threads4, BLAS1, chunk64. Source G hash before/after
совпадает:
`cfd144796094dacb28a42cbae5a834db173b1fc3e62c9923d9d841305bcfe50b`.

## Timing и allocation пересчитаны независимо

Каждая arm warmed один раз; warm .4961/.3637 s исключены. Пять measured
pairs действительно чередуют dense→structured / structured→dense.
Preallocated outputs одинаковой полной формы; in-place apply включает
scratch allocations и packing. GC перед pair вне timer; timed GC0 во всех
10 samples. Хэширование/чтение HDF5/model/cache construction не входят в
operator timing. Это сравнение consumer operations, не construction benchmark.

| Arm | n | Медиана, s | min…max, s | Allocations на apply, B |
| --- | ---: | ---: | ---: | ---: |
| Native generic dense | 5 | .080033348 | .079872273….081585710 | 2261856 |
| Canonical structured | 5 | .036443100 | .035215603….037111629 | 1490144 |

Независимый stdlib пересчёт median ratio:
`.080033348/.0364431 = 2.1961180031336522`. Все пять pair ratios лежат
в **2.1607…2.2681**; более быстрый исход не зависит от единственного порядка
или выброса. Это n5 локальное измерение, не confidence interval и не
гарантия будущего speedup на другом backend/chunk/field.

Actual cache payload dense **1210000 B**, structured **5136 B**, что совпало
с `(11*25)^2*16` и `(2*25*4+11²)*16`. Последнее включает ComplexF64 radial
matrix; более ранняя design оценка с real radial имела бы другой размер.
Payload savings не являются измеренным total RSS savings. Вход/оба outputs,
полный raw K и full SCBA state в scientific run от этого не исчезают.

Process завершён exit0: wall **21.66 s**, user22.45 s/system1.70 s,
max RSS **1049868 KiB**, swap0. Harness elapsed17.3904 s включает отдельный
streaming source hash9.1762 s. Эти полные времена не подменяют operator medians.
Проспективные external4GiB/CPU400%/300s — controls ведущего, а не свойство
алгоритма; рассмотренный log показывает нормальное завершение внутри окна.

## Arithmetic equivalence и immutable input

Сравнены все complex entries последнего полного output обоих deterministic
arms. Tiny unit/reuse oracles отдельно защищают preceding apply semantics.

```text
absolute Frobenius difference  = 4.1501837181654613e-16
reference Frobenius scale     = 0.7511503947801508
relative Frobenius            = 5.525103557164675e-16
maximum absolute entry error  = 8.265575096131627e-18
maximum reference entry scale = 0.005936808968139032
maximum scaled entry error    = 1.3922588953915046e-15
worst linear index            = 8895
```

Оба quotient независимо пересчитаны; scales не около нуля, hidden floor не
объясняет малые discrepancies. Заданный arithmetic limit2e−12 пройден без
изменения criterion. Output hashes закономерно различаются при иной
floating summation order; это не bitwise-output claim. Stored input G hash
проверен после всех repetitions и не меняется.

## Итог и границы

Implementation + exact-source admission + canonical complex-phase/one-ulp
fallback + independent matrix-unit/CP tests и данный real-field comparison
дают достаточное свидетельство roundoff-equivalent optimization выбранной
discrete IFR map и полезного operator timing на проверенном режиме.
Не требуется повторять unchanged broader suites или запускать полный
SCBA только ради этой optimisation acceptance.

Не принимаются замена physics, новое CP свойство всей nonlinear iteration,
Γ/native convergence, новый warm trajectory, grid/boundary adequacy,
IFR share других channels, end-to-end speedup или автоматическое обновление
активного научного source identity. Сохраняются private raw arrays и
предыдущие baseline outcomes; публикационные projections требуют отдельного
отбора путей.

Артефакты этого result:

```text
1765c19d847a18bfae366bad74a2c1709b786e89faa2dbbeef36abf6474bf9e4  ifr-benchmark-attempt01/summary.toml
c4da927a9c3d9fcce56338752f4b463fa07942ead3d1d90cfd1c174713a39bc6  ifr-benchmark-attempt01.log
```
