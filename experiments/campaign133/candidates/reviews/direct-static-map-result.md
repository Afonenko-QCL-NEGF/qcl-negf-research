# Независимый review actual C static-map control, attempt02

## Verdict

Узкий результат принят: на одном фактическом FD64 Green input native
`direct_green` меняет четыре static retarded self-energies примерно на
.78–.97% их старой unweighted Frobenius нормы. Lesser/greater candidates
пяти каналов и весь LO candidate остались побайтово неизменными.
Полная discrete difference identity с фактическим input defect замыкается
на масштабе arithmetic roundoff; разность retarded практически Hermitian.

Это actual confirmation отдельной closure functional и её отличия от
finite-window nodal Hilbert baseline. Не приняты: SCBA, новый retarded
spectrum/embedding, физическая tail attribution, причины прежней несходимости
или превосходство C для устройства. Одного same-input local-map сравнения
для этих утверждений недостаточно.

v1 attempt01 остаётся timed-out partial outcome. Его три завершённых
channel artifacts воспроизводятся v2 точно по всем исходным полям,
но отсутствие final summary и source-final-hash check не позволяет
переименовать v1 в completed. V2 — отдельный успешно завершённый attempt
с механической дельтой harness и иным resource envelope.

## Frozen identity и что независимо проверено

Пути отсчитываются от `.build/campaign133/`.

- Core `2dc5a26c1a7a6260febe3ab07bf6439fd2f0bf19`, tree
  `30747daa72e86c7f59442c600d456c2a725164a1`.
- Source FD64 state SHA256
  `e6fffb997f36397f877aee1e2c98b1dbad7d97417d45458a042b5f074ce57a6d`;
  source input SHA256
  `0fa0d0bc2fc27602fa99b8fe0797071eae99c66ab2661dfde05bd00ebe19b86e`.
- Archived attempt02 harness SHA256
  `98a17f98e12b50fe6678e73285f1513204724fbe09c8a54a4a0a61ce075982a1`.
- `witnesses.h5` SHA256
  `947c92de69af350710e931e0d0f27aa028ae9cc9ee6b99c02f9d73000c691d73`;
  summary SHA256
  `d57658ff4a70f41318d36dba549313c36d201cd034980bb4f98475e34afbb643`.

Прочитаны exact v2 diff, input/identity/summary, четыре small channel TOML,
selector hashes, runtime logs, API/type preflight и small HDF5 witnesses.
Julia/Pkg и большие arrays не читались/не исполнялись. Полный source SHA
проверялся самим harness, а text/small artifact hashes сверены рецензентом.

Независимый tiny audit:

```
timeout 30s python .build/campaign133/reviews/direct_static_result_audit.py
```

Он читает только8 selected5×5 blocks через системный h5dump и Python
stdlib; цена1.03 s, RSS15,620 KiB. Точные команды и результаты:
`reviews/direct-static-result-audit.json`. Все source/code artifacts
оставлены неизменными.

## Что сравнивается

Оба native calls получают один и тот же GR/G</>/A object из source
unconverged FD64. Нет FD repair, μ-fit, density normalization, mixing,
solver, Poisson, cavity refresh или central Dyson. Source Green hashes
проверяются до/после и final assertion прошёл.

All5 actual kernels строятся один раз с теми же physical/grid parameters
и используются обоими calls. Различается только elastic_retarded selector;
cache validation подтверждает общие kernels и FFT plan. Это корректное
reuse, поскольку selector не меняет representation самого K/H.

Algorithm manifest обоих arms имеет coarse impact=physical_model из-за
уже выбранного finite_chain по приоритету общей taxonomy. Это не означает,
что selector C является новым physical-model change: его собственная дельта
остаётся controlled numerical closure, явно сохранённая в manifest.

## Независимый вывод и проверка signs

Обозначим `B=−iG<+iG>`, `D=A−B`, а H — один и тот же finite nodal
principal-value оператор с source wE. Для energy-independent complex-linear
adjoint-preserving static K:

```
ΣR_new = K[GR] = K[HermGR] − i K[A]/2,
ΣR_old = H K[B] − i K[B]/2,
ΔR = K[HermGR−H A] + H K[D] − i K[D]/2.
```

Переход использует commute H и K, поскольку H действует только по E,
K по k/basis. Это не предположение о continuum Hilbert relation source GR.
Также `ΓR_new−ΓK_new=K[D]`.

V2 сохраняет эти signs: residual равен
`delta−firstterm−Hdefect+.5i*Kdefect`. K/GR/adjoint indices и меры
соответствуют ранее reviewed native C. Greater не восстанавливается из
lesser и input D не зануляется. Сам input defect имеет full weighted
norm4.37e−14, что согласуется с исходной почти машинной pair spectral identity.

В8 HDF5 witness blocks независимо получено:

- `new_R−old_R−saved_delta` в точности0;
- пересчитанная difference identity residual в точности равна saved residual;
- energy index→physical energy совпадает до≈8.4e−17 eV;
- anti-Hermitian delta мала и согласуется с сохранёнными global norm scales.

Worst-delta examples: impurity E≈50.1953 meV, k5,
block delta norm.00111840 solver units; IFR E≈228.711 meV, k6,
block norm.000542180. Это local self-energy block norms, не shifts
определённого физического level или accuracy gain.

L/G blocks не сохранены в small HDF5, поэтому ΓK identity нельзя заново
пересчитать только из этих witness blocks. Для неё review использует
full-field scalar measurements native run и inspected code; не выдаёт
small-block reconstruction за независимый повтор всего массива.

## Primary численные результаты

| Static channel | ||ΔR||F / ||Rold||F | weighted relative | identity residual / ||ΔR||F |
|---|---:|---:|---:|
|acoustic|.00960054|.00960054|3.45e−14|
|impurity|.00966257|.00997911|5.54e−14|
|IFR|.00781353|.00811690|7.71e−14|
|alloy|.00830181|.00830182|4.55e−14|

Независимо восстановлены ratios из numerator/denominator TOML, а standalone
channel TOML точно совпали с соответствующими summary subtrees.
Full/trusted measures — `gs/(2π) wE wk` в squared component norms;
trusted interval[−.22,.82] eV. Relative weighted differences близки к
full unweighted здесь, но это разные definitions и они не смешиваются.

Нормы anti-Hermitian Δ относительно Δ равны примерно5.7–8.6e−14.
Main term K[HermGR−HA] воспроизводит практически весь наблюдаемый Δ;
K[D] и H K[D] имеют tiny arithmetic scale. ΓR_new−ΓK_new−K[D] остаётся
малой абсолютной величиной около1e−16–1e−15 unweighted. При практически
нулевом input defect относительное отношение двух tiny остаточных норм
не является физическим критерием; absolute values сохранены.

Сравнение hash tables independently подтвердило:

- Σ</> старого/нового candidate совпадают для всех пяти channels;
- LO ΣR/Σ</> совпадают полностью;
- source arrays unchanged assertion достигнут до final summary.

Dynamic e-e отсутствует: его invariance в этом actual run не измерялась,
только сохранён source exclusion из code review.

## Независимость и научный предел

Discrete decomposition использует production K/H по обе стороны, поэтому
она — consistency actual component, не независимая физическая валидация
самих K/H. Дополнительная независимость приходит из уже выполненных
CP/rational/scalar oracles C, а не из одного совпадения identity residual.
Вместе они поддерживают корректность выбранной route в рассмотренном scope.

Hermitian разность показывает отсутствие значимого нового broadening при
данном согласованном input. Но `HermGR−H_W A` включает и finite-window/PV
discretization, и возможную ошибку global spectral relation самого source.
Нет измеренных exterior masses или independently qualified full-axis GR.
Поэтому всю разность нельзя объявить physical tail correction или
количественным доказательством старой KK ошибки.

Около1% localΣR difference может иметь как малый, так и усиленный эффект
после fresh cavity/Dyson/nonlinear iteration. Эти responses здесь не
вычислены. Нельзя переносить результат PI frozen control на C: там
изменялись все local reconstruction arms и свежий spectrum, здесь только
четыре static functional на неизменённом G. Current, gain, rK/λ full coupled
map, spectral sum и accepted solution C не измерены.

## V2 mechanical qualification и отрицательный v1

V2 вводит concrete Array{ComplexF64,4} barriers, typed per-channel arguments,
hoists dictionary lookups за element loops и возвращает временные arrays
из ограниченного function scope. Формулы, loop/index order внутри audit,
H operator и selectors сохраняются. Telemetry добавлена отдельно.

Программно проверено: все13 original fields acoustic/impurity/IFR partial
TOMLs v1 совпадают с v2 **точно**, включая nested norm records. Добавлены
лишь audit_seconds/allocated_bytes/gc_seconds. Это конкретная qualification
механической дельты для всех завершённых v1 channel outputs; alloy v1 не
завершился и отдельного v1 numerical oracle для него нет.

V1 был остановлен timeout/SIGINT, final summary отсутствует. Строка GNU
time «Exit status:0» после «terminated by signal2» не является completion:
решение wrapper124 и отсутствие terminal artifact имеют приоритет для
execution classification. Три partial channels сохраняют собственную
область, не принимают полный attempt. Его resource trajectory тоже
отрицательный технический результат, а не physics failure.

V2 completed exit0: script≈34.15 s (summary записан на34.09 s), process46.26 s,
RSS6,831,000 KiB, swap0, outputs≈140 KiB. Per-channel audit≈.78–1.10 s,
allocated≈1.22e9 bytes; allocated bytes не равны peak resident memory.
Actual root envelope8 GiB/high7 GiB/CPU400%, Julia4/BLAS1/heap6G,
timeout600 s отличается от v1. В input-v2 остаётся prospective CPU100%
metadata; ресурс определяет actual external launcher, а не это поле.

Нельзя по этим двум attempts изолированно приписать всё ускорение одной
причине: одновременно изменены type/inference barriers, memory envelope
и compilation state. Обоснованный вывод — v2 завершился в новом bounded
окне с qualified numerical agreement по трём прежним partial channels.
Время до accepted scientific C result этим опытом не измерено.

## Решение

Сохранить completed v2 как actual same-input evidence C и v1 как partial
timed-out baseline. Следующий научный вопрос должен явно выбирать fresh
retarded response либо отдельно identified nonlinear C solve; единственная
неизменная L/G карта на старом GR не заменяет такое свидетельство.
Никакого послабления acceptance критериев и никакого принятия source FD64
задним числом из этого результата не следует.
