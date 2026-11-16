# Core32: independent actual INFO-layout result

**Вывод:** benchmark выполнен по frozen протоколу; INFO padding не дал
наблюдаемого выигрыша. Baseline INFO W26 на этом fixture медленнее W4
в2,45944раз. Specific INFO-to-INFO false-sharing hypothesis данными не
поддержана: оба actual layouts имеют26 непересекающихся cache-line spans
изменяемого INFO word. Причина замедления W26 не установлена, end-to-end
ускорение native solver не измерено.

Primary `core32-info-padding-bench/attempt01/summary.toml` SHA
599026d2bdbcff8c9ed6d71b6720d6dde93b379edbdd1a57229f950dafd32c94.
Archived bench1c64732e… совпадает с frozen source; actual input отличается
только runtime authorization. Reproduction —
`reviews/core32_info_padding_result_audit.py`, результат
`reviews/core32-info-padding-result-audit.json` с pinned summary SHA.
Reviewer читал только малые TOML/source и вычислил статистику stdlib;
Julia/HDF/новых измерений не было.

**Protocol/equality PASS:** ровно14 calls; A26/A4 warmups исключены.
12 timed rows имеют строго заданный порядок W4:A/B/A×2, затем W26:A/B/A×2.
Ничего не отброшено как outlier. Все5 raw metrics и12 witness fields
действительно совпадают по Float64 bits между всеми14 rows; служебные
equality flags также true. INFO sentinel tails/прочие buffers addresses
и before/after input hashes неизменны. Witness сохраняет первый canonical
spectral E1/k1 при ratio0; это не exact-PD утверждение.

| Timed group | n | Median wall,s | Median process CPU,s |
|---|---:|---:|---:|
| W4 A:length1 |4|1,424912559|5,553470|
| W4 B:length8 |2|1,4309225645|5,626357|
| W26 A:length1 |4|3,5044928335|85,710753|
| W26 B:length8 |2|3,5712425295|88,663300|

B/A wall1,00422(W4)/1,01905(W26): point medians не показывают улучшения.
Малое число samples не даёт общего statistical/vendor-causality вывода.
Actual BlasInt8/L1line64; первые изменяемые8bytes каждого worker не
пересекаются с INFO другого worker ни в A, ни в B. Это не исключает
contention/false sharing других данных. Allocation medians367507120B(W4)
и367521856B(W26) практически одинаковы, источник allocation не измерен.

**Scope:** source e4fc9b91…/model3de1f7fc… — deliberate sum8 surrogate,
embedding double-counted, candidate===total. rD≈0,3432/rA≈0,4058/rK≈0,4446
не являются прошедшими physical gates. Измерен лишь public residual suite;
family rΣ/λ/freshΦ/SCBA и другие native phases сюда не входят. Поэтому
отношение2,45944 не есть full residual-phase или solver speedup.
Root execution exit0/69,54s, report RSS4323975168B, cgroup peak4166402048B,
high/max/oom events0; bounds соблюдены. Scientific/false-sharing flags false
корректны. Старый scientific contract не изменён.

**Практический следующий вариант — существующий phase API, без Core edit:**

```julia
phase_request = req -> req.name === :residuals ? min(req.task_width, 4) : req.task_width
```

На e171 `_production_phase_begin` принимает integer grant1..requested и
запрещает изменение fixed(nonelastic) phases (`production_telemetry.jl:130–149`).
Residual phase elastic=true; solver передаёт phase.task_width всему
`_production_fixed_point_measurement` (`production_solver.jl:577–603`).
Allocated workspace26 остаётся26; residual worker resolver допускает
активные4 (`residuals.jl:121`). Каждый chunk всё ещё покрыт, per-block
arithmetic/canonical reductions/witness ties и все raw gates сохраняются;
family norm reductions также canonical, но их стоимость здесь не измерена.
Existing `production_hotpaths.jl:177–202` уже проверяет grant через этот API.
Другие phases возвращают исходный requested width, включая fixed phases.

Такой driver-only execution вариант статически допустим: не менять total,
candidate, tolerances, algorithms, native guard либо global worker_count
на4. Нужна отдельная identity/actual evidence нового driver; этот benchmark
не принимает ещё не показанный driver и не обещает уменьшения полного
времени научного решения. Новых опытов данным review не запускалось.
