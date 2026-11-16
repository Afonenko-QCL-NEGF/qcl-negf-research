# Independent result review: bounded fine2x/history1 stage

**Bounded execution и снижение fine raw residual подтверждены; convergence
не достигнута.** Все14 attempts/enteredΠ/completedΦ/successful согласованы:
1 source coarse replay+12 fine+1 terminal coarse restriction, без failures,
fallback или native exit. Этот результат позволяет рассматривать новый
ограниченный Γ-only continuation после source/resource freeze; сам review
новый расчёт не запускает и не принимает SCBA.

Reviewer использовал только малые TOML/log/code/file sizes, без Julia,
maps или HDF5 payload. Audit128MiB `fine2x-m1-result-audit.py/.json/.log`
PASS,0.03s/17848KiB. Все14 map records проверены отдельно; HDF5 arrays
и большие checkpoint SHA reviewer не пересчитывал.

## Frozen evidence и accounting

| Artifact | SHA256 |
|---|---|
| final summary |9a9a1d5c4f7e0176e68910255d283b072fed039401b22e8593e253f54565c18a|
| actual input |844b092393b4c300d29259d550cee522eb76f8b09682ada1df9b1f438d0c8061|
| actual/copied V3 code |8f73ba2feb07d1a2ea15900745e85e338a897fc1ad81cf6b0597c6c258e611c5|

Source7 hashes соответствуют stableA-cont64 unprojected Γ-last64;
малые source TOML hashes подтверждены, model binding согласован с
pinned frame. Coarse cache exact equality=true и whole source coarse
rawΦ replay exact0. Qualifier c48/stream8df/cache7d unchanged. Actual
runtime input отличается от frozen template только authorization.

`fine_last` exact равен map_fine_012, `coarse_terminal` exact равен своему
map TOML. На диске14 map и14 response files, есть seed/initial/terminal
и terminal-coarse checkpoints, нет .tmp/first-domain-failure files.
Reviewed code обновляет G/Φ только после successful evaluator и пишет
last aligned Γ/rawΦ в finally. Это подтверждает program/metadata alignment;
independent HDF payload equality не заявляется. Failed V1/V2 parse controls
до tests/Φ остаются отдельно от V3 tiny и actual.

11 update events:1 linear,10 AAII history1/rank1, affine l1≤6.79945960
при разрешённых4096, condition1≤1e8. Affine sum/L1/coefficient metadata
проверены независимо. Нет PSD clipping/projector; AA uses Herm(rawΦ)−Γ,
criterion использует rawΦ−Γ и reference rawΦ. Max representation
conversion relative1.65e−16; raw/working norms согласованы с conversion
triangle bound и не подменяют друг друга. Все q negative_count0,
number-fit guards/finite records/domain margins выполнены. Это floating
backward-policy witnesses, **не новый exact-PD certificate terminal Γ**.

## Измеренный progress

| Record | Raw relative | Weighted raw relative |
|---|---:|---:|
| Fine001 |2.62528169225e−3 |2.66090136822e−3|
| Fine012 |3.22179074311e−4 |2.87298508163e−4|
| Fresh coarse terminal raw residual |2.44350551119e−3 |2.45763997353e−3|
| Terminal cross-grid Φ difference |2.52003908000e−3 |2.51608127582e−3|

Fine12 residuals монотонно снизились в этой конкретной trajectory:
8.14851709× unweighted,9.26180016× weighted. Threshold остаётся1e−9
по обеим нормам×3 maps; каждый measured fine residual выше него,
strict_reduced_streak=0. SCBA_accepted=false/native exit отсутствует.
History3→1 и coarse→fine — отдельная execution/discretization identity;
это не isolated acceleration comparison с прежним source pilot.

Cross-grid numerator есть ‖Φc(RΓf_terminal)−RΦf(Γf_terminal)‖ в **coarse
measure**, denominator freshΦc. Это свежий operator-response diagnostic
на terminal state, не fine raw fixed-point residual и не continuum error.
Current μ fit собственный у каждой map; boundary hats h→h/2 различаются.
Cross-grid weighted relative≈8.76× fine weighted raw relative не является
физическим error ratio между решениями или доказательством superiority
coarse/fine. Он показывает, что уменьшение fine residual не сделало
terminal response grid-independent. Достаточность h/2/window/Nk и
причина исходной stagnation этим не установлены.

Fine grid61407, nominalLO1878, own native trapezoidal weights/fingerprint
сохранены; old nodes nesting проверено. Saved ghost bound1.9712e−28
получен **из initial source endpoints**, uniform-in-E combined-k Λ norm.
Это не terminal endpoint bound и не resolvent/global-density certificate.

## Ресурсы/continuation boundary

Root log exit0:213.21s process/206.73650s script,6581748KiB≈6.28GiB,
Swap0. Data1,772,810,930B<4GiB; каждый file<2GiB. OS/numerical cap
13019119616B=12.125GiB,High11.5GiB,J4/BLAS1,CPU400%,900s enforced;
31F+3GiB admission=11,597,140,272B<cap, actual full-native estimate
25,860,597,424B честно refused before allocations. Нет budget inflation
или разрешения native full workload из успеха stream map.

Полезный continuation source — `attempt01/fine-terminal.h5`
(Γ/rawΦ/energy/wE;571424042B) вместе с map_fine_012/summary/input и
pinned original frame. Перед новым actual root должен freeze whole
checkpoint SHA и проверить payload/grid/source binding; здесь SHA
большого HDF5 **not_measured_by_reviewer**. Saved terminal — достаточный
Γ-only warm seed, но последние2 AA states/F не сохранены как restart
history: empty-AA continuation — новая identity, не bit-exact trajectory
continuation. Old sourceΓ64 exactPD не переносится автоматически на
новый terminal Γ или будущий Mats source. Нужны свои domain/causal
qualifications по вопросу следующего опыта; этот review не добавляет
нового oracle/runtime или native acceptance.
