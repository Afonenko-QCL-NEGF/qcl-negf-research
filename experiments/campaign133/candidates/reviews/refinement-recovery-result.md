# Independent result review: uniform2x frozen collocation control

**Завершённые карты показывают существенную sensitivity на новых
midpoints, почти невидимую при old-node restriction. Это результат
одного frozen Γ refinement control, не диагноз stall или SCBA solve.**
Original two-map attempt и его postprocessing failure сохраняются.
Recovery действительно восстановил только diagnostics; newΠ=0.
Reviewer HDF5/Julia не использовал. Bounded stdlib small-artifact audit
128MiB: `refinement-recovery-audit.py/.json/.log`, PASS; исходная цена
audit0.05s/22584KiB. Повтор после добавления проверки admission/code hashes
тоже PASS. Никакие source arrays не открывались.

## Provenance/execution

Recovery runtime input1ff5a965567f807aeed56d723d9127093343ee5ee726c6893d585908ba38acd3
и main a1b2d1c4312cda8b4559a5b55f16314ada019252826d3eea9fada85d2885f9e7
совпали с copied artifacts/summary. Original summary3135bba16e5867d2d89cb3f8112d7c2191fdc767f5665d1aa64377e3949c1a4f
не изменён:2Π complete, known postprocess InexactError Int64(1.2).
Copied original code/input hashes independently совпали с original
summary. Actual root recovery tiny15/15 PASS,exit0,2.76s,436916KiB.
Actual recovery exit0,8.41s process/6.82155s script,618828KiB,Swap0.

Три HDF5 whole-file hashes root **precomputed**, записаны в admission
и runtime input; recovery сверил expected SHA, sentinel не использован:

| Artifact | SHA256 |
|---|---|
| coarse.h5 |27cefe8bed437b2afbfc4e0df4f2998565230516db11e20a17878ffc72d49a8a|
| fine-input.h5 |7917b12d55f49f066f70c21f91295832094b0258f96897bb71f172c11f18ef72|
| fine.h5 |87ad129508f6d4f468cec01abb2f2c892d9226f83bea2358073aac6548df8173|

Эти bindings независимо сверены между admission/input/recovery summary;
reviewer большие files не хешировал и underlying fields не пересчитывал.
Raw coarse/fine6 norms independently match original≤1e−10. Γf−PΓc
reported exact zero; residual identity weighted relative1.1996e−18.
Thus original error не означает незавершённую карту, но recovery
completion также не отменяет ни original failure, ни научные ограничения.

## Norms/partition result

D_f=Φf(PΓc)−PΦc, reference PΦc. Fine uses fresh native current μ,
coarse its own μ; D включает этот refit и остальные grid effects.

| Defect | Weighted absolute | Weighted reference | Relative |
|---|---:|---:|---:|
| Full fine |4.3201684586e−4 |0.16278055616 |2.65398310487e−3|
| All new even nodes |4.3201682653e−4 |0.11509207613 |3.75366264163e−3|
| Old-node restriction, coarse measure |1.8275863996e−7 |0.16279633581 |1.12262133577e−6|

Independent squared partition checks passed for numerator и reference,
weighted и unweighted:

    ‖Df‖² = ‖Dmid‖² + ‖Dc‖²
    ‖Df‖²_wfine = ‖Dmid‖²_wfine + ½‖Dc‖²_wcoarse.

Это использует точные fine odd weights=coarse/2. Midpoints несут
99.99998374% unweighted D² и99.99999105% weighted D². Именно numerator
partition устанавливает локализацию по old/new nodes; различие relative
norms само по себе этого не доказывало бы. Weighted D numerator≈445.44×
source coarse raw numerator norm, с явно разными fine/coarse measures.

Fine raw residual relative weighted2.65385873443e−3 против source
5.95754891326e−6; near agreement fine raw/D norms не означает причинной
тождественности. Old-node D заметно меньше, но и он не convergence
criterion. Source coarse replay exact Φ и HPI old-node difference≈1e−14
сохраняются в original summary; endpoint ghost taper variant h→h/2
остаётся отдельной маленькой model delta, без resolvent amplification bound.

Original fixed-source-μ Nfine=.045003199994375624, N/target−1=
7.11109861249e−5. Эти numbers recovery **не пересчитывал**. Current μf
уменьшился относительно μc; сам successful μ fit не доказывает solution.
Отделение его LO-вклада требует отдельного authorised μ diagnostic,
а не предположения по размеру Δμ.

## CSV support и строгий предел small audit

Selected CSV имеет5445 rows=495E×11k, indices/coordinates finite и
полностью совпадают с declared selection. Проверены nonnegative
squared norms и одинаковый inferred weight для raw/D columns в каждом
block. CSV SHA b949312834c423656faefdfd8378bce0389c61825ed8e0c238fdf04b249d1db6.

**Selected ladder покрывает только5.18917214% полного weighted D²**
(fine raw R² —5.20129476%). Global maximum по all-node streaming
e39394,k7,derived canonical E≈469.820607meV —even midpoint—отсутствует
в selected9785±LO ladder. Поэтому этот CSV не показывает главную
энергетическую локализацию нового defect и не позволяет объявить её
«тем же82.4meV feature». Full stream устанавливает maximum/value и
old/even split; ширина/вся support distribution около maximum не
измерены small CSV. Runtime anchor19697 в отдельном μ diagnostic
мог бы включить этот node, без приписывания нового результата исходной
selection задним числом.

Нет independent re-integration raw HDF5 здесь; small audit проверяет
внутреннюю согласованность primary metrics/selected data и frozen
code/provenance. Не установлены sufficient mesh/window/k resolution,
near-peak linewidth, Jacobian/nullmode, причина nonlinear stagnation,
улучшение будущего fine SCBA, physical/model validation или приёмка
native solver. Вывод ограничен measured frozen two-map collocation
sensitivity с explicitly changed endpoint hats и current-μ convention.
