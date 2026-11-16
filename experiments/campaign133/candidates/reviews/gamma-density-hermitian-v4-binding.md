# Независимый V4 filled-binding review

**APPROVE frozen source binding** `input-stable64-frozen.toml` SHA256 `ab77e333726eec4cf258a74ea2572b3aa2336ae642c01b55181f412660b2db7b`. Блокеров не найдено. `execution_authorized=false` сохранён; actual Matsubara runtime/result **not_measured**, текущая heavy работа не получает дополнительной нагрузки от review. Reviewer использовал только малые text/TOML/log и stdlib Python≤128MiB, без Julia/Pkg/HDF5/source arrays.

Исходный template SHA `916a191a11c704c2a9329985f1f4054545c6170a0cbe8107b7b7d086ae589ff4` побайтно равен `frozen-io01/input-pre-digest.toml`. Oracle814dff… и helper192e6f… совпали с ранее принятым static freeze; copied IO helper также побайтно совпал. Численная/source-gate реализация не менялась; `gamma-density-hermitian-v4-static.md` сохраняет область и ограничения.

Completed digest SHA `b121f71d58cbf909cd4087e6ddae4b0216883a375529e5a179a562d2b381ef40` прочитан напрямую, status=completed; original template/helper identities совпали. Root IO log: exit0,8.97s process/7.80987s internal,RSS702332KiB,swap0,J1/BLAS1,timeout120s. Digest сообщает exact checkpoint↔terminal Gamma equality всех8443600 complex entries; это результат reviewed full-slab сравнения, не extrapolation selected nodes.

Независимые whole-dictionary comparisons: ровно13 array hashes и11 A-slice hashes input=digest, keys A точно1…11, все24 значения64hex. Origin−2.9999999999999996 и spacing0.00039084132055378063 перенесены exact. Raw-h hash8cde8ea0… равен digest/metadata; выбранный **actual saved h_H** hash9748cff3… отдельный от raw h. Metadata shapes/types согласованы: Gamma30704×11×5×5, hcell/rawh11×5×5, Phi/chi96×5, остальные basis operators5×5; grids Float64, остальные ComplexF64. hcell metadata не хранит повторный sha leaf, но авторитетный `array_sha256.hcell` есть и совпадает input/digest.

Все15 digest source/environment hashes согласованы с input. Core e171/tree742, μ=0.08459614262982018, file paths, units/T70, target, estimator/checkpoints, caps и исходные identities не изменены. Малые source input/identity/map/summary/ledger подтверждают unprojected stable64, all9 policy, accepted map64/summary ownership и clean Core. Preserved runtime проверит actual whole-file SHA/HDF/logical arrays повторно при actual.

Parsed template→filled delta содержит **ровно35 leaf changes**, независимо проверено точным whitelist:

- 13 array hashes+11 A hashes+2 geometry+raw-h hash;
- 4 поля certificate/proof path+hash;
- status и3 новых provenance leaves freeze_digest/freeze_digest_sha256/source_template_sha256.

Ни source file hash, ни μ/units/grid/estimator/cap/authorization не изменены; UNSET leaves больше нет. Нет скрытого fit/projection/backend изменения.

Original ALL certificate a703ec4c… и proof-input27da5a97… пересчитаны; cert input hash=proof hash. Все9 source triplets cert/proof/input совпали, включая Gamma-last **cab810…**, terminal6eb3… связан с ним completed IO exact equality. Strict fields complete/all/PD337744/20264640 exact divrem/zero failures проверены независимо; script d239… и math c6b… совпадают. Дополнительная независимая scientific qualification certificate сохранена в `global-gamma-exact-pd-result.md`; здесь не повторялся HDF5 proof scan.

После завершения текущей heavy работы и root resource-green может быть создан **отдельный** actual run-input с единственной parsed delta `execution_authorized: false→true`, сохраняя этот false frozen input. Если изменение сделано ровно заменой единственной строки `execution_authorized = false` на `execution_authorized = true` без прочих byte edits, ожидаемый SHA256 run-input — `f6560732368b3ba1a4ccae52004a9a23eac63c18775ff6d614fa6ed6afbeee2b`. Reviewer такой input не создавал и actual не запускал; при ином SHA проверить точную delta, а не переносить approval на любые изменения.

Это готовность binding к bounded actual, не научная приёмка ещё несуществующего Mats result. Exact no-UHP-poles относится к declared frozen Hermitian/Cauchy-hat model; finite-series/Float64/LU/tails, possible atoms и real-axis comparison требуют собственных witnesses. SCBA/native convergence, grid adequacy, эксперимент и старый Gamma128 HOLD не приняты.
