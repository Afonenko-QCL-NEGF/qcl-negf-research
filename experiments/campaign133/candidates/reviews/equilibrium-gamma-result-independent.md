# Независимый результат Γ pilot16, attempt02

**Ограниченный опыт выполнен; raw невязка FD-restricted Γ-map уменьшилась примерно в95 раз, но сходимость не достигнута.** На16-й Π rΓ=.0203879940, weighted=.0186891421; все16 triggerfalse. Native unconstrained validation не выполнялась. Это перспективная конечная trajectory и честный cap outcome, не доказательство существования fixed point, mathematical contraction или native SCBA решения.

## Источники и выполнение

Primary directory `equilibrium-gamma-solver/attempt02`: archived main SHA256 `bcd20186d8be8b3aa7e9f1b6c76677eb85a4c4af33a80933e70193a23434b37e`, functions `470a8e477129c5f23d23482057a638fa7c057922f5d5b0a54d618553b85b5d0d`, checks `cc0e143bd0c0af8fb295a7992a28313c342c3b115d047cc773d159ce47ab111e`; actual authorized input `c0846d37d7f539b33c861d86ff0a3c0eaef1f688c20dd339b7a9c027de75ae0c`. Малые source/input hashes независимо пересчитаны и совпали identity. Imported accepted Core e171/tree742 clean; Julia1.13.3/4threads/BLAS1; actual cap12.125GiB/high11.5GiB/swap0/CPU400/FSIZE5GiB записаны source gate.

Process log exit0,146.48s wall (2:26.48),271.53s user+17.14s system,197%CPU,peak6662832KiB,swap0. Harness134.714s включает его own phases после imports/source setup; эти времена не смешивать. Summary `execution=complete`, `reduced.status=pilot_call_cap`, `scientific_acceptance=false`, data2429909590B. Это штатный конец16Π, не timeout/OOM и не внутренний scientific pass. Initial79small checks прошли; attempt01 programming failure сохранён отдельно и не исчезает из campaign cost.

Рецензент не запускал Julia/Pkg/новую карту и не читал большие HDF arrays. Сохранены `reviews/gamma_pilot_result_audit.py`, `gamma-pilot-result-audit.json`, `gamma-pilot-result-audit.csv`: stdlib пересчёт immutable TOML/source hashes,0.010s,Python17288KiB. `h5ls -r reduced-terminal.h5` использован только для metadata:Γ/rawΦ/GR/GL/GG/A имеют disk shape(5,5,11,30704), соответствующий Julia(30704,11,5,5); condition/scale disk(11,30704). Значения массивов независимо не реконструированы.

## Initial actual consumer gate

`actual-one-map.toml` completed,passedtrue,worstrelative5.04733716058e−16. Static Γ reduced-vs-native comparison≤5.05e−16, LO comparison0 (общий native LO primitive); native sumΣR vsPI(rawΦ)relative2.35865e−16; SAME-input cavity GR dummy-vs-physical-FD discrepancy0. Per-channel raw KMSrelative7.20e−17–9.10e−17, weighted7.20e−17–8.76e−17, maximumblockrelative1.22e−15. Absolute/scale/locator primaries сохранены. Заданный1e−10 пройден без изменения criterion.

Это подтверждает actual consumer equivalence на **seed A/current μ=−17.7368846meV**, включая static4+LO, explicit PI/Hilbert и finite-chain4. Совпадение shared native primitives не является независимым oracle их physics и не uniform KMS certificate последующих Γ/μ. Stencil minorweight4.54e−12/coefficientbound2.22e−16 согласованы с nearly integer shift; weighted-adjoint defect.00409869 остаётся отдельным quadrature/endpoint witness и не становится conservation pass.

## Π trajectory и guards

Независимый TOML audit нашёл ровно maps001–016 и19attempt records: acceptedseed1, acceptedlinear4 (maps2/4/6/8), acceptedconvexfallback3 (maps3/5/7), acceptedAAII8 (maps9–16); AA cone rejections3. Rejections обнаружены **до actual_call!**, поэтому не являются failed evaluated Π; выполняемые fallbacks counted. ИтогΠ16/precheckoperations3 согласован. Precheck3 включает native operators; `native_accounting=0` обозначает отсутствие **native exit stage**, не отсутствие вообще native contractions.

Отрицательные eigenvalues rejected AA trials −8.79e−8,−1.29e−6,−4.75e−8 существенно больше arithmetic allowance7.11e−14 по модулю: cone guards действовали содержательно, clipping не применялся. Максимальный accepted L1=8.7661<10,QRcondition11.205<1e8. History/streak reset после fallback согласован frozen implementation. **Stdout kind=AAII на3/5/7 обозначает исходный trial, а не принятый шаг**: окончательная классификация берётся из attempt TOML. Это naming limitation, не изменение численного пути.

| map | raw relative | weighted relative | фактически принят |
|---:|---:|---:|---|
|1|1.928505285|2.049486268|seed|
|8|.417546533|.440725633|linear|
|9|.160496173|.151865451|AAII|
|12|.0568918055|.0526466866|AAII|
|16|.0203879940|.0186891421|AAII|

Raw numerator уменьшается30.2448969→.353067206, denominator15.6830770→17.3174078. Улучшение relative94.5902×,absolute85.6633×; weightedrelative109.6619×. Поэтому эффект не объясняется только ростом denominator или малым mixed step. Published raw residual и working Herm residual совпадают доroundoff; maximum representation relative9.39106e−17, намного ниже наблюдаемой невязки. Проекция representation не скрывает текущую ошибку.

Observed raw relative строго убывал по этим16 accepted inputs, хотя алгоритм nonmonotone. Из этого не следует contraction constant, linear asymptotic rate или ошибка расстояния до решения. Финальная raw невязка примерно2.04e7 раз выше1e−9 trigger, weighted1.87e7 раз; требование3 последовательных passes не приблизилось к выполнению. Позднее замедление/.0204 floor нельзя достоверно объявить stagnation за этот короткий отрезок.

## Физический смысл и неизмеренное

Final common μ=8.463714668meV, density-fit relative−1.43250e−13, maximum absolute fit error по16maps1.54506e−13; capacity конечна/feasible, derivativepositive. Это **конструктивное условие restricted map**, не независимый raw Keldysh λ, chemical-potential convergence или evidence thermalization. GL/GG задаются ifA/−ihA на каждойΠ; их pair identity, наблюдаемая≈7.06e−13 absolute в конце, также не доказывает unconstrained FDR.

All input/A/map cones прошли заявленный floating backward criterion; финальные tail node minima Γmap2.87e−17 иA2.53e−15 меньше allowance7.11e−14. Это не formal strictPD certificate, не основание исключить atoms/dark states/дискретизационные ошибки. Maximum sampled Dyson condition вырос193.2→1276.8; finite inverse на этих узлах не доказывает resolution между узлами. Total Γ residual особенно не bound individual weak-channelΣ residual: вΦ norm доминирует LO16.5195, тогда как alloy.1677/acoustic.3286/IFR1.0791/impurity1.3140.

Native stage **not_measured reason=no reduced trigger**; materialization/solver/audit events все0, started/completed bounds0. Поэтому нет новых native rK/rΣ/rawλ/candidateλ/rD/rA, per-component RLG/embedding± fixed-point evidence, collision/power/FDR audit, current или strict inner pass. Это корректное отсутствие измерения, не fail native equations и не pass. Нет Poisson/U update, discretization validation или experiment comparison. Сходство μ с прошлым FD64 не устанавливает совпадение решений: grid m939/PI closure/FD restriction отличаются прежнему baseline.

## Решение

Сохранить attempt02 как bounded **raw residual reduction**, не acceptance. Данные показывают работоспособность actual reduced consumer/control и полезное продвижение от broad seed, но оставляют отдельно nonlinear approach и spectrum/quadrature quality. Дешёвый следующий различающий контроль может фиксировать terminalΓ/A и измерить density/spectral quadrature sensitivity; он не принимает SCBA. Если исследовать другой warm seed/продолжение, нужна отдельная identity и тот же raw/native acceptance contract, без ретроактивной смены cap16 или допусков этого attempt. Просто увеличить iterations без различающего вопроса эти16 результатов не обязывают.
