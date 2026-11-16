# Matsubara V4: frozen UNPROJECTED Γ-last64 + actual saved h_H

Статус: write-only, actual запрещён. Новый каталог; V3 e45ac490…, parent5249…, прежний Γ128-HOLD и все старые templates/frozen files неизменны. Новый источник — stableA-cont64: terminal6eb3…, checkpointΓ-last cab810…, spectral-selected06442fe…, raw model700ae…, writer50bc…, stable_spectraleeca…. Source file hashes взяты из root `stableA-cont64-source-digests.json`, не вычислялись через чтение HDF5 этим агентом.

**Единственное изменение модели V4 относительно V3:** h для всех9 cells берётся из ACTUAL сохранённого `spectral-selected-064.h5/h_H`, а не raw model h. Loader проверяет saved h_raw==model h_raw побайтно, spectral grid/E0/Tplus/ledger identity, all11×5×5 h_H finite и EXACT Hermitian без norm tolerance/проекции. Exact equality относится к точным бинарным числам; ±0 — один точный ноль, как в PD proof. Исходные h_H numbers не усредняются. Tminus обязан быть EXACT adjoint Tplus. Γ terminal и Γ-last сравниваются по всем real/imag entries slabs≤1024E через `isequal`, включая signed-zero provenance.

Численный backend не менялся: Cauchyhat coefficients/BigFloat controls, full45 LU central5 columns, Herm(G), moments/M1/M3 tails, N512/1024, Γ=0 controls, sum/checkpoint арифметика скопированы из e45/5249. Hash prefix1:177 равен parent; suffix начиная с independent controls побайтно равен V3. Выбор h_H и ограничения источника — новая явно названная scientific identity, не приёмка Γ128/native старого raw-h пути.

## Freeze-IO (можно до PD certificate, только Root после review)

`input-stable64.toml` содержит actual file hashes/Core/env/μ/units и `execution_authorized=false`; 13 logical-array SHA (hcell теперь saved h_H), 11 per-k A SHA, geometry и raw-h SHA пока UNSET. Новый `freeze-stable64-io.jl` может извлечь их, не читая certificate и не снимая HOLD. Он проверяет whole-file SHA до HDF5, source writer/input/identity/map64/ledger/summary ownership, strict h_H policy и full checkpoint-terminal Γ equality. Один большой массив за раз; A по одному k; только small digest output. `h_raw_sha256` отдельный от 13 array fields. Никакой PD/плотности/SCBA приёмки helper не выполняет.

Root scope: ≤2GiB, Julia1/BLAS1, external120s, output32MiB (internal105s включая SHA):

    julia --startup-file=no --threads=1 --project=.build/campaign133/environment/accepted \
      .build/campaign133/adjacent/gamma-density-hermitian-v4/freeze-stable64-io.jl \
      .build/campaign133/adjacent/gamma-density-hermitian-v4/input-stable64.toml NEW_DIGEST_DIRECTORY

После completed digest переносить только поля из совпадающего small digest, freeze нового input SHA и independent binding review. Никакой автогенерации actual разрешения.

## Жёсткий PD gate actual

Принимается только ORIGINAL native **TOML** certificate b0 Bareiss, плюс exact hashed proof input. JSON не требуется. Пути/SHA остаются UNSET до настоящего completed ALL. V4 actual main останавливается на execution_authorized=false ещё до HDF; после будущего Root green дополнительно требует:

- `execution=complete`, `mode=all`, `all_blocks_strict_PD_certified=true`, expected/requested/checked/PD=337744 и shape30704×11×5×5 в proof input; call64;
- zero nonpositive/nonhermitian/nonfinite counts, empty failures, отсутствие failure поля, 20264640 exact component divrem checks;
- checkpoint/model/input/identity/operator/spectral/ledger/summary/map SHA должны совпасть одновременно в certificate, proof-input bindings и V4 source packet; проверяются runner/math hashes и proof input SHA. Source checkpoint — ровно cab810…;
- после certificate Γ terminal всё равно сравнивается с certified checkpoint по ВСЕМ entries. Ни benchmark5000, ни partial pass, ни selected fixtures не принимаются вместо ALL.

`certificate-gate-controls.jl` определяет source-free отрицательные metadata tests (executionfalse, benchmark, incomplete/failed, requested/checked subset), не настоящий proof. Они здесь НЕ запускались. Root может include V4 + controls и вызвать `certificate_gate_controls(require_all_PD_certificate)` без HDF5 reads.

## Доказательная область и предел

Если ALL доказал PD exact stored Γ_j, convex full hats дают Γ̃(E)⪰0 (внутри support PD), а для Im z=y>0:

    −Im_H Σ(z)=∫ y Γ̃(E)/(2π|z−E|²)dE ⪰0,
    Im_H D45(z)=y I45−I9⊗Im_HΣ(z) ⪰ y I45 >0.

При EXACT Hermitian H45 из saved h_H и adjoint couplings D invertible во всей UHP; верхних полюсов нет. Это опора причинной Cauchy-hat frozen numerical model, не утверждение о SCBA или forward error программы. Возможные real-axis atoms формула Matsubara учитывает; их отсутствие этим PD gate отдельно не заявляется. Нормировка G~I/z и finite μ/T — условия equilibrium spectral functional. Реальные finite-series/Float64/LU/tail ошибки требуют сохранённых witnesses и refinement; agreement N512/N1024 не formal total-error certificate.

Будущий actual scope после ALL+review+Root green: ≤2GiB / Julia1-BLAS1 / CPU200% / external180s / output64MiB; internal160s; выбран M3 N1024, N512 — witness. Никаких новых Φ, μ fit, Γ projector, h replacement beyond declared saved-h_H policy, Pkg или solver campaign. Сейчас ни helper, ни actual не выполнялись автором.
