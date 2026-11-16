# Independent result: frozen Γ64 Hermitian-all9 Matsubara V4

**Вывод:** при одном fixed Γ-last64 и неизменённом μ V4 воспроизводимо
показывает расхождение Matsubara M3 density и сохранённой finite-grid
`∫fA/(2π)`: +9,3376849e−5 target по числу,1,1832969e−4 по Frobenius матрицы.
N512→1024 существенно стабильнее этого расхождения. Принимается это узкое
численное свидетельство для intended exact-Hermitian-all9 frozen модели;
не принимаются SCBA convergence, точная величина continuum quadrature error
или утверждение об источнике расхождения исключительно из ghosts/atoms.
Blocker для такой формулировки не найден.

Первичные данные: `adjacent/gamma-density-hermitian-v4/attempt01/`
summary SHA323cc1773423a3b5a97c1af383acffff229a7c5fdb3b72c3e211ad1da7047644,
input f6560732… и archived script814dff50… совпадают с recorded hashes.
Самостоятельный малый audit `reviews/gamma_density_v4_audit.py` читает
TOML/certificate и только small density HDF (875KiB), не большие source HDF.
JSON `gamma-density-v4-audit.json` и точные h5dump commands/data в
`gamma-density-v4-small-h5dump.txt` сохраняют пересчёт. Выполнение audit:
0,25s/18912KiB/exit0; Julia, Φ и новые solver maps не запускались.

**Source/causal prerequisite.** ALL exact Bareiss certificate a703ec4c…
с proof input27da5a97… подтверждает337744 PD blocks,0 failures,
20264640 exact component divrem. Все9 source bindings certificate/proof
совпадают с V4 packet; checkpoint cab81029… — именно unprojected
stableA-cont64 Γ-last64. V4 требует exact equality всех real/imag Γ entries
между certified checkpoint и terminal, а freeze-IO completed digest также
содержит эту проверку8443600 complex entries. Сам аудит не повторял большой
pass; проверены малые proof/digest hashes и соответствующий loader code.
Прежний Γ128 и projected64 в этот вывод не входят.

Loader берёт **actual saved h_H**, проверяя exact Hermiticity всех11×5×5
entries без нового усреднения; saved h_raw связывается с model только как
контроль. `chain_hamiltonian` ставит этот h_H во все9 cells, T/T† — в
обе стороны; loader требует exact Tminus=Tplus†. Γ содержит5 local channels,
embedding исключён из повторяемого локального источника и реализован chain.
Для точного Cauchy-hat Γ-интерполянта с неотрицательными hats PD узлы дают
PSD spectral measure; exact-Hermitian H45 тогда исключает UHP poles
посредством Im(zI−H−Σ)>0. Это причинная **intended numerical model**, а не
сертификат точности вычисленного Σ/G или физической модели QCL.

**Нормировка и числа, независимо восстановленные из small HDF.**
ρ(k)=I/2+(2/β)Σ_{n≥0}Herm Gc(μ+iω_n), с M1/M3 asymptotic tails;
Herm — `(G+G†)/2`, сохраняющая комплексные offdiagonal элементы.
β из70K/E0=0,1eV, μ=0,08459614262982018 (8,459614meV) взят из source,
никакого fit/renormalization. Число `gs·Σk wk·trρ`, gs2, source wk;
target0,045. Frobenius comparison использует все11×5×5 entries без wk;
его нельзя приравнивать weighted density-number criterion.

| Свидетельство | N512 | N1024 |
|---|---:|---:|
| M3 N/target |1,000093421090491|1,0000933768485263|
| M1 N/target |0,9999097140597122|1,0000704134368297|
| `||ρM3−ρgrid||F` |3,34687095382e−5|3,34687094021e−5|
| `||ρM1−ρM3||F` |3,11449182276e−6|3,89312034759e−7|

Сохранённая grid density даёт0,04499999999999452 и согласуется с source
number fit; это binding/normalization control, не quadrature acceptance.
M3 N1024−N512 =−4,42419647e−8 target;
`||Δρ||F=9,71240560e−10`, relative3,43345111e−9.
Grid discrepancy по числу примерно2110× этого N-step witness.
Его dominant matrix contributions приходятся на low-k slices; обе нормы
сопоставляют один Γ, один μ и один h_H policy, не два SCBA решения.

**Ограничения точности.** M3 refinement и M1/M3 separation — witnesses,
не formal total tail bound. Maximum solve backward error4,86319e−17 не
даёт forward density error без conditioning/error propagation. Нулевая
Herm-part defect конструктивна, selected ImΣ eigenvalue maximum0 лишь
проверяет выбранные frequencies. Аналитический omitted-series Σ bound
для трёхчленного hat expansion не включает Float64 coefficient evaluation,
BLAS summation, LU и дальнейшую density conditioning; BigFloat controls
проверяют отдельные coefficients. Source canonical geometry отличается
от saved energies на1,78e−15 scaled units — это также не formal zero-error
тождество. Scalar/complex/source-bare controls passed поддерживают
конвенции/знаки/tails, не заменяя полную оценку ошибки interacting case.

Поэтому расхождение — сильное evidence недостаточности **данной saved
real-axis density quadrature как точного представления frozen analytic
functional**, с указанными numerical witnesses. Его нельзя без total-error
bound назвать доказанной чистой discretization error известной величины.
Matsubara functional учитывает возможные real-axis spectral atoms; отличие
от finite-window grid само по себе не локализует их вклад. Все ghost/atom
данные V4 — conditional floating witnesses, не interval certificate;
ALL Γ-PD не сертифицирует их endpoint45 bounds или отсутствие внешних атомов.

Actual execution completed:24,7846s script/28,54s process,917176KiB,
swap0/exit0. SCBA convergence/native acceptance, Poisson, energy/basis/window
convergence и эксперимент не измерены. Результат касается V4 frozen Γ64
Hermitian-all9 source; старую Γ128 density задачу он задним числом не принимает.
