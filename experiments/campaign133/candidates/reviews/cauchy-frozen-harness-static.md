**Независимое static review frozen Cauchy harness, 2026-11-16**

Решение: **блокирующей source/index/sign/API-binding ошибки по чтению не
найдено**. Wrapper подходит для source-free preflight, затем ограниченного
actual transform-only probe по решению root. Это не разрешение от reviewer
на запуск и не утверждение о measured accuracy/RSS/time. Julia, Pkg,
большой HDF5, scientific compute и сторонний код reviewer не запускал;
исходники harness/backend не менял. Записан только этот review.

Рассмотрен ровно `adjacent/cauchy-tree-prototype/run_frozen.jl` SHA256
`2e5777af5e21657e793cc5fe529edb2981b1d85816b870265af5006958137765`,
156 строк; `frozen-input.toml`
`546a2f57c7bfa6e51f25c69c05cfd62c48117f9073fef42a14cc58647a5da63a`,
`frozen-protocol.md`
`96120c301d2e86668a25cb43b2fa06b38c7abfde115d75bdc15a993c47cb47b1`.
Эти малые SHA перепроверены. Backend `a74ba641…` и tiny tests `3be8f0d8…`
также совпадают с FROZEN_SHA256SUMS и предыдущим independent math review
`reviews/cauchy-tree-prototype-static.md`. Математика всего tree повторно
не доказывалась; предыдущие 67 tiny tests относятся к root run.

Superproject при review — b515c42e…, существующие dirty Core/Runner и
research gitlink 69f1ef4 сохранены. Accepted Core worktree фактически clean
`e171caa60cce676c6c4238dd671e7c527f2bd977`; код comparator прочитан на нём.

**Проверенные границы wrapper**

| Утверждение | Первичная привязка и вывод |
|---|---|
| Source-free preflight предшествует source | run_frozen:2–7 включает source-free tiny tests и связывает CT/PI aliases; `--preflight` на 51–52 возвращается до cfg/output/source reads. Guard 154–156 — явный if, без прежнего macro `&&` parse trap. Actual импорт/API всё равно not_measured до root preflight. |
| Writer/input/source provenance | 87–95 сверяют SHA prior summary, input, identity, writer и marker; identity input/script и fixed μ проверены отдельно. На 96 source H5 SHA проверяется ровно один раз, до единственного read-only `h5open` на 100, а не трижды для k1/6/11. Reviewer этот 3.8GB файл заново не читал/не хешировал. |
| Γ из пяти LOCAL families | 12–22 явно суммируют IFR/LO/acoustic/alloy/impurity в порядке прежнего PI control; shapes `(30721,11,5,5)` проверяются для real/imag retarded datasets, hyperslab выбирает один k. Формула `i(R[e,a,b]−conj(R[e,b,a]))` — matrix anti-Hermitian part с правильным знаком. Embedding/± не читаются; никакого Hermitize/clip/repair. |
| Axis/scales/basis | 101–110 проверяют v3/E0=.1eV, U=0, NE30721, семь source basis operators по dimension и persisted SHA из закреплённого prior summary. χ/Φ действительно Nz×Nb: accepted basis builder задаёт `χ=Φ/sqrt(Δx)`, types.jl:556–577. E finite/order/ghost order проверяет tree; accepted PI дополнительно требует exact canonical-range values. Positive finite length-matching wE проверяет native transform. Полные k/x/scales не rebuild: для этого transform они не интегрируются; source-file SHA связывает их с known artifact. |
| Aggregate и slice identities различены | 97–98 проверяют aggregate Γ reference из prior summary и прямо не называют его recomputed digest selected k. Actual gamma_k_sha256 и energy_sha256 пишутся на 118. Это корректная ограниченная provenance. |
| Accepted Core aliases/layout | `QCLNEGF.QCLNumerics` действительно включает fft_parallel.jl на QCLNEGF.jl:652. `_ProductionPIFFTHilbertPlan` / `_production_pi_fft_hilbert_transform` определены на fft_parallel.jl:333–399. `(NE,1,5,5)` имеет 25 columns, plan(E,25) совпадает. wE нужны validation, не дополнительный PI factor. Returned Λ — без boundary imaginary term; wrapper на 131 ровно один раз вычитает `iΓ/2`. |
| Literal и canonical geometry разделены | Tree использует stored knots и literal adjacent ghost gaps. Core PI canonical origin/spacing и Toeplitz endpoint hats. 117–123 сохраняют обе geometry metadata, actual/ideal ghosts и их max difference; protocol прямо запрещает bitwise/continuous-interpolant equivalence claim. |

**Независимый off-grid oracle**

Строки 25–48 используют second divided difference `F(t)=t log t`:

`w_j=[(F(z−a)−F(z−c))/(c−a)−(F(z−c)−F(z−b))/(b−c)]/(2π)`.

Это правильный integral положительной hat против `(z−x)⁻¹`, включая обе
endpoint hats с нулевыми ghost nodes. Для real PV используется
`t log|t|`, в t=0 берётся непрерывный предел F(0)=0. Для Im z>0 principal
complex log имеет retarded branch. Big192 применяется к scalar coefficients
на всех source+ghost knots, затем weights округляются до ComplexF64;
сама matrix contraction обычная Float64. Поэтому это независимый
formula/arithmetic witness, а не rigorous Big192 matrix result.

`transpose(weights)*reshape(G,NE,25)` на 37–38 принципиально не conjugates
weights; reshape сохраняет column-major (a,b). Reverse contraction на 39
меняет accumulation order, не оператор. Cache ключ `(Re z,Im z)` допустим:
E/ghosts неизменны между k, τ; real targets и complex target не смешиваются.

Real off-grid interpolation: x=[leftghost;E;rightghost], interval j имеет
left=G[j−1] (или zero при j=1), right=G[j] (или zero у правого ghost).
Это соответствует actual hats; boundary `−iΓ(z)/2` добавляется один раз.
Complex target на 127 получает ordinary Cauchy без лишнего boundary term.
От source-free tiny oracle здесь не используются backend panel/moment
functions, native Core Hilbert или adaptive quadrature.

**Scope, checkpoints и пределы**

При полном завершении ровно 3 k × 2 τ × (16 source nodes+4 off-grid)=120
matrix comparisons на одном и том же selected source field для каждого k.
Четыре off-grid targets: два объявленных midpoint, left ghost-strip
midpoint и μ/E0+0.2i. Native FFT вычисляется на полной original NE для
одного k, tree queries — только declared samples. Обе τ используют те же
G/tree/reference; изменения источника между ними нет.

Checkpoint перед/после queries и между фазами сохраняет full5×5 value/ref,
raw difference, allowance, ledger, roundoff witnesses, timing/RSS. Gate
far_estimate+max(1e−11,1e−10||ref||) — diagnostic tolerance с размерностью
source Σ, не formal whole certificate. NaN difference не проходит `<=`;
backend проверяет finite evaluator values, tree проверяет raw G. Наличие
reverse-order witness не превращает shared arithmetic errors в bound.
Ни один green transform gate не принимает density/SCBA/KMS/adaptive economy.

Фаза build_tree/full FFT/Big192 reference не имеет внутреннего wall-clock
interrupt: soft110s проверяется на checkpoint, а не каждая arithmetic op.
Это явно отражено в protocol; **внешний timeout120s и MemoryMax/Swap0 нужны
обязательно**, root выбирает 1.5 или 2GiB по preflight и пишет фактический cap.
Top-level tiny tests/imports происходят до internal timer, но покрываются
внешним timeout. RSS включает process lifetime. Подавать actual probe до
source-free root preflight оснований нет; переносить только потому, что
ресурсная стоимость ещё not_measured, тоже не требуется.

Малые metadata/code SHA проверены, source axes/basis/operator hashes
проверены **статически как guards**, не как новое HDF5 measurement.
Источник Γ и фактическая source API compatibility будут установлены только
root actual execution. Новых требований к платформе, solver или acceptance
review не вводит; конкретных обязательных кодовых исправлений до первого
bounded transform probe в рассмотренной frozen версии не найдено.
