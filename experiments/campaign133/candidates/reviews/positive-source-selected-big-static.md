**Независимый static review selected split Big oracle**

Решение: **одобряется root source-free tiny, затем один bounded selected
run после tiny/resource admission**. В frozen233LOC
`positive-source-selected-big/oracle.jl`
`2d946fd46e35bdb73d9347d8a5487cbaa7a6491072eb94ebd2e86d4ab0df87f1`
блокирующей math/sign/left-right/source-frame/serialization ошибки не
найдено. Input `a91219212ecd0408376ffeebbbf7bb9e78709e18250fb842064ce94048bfea2f`,
protocol `278e3e2541af4f93667ecbb32c41dacfccd3ab2bafb86490f5a0f7510f62556b`.
Reviewer Julia/Pkg/HDF5/numerical oracle не запускал; код автора не менял.
Проверены тексты, малые SHA и actual Julia Base rational conversion source.

**Frame и intended operator**

Source selected1.3MB, input, summary и preflight script связаны четырьмя
SHA guards до read-only HDF5. Prior должен иметь execution complete,
3started/completed и совпадающие code/input bindings. Expected energies
7893/8038/30696/30704, shapes4×11×5×5, H11×5×5 и T5×5 фиксированы;
finite fields и positive source quadrature metadata проверяются. Это
ровно44 selected blocks, не новый full-grid retarded pass.

В intended_oracle:66–85 saved hH, **preassembly** ΛH, Γ и T сначала
переводятся в Complex BigFloat указанной precision. Только затем
собирается D0=EI−hH−ΛH+iΓ/2 и full45 open-chain D с9 одинаковыми cells,
upper−T/lower−T†. E promote также до addition. h/Λ/Γ повторно не
Herm-проектируются и не clipping; exact Hermiticity отдельно проверяется.
Это соответствует positive-source preflight all9 hH, а не старой модели
с central hraw/eight side hH. LocalR/totalR не используются вместо split
inputs, что сохраняет tiny Γ, потерянную ранее при Double assembly.

**Central row и source covariance**

При P, выбирающем центральные5 columns, solve `D†Y=P` даёт
`Y†=P†D⁻¹`: именно **central ROW**, нужный для repeated-source covariance.
`G=row[:,center]` — правильный central retarded block. Использовать
напрямую centralCOLUMN `D\P` для source sum было бы неверно в complex
noncommuting задаче; здесь такой подмены нет.

`Aprod=Σcell row_cell Γ row_cell†` на82–84 равен
`P†D⁻¹(I9⊗Γ)D⁻†P`. `Asub=i(G−G†)` и adjoint-LU residual относятся к
тому же assembled Big D. Знак источника и i/2 correct; лишнего2π в
спектральном A нет. LU residual scale `||D||||Y||+||P||` безразмерно
нормирует inverse equation и не делит на vanishing A.

Double Acov/Asub/GR comparisons выполняются внутри active precision,
против intended Big model. Они включают split-vs-Double-assembled Σ,
cavity/solve и A arithmetic, а не только ошибку LU. F64 cond(Ddowncast)
явно proxy, не certified Big/exact condition bound. Execution completion
не является numeric accuracy gate: actual residual/consistency/error
records должны быть рассмотрены отдельно после run.

**Exact rational Γ qualification**

`Rational{BigInt}(Float64)` действительно exact binary conversion:
прочитан установленный `/usr/share/julia/base/rational.jl:175–178`,
где вызывается `rationalize(T,x,tol=0)` и проверяется обратное преобразование.
Tiny дополнительно проверяет exact rational0.1. Нет SVD cutoff/jitter/
eigenvalue repair; комплексные real/imag части переводятся раздельно.

No-pivot Hermitian LDL:24–58 использует
`d_j=Qjj−Σ L_jk d_k conj(L_jk)` и
`L_ij=(Qij−Σ L_ik d_k conj(L_jk))/d_j` в exact rationals. У imaginary
pivot проверяется точный ноль. Все пять строго positive pivots вместе
с exact reconstruction `Q=L diag(d)L†` доказывают strictPD данной
**exact binary selected matrix**: L unit-triangular invertible, d>0.
Проверка reconstructed==Q действительно выполняется, не inferred из
eigenvalue знака. Negative/zero pivot сохраняет no-certificate, не
объявляет автоматически indefinite matrix; non-Hermitian input не чинится.

Numerator/denominator/sign каждого pivot — строки/Int, допустимые для
nested TOML; decimal256 — лишь presentation, не источник exact sign.
Ни Rational/Complex BigInt объект в отчёт не протекает. Nonpositive
certificate может coexist с execution complete; смотреть нужно
`exact_PD_certified_count`, `all_selected_exact_PD` и отдельные records.
PD certificate относится только к44 retained Γ, не всему NE grid или
uncertainty исходного physical model. При exact-Herm h/Λ и PDΓ strict
Im D=(I9⊗Γ)/2≻0 гарантирует nonsingular intended D и PSD exact A;
approximate Big inverse всё равно требует residual/precision evidence.

**Precision, q и доказательная граница**

Все44 responses на192bit, primary trace/residual numbers сохранены decimal
strings. По одному заранее заданному worst **Double Acov relative error**
узлу выполняется384 repeat. Точки не выбираются по красивому sign/gain;
worstSub сохраняется отдельно, но не притворяется384 repeated. Сравнение
traces192/384 парсит строки внутри384precision. Это single-point precision
witness, не all44/full-grid384 certificate.

Per-selected-E q192 агрегирует все11 k из тех же trace products,
promoted source Double `gs/(2π)`, wE и wk. Source prefactor rounded
как в native number_weights_A, затем Big accumulation; measure явно
объявлена. Это q_E, не integral по всем E и не fitted density. CSV44
строки позволят независимый aggregation audit при доступе к source weights.
Critical HDF matrices intentionally F64 downcasts, не primary192 precision
arrays; Big traces/metrics strings остаются отдельными evidence.

**Tiny и лимиты**

Tiny содержит exact binary rational, известный complex2×2 PD pivot19/8,
negative и non-Hermitian controls; complex/non-normal three-cell Schur
oracle для centralGR/sourceA и phase-sensitive off-diagonal, weakΓ~1e−27
product/subtraction consistency, nested proof TOML roundtrip. Это
содержательное покрытие row/source/precision/serialization, не только
scalar commuting case. Его результат пока должен дать root; reviewer
не заявляет pass до actual tiny log.

Protocol задаёт J1/BLAS1/heap1G, hard2GiB/Swap0/CPU100%/120s/file8MiB;
internal110s checks между nodes не заменяют external timeout, который
покрывает imports/tiny. Script не импортирует Core/Runner, не строит
solver cache, не вызывает Π/SCBA/Pkg. Resource enforcement остаётся root.
Не найдено причины расширять код/проверки перед первым bounded run.
Scientific_acceptance=false и selected-only scope корректны; дальнейший
result review должен различать execution, exact input PD и accurate
covariance, не объединяя их в один pass.
