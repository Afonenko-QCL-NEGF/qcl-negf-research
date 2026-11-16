**Независимый result review selected intended split Big oracle**

Решение: **selected44 input-PD proof и numerical covariance accuracy
подтверждены**. Для этих4 retained energies×11 k positive-source Acov
соответствует intended all9 hH/splitΛH/Γ model существенно точнее Asub
в слабых хвостах. Это не full-grid causality/SCBA/density/experimental
certificate и не nonlinear acceleration result. Активный cont64 данным
review не изменяется.

Evidence: `positive-source-selected-big/attempt01/summary.toml`, traces.csv,
critical.h5 и actual log; source2d946fd46e35bdb73d9347d8a5487cbaa7a6491072eb94ebd2e86d4ab0df87f1,
summary SHAe0f4622495195f500763b4479f4daffe0477ed842a5f4dbb9c0bd10d68f7f405,
traces CSV SHA2f1af5b2613957d9395b1e748095bb84a8d7666a9a327d2d6eb26e92387a9b8e,
runtime inputaf94ef28…f3e. Source preflight selected1.3MB SHAa74dd7cb…079a
связан с original preassembly hH/ΛH/Γ, не со старым raw-h obstruction case.
Root tiny13/13 pass,4.21s/514744KiB/Swap0; actual44 blocks завершены,
3.89866s script/9.30s process,629816KiB/Swap0/exit0, declared2GiB/120s.
Reviewer Julia/Pkg/Big inverse не запускал. Для независимого proof audit
прочитан только selected1.3MB через existing h5dump, малые TOML/CSV/log;
whole-grid Γ/model/SCBA arrays не открывались.

**Независимая проверка exact premise**

`reviews/positive_source_big_result_audit.py` — stdlib-only Python с
RLIMIT_AS128MiB. Через h5dump17 significant digits получены source Double
real/imag entries; Float roundtrip затем переводится в exact
`fractions.Fraction`, без промежуточных complex arithmetic округлений.
Проверены axis order selectedE/k/a/b и expected indices7893/8038/30696/30704.
Все44 Γ exact Hermitian; независимо выполнены no-pivot rational LDL,
**все220 pivots строго положительны**, каждый равен сохранённому rational
numerator/denominator/sign. Полные reconstruction всех25 entries в каждом
block точно равны source Γ. hH и preassemblyΛH exact Hermitian также
сверены непосредственно по выбранным данным.

Это strictPD proof **этих exact binary input matrices**, сильнее
floating eigensign/allowance. Он не ограничивает uncertainty original
physical Γ или остальные NE×Nk blocks. Для selected intended finite-chain
D с Hermitian hH/ΛH/T links и ImD=I9⊗Γ/2≻0 имеем nonsingular D, а exact
central-row source congruence положительно определена. Actual inverse
по-прежнему численный, поэтому precision/backward evidence рассмотрено
отдельно. Independent audit JSON:
`reviews/positive-source-selected-big-result-audit.json` содержит44 proof
checks и SHA всех использованных малых источников.

**Big192 response и Double comparison**

Все44 response records присутствуют с правильными e/k. Скрипт использует
adjoint5-RHS solve D†Y=P → central rowY† и source sum по всем9 cells,
то есть проверяется intended full45 model без native cavity recursion.
Big192 metric absolute/scale/relative decimal ratios независимо
пересчитаны (176 проверок); trace CSV44 строк совпадает с primary TOML.

| Величина, максимум44 | Relative error/witness |
|---|---:|
| Double Acov vs Big intended Aprod | 6.8368227e−15 |
| Double GR vs Big intended GR | 3.4542522e−15 |
| Double Asub vs Big intended Aprod | 6.5662380e10 |
| Big192 Asub vs Big192 source product | 2.7888870e−31 |
| Big192 adjoint-LU backward residual | 8.5411711e−59 |

Worst covariance node e8038/k1: intended ||A||F=4.8348184e−12,
Double Acov absolute difference3.30548e−26, Asub absolute9.52465e−15
(relative0.0019700). F64 condition proxy for downcast intended D249.665;
он не сертифицированный Big forward-error bound.

Worst subtraction node e30696/k4: intended ||A||F=1.1553521e−28,
Asub absolute difference7.5863169e−18, поэтому relative6.566e10.
Covariance relative там1.71318e−16. Огромное relative value относится к
ультраслабому spectral matrix, не к большой sheet-density/current error.
Double-vs-intended comparison включает split-vs-assembled ΣR, native
cavity/inverse и A arithmetic; оно не изолирует одну LU или финальное
вычитание. Strict floating Gram property materialized Acov не заявляется.

**Один384 repeat и q aggregates**

Повтор выполнен ровно на e8038/k1 — **worstCov** из44 по заранее заданному
selection rule, а не на всех44 или на worstSub. Independent Decimal
пересчёт trace192/384 differences подтверждает relative8.29134e−57 для
product trace и4.28885e−56 для subtraction trace. Big384 source identity
relative6.37252e−103. Это сильное numerical precision evidence в этой
точке, не interval certificate и не384 qualification остальных blocks.

q_E заново агрегированы независимо в Decimal110 по всем11 k для каждой
selected E, с exact promoted **source Double** gs/(2π),wE,wk и сохранёнными
Big traces. Stored192 aggregates воспроизведены relative≤3.981e−58.
Это4 spectral number weights, а не density integral по всем E и не fit μ.

| E-index | Relative error stored Double qcov vs Big q | qsub vs Big q |
|---|---:|---:|
| 7893 | 7.03e−16 | 4.96e−16 |
| 8038 | 1.66e−15 | 2.16e−14 |
| 30696 | 1.79e−17 | 8.65e−8 |
| 30704 | 5.48e−17 | 4.93e−7 |

Последние q всего~9.22e−32 и4.59e−32 в native number-weight units;
их lost relative precision не доказывает заметной ошибки total density.
Ни reconstructed q, ни source Γ/A не нормировались/не clipping.

**Вывод для кампании**

Qualified narrow claim: covariance construction сохраняет intended
selected weak-tail matrix/number response в Double на проверенном seed,
тогда как Asub может терять relative accuracy даже при хорошем GR.
Точные selected Γ PD, Big backward/source identity и single384 repeat
разделены как разные evidence. Whole-grid/future Γ acceptance из них
не следует; raw Φ kernels, nonlinear residual, discretization и original
raw-h native SCBA gates остаются независимыми. StableA16 уже показал
практически прежнюю16-map residual trajectory, поэтому этот oracle не
добавляет заявления об ускорении или принятом fixed point.
