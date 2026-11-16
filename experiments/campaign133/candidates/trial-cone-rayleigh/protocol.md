# Selected rejected-trial Rayleigh diagnostic

Статус: write-only proposal1. Автор Julia/Pkg/HDF5 runtime не запускал;
Core/solver/trial arrays не менялись. До запуска Root проверяет код и tiny
controls. Exact script/input SHA — в SHA256SUMS. Cold/warm/affine outcomes
сохраняются; этот diagnostic не расходует Π, не принимает convergence.

Вопрос: существует ли **точно отрицательная quadratic form** у сохранённого
вычисленного AA proposal, или знак Float64 eigensolve не подтверждается
независимой arithmetic? Это не проверка точности всей Φ, не доказательство
ошибки eigenvalue для continuum model и не разрешение ослабить cone guard.

Default blocks: rejected-trial5/e5090/k1 и rejected-trial16/e9771/k1.
Input pins small attempt TOMLs и parent input. Root предварительно streams
whole-file SHA каждого trial и фиксирует их вместо UNSET в authorized input
copy; diagnostic сам повторяет SHA до HDF5 field read. Default input имеет
runtime_authorized=false и intentionally UNSET large-file hashes, поэтому
его нельзя запустить как actual calculation без concrete source freeze.
При необходимости Root может отдельным input добавить остальные5 known
rejections: максимум7 blocks. Никакого auto-scan/current source selection.

HDF5 читает только logical(1,1,5,5) slab Gamma для заданных e,k. Проверяются
checkpoint schema/identity/call и logical shape(30704,11,5,5), source attempt
negative-PSD reason/indices/call. Дисковая C orientation(b,a,k,E) не выдаётся
за Julia logical(E,k,a,b). Full-file hashing streaming, без full Γ allocation.

Float guard matrix B[a,b]=(X[a,b]+conj(X[b,a]))/2 строится тем же выражением
и b,a loop order, что frozen cone. Float eigenvector minimum eigenpair берётся
через LinearAlgebra Hermitian/LAPACK. Это candidate vector, не oracle sign.
Raw Frobenius/anti-Hermitian norms и allowance вычисляются literal formula
accepted Core e171 metricv3:64*n*eps(Float64)*max(Frobenius(X),1); allowance
должен exactly совпасть записанному rejection reason. Выводятся Float λ,
Rayleigh quotient, eigenpair absolute/relative residual/scales.

Независимый sign witness — Complex{Rational{BigInt}}. IEEE754 bits всех Float64
entries и Float eigenvector декодируются в **точные binary rationals**, без
rationalize/decimal approximation. Вычисляются два нормированных quotient:

1. v† B_float v / v†v — для точных сохранённых значений Float-averaged guard B.
2. v† Herm_exact(X_float) v / v†v — точное усреднение exact binary raw entries.

Разность отделяет averaging roundoff от выбранного vector witness. У exact
Hermitian формы imaginary part обязана быть exactly0. Negative quotient
доказывает non-PSD **этой exact matrix**, без необходимости exact eigenvector;
он не является eigenvalue enclosure или полным forward-error bound Φ.
Сравнение quotient с−declaredallowance — только diagnostic same cone scale.
Предполагаемый negative sign не используется как assertion: nonnegative
witness сохраняется как самостоятельный результат, а не скрывается failure.
Exact numerator/denominator/sign и 256-bit decimal display сохраняются вTOML;
малыйHDF содержит raw block/guard B/Float eigenvector. Decimal не заменяет
exact rational первичное свидетельство.

Root controls: один Julia1/BLAS1 process, heap hint1G, OS memory2GiB/swap0,
CPU100%, external timeout120s, file8MiB. Внутренний110s check между blocks не
заменяет external timeout. Ожидаемый output менее1MiB; максимум7×5×5 blocks,
никакой full-state/Green storage. Existing accepted environment предоставляет
HDF5 и stdlib; нет Pkg/install/Core/Runner imports. Script wall стартует после
imports; total process wall/RSS Root измеряет отдельно.

Tiny command (не исполнен автором):

    julia --startup-file=no --threads=1 --heap-size-hint=1G \
      --project=.build/campaign133/environment/accepted \
      .build/campaign133/trial-cone-rayleigh/rayleigh.jl --tiny

Actual command после Root review + authorized input SHA freeze:

    julia --startup-file=no --threads=1 --heap-size-hint=1G \
      --project=.build/campaign133/environment/accepted \
      .build/campaign133/trial-cone-rayleigh/rayleigh.jl INPUT_RUN.toml NEW_OUTPUT

Tiny controls: exact binary0.5/nextfloat1/minsubnormal/signedzero/nonfinite,
independent complex2×2 negative vector, guard Hermitian expression, exact raw
Hermitian part, nested TOML rational fields и selected HDF slab roundtrip.
Checks written, not pass. Root actual outcome/review нужны до научного вывода.

Test-only V2: исходный tiny run V1 дал11pass/1fail,4.09s,477860KiB,exit1.
Тест ошибочно требовал numerator−1 для quotient **LAPACK minimum eigenvector**;
известный−1/2 относится к другому, вручную заданному v=[1,i]. Manual exact
oracle остаётся unchanged. Исправлена только TOML roundtrip проверка: exact
numerator/denominator должны совпадать исходному rec. Все math/IO/runtime
functions unchanged. V1 source/input/SHA/log/admission/authorized input
сохранены до правки вunexecuted-v1/. Actual selected trial HDF5 run не был
выполнен. Root повторяет tiny после нового source freeze, затем reuse того же
authorized input с уже frozen whole-trial SHAs. Автор runtime не запускал.
