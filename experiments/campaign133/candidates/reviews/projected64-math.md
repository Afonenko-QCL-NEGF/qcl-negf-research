# Projected64: независимая математическая проверка frozen реализации

Область: статическое чтение исходников и доказательство формулы; Julia, HDF5,
численные опыты и изменения frozen файлов не выполнялись. Рецензия относится к
`equilibrium-gamma-projected64/trial_projection.jl`, а не к отсутствующему
`projector.jl`. **Блокирующих математических ошибок в проверенной дельте не
найдено.** Это не результат actual pilot и не доказательство его сходимости.

Проверены все 18 записей `SHA256SUMS`. Ключевые SHA256:

| Файл | SHA256 |
| --- | --- |
| `trial_projection.jl` | `cf9dc3c14e2d640884cfd7a7fef49c6137b1ad7681901b15457388677d905e7e` |
| `gamma_functions.jl` | `1c1eb3c3b5de030e21be5d23af0d74bf446781201e0cd12069a6f77d028027d4` |
| `stable_spectral.jl` | `eeca1f29710d9bb60bf6ce70ce0ebd1620fb86d4c77314b4bd8ffa9a5bcd5aa5` |
| `projection_checks.jl` | `ecaff8af11ef75da78d3e70fc6ad37b15e8f29c5dd0975dbcfc8c17725a73c65` |
| `protocol.md` | `b3d3e3e93398260af4ec5fa14aadc34e3256a5b4858665d3afc9987697e69e45` |
| `cont64-to-projected64.patch` | `39ae7a26bf3e489af79a08a4401300000b2ee651acc62c46a048b3f5342b18f4` |

Core задан `e171caa60cce676c6c4238dd671e7c527f2bd977`, tree
`742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8`. Проверяемый контроль начинает с
того же pinned positive Γlast16, empty AA history и fresh μ, что unprojected
cont64; таймеры явно увеличены 600→900 s и overall1800→2100 s. Parent result
и активный расчёт рецензия не меняет.

## 1. Какую проекцию реализует код

Пространство комплексных матриц рассматривается как вещественное евклидово
пространство с `⟨X,Y⟩=Re tr(X†Y)` и полной Frobenius norm. Положим
`H=(X+X†)/2`, `S=(X−X†)/2`. Для любого Hermitian `Z` имеем
`⟨H−Z,S⟩=0`, следовательно

`||X−Z||F² = ||H−Z||F² + ||X−X†||F²/4`.

Если `H=U diag(λ) U†`, унитарная инвариантность нормы и разложение по
диагональным/внедиагональным элементам дают единственный ближайший PSD элемент

`P+(H)=U diag(max(λ,0)) U†`.

Это также следует из ортогональности положительной и отрицательной спектральных
частей и вариационного условия проекции на замкнутый выпуклый конус. Поэтому

`dist(X,PSD)² = Σ(λj<0) λj² + ||X−X†||F²/4`.

Это комплексное Hermitian продолжение nearest-PSD конструкции
[Higham (1988), авторский препринт](https://www.maths.manchester.ac.uk/~higham/narep/narep126.pdf),
[DOI](https://doi.org/10.1016/0024-3795(88)90223-6); выше приведён необходимый
вывод, а не перенос real-symmetric формулировки без проверки.

`projection_worker!`, строки20–57, реализует именно эти операции. LAPACK `zheev`
возвращает собственные векторы столбцами: формула строки50
`U[a,j]*max(values[j],0)*conj(U[b,j])` имеет правильные complex фазы и индексы.
Полная сумма по `(a,b)` правильно учитывает обе комплексные off-diagonal
компоненты Frobenius norm. В `expected2` стоит ровно `neg2+anti2/4`.

При положительных скалярных весах `wE[e]*wk[k]` независимые blockwise проекции
также минимизируют взвешенное расстояние поля: вес не меняет minimizer блока.
При нулевом весе это seminorm, и minimizer такого блока не уникален; код всё
равно выбирает канонический unweighted `P+`. Ledger отдельно сохраняет обе
невязки, их шкалы и trace; размерность поля Γ и δ — энергия в solver units,
weighted norm использует dimensionless solver-grid `wE*wk`.

## 2. LAPACK, immutable ветка и измеряемая погрешность

ABI `zheev_` совпадает с используемым accepted Core eigenvalue helper,
за исключением нужного здесь `JOBZ='V'`; `UPLO='U'`. Workspace `3Nb` complex и
`5Nb` real не меньше требований `max(1,2Nb−1)` и `max(1,3Nb−2)`. Четыре
logical workers имеют собственные matrices, eigensystem и work buffers;
native worker ranges разделяют `(E,k)` блоки без пересечения записей Y.
Ненулевой LAPACK info и nonfinite eigenpairs отвергаются, не исправляются.

Ветвь `exact && low>=0` копирует raw block без реконструкции; если таковы все
блоки, возвращается сам исходный объект X. Это заявленная bitwise identity
для **computed** feasible exact-Hermitian поля. В прочих блоках порог равен
ровно нулю; положительный floor отсутствует. После реконструкции парная
Hermitian representation не вводит нового спектрального порога. X не меняется.

Ограничение существенно: `low>=0` от Float64 eigensolver не доказывает PSD
точной input matrix, а обратная реконструкция и явная symmetrization могут
создать крошечные отрицательные λ около нуля. Измеряемые eigen backward residual
`||HU−Udiagλ||F/||X||F`, `||U†U−I||F` и defect формулы расстояния — полезные
свидетели, но не interval certificate. `output_positive_rank_counts` означает
число положительных **computed λ до reconstruction**; это не доказанный ранг
Float64 output. Протокол и `arithmetic_scope` это ограничение сохраняют.

## 3. Проекция не меняет проверяемое уравнение на projected equation

Единственный runtime call site — `gamma_functions.jl:718`, только для
`info.kind=="AAII"`, до настоящего `evaluate(c,trial,...)`. Linear/rank reset
и единственный convex fallback не проецируются. После принятия `gamma=trial`
и AA history сохраняют испытанный Y, не исходный X.

В `_evaluate`, строки602–650, проверяется cone реального Y, затем неизменный
positive-source spectral operator, raw q и fresh μ. `fit_mu` по-прежнему требует
все raw `q[e]>=0`, конечность, `0<Ntarget<capacity`, положительную производную и
прежнюю точность `2e−13*Ntarget`. Нет q repair/μ normalization/floor.
Raw Φ вычисляется всеми прежними каналами, проходит старый cone guard и
**не проецируется**. `gamma_residuals` отдельно создаёт Herm working Φ для AA,
но trigger/best score используют исходное `norms(rawphi,gamma,...)`.
Старый raw unweighted+weighted `1e−9`×3 trigger сохранён.

Следовательно, потенциальный ложный stationary point `Γ=P+(Φ(Γ))≠Φ(Γ)`
не проходит raw gate. Скалярный контрпример `Φ(x)=x−1`, `x>=0`, имеет
projected stationary0, но raw residual−1; отдельный tiny check воспроизводит
эту границу. У истинного PSD raw root проекция не меняет допустимость root,
однако это не доказывает, что AA proposals найдут его.

`stable_spectral.jl` byte-identical eeca…: тот же hH/ΛH ledger и covariance A,
без Φ clipping. Public native materialization, исходный raw-h Context,
непроецируемые native maps и строгие consecutive gates в дельте не изменены;
`audit_native` по-прежнему строит fresh scattering/embedding, raw charge λ,
candidate residuals и per-channel component audit. Reduced trigger сам по
себе не объявляет native convergence.

Важно: прежний cone допускает арифметический отрицательный allowance; даже
после trial projection generic Float64 `GΓG†` не является формальным Gram/PSD
certificate. Поэтому отрицательный raw q всё ещё может честно отвергнуть
fresh Π. Проекция trial не должна скрывать такой исход, и код его не скрывает.
Entered domain failure расходует Π; input cone rejection до входа — нет.
Projection eigensolver failure находится до Π и, как прямо заявлено в
протоколе, останавливает опыт; существующий fallback относится к map-domain
failure, не к неисправной eigensystem или недостатку места для свидетельств.

## 4. Приёмка ограниченного опыта и остающиеся границы

Tiny tests содержат независимую complex 2×2 матрицу с λ={−2,3}, явную nearest
PSD matrix с правильной imaginary off-diagonal фазой, skew-Hermitian distance
term, idempotence в Float64, bitwise unchanged branch, zero eigenvalue без
floor и false-projected-fixed-point raw residual. Они уместны; их выполнение
**здесь not_measured**. Первые четыре полные raw-X/Y пары и дальнейшие selected
worst blocks позволяют проверить собственно trial correction независимо от
Φ; это не замена raw Φ snapshots или native audit.

При unit pass и внешних admission/guard проверках данная дельта математически
пригодна для одного заявленного matched pilot64. Нельзя заранее заключать,
что projection уменьшит raw residual, устранит q failures или даст сходимость:
ни contractivity QCL map, ни boundedness AA extrapolation+projection для неё
не доказаны. Неexpansiveness точной ортогональной проекции сама по себе этого
не обеспечивает. Новая identity честно задаёт алгоритмический вариант,
сохраняющий исходные raw/native критерии и не меняющий прежний исход cont64.
