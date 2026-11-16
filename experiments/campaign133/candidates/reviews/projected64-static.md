# Projected64: независимая статическая проверка

**Решение: approve для предусмотренного bounded опыта после нового Root tiny и отдельной math review.** Блокирующих ошибок в зафиксированной дельте не найдено. Это одобрение реализации алгоритмического варианта, не доказательство пользы проекции, сходимости или native acceptance. Рецензент Julia/Pkg/H5 и научные расчёты не запускал.

Проверены `gamma_functions.jl` SHA1c1eb3c3b5de030e21be5d23af0d74bf446781201e0cd12069a6f77d028027d4, `trial_projection.jl` SHA cf9dc3c14e2d640884cfd7a7fef49c6137b1ad7681901b15457388677d905e7e, `run_gamma.jl` SHA a1b49e4a36481da9c363dd449839a6b37d6854ac1eaf4b089b8774a67458936e; input b30102163f569b8051ac52b36c31e698ca9db2a00482c1bf306fc485d155420e. Все18 SHA256SUMS совпали. Стандартный Python-аудит ограничен128MiB; машинные свидетельства — `projected64-static-audit.json`.

## Объект сравнения и исходные уравнения

Input против cont64 меняет только identity, project_AA_trial false→true, pilot600→900s и overall1800→2100s. `stable_spectral.jl` eeca1f29… и прежние physical checks byte-identical. Тела `aa_trial`, `_evaluate`, `evaluate`, `native_warm`, `native_stage` byte-identical; поэтому сохраняются исходная Φ, same-A→q→freshμ→FD, все5 каналов, исходная raw Φ−Γ невязка/denominator/weighted criterion и trigger1e−9×3, исходный raw-h native exit. Projection не входит в Φ или residual.

`gamma_functions.jl:715–719` проектирует только AAII, после получения исходного предложения и **перед** evaluate. `:725` проверяет фактический Y; `:746–750` присваивает gamma=trial и добавляет именно принятое Y вместе с его свежим F в history. Linear/rank reset и единственный convex fallback не проектируются. Неудачный Y не заменяет последнюю принятую Γ, finally сохраняет last/best. Pre-entry cone не расходует Π; уже вошедшая map-domain failure расходует Π, как раньше. Проекция/eigensolver/output failure вне evaluate try останавливает опыт без fallback и без вымышленного Π.

## Матрицы, LAPACK, арифметика

`trial_projection.jl:20–26` правильно восстанавливает малую Nb×Nb матрицу из Field(E,k,a,b); линейные ranges разделяют E×k без пересечения записей Y. При all-finite X вычисляется H=(X+X†)/2; zheev('V','U') перезаписывает U eigenvectors по столбцам. ABI совпадает с accepted Core e171 `src/numerics/optimized/residuals.jl:_production_hermitian_eigvals!`, кроме jobV. Work3Nb≥2Nb−1 и rwork5Nb≥3Nb−2 достаточны. Accepted `_production_worker_ranges!` поддерживает do-аргумент и вызывает именно(worker,first,last); buffers принадлежат worker, не ThreadID.

Реконструкция Σ_j U[a,j]max(λ_j,0)conj(U[b,j]) правильна для комплексной Hermitian матрицы. Следующая симметризация явно сохраняет Hermitian representation. Положительного floor нет. При exact Hermitian и computed minλ≥0 блок копируется bit-for-bit; если changed=0, возвращается исходный X object. Это сохраняет такой feasible root/нулевой correction prefix без reconstruction. Математическое root preservation относится к Γ*=Φ(Γ*)≽0; finite zheev/reconstruction не даёт строгого PSD-сертификата для произвольного очень слабо положительного корня. Код и protocol явно сохраняют эту оговорку и проверяют реальный Y прежним guard/freshΠ.

Измеряются eigen backward residual, U†U−I, расстояние и его формула ||X−P||²=||X−H||²+Σ_{λ<0}λ², в том числе weighted norm с прежними we·wk. На нулевых weights это seminorm. Rank histogram обозначает положительные **computed** λ; это не точный ранг физического оператора. Независимый complex2×2 analytic fixture охватывает знак/сопряжение, Hermitian projection, weighted distance, unchanged-object, zero-floor, idempotence и ложный projected fixed point с ненулевой raw residual. Tiny не заменяет сравнение actual Y/history и common prefix по результатам.

## Происхождение, вывод и ресурсы

Restart строго импортирует тот же Γlast stableA16(call16), emptyAA/freshμ; Γbest/старый μ/Green/history не импортируются. Сохранены SHA parent input/identity/Gamma-last/map016/terminal/operator/spectral/model, frame verification и first whole rawΦ+9 scalar replay по256eps с отдельными bitwise flags. Дельта whitelist допускает ровно заявленную projection/timer identity, остальные physics/algorithm/grid поля строгие.

Первые4 AA proposals, включая unchanged, до опыта объявлены полными парами X/Y; следующие содержат ≤4 selected worst blocks плюс scalar ledger. Projection-attempt index независим от Πcounter, поэтому pre-entry rejection не затирает pair. Полные4 пары требуют1,080,780,800 logical bytes /1,207,958,400 padded dataset bytes, ещё небольшие HDF metadata; существующий admission учитывает именно padded chunk policy+8MiB на pair. Каждый файл проходит pre/post5GiB и общий8GiB admission; raw tmp при отказе сохраняется. Новый Y — один Field135,097,600B; малые4 worker buffers не масштабируются с NE/Nk. Фактическую RSS проверяет Root, не эта статическая оценка. J4/BLAS1,64TOTAL/pre180/pilot900/native900/overall2100 и прежние cgroup/file caps сохраняются; source gate сверяет внешний cgroup с input до model allocation. Input пока runtime_authorized=false.

## Неблокирующее замечание P2

`trial_projection.jl:9` на любой nonzero LAPACK INFO бросает Infeasible; INFO<0 означает ошибку API/programming, а не numerical domain. Здесь это **не** становится физическим fallback: projection находится вне evaluate/fallback try, ошибка всё равно fatal до Π. Поэтому минимальный frozen опыт допустим без churn. Если такой отказ случится, классифицировать INFO<0 как API error и не включать его в число физических domain rejections; будущая поддерживаемая реализация должна разделить INFO<0 и INFO>0. Author сам сообщил эту оговорку до review.

Нужное actual свидетельство: root tiny + measured projection defects; первые4 пары/общий prefix относительно cont64; затем raw residual/число entered failedΠ/время/RSS. Стационарность projected шага сама не принимается. Улучшение относительно старого AA10 нельзя приписывать одной covariance/projection; ближайший matched control — stableA-cont64 AA4096 с тем же seed.
