# Streaming covariance helper: независимая static review

**Решение: approve helper8df7885c… для нового source-free tiny.** Физическая/whole-grid эквивалентность ещё не измерена; actual wrapper/source admission требует отдельной проверки. Рецензент Julia/Pkg/H5 и научные расчёты не запускал. Проверены все4 SHA256SUMS; helper8df7885cbf547ffd4077e611f6a88b933dc75c154053f29e420230de538db898 (118LOC), tiny488ba091cba029befea70a31848aced7d5e0af03d1d8587f058309681edf177d (69LOC). Малый stdlibPython-аудит ограничен128MiB, свидетельство `streaming-spectral-static-audit.json`.

## Математический объект и точные примитивы

Сопоставлено с accepted Core e171caa60cce676c6c4238dd671e7c527f2bd977: `production_operators.jl:_cavity_side!`235–292, `_cavity_production_plan`300–327, `_cavity_embedding_production`330–400, `_retarded_green_production`1266–1313; с qualified spectral eeca1f29… (`stable_spectral.jl`). V0/P4/sparse_plan/uniformPI guards и исходный hH=copy(c.cavity.hlocal) сохранены. Повторная сборка плана требует isequal(plan.hlocal,hH), поэтому центральный и все8 sidecell используют тот же точный hH. Проверяются original h correction и HermΛ correction relative+weighted≤1e−10. Λraw/ΛH для selected сохраняются **до** −iΓ/2, обратного извлечения слабой Γ/Λ нет.

Cavity recursion не переписана: `_cavity_side!` получает тот же localR и source≤=iΓ, source>=0, Tplus/Tplus†, remaps и direction. Все промежуточные источники/Schur matrices не симметризуются/не clip. Сохранены native Nb5 `_production_inverse_small!(…,Val(5))` и `_cavity_sandwich_small!` ветви; central inverse остаётся другим native `_production_inverse!` LAPACK путём. Native row-max-sum scaling D/ds→inverse/ds и cond(D) повторяют порядок old central. Физическая total source K=Γ−i(plusL+minusL), totalR=S+(plusR+minusR) сохраняют старые группировки. A=GR·K·GR† вычисляется теми же двумя generic BLAS mul!, не новым Gram/PSD алгоритмом. Generic floating covariance по-прежнему не строгий computed-PSD сертификат.

## Lazy read/write и параллелизм

`_cavity_local!` читает источник только через field[source,m,a,b], не требует AbstractArray/size/axes. LazyL возвращает именно im*Γ[source,m,a,b], LazyZero — ComplexF64 zero. Всё остальное берётся из настоящего S. `_cavity_side!` записывает только field[e,m,a,b] в конце каждой стороны; BlockSink получает именно соответствующий element Nb×Nb и не требует full destination array. Plus/minus sinks и native workspace принадлежат одному worker и переиспользуются последовательно. Minus/plus полностью перезаписывают свои6 матриц в каждой точке; старые side values не протекают.

Ranges по E-fast linear index не пересекаются; GR/A/κ/scale и samples имеют уникальные(e,m) записи. partials[:,worker] локальны worker. **Гонки q по k нет:** q вычисляется после завершённого sweep в отдельном serial postpass, с тем же qualified e→m→a order, gs/(2π) и we·wk. q не floor; negative count измеряется, caller должен применить прежнюю domain policy. Вычислять μ/Φ helper не пытается.

Source/cache guards проверяют grid object/kernel identity/shift contract; stale original cavity проверяется по source_grid/h/energy/weights/P4/sparse/drop0. PI FFT transform дополнительно проверяет размеры/columns/положительные веса/точную uniform geometry (`fft_parallel.jl:383–399`). Это scoped private helper, не универсальная замена `_validated_production_cache` и не разрешение использовать произвольный cache/physics. Wrapper обязан pin Core/GP/source и исходный Context; конкурирующие mutations c/Γ исключены protocol.

## Память и lifetime

Sweep владеет полными Γ/S/GR/A, то есть4F; coarse F=135,097,600B,4F=540,390,400B, κ/scale ещё5,403,904B. Λraw после selected copy освобождается перед GR/A allocation, HermΛ storage переиспользуется как S. В PI phase есть native verification output/transient FFT work; unreachable массив может остаться resident до GC. Поэтому3F preassembly и4F sweep — оценки живых full arrays, **не** RSS bounds. Explicit GC перед sweep не является общей гарантией RSS.

10selected fields≤128E×11k×25 дают≤5,632,000B данных; position O(NE), q O(NE), новый remap plan O(PNE), worker scratch O(Nb²). Original cp/caches/model/PI plans, wrapper reference results, histories/LO/Φ/GC/JIT остаются вне4F. O(Nb²) описывает reusable scratch; cond(D) может создавать bounded per-point LAPACK/SVD temporaries как старый native path. Полные9 side/total selfenergy fields здесь не создаются. Это оправдывает локальный memory mechanism; fine2x/fine4x или end-to-end memory/speed пока не квалифицированы.

## Tiny и следующие свидетельства

Tiny hash-binds old cf923 GP definitions и Coree171; сравнивает Nb2 и специализированный Nb5, комплексный noncommuting PD source при1e−24 и.03, весь GR/A/κ/scale/q и selected preassembly/localR/totalR/K. Strict V0 и out-of-range selected rejection покрыты. Reference использует старые полные native families/central inverse и generic covariance: это хороший arithmetic equivalence oracle; независимый физический oracle остаётся ранее выполненный selected Big, не этот tiny. Tiny workers=1; безопасность многопоточной записи здесь установлена по index ownership, а actual qualification должна запускаться при планируемом J4/BLAS1. Нулевая норма сравнивается без additive floor; обе нулевые дают0≤0.

Actual wrapper должен выполнить один old full и один stream на одной фиксированной Γlast64, без Φ/μ/SCBA; хранить all-entry scale-aware errors GR/A/q/κ/scale и selected split fields. Whole-source/code/manifest hashes до больших allocations, original Context/native gates и measured caps обязательны. Порядок per-worker scalar diagnostic reductions отличается от старого full norm; bitwise ledger identity не ожидается, а GR/A/q следует проверить отдельно. Блокирующих статических замечаний к helper не найдено.
