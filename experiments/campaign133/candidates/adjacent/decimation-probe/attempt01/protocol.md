**Frozen FD64 boundary-limit diagnostic v2 — до root green не запускать**

Задача: установить чувствительность frozen retarded спектра к числу
одинаковых клеток finite4/8/16/32 и независимому бесконечному surface limit.
Это разрешённый numerical boundary diagnostic. `∞` — отдельная identity,
а не изменение смысла прежнего baseline и не новая принятая SCBA.

Source — `raw/S01-warm-FD64/state.h5`, SHA256 e6ff…a6d (полный в input.toml),
исходная ещё несошедшаяся FD64 итерация, E=[−.3,.9]eV, NE30721,Nk11,Nb5,
70K, U=0, finite_chain4, source μ=0.008495393439818052eV. Пять LOCAL
retarded families суммируются ровно один раз. Архивные embedding и Green
arrays не участвуют в frozen Γ. Γ не clipping/repair/fit. Source read-only.

Шаги исполнения после независимого static review и явного root resource-green:

1. Tiny `preflight.jl` проверяет scalar analytic branch и complex 2×2 finite /
   Bloch oracles. Он не открывает HDF5 и не импортирует Core. Отдельная tiny
   define/bindings проверка `include(run_probe.jl)` проверяет загрузку функций,
   но не вызывает main (guard PROGRAM_FILE). Не выполнять Pkg/установки.
2. Main проверяет actual import Core HEAD e171 и clean status, hashes source,
   input, existing PI-control artifact+summary, saved source writer, fixed μ,
   source axes/scales/Hartree/basis, exact projected h and Γ SHA.
3. Восстанавливает local ΣR native `_ProductionPIFFTHilbertPlan` на принятом
   e171 с тем же ghost-hat uniform PI и source measure. Фиксирует Γ/R/hash,
   uniform geometry и arithmetic witness. Никаких microscopic kernels.
4. Строит свежую native finite_chain4 retarded embedding через штатный Core
   с dummy zero Keldysh buffers (они не физические состояния). Native Dyson
   без η даёт A. До sweep сравнивает per-k M/F, rsum и N/target с сохранённым
   `pi-frozen-control/attempt01` fresh_PI. При mismatch прекращает опыт;
   прошлые планы/history не заменяют это фактическое alignment.
5. Независимый потоковый sweep по E/k: finite4/8/16/32 одним проходом обычной
   рекурсии до L32; обе decimation поверхности одной факторизацией/2Nb RHS.
   Накопление только small per-k M/F, norms, witnesses, iteration histograms.
   Сохраняет independent finite4−native A weighted norm; mismatch — failure.
6. Actual occupied peaks k1/k6/k11: full45 finite oracle и independent Bloch
   θ refinement. Все selected matrices/witnesses сохраняются отдельно.

Соглашения/recurrence, единицы, нормы, stop и numerical thresholds подробно
в `algorithm.md` и `input.toml`. Менять их после просмотра результата нельзя
без новой probe identity. μ, Γ, h, T, Ntarget и quadrature weights фиксированы.
Смена boundary length — предмет опыта, не скрытая optimization.

Ресурсное предложение: один Julia4/BLAS1 процесс, MemoryMax6GiB,
MemorySwapMax0, CPUQuota400%, timeout600s с SIGINT и kill grace30s,
file-size/output32MiB; Julia heap hint4GiB. OS enforcement задаётся root
штатным scope/timeout снаружи; текст и soft checks сами cap не обеспечивают.
До тяжёлого запуска root проверяет общее окно с другими процессами.
Script soft time-stop565s на фазовых границах и каждые256 E/worker.
Внутренний numerical memory-budget Core оставлен исходным; он не заменяет
OS6GiB. Никакой новой platform/scheduler/cluster инфраструктуры.

Проспективные команды, сейчас **не выполнены**:

```sh
OPENBLAS_NUM_THREADS=1 julia --startup-file=no preflight.jl
OPENBLAS_NUM_THREADS=1 julia --startup-file=no --threads=4 --project=.build/campaign133/environment/accepted -e 'include(".build/campaign133/adjacent/decimation-probe/run_probe.jl")'
OPENBLAS_NUM_THREADS=1 /usr/bin/time -v timeout --signal=INT --kill-after=30s 600s julia --startup-file=no --check-bounds=yes --threads=4 --heap-size-hint=4G --project=.build/campaign133/environment/accepted .build/campaign133/adjacent/decimation-probe/run_probe.jl .build/campaign133/adjacent/decimation-probe/attempt01
```

Первая команда предполагает cwd decimation-probe; две другие — repo root.
Не создавать/модифицировать environment; `environment/accepted` уже существует,
но фактическую готовность и import HEAD всё равно проверяет guard.
Main принимает только новый output directory, сохраняет неизменные копии
кода/input/protocol/algorithm. `summary.toml` различает execution status и
surface/Bloch witness status. Failure сохраняет stage/error и завершённые k;
незавершённые M/F остаются NaN с completed_k=false, а не нулями.

Outputs: `alignment.h5`, `boundary_integrals.h5` (5×11 small M/F matrices и
per-k sensitivities), `selected_oracles.h5`, `summary.toml`, exact code/input
copies, внешние stdout/time/enforced-scope. Full Γ/ΣR/GR/A states не пишутся.
`boundary_integrals.h5` сохраняется также в finally, если sweep прерван.
Полный file output должен оставаться <=32MiB; внешний file-size cap нужен
для надёжной защиты даже при неожиданной ошибке.

Приёмка этого пакета: сначала preflight и source alignment; затем residual /
causality свидетели decimation и agreement с converged Bloch на selected
blocks. Сохраняются действительные результаты всех boundary variants.
Изменение M/F/r_sum с L — свидетельство boundary sensitivity на этой frozen
Σ; оно не доказывает interacting infinite fixed point, не объясняет raw
residual в прежней finite4 карте и не принимает energy/k/basis discretization,
модель рассеяния, Poisson или экспериментальную валидацию.
