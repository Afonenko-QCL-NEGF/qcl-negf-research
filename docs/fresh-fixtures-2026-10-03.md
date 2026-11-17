# Малые свежие проверки и следующий научный шаг

Постановка 2026-10-03. Старые матрицы не являются входом. Пользователь разрешил
подготовку реализации и локальные проверки, затем увеличил лимиты времени/RAM.
Новая большая SCBA/Poisson кампания этим планом не запускается.

## Зафиксированные операторные fixtures

Исходники: `QCLNEGF.jl/test/numerics/frozen_scientific_fixtures.jl`,
`test/physics/full_window_basis_preserves_the_spectral_subspace.jl`,
`test/numerics/analytic_product_integration.jl`.
Julia 1.13.0, существующий Manifest, CPU=1, BLAS=1. Конечный общий лимит
операторных запусков с JIT и исправлениями: 30 минут, RAM 3 GiB,
вывод 10 MiB. Численные массивы в Git не сохраняются.

| Вопрос | Inputs и наблюдаемая | Критерий и граница вывода |
|---|---|---|
| Меняет ли дробный energy shift связность? | 17 E nodes [-8,8], trapezoidal widths, положительный фиксированный scalar spectrum, LO shift=2, field shifts=0/0.5; actual remap/contraction support и Sigma lesser/greater | 2/1 connected sectors; weighted integral of Sigma lesser × G greater − Sigma greater × G lesser в каждом секторе <=64 eps times positive weighted collision scale. Это операторная fixture, не уникальность нелинейного решения и не полная QCL transition matrix |
| Точна ли occupation-weighted impurity contraction? | 70 K, Nz=18,Nb=1,P=1,NE=17,Nk=2; geometry/screening/basis неизменны; adaptive tolerance=1e-5; literal positive midpoint oracle 2048/4096 angles | Oracle change <1e-5 и relative contraction difference <5e-5, нормировка norm(oracle4096)>0. Не вывод о точности integrated current |
| Что теряется при basis cropping? | Nz=18,Nb=3,P=1,NE=31,Nk=3; полный 3-period projected reference и spatial operator одной represented geometry, overlaps/tails | Existing tolerances 2e-12 и subspace residual<1e-12. Не меняем seam/clipping, не сравниваем разные spatial grids |
| Проверен ли finite-window PV независимо? | Analytic triangular Gamma на 3 E nodes; Lorentzian 401/801 nodes, window [-20,20] | Triangle error <=2e-15; existing interior error <5e-4 и refinement improvement. Finite-window tails остаются ограничением; FFT/direct agreement отдельно не oracle |

Состояния инфраструктурных fixtures сохраняются новым Runner contract; они
синтетические и не получают scientific acceptance. Операторные checks выше
дают отдельное code/model evidence, без экспериментальной валидации.

## Следующий шаг после анализа fixtures

Сначала проверить один новый stationary state с достоверным архивом и recovery
на постоянном узле. Требуется отдельный frozen scientific plan до запуска:
цель, точный model/numerical/algorithm fingerprint, seed, конечный cumulative
SCBA/outer budget, admission и next-step rule. Сохранить paired_convex,
raw convergence gates и фактические algorithms; nonconvergence сохраняется
как результат и не вызывает scientific retry.

Gates последовательны: integrity/provenance → iterative convergence → physical
checks → discretization evidence → experimental comparison. Отсутствующие
измерения обозначаются not_measured. Только после просмотра этой одной точки
пользователь выбирает уточнение сетки либо отдельное изменение физической
модели. 70 K нельзя напрямую валидировать экспериментом при 200 K.
