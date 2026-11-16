# C: FD64 same-input static retarded operator comparison — WRITE ONLY

Статус: скрипт написан и прочитан статически, **не запускался**. Julia/Pkg,
подготовка новой исполняемой среды и расчёт не выполнялись. Baseline128
остаётся единственным активным научным процессом. Следующий запуск требует
отдельного решения ведущего после этого design/static checkpoint.

Вопрос: как на одном фактическом несошедшем FD64 input изменяются четыре
static retarded candidates при переключении native
`elastic_retarded=:hilbert` → `:direct_green`? Оба native candidate calls
получают один и тот же объект GreenState с исходными GR/G</>/A, без FD
repair, density normalization, mixing, solve_scba, Poisson, cavity/embedding
refresh и central Dyson. LO route и lesser/greater channels не меняются.
Dynamic e−e отсутствует в этом all5 input; script проверяет exact channel set.

Источники фиксированы в input.toml: CoreC
`2dc5a26c1a7a6260febe3ab07bf6439fd2f0bf19` в отдельном чистом worktree
`.build/campaign133/core-direct-static-retarded`, FD64 stateSHA
`e6fffb997f36397f877aee1e2c98b1dbad7d97417d45458a042b5f074ce57a6d`.
Source input/identity остаются исходными VM файлами. SHA большого state
будет проверен только внутри разрешённого будущего процесса; сейчас HDF5
не открывался. Source input определяет zero-bias70K, all5, Nz96/Nb5/Nk11,
NE30721, E=[−.3,.9]eV, basisP3, embedding finite_chainP4 и U=0.

Пути actual implementation на CoreC:

- `production_operators.jl::_scattering_candidate_production`: для обоих
  options выполняется actual native selector целиком, включая LO. Ни
  одна candidate family не заменяется ручной формулой.
- `_selected_static_contraction`: те же Khat/qK/wk для остальных terms.
- `production_fft_hilbert_transform`: тот же primary FFT plan и wE, что
  при старой finite nodal PV reconstruction; no wrap, omitted diagonal
  и finite domain остаются исходными.
- `_validated_production_cache`: source/shift/Hilbert/contraction/FFT
  contracts проверяются для обоих options. Elastic selector не меняет
  kernel или FFT representation. Script проверяет object identity
  kernel dictionary и FFT plan, поэтому дополнительный cache не строится.

Различающие свидетельства — per-channel full complex ΔR, Hermitian и
anti-Hermitian parts, старые/новые R scales; ΓR_new−ΓK_new и остаток после
вычитания K[A−B]; отдельно старая/новая SHA всех R/L/G. Все независимо
вычисленные L/G должны совпасть побайтово, а LO — во всех трёх компонентах.
Source Green hashes проверяются перед/после опыта; неизменность не
предполагается только из API.

Полная discrete identity сохраняет фактический input defect:

    n=−iGL, p=iGG, B=n+p, D=A−B,
    X=Herm(GR)−H_W^q[A],
    ΔR=K[X]+H_W^q[K[D]]−(i/2)K[D].

Script вычисляет каждый из трёх terms, residual identity и
ΓR_new−ΓK_new−K[D]. D не зануляется. Поэтому anti-Hermitian Δ нельзя
автоматически приписать tail или багу selector. Если D≈0, reduction к
Hermitian Δ проверяется измерением; она не навязывается input repair.
Научная опора — section3 static-retarded-assessment.md и complex-linear,
adjoint-preserving actual K. Использование production K/H с обеих сторон
проверяет consistency реального компонента, но само по себе не является
независимым oracle его реализации. Для независимого вывода нужны CoreC
rational/covariance fixtures и их отдельный review.

Каждая норма сохраняет unweighted Frobenius и quadrature-weighted L2
по полной/trusted области. Последняя использует
g_s/(2π) wE_e wk_m в квадрате component amplitudes, domain trusted
[−.22,.82]eV. Arrays и residuals dimensionless; это не интеграл signed
collision и не full SCBA r_K/r_Σ. Absolute norm и scale сохраняются
раздельно; relative к нулевому scale возвращает Inf при ненулевой
разнице. Тут нет новых scientific acceptance tolerances или подгонки μ.

Память и исполнение будущего опыта:

- Один процесс≤300s, OS MemoryMax6GiB/MemorySwapMax0/CPUQuota100%,
  Julia4/BLAS1, heap hint около4GiB, stdout/files≤32MiB. Это prospective
  limits, сейчас они не объявляются реально enforced.
- Native full-solver estimate и guard остаются честными с budget
  12992276070bytes≈12.1GiB. Это оценка более широкой solver trajectory,
  а не новый OS cap и не обещание, что весь solver помещается в6GiB.
- Native old candidates считаются целиком; затем сохраняются только их
  R arrays (5×135.2MB≈676MB) и hashes L/G, old family dictionary
  освобождается с GC до нового native candidate call. New all5 families
  занимают≈2.03GB, source Green≈.54GB, два shared audit fields≈.27GB.
  Per-channel audit временно использует девять135MB fields≈1.22GB;
  они освобождаются перед следующим channel. Без JIT/GC/FFT scratch
  это≈4.74GB. Поэтому6GiB cap — ограниченный опыт с риском OOM, не
  доказанный memory acceptance. Никаких full candidate arrays на диск.
- Результаты пишутся поканально: selector-hashes.toml, channel-*.toml,
  компактные worst-delta/worst-identity5×5 witness blocks, identity и
  summary. При timeout/OOM частичные файлы сохраняются как partial;
  отсутствие summary не означает завершения.

До запуска ведущему необходимо подготовить native Pkg среду, которая
действительно импортирует отдельный CoreC worktree, сохранив её
Project/Manifest и hashes. Baseline environment и live Core/Runner не
менять. Script до чтения state проверяет импортируемые Core HEAD и
clean status. Затем штатный scope+timeout должен быть выбран с новым
output directory. Имеющееся численное guard не обходить и не переписывать
под6GiB. Installation/framework wrapper не требуется.

Optional frozen-embedding retarded response не реализован: первичный
операторный Δ достаточно различает вопрос без дополнительных Green
arrays, sourceLO/embedding copying и трактовки occupation. Даже успешная
identity/LO invariant не установит SCBA convergence, физический tail
или преимущество новой closure в полном mixed LO+static nonlinear solve.
