**Независимое static review positive-source operator qualifier**

Решение: **одобряется один bounded operator-only probe по root admission**
на замороженном `positive-source-preflight/preflight.jl`
`a34b1da110293ba289a78a437354436970ed8362a184779eb99e1f87690ade18`,
input `45b2857d2fc4cbcd01a8765a4745e513b9539c87c481aac044318a95d0cf3fbe`.
Блокирующей math/sign/shape/API ошибки не найдено. Есть P2 caveat reporting
matched sub-control ниже; root явно принял его без code-change для этого
finite probe и не будет определять фактическое исполнение только по label.
Нового tiny churn или обязательной правки до actual run review не требует.

Рассмотрены frozen qualifier, input/protocol/handoff и нужные
функции helper `cf923bfb7e5dc56d5aa296c6863e8615a92984f82e125433ea1002bac9c9ee63`.
Scope основан на `adjacent/positive-source-operator-plan.md` и original
`stable-spectral-plan.md`: здесь первая минимальная covariance-ветвь,
не предложенная более широкая layer-by-layer Gram/reclosed Schur схема.
Accepted Core e171caa60cce676c6c4238dd671e7c527f2bd977 не менялся.
Reviewer Julia/Pkg/HDF5/массивы/научные вычисления не запускал; код автора
не редактировал. Проверены малые SHA и первичные функции исходников.

**Сохранённые контракты**

1. Source guards на105–146 связывают original warm strictly-positive-
   measured seed SHA4142e06d…, исходный input и model SHA6b41…,
   actual Core HEAD/tree/clean status/Julia, schema, все axes/scales/spin.
   Это именно original warm seed, не поздний Γ checkpoint с tolerated
   tiny negatives. `seedcone.minimum_eigenvalue>0` — floating eigenvalue
   witness, не strict interval/Gram certificate; source не проектируется.
   Source/model reads выполняются только после root runtime-authorized
   input и explicit resource guards. Input-runtime должен остаться в той
   же папке (relative source paths) и отличаться только authorization.
2. `hH=copy(c0.cavity.hlocal)` на148 — **точный** native Float64
   Herm(hraw), не independently reconstructed matrix. Новая cavity plan
   построена именно с hH, поэтому guard `plan.source_h==passed h`
   согласован. Центральный Dyson также получает c.h=hH. Side-repeat defect,
   Δh norm/scale и raw/hH сохраняются. Формально all9 cells имеют одну
   Hermitian h-политику; raw-central старого native baseline здесь меняется
   и не объявляется unchanged map. Tplus и Tplus† остаются штатными.
3. На44–48 native uniform PI Λraw переводится в explicit Hermitian copy,
   с прежним1e−10 arithmetic correction target и absolute/weighted ledger.
   `LH` не мутируется при assembly: `S=copy(LH)` перед −iΓ/2. Значит
   selected `Lambda_H_preassembly` действительно сохраняет split input,
   не восстановлен из уже округлённой ΣR. Это позволяет строить intended-
   model Big oracle без утраченных при assembly off-diagonal tiny Γ.
4. Native all-component cavity получает `(ΣR=S,Σ<=iΓ,Σ>=0)` на49–50.
   Это **spectral f=1 test source**, не bath occupation или physical density.
   При соглашении Σ<=iA physical embedding broadening есть −iΣemb<,
   поэтому `physical=Γ−iΣemb<` на51 имеет правильный sign. Σ>=0 не
   подмешивает отрицательный/физический hole source. Native retarded
   recursion не переписана, extra η/clipping/floor не добавлены.
5. Центральный GR от totalR и того же hH; новый A на10–21 — обычная
   `GR physical GR†` через два `mul!` в собственных worker buffers.
   Worker do-block API существует в Core, включая function-first adapter.
   Asub на53 строится из **этого же GR**, без второго inverse. Effective
   source на56 равен `i[(h+ΣR)−(h+ΣR)†]=i(D†−D)` при real E;
   E cancels. Physical/effective и их исходные компоненты сохраняются.
6. На23–34 q вычисляется из того же A_cov:
   `qE=gs wE/(2π) Σk wk ReTr A_cov`. `_seed_number_weights_production(GR)`
   не вызывается. Затем fit_mu→fd_green(GR,A_cov)→phi(c,green) на166–171
   используют один A; static4 потребляет A, LO — FD lesser/greater.
   Нет q correction, density rescaling, AA, trial/Φ PSD projection.
   `same_A_number_check` действительно использует существующий native
   number-functional API и исходную positive quadrature measure.
7. Matched subtraction control на173–180 имеет те же hH/Γ/GR и A_sub;
   свой fresh μ допустим, поскольку его q изменён вместе с A. Отрицательный
   qsub останавливает fit без clipping. Это различающий контроль A
   arithmetic при fixed h-policy; старый central-raw-h trajectory не может
   быть приписан одному изменению A.

**Доказательство и его floating предел**

В exact arithmetic для Hermitian h/Λ, Γ≽0 и nonsingular finite-chain
inverse: `i(G−G†)=G Γtot G†`. Индукция native cavity при Σ<=iΓ даёт
`Σemb<=i T A_surface T†`, поэтому central physical source
`Γ−iΣemb<` воспроизводит full9-cell source covariance. Это оправдывает
qualifier как отдельную arithmetic identity той же exact underlying model.

В Double `C g C†`, assembly local Σ и inverse могут создать ΓR/ΓK mismatch.
Generic `G K G†` не гарантирует строго положительную materialized matrix
в FP и не доказывает accuracy GR. Здесь нет source factor/Gram diagonal
sum(abs2) proof: протокол и названия это честно оговаривают. Cone/condition/
same-G differences сохраняются как witnesses. Если input/covariance/q/rawΦ
не проходит прежний gate, опыт должен остановиться/показать unresolved,
а не разрешать floor или скрытую смену acceptance. Selected h/ΛH/Γ/T/Σ/GR
достаточны для последующего promoted split-input full45 oracle; сам этот
Big actual-oracle в qualifier ещё **не выполнен**, full scientific operator
accuracy не объявляется принятой одним положительным q.

**Tiny controls и actual evidence**

Fixtures78–103 содержат scalar Γ1e−27 analytic response и quadrature q,
complex Hermitian2×2 H с noncommuting PD K и обычной spectral identity,
а также weak `1e−27 K` с independent Big192 inverse/explicit four-index
covariance accumulation. Это существенный matrix-phase control: H,K не
коммутируют, не только scalar/diagonal пример. Big192 здесь проверяет
малый manufactured operator, не actual S01/cavity/global seed.

Root прислал и reviewer прочитал `positive-source-preflight/tiny-root.log`:
8/8 pass, exit0, wall6.14s, peak526984KiB, Swap0, J4, BLAS1. Это root
runtime evidence именно tiny primitives; reviewer его не повторял.
Оно не принимает actual native spectral-cavity consumer или full map.
SHA tiny-root.log:62d9624642bfa27ceb53c41c3908574f785d169bf06546bd0a536363f2160d2d.

**Native comparability, сохранение и замечание P2**

`native_compare` на61–76 вызывает native пятисемейный candidate на
**переданном stable-A Green**. Per-channel kinetic difference/KMS и
retarded sum vs HermΠ[rawΦ] сохраняются. Поле
`kinetic_KMS_control_pass` относится только к своему declared1e−10
kinetic/KMS maximum, не к weighted/retarded/full native acceptance.
Native fresh GR и unconstrained fixed point не вычисляются: новый variant
не должен приниматься как старый raw-h native baseline по одному consumer
pass. Это различие в returned scope явно имеется.

P2: catch173–180 охватывает не только fit_mu, но и вычисление Φsub и cone.
Если Φsub уже вычислена, а её cone бросает Infeasible, label
`matched_sub_control.not_measured=true` слишком широк. Primary truth
сохраняют maps_started/maps_completed, payload `Phi_sub_raw`, μ_sub и
matched norms, добавленные до cone. Root принял caveat для данного finite
probe: result reviewer обязан восстановить stage по ним и не называть
выполненную карту невычисленной только из label. Для последующей reusable
реализации catch следует сузить/разделить статусы, но это не требуется
как новый цикл правок перед уже выбранным operator-only измерением.

Конечный budget по control flow: один spectral-cavity response, максимум
две raw reduced Φ и один native candidate; это не итерационный solver.
Started/completed counters стоят до/после actual map calls. Checkwall
на границах не заменяет внешний300s cap. MemoryMax12.125GiB/High11.5GiB,
Swap0, CPU≤400%, LimitFSIZE≤32MiB проверяются по cgroup/proc; root должен
задать их до запуска. Full fields находятся в памяти, output выбирает
несколько E/all-k matrices, whole q vectors и small h/T/axes. Σ/GR/Γ
full arrays не пишутся; selected data достаточны для различения механизмов.
Пиковый actual RSS/cost до полного run всё ещё not_measured.

Итог: frozen qualifier и root tiny evidence достаточны для одного
ограниченного различающего опыта с этим caveat. Ни прежний STOP, ни raw-h
native acceptance, ни SCBA/Poisson/grid/экспериментальная валидация этим
review не пересматриваются.
