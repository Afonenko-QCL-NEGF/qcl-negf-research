# Green-variable full SCBA: source-equivalence assessment

**FUTURE / UNRUN.** Только source/math assessment; никакого нового запуска,
реализации или изменения Core. По последнему указанию пользователя после
текущего FDseed128R4 — сохранить результаты и завершить. Описанный ниже
pilot является возможным будущим протоколом, а не разрешением его выполнять.

Проверен Core `e171caa60cce676c6c4238dd671e7c527f2bd977`, tree
`742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8`. Область: fixed U/h, fixed V,
all5/noEE, те же kernels, quadratures, shifts, full matrix options и native
`finite_chain`/4. **Да: три Green поля достаточны для root reparameterization.
Дополнительная старая Σ/embedding переменная не нужна.** Это не FD/Γ reduction,
не добавление e–e, не реализованный refactor issue32 и не доказательство сходимости.

## Что именно читает source

В `production_solver.jl:481–516` порядок однозначен:
`candidate_sc = _scattering_candidate_production(problem, green, cache, …)`,
затем `_embedding_self_energy_selected(…; scattering=candidate_sc,h,cavity_plan)`.
Старая `scattering` в candidate embedding не передаётся.

`production_operators.jl:1187–1268`: four static lesser/greater берут GL/GG;
`elastic_retarded=:direct_green` берёт GR для static R; иначе R строится из
свежих L/G выбранным native Hilbert operator. LO (`:974–1093`) читает GL/GG,
те же fixed Bose factors, κ/qK и native energy shifts. Model switches
применяются тем же `_apply_scattering_model!`. Для noEE никакая старая family
не читается. Cache содержит operators/geometry, не self-consistent history.

В finite-chain branch `:872–899` fresh locals суммируются native
`_sum_selfenergies`; `green` используется здесь только для shape. Cavity
`:235–393` стартует на каждой внешней границе с **нулевого embedding**, затем
детерминированно строит четыре side cells из этой fresh local R/L/G с remaps
E+jEp. Plan содержит копии h/energy/weights, remaps и drop; iteration outputs
не хранит и проверяет stale geometry. Сохраняются actual left `adjoint(Tplus)` /
right `Tplus`, side `hlocal=Herm(h)` и central original h. Нельзя незаметно
заменить это на иной all9-h policy. В других embedding branches используется
Green напрямую, но доказательство ниже заявлено для проверенной finite-chain ветви.

## Точная композиция и fixed points

Пусть X=(GR,GL,GG), F(X) — пять fresh local families; Y(X) — эти families
вместе со fresh finite-chain plus/minus. Их сумма embedding также materialized
native. Пусть D0(Y) — raw Dyson и два независимых products
`GR ΣL GA`, `GR ΣG GA`; D(Y) — тот же native paired-convex normalizer.
Тогда Green map **T=D∘Y**. Обе стадии обязаны использовать прежний source,
включая порядок суммирования Dict/enabled families и rounding policy.

Полные raw equations: `Y*=Y(X*)`, `X*=D0(Y*)`, `N(GL*)=Ntarget`.
Если raw domains допустимы, они дают T(X*)=X* и λraw=1. Обратно, T(X*)=X*
**вместе с λraw=1** означает `a=1,c=0`; paired map identity, поэтому
`X*=D0(Y(X*))`. Materialization Y(X*) даёт все полные raw equations.
Это алгебраическое соответствие exact roots; Float tolerances и физическая
приемлемость устанавливаются только полным native audit. T(X)=X без λraw=1
недостаточно: number correction способна возвращать normalized false root.

`_green_state_production(:1316)` всегда вычисляет A=i(GR−GR†). Поэтому A —
deterministic cache, не четвёртая независимая AA coordinate. **GG остаётся
независимым:** native paired branch передаёт explicitly computed greater;
на unresolved state B=−iGL+iGG может отличаться от A. Подмена
`GG=GL−iA` меняет проверяемый map. κD/scaleD также пересчитываются Dyson,
не служат independent state. FD/общий μ нигде не вводятся.

Native iteration сначала делает `X_i=D(Y_i)`, затем fresh `Y(X_i)`, residuals,
и лишь потом mixes self-energies. T(X_i) совпадает со **следующим** Dyson после
немешанного обновления αΣ=1; это не та же trajectory, что нынешнее αΣ=.2 или
native Σ-Anderson. Старт от Green unresolved snapshot допустим как новый warm
seed, но не называется unmodified continuation. На гладкой ветви derivatives
двух композиций имеют вид AB/BA и одинаковые nonzero eigenvalues
(`det(I−zAB)=det(I−zBA)`); поэтому устранение redundant coordinates не
устраняет физические slow modes автоматически. Paired map на Nraw=Ntarget
имеет branch kink: здесь безусловный Fréchet/spectral-radius вывод не заявлен.

## Trial cones не заменяют raw image

AA использует одни real affine coefficients для всех трёх полей; линейные
R/L/G identities histories тогда переносятся как те же affine combinations.
Нельзя независимо нормировать/миксовать три residuals и считать identity
сохранённой. Нормы history должны быть заранее определены для всех полей.

Возможный **только trial** safeguard: оставить Herm(GRtrial), вычислить
Herm(n)=Herm(−iGL), Herm(p)=Herm(iGG), применить nearest-PSD projection к n,p
без floor; затем задать `GRtrial=Herm(GRtrial)−i(n+p)/2`, GLtrial=in,
GGtrial=−ip. Optional native paired normalization при допустимой capacity
сохраняет B=n+p и global N. Это pointwise Pauli/causality/Keldysh working
proposal; оно не доказывает analytic Cauchy causality по всему energy axis.
Все representation/projection/number changes требуют ledger. Нельзя применить
эту reconstruction к **raw T image** или заменить ею source GG.

Утверждение о root correspondence относится к T, не к P∘T: даже scalar
`T(x)=−1`, `P(x)=max(x,0)` имеет projected fixed point0 с raw residual−1.
Поэтому сохраняется unprojected T(X)−X и raw λ. Physical exact roots с B=A
остаются identity points такого trial P; все floating native accepted states
с small nonzero B−A такими points автоматически не являются. Trial strategy
не гарантирует доступность каждого root или сходимость.

## Alignment, counters и обязательный выход

Сокращённая Σ residual `F(X)−F(X)` тождественно нулевая и не имеет статуса
native rΣ. Один T evaluation выполняет по одному fresh scattering build,
fresh cavity embedding и current-Green solve. Сохраняются input X, raw image,
raw λ/N/holes/capacity, paired a/c, norms по каждому Green field,
B−A/Hermitian/PSD/causality witnesses и counters **started/completed** по
трём стадиям, включая failures. На AA trial нет extra map, пока он не оценён.

На предполагаемом выходе materialize Y=Y(X), получить native Xnext=D(Y),
затем свежий Ynext=Y(Xnext). Сравнение Y с Ynext — настоящий fresh rΣ;
rD/rA используют aligned Y/Xnext, rK — unmixed Ynext, как
`production_iteration.jl:1–66`, `residuals.jl:539–668`. После этого нужны
три полных **unprojected, unmixed** native passes с первоначальными
per-family/embedding-total/±, raw products, λ, PSD/causality и остальными
acceptance gates. Green residual, projection step или exact N не принимают
решение о SCBA. Эти passes и отдельный final fresh audit считаются дополнительно;
если первый Dyson повторяет последний T output, повтор также честно counted.

```julia
# Псевдокод: только композиция существующих Core primitives.
function green_map(X, problem, h, cache, options, cavity_plan)
    G = native_green_state(X.GR, X.GL; greater=X.GG) # A from GR, no GG repair
    S, roundoff = native_scattering_candidate(G; return_roundoff=true)
    E, Eplus, Eminus = native_embedding_selected(G; scattering=S, h, cavity_plan)
    image, lambda_raw, total = native_current_green(h, S, E)
    persist_raw_before_any_trial_guard(image, lambda_raw, total, S, Eplus, Eminus)
    return image, lambda_raw                       # full graphs may now be released
end
# AA history stores only X and unprojected image-X; same coefficients for R/L/G.
# Any projection modifies the NEXT trial only. Native final audit materializes Y.
```

## Память и один будущий bounded control

F=`16*NE*Nk*Nb²` bytes. Для NE24001,Nk11,Nb5: F=105604400B.
Native `_ordered_iteration_families(:121)` хранит 5 locals+plus+minus:
**21F**, не24F; derived summed embedding добавляет3F transiently.
Depth5 **points**, state+residual: Green30F≈2.951GiB против native210F≈20.654GiB.
Если «AA5» означает пять secants/шесть points, это36F≈3.541GiB Green history.
Это только history: нужны current/trial/raw Green buffers, A, full fresh
families/embedding/total на один evaluation, kernels, FFT workspaces, residual
workspace и save/audit reserve. Σ graphs не исчезают из peak memory; их нельзя
удерживать во всех history records. RAM feasibility требует нового честного
estimate/OS cap, не обхода прежнего guard. GMRES5 также требует отдельных
Jv/basis buffers; он не предоставляется бесплатно этим выводом.

**Один tiny source oracle до любого будущего pilot:** Nb2/Nk1, V0,
finite_chain4, complex noncommuting T и positive raw local sources. При одном
и том же X два заведомо разных stale Σ/embedding records должны дать identical
fresh F/E/T; stale records не должны читаться. Для того же fixture полный
18×18 block-chain inverse с block-diagonal fresh local sources должен дать
central raw R/L/G как native cavity+central Dyson; затем сравнить native paired
output/λ. Это один fixture для hidden-state, orientation и alignment, не
запущенный тест и не QCL convergence proof.

Если когда-либо отдельно разрешено: один pinned 1mV/noEE fixedU checkpoint,
неизменные grid/h/Core, пустая history, максимум16 Green T calls с AA depth5
points и α=.2 fallback; первоначально guard-only либо явно новая trial-P
identity. Source hashes/initial independent GG сохраняются. Full raw gates
и extra final passes выше обязательны; no trigger → unresolved. Cost одного
T примерно содержит те же основные stages, что native iteration, но wall и
allocation gain не измерены. **В текущей кампании этот control не выполняется.**

Source SHA256: solver `ae776c29b611a84553116d9b25be8e04e2c16b7ff258c5ca5b11e03f555d8eaa`;
state `17dfe45e2cb128af0b7a27b22e7a9a7f9621f2784369c7dcd283f59cf07a29e1`;
operators `e28848511c295ef8741ca7f7baff1f1c8eb30566412e2212db4d5c4a1ed5c20b`;
mixing `5b608db96fb52b75c22b0d7f30da1cf743f872ae872b96ea9aa8124495d84e49`;
iteration `ff4662634f7db6cea8f8b578ab22ec1de10713aa522b6151887ca4cbc3276403`.
