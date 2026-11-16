# Write-only handoff; НЕ разрешение на запуск

Артефакты: `analytic_jv.jl` (finite module, no I/O/model/automatic execution),
`tiny_checks.jl` (explicit callable small tests), `input.toml`, `proof.md`.
Это bounded prototype+API, **не готовый actual runner**. Чтобы не вводить ещё
один general framework, существующий pinned source/model/Gamma loader и raw
output helpers должен вызвать короткий owning adapter после отдельной review.
Автор не запускал Julia/Pkg/H5, не читал большие arrays и не менял Core.

1. Source gate до actual allocation: imported Core HEAD/tree/clean tracked,
   Julia version/J4BLAS1; SHA всех транзитивно included parent sources по
   parent identity, helper8df и этого prototype. Parent directory/model/input/
   identity/Γlast64/terminal/map64 hashes зафиксированы в input из root JSON.
   Γ читается **без** Herm repair/rescale/interpolation; own model rebuild
   обязан пройти existing `GP.verify_problem_model` against model700a….
2. Сначала отдельный tiny: include pinned GP, this module and tiny definitions,
   затем `tiny_checks(GP)`. Это small synthetic full45 и small native kernels;
   не actual S01 map. Time/memory root назначает отдельно. Любой fail → noactual.
3. Сохранить actual source-PD witnesses и статус pending exact proof. Старый
   cone guard и sampled minima не объявлять global exact-PD proof. Для v=Γ
   PSD safety exact следует из scaling **при** PSD base; strictness не нужна
   этому направлению. Без interior proof full physical Fréchet/Arnoldi
   interpretation остаётся условной: smooth algebraic response при finite
   inverses/χ не становится доказательством допустимости finite arbitrary
   direction. Это ограничение не лечится floor. При установленной true source
   negativity qualification physical causal base отвергается.
4. Actual OS envelope: MemoryMax13019119616B=12.125GiB,
   High12348030976B=11.5GiB, Swap0, J4/BLAS1, CPU400%, timeout600s, heap8G.
   Native cache numerical guard не уменьшать/не обходить. Время включает
   context/JIT/hash/model comparison. ≤5 enteredΦ, ≤4 enteredJv,
   ≤3GiB output, ≤2GiB file. Wrapper учитывает wall от process entry, а не
   только от создания `Limits`; перед каждым вызовом и записью проверяет cap.
5. `base=fresh_phi(GP,SS,c,Gamma,limits;keep_base=true)`; primal base должен
   повторить saved rawΦ64/μ/nine scalars within existing256eps replay scope.
   SS reduction diagnostics могут иметь иной порядок сумм; это qualification
   ledger, не замена whole-field rawΦ replay. Сохранить base rawΦ до probes.
6. Четыре freshΦ для `(1±ε)Gamma`, ε=1e−3,5e−4: rawΦ и μ/channel/domain
   summaries сохранять по мере получения, затем Dε и их comparison. Готовый
   `central_control` возвращает arrays только при успехе всех вызовов: actual
   adapter предпочтительно вызывает `fresh_phi` по одному и сохраняет каждый,
   чтобы поздний fail не потерял предыдущие raw witnesses. Это та же формула.
7. Один `jvp(...,Gamma,limits)` сравнить с Dε. До Arnoldi reviewer/root принимает
   derivative consistency по absolute+relative+weighted signal/noise/curvature,
   не по одному guessed threshold. `primal_GR/A replay`, numberδN и raw/Herm
   Jv conversion должны быть finite и <=1e−10 inherited arithmetic target;
   абсолютные числители/шкалы сохранять, zero/zero=0, nonzero/zero=Inf.
8. При квалифицированном Jv вызвать `arnoldi3(...;qualify=record->...)`.
   Сохранить Hbar, Ritz real/imag/residuals, restricted-B SVD, seed conversion,
   each Jv ledger и Q directions (небольшие m≤4 full fields), но не Γ update.
   `qualify` проверяет перечисленные numeric gates; default always-true не
   использовать. Exact breakdown честно останавливает span; near breakdown
   требует noise qualification, не auto-normalization of noise.
9. На любом cone/q/χ/inverse/replay/qualification/wall/memory/output failure:
   stop и save started/completed counters + исходную причину; programming
   errors не переводить в physical mode или fallback. Никаких дополнительных
   Φ/Jv, SCBA/native/materialization calls, adaptive ε или Newton retries.

Limits учитывает entered/completed Φ/Jv, но не OS enforcement и не build cost.
Один full complex field135097600B; baseΓ/R/GR/A/Φ≈0.63GiB, Q≤4≈0.50GiB,
temporary signed Green/scatter/PI fields и cache/JIT дополнительно. Весь
process peak не измерен; fit12.125GiB **не обещается**. Нет full45×NE storage
и eigendecomposition глобального Jacobian. Полезный исход — qualified local
response/near-null witness либо конкретный отказ arithmetic/domain, не SCBA pass.

Минимальный owning adapter (псевдокод, существующие IO helpers переиспользовать):
```julia
source_gate_and_archive(); c = GP.build_context(parent_input)
Gamma = pinned_Gamma64_and_verify_model(c); budget = Limits(seconds=remaining())
base = fresh_phi(GP,SS,c,Gamma,budget; keep_base=true); save_raw(base.raw)
require_parent_replay(base)
for eps in (1e-3,5e-4), sign in (+1,-1)
    raw,meta = fresh_phi(GP,SS,c,(1+sign*eps).*Gamma,budget)
    save_raw_and_meta(raw,meta)          # preserve before next actual call
end
d,meta = jvp(GP,c,base,Gamma,budget); save_raw_and_meta(d,meta)
require_derivative_validation(d,saved_central_differences)
k = arnoldi3(GP,c,base,budget; qualify=arithmetic_gates)
save_small_krylov_and_Q(k); save_counters(budget)
stop_without_Gamma_update()             # scientific_acceptance=false
```
