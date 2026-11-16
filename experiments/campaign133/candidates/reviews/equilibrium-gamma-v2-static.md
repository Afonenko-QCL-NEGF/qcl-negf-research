# Независимый review Γ pilot16 V2: дельта и small runtime

**V1 blockers устранены в проверенном V2 freeze; нового launch blocker не найдено.** Допустим ровно предложенный bounded pilot16 после resource gate с mandatory initial actual consumer/KMS control. Actual QCL Π/native convergence остаются `not_measured`:75 small checks их не заменяют. Рецензент Julia/Pkg/compute не запускал, source не редактировал.

## Freeze и свидетельство

`run_gamma.jl` SHA256 `abda3fb908137eba303ac92aac35dcbf6e312a0dc6edb0081c39be72d7c0b46e`; `gamma_functions.jl` `4475101718dab32cb5f4c4446390a496caaf2731366ca6e8e503c76232970f46`; `checks.jl` `0942e463bc5ff5e26e721d17fb6dc9e72abec6f416d3b3187b029857f7c0f494`; input `5dadb59f19fddbf0aa483687b7addb8a792d0e97681456022585b47796be3498`; patch `01d122e50c02bde91029223c6c807a86b1dd02e13cf1797174048157aae843a3`. `sha256sum -c SHA256SUMS`: все8 entries OK. V1, его parse failure и прежние reviews сохранены `unexecuted-v1/`; не переписывают отрицательный исход.

Input diff с V1 — только identity suffix v2; physics/grid/caps прежние, imported Core e171/tree742. Main дополнительно фиксирует native cap32. Root `environment/accepted/gamma-small-checks-v2.log`:75/75 pass,7.2s testset,14.03s process wall, peak646448KiB, exit0. Command проверяет actual Core path, затем include+small_checks без S01/model build. Marker `GammaSmallV2ChecksOK_NoActualS01Maps` согласован с кодом. Тест непосредственно вызывает helpers при runtime_authorized=false; это намеренно не полный entrypoint/source/resource gate.

## Проверка findings V1

1. Parse исправлен явным `0.2 .* x`; inclusion и tests прошли. Это подтверждённое устранение прежнего execution blocker.

2. `gamma_residuals:170` раздельно вычисляет raw `norms(rawphi,gamma)`, working `norms(Herm(rawphi),gamma)` и representation correction. `evaluate:583` использует **raw** residual для streak/trigger, публикуемого rΓ и best score; возвращаемый AA residual остаётся Hermitian working. Tiny test допускает asymmetric field по native cone allowance, показывает raw relative вышеtrigger и working0. Норма denominatorrawΦ/нулевая политика сохранены, raw output не изменяется.

3. `known_inverse_failure:18` принимает только SingularException с неотрицательным info и LAPACKException с положительным info, плюс рекурсивные TaskFailed/Composite/Captured wrappers, если все причины известные численные. Negative LAPACK argument, MethodError, BoundsError, произвольный/budget ErrorException отвергаются. `retarded_state` обёрнут `inverse_domain`; linear/fallback не начинают бесконечный recovery. Tests включают реально failed task, mixed numerical+programming Composite и propagation programming error;75pass подтверждают доступность Base task exception API. Predicate не является универсальным catcher всех возможных numerical errors: неизвестная причина должна остановить attempt, это честный отказ, не hidden fallback.

4. `write_metadata:389` материализует BitArray в dense Int8, Symbol в String; model/stencil loops используют его. Actual HDF tiny roundtrip trusted/valid masks и eltypeInt8 прошёл. Complex hyperslab2/3/4/6D checks также прошли. Whole model consumer ещё не измерен, но найденный packed-pointer путь устранён; schema/mask storage явно названы. Actual AlgorithmOptions архивируются, Core API не менялся.

5. Protocol и algorithm теперь явно объявляют **nonmonotone, feasibility-safeguarded AA**, без отсутствовавшего merit/regression criterion. Это разрешает расхождение promise/implementation без выбора thresholds по trajectory. Fixed16/onefallback/QR/L1/PSD guards прежние. Retained firstpair/Rdiag/initial scales добавлены как small evidence; formula неизменна.

## Native diagnostics и accounting

Per-channel KMS сохраняет absolute/scale, weighted wE*wk counterparts, max block Frobenius error/locator и max block relative. Initial gate по прежнему относится к actual **seed spectrum/current μ**, не uniform всем будущим Γ. Scope прямо назван. Fresh-candidate embedding в precheck больше не обещана; same-input dummy/physical-FD cavity comparison сохранён, full fresh embedding остаётся в native materialization/final audit.

Final audit теперь включает пять local families и embedding total/plus/minus: stored/fresh collision, ΓR−ΓK identity, streaming separate R/L/G residuals с fresh denominator и weighted counterparts. Это полезная расшифровка, отдельно от unchanged owning native max rΣ. `embedding` здесь — сумма двух boundary contributions, не total local+embedding; общий total используется в owning native rK/λ/rD/rA, отдельного total identity metric не появляется. Формулировка actual report не выдаёт это за большее coverage.

FDR остаётся diagnostic-only common μ, arrays не переписываются. V2 правильно разделяет unweighted joint FDpair denominator, weighted wE*wk, occupied lesser error/denominator, trusted weighted pair и pair/A. Raw pair строится из **stored total Σ</>**, normalized — published native Green; fresh candidate consistency отдельно owning rK/candidateλ. Это naming предотвращает прежнее смешение occupied FDR с pair/A. Measures/solver units и отсутствие spin/2π в Frobenius norms названы. Collisions остаются signed scaled moments с явно указанной мерой без SI power conversion.

Materialization/audit start отмечены до вызова, completion/report отдельно; observer назван count завершённых public maps. `native_accounting` даёт conservative started lower/upper и completed lower bound; при interrupted solver не выдумывает exact started count. Public max_scba30 фиксирован, reserve1+30+1 и main cap32 сохранены. Upper bound может оставаться31 даже после раннего нормального выхода solver; это консервативная граница, не exact telemetry. Unit bounds проверены. Budget wall check между maps не заменяет external timeout/cgroup; версия не меняет эту ограниченность.

## Научный scope и решение

Physical Π/μ/LO/static formulas не менялись относительно independently reviewed V1; новый freeze меняет correctness/accounting/output, а не source Core, допуски, energy grid или closures. Φ restriction остаётся common-μ FD, explicit Hilbert PI, без Cdirect; ηseed — только initial Γ. Algorithm.md делает эту схему понятной отдельно от I/O guard scaffolding. После trigger history/old Green освобождаются barrier/GC, fresh native warm далее не подвергается FD projection. `scientific_acceptance=false` остаётся неизменным.

Small checks снимают конкретные compatibility/programming blockers, но actual initial Π+native-five/KMS+cavity equality, whole model HDF, runtime/RSS и convergence ещё предстоят. Сначала применить сохранённый threshold1e−10 к actual initial control, только затем до16Π. При fail/domain/cap сохранять отрицательный исход, не увеличивать бюджет и не менять tolerance внутри attempt. Native32 только условно по raw reduced trigger; даже native strict fixed-U pass потребует самостоятельной physical/discretization/Poisson/experiment интерпретации.
