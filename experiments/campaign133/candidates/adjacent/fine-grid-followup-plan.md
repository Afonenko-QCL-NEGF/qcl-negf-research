# Следующий научный шаг после stableA64: ограниченная fine2x Γ-итерация

2026-11-16. **Read-only design; active stableA-cont64 не затрагивается.**
Core `e171caa60cce676c6c4238dd671e7c527f2bd977`; новые расчёты, HDF5 reads,
изменения frozen code и реализация не выполнялись. Superproject HEAD
`b515c42ef14c1f61151626280da0bcf100d116be`, чужие gitlinks сохранены.

**Рекомендация: сначала один bounded fine2x stage с private Γ-only cache,
qualified stableA и AA history1.** LO-closed phasegrid/Cauchy-tree оставить
следующим условным вариантом: сейчас uniform2x отделяет nonlinear accommodation
от ошибки coarse sampling с меньшим числом новых операторных предположений.
Даже near-root coarse и fine не означают continuum accuracy или native pass.

## Что уже установлено и что ещё не различено

На **одном Γ128 input** two-map control дал weighted raw residual
`5.95755e−6 → 0.00265386` при NE30704→61407. При source μ число изменилось
на `+7.1111e−5 Ntarget`, refit μ — на `−0.4946 µeV`.
HPI(PΓ) на old nodes изменился лишь относительно `1.02e−14`, ghost-bound
`4.75e−20` native energy. Источник:
[summary.toml](../refinement-control/attempt01/summary.toml), §§fine_map,
HPI_source_node_difference, fine_cache_admission; обе карты completed, хотя
исходный postprocess затем остановился с `InexactError`.

Это сильная **чувствительность discrete raw map** к h/2 на данном input,
не доказательство конкретной linewidth или единственной причины стагнации.
Новые midpoint resolvents, density quadrature/current-μ и LO/static contractions
участвуют совместно. Малое изменение occupied number не позволяет снять
residual в пустой области: LO связывает её с occupied ladder.
StableA qualification44 Big/ΓPD и same-trajectory16 сообщены ведущим;
они улучшают arithmetic tails, но не доказывают разрешённость сетки.

## Единственный follow-up и различающее свидетельство

После остановки active64 зафиксировать **terminal либо last-admissible Γ**,
его SHA/μ/source/operator identity. Правило выбора объявить до чтения outcome;
не выбирать наиболее красивый checkpoint. Новое исследование — та же all-five
модель, common70K, U0, basis/k, finite4 central-driven repeated Σ, η=0;
меняется energy discretization, history запускается пустой. Исходную positive
Γ prolongировать convex piecewise-linear в actual coordinates; не переносить
Gr/A/μ/history. Если source только допущен слабой cone tolerance, это не strict
PSD: сохранить eigen/backward witnesses, не clip и не объявлять proof.

Один процесс J1/BLAS1, до **14 map starts всего**, включая failures:

1. Один fresh stableA coarse replay выбранного source. Зафиксировать исходные
   raw norms, spectral/cone/number ledger. Rebuild fine61407 с прежними endpoints,
   bitwise old-node nesting, ΔE/2 и native LO plans (nominal939→1878).
2. Максимум12 fine Π, **включая первый fine seed replay**, AA-II history1,
   прежние α/β/coefficient/cone controls. Один штатный convex fallback при
   infeasible AA тоже расходует этот cap. Raw Φ−Γ и original weighted/unweighted
   denominators остаются stopping evidence; checkpoints 1/4/8/12, без keep_channels,
   Green dumps, hidden best-field copies и full-native precheck.
3. Один terminal coarse map на RΓfine (odd-node restriction), если не было
   resource/domain stop. Отдельно измерить `Γfine−P(RΓfine)`: этот coarsening
   меняет continuous source, поэтому итоговая coarse/fine разность включает
   потерю midpoint structure и не называется same-source map defect.

Primary question: способен ли fine map **самосогласованно поглотить** начальный
~0.00265 defect за ограниченное число evaluations, сохранив causal/number/stencil
witnesses? Сохранить raw residual reduction, changes Γ и μ, support около
прежнего anchor9785 и его ±Ω ladders, включая even-node contribution; полные
нормы обязательны. Если fine residual уменьшается — coarse solution нуждается
в сеточной коррекции, но fine accuracy не доказана. Если fine stagnation/PSD
barrier воспроизводится — гипотеза «достаточно удвоить N» этим budget не
поддержана; это не доказательство отсутствия fine root. Граница PSD может быть
свойством AA/state или arithmetic, не обязательно физическим запретом.
Wall cap900s/soft840s; 12 fine maps по прежнему22s map-only — лишь scale для
admission, **не прогноз** stableA nonlinear runtime. Неудачные starts/retries
не бесплатны. Сохранить outcome `budget_reached/domain_failed`, если нужно.

## Честный cache и live-set, без обещания native exit

Переиспользовать owning constructor recipe
[`operators.jl:4–23,39–61`](../refinement-control/attempt01/operators.jl), SHA
`7d71dadce03edc17bcd7c68bae2fd5d1057e17cf95aae66432859ce4caa4a4dd`:
exact native kernels, compact shifts, one-worker/one-column FFTPI plan;
cache.estimate хранит **настоящий** full-native estimate. Старый envelope
`24F+3GiB` относится к history0 two-map scope; нелинейный stage требует новой
явной ledger, не переименования прежнего admission.

F=`270190800 B` (full complex fine field); packed-Hermitian storage не
предполагается. Actual stable_spectral state :53–103 удерживает примерно21F:
Γ; Λraw/ΛH/localR; probe L/G; emb/plus/minus по3; totalR/physical; Gr/A/Asub/
effective. Это крупнее одного положительного Gram output. Actual AA хранит
states/residuals по m+1, QR columns m и отдельные mapped/best/trial:
[`gamma_functions.jl:292–365,695–751`](../equilibrium-gamma-stableA-cont64/attempt01/gamma_functions.jl).
Conservative accounting поверх прежнего24F envelope: `2(m+1)F` history,
`2F` mapped/best и `mF` unreclaimed QR; trial-input уже в spectral envelope.
Без уменьшения frozen lifetimes это даёт:

| Вариант | Envelope с3GiB cache/JIT/GC allowance | Решение относительно hard12.125GiB |
|---|---:|---|
| history0 control | `24F+3GiB = 9,705,804,672 B` | прежний scope; не nonlinear admission |
| history1 | `31F+3GiB = 11,597,140,272 B` | ниже cap; ещё нужна actual lifetime/RSS проверка |
| history3 | `37F+3GiB = 13,218,285,072 B` | **выше cap 13,019,119,616 B; отказаться** |

Это conservative allocation estimate, **не измеренный bound RSS**. Удаление
best-copy/QR retention может уменьшить peak, но требует отдельного reviewed
lifetime design; экономию здесь не предполагаем. Перед запуском — доступная
память всей VM20.6GB, отсутствие concurrent heavy work, реальный MemAvailable
≥cap+1400MiB. Реальное enforcement: MemoryMax12.125GiB, High11.5GiB, Swap0;
scoped output≤3GiB, один файл≤1GiB, checkpoints Γ-only. Coarse context и крупные
buffers не должны переживать fine stage barrier. Не резервировать всю VM.

**Full-native estimate `25,839,982,984 B` превышает и cap, и VM.** Сохранить
`native_validation_refused_for_full_budget=true`; public guard/Core budget
не менять. Этот опыт вообще не содержит materialization/native exit. Даже
rΓ→0 означает только reduced fine fixed point; при resource refusal native
rK/rΣ/λ/free release остаются `not_measured`, не pass. Последующий native
endpoint потребует отдельно разрешённого streaming verification либо другого
ресурсного окна, а не повышения числа в options. Restrict на coarse и принять
coarse native нельзя выдавать за fine acceptance.

## Почему phasegrid/Cauchy — пока второй вариант

[Существующий phase-ladder plan](residual-support-refinement-plan.md#future-alternative-lo-closed-local-midpoint-ladders)
сохраняет все coarse узлы и добавляет целые LO ladders выбранных fine midpoint
phases: порядка33 nodes на phase в этом окне. Это может экономить mesh, но
число нужных phases/ширина features пока неизвестны. Shared all-k union,
новые положительные weights и actual Float64 LO matching обязательны;
pointwise KMS не доказывает weighted-adjoint/particle conservation.

[Cauchy-tree prototype](cauchy-tree-prototype/frozen-attempt01/summary.toml)
прошёл transform-only selected-target diagnostics, **density_SCBA_adaptive_economy
=not_measured**. Nonuniform caller должен сохранять explicit ghosts, near-panel
PV cancellation, far-error budget и analytic input contract. PSD hats дают
causal **точный** Cauchy integral; truncated moments сами по себе не гарантируют
каузальность или малую nonlinear map error около narrow poles. Готовая tree
не устраняет mid-energy sampling error и не является completed adaptive SCBA.
Если bounded fine1 не помещается либо новые even-node features доминируют без
достаточного прогресса, следующий вопрос — frozen all-k LO-closed phasegrid
map accuracy против fine uniform, до новой nonlinear/adaptive кампании.
Direct-static hybrid остаётся прежним отдельным future protocol.
