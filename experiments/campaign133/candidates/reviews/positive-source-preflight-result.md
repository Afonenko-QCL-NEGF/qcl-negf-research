**Независимый result review positive-source qualifier, 2026-11-16**

Решение: **узкое исполнение оператора и supplied-G consumer control
подтверждены; дальнейший stable-A-only pilot допустим как отдельный
ограниченный арифметический вариант по root admission**. Это provisional
qualification, а не доказательство intended causal full45 accuracy или
native fixed-point acceptance. Selected split-input Big oracle готовит
отдельный автор b0_review; код/результаты будут проверены независимо.
Нелинейного преимущества covariance этот один опыт не установил.

Первичные evidence: `positive-source-preflight/attempt01/summary.toml`,
SHA256 `06d6fa46c058187a8368663c9d19ef3f3edf3ab3582a5cf629286a2b9e231367`;
неизменный script `a34b1da110293ba289a78a437354436970ed8362a184779eb99e1f87690ade18`;
runtime input `de3ef182dd039cffb46a1cedf918a87ba2dfec212e98c703a7619952206c1d18`;
selected.h5 `a74dd7cb0afdfdca6f37af61f0ba1b43201306249816f3acc87fb080f1c8079a`;
`attempt01.log` `0fb465480f1b85475b278a192e8f1311f470b3cab1df3c117a4209d08938638f`.

Reviewer обработал только малые TOML/log и inventory/indices selected1.3MB
через existing h5ls/h5dump для подготовки передачи Big oracle. Whole Γ,
source/model HDF5 и Julia не открывались/не запускались; пакет/код не
устанавливался. stdlib Python-аудит с RLIMIT_AS128MiB сохранил
`reviews/positive-source-result-audit.json` и воспроизводимый
`reviews/positive_source_result_audit.py`. Проверены31 norm ratios и
scalar stage/status identities; **не** пересчитаны full-field нормы и
глобальные eigenvalue minima по отсутствующим полным массивам.

**Что действительно произошло**

`execution=complete`, maps_started=maps_completed=3: два reduced raw Φ
(covariance и matched subtraction) плюс один native candidate consumer.
Оба fit μ и оба Φ завершены; имеются μ_sub/Phi_sub_cone/matched norms и
selected обеих raw Φ. P2 broad-catch caveat из static review **не сработал**:
`matched_sub_control.not_measured` в результате отсутствует. Здесь нельзя
перепутать uncomputed control с computed-but-unaccepted control.

Actual MemoryMax12.125GiB/High11.5GiB, J4/BLAS1/CPU400%,300s/Swap0/
output32MiB заданы admission; программа завершилась за26.8492s internal,
34.47s whole process, peak6227284KiB, Swap0, exit0. Data около1.3MiB.
Runtime input отличается от frozen только `runtime_authorized=true`;
параметры/модель/seed и thresholds не ослаблены.

Exact h-policy проверяется результатом: hH взят из прежнего native plan,
side-repeat defect0; Δh_F=1.6335754e−15 при scale16.4803651,
relative9.9122524e−17. Λ projection relative2.5072768e−16,
weighted2.5072153e−16; PI padding witness3.4934751e−16.
Это не bound эффекта на самые слабые linewidths. Selected preassembly ΛH
сохранена отдельно от assembled localR; tiny Γ могла потеряться при
assembly/retarded propagation, поэтому Big oracle должен promote split
hH/ΛH/Γ/T **до** сложения.

**Matched результат при фиксированных Γ/hH/GR**

| Witness | Covariance A / Φ | Subtraction A / Φ |
|---|---:|---:|
| Minimum computed spectral eigenvalue | +2.5949557e−29 | −1.4648561e−15 |
| Minimum computed raw Φ eigenvalue | +2.4945639e−31 | −6.6180748e−20 |
| Minimum q | +4.5866979e−32 | +4.5866956e−32 |
| Negative/nonfinite q | 0 / none | 0 / none |
| Current μ (native E0 units) | 0.08464025621945251 | 0.08464025621945251 |

Оба q допустимы и обе μ **точно совпали** в saved Float64 representation;
number/capacity тоже совпали (capacity28.64390675, target0.045).
Следовательно этот matched контроль **не показывает, что covariance одна
устранила отрицательные q**. В сравнении со старым HermΛ trajectory вместе
изменена central h-политика: теперь hH во всех9 cells, раньше центральный
hraw. Прежний STOP нельзя задним числом приписать последнему вычитанию A
или считать повторённым одним matched Asub.

Глобальная разница A_cov−A_sub: absolute9.1335109e−12,
relative1.1545006e−15, weighted1.0103576e−15. Разница raw Φ:
absolute1.5477284e−14, relative8.9306847e−16, weighted8.6100937e−16.
То есть source product улучшил знаки очень слабых хвостовых witnesses
при изменениях всего поля порядка машинной точности. Обе старые cone
проверки допускают свои tiny negatives в backward-error budget; это не
новый строгий positivity criterion. У covariance A максимальная absolute
anti-Hermitian defect4.7998e−13: generic multiplication не является
strict FP Gram primitive. Малые положительные computed eigenvalues —
floating witnesses, а не рациональный/interval certificate PSD.

Фitted number relative error−1.1611082e−13; отдельный native number
functional на том же G<=ifA даёт−2.0805579e−13 relative к target.
Разность порядка суммирования не скрыта под normalization; A/Φ не
рескалировались. Не нужно превращать этот roundoff-level witness в
отдельное разрешение менять scientific tolerance.

**Native consumer и исходные уравнения**

Пять enabled families сохранены. Независимо recomputed maximum из
per-channel kinetic/KMS scalar records5.6869146e−16, меньше declared1e−10;
максимум weighted kinetic errors5.0616412e−16. KMS joint global norms
порядка8e−17–1e−16; maximum-block witness LO8.7353e−15 также записан,
но не смешивается с выбранной global нормой. LO kinetic difference0
здесь от общего call path, не независимое подтверждение LO формулы.
Retarded sum vs HermΠ[rawΦ] relative2.7616483e−16,
weighted2.7410997e−16. Это хороший **same supplied stable-A Green**
consumer control, не fresh unconstrained Dyson/Keldysh/native exit audit.

Raw fixed-point residual одного cov Φ−input Γ остаётся0.0088999078
(weighted0.0088372481), существенно выше root target1e−9. Solver не
выполнялся, никакой нелинейной сходимости этот result не предъявляет.
Ни positive q, ни маленькие matched differences не принимают SCBA,
discretization, window/basis, Poisson или экспериментальную валидацию.

**Primary selected scope и следующее различающее свидетельство**

Сохранены E-indices7893/8038/30696/30704 и все11 k, full5×5 split/assembled
операторы, обе A и raw Φ. Это покрывает worst Asub e8038/k1 и minimum
Acov e30704/k1. Worst Φsub matrix находится e4279/k1 и **не попала**
в selected energies: его глобальный negative minimum сейчас только
summary cone witness, не independently inspectable primary matrix.
Ограничение сохранено; повтор wide probe ради него здесь не требуется.

Sole Big-oracle author — b0_review. Правильный независимый input —
promoted split hH/ΛH_preassembly/Γ/T во всех9 cells, а не только promoted
assembled ΣR, где часть Γ могла исчезнуть. Он проверит selected intended
full45 response/source identity и precision repeat, с explicit Γ positivity
qualification. До этого exact-causal-model accuracy остаётся частично
unresolved даже при улучшении знаков. Reviewer не написал дублирующую
реализацию oracle.

Следующий ограниченный stable-A-only16 pilot должен иметь отдельную
identity, тот же original seed/hH/source construction и пустую AA history;
PSD trial projection не объединять с ним. Raw Φ не проектировать, native
handoff/acceptance не заменять supplied-G consumer pass. Root admission
ресурса и Big-oracle evidence определяют дальнейшее исполнение.
