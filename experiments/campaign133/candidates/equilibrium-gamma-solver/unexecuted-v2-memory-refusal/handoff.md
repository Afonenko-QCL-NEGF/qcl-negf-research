# Γ pilot16 — V2 freeze после independent reviews

V2 реализована только в этом каталоге; Core/Runner/source environment не
редактировались. Автор Julia/Pkg/compute не запускал. Runtime V2 и RSS —
`not_measured`; input `runtime_authorized=false`. Root повторяет small checks
на frozen V2 до любого actual S01 map. Физическая схема и caps не менялись;
input отличается только metadata identity `...pilot16-v2`.

V1 сохранена **до правок**, исходные SHA проверены: `unexecuted-v1/` содержит
все6 source/input/doc files, неизменный original SHA256SUMS, archive-relative
SHA256SUMS.archive, оба final review и копию Root gamma-small-checks.log.
Root V1 outcome: ParseError checks.jl:86 (`.2.*x`), exit1,21.98s,
peak1029512KiB; **actual S01 maps0, small-check pass отсутствует**.
Исходный log остался environment/accepted/gamma-small-checks.log.

Accepted Core PR30: HEAD `e171caa60cce676c6c4238dd671e7c527f2bd977`, tree
`742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8`, source worktree
`core-accepted-science`, clean. Root native Pkg environment `environment/accepted`
подключает его relative develop path; unused Runner удалён Root, harness не
имеет Runner/Pkg dependency. Gate проверяет actual source/version/resources;
стандартные Project/Manifest архивируются при actual run.

Одна V2-дельта (`v1-to-v2.patch`, полный SHA в SHA256SUMS):

- Однозначное `0.2 .* x`, tiny asymmetric control: raw Φ−Γ управляет trigger
  и best; Hermitian working AA residual и representation correction отдельно.
- Узкий known inverse predicate: Singular/LAPACK, nested TaskFailed/Composite/
  CapturedException только если **все** причины численные. Method/Bounds/
  budget/произвольные ошибки не превращаются в fallback. Одна convex fallback
  α=.2 остаётся counted; failed fallback останавливает attempt.
- Packed BitArray model/stencil masks пишутся dense Int8; tiny HDF roundtrip.
  Actual AlgorithmOptions добавлены в model metadata.
- Initial per-channel KMS absolute/scale/relative, weighted wE*wk и maximum
  block error/locator. Precheck fresh-candidate embedding не обещана/не
  добавлена: она только в conditional native materialization/final audit.
- Final per-family R/L/G fresh−stored norms/scales — typed streaming, пять
  channels+embedding total/±. Collision и ΓR−ΓK identity тоже включают embedding.
  Native max rSigma с прежними definitions сохранён отдельно.
- Final FDR явно именует full-domain joint pair norm/FDpair denominator;
  добавлены weighted, occupied lesser, trusted weighted и pair/A показатели
  на существующих arrays. Это posthoc observation без FD projection/maps.
- Materialization/final-audit start отмечаются до call, completion отдельно.
  Public solver observer означает completed maps; reserved cap30 даёт честный
  started lower/upper bound. Структурный максимум1+30+1=32, не утверждение,
  что observer знает начавшийся и упавший внутри map. События остаются в summary.
- AA объявлен feasibility-safeguarded **NONMONOTONE**, без отсутствовавшего
  regression threshold; fixed16/caps прежние. QR сохраняет firstpair/Rdiag/scales.
- `algorithm.md` объясняет unknown/map/FD restriction/AA/native exit с25
  строками Julia-style pseudocode; I/O/resource policy выделена отдельно.

API impact: **0 Core/Runner files, 0 public keywords, 0 dependency changes**.
Три локальных Julia files: main `run_pilot`, typed operators/stages,
small controls. Native source pinned e171 и explicit PI/Hilbert; directGreen
не включается. Γ-stage barrier+GC освобождает history/QR/old Green перед
native; consumed warm Green освобождается перед fresh final audit.

Root first command после static V2 review в accepted environment:

```julia
include(".build/campaign133/equilibrium-gamma-solver/run_gamma.jl")
small_checks(TOML.parsefile(".build/campaign133/equilibrium-gamma-solver/pilot16.toml"))
```

Эти controls включают prior CP/QR/nonlinear/PSD/μ/writer controls и новые
raw-vs-working, wrapped numerical-vs-programming exceptions, packed-mask
roundtrip, accounting bounds. Автор их не исполнял. После их pass и resource
review Root создаёт authorized copy input с новым SHA; normal entrypoint
`run_gamma.jl INPUT NEW_OUTPUT_DIRECTORY` сначала actual one-map native
comparison, только затем loop. До этого никаких actual S01 execution claims.

Ресурсы: hard12.375GiB/high11.5/Swap0/CPU400/heap8G; globalΠ≤16 включая seed,
failed maps и fallback;180/600/900/1800s; file≤5GiB/data≤6GiB. Γ-only checkpoints,
initial rawΦ/native sums, terminalΓ+GRA/GFD/rawΦ, model/kernels/weights/stencils;
conditional full native state≈3.52GiB. Data admission не удаляет evidence.
При stop между maps lastG может отсутствовать; last validΓ/best/failedΓ
сохраняются, дополнительная Π сверх16 ради output не выполняется. Fresh native
candidate arrays не дублируют второй3.52GiB файл: full terminal G/Σ+source/model
и scalar audit позволяют воспроизведение. `scientific_acceptance=false`;
fixed-U native strict pass не принимает Poisson/discretization/experiment.
