# Конечный adapter к frozen analytic Jv: новая дельта

Этот документ дополняет исходный `protocol.md`; прежние module0271/input/proof
не изменены. Новый `run_diagnostic.jl` — finite adapter, не новый solver.
`adapter-input.toml` false; actual допускается только после root source review,
tiny и independent review. Автор не запускал Julia/Pkg/H5/численные тесты.

Два входа CLI: `--tiny` либо `AUTHORIZED_INPUT NEW_OUTPUT`. Пути frozen frame
в config разрешаются относительно каталога adapter, независимо от места
authorized copy. Actual copy может отличаться от frozen `adapter-input.toml`
**только** `runtime_authorized=true`. Вход hash сохраняется отдельно. Нельзя
менять amplitudes/χ/μ/core/geometry/критерии ради запуска.

До model arrays проверяются: Core e171/tree742 clean tracked; Julia1.13.3;
manifestda078; все транзитивные parent sources по pinned identity973d;
literal module0271/stream8df/parent50bc/spectralEECA; новый adapter manifest;
Γcab810/model700a/terminal6eb3/input8d5c/map073d/summary760c hashes. Перед build
проверяются checkpoint schema/parent identity/call64/axes и model schema.
После model build используется прежний `GP.verify_problem_model` (whole basis,
grid/profiles/H/kernels/qK). Γ загружается без repair/interpolation/rescale.
Новая identity hashes **архивированных** sources после copy, до model build.
Native Context/SCBA/acceptance не вызываются и не меняются.

`--tiny` source/manifest gate сохраняет, actual OS envelope gate пропускает:
это только small synthetic full45+number+native-channel linearity tests.
Root выбирает один отдельный bounded tiny process, например J1/BLAS1,
MemoryMax2GiB/Swap0/CPU100%/wall60s; это allocation suggestion, не измеренный peak.
Source import выполняется через world-age-safe latest invocation. При tiny
ошибке actual запрещён. Tiny не трогает actual Γ/H5/model frame.

Actual OS caps ровно из задания:12.125GiB13019119616B,
High11.5GiB12348030976B, Swap0, J4BLAS1/CPU<=400%, external timeout600s,
LimitFSIZE<=2GiB, data<=3GiB. Guard проверяется через active cgroup и proc;
native numerical memory guard остаётся parent12.125GiB. Wall START перед
package imports; internal Limits включает hash/model/JIT и проверяется перед
новой фазой. Hard timeout внешне обязателен. Peak/time не измерены.

Сначала base fresh Φ и fresh current μ. Сохраняются raw baseΦ и condition/
spectral/μ summaries, затем **прежний** whole-field rawΦ64/nine-scalar256eps
restart replay с bitwise flags. Replay fail останавливает дальнейшие maps.
Далее четыре fresh Φ для `(1±1e−3)Γ`, `(1±5e−4)Γ`; каждый raw Φ+summary
сохраняется до следующего вызова. Dε сохраняются отдельно. Γ scaling PSD-safe
при PSD source; global exact-PD proof остаётся pending и так помечен в identity.

Один analytic JΓ сохраняется full raw до qualification. Numeric gates каждого
Jv: finite primal G/A replay, PI roundoff, raw/Herm conversion и δN relative
<=1e−10 (прежняя arithmetic target, не SCBA criterion). FD gate отдельно
**эмпирическая**: fine error не больше coarse/fine difference в absolute и
weighted norm (или уже <=arithmetic relative); fine error не растёт против
coarse. Leading central truncation errors имеют ratio4 и fine error=diff/3;
без high-derivative/roundoff bounds эта gate **не rigorous error enclosure**.
Actual raw comparisons сохраняются независимо от gate. Fail → no Arnoldi.

При gates проходят максимум три дополнительных analytic Jv (Arnoldi3), из
Herm(raw residual), без trial Φ вне PSD cone. Каждый завершённый Jv ledger
пишется callback **до** следующего Jv/guard; small Hbar/Ritz/residual/SVD и
Q≤4 сохраняются при завершении. GMRES только coordinates/predicted restricted
residual, Γ correction не формируется/не применяется. Short Ritz spectrum не
сертифицирует full spectral radius/nonnormal stability/grid accuracy/root.

Перед каждой Φ/Jv фазой есть immutable request TOML с entered lower/upper и
completed counters. При hard kill finally может не выполниться: request
структурный upper bound, не выдуманное exact completed count. При обычной
ошибке summary записывает actual entered/completed, exception и elapsed;
причины не маскируются fallback. Programming errors также останавливают процесс.
Φ pre-cone rejects до entered остаются requests, не называются начатыми maps.

Raw output ожидается ~12 full fields≈1.51GiB + metadata/source archive; каждый
file≈0.126GiB, padding учитывает прежний atomic H5 admission. Это оценка по
заданной форме, не actual data/peak measurement. Partial raw/tmp сохраняются.
На source gate fail до output allocation причина остаётся в captured root log;
никакого scientific pass. В каждом summary `scientific_acceptance=false`,
`native_calls=0`, `Gamma_updates=0`, даже при успешном diagnostic.
