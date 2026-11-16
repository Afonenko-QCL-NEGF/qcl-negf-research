# V2: serialization и raw preservation, без изменения derivative

V1 целиком сохранён в `unexecuted-v1/`: module0271, adapter8501, все входы,
proof/protocol и оба первоначальных manifests — побайтные копии. Дополнительный
`archived-files-SHA256SUMS` имеет archive-local пути. V1 actual не запускался;
HOLD independent review P1/P2 не является numerical failure/outcome.

**P1.** Удалён `plain(AbstractArray)=map(...)`, сохранявший Julia Matrix как
Matrix, которую TOML не поддерживает. Новый source-free `serialization.jl`
кодирует real Matrix как `{format,shape,rows}` с вложенными row vectors;
complex Matrix — `{format,shape,real_rows,imag_rows}`. Порядок строк и исходная
форма явны. Vectors/NamedTuple/Dict рекурсивны; Complex scalar имеет real/imag;
необъявленный tensor ndim>2 отвергается и должен идти через H5 raw writer.
Таким образом actual Hbar сериализуется с shape, без transpose/conjugation.

`serialization_checks.jl --tiny-serialization` — отдельный stdlib-only fixture:
rectangular Float64 и complex2×2 matrices, bitwise UInt64 round-trip (включая
signedzero,1e−300,2e200), actual-like vector/records fields и nonoverwrite.
Не импортирует Core/Unitful/HDF5 и не строит model. Adapter `--tiny` вызывает
тот же fixture до mathematical tiny checks. Автор ничего не запускал.

**P2.** `fresh_phi(...;persist_raw=nothing)` вызывает optional persistence сразу
после native raw Φ contraction/summary, **до** `GP.cone(raw)` и до `keep_base`
extra PI reconstruction. Adapter задаёт read-only callback для **всех пяти**
base/probe Φ. Callback сохраняет full raw H5 и metadata до guard; после return
нет duplicate raw write. При cone failure или позднем extraPI failure raw
остаётся; если callback IO сам неуспешен, cone/extraPI не вызывается и причина
не скрывается. Actual tested Γ/raw Φ callback не модифицирует.

`Phi_completed` теперь точно означает «raw contraction произведена», не
«cone/IO/return/extra-PI admitted». Инкремент стоит до callback/guard; metadata
и summary содержат этот scope. `Phi_entered` по-прежнему после input cone,
failed calls остаются расходом cap. Отдельный tiny behavioural fixture
проверяет snapshot-before-reject, порядок persist→guard, немодификацию raw,
completed scope и IO failure→no guard. Ни root, ни автор не объявляют
callback/machine cone строгим PSD certificate.

Другие mathematical code paths `spectral_tangent`, covariance3terms,
implicit μ, native5 channels, central differences, Arnoldi/reorth/SVD/GMRES
coordinates **не изменены**. Γ update/native calls всё ещё0. Constants для
module/serializer обновлены, manifest включает новые sources/tests. Новая
identity суффиксv2-IO;12.125GiB/600s/5Φ/4Jv/3GiBdata/2GiBfile неизменны.
Матрица Hbar в TOML теперь dictionary: читатель использует `Hbar.shape` и
`Hbar.rows`, не предполагает top-level2D native TOML array.

Весь V2 остаётся unexecuted-author, runtimefalse. Root source+tiny+independent
delta review обязательны до actual. Научная acceptance всегдаfalse; прежние
conditional-PD/finite-difference/nonnormality limitations полностью действуют.
