# Независимый static review: Hermitian density V4

**APPROVE implementation/source-gate delta и отдельный root freeze-IO запуск** в заявленном внешнем лимите ≤2GiB, Julia1/BLAS1, timeout120s, output32MiB. Блокеров на этих frozen текстовых исходниках не найдено. **Actual Matsubara остаётся HOLD**: текущий input не заполнен, `execution_authorized=false`; данное решение не принимает плотность или SCBA. Reviewer не запускал Julia/Pkg, не открывал HDF5/массивы и не менял проверяемые файлы.

## Frozen identities

Каталог `adjacent/gamma-density-hermitian-v4`, независимо пересчитанные SHA256:

| Файл | SHA256 |
|---|---|
| matsubara_oracle.jl | `814dff509cef6cdb1e0e1d1a762cabbd885a761c2a5b4964f80373df2b01b21a` |
| freeze-stable64-io.jl | `192e6f83bd808179db98beecb375642c255480b41a9349d0644ca4c8d32165c6` |
| input-stable64.toml, pending | `916a191a11c704c2a9329985f1f4054545c6170a0cbe8107b7b7d086ae589ff4` |
| protocol.md | `bf451daeeb9728e05c35481f76e20850b60328fc64dabffa7f7e8945de99a3fc` |
| certificate-gate-controls.jl | `93297c5e77659c6480b1c8b6791534b6a2268ca2aee4bf10a64db18473ef8268` |
| v3-to-v4.diff | `4727c72089d31b271eeec28dece3f991b3e4a2dda98e46df89ade1ea0df50e36` |
| freeze.toml | `42a9f2c52983bc134531167c25bfc64e5d5e5d63ef2aad1cc73739278d63a3aa` |

Численный parent5249 и V3 `e45ac4908f5ea3487aff8e1a5580d1ccc9277fb8bc16704d15bc143b3428c121` остались неизменны. Независимое byte comparison подтвердило:

- первые177 строк V4 = V3 = parent, SHA `c0af40bc7b1ac14773201a848eb1fb5ddf28541d2794f3a045f3ec1a0357205e`;
- весь suffix от `checkstage("independent_kernel_and_bare_controls")` до EOF V4 = V3, SHA `5dd646a0d3dbc1e787d4f91b2c77f98119f87d36679cd9bd0f30eccd2d51a83a`.

Таким образом, hat series/BigFloat/remainder, moments/tails, complex bare controls, full45 LU/central5 RHS, Herm(G), accumulation/refinement/output arithmetic не получили новой implementation delta. Прежняя численная qualification сохраняет прежние ограничения; новая h/source identity требует собственных binding/runtime witnesses.

## Source policy и all-nine h

`matsubara_oracle.jl:225–327` и helper:45–118 используют known writer-owned schemas/datasets, а не произвольный configurable export. Whole-file hashes проверяются до HDF; realpath ownership, input/identity/Coree171-clean, call64, unprojected source, spectral operator, saved ledger и completed source summary связывают packet с одним attempt. Малые input/identity/writer/map/ledger/summary hashes и accepted Project791aec/Manifestda078 пересчитаны reviewer и совпали. Большие HDF5 hashes reviewer не пересчитывал; их проверяет actual helper/oracle.

h для всех nine cells берётся ровно из **saved `spectral-selected-064.h5/h_H`**, finite ComplexF64 shape11×5×5 и elementwise exact Hermitian. Raw model h только crosscheck с saved h_raw; новая symmetrization не выполняется. Tminus проверяется elementwise как exact adjoint Tplus. `chain_hamiltonian` повторяет saved h_H в9 diagonal cells и задаёт сопряжённые offdiagonal couplings; tolerance guard H0 не является основанием exact-Hermitian доказательства, поскольку H0 в эту сборку H45 не подставляется.

Это соответствует проверенному source `stable_spectral.jl::stable_context`: копия native `c.cavity.hlocal` применяется центрально и к sidecells, а `spectral-selected` сохраняет именно s.h. Исходный native raw-h путь остаётся отдельным контрактом, не получает приёмку задним числом.

Known `gamma_functions.jl:649–661,737–763` пишет terminal Gamma=input trial, μ и Green одного evaluate; после принятия тот же trial становится last gamma/checkpoint. Дополнительно helper:104–114 и oracle `checkpoint_terminal_equal` сравнивают **все real и imag entries** terminal6eb3… и checkpointcab810… slabs≤1024E, с exact shape/type/schema/call/identity checks и `isequal` (включая signed-zero provenance). Это не selected-subset comparison и не сравнение только norms.

## Freeze-IO: read-only и ресурсная область

Все HDF5 opens helper имеют mode `"r"`. Пишутся только новый output directory, копии малого input/helper и digest TOML. Source input не изменяется; `execution_authorized` не включается. Существующий output отвергается. Никаких solve, Φ, fit μ, density, PD inference, projection или нового H не выполняется.

Logical digest expression в точности совпадает с V3 `sha256(reinterpret(UInt8,vec(x)))`; чтение через HDF5.jl сохраняет native (E,k,a,b), без h5py axis reinterpretation.13 hashes =4 grids+saved hcell+7 basis+Gamma, raw-h отдельный;11 A hashes считаются по одному k. PI geometry проверяет native `range` equality/positive steps и сохраняет origin/spacing/witnesses.

Один полный complex Γ имеет8443600 entries /135097600 bytes; construction с двумя real temporaries ≈270195200 bytes payload. `digest_one` снимает ссылку и GC перед следующим массивом; A по одному k, checkpoint comparison — slabs. Это не измеренный Julia RSS; import/JIT/HDF buffers/allocator требуют внешнего MemoryMax2GiB и swap policy. Internal105s включает whole-file hashing, но не прерывает blocking I/O; внешний120s нужен обязательно. Output маленький структурно; финальная32MiB проверка не заменяет внешнего лимита. Helper runtime/RSS/API **not_measured** на момент этой static проверки.

## Original global ALL certificate и strict actual gate

Прочитаны original native `global-gamma-exact-pd/all-attempt01/summary.toml`, exact hashed proof `all-run.toml` и exit log. Certificate SHA `a703ec4c3e2feda1ae9571ab3d1d20b69dd9de48898f3d620ada008da2583335`; proof SHA `27da5a97ab61e9f290f5f8c7a2979d5873ec55e81f45e6b336fba99c495711d4`. Native process exit0,17.45s,RSS1180640KiB,swap0; internal16.141s. Это original TOML, не JSON projection и не benchmark5000.

Cert сообщает complete/all, expected=requested=checked=PD337744, all_blocks_strict_PD_certified=true,20264640 exact divrem components=60×337744, нулевые nonpositive/nonhermitian/nonfinite counts, empty failures и отсутствие failure field. Runner `d2395572…` и math `c6b14339…` hashes пересчитаны и совпали с cert и pending V4 config. Cert input hash равен exact proof hash. Все9source hash triplets cert/proof/V4 равны, включая checkpoint **cab810…**, model700ae…, input/identity/writers/ledger/summary/map. Mathematical PD implementation reviewed отдельно в `global-gamma-exact-pd-static.md`; новый review не переисполняет337744 proofs.

`require_all_PD_certificate` (oracle:188–211) действительно проверяет именно эти поля, фиксированную shape30704×11×5×5 и count337744, call64/source identity, exact proof/cert file hashes и все9 bindings. Partial, failed, benchmark, requested/checked subset не проходят. Main проверяет `execution_authorized===true` до source HDF чтения; gate повторяет его. Source-free negative controls статически соответствуют gate, но reviewer их не запускал. После certificate oracle всё равно требует full checkpoint↔terminal equality, а не переносит proof на другой Γ по декларации.

## Следующий допустимый переход и предел вывода

После completed freeze-IO нужны точный digest→input перенос13+11 hash dictionaries, raw-h hash, origin/spacing, binding source/env/μ, original cert/proof path+SHA; новый input hash и независимая small-text binding проверка. Текущие UNSET values/false flag — корректный HOLD. Запуск helper не разрешает actual автоматически. Actual Mats также требует отдельного root resource-green; текущая static проверка не принимает ещё не выполненный numerical result.

Математический no-UHP-poles вывод относится к exact frozen model: PD stored Γ nodes → PSD nonnegative full-hat interpolation → −ImΣ(z)⪰0; exact Hermitian H45 даёт ImD45(z)⪰Im(z)I>0. Truncated-series/Float64/LU errors не становятся сертифицированными этим доказательством. Возможные real-axis atoms отдельно не исключены; N512/N1024 refinement и M1/M3 — witnesses, не formal total error. Нет SCBA convergence, grid adequacy, native acceptance или experimental validation.

## Отдельная docs-only closure final pair

По steering root дополнительно прочитан исправленный `environment/final-integration/HANDOFF.md` SHA `f0b333e2d826d3684db700a93487622197ad279a71d4c0ab52fdc45dbbb94ded`. Он теперь явно описывает actual45.86s/exit0, direct Unitful1.28/HDF50.17.4 Pkg.add, новый output path и narrow scope. Единственное handoff замечание из `final-root-pair.md` закрыто на уровне документации; runtime/source assertions повторять не требуется. Reviewer не менял handoff.
