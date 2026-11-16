# Независимый результат: AAtrialPSD projected64 против same-seed control

**Вердикт:** actual projected64 завершил64 начатых и64 завершённых reduced Π
без domain failure. Конечная **исходная**, не projected/mixed, невязка
меньше same-seed stableA-cont64 control на39,9354% (weighted40,1581%).
Это принимается как ограниченный результат данного опыта с одинаковым
лимитом64 начатых Π. Решение SCBA, native acceptance, превосходство алгоритма
вообще или достаточность сетки **не установлены**. Scientific blocker для
этого узкого сравнения не найден. Exact PSD всех вычисленных proposals также
не устанавливается: измерены eigensolver/reconstruction errors и прохождение
существующего cone guard, а не знаки точных dyadic eigenvalues.

## Первичные свидетельства и воспроизведение

Источники: `equilibrium-gamma-projected64/attempt01/` и
`equilibrium-gamma-stableA-cont64/attempt01/`, их archived input/identity/
operator/stable_spectral, map/attempt/projection TOML, boundary replay,
summary, actual-one-map, log и небольшие HDF selections. Авторский диагноз
для вывода не использовался. Raw Φ и полные Γ/Green arrays повторно не считались.

Проверка: `python3 .build/campaign133/reviews/projected64_audit.py`.
Сохранены `projected64-audit.json`, `projected64-{control,projected}-maps.csv`,
`projected64-h5-selected.txt` с точными h5dump командами/выборками.
Фактическая цена окончательного аудита `/usr/bin/time`:0,48s,21944KiB peak, exit0;
Julia/Pkg/новые карты не запускались.

Все archived source SHA и input SHA совпадают с identity обоих runs.
Core e171caa60cce676c6c4238dd671e7c527f2bd977, tree742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8,
Julia1.13.3, Manifest da0780fb… одинаковы.
Projected input SHA2998d854ff2e734b2a5357f3f6e26f888c9b5fe02d46027d14806d5090596a06;
projected `gamma_functions.jl`1c1eb3c3…,
`trial_projection.jl`cf9dc3c1…;
stable_spectral eeca1f29… одинаков с control.
Полные SHA малых первичных входов/summary/map64 перечислены в audit JSON.

## Сопоставимость и остаток исходного уравнения

Оба run импортируют тот же GammaLast16 SHAab44aa53… и тот же model/physical
input, запускают emptyAA/freshμ; history3, affineL1≤4096, condition≤1e8,
linearα0,2, AAβ1, thermal70K, fixed U0 и пять каналов одинаковы.
Все input deltas: identity, `project_AA_trial=false→true`, pilot600→900s,
overall1800→2100s. Обе программы дошли до call cap, поэтому разные wall caps
не объясняют отсечение одной траектории. Это одинаковый **Π call budget**,
не одинаковая стоимость или число завершённых карт.

Обе boundary replay: whole raw Φ difference0, все9 проверяемых scalar fields
bitwise equal parent map16; новый call1 входит в64TOTAL.
Между runs карты1–4 имеют побитовое числовое равенство всех записанных
residual/working-residual/μ/input-cone/rawΦ-cone/channel-norm/whole spectral
ledger полей. У projection1/2 перед call3/4 correction0 и возвращён исходный
объект proposal. Первое изменение — projection3 перед call5:22 блока,
absolute1,22147e−32, relative7,05579e−34. Различие траекторий начинается
именно с call5. Это подчёркивает чувствительность дальнейшего AA к очень
малым исходным возмущениям; не доказывает, что вся конечная разница вызвана
устранением одного крупного cone obstruction.

В archived `gamma_functions.jl::norms/gamma_residuals/_evaluate`:

`r=||raw Φ(Γ)−Γ||F / ||raw Φ(Γ)||F`;
weighted аналог использует `Σ_E,k wE·wk·Σ_a,b |…|²` под корнем.
Проверяется полная сетка30704×11 и все5×5 элементы, без occupied mask.
Γ и Φ заданы в self-energy units E0=0,1eV; unweighted absolute нормы можно
перевести в eV умножением на E0. Это агрегированные нормы полей, не ошибка
одного матричного элемента. Нормировка без additive floor:0/0→0,
nonzero/0→Inf. AA использует `Herm(raw Φ)−Γ`, но это **другая** working
норма и она не подменяет критерий исходного уравнения.

| Конечная карта64 | Control | Projected |
|---|---:|---:|
| raw relative |4,125583337831235e−6|2,4780169900783377e−6|
| weighted raw relative |3,412502855963494e−6|2,042106454929957e−6|
| raw absolute (scaled E0) |7,142007180515379e−5|4,2898225631428256e−5|
| weighted raw absolute |5,555432467297159e−7|3,324476925370057e−7|
| μ (scaled E0) |0,08459614262982018|0,08459614535936885|

У projected и control лучшие записанные raw residual достигнуты на карте64.
Требование handoff — **обе** relative нормы<1e−9 три карты подряд — не
менялось. В projected последние нормы превышают этот порог в2478 и2042раз;
ни одна записанная карта не вызвала trigger. Number fit projected
relative error4,63e−15 при target0,045 — выполнение density constraint,
не решение остальных нелинейных уравнений.

## Счётчики, проекции и точность утверждения PSD

Projected:64 attempt records = seed1+linear1+AAII62, все accepted;
64 map records;0 domain-failure records;0 fallback;62 projection ledgers.
Control:65 attempt records, из них62 accepted = seed1+linear4+AAII54+
convex_fallback3; два q-domain failure (calls5,23) потребляют Π call, один precone
rejection не потребляет Π. Поэтому64 started Π у control дали только62
завершённых map records. Нельзя назвать это сравнением64 успешных updates
в каждом алгоритме или63 одинаково определённых успешных новых шагов.
У projected63 новых завершённых карт после общей boundary; у control61.

Из62 projection attempts12 действительно меняют proposal (next calls
5,21–25,38–43);50 оставляют exact-Hermitian, numerically nonnegative block
побитово неизменным. Every ledger coverage/rank-count =337744 blocks;
inline projection metadata совпадает с отдельными projection TOML.
Проекция применяется **только** к AAII proposal до fresh Π. В raw Φ,
local Γ components, spectral map или критерии проектора нет.
После freshΦ выполняется обычный cone guard; затем Herm representation
для AA, без eigenvalue clipping rawΦ.

Максимальная field correction на call23:
absolute7,136751835583073e−7; **relative4,122547515804853e−8**;
weighted absolute7,366941588181813e−9. Число4,12e−8 нельзя обозначать
absolute delta. Исходные proposals exact-Hermitian (max antiHerm norm0),
так что изменение здесь связано с eigenspectrum/reconstruction, не
устранением исходной antiHerm части. Max relative eigen backward residual
1,98487e−15, orthogonality defect4,36402e−15; суммарно projection18,862s.

Замена negative eigenvalues на0 не создаёт positive floor. Однако Float
reconstruction не гарантирует exact PSD: minimum computed input eigenvalue
среди projected maps −2,06875e−23; существующий absolute backward allowance
допускает это. Поэтому здесь принимается guard-admissibility, не global
strictPD. Отдельный exactPD audit **unprojected Γ-last64** не переносится
на projected Γ или на его промежуточные proposals.

## Terminal alignment и остальные границы

Map64 точно совпадает с `summary.reduced.last_map`, standalone spectral
ledger — с ledger внутри map64. Code path записывает reduced-terminal
Γ/Green/rawΦ из одного fresh evaluate64, затем accepted Γ становится
`Gamma-last`; best score также достигается на64.
HDF Gamma-last имеет нужную identity и pi_calls_total64; terminal flags
native_converged=false/handoff_trigger=false. Terminal μ совпадает с
map64 побитово при17-digit h5dump. Дополнительно **только selected** Γ
блоки (E382/k11,E5090/k1,E30704/k1) одинаковы побитово во всех трёх
Gamma-last/Gamma-best/reduced-terminal real+imag payloads. Полного
сравнения массивов этот audit не выполнял; alignment обоснован также
прямым writer/control-flow на frozen source, не выдан за новый whole-field
hash/equality test.

Начальный actual-one-map passed,worst5,64977e−16 — matched hH-all9 native
consumer control и начальная KMS проверка при currentμ. Это не исходный
raw-h native solve и не uniform KMS certificate всех последующих64 карт.
Summary native started/completed bounds0; все native accounting events0;
`not_measured=true`, `scientific_acceptance=false` корректны. Execution
complete/exit0 подтверждается log:547,916s harness/557,76s process,
peak7141588KiB,swap0; control519,697s/528,69s,peak7616820KiB. Не заявляется
ускорение: projected получил меньшую невязку при несколько большей wall
стоимости и отличающейся успешности map attempts.

Физически это уменьшение fixed-Hartree reduced Γ-map defect выбранного
intended hH(all9) оператора. Poisson не решался; original raw-h native
SCBA не запускался; точность real-axis quadrature/Γ-continuum, fullΠ
forward error, energy/basis/window convergence и эксперимент не измерены.
Следующая полезная проверка — один fresh unprojected rawΦ из сохранённого
terminal Γ (boundary replay) перед продолжением/изменением алгоритма;
projection trace остаётся отдельным алгоритмическим свидетельством,
а прежние native/physical gates сохраняются для любого handoff.
