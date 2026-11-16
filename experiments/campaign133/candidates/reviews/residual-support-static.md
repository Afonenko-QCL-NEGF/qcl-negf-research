# Static review: cont128 residual-support snapshot

Статус: **пригоден для bounded fixture и затем fixed-source diagnostic**,
при сохранении указанной области. Один крайний случай replay comparator
ниже не затрагивает фактический input map128. Ни причина stalled residual,
ни сходимость/accuracy решения этим review не установлены.

Reviewer прочитал только code и малые TOML; Julia/Pkg, HDF5 datasets,
fixture и science compute не выполнялись. Source не редактировался.
Применён qcl-physics-audit. HEAD superproject b515c42ef14c1f61151626280da0bcf100d116be;
live Core f50fa9b, Runner18f45ed, research69f1ef4; чужие изменения сохранены.
Источником diagnostic служит archived accepted Core e171, не live Core.

## Frozen provenance и source relation

| Artifact | SHA256 |
|---|---|
| probe.jl,118 lines |03a1fdd1a43023091fe067fbfb9370f7849b752c696e4a02ce7cc87192512148|
| input.toml |46a74dd52a61b81de7d97c994f6ecb76ae2916d07a07208d3a8d31096c5fcc18|
| protocol.md |46e8369cc03bd5f00615747769de7c6e132d8b1a733f95ba92fbcbf38715edfa|
| archived gamma_functions.jl |57f2c277255ae1f8265659e65176cc88f97d5a40cbf8d6e32b7ea9257ccfb032|

Первые три hashes независимо пересчитаны и совпали с freeze. Малые
source map/input/identity/writer hashes совпали с input. Все шесть paths
ведут в cont128/attempt01; код проверяет realpath ownership, exact file
SHA, identity input/writer/Core/tree/clean и map call128. Фиксированы
Core e171caa60cce676c6c4238dd671e7c527f2bd977 и tree742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8.
Большие terminal/model files reviewer не открывал и их hashes не
пересчитывал: указанные fb7cea… и47682a… совпали с root-owned
affine4096/pilot32.toml parent_terminal/target_model. Probe проверит эти
полные hashes до datasets; это остаётся будущим actual input check.

Сохранённый writer `gamma_functions.jl:170` определяет raw residual через
`norms(rawphi,gamma,we,wk)` до Hermitian working-map conversion.
`norms:106` суммирует |rawphi−gamma|², reference |rawphi|² и их wE·wk
варианты. `evaluate:599–633` сохраняет эти raw norms в map и те же input
Gamma/raw_Phi, fitted μ и matched Green в reduced terminal. Поэтому
snapshot relation не смешивает working residual или mixed AA step.
Native trigger в map128=false; source остаётся несошедшимся.

## Массивы, меры и units

`readrows` сохраняет Julia logical (E,k,a,b) axes; full shape guard
(NE,Nk,Nb,Nb) применяется ко всем real/imag Gamma/raw_Phi datasets.
Slabs1024E читают оба complex поля; `blockstats` использует complex abs2
и все orbital entries без transpose, потери imaginary part или удвоения.
Header/model-v2 grid_ε/grid_wᴱ/grid_κ/grid_wᵏ и terminal-v1 paths
соответствуют archived writer. Для HDF5.jl это logical axes; Python
physical reversed axes здесь не предполагаются и не угадываются.

`W=we*transpose(wk)` даёт положительную energy и неотрицательную k меру;
суммирование r²/s² по orbital block перед умножением на W математически
равно source full-field нормам. Отличается порядок floating reductions,
поэтому replay_relative1e−10 — arithmetic/source check, а не scientific
acceptance. Нет spin,2π, quadrature renormalization или occupation mask в
primary norms. Actual map scales17.311496519918023 и0.16279633580600295
конечны и положительны. Relative raw≈7.192305218634285e−6 и weighted≈
5.957548913257313e−6 не заменяются новой f² нормой.

E0 берётся из matched map (.1eV), L0=1e−8m проверяется через physical
kmax=.6/nm и model κmax; weights остаются dimensionless source measure.
E_eV=E·E0 относительно source Eref; κ dimensionless. Uniform spacing
проверяется по actual model E, не подставляется старый30721 grid.
Stored μ сверяется с map μ=.08459613668783181; T=70K и exact-SI kBT/E0
используются без нового fit. Никакие Φ, solver, regrid, inverse, clipping,
PSD repair или Hermitization input не вызываются.

## Зоны, f² и near-zero поведение

Зоны E<μ−5kBT, |E−μ|≤5kBT, E>μ+5kBT являются disjoint и покрывают finite
source window. Boundary equality относится к средней зоне. Доли — доли
полного **squared weighted numerator**, не доли linear norm. top20 тоже
выбирается по squared weighted numerator, не по локальному relative
error; perE/perK сохраняют также raw reference/numerator.

Stable Fermi branch вычисляет только exp неположительного аргумента.
Большие ±arguments, в том числе overflow division к ±Inf, дают корректные
лимиты0/1; underflow exp к0 не создаёт overflow или 0/0. kBT проверяется
на finite и >0 до measures. Для actual μ/E/T все thresholds конечны.
Broadcast R.*f.^2 применяет energy f к каждому k. Поэтому эта норма
измеряет matrix residual, умноженный на f, в той же мере; она не равна
ошибке density/current, collision balance или невязке original SCBA.
Малость f² diagnostic не разрешает игнорировать unoccupied residual.

ratio(0,0)=0 и ratio(nonzero,0)=Inf соответствуют source convention.
Для нулевого полного weighted numerator fractions=NaN и отдельный
zone_fraction_status="unavailable_zero_numerator"; нет подстановки нуля
в недоступные доли. Source-free fixture проверяет complex blocks,
неравномерные weights, зоны и f² отдельно; stated128/42,234/49 и
14/234,52/234,168/234 следуют из его явных2×2 inputs. Runtime fixture
по-прежнему not_measured.

**R1, low, fixed-source nonblocking:** replay comparator при valid
reference=Inf и value=Inf вычисляет abs(Inf−Inf)=NaN и ошибочно останавливает
replay. Helper ratio поддерживает nonzero/zero, но fixture проверяет
только helper, а не comparator. Actual map128 все шесть references finite
и positive, так что frozen input не затронут. Minimal general fix —
сначала `isequal(value,ref)`, затем finite relative arithmetic; либо явно
ограничить replay на finite source norms. Допуск менять не требуется.
Reviewer передал root замечание, source не правил.

## Selected annotations и допустимый вывод

traceΓ/Nb·E0 и ratio к h·E0 являются broadening proxy, не linewidth
certificate: отсутствуют spectral pole/derivative/overlap и boundary
analysis. Selected HermΓ min-eigenvalue и norm(Γ−Γ†) только измеряют input;
они не исправляют его. Anti-Hermitian Frobenius column остаётся в source
dimensionless self-energy units, min-eigenvalue и trace proxy — eV.
Condition и complex Atrace берутся с того же selected(e,k); при отсутствии
dataset записываются NaN и explicit availability flags. Atrace — scaled
Green units, не физическая DOS без остальных мер.

Можно будет утверждать, какие сохранённые E/k nodes и occupation zones
несут raw residual и коррелируют ли top nodes с узким proxy/conditioning.
Один snapshot не определяет Jacobian, null/slow mode, причинность stall,
discretization sufficiency, baseline acceptance или native SCBA convergence.
Корреляция ширины/condition с residual не является причинным доказательством.

## Prospective resources и execution boundary

Внутренний guard160s начинается после full-file hashes; elapsed SHA
сохраняется отдельно. Общий timeout180s должен быть внешним и охватывать
hash/JIT/всё выполнение. MemoryMax2GiB/Swap0,Julia1/BLAS1 и output32MiB
должны быть установлены root scope: input сам их не применяет. Slab memory
порядка единиц/десятков MiB плюс O(NE·Nk) scalar matrices; это статическая
оценка, не измеренный peak. Мелкие block Matrix allocations и full-file
SHA могут определять цену; runtime/RSS пока not_measured. Output cap
проверяется после небольших CSV, timeout проверяется между slabs и в
конце; это не штатная hard preemption отдельного I/O.

Нет новых execution полномочий от review. Root сначала выполняет
source-free fixture и отдельно выбирает ресурсное окно после активного
AA4096; только затем bounded frozen-source run даёт actual evidence.

## R1 delta closure: b77d freeze

**R1 закрыт. Новых блокеров нет; bounded actual fixed-source diagnostic
может выполняться root в заявленном2GiB/180s окне.** Это review diagnostic,
не научная приёмка source или утверждение причины stalled convergence.

Проверены SHA main b77dac503c3bf35e9bc7bc927242fbcbdeeafe7dec608e91ab5e0c05fb81edee,
R1-comparator.diff1212267c716b2df295b35b0823c5743be192db4cbce9a62a83121804eff29ab7,
unchanged input46a74dd52a61b81de7d97c994f6ecb76ae2916d07a07208d3a8d31096c5fcc18
и protocol df77034c37ae3f84b2e7b6a553682bc34a54a302f54ef7a649002bd3c9810ec3.
Archived V1 сохранил exact03a1 hash. Independent text reversal только
helper/8 assertions/callsite дал побайтовое равенство V1; source I/O,
численные residual expressions, weights и физическая область не менялись.

`replay_ok` сначала отвергает NaN и invalid tolerance; затем принимает
точное равенство, в том числе одинаковые signed Inf и ±zero. Для
неравных значений разрешена только finite relative arithmetic с прежним
tol·abs(reference); finite/nonfinite mismatch, arithmetic overflow delta
или limit отвергаются. Reference=0 по-прежнему требует точный ноль.
Это устраняет Inf−Inf false mismatch, не ослабляя finite map128 tolerance.

Прочитан root-owned fixture-v2-root.log:12/12 PASS, exit0, process1.73s,
peak387532KiB, swaps0. Это actual source-free execution evidence от root,
а не собственный запуск reviewer и не actual HDF5/source-map replay.
Review не читало source arrays; все прежние scientific limitations
сохраняются до и после будущего bounded diagnostic.
