# Паспорт происхождения reference-2019 / S01

Снимок проверки источников: 2026-10-19. Это паспорт declared inputs и их
происхождения; numerical acceptance и experimental validation не выполнялись.
Сохраняются 0 mV/period, lattice/LO 70 K, LO/acoustic/impurity/IFR/alloy и
electron–electron mode `none`. Значения модели, конфигурации и допуски
этим документом не изменяются.

## Неизменяемые источники и смысл identity


Снимок чтения: root `e0f244d075dd49874ba4f3847f4efa693f3b3622`,
Research `5458212063a41051328ef8b5d14d63d5ac7dc781`,
Core `4af6e94808c6e60affbf1050a5184dbc12083335`,
Runner `0c53c6d0f68e0b25fce79388db831da6fc3ae856`.
Core `4af6e948…` меняет только три документа относительно numerical parent
`3108150dd1b6b1377470b384b15c26e272fe5c5d`; рассматриваемые formulas/functions
прочитаны в самом `4af6e948…`. Это evidence locators, не дополнительный lock-файл:
состав опубликованного release по-прежнему задают gitlinks superproject.

Версии runtime/imported producer здесь **unknown: выполнение не создано**.
Frozen plan bytes/hash/fingerprint, configuration/problem/physics fingerprints,
run/WorkChain/CalcJob UUID, execution/point/attempt/state identity отсутствуют:
**run = not_created**, остальные поля `null` с причиной «source-only design».
Паспорт не подменяет будущий resolved plan, runtime identity или artifact receipt.

Адресные исходники (пути далее относительно owning repository):

| Репозиторий | Путь / symbol | Роль свидетельства |
| --- | --- | --- |
| Research | `literature/reference-2019/{source.yaml,structure.csv}` | declared publication/data identity, исходные units и атрибуция |
| Research | `config/equilibrium.yaml`, `config/studies/equilibrium-0-70.yaml` | выбранная ordinary S01, axes/policies/output и explicit overrides |
| Research | `config/model/reference-2019-70k.yaml` | явные material/scattering/numerical inputs до override |
| Research | `config/policies/analytical-diagnostic.yaml` | численная политика и явные alloy values/channel |
| Core | `src/domain/{types,presets,scaling,electron_electron}.jl` | typed physical inputs/default closures/scales и преобразования |
| Core | `src/numerics/reference/grids.jl::{build_grids,build_profiles}` | выбор слоёв на центрах ячеек и дискретная donor normalization |
| Core | `src/numerics/reference/{hamiltonian,basis,kernels}.jl` | конкретные операторы, использующие числа |
| Runner | `src/infrastructure/config/configuration.jl::{_deep_merge!,_parse_physical,_parse_physical_models,_parse_domain_adaptation}` | порядок присваивания, decoding и typed defaults |
| Runner | `src/infrastructure/scientific/definitions.jl::resolve_scientific_configuration` | study/variant/inclusion overrides и generated assignments |
| Runner | `src/application/scientific/planning.jl::{_point_resolved_raw,_point_passport,scientific_plan_dict}` | effective point axes, hashes и сохранённая provenance |

## Проверенный первичный источник

Авторский текст: L. Bosco, M. Franckié, G. Scalari, M. Beck, A. Wacker, J. Faist,
*Thermoelectrically cooled THz quantum cascade laser operating up to 210 K*.
Идентичность: [arXiv:1911.06582v1](https://arxiv.org/abs/1911.06582v1),
[DOI 10.1063/1.5110305](https://doi.org/10.1063/1.5110305),
*Applied Physics Letters* 115, 010601 (2019). Title/authors/version прочитаны
на первой странице PDF; связка DOI/journal сверена в primary audit по actual
arXiv abstract metadata. Прямое чтение publisher version через DOI было недоступно;
это не утверждение о тождестве всех версий статьи.

Проверенные bytes: [author PDF v1](https://arxiv.org/pdf/1911.06582v1),
2 760 444 bytes, SHA256
`b103b3da9fa977d1dc682d1e2caf96d74b5f081731e4db0da944af675a58dc70`.
Это точное совпадение с `source.yaml::source_artifacts[1911.06582.pdf]`.
PDF имеет 5 страниц; relevant страницы 3 и 4 просмотрены полностью визуально,
а стр. 1–4 адресно извлечены в текст для поиска локаторов. PDF, уже полученный
координатором, повторно не загружался. Source tarball и отдельная figure PDF,
перечисленные в `source.yaml`, данным task не получались и их bytes не проверены.

| Величина / наблюдение | Primary locator в PDF с указанным hash | Что устанавливается |
| --- | --- | --- |
| Последовательность `32.6/79.9/19.0/84.0/29.0/51.6 Å` | printed p. 3, левая колонка, первый абзац; p. 4, Fig. 3(a) caption | слой 29.0 Å подчёркнут; 32.6 и 19.0 Å выделены как barriers |
| GaAs/Al₀.₂₅Ga₀.₇₅As | printed p. 4, начало Fig. 3(a) caption | composition/geometrical attribution, не набор всех material constants |
| `4.5 × 10¹⁰ cm⁻²` | printed p. 3, левая колонка, первый абзац, строка о подчёркнутом well section | primary main-text sheet quantity |
| `4.5 × 10¹⁰ cm⁻³` | printed p. 4, Fig. 3(a) caption, первая/вторая строки | внутренняя dimensional discrepancy того же авторского PDF |
| `56 mV/period`, `200 K` | printed p. 4, Fig. 3(a,b) caption; bias также p. 3, правая колонка, первый абзац | условия показанного finite-bias расчёта, отличные от S01 |
| `16.5 meV`, extraction energy около `41 meV`, `fosc = 0.20` | p. 3, первые абзацы колонок; p. 4, Fig. 3(a) caption | опубликованные характеристики design; это не заданные acceptance targets S01 |
| Обсуждение rudimentary e–e extension при 200 K | p. 3, правая колонка, первый абзац | дополнительное обсуждение статьи; не разрешение добавить e–e в S1 |

**Наблюдение:** конфликт cm⁻²/cm⁻³ подтверждён primary bytes и визуальными
локаторами; это уже не secondary report или failed screenshot.
**Вывод:** обычная dimensional conversion не устраняет конфликт текста и caption.
Не доказано, что именно авторы intended в discrepant caption; «опечатка» остаётся
интерпретацией, а не подтверждённым erratum. **Сохранённая assumption модели:**
используется sheet dose `4.5e14 m⁻²`, совместимая с main text p. 3, с явным caveat.
Никакого silent scientific выбора нового числа/размерности здесь нет.

## Полная цепочка дозы и нормировки

Печатная величина p. 3 интерпретируется как число доноров на площадь периода:

```math
N_{D,\mathrm{declared}}^{2D}=4.5\times10^{10}\;\mathrm{cm^{-2}}
=4.5\times10^{14}\;\mathrm{m^{-2}},\qquad
N_{D,+}^{2D}=f_{ion}N_{D,\mathrm{declared}}^{2D}.
```

Преобразование использует `1 cm² = 10⁻⁴ m²`; оно корректно **при sheet
interpretation**. Research model явно задаёт `donor_sheet_density: 4.5e14 m^-2`
и `ionization_fraction: 1.0`. Runner `_parse_physical` проверяет положительную
sheet dimension, переводит в `PhysicalParameters.N_dop²ᴰ` и `f_ion`; Core
`src/domain/presets.jl::ionized_sheet_density` возвращает `f_ion * N_dop²ᴰ`.
Полная ионизация — выбранная assumption, не установленная процитированным
печатным числом степень ионизации при 70 K.

Непрерывный reference профиль равномерного распределения в объявленном doped
segment `d_D = 2.90 nm` имел бы

```math
N_D^+(z)=N_{D,+}^{2D}\,\chi_D(z)/d_D,
\qquad N_D^+=1.5517241\times10^{23}\;\mathrm{m^{-3}}
\quad(f_{ion}=1).
```

Это dimensional reference, **не фактическое дискретное значение**.
Равномерное распределение по помеченному segment и выбранная ionization —
предпосылки модели; печатная доза сама их не выводит. Границы периода из сумм
слоёв: `0/3.26/11.25/13.15/21.55/24.45/29.61 nm`. Доза относится к
`[21.55,24.45) nm`. IFR interfaces модели — `0/3.26/11.25/13.15 nm`;
края doped segment внутри GaAs не добавляются как material interfaces.

Core `grids.jl::build_grids` задаёт `x_g=(g−1/2)Δx`, `Δx=L_p/(L_0 N_z)` и
`w_g^x=Δx`. `build_profiles` назначает слои по half-open rule (совпавшая
граница принадлежит слою справа), затем берёт `χ_{D,g}=1` только в `doped=true`.
`src/domain/scaling.jl::scaled_value` задаёт
`\bar N_{D,+}^{2D}=L_0² N_{D,+}^{2D}`, `\bar N_D=L_0³ N_D`.
Реализованная нормировка:

```math
d_{D,h}=L_0\sum_g w_g^x\chi_{D,g},\quad
\bar N_{D,g}=\frac{\bar N_{D,+}^{2D}\chi_{D,g}}
{\sum_{g'}w_{g'}^x\chi_{D,g'}},\quad
\sum_g w_g^x\bar N_{D,g}=\bar N_{D,+}^{2D},\quad
N_{D,g}=\bar N_{D,g}/L_0^3,\quad
\sum_g(L_0 w_g^x)N_{D,g}=N_{D,+}^{2D}\quad[\mathrm{m^{-2}}].
```

Нормировка имеет положительный denominator: `normalization > 0`, иначе Core
выбрасывает `ArgumentError("no doped grid cells were selected")`. Нулевая
target density не принимается конструктором `PhysicalParameters`; near-zero
electron observable здесь не используется. Тождество относится к donor input,
не к computed electron density, neutrality, FD/FDR или сходимости SCBA/Poisson.

Для **ожидаемого source candidate S01** `N_z=96`, `L_0=10 nm` независимая
арифметическая подстановка в эту формулу даёт:

| Величина | Source-derived expected значение | Статус |
| --- | --- | --- |
| `L_p`, `Δz` | `29.61 nm`, `0.3084375 nm` | сумма объявленных толщин / деление, не measurement |
| Число ячеек шести слоёв | `[11,25,7,27,9,17]` | half-open cell-centre assignment |
| Doped indices Julia | `71:79`, 9 cells | source-derived ожидаемая выборка |
| `d_{D,h}` | `2.7759375 nm` | отличается от непрерывных `2.90 nm` |
| `\bar N_{D,+}^{2D}`, `\bar N_{D,g}` в doped cells | `0.045`, около `0.16210739615` | безразмерные sheet и volume representations |
| `N_{D,g}` в doped cells | около `1.6210739615e23 m⁻³` | discrete normalization, вне segment 0 |
| Интегрированная donor sheet density | `4.5e14 m⁻²` | алгебраическое тождество при `f_ion=1` |

Это **не запуск Julia, не готовый профиль из run и не test pass**. Проверка
исполненного профиля остаётся будущему bounded preflight/frozen-plan audit.
Core docs `02_reference2019.md` приводят другой `N_z=198` с indices `145:163`
и около `1.5837466e23 m⁻³`; их нельзя приписывать S01 `N_z=96`. Существующий
`test/physics/reference2019_units_geometry_and_quadratures.jl` содержит такой
198-node test и `rtol=2e-15` для sheet sum, но в SCI02 он не исполнялся.
Это test-source contract roundoff normalization, не порог физической приёмки.

Как dimensional negative control, **не альтернативный input**: буквальные
`4.5e10 cm⁻³` соответствуют `4.5e16 m⁻³`; их интеграл по `2.90 nm` был бы
`1.305e8 m⁻²`. Это иная доза, а не эквивалентная запись `4.5e14 m⁻²`.
Нельзя молча передать volume value через sheet field, поменять единицы или
объявить конфликт исправленным одним преобразованием.

Один и тот же donor profile входит в `impurity_kernel` (`weights_donor =
wˣ .* Nᴰ / L₀²`) и источник Poisson; тождество sheet sum не доказывает,
что электронный final state имеет правильную форму или принят по charge gate.

## Formula, value и assigning origins

**Formula origin** — конкретное математическое предположение/оператор,
**value origin** — первичное свидетельство или явно выбранное число,
**assigning origin** — последний source/override/default, который реально
присваивает leaf. Одно не заменяет другое. Ниже formula locators — читаемые
уравнения и symbol текущего Core, а не новый независимый аудит вывода всех
микроскопических формул. Первичная независимая проверка FD/FDR/error budget и
допусков относится к SCI03; generic bibliography не доказывает происхождение числа.

Все `M` ниже — `config/model/reference-2019-70k.yaml`, `P` —
`config/policies/analytical-diagnostic.yaml`, `S` —
`config/studies/equilibrium-0-70.yaml`. Значение в колонке effective — ожидаемое
после прочитанной композиции, **не наблюдение готового frozen plan/runtime**.

| Вход / expected effective S01 | Formula origin / использование | Value origin и предел attribution | Assigning origin |
| --- | --- | --- | --- |
| Layers `3.26/7.99/1.90/8.40/2.90/5.16 nm`, barriers Al₀.₂₅Ga₀.₇₅As / wells GaAs | Core `grids.jl::build_profiles`, `02_reference2019.md`; Å→nm деление на 10 | Primary p. 3 и Fig. 3(a), caveat дозы отдельно | M `physical.layers` (sequence заменяет sequence) |
| Doped layer 2.90 nm, sheet `4.5e14 m⁻²` | Цепочка нормировки выше | Primary p. 3 sheet; caption conflict сохранён | M `layers[5].doped`, `physical.donor_sheet_density` |
| Ionization `1.0` | `ionized_sheet_density = f_ion*N_dop²ᴰ` | Явная model assumption; полнота ионизации при 70 K первично не установлена | M `physical.ionization_fraction` |
| IFR interfaces `0/3.26/11.25/13.15 nm` | `kernels.jl::interface_roughness_kernel`, `_band_jump` | Derived from declared material boundaries; независимые rough interfaces — model assumption | M `physical.interfaces` |
| Barrier edge `0.20775 eV`, well edge `0 eV` | `hamiltonian.jl::build_bdd_hamiltonian`, basis projection; theory `05_basis.md` | Явные M/Core preset числа; primary value citation/temperature/composition interpolation не установлены. Не выводить автоматически из Al fraction | M layer `conduction_band_edge` |
| Longitudinal/transverse relative masses `0.08775` barrier, `0.067` well | BenDaniel–Duke и projected transverse kinetic operator, Core `05_basis.md`, `hamiltonian.jl`/`basis.jl::project_hamiltonians` | Явные model assumptions; первичный origin точных чисел не установлен | M per-layer `longitudinal_relative_mass`, `transverse_relative_mass` |
| Layer/Poisson `εr=12.9`, static `εs=12.9`, high-frequency `ε∞=10.89` | Core `10_poisson.md`; impurity `EQ-KIMP-001`; LO `EQ-KLO-001` | Явные model inputs; ни совпадение с preset, ни единый DOI не доказывают material dataset | M layer `relative_permittivity`; M physical permittivity leaves |
| `ħωLO=36.7 meV` | LO kernel и сдвиги `EQ-KLO-001`, `EQ-CONTRACT-LO-001` | Явное M/Core preset число; первичный value origin не установлен | M `physical.lo_phonon_energy` |
| Impurity/LO screening `0.20/0.20 nm⁻¹` | `impurity_kernel`: `q+qs`; `lo_phonon_kernel`: `q²+qz²+qLO,s²` | Фиксированные interaction assumptions, не численные регуляторы для pass; value citation не установлена | M two screening leaves; typed `screening=fixed` |
| IFR height `0.10 nm`, correlation length `9.0 nm` | Gaussian `S(q)=πΔ²Λ² exp(−q²Λ²/4)`, `EQ-KIFR-001`, `interface_roughness_kernel` | Явные model assumptions; измерение именно этой структуры не предъявлено | M roughness leaves |
| Acoustic `Ξ=7.0 eV`, `ρ=5317 kg/m³`, `vs=5240 m/s` | Elastic equipartition `Ξ² kB TL/(ρvs²)`, `EQ-KAC-001`, `acoustic_kernel` | Явные model constants; primary value origins и error от elastic approximation при 70 K не установлены | M acoustic/mass_density/sound_velocity leaves |
| Alloy `ΔV=0.60 eV`, `Ω0=4.5168e−29 m³` | `Ω0 ΔV² x(1−x)` integral, `EQ-KALLOY-001`, `alloy_kernel` | Явные P assumptions; первичный origin/primitive-cell convention числа не установлен | P overrides M `null/null`; P enables alloy |
| `TL=TLO=70 K` | Acoustic prefactor, thermal Bose `NLO=[exp(ħωLO/kB TLO)−1]⁻¹`; Core `08_kernels.md` | Явная исследовательская постановка S01, не 200 K Fig. 3 validation | M 70 K → S axes → generated temperature-axis → point resolution |
| `Vperiod=0 mV`, `F=0 V/m` | `F=Vperiod/Lp`, `UF=−eF(z−z0)`, theory `01_model.md` | Явный выбор S01; primary 56 mV относится к иной точке | M 56 mV → S physical override 0 → point voltage 0 |
| `z0=0 nm`, `Eref=0 eV` | Coordinate/energy origin and scaling, Core `01_model.md`, `scaling.jl` | Выбор нулей, не измеренные числа | M position/reference leaves |
| Spin degeneracy `2` | Electron density integration `observables.jl::_electron_density_bar`, Core `11_observables.md` | Выбранная degeneracy convention; не удваивать donor sheet при conversion | M `physical.spin_degeneracy` |
| Scales `E0=0.1 eV`, `L0=10 nm`, reference mass `0.067m0` | `scaling.jl::ScaleSystem/scaled_value/physical_value` | Выбор единиц; reference mass не ещё одна fitted layer mass | M `scales` |

### Typed defaults и closure identity

В M/P/S нет `physical_models` или `domain_adaptation`. Runner
`configuration.jl::_parse_physical_models/_parse_domain_adaptation` явно создаёт
соответствующие defaults; Core `types.jl::PhysicalModelOptions` и
`electron_electron.jl::ElectronElectronOptions` подтверждают baseline.
Отсутствие YAML leaf здесь **не равно** неизвестному отсутствию physics.

| Closure / expected mode | Assigning origin | Применимость и ограничение |
| --- | --- | --- |
| `screening=fixed` | Runner parser default → Core typed options | Два заданных screening wavenumbers сохраняются; density-temperature update не активирован |
| `dispersion=parabolic` | Runner parser default → Core typed options | Однозонная effective-mass модель; Kane не активирован |
| `lo_population=thermal` | Runner parser default → Core typed options | Thermal LO bath при `TLO=70 K`; fixed/hot population modes не активированы |
| `electron_electron.mode=none` | Runner default `get(ee,"mode","none")` → Core `ElectronElectronOptions(mode=:none)` | Сохраняет явное решение оператора S1. Нет e–e collision/exchange self-energy; Hartree/Poisson остаётся |
| `domain_adaptation.mode=none` | Runner parser default | Нет автономного изменения energy window |

Другие нулевые numeric fields этих неактивных closure branches — typed sentinels,
а не измеренные `0 K`, нулевая carrier density или physical absence evidence.
`include_exchange=true` в options при `mode=none` не включает обменный канал;
mode отсутствия должен быть явным в future resolved/runtime identity. Нельзя
добавлять e–e потому, что статья обсуждает extension при 200 K. Перед любым
S1 submit effective mode `none` проверяется в frozen plan/runtime model identity;
missing или необъяснённый mode блокирует submit, а не считается pass.

### Политика, численные inputs и конечный point

| Группа / expected effective S01 | Назначение и assigning history | Предел свидетельства |
| --- | --- | --- |
| LO/acoustic/impurity/IFR enabled; alloy enabled | Первые четыре M=true; alloy M=false → P=true; Core `ScatteringOptions/build_kernels` | Пять каналов сохраняются. Enabled не доказывает корректность/разрешение kernel quadrature |
| `Nz=96`, `Nb=5`, `Pbasis=3` | Nz M198 → P96; Nb/Pbasis M | Базисное окно `2P+1=7` periods, 35 selected states до центрального усечения; не путать с embedding4 |
| Energy `NE=30721`, `[−0.30,0.90] eV`, margin `80 meV`, seed `1.25 meV` | P overrides NE2401/window/seed0.5; margin M | Seed не физическое уширение, `ηwork=0`; окна и их tails не сертифицированы |
| `Nk=11`, `kmax=0.60 nm⁻¹`, angular32, `Nqz=81`, qzmax10 nm⁻¹ | Node counts M61/64/401 → P11/32/81; cutoffs M | Сохранённые quadratures/cutoffs, не доказанная достаточность |
| sparse_plan / FFT / paired_convex / finite_chain4 | energy_shift M sparse_plan и P same; hilbert M direct → P fft; normalization/embedding P | Нормировка и boundary — объявленные operator choices. Embedding4→8 не «автоматический refinement» |
| localization `pzp`, full self-energy, resolved transverse momentum, Kramers–Kronig | M algorithms leaves | Ограниченный projected basis; нет переноса решения о достаточности из tutorial/source tests |
| dense_blas contraction, direct kernel_build, linear mixing, low_rank tolerance0/rank0 | M algorithms leaves | Source-level identity; auto_exact не меняет scientific approximations |
| auto_exact / production backend / BLAS1 | M execution; P повторяет auto_exact | Hardware chunk/worker resolution и фактические threads требуют fresh runtime/admission |
| SCBA cap2048 / Poisson cap128 / mixing0.2/0.15 | M2000/100/0.1/0.15; P SCBA64; S2048/128/0.2; Poisson mixing M0.15 | Iteration caps не равны scientific acceptance и не являются compute grant |
| adaptive_working, stagnation window0/improvement0 | P convergence → S policy → generated convergence-policy; P stagnation overrides | Working return не strict final pass |
| Minimum SCBA/Poisson iterations2/2, consecutive passes3/2; прочие strict tolerances и validation windows | M `solver.convergence/tolerances/validation_windows`; соответствующие leaves P/S не меняют | Полные numeric thresholds остаются в owning config/Core; SCI02 не назначает/не ослабляет их и не выводит accuracy из их наличия |
| Kernel construction impurity tolerance1e−5, lookup65..513, relative2e−4, power2, validation fractions0.25/0.5/0.75, pairs3, strict=true | M `kernel_construction` | Эти controls не доказательство достаточности angular/qz discretization всех каналов |
| Checkpoint every SCBA32, outer0; phase timing true | M SCBA0/outer0/timingfalse → P timingtrue → S SCBA32 | Persistence/telemetry policy, не convergence observable |
| Archive full_final/projections true, optical false, intermediate history64 | S scientific outputs → generated aliases | M optical-study true и light snapshots64/затем P16 не final output policy; generated study отключает optical, output sets history64 |
| Recovery enabled / interval1800 s / retain2 / budget8 GiB / reserve64 MiB | S outputs.recovery | Не compute grant и не доказательство достаточного storage; SCI07/08 отдельно |
| Telemetry enabled / buffer256 | S outputs.telemetry | Не measured RSS/I/O или science gate |
| Independent single point; one baseline variant, operation stationary, purpose research | `config/equilibrium.yaml` includes only S; baseline has no variant overrides; S axes0 mV/70 K | Generated labels — административная identity. 48 mV bundle не включён |

M также содержит output device geometry `405 periods /150 μm /1.8 mm` и photon
grid `1..80 meV/321`. Они не становятся physical S01 acceptance targets:
scientific output/`generated:single-study` задают optical=false, а device periods
не входят в `PhysicalParameters` периода и не умножают donor sheet. `source.yaml`
прямо отмечает, что 405 periods не установлены публикацией. Native output
identity и перечень сохранённых arrays всё равно требуют будущей проверки.

Порядок присваивания читается в Runner `resolve_scientific_configuration`:
ordered M → P, затем study overrides → variant overrides → inherited inclusion
overrides; затем generated run identity/convergence/temperature/output/single-study
assignments. Mapping merge recursive; sequences заменяются; histories leaf paths
дополняются `_deep_merge!`. В выбранной baseline variant overrides пусты,
meta equilibrium не добавляет inclusion overrides. `_point_resolved_raw` вносит
point-specific voltage/TL/TLO и study axes, `_point_passport` хэширует эти bytes.

Поэтому ожидаемые `0 mV/70 K` должны совпасть в physical override, generated
study axes и point passport. Raw M preset `56 mV` и Core preset `200 K` —
история/default context, а не effective S01. `provenance.sources` показывает,
кто присвоил значение; ссылка на generated label не является литературным
обоснованием этого числа. `resolved_configuration` и typed effective closures
должны сверяться вместе: сохранённый raw dictionary сам не перечисляет все defaults.

## Вердикты и пробелы

Уровни свидетельств разделены: публикационная атрибуция геометрии установлена
для проверенной версии PDF; donor sheet interpretation сохранена с caveat;
equation convergence, physical gates, discretization adequacy и experimental
validation S01 здесь **unknown / не измерялись**. Нет final U/G, attempt или
численного результата, из которого можно было бы вывести FD/FDR, J=0,
нейтральность или uncertainty при 70 K. Показанная в статье finite-bias NEGF
кривая200 K и её digitization не являются таким свидетельством.

| Утверждение | Категория / свидетельство | Ограничение / дешёвая различающая проверка | Решение |
| --- | --- | --- | --- |
| Publication identity и geometry соответствуют declared reference | Наблюдение: проверенный PDF hash; p. 1 title/authors, p. 3 sequence, Fig. 3(a) | Publisher bytes не прочитаны; source tarball/figure PDF не проверялись | сохранить |
| Main text содержит sheet dose, caption содержит volume unit | Наблюдение: тот же PDF, p. 3 left first paragraph vs p. 4 Fig. 3(a) caption | Intended caption correction/author erratum неизвестны; при необходимости проверить official erratum/author clarification, не менять модель | сохранить |
| `4.5e14 m⁻²` сохраняет выбранную sheet interpretation | Вывод/assumption: p. 3 + exact cm⁻²→m⁻² conversion + M | Не окончательное разрешение внутренней ошибки источника; иной dose choice — отдельное научное решение | сохранить |
| Donor profile нормирован на `fion*sheet` на реальной z-сетке | Source-derived вывод: Core `build_grids/build_profiles`, `scaling.jl`; алгебра выше | Runtime/profile not_created. Проверить bounded profile slice/weights/source identity в готовом preflight; не запускать solver ради паспорта | сохранить |
| Все material/scattering constants измерены/дословно взяты из Bosco2019 | Такой evidence не предъявлен; M/P фиксируют assumptions, formula locators известны | Первоисточник точных чисел/температурных law/uncertainties неизвестен; отдельный адресный literature audit, если требуется сильная attribution | данных недостаточно |
| S01 сохраняет пять каналов и e–e none | Source-derived expected composition M/P/S + Runner/Core typed defaults + прямой S1 выбор | Exact frozen/runtime identity ещё не создана; проверить до submit | сохранить |
| Published200 K curve валидирует equilibrium0 mV/70 K | Данных такого вида нет; несовпадающие conditions | Нужен отдельный model-verification evidence и sensitivities; новая finite-bias/experimental study вне S1 | данных недостаточно |
| Наличие YAML/normalization/source test принимает stationary solution | Источник не даёт такого следствия | Все четыре уровня: convergence/physical/discretization/experiment неизвестны; SCI03+future attempts отдельно | данных недостаточно |
