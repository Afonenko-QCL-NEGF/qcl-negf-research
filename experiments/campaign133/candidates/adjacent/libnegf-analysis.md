# Стационарные схемы libNEGF и DFTB+: переносимость на S01

Дата: 2026-11-15. Пакет: read-only разбор открытых исходников; сборок, установок,
Julia, новых физических расчётов, изменений Core, push и серверных действий не было.

## Алгоритм и решение в первую очередь

libNEGF HEAD действительно реализует стационарный optical-phonon SCBA, включая
обновление retarded и electron-correlation self-energy. Его главная полезная
структура: **неупругий глобальный цикл по E,k с внутренним упругим SCBA в каждой
точке**, блоковое исключение principal layers и заранее вычисляемое embedding
электродов. Это не готовый oracle для нашей QCL: другая геометрия, другие ядра,
другая retarded closure и иной критерий завершения.

**DFTB+ не подтверждает полный optical-SCBA libNEGF HEAD.** Проверенный DFTB+
использует gitlink libNEGF 2023 года; этот tree содержит только упругие
local/block/overlap-mask dephasing-модели. Его SCC/Poisson+Broyden/Anderson —
внешний зарядовый цикл, не ускоритель полного phonon-SCBA. Для S01 с фиксированным
U0 этот внешний цикл вообще отсутствует.

Рекомендую использовать три проверяемые идеи: (1) точная KMS-редукция на
комменсурной LO-сетке как способ убрать медленное кинетическое направление,
(2) ограниченное вложенное решение упругого подблока как preconditioner полного
native map, (3) отдельная frozen-Γ проверка краёв конечного Hilbert-окна до
удлинения SCBA. Первая уже соответствует плану ведущего; ниже дан независимый
дискретный вывод. Не переносить current-only остановку или контурное
игнорирование inelastic Σr в научную приёмку S01.

## Фиксированные источники и лицензии

| Источник | Точная ревизия / происхождение | Лицензия и локальный checkout |
|---|---|---|
| [Официальный libNEGF](https://github.com/libnegf/libnegf) | `3b4120f65b19107149a46034ecd268be00622042`, 2026-11-12, shallow master | [LGPL-3.0-or-later](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/LICENSE); `src-libnegf/`, около 3.1 MiB |
| [Официальный DFTB+](https://github.com/dftbplus/dftbplus) | `f866d98910798356b738bcf2dde92dcd81668e25`, 2026-11-13, sparse shallow main | [LGPL-3.0-or-later, docs CC-BY-SA-4.0](https://github.com/dftbplus/dftbplus/blob/f866d98910798356b738bcf2dde92dcd81668e25/LICENSE); `src-dftbplus/`, около 26 MiB |
| libNEGF, выбранный самим DFTB+ | [gitlink `34db682aa7135e4ac98ddd25013bfbf0784d7ed2`](https://github.com/dftbplus/dftbplus/tree/f866d98910798356b738bcf2dde92dcd81668e25/external/libnegf), 2023-05-31 | Отдельно shallow-fetched в object database `src-libnegf`, checkout HEAD не менялся |
| Наш Core, узкое сопоставление | `e171caa60cce676c6c4238dd671e7c527f2bd977`, `.build/campaign133/core-accepted-science` | Локальный frozen accepted science tree; не изменён |

Суммарно новые два checkout значительно меньше лимита 300 MiB. DFTB+ submodules
не инициализировались; provenance архива и исходные gitlinks рабочего superproject
не менялись. AGENTS.md в деревьях этих двух внешних репозиториев не найден.
Первоначально зафиксирован superproject `b515c42ef14c1f61151626280da0bcf100d116be`
с уже существовавшими изменениями gitlinks Core/Runner и untracked research;
они оставлены как есть.

Все нижеуказанные строки относятся к **локальному сырому исходнику**, до fypp.
Браузерный raw-render некоторых Fortran файлов показал другие номера строк;
точная ревизия + имя процедуры — основной locator.

## Уравнения и циклы, действительно присутствующие в коде

Обозначения libNEGF: `Gn=-iG<`, `Gp=iG>=A-Gn`, `A=i(Gr-Ga)`,
`Sigma_n=-iSigma<`; размерности Gr/Gn — 1/energy, self-energy — energy.
В DFTB+ базис неортогональный: `M(z)=zS-H-SigmaR`, в нашем Core S=I.

### Dyson, Keldysh и блоки

\[
G^r(z)=[zS-H-\Sigma_L^r-\Sigma_R^r-\Sigma_{sc}^r]^{-1},\qquad
G^n=G^r(\sum_\alpha f_\alpha\Gamma_\alpha+\Sigma^n_{sc})G^a.
\]

`calculate_gsmr_blocks` делает правое Schur-исключение
\(g_i^R=[M_{ii}-M_{i,i+1}g_{i+1}^RM_{i+1,i}]^{-1}\);
затем `calculate_Gr_tridiag_blocks` восстанавливает необходимые диагональные и
соседние блоки. `calculate_Gn_tridiag_blocks` переносит также injection с обеих
сторон, включая соседние блоки Sigma_n. Это не приближение полной инверсии,
**если реальный оператор уже block tridiagonal**. Ограничение взаимодействия
диагональными/соседними блоками — отдельное приближение, если исходная self-energy
имеет дальние блоки.

Source: [iterative_cpu.F90:51–224,230–372](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/iterative_cpu.F90#L51).
Наш [finite_chain_green:172](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/periodicity.jl#L172)
уже устраняет левые/правые полости, включая lesser/greater injection. Поэтому
«внедрить RGF» само по себе не является новым лекарством S01 Nb=5. Нужна точная
проверка, какая часть времени ещё расходуется на повторные одинаковые полости.

### Полубесконечный lead embedding

`surface_green` строит `Ao=zS00-H00`, `Bo=zS01-H01`, `Co=zS10-H10`, затем
`decimation` удваивает исключённую область. Контактная self-energy равна
\(\Sigma_c^r=(H_{dc}-zS_{dc})g_s^r(H_{cd}-zS_{cd})\), с согласованными знаками
двух coupling factors. Временная независимость lead Hamiltonian позволяет
кэшировать g_s между SCBA passes. Новый cache label включает фактические complex
E, k-vector, spin, не только индекс.

Source: [contselfenergy.F90:82–266,339–410,571](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/contselfenergy.F90#L82),
[SelfEnergy_dns:571](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/contselfenergy.F90#L571).
CPU decimation имеет лимит 300 и два последовательных малых `C1`; это критерий
исключённого coupling, а не тождество спектральной суммы. В нашей finite_chain
граница — конечная цепь с локальным scattering, который меняется на SCBA: кэш
полостей между разными scattering states нельзя переносить как lead cache.
Замена finite_chain(4+4) полубесконечным bulk lead — **новая boundary identity**.
Независимый scalar oracle для lead: `g(z)=(z-e0-sqrt((z-e0)^2-4t^2))/(2t^2)`,
с ветвью `g~1/z`, Im g<0 при Im z>0. Он различает bulk и surface Green functions
без сравнения двух реализаций decimation. Для конечной цепи oracle — полная
инверсия той же объявленной цепи, а не этот бесконечный lead.

### Подлинный phonon SCBA HEAD

Публичные setters `set_elph_polaroptical` / `set_elph_nonpolaroptical` создают
отдельные interactions с `wq`, T, screened Fröhlich либо deformation-potential
ядром и опцией `tridiag`.
Source: [lib_param.F90:467–522](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/lib_param.F90#L467).
Здесь `wq` употребляется как энергия \(\hbar\omega\), поскольку делится на dE.

Пусть \(\mathcal K\) означает суммирование по поперечному q с уже включёнными
весами. По фактическим coefficients `compute_Sigma_r/n` реализуют

\[
\Sigma^n(E)=\mathcal K\{N G^n(E-\omega)+(N+1)G^n(E+\omega)\},
\]
\[
\Sigma^r(E)=\mathcal K\{(N+1)G^r(E-\omega)+NG^r(E+\omega)
 +\tfrac{i}{2}[G^n(E-\omega)-G^n(E+\omega)]\}.
\]

Source: [elphinel.F90::compute_Sigma_r:1080–1210, compute_Sigma_n:1293–1384](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/elphinel.F90#L1080),
[общий shifted contraction selfenergy:1539–1754](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/elphinel.F90#L1539).
Это finite-temperature absorption/emission, не просто elastic broadening.
Но отдельного principal-value/Hilbert вклада от Gn в этом пути не обнаружено;
показанная формула не равна автоматически нашему
\(\Sigma^r=\mathcal H[\Gamma]-i\Gamma/2\) на конечном окне.
Утверждать, что чужая реализация принимает нашу полную retarded closure, нельзя.
В частности, если Re Gr=H[A] на полной прямой, показанной формуле при thermal
Gn=fA для Kramers–Kronig consistency не хватает
`H[K(-f_minus*A_minus + f_plus*A_plus)]`. Это алгебраический вывод из двух
формул, не утверждение о намерениях авторов; на конечном окне добавляются
собственные tail errors.

`real_axis_int_inel`:

1. Готовит K и равномерную E-сетку; начинает без inelastic Sigma на iteration 0.
2. По всем k,E вызывает `compute_Gn`, который внутри вызывает
   `calculate_elastic_scba`: обновления Gr и Gn для elastic channels в этой точке.
3. Сохраняет Gr/Gn, интегрирует плотность и layer current.
4. По накопленным Gr/Gn всего окна вычисляет обе inelastic Sigma.
5. Проверяет изменение интегральных layer currents, повторяет внешний цикл.

Source: [integrations.F90:1202–1397](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/integrations.F90#L1202),
[iterative.F90:292–370](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/iterative.F90#L292).
В просмотренном пути нет Broyden/Anderson/homotopy внутри phonon-SCBA.

### Энергетическая интеграция и KMS

В coherent/elastic density path libNEGF разделяет заполненную часть, вычисляемую
контуром Gr, и nonequilibrium real-axis correction. Контур содержит круговую
дугу, линию над real axis и явно учитываемые Fermi poles
\(z_m=\mu+i(2m+1)\pi k_BT\). Число Gaussian points пользователь задаёт вручную.
Source: [integrations.F90::contour_int_def:693–855](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/integrations.F90#L693),
[DFTB+ manual transport.tex:534–584](https://github.com/dftbplus/dftbplus/blob/f866d98910798356b738bcf2dde92dcd81668e25/doc/dftb%2B/manual/transport.tex#L534).

Для inelastic путь принципиально иной: real interval расширяется на
`n_kT*kBT+wq_max`; применяется uniform trapezoidal grid, а не Gaussian.
`iEshift=nint(wq/dE)` округляет LO-сдвиг; вне окна подставляется G=0.
Source: [integrations.F90::real_axis_int_def:1011–1096](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/integrations.F90#L1011),
[elphinel.F90:1101,1314,1608–1638](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/elphinel.F90#L1608).
Точная KMS-ковариантность сдвигов кодом не обеспечивается при нецелом wq/dE;
это следует из несовпадения Bose-фактора исходной ω с округлённой энергией.
Ошибка detailed-balance exponent ограничена `beta*dE/2`, даже до ошибки
решения SCBA. Комменсурность устраняет эту причину, но не ошибку окна.

При inelastic contour обрезается ниже `mu_min-n_kT*kBT-wq_max`, Fermi poles
убираются; сам комментарий говорит, что inelastic Sigma_r на контуре **не
вычисляется**. Source: [integrations.F90:715–723,920–925](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/integrations.F90#L920).
Это approximation глубоко заполненной области, а не универсальный analytic
continuation полного SCBA. У нашего semiconductor density небольшой occupied
weight; нельзя считать эту область безразличной по аналогии с молекулярным DFTB.

В `energy_mesh.F90` существуют midpoint-refinement data structures, но поиск
production calls нашёл `create_mesh` лишь в `test_mesh`; active density path
пользуется описанными фиксированными quadratures. Поэтому наличие этого файла
не является свидетельством рабочего adaptive energy SCBA.
Source: [energy_mesh.F90:69,159,278](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/energy_mesh.F90#L69).

### SCC, Poisson, mixing и сохранение

DFTB+ `negf_density` передаёт H,S, μ в `compute_density_dft`; затем внешняя SCC
обновляет заряды → Poisson potentials → Hamiltonian. Poisson charge/potential
hooks находятся в [scc.F90:474,521](https://github.com/dftbplus/dftbplus/blob/f866d98910798356b738bcf2dde92dcd81668e25/src/dftbp/dftb/scc.F90#L474).
На зарядовом цикле берётся **сырой** `qDiff=qOut-qIn`,
`max(abs(qDiff))<sccTol` до смешивания. Source:
[main.F90:4971–4991](https://github.com/dftbplus/dftbplus/blob/f866d98910798356b738bcf2dde92dcd81668e25/src/dftbp/dftbplus/main.F90#L4971).

Modified Broyden ссылается на Johnson PRB 38,12807 (1988);
Anderson — на Eyert JCP 124,271 (1996), history, initial mixing и diagonal
rescaling. Anderson использует real inner products для flattened Hermitian
матриц и добавляет смешанный residual к комбинации inputs. Source:
[broydenmixer.F90:11–14,204–342](https://github.com/dftbplus/dftbplus/blob/f866d98910798356b738bcf2dde92dcd81668e25/src/dftbp/mixer/broydenmixer.F90#L204),
[andersonmixer.F90:192–351](https://github.com/dftbplus/dftbplus/blob/f866d98910798356b738bcf2dde92dcd81668e25/src/dftbp/mixer/andersonmixer.F90#L192).
Сам алгоритм mixing переносим; его charge-space tuning не обоснован для пяти
QCL-механизмов, а его память при хранении каждого SCC step не масштабируется
автоматически на наши NE=30704/30721 arrays.

Контраст для научной приёмки: libNEGF `check_J_convergence` использует
`max(abs(Jold-J))`, а не `max(J)-min(J)` и не raw self-energy residual.
`check_Mat_convergence` использует абсолютное max изменения элементов Gn.
Source: [scba.F90:119–173](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/scba.F90#L119).
При zero-V common-bath equilibrium J≈0 на многих несамосогласованных states;
этот stopping rule особенно слаб для S01. Ток нельзя нормировать на почти ноль:
нужен заранее объявленный абсолютный scale и полный Dyson/Keldysh/Σ residual.
Одинаковость двух последовательных токов не доказывает ни conservation, ни KMS.

## Первичная опора и ограничения библиографического вывода

1. [Pecchia, Penazzi, Salvucci, Di Carlo, NJP 10,065022 (2008), DOI](https://doi.org/10.1088/1367-2630/10/6/065022)
   — первичная методическая статья, прямо указанная libNEGF. Реализация
   block-iterative density/open boundaries/Poisson связывается с ней в официальной
   документации. Publisher PDF в браузерном инструменте недоступен; в этом пакете
   статья **не использована как построчный oracle новой optical-SCBA формулы**.
2. [Официальный DFTB+ recipe](https://dftbplus-recipes.readthedocs.io/en/latest/transport/carbon2d-trans.html)
   независимо от README показывает реальный user workflow PL partition,
   real-space Poisson, Green density и charge Broyden. Это workflow coherent
   transport; не свидетельство всех новых inelastic возможностей HEAD.
3. [Conrads, Di Napoli, arXiv:2605.21334v2, §4.3](https://arxiv.org/html/2605.21334v2#S4.SS3)
   — первичный отчёт разработчиков о дефектах boundary handling libNEGF.
   Применимый вывод узок: нужен независимый surface/resolvent oracle; известность
   кода и прежние regression outputs не заменяют физическую проверку.
4. Ниже KMS identity выведено непосредственно из Fermi/Bose distributions;
   для scalar lead независимый oracle — аналитический quadratic resolvent, а
   не второй вызов libNEGF. Дополнительный inelastic repository не требовался:
   нужная реальная optical-SCBA capability уже установлена по открытому HEAD.

## Три различающие гипотезы для нашей S01

Исходная постановка из поручения: zero V, 70 K, common baths, LO плюс
acoustic-elastic/impurity/IFR/alloy, Nb=5 Nk=11 Nz=96 Pbasis=3,
finite_chain 4 клетки с каждой стороны, fixed U0. Старый NE=30721 map после
190 free updates имел rK≈8.24e-5, rΣ≈1.32e-4. Frozen-Γ переход nodal PV →
linear-hat PI дал rsum .004922→.000406 и N(μ fixed) +.434%.
Это сообщённые ведущим результаты, не заново воспроизведённые здесь.

### 1. **Медленный кинетический модуль можно удалить точным thermal Γ-map, но только дискретно KMS-совместимым**

При общей μ,T для всех baths:
\(G^n=fA\), \(G^p=(1-f)A\), \(\Sigma^n=f\Gamma\).
Из показанного phonon collision map следует
\[
\Gamma_{LO}(E)=\mathcal K\{[N+1-f(E-\omega)]A(E-\omega)
                        +[N+f(E+\omega)]A(E+\omega)\}.
\]
Тождества
\(Nf_-=f(E)(N+1-f_-)\),
\((N+1)f_+=f(E)(N+f_+)\)
доказываются подстановкой \(f=(1+e^{\beta(E-\mu)})^{-1}\),
\(N=(e^{\beta\omega}-1)^{-1}\). Для static channels Γ=K[A].
Округлённые или интерполированные сдвиги могут нарушить эти тождества:
`W*(f*A)` вообще не равно `(W*f)*(W*A)`.

На m=939 commensurate grid (план ведущего, NE=30704) LO shift целочисленный,
поэтому в interior эта причина исчезает с точностью roundoff. Γ-only unknown
Hermitian; real packing Nb² действительных степеней на E,k и real Anderson-II
разумны. Это **новая solver parameterization того же equilibrium fixed point**,
только если исходный finite-window map действительно сохраняет KMS. Изменение
PI относительно nodal PV — отдельная discretization identity.

Различающий опыт: взять один Γ-state; собрать thermal lesser/greater без
нормировочной проекции, выполнить **исходный native полный map**, разделить
KMS defect по LO/static/embedding и interior/edge shells шириной ω. Далее 2–3
free full-map updates без FDR-projection и с неизменным U0. Фальсификация:
raw rK/rΣ растут или нарушается FDR/number; тогда thermal solve не принимает S01.
Размер допуска — существующие native tolerances; для identity defect записать
absolute max плюс относительный norm с denominator `max(norm(Sigma), scale_floor)`
в native energy units. Не объявлять произвольный 1e-12 физическим gate.

Стоимость: 1 collision audit + 2–3 raw sweeps; cap 15 min / 4 sweeps,
checkpoint после каждого, без параллельного второго SCBA. Один полный complex
Nb²*Nk*NE field ≈128.8 MiB; Hermitian real packing ≈64.4 MiB. History m=5,
пары Δx/Δr требуют ≈644 MiB сверх текущих states (без temporary arrays).
Это оценка storage, не измеренный RSS или runtime.

### 2. **Внешний Jacobi map страдает от смешения быстрого упругого блока и медленной LO-передачи; 2–3 упругих подшага могут дать лучший preconditioner**

libNEGF действительно решает elastic local-energy feedback внутри общего
inelastic iteration. Наш [Core _scattering_candidate:625](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/kernels.jl#L625)
строит все механизмы Jacobi из общего G. При сильном elastic redistributing
feedback outer Anderson вынужден аппроксимировать весь coupled Jacobian.
Проверяемая гипотеза — block preconditioning уменьшит raw full-map residual
на единицу стоимости, без смены механизмов/сетки/границ.

Опыт: из одного frozen state сделать A: 6 обычных full updates; B: 2 macrosteps,
в каждом 2 elastic-only inner updates с фиксированной LO self-energy, затем
один full update. Всё в одной общей budget 12 equivalent full sweeps / 20 min;
проверять full candidate и native residual после каждого macrostep. Объявить
inner stopping по elastic raw residual или максимум 2 шага; не по mixed step.
Сравнение по затраченным contraction/Dyson evaluations, а не только macrosteps.
Если B не лучше A по rΣ/rK на тот же wall cost, этот preconditioner отклонить.

Важно: finite-chain embedding должен обновляться согласованно с inner elastic
self-energy. Hilbert closure каждой elastic family остаётся прежней. Замена
на `K[Gr]` здесь запрещена как скрытая оптимизация: Core сам помечает
`elastic_retarded=:direct_green` отдельной finite-window closure.
Вложенное конечное число шагов меняет iteration map, но его принятие возможно
лишь по тому же полному endpoint residual; истинные fixed points совпадают,
когда все блоки действительно решены. Не подменять исходные данные LO=0.

### 3. **Оставшийся rsum после PI может ограничиваться краевым retarded tail, а не скоростью нелинейного solver**

Число .004922→.000406 на одном frozen Γ показывает чувствительность к real-part
quadrature. Оно не устанавливает отсутствие truncation error. Для причинной
closure вклад tail в Λ(E) не локален:
\[
\delta\Lambda(E)=\frac1{2\pi}\int_{\mathbb R\setminus[E_a,E_b]}
                       \frac{\Gamma(E')}{E-E'}\,dE'.
\]
Ошибка Λ меняет положения narrow resonances и N даже при фиксированной μ.
LibNEGF отрезает deep filled contribution иначе и без inelastic Σr; копирование
его contour split не является решением этой проблемы.

Дешёвый различающий опыт: на сохранённом Γ выполнить только retarded rebuild
и плотность при фиксированной μ для текущего окна и диагностического tail
continuation (например, fitted positive 1/E² wings, fit по двум разным edge
shells), отдельно от native baseline. Не прогонять новый SCBA. Проверить
чувствительность interior Λ/A/N к способу continuation и отчётный вклад краёв
в spectral sum. Дополнительный oracle: scalar Lorentzian с известным
полным-line Gr и finite-window log PV; это отделяет ошибку PI weights от
отсечённых tails. Фальсификация: обе разные разумные wings дают вклад существенно
меньше наблюдаемого остатка; тогда искать nonlinear fixed point/другую причину.

Это **диагностический вариант retarded closure**, не optimization и не принятие
S01. Если tail effect значим, следующий настоящий опыт — wider-window state с
tем же ΔE=ω/939, LO bath и физической identity, отдельная новая discretization
identity, baseline сохранён. Cap frozen diagnostic: 2 retarded rebuilds,
10 min, +≤1 GiB; никакой dense NE×NE Hilbert matrix (30704² doubles ≈7 GiB).
Для нового окна сначала считать число добавленных узлов и стоимость; здесь
полномочий на его запуск нет.

## Псевдокод одной идеи: thermal Γ solve с обязательным unprojected release

Иллюстрация H1; это не существующий API и не запущенный код. E и ω в одних
native units; `shift` — нулевое продолжение на **заявленном** окне, точный
целочисленный сдвиг, включая явный edge mask. `mu_for_number` находит μ для
исходной плотности без умножения G на normalization factor.

```julia
function thermal_map(Γ, problem, μ)
    ΣR = hilbert_PI(Γ, problem.E) - 0.5im * Γ
    GR = finite_chain_Dyson(problem, ΣR) # прежние 4+4, fixed U0
    A = im * (GR - adjoint_blocks(GR))
    f = fermi.(problem.E, μ, problem.kBT)
    μ = mu_for_number(A, problem.target_N, problem.weights)
    f = fermi.(problem.E, μ, problem.kBT)
    Γnew = sum(K[A] for K in problem.elastic_kernels)
    Am, Ap = shift(A, -939), shift(A, +939)
    fm, fp = shift(f, -939), shift(f, +939)
    Γnew += KLO[(Nlo + 1 .- fm) .* Am + (Nlo .+ fp) .* Ap]
    return Γnew, μ
end
for iteration in 1:maxiter
    Γnew, μ = thermal_map(Γ, problem, μ)
    r = raw_weighted_residual(Γnew, Γ)  # исходная норма, не AA step
    save_checkpoint(Γ, μ, r)
    r < tol && break
    Γ = real_hermitian_AndersonII(Γ, Γnew; history=5)
end
state = thermal_family_without_number_projection(Γ, μ, problem)
audit = full_native_unprojected_map(state, problem)
state = free_native_updates(state, problem; count=3)
report(audit, full_native_audit(state, problem))
```

На реальном Core осторожно: paired number projection переносит occupied/empty
weight, не просто умножает G<. Source:
[greens.jl::_normalize_keldysh_pair!:185](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/greens.jl#L185),
[scba.jl::_current_green:102](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/scba.jl#L102).
Следовательно, фиксированный GR ещё не делает весь current normalized map
линейной kinetic equation. Предложение «решить (I-K)G< GMRES» без явной обработки
number constraint было бы необоснованно.

## Итог области доказанного

Доказано чтением исходников: настоящий inelastic HEAD path, вложенность циклов,
формулы используемых shifted terms, uniform/rounded shifts, contour limitation,
DFTB+ старый gitlink и elastic-only interface, различие charge residual и
current-change stopping. Ни сходимость S01, ни физическая адекватность libNEGF
optical approximation, ни скорость предложенных способов здесь не измерялись.
Следующий предпочтительный опыт ведущего — H1 native release на m=939; H2 нужен
при провале/слишком медленном release, H3 — при устойчивом spectral/number defect.
