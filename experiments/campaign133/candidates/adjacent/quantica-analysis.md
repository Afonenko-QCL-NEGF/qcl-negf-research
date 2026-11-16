# Quantica.jl: алгоритмы, пригодные для QCL-NEGF

Дата проверки: 2026-11-15. Это read-only анализ исходников и первичных публикаций; Quantica не устанавливалась и не запускалась, новых научных расчётов здесь нет.

**Алгоритмический вывод.** Главная переносимая идея Quantica — задавать Green-функцию как вычисляемый объект и выбирать алгоритм **под наблюдаемую**: плотность интегрировать по гладкому комплексному контуру, пространственную границу устранять через surface Green/QZ, требуемые блоки вычислять лениво. Для нашей задачи Nb=5 общий sparse LU — слабая цель оптимизации. Наиболее полезны (1) причинное спектральное представление локальной Σ и независимая контурная плотность; (2) точная бесконечная cavity при нулевом смещении как отдельный вариант границы; (3) функциональный равновесный state, из которого временно материализуются только нужные поля. Ни один из этих выводов не свидетельствует о сходимости текущего SCBA.

## Зафиксированные исходники и границы

- Официальный репозиторий: [pablosanjose/Quantica.jl](https://github.com/pablosanjose/Quantica.jl). Read-only shallow clone: `.build/campaign133/adjacent/src-quantica`, размер 47 MiB.
- Проверенный commit **`ed2bb1ce1598cf5e03f721e70fe1b547f8103473`**, 2026-11-02, `master`; [Project.toml](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/Project.toml#L1-L4) объявляет 1.4.0. Это snapshot исходников, не утверждение о последнем опубликованном релизе.
- [LICENSE.md](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/LICENSE.md): MIT, copyright Pablo San-Jose and contributors (2020); при переносе существенного кода сохраняется уведомление. [CITATION.cff](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/CITATION.cff): software DOI [10.5281/zenodo.4762963](https://doi.org/10.5281/zenodo.4762963).
- QCL сравнивался по frozen `.build/campaign133/core-accepted-science`, **`e171caa60cce676c6c4238dd671e7c527f2bd977`**. Состояние superproject на входе `b515c42ef14c1f61151626280da0bcf100d116be`; live Core был `f50fa9b…` с внешними изменениями. Ничего из них не менялось.
- Внешние docs [GreenFunctions](https://pablosanjose.github.io/Quantica.jl/stable/tutorial/greenfunctions/) и [Observables](https://pablosanjose.github.io/Quantica.jl/stable/tutorial/observables/) помогли навигации. Научные утверждения ниже проверены по pinned source: moving `stable` и comments могут отставать. Например, actual `retarded_omega` добавляет `5eps`, а не встречающееся в docs `sqrt(eps)`.
- Проверенный scope Quantica: lattice tight-binding Green solvers, leads/self-energies, статический Hartree–Fock–Bogoliubov mean field. В просмотренном `src/solvers`, `observables.jl`, `meanfield.jl` **нет реализованного QCL electron–phonon SCBA или набора LO/acoustic/impurity/IFR/alloy**. Общая возможность задать Σ(ω) контакта не является таким SCBA.

## 1. Что именно вычисляет Quantica

### 1.1 Partial evaluation и выбор алгоритма

`greenfunction(oh, solver)` применяет solver к Hamiltonian и контактам. `g[selection]` фиксирует пространственные индексы (`GreenSlice`), `g(z)` — энергию (`GreenSolution`). Вызов slice возвращает копию; внутренний `call!` переиспользует буфер. Источники: [greenfunction.jl:5–20, 45–86](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/greenfunction.jl#L5-L86), [docstrings.jl:1660–1739](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/docstrings.jl#L1660-L1739).

Это substantive Julia abstraction: dispatch выбирает математический метод по типу solver/geometry/observable, а не один универсальный dense tensor pipeline. Малые орбитальные блоки представлены StaticArrays/HybridSparseMatrix; есть отдельные workspace-типы и `FastLapackInterface`. Не следует смешивать это с обещанием автоматической thread safety: lazy поля и aliases требуют независимых call-safe workspace на worker. В [SchurGreenSlicer.getproperty:782–798](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L782-L798) прямо описана опасность изменения aliased `h0` между вычислением и отложенным indexing.

### 1.2 SparseLU: factorization + нужные RHS, а не selected inversion

Для конечной системы строится

\[
D(z)=zI-H-\Sigma(z),\qquad D X=P_J,\qquad G_{IJ}=P_I^\dagger X.
\]

`apply(::GS.SparseLU,...)` создаёт sparse inverse-Green; `build_slicer` обновляет численные entries и вызывает **`lu(igmat)` заново для данного z**. `compute_or_retrieve_green` записывает выбранные единичные RHS и делает `ldiv!`. Все это реализовано в [sparselu.jl:54–86, 97–150](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/sparselu.jl#L54-L150).

Сохраняются LU и source buffers на одну energy solution. Если запрошены все нерасширенные столбцы, `unitg` становится плотным full-cell cache (стр. 136–138). Если нужны только q столбцов, storage RHS — O(Nq), factor storage — O(nnz(L+U)); время — factorization плюс q triangular solves. Fill зависит от geometry/order, общего O(N) обещать нельзя. При full request снова появляется O(N²) память. Это **не Takahashi/PEXSI-style selected inversion** и не доказанное переиспользование symbolic LU между энергиями: в проверенном пути виден новый `lu`, а не `lu!(oldfactor,newvalues;reuse_symbolic=true)`.

Для energy-dependent контактных Σ другие solvers используют T-matrix:

\[
T=(I-\Sigma G_0)^{-1}\Sigma,\quad G=G_0+G_0TG_0.
\]

Extended self-energy хранится блоками `(V†,g⁻¹,V)` и входит в расширенную систему без предварительного вычисления потенциально плохо обусловленного inverse. Реальные LU/gluing: [greenfunction.jl:609–627](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/greenfunction.jl#L609-L627).

### 1.3 Schur: исключение внутренних степеней свободы до QZ

Для nearest-cell lattice с `h₋=h₊†`, Bloch ansatz `ψ_j=λ^j φ` даёт quadratic pencil

\[
[h_+\lambda^2+(h_0-zI)\lambda+h_-]\phi=0
\]

(перестановка `h±` зависит от выбранного направления). Устойчивый surface solver выбирает decaying/outgoing subspace, а не обращает hopping. Для Im z>0 retarded modes затухают в сторону lead; на действительной оси propagating modes требуют знака скорости.

В актуальном коде `left_right_projectors` раскладывает `h₊=LR†` и выбирает **d=min(число поддерживаемых левых столбцов, правых столбцов)**, не SVD-rank с численным cutoff. Затем:

1. Формируется `K=zI-h₀+iΩP†P`; symbolic sparse structure и pointers сохранены.
2. `lu(K)`, затем `GL=K⁻¹L`, `GR=K⁻¹R` через `ldiv!`.
3. Из contractions `[R†;-L†]GL/GR` строится dense pencil `A−λB` размера **2d×2d**, затем `schur!(A,B)`.
4. `ordschur_retarded!` выбирает |λ|<1, а около unit circle использует projected velocity matrix и знак unpivoted LDL† pivots. Если классификация не даёт ровно d modes, код увеличивает Im z и повторяет QZ.
5. Из reordered Schur vectors хранятся surface factors; например right Σ имеет вид `P_R R Z₂₁ Z₁₁⁻¹R†P_R†`, но представляется extended blocks, а не обязательно материализуется inverse.

Здесь Ω — auxiliary shift pencil, который должен исключаться из окончательного surface result, **не физическая ширина уровня**. Проверка invariance при нескольких Ω нужна отдельно. Напротив, fallback, увеличивающий Im z для разделения modes, действительно меняет точку вычисления; его ошибку нельзя прятать в физической Γ.

Точные loci: [nearest_cell_harmonics/left_right_projectors:122–154](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L122-L154), [call!/pencilAB!:203–309](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L203-L309), [ordering:344–418](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L344-L418), [surface factors:22–26](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L22-L26).

Стоимость на z: sparse factorization + 2d solves + O(Nd²) contractions + O(d³) QZ. Workspace содержит O(Nd+d²) обычных массивов, **но также EE размера (N+d)² и EN размера (N+d)N**, поэтому текущий общий SchurWorkspace нельзя целиком рекламировать как O(Nd) память ([allocation:92–119](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L92-L119)).

QZ выполняется один раз на energy solution; `gL/gR/g∞`, surface Green и products создаются только при необходимости. Все они разделяют один `SchurFactorsSolver`, `skipsolve_internal=true` предотвращает повтор QZ: [lazy fields:606–655](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L606-L655), [slicer:745–821](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L745-L821).

**Первичная опора и её предел.** [San-Jose, *Deflation algorithm for the Green function of quasi-1D lattices*, arXiv:2105.08038v1 (2021), §§I–III](https://arxiv.org/pdf/2105.08038) выводит QZ-deflation и reconstruction generalized eigenvectors/Jordan blocks, необходимую при недиагонализируемом transfer matrix. Это серьёзная опора идеи reduction before QZ, но **не доказательство тождества с текущим кодом**: первая строка `schur.jl` ссылается на unpublished `scattering.pdf`, которого в clone нет. Актуальный shift-regularized pencil и velocity LDL† надо проверять по коду отдельно. Родственная независимая опора — [Wimmer, диссертация, гл. 3, особенно §§3.2–3.3](https://epub.uni-regensburg.de/12142/1/wimmer.pdf): стабильное использование Schur subspaces для singular hoppings.

В §IV статьи San-Jose implementation действительно называется Quantica.jl, но benchmark фиксирует **v0.5.0 / Julia 1.7**, а не проверенный здесь 1.4.0. Исторические graphene speedups нельзя переносить на Nb=5 или на нынешний workspace без измерения.

### 1.4 Две разные densitymatrix схемы

**Общая комплексная квадратура.**

\[
\rho=-\frac1\pi\operatorname{Im_H}\int_C f_\beta(z-\mu)G^R(z)\,dz,
\quad \operatorname{Im_H}X=(X-X^\dagger)/(2i).
\]

`densitymatrix(gs,path)` создаёт callable integrand. Для каждого scalar parameter x оно получает complex z=point(x), Jacobian dz/dx, f(z−µ), затем GreenSlice с Hermitian symmetrization. Array result сериализуется в flat buffer и интегрируется `QuadGK.quadgk!`; предусмотрен callback, `atol=1e−7` по умолчанию. Возвращаемый QuadGK `err` **discarded**; для научной qualification QCL его следует сохранять плюс проводить tolerance refinement. Источники: [observables.jl:398–449](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/observables.jl#L398-L449), [integrator.jl:158–179](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/integrator.jl#L158-L179).

Реализованы polygon, sawtooth и radial contour. Sawtooth вставляет µ в nodes и отрезает незанятый tail при T=0. Radial с θ=π/4 по умолчанию идёт с −∞ сверху к µ, затем к +∞ сверху при T>0; из-за G~I/z добавляется arc contribution `(θ/π) I` на выбранном срезе. См. [path construction/Jacobians:32–87](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/integrator.jl#L32-L87), [arc postprocess:414–425](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/observables.jl#L414-L425).

Это **adaptive quadrature, не Padé/Matsubara/Ozaki pole summation**. Ферми-полюса на `µ+iπ(2n+1)kBT` нельзя пересекать без residues. Radial wedges и sawtooth, разбитый в µ, обходят вертикальную линию полюсов. Произвольный пользовательский polygon этого не гарантирует. Docs предупреждают, что Cauchy deformation требует аналитичности G/пользовательской Σ; реализация её не доказывает. Контур не легализует подстановку неаналитической интерполяции Re/Im Σ.

Родственный, но иной алгоритм: [Ozaki, PRB 75, 035123 (2007), §II](https://www.openmx-square.org/tech_notes/CF_Fermi.pdf) получает неравномерные imaginary poles Fermi function из continued fraction и суммирует G в этих точках. [Ozaki–Nishio–Kino, PRB 81, 035116 (2010), §II.C](https://www.openmx-square.org/tech_notes/NEGF.pdf) разделяет equilibrium contour density и nonequilibrium real-axis correction. Это возможная следующая реализация density backend после квалификации аналитической Σ; наличие её в Quantica не утверждается.

**Специальная Schur density.** Вопреки названию solver это не интеграл scattering Green на реальных энергиях. Для чистого бесконечного Hermitian lattice вычисляется

\[
\rho_{mn}=\int_{BZ}\frac{d\phi}{2\pi}\,e^{i(m-n)\phi} f_\beta(H(\phi)-\mu).
\]

Fermi points получаются из propagating λ при µ; quadrature interval делится в этих точках. На каждом k выполняется `eigen_ws!`, умножение eigenvectors на occupations, затем извлечение slice. В 2D — вложенная quadrature. Имеется явный запрет boundaries и nonempty selfenergies: [schur.jl:1075–1180](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L1075-L1180). Преимущество — не разрешать узкие delta-like energy peaks, а интегрировать occupations собственных состояний. Прямого SCBA-обобщения нет: Σ(z) даёт spectral broadening/спутники, которые f(H) не содержит.

`GS.Spectrum` один раз diagonalizes nonparametric H и использует Lehmann sums; fixed params можно вычислить заранее ([spectrum.jl:38–82](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/spectrum.jl#L38-L82)). KPM density сейчас запрещает nonempty Σ и finite T ([kpm.jl:245–290](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/kpm.jl#L245-L290)). Поэтому ни KPM density, ни bare spectral sum не принимают наш all-five 70 K SCBA.

## 2. Связь с текущим QCL Core

В frozen `e171caa…`:

- `Reference.finite_chain_green`, `src/numerics/reference/periodicity.jl:172–264`, делает точную block cavity elimination для **конечной** открытой цепи с локальными ΣR/Σ</Σ>. `cavity_embedding_self_energy:293–376` строит shifted copies ячеек и явно обнуляет Σ за energy window. Докстрока требует отдельного предела по числу ячеек и окну; это не semi-infinite solver.
- `_embedding_self_energy_selected`, `src/numerics/optimized/production_operators.jl:872–903`, направляет `:finite_chain` в эту cavity концепцию с irreducible local scattering; задача сохраняет finite_chain4 на каждой стороне как baseline.
- `_retarded_green_production:1266–1314` уже использует per-worker 5×5 inverse workspace; `_green_state_production:1316–1343` материализует G>, A рядом с GR/G<. Переход на общий `SparseLU` не убирает NE×Nk workload.

Размер одного полного ComplexF64 поля `(30721,11,5,5)` — **135,172,400 bytes ≈128.9 MiB** без metadata; четыре Green поля ≈515.6 MiB, три Σ компоненты для пяти mechanisms ≈1.89 GiB. Это аналитический размер, не измеренный RSS. При Nb=5 общая sparse machinery может стоить дороже плотной; надо сокращать математическое состояние/число evaluation, а не только ускорять 5×5 loops.

Входное свидетельство кампании (передано ведущим, **не пересчитано этим анализом**): common baths, zero V, T=70 K; после 190 ordinary updates raw rK≈8.24e−5 и rΣ≈1.32e−4 остаются нерешёнными. Frozen-identical Γ nodal PV→linear-hat PI уменьшил rsum .004922→.000406; это отделяет ошибку retarded reconstruction от ещё открытой fixed-point convergence. Точное происхождение чисел должно ссылаться на родительский campaign artifact, а не на Quantica.

## 3. Три переносимые гипотезы по научной пригодности

### H1 — единая причинная Σ из Γ + отдельная contour density (наиболее близко)

В thermal equilibrium с общими µ,T/baths:

\[
\Gamma=i(\Sigma^R-\Sigma^A),\quad
\Sigma^<=if\Gamma,\quad\Sigma^>=-i(1-f)\Gamma,
\quad G^<=ifA.
\]

Вместо независимой интерполяции ReΣ и ImΣ задавать

\[
\Sigma^R(z)=\Sigma_\infty+\frac1{2\pi}\int_{E_a}^{E_b}\frac{\Gamma(E')}{z-E'}dE'.
\]

**Краткое доказательство причинности для discretized representation.** При Hermitian PSD Γ узлов линейная hat-interpolation сохраняет PSD на интервале. Для y=Im z>0,
`−Im_H Σ(z)=∫ y Γ(E')/[2π((Re z−E')²+y²)]dE' ⪰0`.
При Hermitian H, `Im_H[zI−H−Σ(z)]⪰yI`, поэтому resolvent существует в upper half-plane. `Σ(z)=Σ∞+O(1/z)` обеспечивает `G(z)=I/z+O(1/z²)`. Это даёт independent zeroth spectral normalization **при интегрировании всей спектральной меры, включая возможные реальные bound-state delta poles**. Оно не доказывает правильность Γ SCBA или физическую полноту конечного окна.

На interval [a,b], h=b−a, ΔΓ=Γb−Γa exact hat integral:

\[
\int_a^b\frac{\Gamma_a+(E-a)\Delta\Gamma/h}{z-E}dE
=[\Gamma_a+(z-a)\Delta\Gamma/h]\,[\log(z-a)-\log(z-b)]-\Delta\Gamma.
\]

Boundary value даёт согласованные PV и −iΓ/2; off-axis формула позволяет Quantica-style density без реального energy mesh для **наблюдаемой**. Нужны stable `log1p`/asymptotic branches вдали от interval, иначе cancellation уничтожит хвост. Это математическое продолжение выбранной Γ-representation, а не нестабильное inverse analytic continuation из noisy imaginary data.

**Различающий минимальный опыт, не запущен:** freeze одну Γ с каждого mechanism/один µ и тот же finite-chain4. Проверить несколько k, затем все Nk: rho из radial θ=π/6 и π/4, два tolerance levels; отдельно исходная real-axis hat integration. Сначала analytic scalar Lorentzian и finite Hermitian chain через точный eigenvector Fermi sum, затем frozen QCL. Сохранить norm matrix rho differences, per-orbital absolute errors, QuadGK error estimate/count, tail/bound-state accounting, `zG→I`, min eigenvalue Γ. Proposed oracle tolerances: manufactured cases `||Δrho||F≤1e−8 max(1,||rho||F)`; QCL не ослаблять существующий scientific threshold, quadrature contribution просить ≤0.1 его абсолютной error budget. Если positivity, contour agreement или tails не контролируются — результат `not_measured/failed`, не projection Γ на PSD без отдельной identity.

**Лимит и смысл:** contour может доказать, что наблюдаемая интегрируется дешевле/надёжнее, и отделить missed poles/window от PV error. Она **не заменяет** real-energy shifted SCBA map Γ→Γ: electron–phonon retarded Σ не получается наивной заменой E→z в реальных Fermi factors. Γ-only equilibrium pilot и Anderson-II уже принадлежат ведущей кампании; эта гипотеза даёт независимую density qualification, а не дублирующий запуск.

Особый случай: compact-support Γ может оставить связанные состояния вне своего support, где Γ=0. Contour density содержит их equilibrium occupation; чистая injection formula `G<=GR Σ< GA` без bound-state/initial-state contribution может эту delta weight не определить. Такое расхождение требует отдельной классификации и не является автоматически ошибкой PI или свидетельством несходимости SCBA.

### H2 — точная semi-infinite cavity при нулевом bias (сильное изменение границы)

В zero V все cell onsite/local Σ одинаковы. Для frozen causal Σ(z) можно решить

\[
g_R=[zI-H_0-\Sigma_{loc}(z)-Tg_RT^\dagger]^{-1},\quad
G_{00}=[zI-H_0-\Sigma_{loc}-Tg_RT^\dagger-T^\dagger g_LT]^{-1}
\]

через ordered QZ либо stable decimation. Это устраняет finite-chain standing-wave poles и reflections; не требует 4→8→16 последовательных neighbours внутри каждого energy evaluation. В нашем маленьком basis QZ размера ~10 дешёвый, но если T dense, d=Nb=5 и deflation не даёт dimension reduction. Фактический numerical rank T надо измерить отдельно; truncation по singular values была бы новой approximation.

**Независимый oracle:** для того же frozen Σ и комплексных z интеграл
`G00(z)=∫₋π^π [zI−H0−Σloc(z)−T exp(ik)−T† exp(−ik)]⁻¹ dk/(2π)`.
На scalar chain есть analytic square-root surface Green (`Reference.chain_surface_green`, periodicity.jl:268–275), не reference/optimized copy одной Riccati recursion. QZ и dense Bloch integral различают неверную ветвь/ориентацию и mistaken bulk-as-surface.

**Минимальный опыт, не запущен:** analytic scalar, затем frozen 5×5 Σ: 8–16 upper-plane z и несколько реальных nodes с физической Γ; сравнить QZ, Bloch quadrature, c4/c8/c16. Сохранять Riccati raw residual `||g⁻¹−(zI−H−Σ−TgT†)||` с energy units, causality, density, energy-resolved A. Для околонулевой матрицы использовать абсолютный масштаб из declared energy unit, не деление на нулевой g/Σ. Если finite-chain trend отсутствует, не объявлять c4 бесконечным. Потом только отдельная SCBA identity `semi_infinite_zero_bias`; старый finite-chain run остаётся не принятым этой заменой.

**Лимит:** текущий Quantica velocity-ordering предполагает Hermitian lead structure; произвольную dissipative energy-dependent Σlocal нельзя считать проверенно поддержанной простой подстановкой в `h0`. Для Im z>0 и passive Σ decay-subspace formulation математически естественна, но требует вышеописанной qualification. При finite field ячейки energy-shifted, Bloch reduction неприменим. Nonlocal intercell scattering также меняет pencil.

### H3 — quantity-oriented equilibrium state и ленивые contractions (инженерный второй этап)

Перенести **структуру API**, а не весь Quantica DSL: typed `EquilibriumSpectralState(Γ,μ,T,representation)` → callable `LocalResolvent(state,H,embedding)` → `density/SCBAResidual`-specific materialization. Γ, FDR и аналитическая retarded reconstruction должны быть единственным научным source, а G</G>/A — views/временные chunks. Для SCBA нужны реальные E±ħωLO и momentum contractions; их plan хранится отдельно и проверяет conservation/adjointness. Density backend может быть contour, SCBA map пока uniform real mesh.

Хранение одного полного Γ вместо трёх ComplexF64 Σ полей даёт factor 3 по полям механизма; packed Hermitian Γ с Nb² real numbers теоретически даёт factor 6, но packing может замедлить contractions. Это расчёт storage, не benchmark. Общий raw residual каждого mechanism проверять на исходной Γ→Γ map; FDR-конструкция не является самостоятельным convergence oracle.

**Минимальная qualification:** на frozen Γ проверить полные mechanism candidates и исходные equation residuals перед/после lazy materialization; достаточно одного неизменного map evaluation без новой SCBA campaign. Отдельно measured peak RSS/allocations/time и consumer API. Oracle — analytic detailed balance/FDR и независимая contraction на маленьком domain; agreement старого/нового code path — supplementary. Существенный риск — aliasing immutable-looking closures на shared mutable workspace: worker-local workspace обязателен, результаты checkpoint должны владеть своими arrays. Core уже имеет per-worker inverse buffers; очередное ручное переписывание loops — не новая научная схема.

## 4. Julia-псевдокод (23 строки; новая схема, не готовый API Quantica)

```julia
struct EquilibriumSpectralState{R,B}
    gamma::R                 # one Hermitian Γ per mechanism, E and k
    bath::B                  # common μ and kBT, in one declared energy unit
end
struct LocalResolvent{S,H,C}
    state::S; hamiltonian::H; cavity::C
end
function (g::LocalResolvent)(z, k, ws)
    Σ = causal_hat_integral(g.state.gamma, z, k, ws) # same PI on/off axis
    embed = surface_selfenergy(g.cavity, z, k, Σ, ws)
    return solve_local!(ws, z, g.hamiltonian[k], Σ, embed)
end
function density(g, k, ws; atol)
    f(z) = stable_fermi(z - g.state.bath.μ, g.state.bath.kBT)
    ρ, err = contour_quadgk(z -> hermitian_integrand(f(z), g(z,k,ws)),
                          radial_path(g.state.bath.μ, π/4); atol)
    return (; rho = ρ + arc_identity(π/4), quadrature_error = err)
end
function raw_map(state, real_grid, worker_workspaces)
    localG = real_axis_green_chunks(state, real_grid, worker_workspaces)
    Γnew = all_five_scba_gamma(localG, state.bath, conservative_shift_plan)
    return Γnew, physical_raw_residual(Γnew, state.gamma)
end
```

Для H1 сначала оставить `cavity=finite_chain4`, чтобы не смешивать контур и новую границу. `surface_selfenergy` в H2 вводится отдельной identity. `all_five_scba_gamma` не вычисляет LO/acoustic/etc при complex z: это явный real-axis map. Псевдокод не означает, что действующие Core checks/units/bath conventions можно пропустить.

## Решение для текущей кампании

**Сейчас:** использовать H1 как ограниченный независимый frozen-Γ diagnostic после/рядом с equilibrium-PI pilot ведущего; не делать package migration. **Следующий научно различающий вариант:** H2, если окно позволяет, с фиксированными Σ и явным boundary identity. **После научного исхода:** H3 для памяти/API; оценка должна включать time-to-accepted-result, а не isolated inverse throughput. Ни quantica source, ни перечисленные публикации не дают основания принять нерешённый all-five raw SCBA по малому mixed step или по одному выполнению программы.
