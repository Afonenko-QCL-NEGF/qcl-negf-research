**Переносимые стационарные GF-методы: Kwant и TranSIESTA, разбор исходников 2026-11-15**

Короткое решение: для нынешних `Nb=5`, `finite_chain=4+1+4` не переносить большой sparse/scattering framework. Core уже выполняет точное исключение полостей. Наиболее полезны два небольших независимых оператора: **doubling/decimation для проверки предела числа клеток при zero V** и **контурный интеграл плотности из причинного продолжения замороженной Γ**. Третий перенос — **масштаб ошибки около закрытого канала, заданный полной Σ/оператором, а не исчезающей Γ или током**. Ни один из них не принимает SCBA с ненулевыми исходными `r_K/r_Σ`.

Это анализ кода и постановка различающих опытов. Julia, Kwant, SIESTA и сторонние тесты не запускались; пакеты не устанавливались. Числа текущей кампании ниже получены из задания ведущего, а не перепроверены по массивам.

**1. Алгоритмы, которые действительно можно перенести**

Пусть `D(z)=zI-h-Σlocᴿ(z)`, а `T=H[j,j+1]`. Для правой однородной полубесконечной цепи поверхностный GF удовлетворяет

\[
 g_s(z)=[D(z)-Tg_s(z)T^\dagger]^{-1}.
\]

Это **cavity**, а не центральный bulk GF. Центральный GF бесконечной однородной цепи имеет две поверхности: `G=[D-TgR T†-T†gL T]⁻¹`. Именно это различие существенно для QCL. Для конечного числа одинаковых клеток применяется та же Schur-рекурсия с нулевым внешним embedding; она точно решает конечную задачу и не обязана совпадать с бесконечной.

**Doubling López-Sancho.** Начинаем с `A₀=T`, `B₀=T†`, `D₀=D`, `Ds₀=D`. На каждом шаге исключаем каждую вторую клетку:

\[
 X_n=D_n^{-1}A_n,\quad Y_n=D_n^{-1}B_n,
\]
\[
 D^s_{n+1}=D^s_n-A_nY_n,\quad
 D_{n+1}=D_n-A_nY_n-B_nX_n,
\]
\[
 A_{n+1}=A_nX_n,\quad B_{n+1}=B_nY_n.
\]

После сходимости `g_s=(Ds)⁻¹`. Один шаг удваивает расстояние между оставшимися клетками. Один LU `D_n` и два набора RHS дают `X,Y`; не нужно явно строить `inv(D_n)` для этих произведений. Работа — `O(n_iter Nb³)`, рабочая память — `O(Nb²)` на энергию/k. Число шагов логарифмично эффективной длине при затухании; у реальной оси с очень слабой диссипацией быстрый предел не гарантирован. Нулевой hopping — отдельный тривиальный случай, а не ошибка физики.

В SIESTA эта схема находится в [`SSR_sGreen_NoDOS`, строки 428–637](https://gitlab.com/siesta-project/siesta/-/blob/71c8602911fd632fa3331f2588ea55b4d6ed684e/Src/m_ts_electrode.F90#L428): `zgesv(no,2*no,...)` на строке 561, обновления на 571–583, `8*no²` complex workspace на 492. Стоп на 585–588 задан абсолютным `maxabs` поправки поверхности, то есть в единицах энергии, не безразмерным residual уравнения. При переносе нужен также независимый residual `||(D-TgT†)g-I||/(||D-TgT†||·||g||+√Nb)` и проверка причинности. Копировать только их stop как научную приёмку нельзя. Первичная публикация: [López Sancho et al., J. Phys. F 15, 851 (1985), DOI](https://doi.org/10.1088/0305-4608/15/4/009); полный текст DOI роботам недоступен, формулы здесь проверены по указанному коду как Schur elimination.

Применимость к S01: при `U=0` клетки имеют одинаковые `h` и замороженную локальную `Σᴿ`; Schur elimination остаётся допустимым для диссипативной причинной Σ. Никакой unitary scattering matrix для этого не требуется. При ненулевом поле `D_j(E)` различаются, поэтому простой однородный doubling неприменим. При межклеточной scattering self-energy блоки `T,T†` тоже перестают задавать всю связь. Для lesser/greater равновесный shortcut разрешён только после проверки общих `μ,T` и thermal baths: `Σembˡ=ifΓemb`, `Σembᵍ=-i(1-f)Γemb`. Это **новая бесконечная постановка/оракул области**, не эквивалентная оптимизация принятой конечной цепи.

Высокоуровневый Julia-псевдокод независимого retarded оператора (20 строк; это не исполнявшийся опыт):

```julia
function surface_doubling(D0, T; atol_E, rtol, maxiter=80)
    n = size(D0, 1); eye = Matrix{ComplexF64}(I, n, n)
    iszero(norm(T)) && return D0 \ eye
    D, Ds, A, B = copy(D0), copy(D0), copy(T), Matrix(T')
    for iteration in 1:maxiter
        F = lu(D)                         # pivoted, one factorization
        XY = F \ hcat(A, B)                # two block right-hand sides
        X, Y = XY[:, 1:n], XY[:, n+1:2n]
        cR, cL = A * Y, B * X
        Ds_next, D_next = Ds - cR, D - cR - cL
        A_next, B_next = A * X, B * Y
        done = max(norm(A_next), norm(B_next), norm(cR)) <=
               atol_E + rtol * max(norm(Ds_next), norm(D0))
        D, Ds, A, B = D_next, Ds_next, A_next, B_next
        if done
            g = lu(Ds) \ eye
            certify_surface_equation_and_causality(D0, T, g)
            return g
        end
    end
    error("surface limit not converged")
end
```

Не добавлять скрытое `iη`. Для полубесконечного Hermitian предела сначала законно брать `z=E+iη`, а затем отдельно исследовать `η→0⁺`; `η` должен остаться частью identity. Для замороженной диссипативной Σ на real E строгое затухание также надо проверить, особенно на краях окна и в подпространствах с почти нулевой Γ.

**Причинный контур из linear-hat Γ.** В равновесии плотность матрицы задаётся `ρ=∫f(E) A(E)dE/(2π)`, `A=i(Gᴿ-Gᴬ)`. Тот же интеграл можно вычислить по комплексному контуру или конечным полюсам рационального приближения Fermi. Вдали от реальной оси узкие резонансы сглажены аналитичностью; это переносится как независимая проверка интеграла, а не как замена energy-resolved SCBA.

Ключевой bridge к текущему PI-оператору:

\[
 \Sigma^R(z)=\Sigma_\infty+\frac{1}{2\pi}
 \int_{E_{min}}^{E_{max}}\frac{\Gamma(x)}{z-x}\,dx,\qquad \Im z>0.
\]

Для линейной Γ на `[a,b]`, `Γ(x)=Γa+M(x-a)`, матричный вклад **аналитический**:

\[
 [\Gamma_a+M(z-a)]\,[\log(z-a)-\log(z-b)]-M(b-a).
\]

Граничное значение даёт одновременно ту же PV real part и `-iΓ(E)/2`, если соглашения по Γ и краевым функциям совпадают. Если узловые Γ положительно полуопределены, их линейная интерполяция также PSD, и `-Im Σᴿ(z)≽0` для `Im z>0`. Это независимое математическое ограничение на продолжение; его нельзя получить произвольным spline отдельно для Re/Im Σ. На больших `|z|` нужна устойчивая форма через `log1p`/моменты, иначе разность больших слагаемых теряет точность.

Так можно вычислять `Gᴿ(z)` на нескольких десятках/сотнях комплексных точек, применяя **ту же конечную** `4+1+4` цепь, и сравнивать `ρcontour` с `ρreal` при фиксированной Γ. Это сохраняет область и локальную модель; смена схемы интеграла остаётся отдельной численной identity. `Σ∞` и high-frequency moment должны быть явно заданы из модели; ноль нельзя молча подставить вместо неизвестной статической части. Конечное support Γ — предположение, которое нужно обозначить, особенно при исследовании хвостов окна.

Фактическая опора SIESTA: [`contour_poles`, 1310–1336](https://gitlab.com/siesta-project/siesta/-/blob/71c8602911fd632fa3331f2588ea55b4d6ed684e/Src/m_ts_contour_eq.f90#L1310) задаёт Matsubara-полюса `μ+iπkT(2ν+1)`; [`contour_continued_fraction` / `Ozaki_residue`, 1339–1473](https://gitlab.com/siesta-project/siesta/-/blob/71c8602911fd632fa3331f2588ea55b4d6ed684e/Src/m_ts_contour_eq.f90#L1339) строит continued-fraction полюса через generalized symmetric eigenproblem и отдельно учитывает zeroth moment. Это существенно: простая обрезанная Matsubara-сумма без moment/tail correction не равна плотности. В [`ts_trig`, 359–449](https://gitlab.com/siesta-project/siesta/-/blob/71c8602911fd632fa3331f2588ea55b4d6ed684e/Src/m_ts_trig.F90#L359) реально выполняются complex-energy GF и накопление DM, а при нулевом bias обработка equilibrium заканчивается до nonequilibrium части. Первичная методическая статья: [Ozaki, Phys. Rev. B 75, 035123 (2007)](https://doi.org/10.1103/PhysRevB.75.035123); [авторская копия](https://www.openmx-square.org/tech_notes/CF_Fermi.pdf). Обзор соответствующей реализации: [Papior et al., CPC 212, 8–24 (2017), §§IIA, IIIA, IIID](https://arxiv.org/pdf/1607.04464).

Стоимость предложенного frozen-source audit: `O(Nc NE Nk Nb²)` для прямого суммирования hat-вкладов плюс `O(Nc Nk P Nb³)` для GF. Это не бесплатное получение Σ(z): при `Nc=128`, `NE=30721,Nk=11,Nb=5` получается около `1.08·10⁹` комплексных весовых вкладов для одной матричной Γ. Потоковая обработка не требует `Nc` полных массивов; одна полная complex Γ занимает примерно 129 MiB, scratch — `O(Nb²)`. Сначала достаточно одного физически выбранного k и `Nc=16,32,64`; оценку wall time следует измерять отдельно. Этот аудит не требует нового полного SCBA.

Контур не переносится автоматически на исходные LO collision convolutions, `Gˡ` и `Σˡ`: они содержат real-energy occupations и shifted products, а не одну retarded аналитическую функцию. Даже точная контурная плотность не доказывает `r_K`, `r_Σ`, дискретизацию рассеяния или отсутствие FDR-дефекта. При конечном bias остаётся неравновесная real-axis поправка; thermal FDR-подстановка уже недопустима.

**2. Что именно делает Kwant**

Официальный сайт [указывает основной GitLab-репозиторий](https://kwant-project.org/community). Исследован его `main` на `e1603f8817e1459f2544ae84831e78a3e4b79430` (commit 2026-10-22), локально `adjacent/src-kwant`, 2.7 MiB. [LICENSE.rst](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/LICENSE.rst) содержит двухусловную BSD-лицензию. SIESTA: официальный [`siesta-project/siesta`](https://gitlab.com/siesta-project/siesta), `71c8602911fd632fa3331f2588ea55b4d6ed684e` (`dev`, commit 2026-08-05), `adjacent/src-siesta`, 67 MiB; [COPYING](https://gitlab.com/siesta-project/siesta/-/blob/71c8602911fd632fa3331f2588ea55b4d6ed684e/COPYING) — GPL version 3, заголовки исследованных Fortran-файлов отсылают к нему. Это provenance анализа, не предложение переносить GPL-код в Core. Все предлагаемые формулы/псевдокод сформулированы заново.

**SVD + Schur/QZ вместо наивных transfer-eigenvectors.** [`setup_linsys`, 200–391](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/physics/leads.py#L200) делает SVD hopping, отбрасывая singular values ниже `tol·epsmach·smax` (260–261). В невырожденном случае возможна обычная transfer matrix; в другом случае строится pencil размером `2r×2r`, где `r=rank(T)`. Используется basis с `sqrt(s)` на обеих половинах, что избегает деления на малые singular values. При необходимости алгебраически добавляется anti-Hermitian stabilizer `i(uu†+vv†)` и одновременно меняется pencil/reconstruction (318–370). Это **компенсированное преобразование задачи**, а не разрешение добавить физическую imaginary broadening в QCL.

[`unified_eigenproblem`, 394–466](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/physics/leads.py#L394) строит Schur либо generalized Schur (QZ), а [`compute_block_modes`, 872–909](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/physics/leads.py#L872) сохраняет ортонормированный Schur basis затухающего подпространства вместо индивидуальных плохо обусловленных evanescent eigenvectors. Для propagating modes выполняется current normalization. [`StabilizedModes.selfenergy`, 179–196](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/physics/leads.py#L179) берёт только outgoing propagating+evanescent часть и использует solve. Стоимость lead-проблемы плотная, `O(Nb³)` времени и `O(Nb²)` памяти, с меньшим pencil при `r<Nb`. Это полезнее decimation в отдельных почти вырожденных когерентных контактах, но требует самостоятельного аккуратного порта, не одной функции `eig(T)`.

Ограничение: стандартные propagating/current-paired modes основаны на Hermitian периодическом lead Hamiltonian. Подстановка `h+Σlocᴿ(E)` с поглощением меняет задачу; деление мод пополам и normalization в `make_proper_modes` уже нельзя считать доказанными. Поэтому для нынешнего frozen dissipative chain безопаснее перенос Schur-decimation, а Kwant — независимый coherent-limit oracle (`Σloc=0`, scalar-chain analytic answer или заданный Hermitian tight-binding lead).

Мелкое, но важное наблюдение исходников: docstring второго `stabilization` флага обещает generalized problem, тогда как фактические строки 380–390 при true возвращают `B⁻¹A`. Поэтому ручные флаги нельзя переносить из документации без теста выбранного pencil. Это не заявленная ошибка результата Kwant: здесь не выполнялись ни тесты, ни трассировка публичного вызова.

**Разреженный LU и selected solves; не selected inversion.** [`greens_function`, 403–497](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/solvers/common.py#L403) формирует RHS только для интерфейсов заданных incoming leads, факторизует всю sparse LHS один раз и оставляет запрошенные output rows. [`mumps.Solver`, 107–126](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/solvers/mumps.py#L107) решает RHS пакетами (`nrhs=6` default) и сохраняет `kept_vars`. Это не алгоритм selected inversion типа Takahashi/PEXSI; доступ к нужным столбцам обратной матрицы получается треугольными решениями. Рабочая память включает LU fill-in плюс `O(N·nrhs)` RHS, а не только возвращаемую маленькую матрицу. Nested-dissection ordering выбирается из METIS/SCOTCH/auto (94–97). Реальное преимущество зависит от графа и fill-in, а не от sparse-метки.

Для квазиизмерной block chain с `P` одинаковыми по размеру блоками RGF даёт `O(P Nb³)` и `O(Nb²)` потоковой памяти для одной центральной клетки; полный inverse имел бы `O(P² Nb²)` памяти. У нас `P=9,Nb=5`, центральные GF нужны для рассеяния на всех `E,k`, и уже есть cavity R/</> recursion. Большой sparse LU не имеет очевидного основания: он увеличивает инфраструктуру и не меняет stiffness/нелинейность карты Σ. Kwant `greens_function` сам описывает себя как обычно более медленный/менее устойчивый путь, чем scattering matrix, для testing/RGF compatibility (444–448), но это сравнение относится к coherent observables, а не к полноценной SCBA.

Для сопоставления: [`SIESTA invert_TriMat`, 87–180](https://gitlab.com/siesta-project/siesta/-/blob/71c8602911fd632fa3331f2588ea55b4d6ed684e/Src/m_trimat_invert.F90#L87) рекурсивно вычисляет выбранные diagonal/nearest-neighbour blocks (`calc_parts`), а [`invert_BiasTriMat_rgn`, 453–701](https://gitlab.com/siesta-project/siesta/-/blob/71c8602911fd632fa3331f2588ea55b4d6ed684e/Src/m_ts_trimat_invert.F90#L453) дополнительно выбирает регионы/столбцы. Таким образом, даже соседняя большая transport-программа использует структуру chain, когда она доступна.

**3. Равновесие и conditioning: что является свидетельством**

Coherent Landauer equilibrium `Iα=(q/h)∫Σβ Tαβ(E)[fα-fβ]dE=0` при общих Fermi functions встроено в формулу. Оно не является независимым тестом collision SCBA. Более того, Kwant [`BlockResult.transmission`, 703–739](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/solvers/common.py#L703) может **достраивать отсутствующие transmissions из предполагаемого current conservation**, когда Hamiltonian был проверен на Hermitian. Это удобно для расчёта, но такое достроенное равенство нельзя затем считать независимой проверкой conservation.

Для interacting equilibrium остаётся точечное FDR-условие `Gˡ=ifA`, `Σαˡ=ifΓα`, `Σαᵍ=-i(1-f)Γα` для всех согласованных источников. Тогда collision integrand `Tr[ΣαˡGᵍ-ΣαᵍGˡ]` равен нулю. В equilibrium-per-reservoir тесте нужно вычислить обе стороны независимо по исходной unconstrained map, а не сконструировать вторую из первой. При LO detailed balance нужны правильные shift, Bose occupation и thermal baths; нулевой bias сам по себе этого не обеспечивает. В системе с магнитным полем возможен equilibrium circulating bond current; нулевой bath flux не равен требованию занулить всякий bond current. У нынешнего S01 из задания предпосылки thermal/common baths согласованы, но native map audit всё равно необходим.

Перенос **near-zero scale** виден в Kwant [`GreensFunction._transmission`, 965–985](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/solvers/common.py#L965): при числе открытых каналов tolerance масштабируется `||Σ||∞`, поскольку `||Γ||≈0` в закрытом канале. Значение `10⁶ eps` — их конкретная эвристика, не физический допуск QCL. Для QCL переносима логика: ошибка anti-Hermitian разности `Γ=i(ΣR-ΣR†)` должна соотноситься с масштабом обоих вычитаемых операторов/конструкции; residual тока около нуля требует абсолютной части в A/cm², а не деления на нулевой expected current. Отдельно сообщать relative error там, где сигнал разрешён.

Core на принятой ревизии уже содержит backward Dyson residual, [`validation.jl:33–44`](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/validation.jl#L33), и production хранит `cond(D)` на каждой `(E,k)` [`production_operators.jl:1266–1313`](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/production_operators.jl#L1266). Малый backward residual при большой κ не исключает forward error; нужно различать `κ·eps` и физические SCBA `r_K/r_Σ`. С другой стороны, κ без большого residual/разницы с независимым solve не доказывает, что проблема именно в linear algebra.

**4. Три ранжированные гипотезы и минимальные опыты**

| Приоритет | Гипотеза / перенос | Минимальное различающее свидетельство | Что результат не установит |
|---|---|---|---|
| 1 | **Часть спектральной ошибки/узких пиков связана с конечной длиной 4+1+4; doubling позволяет дешёво проверить предел при frozen Σ.** | На одной и той же Γ/Σ, `U=0`, выбранных k и энергиях с большим spectral-sum вкладом сравнить существующие конечные `Pside=4,8,16,32` с независимым decimation surface. Проверить surface equation, `-Im g≽0`, central spectral matrix и последовательно интегральные `r_sum`, charge. Сначала scalar-chain analytic oracle, затем реальный Nb=5. Разница должна уменьшаться с длиной; если не уменьшается, различить branch/conditioning и истинное отсутствие затухания. | Даже большой finite-size effect не объясняет ненулевой raw residual **внутри фиксированной finite_chain4 карты**. Переход на бесконечную цепь — новая identity и собственный SCBA. |
| 2 | **Linear-hat причинное продолжение + complex contour даст независимую оценку равновесной charge и отделит real-axis quadrature от недостатка окна.** | При frozen Γ и той же finite_chain4 вычислить `ρcontour` для 16/32/64 полюсов или fixed contour refinement на одном выбранном k. Первично сравнить с finite Hermitian diagonalization/analytic scalar resolvent, где answer известен. Затем сравнить с real-axis плотностью на native Γ interpolation; отдельно оценить tails и проверить zeroth moment. Если контур стабилен, а real integral меняется с quadrature, причина локализована в интегрировании; если оба меняются с support Γ, окно ещё не принято. | Не принимает collision convolution, pointwise Γ, FDR или SCBA. Новая неаналитическая интерполяция Σ(z) уничтожила бы независимость. |
| 3 | **Околонулевые current/Γ witnesses могут быть ограничены cancellation/conditioning, а не физическим нарушением.** | Без нового SCBA выбрать блоки с max κ и min Γ; независимо собрать `D`, pivoted LU solve и `A=i(G-G†)`, source/current terms с compensated summation. Сопоставить абсолютные дефекты с `eps·construction_scale` и backward error; проверить энергетически интегрированный дефект с исходными весами. Сравнить с блоками, где сигнал велик. При необходимости только эти 5×5 блоки — high precision oracle. | Если raw `r_K/r_Σ≈10⁻⁴` остаются намного выше арифметического бюджета, это не roundoff и снижать их допуск нельзя. |

Все три — будущие ограниченные опыты ведущего; в этом пакете они **не выполнялись**. Предложенное начальное окно для каждого — один frozen artifact, не более 1 CPU / 1 GiB дополнительной памяти, до 2 минут на пробу; если первая проба не укладывается, сохранить причины/стоимость и пересмотреть план до запуска полной сетки. Выход — компактные scalar/matrix diagnostics с `(E,k)` локаторами, не дублирующие исходные массивы. Ресурсные числа — предложение для отделения механизма, не обещание фактической скорости.

Для текущего ведущего решения `Γ-only thermal PI + Anderson II` эти методы — следующие независимые проверки. Приоритет самого native unconstrained full-map audit выше запуска новой схемы контактов: equilibrium-constrained решение обязано выдержать исходные уравнения и исходный набор bath/scattering assumptions.

**5. Локальная привязка и границы вывода**

Superproject на старте `b515c42ef14c1f61151626280da0bcf100d116be`; рабочее дерево имело чужие изменения Core/Runner/research, они не изменялись. Прочитан только нужный Core accepted science worktree `e171caa60cce676c6c4238dd671e7c527f2bd977` (`.build/campaign133/core-accepted-science`): finite chain reference [`periodicity.jl:172–264`](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/periodicity.jl#L172), production cavity [`production_operators.jl:235–292`](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/production_operators.jl#L235). В reference явно сказано, что finite-domain требует size convergence и не допускает межклеточную scattering self-energy (155–170). Сам новый пакет анализа меняет только этот Markdown и две read-only source copies.

Входы из задания: S01 zero V, 70 K, общие thermal baths, 5 LO + acoustic-elastic + impurity + IFR + alloy; `Nb=5,Nk=11,NE=30721,ΔE≈39 μeV,Nz=96,Pbasis=3`, finite chain по четыре клетки с каждой стороны, fixed U=0. Переданный baseline после 190 ordinary updates от одного FD-init имеет `r_K≈8.24e-5,r_Σ≈1.32e-4`, то есть SCBA не принята. Переданный frozen-local-Γ эксперимент nodal PV → linear-hat PI снизил `r_sum≈0.004922→0.000406`; это поддерживает изучение причинной product integration, но не доказывает принятый бесконечный предел, full-map solution, grid adequacy или экспериментальную валидацию. Архивы и числа не открывались/не менялись.

Первичная статья Kwant: [Groth et al., New J. Phys. 16, 063065 (2014), §§II, VII, Appendix C](https://arxiv.org/html/1309.2926v2). Её выигрыш sparse solver над RGF относится к изученным когерентным геометриям; перенос этого benchmark на 5×5 SCBA нельзя обосновать одним названием NEGF. Установка SIESTA/Kwant, импорт большого framework, изменение Core, новый научный расчёт или публикация в данном пакете не производились.
