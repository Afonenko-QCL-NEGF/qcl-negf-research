# Неравномерная энергия: exact panel PI с матричным Cauchy treecode

Предлагаю один следующий алгоритм: **оставить положительный кусочно-линейный
Γ как определение self-energy, а плотную nonuniform Hilbert matrix заменить
деревом аналитических panel moments с точным ближним PV**. Это открывает
неравномерные target energies без `N_E²` памяти. Само дерево не уменьшает
число Dyson solves и не доказывает, что adaptive grid будет экономнее
существующего uniform PI FFT. Такое уменьшение надо отдельно измерить на
одной frozen PI модели относительно уже полученной Matsubara density.

Статус: проект алгоритма и bounded протокол, **не реализация и не запуск**.
Новые Julia/Pkg/SCBA/third-party build не выполнялись, Core/Runner не менялись.
Все performance, accuracy и convergence outcomes ниже `not_measured`.
Применён `qcl-physics-audit`; источники — actual Core и первичные статьи.

## Фиксированные исходники и предел вывода

На момент чтения superproject HEAD
`b515c42ef14c1f61151626280da0bcf100d116be`; рабочие gitlinks Core/Runner
изменены, research имеет untracked content. Чужие изменения не затрагивались.
Рассмотрен отдельный accepted checkout
`.build/campaign133/core-accepted-science`, HEAD
`e171caa60cce676c6c4238dd671e7c527f2bd977`.

Локаторы на этой ревизии:

- `src/numerics/reference/kernels.jl:531`,
  `product_integration_hilbert_matrix`: exact linear-panel PV и два zero
  ghost nodes; nonuniform geometry допустима, но строится dense `H`.
- Там же `:571`, `product_integration_hilbert_transform`: `wᴱ` проверяется,
  но **не умножается повторно** на PI operator.
- `src/numerics/optimized/production_cache.jl:765`: uniform branch выбирает
  `_ProductionPIFFTHilbertPlan`, nonuniform branch — dense PI matrix.
- `src/numerics/reference/grids.jl:24`: builder сейчас создаёт uniform ε.
  Наличие nonuniform transform не означает готовый adaptive solver.
- `src/numerics/optimized/production_cache.jl:266`, `build_shift_plan`:
  positive two-point interpolation для `X(E+δ)`, zero outside, без wrap.
- `src/numerics/reference/periodicity.jl:49`,
  `build_conservative_shift_pair`: weighted-adjoint finite-volume pair,
  **piecewise constant**, другой discretization identity. Нельзя использовать
  его conservation theorem для nodal-linear shifts.

Прочитаны `nessi-analysis.md` и
`../pi-causality-window-certificate.md`. Повторять PI-versus-nodal pilot и
Matsubara implementation не предлагается. Existing numbers из первичных
`pi-frozen-control/attempt01/summary.toml` и
`adjacent/contour-oracle/attempt01/summary.toml`:
real PI `N/target=1.0043381107227853`, off-axis M3/1024
`1.0045648605554407`; разность `+2.26749832655e−4`.
512→1024 off-axis изменение `−4.4221953e−8` — refinement witness,
не строгий tail bound. Рецензия — `../reviews/contour-oracle-result.md`.

Отрицательный QuadGK опыт `../adaptive-frozen/report.md` относится к
**другому interpolant: total stored ΣR вместе с embedding**, source S01
iterate32 и μ=8.519932576353691 meV. Его 16369 общих E/180059 inverses,
estimated error 0.00434553 и недостигнутая цель 1e−5 сохраняются.
Нельзя сравнить его charge непосредственно с FD64 frozen PI и объявить
улучшение. Полезные уроки: считать cumulative solves; shared-k feature union
дорога; GK-панели пересекают множество изломов source interpolant.

## Что берётся из первичных работ

1. Carrier–Greengard–Rokhlin, *A Fast Adaptive Multipole Algorithm for
   Particle Simulations*, SIAM J. Sci. Stat. Comput. **9**, 669–686 (1988),
   [авторский PDF](https://math.nyu.edu/~greengar/cgr_88.pdf),
   Remark 1.1 и §§2–3: Cauchy/Hilbert matrix application через разделение
   near/far и multipole expansions; adaptive FMM допускает clustered points.
   Наша continuous linear-panel версия и приведённый ниже matrix bound —
   самостоятельная специализация, не готовый QCL routine из статьи.
2. Greengard–Rokhlin, *A Fast Algorithm for Particle Simulations*, JCP
   **73**, 325–348 (1987), [DOI](https://doi.org/10.1016/0021-9991(87)90140-9),
   [первичный текст](https://statistics.uchicago.edu/~lekheng/courses/302/classics/greengard-rokhlin.pdf),
   §2: аналитические multipoles и translations. Для первого Julia prototype
   достаточно более простого tree traversal `O(M log N)`, без полного FMM.
3. Helsing–Ojala, *On the evaluation of layer potentials close to their
   sources*, JCP **227**, 2899–2921 (2008),
   [DOI](https://doi.org/10.1016/j.jcp.2007.11.024),
   [первичный preprint](https://citeseerx.ist.psu.edu/document?doi=af56486c8a8a727cc522590913625b4667677f3f&repid=rep1&type=pdf):
   fast far-field evaluation не заменяет специальную near-singular quadrature.
   Здесь near field проще: Γ точно linear, и PV интегрируется аналитически.
4. Barnett–Magland–af Klinteberg, SISC **41**, C479–C504 (2019),
   [1808.06736v2](https://arxiv.org/abs/1808.06736v2), §1 equations (1)–(4):
   NUFFT ускоряет exponential sums. Это не формула nonperiodic PV с нашими
   ghost ramps. NUFFT плюс Fourier multiplier потребовал бы отдельного
   free-space kernel split, padding/cutoff и near correction. Его преимущество
   перед уже exact uniform PI FFT здесь не установлено.

Wavelet compression также математически допустима: Beylkin–Coifman–Rokhlin,
CPAM **44**, 141–183 (1991),
[первичный PDF](https://amath.colorado.edu/faculty/beylkin/papers/BE-CO-RO-1991.pdf),
§§2,4–5, рассматривают быстрые Calderón–Zygmund operators. Для текущих
linear hats явное Cauchy expansion даёт более простой error budget и не
требует менять положительное nodal representation. Это выбор сложности
реализации, не измеренное сравнение скорости методов.

## Exact causal model и численно приближённый evaluator

Пусть `x₀ < x₁ < … < x_N < x_{N+1}`, `Γ₀=Γ_{N+1}=0`,
`Γ_j=Γ_j†≽0`. На каждом panel `[x_j,x_{j+1}]` Γ̃ — convex linear
interpolant, снаружи support — ноль. В одной energy unit:

    Σ(z) = ∫ Γ̃(x) / [2π(z−x)] dx,      Im z>0,
    Im Σ(z) = −Im z ∫ Γ̃(x)/[2π|z−x|²] dx ≼0,
    ΣR(E) = Λ(E)−i Γ̃(E)/2.

Γ имеет размерность энергии, moments ниже — energy², Σ/Λ — energy.
Обозначение `Im` матричное, `(Σ−Σ†)/(2i)`. Для неортогонального basis
потребовалась бы другая resolvent normalization; actual basis здесь
ортонормирован. Finite_chain4 и fixed h/T± сохраняются.

**Causal объект — этот интеграл.** Truncated far expansion является его
приближённым вычислением; оно само по себе не является новой доказанной
Herglotz representation. На real axis Γ̃ и `−iΓ̃/2` задаются точно по hats,
Λ собирается Hermitian при Hermitian input. Это гарантирует sampled
anti-Hermitian sign, но global analyticity/coherent KK требует контроля
ошибки evaluator. На complex targets можно утверждать sign лишь с учётом
bound, например если `λmin(−Im Σcomputed)>error_bound`.
Нельзя назвать нулевой sampled causality defect доказательством exact KK.

Ghost nodes фиксируются **как часть модели**, отдельно от target grid.
Refining target energies не сдвигает ghost strips. Перестроение Core default
ghost distance из нового adjacent spacing меняет support/tail model — это
новая identity, а не эквивалентная оптимизация.

## Алгоритм без плотной H

Source panels организуются в binary interval tree; leaf содержит ≤16 panels.
На сильно graded mesh нужен геометрически сбалансированный interval tree
и контроль соседних scale ratios. Случайная median tree и произвольные
clusters не дают автоматически обещанную `O(M log N)` оценку. Первый frozen
source uniform; будущий shared adaptive mesh можно делать 2:1 graded.

Для node `B=[c−r,c+r]` хранятся normalized matrix moments

    μ_q(B) = ∫_B Γ̃(x) [(x−c)/r]^q dx/(2π), q=0,…,p−1,
    m_B = ∫_B tr Γ̃(x) dx/(2π),
    F_B(z) ≈ [1/(z−c)] Σ_q μ_q(B) [r/(z−c)]^q.

Leaf moments интегрируются аналитически для linear Γ; parent moments
переносятся binomial translation с нормированными coordinates. Geometry и
translation coefficients reusable; matrix moments обновляются при каждом
изменении Γ. Нельзя переиспользовать старые moments после SCBA update.

При `d=|z−c|>r`, `q_B=r/d` geometric series даёт собственный явный bound:

    ‖Σ_B(z)−F_B(z)‖₂ ≤ m_B q_B^p/(d−r).

Доказательство: remainder ядра ограничен
`q_B^p/(d−r)`; для PSD density `‖Γ̃‖₂≤trΓ̃`; затем triangle inequality.
В bound уже присутствует `1/(2π)` через m_B. Если raw Γ не PSD, заменить
m_B на независимую upper bound `∫‖Γ̃‖₂/(2π)`; PSD claim тогда unavailable.
Float64 moment/translation/near-sum errors этим exact-arithmetic bound не
покрываются. Их хранить отдельно, а не приписывать multipole theorem.

Traversal принимает far node только при `q_B≤1/2` и достаточном order;
иначе открывает children. Можно выбрать p из локальной absolute allowance
`τ_B=τ_far m_B/m_root`: сумма allowances по любому disjoint accepted far
partition ≤τ_far. Если `p_required>pmax`, node открывается; leaf считается
direct, а не повышается tolerance. Cancellation требует **absolute** budget:
деление на почти нулевую Λ недопустимо.

Near panels вокруг real target E сначала объединяются в непрерывный
интервал `[a,b]`, содержащий E внутри, и считаются subtraction формулой

    PV∫_a^b Γ̃(x)/(E−x) dx
      = Γ̃(E) log[(E−a)/(b−E)]
        + ∫_a^b [Γ̃(x)−Γ̃(E)]/(E−x) dx.

На touching panels residual linear и removable. В частности, для panel
`[l,r]`, `B=(Γ_r−Γ_l)/(r−l)` при нетouching target:

    ∫_l^r Γ̃(x)/(z−x) dx
      = [Γ_l+B(z−l)] log[(z−l)/(z−r)] − B(r−l).

На node нельзя независимо складывать два бесконечных endpoint logarithms;
их отменяет grouped subtraction. Для exterior support edge с Γ̃=0 нужен
аналитический limit `u log u→0`. Для complex z использовать principal
`log(z−l)−log(z−r)` в upper half-plane, контролируя branch, а при больших
distance — moments/`log1p` и compensated summation. Вставка η ради обхода
PV меняет задачу и здесь не предусмотрена.

Псевдокод Julia, **24 строки**, names helpers обозначают новый proposed
контракт, а не существующий Core API:

```julia
function causal_panel_pi(Esrc, Γ, Etgt; ghosts, τΛ, pmax=48)
    check_sorted_finite(Esrc, Etgt); check_hermitian_psd(Γ)
    panels = linear_panels(Esrc, Γ; ghosts) # no extra wᴱ multiplier
    tree = balanced_interval_tree(panels; leaf_panels=16)
    moments!(tree, panels, pmax)            # exact panel integrals
    mass!(tree, panels)                    # ∫trΓ/(2π), not ‖∫Γ‖
    out, bound = similar_targets(Γ, Etgt), zeros(length(Etgt))
    for (i, E) in pairs(Etgt)
        near, queue = grouped_near_and_root(tree, E)
        Λ = analytic_subtracted_pv(panels, near, E) / (2π)
        err = 0.0
        while !isempty(queue)
            B = pop!(queue); overlaps(B, near) && (split_except!(queue,B,near); continue)
            d, q = abs(E-B.center), B.radius/abs(E-B.center)
            p = required_order(B.mass, q, d-B.radius, τΛ*B.mass/tree.mass)
            if q ≤ 0.5 && p ≤ pmax
                Λ += evaluate_moments(B, E, p); err += B.mass*q^p/(d-B.radius)
            elseif isleaf(B)
                Λ += analytic_regular_panels(B, E) / (2π)
            else; append!(queue, children(B)); end
        end
        out[i], bound[i] = Λ - 0.5im*linear_value(panels,E), err
    end; return out, bound
end
```

`required_order` должен возвращать sentinel/open для `d≤r`, q≥1,
zero mass и underflow/overflow cases, не считать `log(q)` без проверки.
Нулевая PSD mass означает нулевой node и может быть пропущена. Returned
bound — только far truncation; near/roundoff allowance добавляется отдельно.
Для production — batches matrix columns, reuse plan и workspace, никаких
`Vector{Matrix}` allocations в target loop. BigFloat проверяет маленькие
near formulas, а не выполняет весь actual transform.

Tolerance transform связывается с observable через resolvent, а не с
желаемым residual по названию. Если `Dhat=zI−Hchain−I9⊗Σhat`,
`g=‖Dhat⁻¹‖₂`, `δ=‖Σhat−Σ‖₂` и `gδ<1`, то resolvent identity даёт

    ‖Ghat_c−G_c‖₂ ≤ g² δ/(1−gδ),
    ‖Ahat_c−A_c‖₂ ≤ 2g² δ/(1−gδ).

Здесь **нет factor9**: norm repeated block diagonal δΣ равна norm δΣ.
Нужен full-chain inverse norm, не только norm central5. Floating smin/LU
estimate — practical witness с отдельным backward-error allowance,
не interval certificate. Для положительных observable weights wα сумма
`Σα wα fα gα² δα/[π(1−gαδα)]` ограничивает transform вклад в
`‖δF‖F/√5`. В frozen pilot потребовать этот вклад ≤1e−7 либо сообщить,
что bound недостаточен. При слишком большом bound переоценить Λ с меньшим
τΛ; это не новые Green abscissas, но wall cost учитывается. Малый local
Γ eigenvalue или узкий resonance могут сделать исходный τΛ недостаточным.

## Как будущая адаптация связана с KMS и LO

Не вводить cubic/high-order Γ interpolation ради внешнего вида точности:
положительность matrix hats имеет прямую causal опору. Refinement добавляет
узел с barycentric Γ от прежних neighbours: **до следующего физического
map update** continuous Γ̃ и Σ не меняются. Coarsening обычно их меняет;
нужен отдельный representation error, а не label exact optimization.

Для shared μ/T equilibrium на nonuniform mesh shift должен вычислять

    A_s = positive_linear(A, E+sΩ),
    Gn_s = f(E+sΩ) A_s + interpolate(δGn, E+sΩ),
    Gp_s = [1−f(E+sΩ)] A_s + interpolate(δGp, E+sΩ),
    δGn = Gn−fA,     δGp = Gp−(1−f)A.

В equilibrium residuals нулевые; Fermi/Bose identities выполняются в
**physical shifted energies**, поэтому fractional Ω/ΔE само по себе не
нарушает LO KMS algebra. Это новая shift variant. Вне equilibrium residuals
могут нарушать PSD суммы; clipping и forced-FD projection не являются частью
метода. Zero-outside policy применяется одинаково к A и residuals; edge loss
остаётся измеряемой ошибкой domain coverage. Нельзя превращать PI ghost
taper local Γ в ghost extrapolation электронного A без отдельного решения.

Это решает один вопрос согласованности, но **не** доказывает discrete
particle conservation: на произвольной nonuniform сетке positive nodal
`S+`/`S−` обычно не удовлетворяют `W S+ = S−† W`. Native conservative
finite-volume pair решает weighted-adjoint вопрос для другого representation
и не даёт одновременно pointwise KMS для произвольного analytic f.

Перед full adaptive SCBA потребуются на том же варианте:
KMS residual с occupied absolute floor; min eigenvalues Σin/out; weighted
shift-adjoint defect; signed particle collision integral с absolute event
scale; LO energy flow к bath. Если nodal KMS-aware map не проходит balance,
следующий **отдельный** вариант — common-refinement event quadrature на
`mesh ∪ (mesh±Ω)` с одной положительной quadrature measure для forward/reverse
transition pair. Detailed balance там проверяется в event energies:
`NΩ f(E)[1−f(E+Ω)] = (NΩ+1) f(E+Ω)[1−f(E)]`.
Conservative deposition и pointwise nodal reconstruction при этом надо
вывести совместно; здесь готового theorem/API для них не предлагается.
Таким образом план не обещает, что один fast Hilbert operator автоматически
сделает full nonuniform LO solver физически согласованным.

## Стоимость и честные ограничения

N — source knots, M — cumulative distinct targets, C — число matrix columns,
b=16 leaf panels, p≈32–48. Geometry memory `O(N/b)`;
moments `O(p C N/b)`, Γ/output `O(C(N+M))`. Build имеет
`O(p C N + p² C N/b)` arithmetic с обычными translations; balanced tree
traversal — ориентировочно `O(p C M log N)` плюс analytic near panels.
Это treecode, поэтому `O(N)` FMM theorem не приписывается реализации.
При degenerate grading/admissibility cost может ухудшиться; логировать
accepted nodes, opened leaves и direct panel evaluations.

Для N=30721, Nb=5, Nk=11, C=275: один ComplexF64 Γ ≈129 MiB;
все moments при p=40, b=16, ~3840 nodes ≈645 MiB. Column chunk=16
уменьшает moments до ≈37.5 MiB; retained full Γ+Λ тогда около258 MiB
плюс buffers/geometry. Это **transform-only estimate**, без SCBA histories,
kernel tensors, HDF5 buffers и Julia/JIT; не обещание RSS всей программы.
Дерево не требует dense H, которая здесь занимала бы ≈7.0 GiB на одну
Float64 matrix, и не требует такой же verification temporary.

На uniform source/targets exact PI FFT почти наверняка проще и должен
оставаться штатным путём; treecode crossover не измерен. Для nonuniform
M~N tree traversal имеет большие constants. Экономия возможна, если
adaptive M существенно меньше uniform doubled N и/или grid representation
потребует меньше всех family arrays. Снижение transform memory отдельно
от снижения time-to-accepted-SCBA. Wall time в секундах без prototype
предсказывать нельзя. Обоснованный operational ceiling первого frozen
прототипа — 120 s CPU/wall stop, 1.5 GiB RSS, а не прогноз его скорости.

## Один bounded frozen oracle: разрешить discrepancy на одном k

Предлагаемый experiment identity:
`FD64-PI-frozen-paneltree-singlek-v1`. Не запускался этой задачей.

**Модель фиксирована.** Входы берутся из существующего
`adjacent/contour-oracle/attempt01/input.toml`: source
`e6fffb997f36397f877aee1e2c98b1dbad7d97417d45458a042b5f074ce57a6d`,
aggregate local Γ SHA
`edb8530e9a59944afe1a447d999bf0c6c42d29d726e16144a5d85bdccdaecc55`,
h SHA `8cde8ea024d9e690bf0396428574496734d7bad5962fc00617b9b380ae1194a7`,
μ=0.008495393439818052 eV, T=70 K, original30721 hats и original ghosts.
Никакого thinning/coarsening Γ, μ fit, fresh scattering update, mixed
embedding или η. Новые adaptive **targets** не меняют source representation.

1. Из saved full matrices выбрать один k с максимальным
   `g_s w_k |tr(ρM3_1024−FrealPI)|/target`; при tie меньший index.
   Сохранить также matrix discrepancy и occupation scale. Selection по
   existing evidence, не поиск красивого результата. Использовать именно
   saved central5 density на этом k, без factor9. Это one-k question,
   не воспроизведение sheet density всех11 k.
2. До Green evaluations проверить tree Λ на ≤64 заранее заданных targets:
   endpoints, knots/midpoints в occupied features, near-null Λ cancellation
   и ghost edges. Oracle — independent direct analytic **panel integral**
   streamed через все original panels, не dense H и не another tree.
   Цель `max‖ΔΛ‖F ≤ 1e−10 E0 + 1e−9 max‖Λdirect‖F`; near-zero часть
   absolute. Это evaluator diagnostic, не новый SCBA tolerance. Если
   arithmetic allowance его не поддерживает, stop/negative result.
3. На том же k интегрировать central A/(2π) и fA/(2π) с adaptive target
   panels и cached queries, первоначально ≤128 panels по saved Γ/A features
   и равномерному coverage; boundaries coalesce по явному energy tolerance.
   Разделять panels у крупных source slope changes; использовать local
   Γ curvature **и** sampled smallest singular value Dyson denominator
   как refinement hints. Не объявлять эти hints доказательством отсутствия
   скрытого резонанса. Source breakpoints не разрешено молча выкинуть.
   Конкретная integration policy: GK7/15 на каждой target panel, норма
   local embedded difference `max(‖δM‖F,‖δF‖F)/√5`; делить пополам panel
   с наибольшей оценкой, переиспользуя cached energies. До bisection
   предпочесть ближайший source knot с крупным slope jump, если он лежит
   в средней трети panel; ties по меньшей энергии. Source knots с малыми
   jumps остаются внутри panel, а их сохранённая density участвует точно
   в Λ queries. Сумма local estimates должна стать ≤1e−5; это отдельный
   практический gate, не rigorous quadrature bound. Заранее заданные
   physical feature hints не меняются по желанию подогнать Matsubara match.
4. Cumulative Green target budget **8192** на выбранный k, including
   discarded panels и проверочные повторения; wall120 s, OS MemoryMax1.5GiB,
   swap0, Julia1/BLAS1, output≤16MiB. Остановка штатным timeout/scope с
   checkpoints и итогом `budget_reached`, если задача не решена. Столько же
   inverse calls на all11 k не запускать по этому протоколу.
5. Раздельно сохранить `M=∫A/(2π)`, `F=∫fA/(2π)`, errors embedded rule,
   cumulative target count, transform bounds и solves. Проверить
   `‖F−ρM3_1024‖F/√5 ≤ 1e−5`, а также occupied
   `‖F−ρM3_1024‖F / max(‖ρM3_1024‖F,1e−8) ≤ 1e−4`;
   сообщить обе нормы, не только trace. Matsubara512→1024 matrix difference
   и её вклад в эти числа привести отдельно. Эти diagnostic thresholds
   различают known ~1e−4 discrepancy и не принимают baseline.

Если feature seeds уже требуют больше8192 cumulative targets, опыт
заканчивается отрицательно до попытки полного adaptive solve. При budget
stop сохранить обе причины отдельно: embedded quadrature estimate и
независимую measured matrix difference. Agreement с Matsubara при большой
embedded estimate — полезное наблюдение, но по этому протоколу accuracy
gate ещё не достигнут. Это предотвращает принятие случайной cancellation.

Matsubara oracle включает атомарную меру; real A integration — нет.
Поэтому отсутствие exterior/interior atoms и ghost contribution для
выбранного k остаётся явной предпосылкой интерпретации. Переиспользовать
existing edge/ghost witnesses из contour artifact и проверить применимость
interior-PD условия из certificate. Если оно недоступно, failure matching
ρ нельзя автоматически назвать quadrature error; outcome
`density_difference_unattributed`. Никакого forced mass normalization.

Главное различающее свидетельство: совпали ли **полные F matrices** одной
continuous causal модели с off-axis result при ≤8192 cumulative targets,
когда original real PI давал видимую разность. Это информативнее сравнения
двух uniform quadratures. Успех показывает только возможность уменьшить
число **одно-k frozen observable** solves; перенос на shared all-k grid и
SCBA не доказан. Недостижение bound/cap или Matsubara match — содержательный
отрицательный итог; uniform doubling не объявляется побеждённым.

## Решение

Этот кандидат заслуживает маленького frozen прототипа, если кампания
разрешает следующий опыт: явный bound, positive-hat source и независимая
плотность уже имеются. Сейчас техническое решение — сохранить exact FFT
для uniform grids; рассматривать panel tree как замену **nonuniform dense
fallback**, а adaptive SCBA — как последующую отдельную научную задачу.
Ни новая модель, ни сходимость, ни выигрыш времени этой запиской не приняты.
