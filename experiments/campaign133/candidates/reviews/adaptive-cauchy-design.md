# Independent math review: adaptive Cauchy design

Проверен `adjacent/adaptive-causal-energy.md`, SHA256
`0eccec800e139fe321f64f92edeb418c3aa37339883c5eed88b5d98e11cbea7d`.
Accepted Core HEAD `e171caa60cce676c6c4238dd671e7c527f2bd977`, tree
`742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8`. Здесь только independent
proof/design review: Julia/Pkg/массивы/Green/SCBA не запускались,
candidate/Core не менялись. Использован physics-audit.

**Вывод:** exact panel/PV algebra, matrix far remainder и resolvent
propagation правильны. Кандидат может стать honest error-budgeted
nonuniform evaluator **той же positive-hat модели**. Но returned `bound`
сейчас по design ограничивает только far truncation; полного floating
error certificate ещё нет. До prototype contract надо явно разделить
real PV и complex Cauchy paths и определить область density integration
для сравнения с Matsubara. Это содержательные qualifications, не обнаружение
неверного far-kernel theorem и не acceptance нового adaptive SCBA.

## 1. Что математически подтверждено

Для ξ=x−c и |ξ|≤r<d=|z−c| точный kernel remainder после p terms:

```
1/(z−x) − sum_{q=0}^{p−1} ξ^q/(z−c)^(q+1)
    = [ξ/(z−c)]^p /(z−x).
```

Поэтому его модуль ≤(r/d)^p/(d−r). Для matrix PSD Γ̃,
||Γ̃||2≤trΓ̃; integration/triangle inequality дают ровно

```
||ΣB−FB||2 <= mB (r/d)^p/(d−r),
mB=integral_B trΓ̃/(2π).
```

Phase — positive powers r/(z−c), не conjugate и не |z−c| в самом
evaluator. Абсолютная величина d нужна только bound/admissibility.
Moments имеют energy², F — energy; дополнительного h или wE нет.
Noncommuting matrix Γ не мешают: scalar kernel не меняет matrix ordering.

Для disjoint accepted far nodes B сумма mB≤mroot, следовательно
allocation τB=τfar mB/mroot действительно даёт суммарный truncation
budget≤τfar. Это требует настоящего partition: near, accepted far и
direct leaves должны покрывать каждый source-panel ровно один раз.
Zero mass при true PSD означает zero matrix density a.e.; node можно
пропустить. При non-PSD input этот shortcut/trace majorant недопустим.

Для near panel с slope B формула candidate
`[Γl+B(z−l)] [log(z−l)−log(z−r)]−B(r−l)` имеет правильный sign.
На real near interval [a,b] subtraction term
`Γ(E) log[(E−a)/(b−E)]` также правильный. Оба singular endpoint log при
target knot должны отменяться **аналитически до evaluation**, а не как
Inf−Inf. На outer zero-Γ ghosts использовать finite x log|x| limit.

## 2. Panel-moment и PV oracles для будущей реализации

Parent translation стоит зафиксировать явно. Для child center cs/radius rs
и parent center c/radius r:

```
mu_q(parent)=sum_children sum_{j=0}^q binomial(q,j)
    *((cs−c)/r)^(q−j)*(rs/r)^j * mu_j(child).
```

Leaf panel [l,u], L=u−l, Γ(x)=Γl+θΔΓ, a=(l−c)/r, b=L/r:

```
mu_q(panel)=L/(2π) sum_{j=0}^q binomial(q,j) a^(q−j)b^j
    *[Γl/(j+1)+ΔΓ/(j+2)].
```

Это identity, **не рекомендация** наивно суммировать cancellation-prone
binomial expression при p48. Stable power-difference/recurrence и его
roundoff allowance нужны отдельно. Для PSD source ||mu_q||2≤mB —
дешёвый sanity witness, не полный forward-error proof.

Независимые small oracles до actual frozen queries:

* asymmetric single hat с peak Γ и widths hL,hR:
  `Lambda(peak)=Γ log(hL/hR)/(2π)`. Nonuniform diagonal PV **обычно не0**;
  symmetric hat даёт0. Это сильная защита от заимствования uniform k0.
* panel zeroth moment `L(Γl+Γu)/(4π)` и first moment
  `L[(2l+u)Γl+(l+2u)Γu]/(12π)`; sums должны совпасть с root moments.
  Проверить parent-vs-direct moments и partition mass ledger.
* complex-Hermitian matrix hats, target knot±tiny offset, midpoint,
  exterior ghosts, strongly graded gaps; reference — BigFloat analytic
  integration с независимым expression/series и manufactured scalar cases.
* high-z `zΣ→S0`, `z²Σ−zS0→S1`; Schwarz adjoint relation при extension
  и знак −ImΣ на upper-half-plane points с явным error allowance.
* uniform limit против existing exact uniform PI operator. Само согласие
  двух paths, использующих одну panel formula, не исключает shared error.

## 3. Что обязательно уточнить в error contract

**A. Полный budget.** Строки144–160/212–218 candidate честно оговаривают,
что bound только analytic truncation. Его нельзя затем подставить как
полную δ в resolvent theorem. Нужен ledger для каждого target:

```
δΣ <= sum_far[trunc_B + (1/dB)sum_q(rB/dB)^q e_muBq + e_evalB]
      + e_near + e_final_sum + e_geometry/input.
```

e_mu включает leaf construction и **накопленные translations**. Обычная
compensated summation уменьшает ошибку, но сама не даёт verified upper
bound. Floating m/r/d требуют upper/lower one-sided allowances, иначе
округлённый `m*q^p/(d−r)` не доказанный upper bound. В practical prototype
можно сохранить conservative estimates/empirical direct comparisons и
называть их numerical witnesses; formal interval certificate не обязателен
для exploratory comparison, но его отсутствие должно остаться в выводе.

Fixed Float64 source knots/Γ можно объявить точными числами frozen model.
Их coarsening, source-node coalescing, смена ghosts или reuse stale
moments — уже representation error. Merging target coordinates по energy
tolerance также требует evaluation/geometry allowance; near resonance
«близкие E» не автоматически дают близкие G. Не включать этот эффект
молча в far theorem. Near/far partition меняет evaluation, не source.

**B. Два target contracts.** Псевдокод186–208 соответствует **real** E:
PV Λ, затем ровно один `−iΓ̃(E)/2`. Для complex z нужен ordinary analytic
near-panel Cauchy integral (с principal log difference при Imz>0), без
добавочного `−iΓ̃(z)/2`/PV; Γ̃(z) вообще не определена как spectral density.
В real limit panel log содержит `−iπΓ̃(E)` до division2π. Нельзя взять
этот полный complex-limit integral и затем второй раз добавить imaginary
term. Regular real panels вне near должны действительно не touching E.

**C. Propagation относится к controlled object.** Для полного45 Dhat
g=||Dhat^−1||2, δ=||Σhat−Σ||2, gδ<1 доказательство resolvent identity
даёт `||δGc||2≤g²δ/(1−gδ)` и `||δAc||2≤2g²δ/(1−gδ)`.
Factor9 не нужен: norm I9⊗δΣ равна local operator norm. Central5 inverse
norm нельзя использовать вместо full45 g. При verified γmin(E)>0 можно
дешёво использовать g≤2/γmin(E), но такой bound может оказаться бесполезен;
computed positive eigenvalue без allowance не verified lower bound.

Для **одной и той же positive quadrature sum** веса wα дают
transform contribution к normalized Frobenius density ≤
`sum wα fα gα²δα/[π(1−gαδα)]`, как написано. Это не bound continuous
integral между queried points. Отдельно остаются real-axis quadrature,
inverse/backward/accumulation errors, atoms/ghosts и Matsubara tail.
Embedded GK difference не превращает их сумму в rigorous bound.

## 4. Causality и one-k oracle

Positive-hat **точный интеграл** causal/Herglotz и задаёт normalized
spectral measure finite-chain system. Tree approximant может менять
partition с target, не обязан быть глобально analytic/Herglotz. Real
Hermitian Λhat с exact −iΓ/2 сохраняет sampled dissipative sign и Ahat PSD
при invertibility; это не proof exact KK/spectral mass. Candidate правильно
разделяет эти утверждения; нельзя убрать qualification после удачного run.

One-k comparison имеет корректную model identity: same local Γ/h/T,
same finite4, fixed source μ/T, actual selected k и Γ ghosts. Однако
integration domain в prototype input нужно определить явно:

* интегрировать **весь** continuous support [outer-left-ghost,
  outer-right-ghost], затем отдельно qualified absence of atoms; либо
* интегрировать original physical [a,b] и сравнить с Matsubara только
  после ghost contribution correction/bound, достаточно малого для gate.

Вне full support continuous A=0, но possible delta atoms остаются.
Matsubara включает их; real A quadrature — нет. Если их отсутствие и
ghost accounting не qualified, outcome `density_difference_unattributed`,
а не «quadrature failed». Matsubara512→1024 difference — refinement
witness, не rigorous tail error. Selection по максимальному saved weighted
trace discrepancy допустим для одного конкретного k, без all-k claims.

Таким образом ≤8192-target comparison может различить возможность
экономной integration **frozen one-k observable**, если source/ghost/atom
conditions, transform contribution и independent matrix-density agreement
предъявлены раздельно. 64 direct-transform probes проверяют реализацию,
но не гарантируют все future targets; per-query ledger остаётся нужен.
Совпадение полной F matrix при большом GK estimate не принимает accuracy
gate; negative budget outcome также сохраняется. Performance и SCBA gain
остаются not_measured. Следующая реализация не разрешена этим review.

## Первичная опора

[Carrier–Greengard–Rokhlin, SIAM SISC9 (1988),669–686,
Remark1.1/§2](https://math.nyu.edu/~greengar/cgr_88.pdf) действительно связывает
adaptive multipoles с Hilbert-matrix action и translations. Он не даёт
готовый continuous matrix-hat floating error budget: bound выше выведен
независимо. [Helsing–Ojala, JCP227 (2008),2899–2921,
author preprint](https://lup.lub.lu.se/search/files/5505653/3878573.pdf)
обосновывает специальное обращение с near-source quadrature. Ссылки
проверены primary browsing; paper authority не заменяет наши PV/PSD tests.
