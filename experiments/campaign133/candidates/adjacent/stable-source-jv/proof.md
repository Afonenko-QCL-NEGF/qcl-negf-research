# Stable-source Γ tangent: ограниченное доказательство и квалификация

**Scope.** Write-only prototype, не Newton solver и не новый SCBA run. Base —
unprojected stableA-cont64 Γlast64 `cab81029…db55`, canonical coarse m939 grid,
Core e171/tree742, uniform PI, hH во всех девяти cells, η=0. Проекции Γ/Φ нет.
Предполагаются finite inverses, PSD-interior base, nonnegative raw q,
`0<Ntarget<capacity`, `χ>0`. Global strict-PD proof ещё не получен: прежний
cone allowance сам не доказывает smooth interior. Внутренние derivative
directions могут быть indefinite; они не являются trial physical Γ.

## 1. Exact covariance derivative для фактической finite4 recursion

Рабочее пространство Hermitian Γ — вещественное. `δR=Herm(HPI[v])−iv/2`:
PI и Herm — real-linear, hH и zero-ghost geometry фиксированы. Для стороны
`C=T+` или `T+†`, начиная с `g=L=δg=δL=0`, четыре раза:

```
V = R + C g C†;          Lin = iΓ + C L C†;
gnew = [E I−hH−V]⁻¹;    Lnew = gnew Lin gnew†;
δV = δR + C δg C†;       δLin = i v + C δL C†;
δgnew = gnew δV gnew;
δLnew = δgnew Lin gnew† + gnew δLin gnew† + gnew Lin δgnew†.
```

Это прямая производная native `_cavity_side!`; только local source lesser
`iΓ`, greater0 — spectral probe, не physical occupation. Actual V0 sparse
remaps проверяются как identity. После четырёх cells embedding каждой
стороны есть `C(g,L)C†`, и derivative так же sandwich. В центре:

```
G = [E I−hH−R−ΣembR]⁻¹; K = Γ−iΣembL;
δG = G(δR+δΣembR)G;      δK = v−iδΣembL;
A = G K G†;
δA = δG K G† + G δK G† + G K δG†.
```

Первый inverse derivative содержит G справа, **не G†**. δA берётся именно
от covariance; `i(δG−δG†)` в реализации отсутствует. Native hH/ΛHerm
фиксируются, а их малые roundoff-correction witnesses наследуются от
прежнего оператора. Различие primal G/A recomputation от qualified streaming
base сохраняется в каждом Jv; это необходимая численная gate, не theorem о
bitwise derivative floating program. Перестановка небольших products может
дать roundoff, хотя mathematical map тот же. Per-worker5×5 scratch bounded;
full45 Green/derivative на energy grid не строится.

Независимый tiny oracle собирает complex noncommuting full9-chain (18×18
и45×45), прямо вычисляет `δGfull=Gfull(I9⊗δR)Gfull`, covariance derivative и
full-matrix central difference. Это не shared cavity reference.

## 2. Текущая μ, signed native scattering и raw response

При `w=gs wE wk/(2π)`, `f=f(E−μ)`, `p=1−f`,

`χ=Σ w f p Tr A/kBT`, `δμ=−Σ w f Tr δA / χ`,
`δF=fδA+(fp/kBT)A δμ`, `δGL=iδF`, `δGG=−i(δA−δF)`.

Тогда `δN=Σw TrδF=0`. Другая μ или derivative числа bisection steps дали
бы другую map. Сохраняются χ, числитель, δμ, residualδN и абсолютная шкала
из сумм модулей; χ≈0 означает обусловленность/остановку, не division floor.
`δN_relative<=1e−10` — inherited arithmetic gate, не scientific tolerance.

Four static channels дают `Ks[δA]` прежними full kernels/qK/wk; native
`_production_lo_contraction` применяется к `(δGL,δGG)` с теми же fixed Bose,
zero extension и actual shifts. Поэтому `Jv=Σstatic Ks[δA]+i(δΣLO>−δΣLO<)`.
LO contraction linear в этих полях; никакой differentiation нового kernel,
Bose, screening или Poisson здесь нет. `phi` читает только A/GL/GG; GR slot
signed tangent storage явно dummy, он не заявляется физическим GreenState.
Tiny проверяет δN0 и линейность каждого из **пяти actual native** channels.

## 3. Один различающий bounded diagnostic

Один fresh base Φ и четыре probes `Γ±εΓ`, ε=1e−3,5e−4: при true PSD Γ
scaling положителен и не требует eigen clipping. Для finite ε нужны все
те же raw q/current μ/domain gates. Central differences сравниваются с
analytic `JΓ`; их взаимная разность и ε² scaling — evidence, не rigorous
error bound: μ-fit, FFT/inverse arithmetic и curvature учитываются отдельно.
Не уменьшать ε автоматически до cancellation и не называть truncation
error1e−10 mandatory. Каждая начатая failed Φ/Jv расходует свой cap.

Остаётся **три** Jv (всего≤4) для real weighted Arnoldi из
`Herm(Φraw−Γ)`. Удалённая anti-Hermitian часть seed записывается; raw residual
не заменяется working seed. Two-pass reorthogonalization даёт small Hbar;
Ritz residual, `σ(Brestricted)`, `Brestricted=Ibar−Hbar` и out-of-span term
сохраняются. Small `σ` даёт найденное near-null direction; большой `||Jv||`
не доказывает unstable eigenvalue. Nonnormality, modes вне span и validation
error исключают spectral-radius/full-σmin certificate. Small GMRES coordinates
решают лишь restricted `(I−J)δ≈R`; код **не формирует/не применяет** Γ update.
Ни один Krylov output не принимает SCBA/native convergence и не различает
самостоятельно edge physics от quadrature: нужен matched refinement evidence.

Первичная опора: [Saad1980 Arnoldi](https://doi.org/10.1016/0024-3795(80)90169-X),
[Saad–Schultz1986 GMRES](https://doi.org/10.1137/0907058),
[Knoll–Keyes2004](https://doi.org/10.1016/j.jcp.2003.08.010). Эти источники уже
проверены в `../gamma-linear-response-plan.md`; inverse/number/source proof
выведен выше для actual map. QCL contractivity/сходимость отсюда не следуют.
