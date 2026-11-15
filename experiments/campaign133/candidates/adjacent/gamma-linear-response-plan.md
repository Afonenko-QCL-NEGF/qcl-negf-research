# Γ-map: проверяемый linear-response diagnostic и один warm control

**Решение сейчас:** дождаться numerical outcome исправленного pilot.
Его `Dict`/report failure не является свидетельством математической
неустойчивости. Если constant-seed pilot не даёт решения, следующий
минимально отличающийся control — **одна** Γ initialization из сохранённого
S01-FD64-continuation128, перенесённая на уже выбранную m939 grid, затем
тот же Π/AA-II/native audit. Не запускать одновременно новую homotopy,
Newton/Krylov solver и другую сетку. Ниже design, не реализация/запуск.

Настоящий Jacobian thermal Γ-map существует при interior domain conditions.
Его exact action можно получить без finite differences; короткие secants
из уже имеющейся AA history дают более дешёвый, но менее сильный diagnostic.
**Один frozen-grid spectral radius не отделяет discretization error от
physical near-neutral mode.** Такие attribution требуют независимых
occupied-quadrature свидетельств и matched refinement/сохранённой symmetry.

## 1. Scope и реальные источники

Fixed V0/U0, TL=TLO70K, same basis/k/finite_chain4, static acoustic,
impurity, IFR, alloy + native thermal LO; fixed screening, no e–e.
Local aggregate unknown Γ не содержит embedding. Retarded closure —
uniform PI/full hats/zero ghosts, `elastic_retarded=:hilbert`, ηwork=0.
Это ветвь Appendix A1–A6 в `../pi-causality-window-certificate.md`,
не `direct_green` и не Poisson derivative.

Native source прочитан на accepted Core HEAD
`e171caa60cce676c6c4238dd671e7c527f2bd977`, tree
`742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8`, tracked clean.
Точки привязки: `optimized/production_operators.jl` —
`_cavity_side!`, `_cavity_embedding_production`, `_production_lo_contraction`,
`_selected_static_contraction`; `optimized/seed.jl` — exact q_e measure;
`optimized/fft_parallel.jl` — linear uniform PI operator. Current harness
математически разобран в `../reviews/gamma-math-code-v2.md`; последующие
report-only fixes требуют собственного frozen identity/diff, не меняют
автоматически этот proof. Рецензент здесь не запускал Julia/Pkg/arrays.

## 2. Domain: real Fréchet derivative, а не complex AD через всё подряд

Working space — Hermitian matrix fields как **real** vector space,
inner product Re sum wE wk Tr(v†w). Нужны finite sampled inverses,
finite A/q, 0<Ntarget<C[A], χ=dN/dμ>0 и достаточно interior PSD margin.
При cone boundary или capacity/root/inverse failure smooth Jacobian
полной admitted map может отсутствовать. Cone rejection/QR branch,
clipping и bisection iteration count не следует дифференцировать.

Продифференцировать physical smooth Π на допустимом base, number fit —
через implicit-function formula. Herm/adjoint делают derivative real-linear;
complex-step derivative `imag(Φ(Γ+iεv))/ε` здесь неприменим. Переход к
packed Hermitian coordinates требует sqrt2 weights off-diagonal; полный
matrix storage с real Frobenius inner product уже имеет эти веса.

## 3. Exact retarded/cavity tangent

Для real Hermitian perturbation v=δΓ:

```
δRlocal = HPI[v] − i v/2.
δG45 = G45 (I9⊗δRlocal) G45,
δGc = P δG45 P†,   δA=i(δGc−δGc†).
```

В первом product справа **нет adjoint G45**: это derivative inverse,
δD=−δΣ. Adjoint появляется только при δA/advanced variation.
Fixed V0 означает одинаковую local Σ в девяти cells: никаких derivatives
energy remaps/bias/hcell/T не добавляется. При V≠0 этот shortcut неверен.

Full45 не нужен на всех E,k. Пусть C+=T+, C−=T+†; для каждой стороны
g_s,0=0. От outward boundary к центру, j=1,…,4:

```
g_s,j = [E I−hcell−Rlocal−C_s g_s,j−1 C_s†]^−1,
δg_s,0=0,
δg_s,j = g_s,j [δRlocal+C_s δg_s,j−1 C_s†] g_s,j.
δΣemb = sum_s C_s δg_s,4 C_s†,
δGc = Gc [δRlocal+δΣemb] Gc.
```

Это derivative actual native recursion. В local block loop достаточно
текущих g/δg и bounded5×5 scratch; восемь full cavity Green tensors не
требуются. AD-friendly реализация должна использовать solve/tangent-solve
с теми же primal matrices, а не собирать dense global Jacobian.
При fixed V0 равенство small finite45 derivative и recursion derivative —
независимый algebra oracle на manufactured complex-Hermitian chain.

## 4. Current μ: один global rank-one член

Определим f_e, p_e=1−f_e stable complementary logistic и

```
L_f(δA)=gs/(2π) sum_e,k wE wk f_e Tr δA_e,k,
χ=gs/(2π) sum_e,k wE wk f_e p_e/kBT Tr A_e,k > 0,
δμ=−L_f(δA)/χ,
Q_e,k=(f_e p_e/kBT) A_e,k,
δF=δ(fA)=f δA+Q δμ.
```

Тогда δN=L_f(δA)+χδμ=0 — точный scalar tangent oracle на исходной
number quadrature. χ, δμ, δN absolute и scales обязательны: small χ
усиливает response, даже если raw Γ residual мал. Stale μ или δμ=0
дифференцируют другую map. Отдельный frozen-μ response полезен только
как явно помеченный decomposition.

`δGL=i δF`, `δGG=−i(δA−δF)`, поэтому i(δGG−δGL)=δA.
Four static contributions: δΓs=Ks[δA], прежние qK/wk/kernels.
Для fixed Bose N_B и actual linear zero-extended shifts S± native LO
linearization можно записать без quotient/interpolation rewrite:

```
δΓLO = KLO{ S+[N_B δA+δF] + S−[(N_B+1)δA−δF] }.
JΓ v = sum_static Ks[δA] + δΓLO.
```

Формула следует непосредственно из native Σ</> contractions; она верна
также для fractional linear shifts как derivative **их** discrete map.
Не заменять S(fA) на интерполированное f × интерполированное A.
Для full thermal native fixed-point equivalence дополнительно нужен
actual near-exact KMS на commensurate Ω/T/stencils; derivative formula
сама не устраняет fractional KMS defect.

Если T:δΓ→δA, то

```
JΓ = J_(μ fixed) − u (L_f∘T)/χ,
u = KLO[(S+−S−)Q].
```

Это global rank-one correction. Static Γ не зависит от μ непосредственно,
поскольку Ks[A]; electron–electron, hot-phonon или self-consistent screening
добавили бы terms и сюда не входят.

## 5. Что именно говорят eigenvalues/near-null witnesses

Для F(Γ)=Φ(Γ)−Γ derivative B=I−JΓ меняет знак относительно DF.
λ≈1 у JΓ означает near-null **B**, а не обязательно instability.
λ≈0 — сильное Picard contraction direction, не плохо обусловленный root.
Для linear mixing α:

```
Jmixed = (1−α)I+αJΓ,
stability of eigenmode: |1+α(λ−1)|<1.
```

Если Reλ<1, положительный α должен быть меньше
`2(1−Reλ)/|1−λ|²`; Reλ>1 нельзя исправить положительным scalar damping
в этой локальной eigenmode. Это не theorem о существовании SCBA root или
его доступности AA. Вдали от root JΓ — только local sensitivity; отказ
constant seed/cone ещё не asymptotic instability.

Γ-only FD map **исключает** independent occupation/LO-ladder populations
полной kinetic map. Поэтому её neutral mode нельзя автоматически назвать
slow redistribution по E mod Ω. При ideal LO+elastic ladders могут иметь
дополнительные stationary freedoms; exact discrete conservation при
finite window/trapezoid endpoints следует доказать отдельно. One common μ
выбирает thermal ensemble. Утверждение о full kinetic nullspace потребует
unconstrained native tangent/collision audit, которого этот plan не строит.
Number constraint удаляет одну scalar density variation, но не доказывает,
что у JΓ автоматически есть или нет λ=1.

Независимый scalar tangent oracle: для G=(ξ+iγ/2)^−1,
A=γ/(ξ²+γ²/4) должно быть
`dA/dγ=(ξ²−γ²/4)/(ξ²+γ²/4)²`. Для toy Φ=g²A на resonance root γ*=2g
имеет J=−1, хотя root существует; damping α=.2 даёт eigenvalue .6.
Это аналитическая проверка signs/interpretation derivative, не новый
QCL/PI fixture и не приёмка физических каналов.

Nonnormality важна: matrix [[.5,M],[0,.5]] имеет radius .5, но может
сильно усиливать perturbation за один шаг. Поэтому `||Jv||/||v||>1`
не доказывает ρ(J)>1; ρ<1 не исключает transient growth/cone rejection.
Короткий Arnoldi не видит не возбуждённые modes. Даже small residual
approximate eigenpair на nonnormal J даёт pseudospectral witness, а не
безусловную eigenvalue enclosure. Не сохранять Ritz maximum под именем
`spectral_radius_certified`.

## 6. Дешёвые данные сначала; Krylov только при plateau и отдельном решении

**Без extra Π:** live AA history хранит X=ΔΓ и Y=ΔΦ_work=ΔΓ+ΔF_work.
Weighted real QR **X**, не только имеющийся QR ΔF, даёт X=QR;
Z=Y R^−1, H=Q*Z. Сохранить small H/R/condition, observed Ritz values,
out-of-span defect ||Z−QH||, step/residual scales, representation error.
Это working-Hermitian secant response, усреднённый вдоль разных nonlinear
segments, не exact raw J; representation-conversion error отдельно нужен.
Крупные steps/rank loss/сильная μ variation делают такую оценку только
описательной. Не восстанавливать невыписанную history из scalar CSV.

При analytic Jv real Arnoldi m=4 даёт
JQ_m≈Q_(m+1) Hbar. Для real Schur 2×2 blocks complex conjugate pair
обрабатывается в real Hermitian span, без complex Γ trial coefficients.
Ritz residual и Jv-validation error сохраняются отдельно. Для B=I−J
small tall matrix `Ibar−Hbar` имеет singular values restricted BQ;
small minimum даёт **найденный** near-null witness v и ||Bv||, но отсутствие
малого singular value на этом span не исключает modes снаружи. СVD этой
малой матрицы предпочтительнее normal equations/домысла о полном σmin.

**FD validation не заменяет analytic oracle.** Для одного interior base
и направления v сравнить central D_ε=[Φ(Γ+εv)−Φ(Γ−εv)]/(2ε) при ε и ε/2:
четыре fresh Π, current μ в каждом, без смешивания Green/μ. На первом
control дешёвое безопасное направление v=Γ даёт PSD для |ε|<1, включая
почти нулевые tail blocks; оно не претендует на все eigenmodes. Relative
ε=10^−3,5·10^−4 — заранее объявленные probe amplitudes, не scientific
tolerance. Записать absolute/weighted D difference и analytic error,
actual number-fit error/χ, representation/FFT noise. Разность двух ε
не является rigorously certified derivative-error bound без smoothness/
roundoff bounds; нельзя автоматически уменьшать ε до cancellation.

Для arbitrary v sufficient cone condition:
`|ε| max_blocks ||Γ^−1/2 v Γ^−1/2||2<1` при true PD Γ.
Малый global minimum eigenvalue может сделать FD непригодным: не добавлять
Γ floor. Relative congruence coordinates v=Γ^1/2 W Γ^1/2 контролируют
cone; при invertible transform spectrum exact J подобен исходному,
но norm/conditioning меняются. При rank deficiency restricted-active
diagnostic должен так и называться, не full spectral radius.

Псевдокод design (22 строки; owning tangent helpers пока отсутствуют):

```julia
function Jv(base, v)
    dR = PI(v) - (im/2)*v
    dG = cavity4_tangent(base, dR)       # actual fixed-V0 recursion
    dA = im*(dG - blockadjoint(dG))
    dmu = -density_trace(base.f .* dA)/base.chi
    dF = base.f .* dA + base.Q*dmu
    out = sum(K(dA) for K in static4)
    out += KLO(Splus(NB*dA+dF) + Sminus((NB+1)*dA-dF))
    check_number_tangent(base, dA, dmu)
    return out
end
base = Pi(Gamma_warm)                   # counted fresh raw physical map
require_interior_domain(base)
v = Gamma_warm / normW(Gamma_warm)
D1 = central_probe(Pi, Gamma_warm, v, 1e-3*normW(Gamma_warm)) # 2 maps
D2 = central_probe(Pi, Gamma_warm, v, 5e-4*normW(Gamma_warm)) # 2 maps
d1 = Jv(base,v); save_derivative_witness(d1, D1, D2)
qualified_signal_or_stop()              # no guessed noise/PSD pass
Q, Hbar = real_arnoldi(v -> Jv(base,v), v; m=4, first_action=d1)
save_ritz_and_residuals(Q, Hbar)
save_restricted_null_witness(Ibar(Hbar)-Hbar, Q)
stop_without_updating_Gamma()           # diagnostic, not new solver
```

Normalization v меняет units шага: amplitudes выше относятся к Γ±εΓ;
в pseudocode шаг поэтому умножен на ||Γ||W. Не подставлять dimensionless
ε blindly в normalized-field formula.
`real_arnoldi` здесь включает two-pass reorthogonalization и reuse d1:
validation action не превращается в пятый analytic Jv сверх cap4.

## 7. Underresolution и physical interpretation требуют другого evidence

Для candidate mode сохранять δμ/χ, ||δA||, occupied δF weighted response,
edge/trusted fractions, per-channel δΓs и energy/k locators. Global
Γ/ΔE quantiles не заменяют occupied linewidth/quadrature checks.
Fine oscillatory/endpoint-localized mode — повод проверить discretization,
не автоматический artifact verdict. Shift-adjoint error и KMS defect
различны; λ current μ fit сам не является kinetic convergence metric.

Для attribution понадобился бы **один matched frozen response**, не fine
SCBA: на той же continuous Γ representation проверить occupied density
independent contour/real quadrature при fixed μ. Если проверяется сама
Γ-map/Jv, matched commensurate refinement m939→1878 сохраняет Ω, a,b
и имеет NE61407. New half-step endpoint ghost taper меняет support;
его вклад явно учесть или перенести старую taper как continuous input.
Full nonlinear map на fine grid всё равно отдельная numerical identity.
Только grid-stable mode/response plus physical conservation/symmetry
evidence позволяет говорить о physical near-nullness. Контурный FD64
oracle на другом Γ не сертифицирует quadrature нового warm Γ.

Это **не** поручение запустить refinement/contour/Jv сейчас. Если имеющиеся
pilot artifacts не показывают plateau или вывод уже решён дешевле,
дополнительные derivative maps не нужны.

## 8. Один следующий warm seed, бюджет и stop

Source candidate: `raw/S01-FD64-continuation128/state.h5`, expected SHA из
его summary `a690ba56918d4c9803d70d759f2c228f4e9133fdc7f4bcede132f19229c0042f`.
Это **unresolved** Core f50/nodal/fractional state, не принятое решение.
Source rK/rΣ/λ, old grid/operator/input identity сохранить; actual H5 hash
перед future use сверить отдельно. Не заменять source на FD64 молча.

Из source взять aggregate **local** Γ=i(sum_sΣsR−adjoint), five channels
в объявленном порядке, embedding исключить. Проверить finite/Hermitian/
PSD witnesses исходного и перенесённого Γ; unresolved snapshot сам по
себе не доказывает cone admission. Перенести на m939 через
convex piecewise-linear interpolation Hermitian Γ, zero ghost extension
явно объявить: новый Emax чуть выше old physical endpoint, но находится
внутри old ghost strip. Не extrapolate constant edge/не занулять его
молча, не interpolate старую complex ΣR/GR independently. New PI Γ
определяет новый seed R; new current μ строится заново. Не rescale Γ
по плотности и не переносить full kinetic history/old μ.

Это legitimate warm **initial guess** новой declared numerical variant,
не exact native continuation и не свидетельство convergence source.
Сравнивать с constant-seed control можно по same raw Π/norm definition,
feasibility, maps/wall до same native gates, с явным различием initial
condition. Два похожих unresolved исхода не докажут root/uniqueness.

Первый рекомендуемый control использует прежний finite pilot16 budget,
same cap/time/AA/grid/operators и mandatory native handoff; меняется
только initialization identity. Cold/warm не запускать параллельно.
Не увеличивать total cap ради derivative probes. Homotopy пока не нужна.
Если после данных есть thermal plateau, сначала free secant diagnostics;
следующий **отдельно выбранный** derivative witness может иметь ≤5 Π
(base+4central) и ≤4 analytic Jv, wall≤300s всего процесса с cache/JIT,
single process/Julia4/BLAS1/прежний MemoryMax12.375GiB+Swap0. Внешний
timeout/предварительная оценка memory обязательны; это design allocation,
не уже разрешённый запуск и не гарантия уложиться.

Один full ComplexF64 Γ на m939:135097600 bytes; basis5 vectors около
0.63GiB сверх primal/tangent/cache. Streaming5×5 cavity избегает больших
full45 tensors, но **весь процесс <3GiB не доказан**. Native guard не
понижать/не обходить. Exact Jv требует PI и scatter contractions и может
стоить порядка одной Π; первое измерение определит реальную стоимость.

Stop: cone/capacity/χ/inverse failure, unqualified FD noise, rank loss,
call/wall/memory/output cap или отсутствие различающего вопроса. Сохранить
raw reasons/scales, не создавать eigenvalues из NaN/zero denominator.
При Reλ-like >1 на reliable subspace можно обсуждать small Krylov solve
`(I−J)δ=Φ−Γ` с cone/feasibility line search как **будущий** метод; сейчас
не реализовывать Newton framework и не обещать сходимость. Anchored
homotopy Φ_t=(1−t)Γwarm+tΦ имеет derivative tJ и возвращает исходные
уравнения только при t=1; промежуточный root не physical acceptance.
Не объединять homotopy и Krylov в первый warm control.

## 9. Первичная математическая опора и её границы

Inverse/implicit-number/LO tangent identities выше выведены непосредственно
для actual source, не взяты из generic solver recommendation.
[Saad, LAA34 (1980),269–295](https://doi.org/10.1016/0024-3795(80)90169-X)
задаёт Arnoldi projection для nonsymmetric eigenproblems; title/abstract/
DOI дополнительно проверены по [авторской university record](https://experts.umn.edu/en/publications/variations-on-arnoldis-method-for-computing-eigenelements-of-larg/).
Короткий Krylov projection не даёт полного spectrum certificate.

[Saad–Schultz, SIAM SISC7 (1986),856–869](https://epubs.siam.org/doi/abs/10.1137/0907058)
обосновывает residual-minimizing Krylov solve и связь с Arnoldi; это опора
условного linear correction, не convergence theorem QCL.
[Knoll–Keyes, JCP193 (2004),357–397](https://doi.org/10.1016/j.jcp.2003.08.010),
[primary paper PDF](https://www.inf.ufes.br/~luciac/comcie/JFNK-Knoll.pdf)
рассматривает matrix-free Jacobian actions и существенную роль
preconditioning. Его общий JFNK apparatus здесь не предлагается строить.

[Toth–Kelley, SIAM JNA53 (2015),805–819](https://epubs.siam.org/doi/abs/10.1137/130919398)
доказывает local AA convergence при contraction/bounded coefficients;
этого нельзя предполагать по имени AA-II или одному small manufactured
test. Источники проверены browsing; никаких новых numerical результатов
или полного доказательства sufficient discretization они не предоставляют.
