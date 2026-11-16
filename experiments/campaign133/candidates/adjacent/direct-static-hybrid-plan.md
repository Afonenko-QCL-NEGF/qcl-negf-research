# Резервный hybrid equilibrium map: direct-static C + LO HPI

2026-11-16. **Только протокол, не реализация и не запуск. Активный Γ128 остаётся
без изменений.** Исходник для контракта — accepted Core
`e171caa60cce676c6c4238dd671e7c527f2bd977`. Сохранены S01: V=0, общие μ и T=70 K,
LO+четыре static channels, fixed U0, Nb=5/Nk=11, существующее commensurate
energy window/ΔE, basis P3 и finite_chain **4+4**, η=0. Не предполагается
самостоятельная Σ каждого boundary cell: повторяется central-driven local Σ,
как в actual Core. Это новый contingency experiment identity, не exact resume.

## Неизвестные и физическая опора

Вместо total Γ с Hilbert-восстановлением всех механизмов взять
`X=(ΓLO, Sel)`, где ΓLO Hermitian PSD, Sel — общая complex retarded Σ четырёх
static channels. Map:

\[
\Sigma^r_{loc}=S_{el}+H_{PI}[\Gamma_{LO}]-i\Gamma_{LO}/2,\quad
G^r=\text{finite4}(H,\Sigma^r_{loc}),\quad A=i(G^r-G^a),
\]
\[
\Phi_S=\sum_{s\in\{acoustic,impurity,IFR,alloy\}}K_s[G^r],\quad
\Phi_{\Gamma LO}=i(\Sigma^>_{LO}[ifA,-i(1-f)A]-\Sigma^<_{LO}[ifA,-i(1-f)A]).
\]

Каждый K_s тот же exact static contraction с qK/wk, без дополнительного
energy weight или spin. При scalar f(E) и energy-local K_s:
`Σs<=if K_s[A]`, `Σs>=−i(1−f)K_s[A]` и
`i(ΦS−ΦS†)=ΣK_s[A]`. Поэтому fixed point совпадает с **native C + KMS**, если
finite4/Keldysh identity и discrete LO detailed balance действительно выполнены.
Это не совпадение с baseline HPI-all-static либо нынешним Herm-total-Γ map:
`K[Gr]` сохраняет static off-window dispersive response, которого finite-window
`HPI[K[A]]` не восстанавливает. Даже successful C не принимает старую identity.

Код: [reference kernels:646–654](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/kernels.jl#L646),
[production:1240–1246](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/production_operators.jl#L1240).
Предусловия: exact kernels имеют covariance/Kraus форму `Σρ VXV†`, ρ≥0,
положительные wk и adjoint preservation. Они разобраны в
[static-retarded-assessment](../static-retarded-assessment.md), §§1–3;
первичная SCBA опора — [Kalvová–Velický, PRB65,155329, Eqs.63,66,68](https://arxiv.org/pdf/cond-mat/0202113).
Независимый code precedent: [libNEGF elphdd.F90:227–250](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/elphdd.F90#L227)
считает energy-local static Σr непосредственно из Gr; это diagonal dephasing,
не доказательство эквивалентности наших четырёх kernels.
CP сохраняет cone одного map, не доказывает сходимость. Continuum KK equivalence
требует analytic Gr в upper half-plane, полного spectral measure, 1/z asymptotic
и energy-independent bounded K. LO HPI имеет свой finite-support/taper contract.
Pointwise cone не гарантирует глобальную analyticity произвольного AA state.
[Alt–Erdős–Krüger](https://arxiv.org/abs/1804.07752v5) не даёт автоматически
uniqueness/contraction для LO + нашего central-driven finite4 closure.

## Spectral source, μ и смешивание

Local source `Γloc=ΓLO+i(Sel−Sel†)` должен быть PSD. Переносить его через
**каждую** повторённую cavity cell; в центре суммировать local + left/right
injection. A_Gram строится из этих sources и центрального resolvent; только
lead Γ недостаточно. Положительный Gram-LDOS имеет независимый implementation
прецедент [Kwant common.py:552–561](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/solvers/common.py#L552).
На каждом state также сохранять A_sub=i(Gr−Ga), optical-theorem defect и solve
backward error. Нельзя скрыть несовпадение выбором положительного A.

Для слабого уширения Sel удобно хранить **в эквивалентных координатах**
`Sel=Lel−iΓel/2`, Lel/Γel Hermitian: это ровно столько real degrees of freedom,
сколько одна general complex matrix. ΦL=ΣK[Herm Gr], ΦΓel=ΣK[A] эквивалентны
ΦS=ΣK[Gr] и не теряют малую Γel при вычитании больших complex entries.
Проверить эту эквивалентность с принятым native path; дополнительных hidden
Γ storage не обещать. Cone: ΓLO≽0, Γel≽0, causal total/finite4. Никакого
clipping eigenvalues, η-floor или замены Hilbert для LO.

μ подбирать по **исходной** number functional `N[ifA]=Ntarget`; f общий для
всех k/cells/baths. Не умножать A/G< после этого на число и не объявлять такую
нормировку native convergence. Raw native λ/candidateλ должны стремиться к 1.
AA-II depth=3 использует **real** coefficients и real inner product:
`Re⟨ΔSel,ΔSel⟩ + ⟨ΔΓLO,ΔΓLO⟩`; для split Sel это
`⟨ΔL,ΔL⟩+¼⟨ΔΓel,ΔΓel⟩`. В Hermitian packing учитывать √2 off-diagonals.
Масштабы block norms фиксировать до probe в native energy units; per-block
raw residual публиковать отдельно. Negative AA coefficients разрешены только
если **совместный** trial проходит cone/finite checks. Иначе одна общая convex
backtrack между Picard output и AA trial, затем Picard fallback; если сам
Picard infeasible — остановка с диагнозом, не causal repair. Непринятый trial
не добавляется в историю; малый accepted step не является stopping gate.

## Один bounded probe — максимум 16 map evaluations, не запускать сейчас

Только если Γ128 не достигает native acceptance и ведущий выбирает этот опыт:
один seed, без параметрического sweep. №1 — новая C-seed conversion:
из сохранённого G/A отдельно вычислить native FDLO ΓLO и ΣK_s[Gr]; **total Γ
не переименовывать в ΓLO**. №2–12 — hybrid raw maps/AA3 с checkpoints 4/8/12.
№13 — fresh unprojected native-C audit; №14–16 — три free native-C updates
без FD projection. Все retarded/kinetic candidates, embeddings и per-channel
rΣ/rK, λ, causality/KMS/spectral gates — исходные; time cap 20 min либо
оставшийся явно выделенный ресурс, что меньше. Exhausted budget = unresolved.

На exit восстановить каждую из четырёх families отдельным K_s[Gr], LO —
полным native operator, embedding — свежим finite4. До materialization сохранить
`||ΣK_s[Gr]−Sel||`; новая per-family сборка не стирает aggregate residual.
Restart contract должен явно содержать `elastic_retarded=:direct_green`, LO HPI,
finite4, сетку, projection и новую solver identity, **без старой AA history**.
Core [scba.jl:159–166](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/scba.jl#L159)
отвергает warm state с другим retarded selector;
[reference audit:3–8](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/final_audit.jl#L3)
тоже требует `:direct_green`, но этого **недостаточно**: его reference LO
`retarded_self_energy` использует nodal direct H, не HPI. Accepting fresh
candidate обязан сохранить selected HPI LO + C; nodal reference difference
отчитывать как cross-discretization check, не same-map residual.

Memory estimate, **не RSS claim**: для NE=30704 один Hermitian-real field
S≈64.4 MiB; X=3S≈193 MiB (в 3 раза больше Γ-only). Три пары ΔX/ΔF≈1.13 GiB;
четыре live X-like buffers≈0.75 GiB; вместе≈1.89 GiB, до GR/A/FFT/kernels/native
families. Две копии 5×3-component native families добавляют≈3.77 GiB, поэтому
перед materialization освобождать AA buffers. Реальный peak должен быть
предварительно посчитан из driver allocations и общего процесса в пределах
12.125 GiB; здесь соответствие лимиту **не установлено**. Если не помещается —
не запускать; не сокращать scientific state молча.

Псевдокод (23 строки, новые имена условны; общий бюджет включает exit):

```julia
function hybrid_map(x, ctx)
    (; ΓLO, Sel) = x
    Γel = stored_broadening(Sel)             # split coords, не large subtraction
    Γloc = Γel + ΓLO
    Rloc = Sel + HPI(ΓLO, ctx) - 0.5im * ΓLO
    g = finite4_same_central_repeat(ctx, Rloc, Γloc; eta=0)
    A = source_gram_spectral(g, Γloc, ctx)
    μ = fit_mu_original_number(A, ctx)
    f = fermi(ctx.E, μ, ctx.kBT)
    fd = (lesser=im .* f .* A, greater=-im .* (1 .- f) .* A)
    nextS = pack_static_complex_map(g.retarded, A, ctx) # ΣK_s[Gr], split cone
    nextΓ = native_FDLO_broadening(fd, ctx)
    return (ΓLO=nextΓ, Sel=nextS), (μ=μ, g=g, A=A, fd=fd)
end
x = explicit_new_C_seed(saved_state, ctx)   # provenance conversion, map №1
for probe in 2:12
    y, witness = hybrid_map(x, ctx)
    save_raw_map_residual(x, y, witness)
    trial = real_AAII((ΓLO=x.ΓLO, Sel=x.Sel), y; depth=3)
    x = joint_cone_backtrack_or_picard(y, trial, ctx)
end
native = materialize_five_families_with_C_contract(x, ctx)
audit_and_free_native_C(native, ctx; audit_maps=1, free_maps=3)
```
