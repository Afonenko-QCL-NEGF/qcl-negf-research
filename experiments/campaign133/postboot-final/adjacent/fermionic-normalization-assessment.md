# Fermionic odds normalizer: корректная ограниченная domain, не universal repair

Research-only; Core e171 `production_state.jl:54–74`,
`reference/greens.jl:152` прочитаны. Core/Julia/Pkg/H5 не менялись/не запускались.
**Предложение математически верно на допустимой domain и сохраняет exact raw
SCBA roots. Но замена PairedConvex не универсальна: pinned0/1 могут сделать
target недостижимым.** Это отдельная численная working-map identity, не α tweak
и не добавление e–e. Сходимость/ускорение не доказаны; warm128 остаётся unresolved.

## Support, positivity и единственная scalar root

Для каждого energy/k блока пусть n,p⪰0 и B=n+p. Если v∈kerB, то
`v†nv+v†pv=0` ⇒ nv=pv=0. Поэтому обе матрицы живут на P=ranB. На нём
`C=B^{+/2} n B^{+/2}`, `0⪯C⪯P`. При C=Udiag(c_i)U†,

`φδ(c)=exp(δ)c/[1−c+exp(δ)c]`, `φδ(0)=0`, `φδ(1)=1`,
`nδ=B^{1/2}φδ(C)B^{1/2}`, `pδ=B^{1/2}[P−φδ(C)]B^{1/2}`.

Это доказывает PSD и nδ+pδ=B без commutativity B,C. B=0 вносит0; tiny
положительные eigenvalues B не являются нулевой support. Для точных
weights `w=gs*wE*wk/(2π)≥0` положим `a_i=w u_i†B u_i≥0`. Тогда

`N(δ)=Σ_i a_i φδ(c_i)`, `N′(δ)=Σ_i a_i φδ(c_i)[1−φδ(c_i)]≥0`.

Здесь сумма по всем blocks/modes; **не** Tr C и не unweighted occupation.
Strict monotonicity имеет место, если есть хотя бы один positive-weight
mode с0<c<1. При `Nmin=Σ_(c=1)a`, `Nmax=Σ_(c>0)a` finiteδ имеет range
`(Nmin,Nmax)`; там root единственен. Endpoints требуют δ→±∞, что не входит
в finite proposal. Если active weights нет, N постоянно: target=Nraw →
выбратьδ0/identity; иной target → unavailable. Zero-weight modes не задаются
number constraint, хотя формула выбирает для них тот же commonδ.

Контрпример universal replacement: B=I₂, n=diag(1,0), p=diag(0,1), w=1.
Всеδ даютN=1; target=.5 допустим для Pauli capacity2 и старого PairedConvex,
но недостижим odds shift. Нельзя открыть pinned modes clipping/floor.

## FD, fixed points и ограниченный entropy смысл

При `n=fβ(E−μ)B`, `p=(1−fβ)B` имеем
`φδ(fβ(E−μ))=fβ(E−μ−δ/β)`: global FD shape сохраняется, μ→μ+kBTδ.
Это FDT относительно **B**; физический FDT относительно retarded A требует
также B=A, чего на несошедшемся native state не предполагаем. Под bias C
может быть произвольной matrix occupation: uniform odds shift не навязывает
FD, температуру или e–e thermalization. Sector-FD µ_r получают одинаковый
offset, поэтому их различия этот normalizer также не обязан устранять.
Это свойство **одного normalizer**, не всей native trajectory: linear mixing
ΣL/G families при разных µ обычно также нарушает общий FD factor.

При Nraw=Ntarget следует возвращать исходные raw n,p **без reconstruction**.
Тогда любой exact SCBA root с rawλ1 остаётся root. Обратного утверждения
без raw gates нет: `normalizer(F(x))=x` может выполняться при F(x)≠x.
Даже scalar raw map nraw=a n сa≠1 можно normalization вернуть к target.
Не утверждается сохранение всех fixed points старой normalized iteration.

На fixed eigenbasis исходного C формула — minimizer **классического** binary
relative entropy `Σa_i[x_i log(x_i/c_i)+(1−x_i)log((1−x_i)/(1−c_i))]` при
Σa_i x_i=Ntarget: stationarity даёт logit(x_i)=logit(c_i)+δ; endpoints pinned.
Это не unrestricted quantum-KL projection. При некоммутирующем B derivative
`Tr B h(C)` содержит Fréchet derivative h, а не просто B h′(C); широкой
quantum variational интерпретации здесь не доказано.

В regular interior около root correction имеет ранг1:
`δnnew=δnraw−D·ℓ(δnraw)/χ`, `δpnew=δpraw+D·ℓ(δnraw)/χ`,
`D=B^{1/2}C(P−C)B^{1/2}⪰0`, `χ=ℓ(D)>0`, ℓ — native number functional.
Это гладкая global charge correction, но не общий SCBA preconditioner theorem;
χ≈0 означает плохую conditioning, а charge-neutral modes она не подавляет.

## Float64 и prospective algorithm

Raw pair сохраняется для audit. Любая Hermitian representation correction
должна быть явной и measured. Машинный PSD allowance не доказывает exact
rank/nullspace: support берётся из доказанного structural nullspace либо
resolved positive spectrum; иначе unavailable/заранее объявленный **global**
fallback. Не обнулять tiny b_i и не заменять c вне[0,1] на endpoints.
После whitening errors усиливаются B^{-1/2}; нужны eig/orth/reconstruction
witnesses, проверка kerB n/p, nδ+pδ−B и фактического number error. Local
power-of-two rescaling B,n,p допустима как arithmetic, не eigenvalue floor.

Даже computed c==1 не доказывает exact pin: B=n+p может округлиться к n при
tiny положительном p. C-only logit тогда теряет hole information. Требуется
resolved complementary occupation либо отказ; не скрытое pinning. Считать
logistic и hole совместно, pδ не получать subtractive B−nδ. Альтернативная
exact identity `nδ=t n(p+t n)^+B`, `pδ=p(p+t n)^+B`, t=eδ, полезна как oracle,
но сама не гарантирует stable/PSD Float products и дорога при каждом δ trial.

```julia
raw_number == target && return original_pair, 0.0  # no repair; raw guards remain
for block in pair
    P, b, V = qualified_support(B = n+p)           # unresolved -> unavailable
    c, U = eig_whitened_n(n, V, b)                # no clipping; exact pins only
    W = V * Diagonal(sqrt.(b)) * U
    a = native_weight .* vec(sum(abs2, W; dims=1))
    retain(W, c, a, complementary_witnesses)
end
require_target_in_finite_odds_range()
delta = monotone_root(N_of_delta, target)          # conditioning/precision guard
for block in retained
    f, h = stable_complementary_odds(c, delta)
    nwork = W * Diagonal(f) * W'; pwork = W * Diagonal(h) * W'
end
measure_raw_to_working_and_number_defects()        # no scientific pass here
```

Это high-level design, helpers ещё не реализованы. Стандартный путь требует
двух dense eigensystems **с vectors** на block плюс products; globalδ root
после decompositions scalar-only, но W/c/a надо хранить либо повторять eig.
Хранение W — примерно один дополнительный ComplexF64 field, c/a — ещё два
real mode arrays; native memory guard потребуется пересчитать. Выигрыш wall
не предполагается, особенно на фоне issue32 residual eig cost.

`λ=Nraw/Ntarget` сохраняется **до** correction; raw rK/rΣ/rλ definitions и
fresh unprojected candidates/strict3pass/PSD/causality gates не меняются
(`production_iteration.jl:1–66`). На иной trajectory числа residuals будут
другими, но не их equations. Первый FD-preserving δ или точный Nwork не есть
сходимость: normalizer не принимает ни FD-seed λ, ни reduced step вместо
fresh full native audit. Это проверяемая numerical hypothesis, не результат.
