# Равновесная итерация одного поля Γ

Неизвестная — local scattering broadening Γ(E,k), Hermitian matrix5×5 на
каждом узле. Embedding в неё не входит. Физика прежняя: static4+LO,70K,
V=0,U=0, finite-chain4. Из Γ native PI/Hilbert строит
ΣR=HPIΓ−iΓ/2, cavity и GR при η=0; A=i(GR−GA). По **текущему** A
подбирается common μ для заданной плотности. Явное равновесное ограничение
GL=ifA, GG=−i(1−f)A сокращает unknowns, но само не доказывает thermalization.
Physical map Φ — сумма четырёх K[A] и i(ΣLO>−ΣLO<) от native LO operator.
Вариант Cdirect здесь выключен: `elastic_retarded=:hilbert`.

Ниже псевдокод науки; Π соответствует `evaluate` и его owning helpers,
g,gs,Ntarget,kBT берутся из ctx; счётчик calls увеличивает owning evaluate
до map. three_passes и ограничение history здесь сокращены. One-fallback ветка применима
только к AA trial; failed linear/fallback означает honest stop.

```julia
function Π(Γ)
    cone(Γ)                              # finite/Hermitian/PSD budget
    GR, κ, scale, _ = retarded_state(ctx, Γ) # PI → cavity → Dyson, η=0
    A = spectral_function(GR); cone(A; require_nonzero=false)
    q = N._seed_number_weights_production(GR, g, gs; options=ctx.options)
    μ = fit_mu(g.ε, q, Ntarget, kBT).mu    # current A, capacity/bracket
    GFD, _ = fd_green(GR, A, κ, scale, g.ε, μ, kBT)
    Φraw, _, _ = phi(ctx, GFD); cone(Φraw) # static4 K[A] + native LO Γ
    return (Φraw=Φraw, GFD=GFD)
end
Γ = 2ηseed * I; states, residuals = [], []; ev = Π(Γ) # no persistent η
while true
    Φwork, rraw, rwork, δrep = gamma_residuals(ev.Φraw, Γ, wE, wk)
    trigger = three_passes(rraw.relative < 1e-9, rraw.weighted_relative < 1e-9)
    (trigger || calls == 16) && break     # Γ trigger is not native convergence
    push!(states, Γ); push!(residuals, Φwork-Γ) # retain at most4 pairs
    trial, info = aa_trial(states, residuals, wE, wk, cfg)
    Γ, ev = try (trial, Π(trial))
    catch err
        err isa Infeasible && info.kind == "AAII" && calls < 16 || rethrow()
        empty!(states); empty!(residuals); linear = 0.8Γ+0.2Φwork; (linear, Π(linear))
    end
end
if trigger; empty!(states); empty!(residuals); ev=nothing; GC.gc(true)
    native_stage(ctx, Γ, budget, out, cfg); end # native G; no repeated FD
```

AA-II использует real weighted Frobenius inner product wE*wk, two-pass QR,
rank/condition/L1 guards и cone/domain checks. Он **nonmonotone**: рост raw
residual сам по себе не вызывает отказ. Известные Singular/LAPACK failures,
включая только численные task wrappers, допускают один convex fallbackα=.2;
programming/budget errors не маскируются. Fallback сбрасывает history и
pass-streak; clipping/Γ floor отсутствуют.
Cone допускает объявленную floating backward error, а не rigorous PSD proof.

I/O и ресурсный budget не являются частью physical map:16Π включая seed и
failed evaluated trials, precheck≤180s, pilot≤600s, overall1800s; hard12.375GiB,
high11.5GiB, Swap0. Γ checkpoints и raw witnesses сохраняются отдельно.
Native exit требует unprojected пятиканальных fresh maps, raw rK/rΣ/λ и
physical audits; μ-fit там только diagnostic. Даже native strict fixed-U
pass не принимает Poisson, дискретизацию или эксперимент.
