# Same Gamma functional, new empty-history continuation

Unknown Gamma(E,k) and Π are inherited from Herm16: native uniform PI realpart
is explicitly Hermitian-projected within original1e−10 arithmetic guard,
ΣR=Herm(HPI Gamma)−iGamma/2; native cavity/Dyson η0 gives A, current common mu
fits target density, FD pair feeds static4 K[A]+native thermal LO Γ. Γ-only
real weighted AAII m3, PSD/domain safeguards and one convex fallback unchanged.

```julia
Gamma = load_gamma_restart(context, cfg, output) # ONLY last Gamma, exact bytes
states, residuals, streak = [], [], 0
mapped, F, stats, streak, trigger = evaluate(context, Gamma; precheck=true)
# evaluate also proves raw first Phi repeats parent map016 within machine error
while total_Pi_calls < 128 && !trigger
    trial, info = aa_trial(states, residuals, weights, cfg)
    Gamma, evaluation = fresh_Pi_with_one_counted_feasible_fallback(trial)
    update_empty_started_history!(states, residuals, Gamma, evaluation)
end
if trigger
    release_history_and_old_Green(); GC.gc(true)
    native_stage(context, Gamma) # original unprojected all5 strict exit/audit
end
```

No parent mu/R/L/G/embedding/AA history or interpolation. First boundary counts
in128. Raw residual, not mixed step, governs trigger. New budget1200/2400s,
data8GiB,file5GiB; padded-HDF estimate changes admission only, never arrays or
mathematics. Root OS hard12.125GiB/high11.5/Swap0. Native convergence/Poisson/
discretization/experiment remain separate from this bounded trajectory probe.
