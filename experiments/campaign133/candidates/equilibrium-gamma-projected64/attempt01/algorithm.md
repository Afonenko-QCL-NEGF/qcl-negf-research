# Projected-AA Γ64: прежний physical map, новое предложение

```julia
native = build_context(input)                    # original raw h/cavity
Γ = load_gamma_restart(native, input)             # pinned Γlast16, no old μ/Green
states, residuals = [], []                        # empty actual-state history
for call in 1:64                                  # includes entered failed maps
    cone(Γ)
    GR, A, q, ledger = stable_spectral_state(native, Γ) # unchanged eeca operator
    μ = fit_mu(q, target, 70K)                    # same A, fresh μ
    GFD = fd_green(GR, A, μ, 70K)
    Φraw = static4(A) + native_LO(GFD)             # never projected
    Fraw, Fworking = Φraw - Γ, Herm(Φraw) - Γ
    record(Fraw, Fworking, ledger)
    call == 1 && boundary_replay_full_Phi_and_9_scalars()
    strict_raw_trigger() && break                 # raw relative/weighted1e−9 ×3
    push_actual_history!(Γ, Fworking)
    X, info = AAII(states, residuals; m=3, L1=4096, cond=1e8)
    if info.kind == "AAII"
        Y = PSD_positive_part(Herm(X))             # eigen floor exactly zero
        preserve_projection_witnesses(X, Y)        # full first4, selected later
        X = Y                                    # history will store actual accepted Y
    end
    Γ = feasibility_evaluate_or_single_convex_fallback(X) # no Φ clipping
end
release_history_and_old_green()
trigger && original_raw_h_native_stage(native, Γ) # unchanged unconstrained acceptance
```

Per-worker reusable LAPACK buffers дают `P+(X)` без per-node workspace allocations. Exact feasible Herm blocks копируются, а общий unchanged field возвращается без численной реконструкции. Это mathematical nearest-PSD projection с измеряемой backward error, не строгий floating positivity certificate. Projection step stationary не означает root: принимается только unprojected raw residual и native audit.

IO и ограничение времени не входят в физический map: первые4 full raw/projectedAA pairs заданы до опыта; затем small worst blocks+scalar ledger. Бюджет64TOTAL/900pilot/2100overall, resource/file caps — в protocol/input; никаких автоматических extensions.
