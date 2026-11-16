# Равновесный Γ cont64: тот же stable-A оператор

Неизвестная — total local Hermitian Γ(E,k), общая равновесная μ при70K. Внутренний reduced оператор использует causal PI/Hilbert closure и native finite_chain4; covariance A вычисляется из физического spectral source. Это экспериментальная double арифметика без гарантии PSD через факторы. AA меняет только следующее предложение Γ; Φ и исходная невязка не проектируются.

```julia
native = build_context(input)                    # raw h/cavity preserved
Γ = load_gamma_restart(native, input)             # only last Γ16, exact source/frame guard
history = []                                     # no old μ/Green/AA history
for call in 1:64                                  # includes entered failed maps
    cone(Γ)                                      # no repair
    spectral = stable_context(native)            # exact plan.hlocal all9
    Λraw = PI_fft(Γ)
    ΛH = Herm(Λraw); record_correction(Λraw, ΛH)
    ΣR = ΛH - im * Γ / 2
    emb = native_cavity(spectral, ΣR, im * Γ, 0)   # spectral probe, not occupations
    GR = native_retarded_green(spectral, ΣR + emb.R)
    A = GR * (Γ - im * emb.lesser) * GR'            # blockwise F64 congruence
    q = number_weights_A(A)                       # same A throughout
    μ = fit_mu(q, target, 70K)                    # fresh each call
    GFD = fd_green(GR, A, μ, 70K)
    Φ = static4(A) + native_LO(GFD)               # raw, no clipping
    record(Φ - Γ, Herm(Φ) - Γ, spectral_ledger)
    call == 1 && boundary_replay_full_Phi_and_9_scalars()
    strict_reduced_trigger() && break
    Γ = AAII(history; m=3, L1=4096, cond=1e8)      # existing feasibility/fallback
end
release_history_and_old_green()
trigger && original_raw_h_native_stage(native, Γ) # materialize + solver + fresh audit
```

Scientific map/AA и IO budgets разделены: запись selected ledger не добавляет map. Domain failure после входа в Π считается вызовом; input-cone reject до него — нет. Выход из reduced trigger не является scientific acceptance: original unprojected native graph/gates остаются обязательны. Runtime64 и finite resource limits заданы протоколом, а не формулой.
