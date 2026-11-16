# Same Pi, bounded larger affine coefficients

Load ONLY final Gamma128; start AA/streak empty and recompute current mu.
First fresh Pi repeats parent rawPhi128 within machine replay allowance.
Same Hermitian-PI/cavity/static4/nativeLO/FD-restricted functional; no new physics.

```julia
Gamma = load_gamma_restart(ctx, cfg, output) # exact final128 Gamma only
history = empty_history()
evaluation = Pi(Gamma)                     # boundary/precheck call1 of32
prove_parent_boundary_replay(evaluation)
while calls < 32 && !original_trigger
    proposal, info = aa_trial(history; affine_l1_max=4096)
    save_attempt_TOML(info.suffix_decisions) # rank/condition/L1/reason, scalars
    Gamma, evaluation = original_fresh_domain_guard_and_one_fallback(proposal)
    remember_original_history_and_cached_working_norms(evaluation)
end
if original_trigger
    release_history_and_old_Green()
    original_native_unprojected_exit_and_fresh_audit()
end
```

Only cap10→4096 changes numerical admission. Hist3/cond1e8/α.2/β1/μprecision/
PSD/relative trigger/native gates unchanged. A larger affine span may still
fail feasibility or amplify fullPi noise; FFT roundoff alone is not its bound.
Telemetry and padded-HDF/storage/resource budgets are separate from Pi physics.
