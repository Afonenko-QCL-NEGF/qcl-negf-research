# One-time FD initialization, then unrestricted native SCBA

```julia
parent = verify_and_load_private_full28(same_voltage_coarse32)
check_native_psd_budget(parent.A)             # every block; no repair
for family in all8families(parent)
    check_native_psd_budget(i*(family.R-family.R'))
end
q = number_measure(parent.A, energy_weights, k_weights)
mu = fit_number(q, target, temperature)        # once, before map1
f, hole = stable_fermi_pair(energy, mu, temperature)
for family in all8families(parent)
    gamma = i*(family.R-family.R')            # block workspace, R unchanged
    family.L .= i*f*gamma
    family.G .= -i*hole*gamma
end
parent.GL .= i*f*parent.A
parent.GG .= -i*hole*parent.A
verify_retarded_digests_unchanged()
save_seed_witnesses_and_changes()             # no second full-state seed file
result = native_full_scba(parent; max_maps=128,
    consume_initial=true, resume_iterations=false)
# Native fresh Keldysh/LO/embedding maps, no FD projection thereafter.
# First GR/A replay gated; expected seed L/G comparison is diagnostic only.
# Initial lambda≈1/J≈0 are constructed. Native raw strict3pass decides.
save_final_full_native_state(result)
```

This is one initialization knob on the same matched0/1mV model. A finite-bias common-FD seed is not an equilibrium statement; native shifted embedding subsequently breaks that common factor. IO/resource/UTC limits are separate admission conditions, not changes to the physical fixed-point map.

Resource scheduling: residual phase uses up to4 workers, all other phases keep requested width up to26; unchanged native map/guard/workspace capacity. This is declared separately from the one-time FD scientific initialization knob.
