# Hermitian-PI Gamma-last continuation128

Implementation-only freeze; автор не запускал Julia/Pkg/compute/commit.
Separate copy from completed Herm16, parent sources/raw unchanged. Parent
Gamma operator source SHA3b296973... (Coree171/tree742) inherited exactly for
Π/local_retarded/AA; only Gamma initializer, finite envelope and disk estimate
change. This is not full native exact restart and not a new convergence result.

Parent `../equilibrium-gamma-warm-hermitian/attempt01/`:16Π, raw rGamma
0.0006837267126878922, weighted0.0006250884102415654, no trigger/native.
Import ONLY Gamma-last.h5 SHA8300b63e... . No interpolation, Hermitian repair,
scaling, old R/L/G/embedding/current mu/mixer/streak imported. Current mu is
refit inside every fresh Π; new AA history starts empty. Scientific question:
asymptotic raw-map trajectory from a near-state under the SAME Hermitian-PI
functional, not another cold16 sweep.

Source guard pins parent input6f9813c6..., identity09b29ae8..., operator3b296973...,
checkpoint8300b63e..., model74f9578f..., map01699f539f7..., terminaled638d17... .
Parent identity binds operator/Core; all inherited physics/algorithm options
must match. Allowed changes are only initialization metadata, identity/runtime,
128 call limit,1200/2400s and8GiB data. Checkpoint schema/identity/call16 and
logical shape checked. This checkpoint contains no axes: E/k/weights/scales/
basis/kernels/profiles/zeroU are bound through parent model SHA and current
native problem/model exact array comparison, not invented checkpoint fields.
Gamma cone is checked without modification; no source state arrays are loaded.

Fresh first map, INCLUDING current mu, native five-channel/KMS/cavity precheck,
uses call1 of128. Then boundary-replay.toml compares raw Phi against parent
terminal rawPhi via1024E hyperslabs and rGamma/weighted/absolute/scale plus
mu/number/derivative against parent map016. Replay allowance256eps relative,
exact zero/zero, no additive floor; this is machine-replay guard, not solver
acceptance tolerance. Failure stops before continued AA; raw first-map evidence
is already preserved. No extra Π is spent to check replay.

Caps fixed BEFORE run:128TOTALΠ including boundary/rejected/fallback maps,
precheck≤6/180s,pilot1200s,native≤32=(1+30+1)/900s only by original trigger,
overall2400s. Existing trigger raw+weighted1e−9 for3 valid passes and all native
strict gates unchanged. Hard12.125GiB/high11.5GiB/heap8G/J4/BLAS1/Swap0,
Root reserve1400MiB; file5GiB, new finite data8GiB. Root global32GiB, original
raw preserved. Intermediate Gamma call8 inherited, best/last finally; no
implicit128→larger extension and no new orchestrator.

HDF storage admission now uses EXACT existing writer policy: row chunk
min(NE,max(1,floor(8MiB/(8*other_dims)))), TWO Float64 datasets. Expected complex
bytes=16*ceil(NE/rows)*rows*other_dims. For30704×11×5×5 rows3813 pad to34317:
150994800B per field,28 native fields4227854400B before condition/scales and
existing32/48MiB metadata margins. All Gamma checkpoints/actual map/reduced
terminal/native save+admission use this; writer/chunk policy unchanged.
New tiny real-HDF fixture forces16MiB logical→24MiB chunks and checks declared
bound>=actual filesize. No deletion/rewriting of prior data; admission refusal
is honest stop/not_measured native exit.

Inherited physical104 tests reuse parent evidence (same microscopic/PI/math
operators). Only new restart/source-contract/budget128/storage tests run here;
parent test definitions are copied for provenance, not rerun. Root tiny command:

```sh
julia --startup-file=no --check-bounds=yes --threads=1 --heap-size-hint=2G --project=.build/campaign133/environment/accepted -e 'include(".build/campaign133/equilibrium-gamma-hermitian-cont128/run_gamma.jl"); restart_small_checks(TOML.parsefile(".build/campaign133/equilibrium-gamma-hermitian-cont128/pilot128.toml"))'
```

After independent source review/admission Root makes authorized input copy and
runs same run_gamma.jl INPUT NEW_OUTPUT under stated OS caps. No author runtime.
SHA256SUMS/parent-to-cont128.patch/handoff.md preserve fixed reviewable inputs.
