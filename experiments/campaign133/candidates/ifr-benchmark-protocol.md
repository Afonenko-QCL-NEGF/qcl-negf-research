# Frozen actual-field IFR operator benchmark — not executed by author

Question: does the exact canonical structured IFR contraction reduce actual
operator wall time at S01 shape while preserving every complex output entry?
No full solver/convergence/LO/physics speed claim. Core sources remain frozen.

Input: S01-warm-FD64/state.h5 SHA e6fffb997f36397f877aee1e2c98b1dbad7d97417d45458a042b5f074ce57a6d,
source v3 green/GR_real,GR_imag, shape30721×11×5×5. Stream H5 hash separately.
Source basis seven operators, k/wk, x/wx, scales, normalized IFR K and qK must
match canonical construction bitwise; no silent basis alignment or rescaling.

The canonical basis/kernel is built with the same source Nz96/Nb5/Pbasis3/
Nk11/Nphi32 but NE17. IFR is energy independent. This selects its cache operator
without constructing a full31k solver/cache; applying it to full savedGR is an
explicit operator-only benchmark, not a bypass of full solver RAM admission.
Source numerical budget is retained for the small builder. Source rawK/qK/wk
are checked, qK and wk applied exactly once identically in both arms.

Dense reference is native generic _production_kernel(source normalizedK),
structured is selected canonical build_production_cache(...).kernels[:IFR].
Both get the same preallocated output shape, input, chunk64, Julia4/BLAS1.
Warm once each; 5 alternating-order measured pairs, GC outside timing before
pair; @timed captures apply wall time, bytes and GC time including scratch/
packing. No third-party benchmark package. Compare all complex entries with
raw Frobenius numerator/scale and max entry numerator/global reference scale;
limit2e-12 is an operator arithmetic check, not changed scientific tolerance.
Record input/output hashes, actualtypes/payloadbytes, median/min/repetitions.
No full field output. Source unchanged hash is checked after measurement.

Prospective external cap:4GiB, CPU400%,300s TERM+10s KILL, output≤32MiB.
Root launches only after current scopes/resource green and static review.
No Pkg/install inside script; root-prepared environment must already contain
QCLNEGF candidate/HDF5/Unitful. Code checks actual import, base e171 and frozen
patch8b47d71a... plus implementation diff00e2e7fa....

```sh
OPENBLAS_NUM_THREADS=1 /usr/bin/time -v timeout --signal=TERM --kill-after=10s 300s julia --startup-file=no --check-bounds=yes --threads=4 --project=.build/campaign133/environment/ifr .build/campaign133/ifr-benchmark.jl .build/campaign133/ifr-benchmark-input.toml .build/campaign133/ifr-benchmark-attempt01
```

Archive stdout/time-v and applied OS cap separately. A positive operator
result still needs actual IFR runtime share before any end-to-end claim.
