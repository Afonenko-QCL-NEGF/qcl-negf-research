# ROOT runtime handoff, unexecuted V2

Author Julia/Pkg/HDF5/solver executions: **0**. V1 module7e847d…191LOC preserved
under `unexecuted-v1/` before reviewer-driven P1/P2 changes. V2 module1434e86…218LOC
and adapter are independently reviewable; runtime-input.toml stays false.
Only runtime_authorized may change in a ROOT-owned actual copy; every other
input field and script source SHA is frozen. No automatic step extension.

External ROOT enforcement must cover **the entire process including imports,
source hashes/model build, controls, both N stages and writes**: wall600s,
MemoryMax2GiB, High1.75GiB, Swap0, Julia1/BLAS1, CPU100–200%, LimitFSIZE32MiB,
output64MiB. Source GP build_context retains its original large native-cache
estimate/options for model reconstruction; no native solve occurs and this is
not a claim that its full-native workload fits2GiB. Actual cache payload is
small, but imported JIT/source build RSS is unmeasured: resource refusal remains
an unresolved outcome, never authorizes cap increase.

Command shape after root review/admission:
`julia --startup-file=no --project=.build/campaign133/environment/accepted run_matsubara.jl AUTHORIZED_COPY NEW_OUTPUT`.
Use absolute adapter/input paths or invoke from its directory. No packages are
installed and no input/source file is overwritten. Target output must be new.
The script checks active cgroup memory/swap/CPU and OS file-size ceiling.

Source binding reuses V4's frozen-io01 logical digest helper192e6f… and
completed digest packetb121f7…; only small h_H/h_raw/T+/axes/units and native
raw5 kernels from reconstructed source model are used. **No Γ/Φ source arrays
or source PD certificate are loaded**. New Mats fixed point is unrelated to
acceptance of old Γ. GP50bc +3 directly included definitions are SHA-bound.
`verify_problem_model` compares every native source grid/profile/basis/kernel,
qK and raw h against model700ae…. Saved h_H comes from actual spectral42fe…,
exact-Hermitian and V4-digest bound, Tminus exactly adjointTplus. β is inv native
kBT; Ω uses native scaled physics. Fixed μ is exact map64 value
0.08459614262982018 native =0.008459614262982018eV, not rounded0.0084596.
There is no μ root/density normalization; fixedμ density is an output.

Before actual inner iterations, adapter records independent17freq native vs
six-index contractions, onelevel raw Mats→exact rational refinement, actual
module two-pole wiring (reference/lo_delta/greens/density), and conditional
CP L/ω0². CP diagnostic is not used as finite-cutoff uniqueness/pass claim.
Fixed α=1 is ROOT's predeclared undamped iteration; significant candidate
noncausality, time/resource failure stops, without clipping/damping adaptation.
N128 andN256 each begin from Σ=0 and get at most16 updates,32 total. No largerN,
extra iterations, warm transfer or μ fitting is silently added.

P1 preservation: module invokes aligned snapshot callback with **input** Σ,
its actual G/ρ, raw Φ and frequencies before guards/update; nonmutating mixing
keeps latest snapshot aligned. Every completed raw map writes full Σ/G/ρ
`N*-i*.h5` plus a small scalar record. `N*-latest.h5` additionally retains Φ;
only this explicitly named own-output latest file may be atomically replaced.
Thus raw per-iteration Σ/G/ρ remain, with last raw Φ perN, under≤64MiB (rough
payload51.6MiB per-map +4.8MiB latest +headers/sources). Admission accounts for
old latest plus new temp; no deletion to rescue an exhausted data cap.
Global hard kill may leave .tmp; request ledger bounds entered/completed calls
honestly. Soft module return contains latest coherent fields, not mixed-S-only.

P2 ratios use explicit0/0→0 and nonzero/0→Inf, no hidden epsilon denominator.
Absolute norms remain alongside relative. Shared low32 frequency/ρ comparison
is only a cutoff witness if both raw Σ residuals are small; iteration-budget
outcome may confound truncation/nonlinearity. Scientific/native acceptance and
DOPE target pass remain false even if diagnostics improve.
