# Isolated INFO-layout benchmark, WRITE-ONLY

Input false, author Julia/HDF5/runtime executions0. Root runs only after
source/binding review and after sole native process has finished. No Core edits,
Pkg, model/cache construction, Phi, SCBA, mu-fit or normalization. Source is
clean e171; accepted postreboot Project/Manifest and GP serialization helpers
are SHA-bound. New output required, read-only parent files unchanged.

Fixture: completed warm case0 NE24001/Nk11/Nb5, full state e4fc9b91…,
model3de1f7fc…, parent128maps/unresolved. Load four saved Green fields,
condition/scale and all24 selfenergy fields with existing GP.read_complex;
load only energy/h_fixedU0 from model. State/model axes, U0, units, ownership,
schema and complete-source SHA checks precede use. Build total ONCE by ordered
sum [LO,acoustic,impurity,IFR,alloy,embedding,embedding_plus,embedding_minus].
This deliberately counts embedding sum AND both sides: **SURROGATE FIXTURE,
not the physical total**. Candidate is that same total object, NOT fresh Phi.
Outputs are kernel witnesses, never a new equation/convergence result.

Only production_residual_suite is timed: rD/rA/rK/PSD/causality, six PSD plus
causality eigensolves on finite path. The8family selfenergy-residual passes in
full native residual phase are excluded. Consequently timings do not directly
reproduce the previously measured3.33s full phase.

Create ONE26-capacity workspace with chunk64. Save its original info vectors
(A,length1) and allocate padded INFO vectors (B,length8, tail sentinel).
Swap only entries of ws.lapack_info, not any H/work/rwork/values/scalar array.
Record pointers/cache-line grouping for both layouts, hardware L1 data line
size and sizeof(BlasInt); length8 is64bytes only for BlasInt8. Check all other
buffer addresses unchanged. Core suite resets its own partial arrays each call.

Excluded warmups: A26 then A4. Timed: W4 A/B/A twice, then W26 A/B/A twice;
12samples,14calls total. Save wall/processCPU/allocated/GC per call with
existing Core counter API (no new ccall). Exact isequal on ALL five metrics
and all12 worst-witness fields, including coordinates and first-tie outcome,
against first warmup. Tail sentinels and layout/buffer addresses must remain
unchanged. Final logical SHA of all loaded inputs/total must equal initial.
No result discarded as an outlier, no automatic repetitions/parameters.

Envelope: proposed6GiB/swap0/J26/BLAS1/OMP1/CPU2600%, file/output32MiB.
Retained data28fields=2.957GB plus total3fields0.317GB; per-field GP read
transients~0.106GB, explicit collection after each field, workspace small.
RSS is unmeasured; cap failure stops without enlargement. Source hashing,
load/JIT/fingerprints are outside sample timing but inside Root external300s.
Compute soft budget180s after setup includes warmups and metadata checks;
before a new kernel admission reserve max(5s,1.5*largest prior width wall).
A single admitted kernel cannot be interrupted by this software check;
Root enforces external300s for the whole process. Soft stop preserves samples
and reports stopped, never silently reduces planned repetition count.

After static approval Root makes an authorized copy changing ONLY
runtime_authorized; establishes OS caps and invokes:

```text
env JULIA_NUM_THREADS=26 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 julia --startup-file=no --project=.build/campaign133/environment/postreboot .build/campaign133/core32-info-padding-bench/bench.jl AUTHORIZED_INPUT NEW_OUTPUT
```

Source include crosses Base.invokelatest(run_loaded,...) after GP definitions.
Metadata only summary/input/script; no arrays/HDF outputs. A speed difference
alone does not diagnose false sharing: inspect actual line sharing and A/B/A
repeatability; unlike two independently constructed workspaces, all other
buffers retain their addresses. Original scientific criteria and active Core
stay immutable whatever the result.
