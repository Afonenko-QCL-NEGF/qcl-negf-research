# Runtime-only V2: BLAS setter/getter

ROOT attempt01 timed out120s, exit124: N128 result preserved, N256 unfinished.
No completion is claimed. Its BLAS getter was **unrecorded/unknown**. Root later
observed getter128 despite environment1 and getter1 after explicit setter;
trace contained zgemm_thread_nn/GOMP_parallel. Old loader.source() did not inherit
old run_actual's BLAS.set_num_threads(1), causing this execution-control gap.

Adapter3fe642/inputc9b451/SHA manifest +accepted environment Manifest archived in
adapter-v1-timeout. Actual attempt01 and logs untouched. V2 calls explicit setter
and asserts getter1 before source/model work, after dynamically imported sources,
and before each moment linear solve. Startup log/report record actual BLAS and
Julia thread getters. Input remains false; only code SHA changes. No math,
source, tolerance, step, μ or model change. MM9182 and prior tiny9 unchanged.

ROOT may authorize one new read-only attempt02, bothN, same2GiB/120s/J1/file32M/
data32M/Swap0, using new output. No new innerΦ or μ fit; runtime limit is not
extended and partial V1 results are retained. World-age warnings are recorded,
not addressed in this runtime-only patch.
