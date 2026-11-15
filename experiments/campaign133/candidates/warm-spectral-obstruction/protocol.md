# Selected warm spectral obstruction diagnostic

Implementation-only freeze; автор Julia/Pkg/compute не запускал, warm/Core/
Runner не редактировал. Root уже выполнил один native retarded/spectral-source
probe (20.34s process, peak2864156KiB, exit0): `../warm-q-quick/attempt01/`.
Snapshot selected.h5 SHA45fa3f6859802aa25ef19b7d16ccfae3e90c297463e3f374ca57ba85b3ca918d.
Все1429 negative q повторяются same-G Γeff product; physical source probe
negative0. Это измерение Root, не диагноз по oracle, который ещёне выполнялся.

Теперь **никакого full retarded pass**: one-file oracles.jl читает только
завершённый1.55MB snapshot и Tminus из hashed cold model6b41... . No fitμ,
Φ/scattering/SCBA. Root задаёт6GiB/300s/J4/BLAS1/heap4G, file/data32MiB.
`input.toml` runtime_authorized=false; Root создаёт authorized copy.

Для first negative E104 и worst q E4145 — **все11k**; дополнительно последний
E30704,k11 и occupied reference E7894,k1. Последние имеют label selected-k-only,
не доказывают whole q. Для каждого узла четыре counterfactual: raw, HermΛ,
Hermh, both. Только selected nodes, никаких global fixes.

45×45 finite chain с9 cells и общей localΣR. **Native side h отличается:**
8 neighbouring cells используют frozenF64 hlocal=(hraw+hraw†)/2 по
production_operators.jl:309–315, central cell используетhraw. Эта F64
side matrix реконструируется до Big conversion и сохраняется. Γeff sum
использует собственный onsite каждого cell. V1 repeatedrawh invalid для
native attribution; source/SHA сохранены invalid-repeated-h-v1/. Accepted native finite_chain
`production_operators.jl:370,384` реально использует left adjoint(Tplus),
right Tplus. Поэтому uniform upperTplus/lowerTplus† совпадает с этой веткой;
Tminus−Tplus† из actual hashed model записывается отдельно, не усредняется.
Примитивы energy,h,Tplus,localΣR,Γ берутся как **точные сохранённые Double**.
Float direct inverse full45 сопоставляется с native selectedGR. Big192 LU
**пять RHS D†Y=P**, row=Y†: это central ROW, а не mistaken central column.

Big raw D наследует frozen input non-Herm defects. Для HermΛ variant сначала
в working precision Λ=localΣR+iΓ/2, затем HermΛ−iΓ/2; Hermh отдельно/both.
Snapshot не содержит original pre-subtraction FFTΛ, поэтому antiHermΛ здесь
явно reconstruction witness (включает rounding локальногоΣR), не чистый
отдельный FFT-error bound. h antiHerm записывается отдельно.

Сравниваются i(Gcc−Gcc†), stable sum_j G_cj Γ G_cj† и stable effective sum с
Γeff_j=i[(h_j+ΣRlocal)−(h_j+ΣRlocal)†]. Same-G native A_sub/Aeff/Aphysical и
Γeff/physical Γtotal probe сохраняются рядом. Root physical probe использует
local lesser=iΓ/greater0 только как технические spectral sources, не occupation.
Γ eigenvalues/Hermiticity/norm и относительный eigensolve scale сохраняются
без floor/clipping. Big exact traces/residuals — strings192bit, signed q sums используют эти
192bit traces и192bit arithmetic, а не downcast matrices. Float q sums
сохранены отдельно. На strongest-negative k при worst q energy выполняется
384bit repeat всех4 variants с trace difference; остальным nodes не
приписывается per-node precision qualification. Float–Big difference включает
projection/assembly rounding и LU, не объявляется только inverse error. Small matrices
дополнительно Float64 для HDF. Все norm ratios имеют0/0=0, nonzero/0=Inf.

Root tiny/API command (no S01 map):
`julia --project=.build/campaign133/environment/accepted .build/campaign133/warm-spectral-obstruction/oracles.jl --tiny`

После static review и OS gate: тот же script `AUTHORIZED_INPUT NEW_OUTPUT`.
Tiny controls — scalar three-cell analytic Schur, noncausalΣ input vs HermΛ,
complex Hermitian two-band chain, mixed centralraw/sideHerm analytic Schur
и Big/Float consistency. No new dependencies.
Output witnesses.h5, traces.csv, summary.toml и exact input/script copies.
Критерии solver не меняются, acceptance не заявляется. Следующий вывод должен
разделить frozen-input antiHerm defect и inverse arithmetic по actual outputs.
