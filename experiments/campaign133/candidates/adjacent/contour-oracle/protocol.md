# Протокол ограниченного frozen-Γ Matsubara oracle

Статус: **write-only, не запускался**. Ведущий должен сначала завершить/приостановить эксклюзивный Γ-pilot, выполнить независимый static review, preflight и явно дать resource-green. Этот каталог не меняет Core/Runner, не устанавливает packages и не включает внешние реализации.

Цель — off-axis equilibrium density точно той же frozen local Γ и PI full-hat модели, source μ=0.008495393439818052 eV, T=70 K, zero bias, finite_chain4 на каждой стороне, все11 k. Source H5 SHA e6fffb…; Γ digest edb853… должен совпасть bitwise с prior control. Нормировка плотности/fitting μ запрещены. `input.toml` фиксирует primary dependencies, hashes и prospectively enforceable limits.

Source Hamiltonian hcell берётся **из сохранённого small reference H5** с exact file+array hash: это тот же h, построенный prior control из source-v3 numerical basis operators, а не новый rebuild. Семь basis matrices берутся из source H5; их digests сверяются с prior frozen-control summary. T−=T+†, нулевые bias/Hartree и все grids/scales сверяются до solve. Для уменьшения shared error harness не импортирует QCLNEGF/QCLNEGFRunner и не вызывает их helpers. Окружение pi-B используется только для уже установленных HDF5/stdlib dependencies; никаких Pkg действий. Поэтому binding preflight здесь относится к HDF5 layout и Julia LinearAlgebra, не Core private helper names.

Перспективная команда, **не исполненная**:

    julia --startup-file=no --threads=1 --project=.build/campaign133/environment/pi-B \
      .build/campaign133/adjacent/contour-oracle/matsubara_oracle.jl \
      .build/campaign133/adjacent/contour-oracle/input.toml NEW_OUTPUT_DIRECTORY

Ведущий оборачивает её штатным scope с MemoryMax=3GiB, CPUQuota=200%, OMP/OPENBLAS=1, timeout ≤600s, writable output cap32MiB. Внутренний phase/chunk stop570s дополнительный, не заменяет внешнего ограничения. Должен быть единственный bounded compute process в разрешённом resource gap. Не использовать cluster/CI/global settings.

Последовательность после green:

1. Hash guards (source, prior script/summary/small H5, active project/Manifest), затем чтение small frame/basis/Γ без kernels. Validate all actual input identities, guard exact Γ hash. В source frame source Γ из пяти retarded fields суммируется один раз в declared order; embedding исключён.
2. Scalar Fermi, complex-Hermitian45 Γ=0 vs exact f(H45), source-H45 Γ=0 controls; BigFloat256 exact coefficient vs series+remainder. При control failure остановить и сохранить primary witnesses.
3. Σ(z) chunks32 → full45 LU central5 RHS → explicit Hermitian accumulation. Сохранить N128/256/512/1024 checkpoints, M1 and M3 tails. Никакой SCBA/Poisson/real-axis full Green.
4. Сравнить каждый checkpoint с `fresh_PI/F` исходной positive real trapezoid measure (original [a,b], ghost strips в той мере отсутствуют). Первичны matrix absolute norms/per-k Δρ и signed ΔN/target. Relative matrix norms сохраняют явный denominator; при нуле denominator relative помечается unavailable, не делится на arbitrary floor.
5. Report execution, controls, truncation bounds, observed N refinement и ограниченный physical conclusion отдельно. Нельзя назвать contour result formal exact из согласия N512/N1024. Нельзя объявить разность density доказанной real-grid quadrature error без ghost/atom/evaluation qualification.

Outputs: copied input/script, summary.toml и small `density_witnesses.h5`, без full Γ/G/Σ tensor. При timeout частичные completed N checkpoints и summary сохраняются. Минимальное различающее решение: если controls проходят и M3 sequence стабилизируется значительно точнее Π−nodal density delta, оценить off-axis density discrepancy; иначе ограничиться отрицательным runtime/convergence выводом. Этот protocol не принимает исходный SCBA и не меняет source tolerance.

Дополнение: code также сохраняет conditional floating ghost-strip margins/bounds из уже загруженной Γ, без дополнительных Green solves. Этот witness не принимает positive computed PD как interval certificate; никакого abs/floor отрицательных margins. Он может быть сопоставлен с контурной разностью по theorem `pi-ghost-mass-bound.md` только с явно указанными PSD/input/evaluation assumptions.
