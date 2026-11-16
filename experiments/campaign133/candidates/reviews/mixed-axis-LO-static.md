# Mixed-axis LO: independent static sign review

Проверен research-only план SHA256
`736aca366fb53b2a8c0397f487426957cb144a8a96e61c3e1bd09efbb5b374e8` и
`matsubara-all5/matsubara_all5.jl` SHA256
`1434e86fe452f02c5831eaeb0fc83dc1dc12710dcc3d7e81f7f36cd985a31447`.
**Знаки C, shifted weights, β и scalar Fock oracle правильны.** Будущий
single-ladder question математически оправдан при оговорках ниже. Это не
реализация/эксперимент/сходимость; Julia/Pkg/H5 и изменения плана не выполнялись.

Независимый вывод: используем matrix spectral **measure** dM, включая atoms,
`g(w)=∫dM(x)/(w−x)`, `Q(w)=∫f(x)dM(x)/(w−x)`. Абсолютно сходящаяся identity

`β⁻¹Σn[(iωn−x)(w−iωn)]⁻¹=[f(x)−f(w)]/(w−x)`

даёт `Q(w)=f(w)g(w)+β⁻¹Σn g(iωn)/(w−iωn)`. При `w±=ζ±Ω` разложение
`L(ζ−iωn)=1/(w+−iωn)−1/(w−−iωn)` даёт
`C=K[Q(w+)−Q(w−)−f(w+)g(w+)+f(w−)g(w−)]`. Подстановка в spectral Fock
`K[(NB+1)g(w−)+NB g(w+)+Q(w+)−Q(w−)]` воспроизводит рамочную формулу.
Для одного ε это
`K[1]·[(NB+1−fε)/(ζ−ε−Ω)+(NB+fε)/(ζ−ε+Ω)]`;
на ζ=iωl shifted coefficients равны0, а L(iν)>0 — actual Mats sign.

Никакая matrix commutativity не использована. Source `applyK/literalK` —
complex-linear raw contractions с qK/wk один раз; отрицательные Mats частоты
дополнительно используют Born property `K[X†]=K[X]†`. Four elastic static
channels дают K[g]; `greens` повторяет central-driven5×5 Σ в9 blocks того же
H45, без factor9/spin2 в self-energy. Для genuinely noncommuting2×2 spectral
oracle нужны≥3 residues: два с суммой I неизбежно коммутируют.

Проверен полный [Margine–Giustino2013, §II.1–II.3](https://arxiv.org/html/1211.3345):
Eq10 normal electron block даёт L=−D; thermal kernel Eq29 и odd continuation
Eq31 дают оба веса, Eq30 — тот же L. Но Eq26–27 уже используют FS/constant-DOS
approximations. Full orbital-matrix result обоснован **самостоятельным
spectral выводом**, а не отменой этих approximations через Δ→0; iπ/NF оттуда
не переносится. [APS abstract MSC1988](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.37.4965)
проверен; полный текст1988 tool недоступен, его equations напрямую не сверялись.

Существенные ограничения: f(w) имеет Fermi poles; отдельно C/shifted terms
имеют artificial poles ζ=±Ω+iωn, отменяющиеся только при согласованных полном
anchor и g. Complex f обязателен. Для первого будущего опыта разумно
`0<η≤ω0/2`, `ω0=π/β`, затем η/2 при fixed qualified Mats μ/anchor; η→0
проверяется отдельно. На real boundary C Hermitian и weights≥0; **при η>0
C не Hermitian**, Herm(C/Bn) недопустима. Точный global spectral functional
causal при positive Born K; fixed-anchor finite-ladder iteration и pointwise
causal eigenvalues этого не доказывают. Mats contraction к малому η не
переносится. Same-H45/μ bare reference и Σ∞=0 оправдывают Δg=O(ω⁻³), paired
anchor tail O(N⁻⁵), но не finite error certificate. Boundary asymptote/M→2M,
anchor cutoff и η — независимые limits. Phases moduloΩ независимы только при
fixed μ/anchor/kernels; одна ladder не задаёт density/sum rule. Старый frozen-Γ
endpoint certificate к новой Mats branch не относится. Новых tests/solvers
этой кампании не предлагается.
