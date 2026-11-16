# Условная contraction-оценка для finite-N reference completion

Научная область: fixed μ, β, Ω и Hermitian finite9/45 geometry; static4 direct C, thermal LO Fock; положительные Mats frequencies n=0,…,N−1 и отрицательные через adjoint. Никаких density/μ-fit, real-axis continuation или исходной PI/native acceptance. Вывод точный **при CP предпосылке на K**, не certified theorem для сохранённых Float kernels. Julia/HDF5/новые тесты не выполнялись.

Пусть `Ks=Σstatic4 Ks`, `KL=KLO`, I=(I_nb) по k, `cB=2NB+1=coth(βΩ/2)`, `ω0=π/β`. Область `C_N={Σ_n(k): ImΣ_n(k)≼0}`, норма `||Σ||∞=max_(n,k)||Σ_n(k)||op`; это замкнутое полное пространство. Предполагаем, что все Ks/KL — complex-linear CP maps на direct sum k matrices, сохраняют adjoint; qK/wk включены один раз и положительны. Это физически мотивировано kernels, но не доказано точно для raw Float tensors этой кампании.

**Reference completion.** В `reference + lo_delta` суммарный LO map алгебраически равен

`ΦLO,n = β⁻¹ Σ_(j inside) D(ωn−ωj) KL[G_j(Σ)] + β⁻¹ Σ_(j outside) D(ωn−ωj) KL[Gref_j]`,

где `D(x)=2Ω/(x²+Ω²)` и outside **фиксирован**. Следовательно его variation/derivative outside равна0. Это утверждение finite reference-tail map, не подстановка нулевого interacting tail и не exact infinite SCBA.

**Causal self-map.** Для каждого positive n full resolvent `R_n=[(μ+iωn)I45−H45−diag9 Σ_n]⁻¹` имеет `Im R_n=−R_n[ωn I45−diag9 ImΣ_n]R_n†≼0`, central compression G_n тоже causal. В LO сгруппировать positive j и его negative partner: imaginary part равна `(d−−d+) KL[Im G_j]/β`, где `d−=D(ωn−ωj)≥d+=D(ωn+ωj)≥0`. Это верно и для inside interacting G, и для outside bare causal Gref. Static maps сохраняют знак по CP. Значит Φ отображает C_N в себя в точной арифметике.

**Resolvent Lipschitz без9.** Accretivity даёт `||R_n||≤1/ωn`; resolvent identity и central compression дают

`||G_n(Σ)−G_n(T)||op ≤ ||Σ_n−T_n||sup-k / ωn² ≤ ||Σ−T||∞/ω0²`.

Причина отсутствия9: `||diag9(δΣ)||op=||δΣ||op`, а central projection имеет norm1. Nine-cell repetition не означает сумму девяти norm bounds.

**Combined row norm.** Расширим входы временно до независимых positive/negative Mats matrices. Static+LO inside row есть CP map T_n,N, для которого

`T_n,N[I]=Ks[I]+W_n,N KL[I]`,
`W_n,N=β⁻¹ Σ_(j=0)^(N−1) [D(ωn−ωj)+D(ωn+ωj)] ≤ β⁻¹ Σ_(l∈Z) D(2πl/β)=cB`.

Последняя сумма — elementary bosonic Mats coth identity. Для CP T, `||T||=||T[I]||`: применить T к PSD block matrix `[[I,X],[X†,I]]` при ||X||≤1. Поскольку KL[I]≽0, Loewner order даёт `||T_n,N[I]||≤L`, где **`L=max_k ||Ks[I](k)+cB KL[I](k)||op`**. Это combined norm из cp_bound, а не более грубая сумма отдельных norms. Adjoint-constrained входы — подмножество расширенной области, поэтому bound сохраняется.

Получаем **`||Φ_N(Σ)−Φ_N(T)||∞ ≤ q ||Σ−T||∞`, `q=L/ω0²`**, для любого N, включая finite128/256 с fixed exact reference completion. При q<1 exact CP map имеет единственный finite-N fixed point и Picard α1 сходится к нему из causal Σ0; mixing0<α≤1 имеет bound1−α+αq. Absolute a-posteriori error bound `||Σ−Σ*||∞≤||Φ_N(Σ)−Σ||∞/(1−q)` требует именно этой абсолютной norm и точного/rigorously enclosed map, не произвольного относительного Frobenius residual.

Root сообщил cp_bound ratio≈0.305293<1. **Это conditional explanatory bound:** CP/adjoint preservation raw Float tensors и interval enclosure L/roundoff не сертифицированы, а native `Hermitian` wrapper/observed floating residual0 не заменяют их. Поэтому нельзя объявлять machine-precision actual theorem или continuum/canonical uniqueness. Тем не менее scope прежнего текста «infinite only» излишне узок: тот же условный bound применим к finite-N frozen-reference map. Он не контролирует N→∞ error, density tail, μ-fit derivative или экспериментальную валидность.

Code anchors: module1434e86… `reference`, `lo_delta`, `greens`, `cp_bound`, `pilot`; accepted Core e171 `production_operators.jl::_production_static_contraction!`; independent native six-index qualification и выполненный onelevel tiny относятся к operator/sign wiring, не CP certificate. Основание — приведённые block-positivity/resolvent/coth derivations; никаких внешних общих MDE предпосылок не импортируется.
