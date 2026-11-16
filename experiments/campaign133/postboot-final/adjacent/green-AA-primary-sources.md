# Green-variable AA: первичные основания и пределы

2026-11-16. **FUTURE / UNRUN.** Только ограниченная проверка литературы и существующего source; Julia/Pkg/HDF5 не запускались, физика/Core не менялись. После текущего FDseed-run новых опытов/реализаций не предлагается запускать без нового поручения. Исход текущего run эта записка не оценивает.

## 1. Что в действительности доказали Walker–Ni

Проверен полный авторский [PDF, SIAM J. Numer. Anal. 49, 1715–1735 (2011)](https://users.wpi.edu/~walker/Papers/Walker-Ni,SINUM,V49,1715-1735.pdf), [DOI 10.1137/10078356X](https://doi.org/10.1137/10078356X). Algorithm AA, (1.1)–(1.2): f_i=g(x_i)−x_i; минимизировать ‖Σα_i f_i‖₂ при Σα_i=1; x_{k+1}=Σα_i g(x_i), β_k=1. Assumption 2.1 и Theorem 2.2 требуют g(x)=Ax+b, несингулярного I−A, полной истории m_k=k и общего x₀ с GMRES для (I−A)x=b. Дополнительно r_{k−1}≠0 и строгое уменьшение GMRES residual на шагах 0<j<k. Тогда

\[
\sum_{i=0}^{k}\alpha_i^{(k)}x_i^{AA}=x_k^{GMRES},\qquad
x_{k+1}^{AA}=g(x_k^{GMRES})=x_k^{GMRES}+r_k^{GMRES}.
\]

Это **не равенство одноимённых iterates**. Remark 2.3 допускает stagnation на самом шаге k; Remark 2.7/Proposition 2.8 показывают последующую остановку AA и rank deficiency. Авторы различают truncated и специально restarted варианты; обычную скользящую историю нельзя называть GMRES(m) без оговорки. Corollary 2.10 относится к линейному left-preconditioned splitting. Конечная история, damping, regularization, nonlinear normalization и cone projection Green-SCBA не удовлетворяют теореме; гарантии сходимости она нам не даёт. Фиксированная положительная inner-product norm допустима; real/imag complex-поля можно realify.

## 2. Проверенный NEGF-пример: смешивается именно G

Первичный полный [Künzel et al., NESSi 2.0, arXiv:2605.02606v1, §3.2 и §5.4.2](https://arxiv.org/html/2605.02606v1) содержит stationary Green iteration. Уравнения (57),(62): Gᴿ=[ω+i0−ε−Σᴿ]⁻¹ и G⁼=GᴿΣ⁼Gᴬ. Последнее использует исчезновение initial-boundary term при decay; exact bound-state poles требуют отдельного treatment. §5.4.2 прямо даёт G_new=αG′+(1−α)G_old; `distance_norm2(G′,G_old)` вычисляется **до** mixing. Пример — scalar Bethe-DMFT, self-consistent perturbation/IPT; он подтверждает выбор Green function как переменной, не QCL-equivalence/AA/cone theorem. Stationary R/< остаются независимыми; equilibrium FD в (55) — отдельное условие, не обязательное ограничение под bias.

Сверен локальный source официального репозитория, HEAD `6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6`: [ness2_bethe.cpp:103–139](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/ness2_bethe.cpp#L103-L139). Порядок: сохранить G_old; построить Σ[G] (при `ipt_flag=0`); K=Σ+G (Bethe Δ=G); Dyson; raw distance; linear mix. SHA256 файла `f300eaee183b6bf5670f6f37077a537d0bc1e9d43206025fe68d1ad5b19e46d4`. Это source evidence, **не** пример точной нашей finite-chain SCBA и не AA-реализация. Никакой перенос их FFT cutoff/bath/нормировки не предполагается.

## 3. Наша гипотеза и отдельные обязательства доказательства

Будущая переменная X=(Gᴿ,G⁼,G⁾), без FD restriction. Требуемая карта
`T(X)=native_paired_normalize(DysonKeldysh(embed(fresh_local_Phi(X), h)))`.
Все пять local channels, interperiod embedding, cavity remaps, energy boundary и исходная number rule сохраняются. Реконструировать GG из A−Gn нельзя, если native pair-normalizer оставляет его независимым. Фиксированный h/bias/model обязателен; дополнительные обновляемые неизвестные нельзя молча исключить.

По отдельной source/math-проверке `equilibrium_oracle` (готовящаяся `green-variable-map-assessment.md`): e171 `production_solver:481–516` строит fresh scattering и embedding именно из него; selected finite-chain operator суммирует fresh local. Точные source anchors и доказательство следует брать из той записки. При raw normalization factor λ=1 full native raw roots и Green roots с λ=1 взаимно восстанавливаются. Малый Green defect или fitted number без λ=1 этого не доказывают. Same nonzero linearized spectra AB/BA применимы только на smooth branch, не автоматически на kink paired normalizer.

История одного state — 3F вместо native 21F (24F включает derived total embedding), F=sizeof(ComplexF64)·NE·Nk·Nb². Для пяти сохранённых state+residual пар — 30F вместо 210F: это оценка **историй**, не полного live-set/RSS. Fresh scattering/embed/Dyson workspaces и строгая native проверка остаются.

**Trial cone projection — наше предложение, источниками не доказанное.** Хранить/считать residual только raw `T(X)−X`; projected AA trial проверять как следующий вход, raw image не eigenclip/project. Даже неподвижность projected шага не принимает raw equations. Нужны отдельные domain/causality/конвенционные проверки и проверка, что допустимый exact root projection не меняет. Финальная приёмка — materialize свежие native поля и полный исходный strict native exit (включая требуемые последовательные fresh unprojected passes); новый residual/малая mixing-step не заменяет его. Здесь это только будущая стратегия, без implementation/convergence claim.
