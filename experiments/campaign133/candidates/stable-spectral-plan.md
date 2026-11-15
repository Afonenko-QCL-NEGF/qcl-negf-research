# Stable spectral A: proof и различающий frozen probe

**Факт:** warm attempt01 остановлен на первом Π в `fit_mu`, до μ fit/
scattering/native validation: `nonpositive/nonfinite raw spectral weights`.
Прочитаны его summary/log; source Γ-first/model сохранены. q-bad indices,
их величины и raw blocks пока **not_measured**. Сообщённый positive seed
minimum1.2585e−27 при block Frobenius3.6286e−27 делает precision hypothesis
разумной, но не устанавливает причину STOP и не доказывает true PSD всех
blocks. Здесь proof/design only; Julia/Pkg/H5 arrays/compute не запускались,
Core/harness не менялись.

**Важная поправка:** native q использует только
`Aaa=real(i*(Gaa−conj(Gaa)))=−2Im Gaa`. Финальная diagonal subtraction
не вычитает два больших почти равных imaginary numbers: они складываются.
Ошибка могла возникнуть раньше — complex inverse/cavity products,
anti-Hermitian Λ/h, inputs — а не обязательно в последней операции
`i(G−G†)`. Off-diagonal extraction действительно subtracts близкие entries.

## 1. Exact identity и положительность

При real E, Hermitian coherent H, retarded Σ=Λ−iΓtot/2 с Hermitian Λ,
Γtot≽0 и nonsingular D=EI−H−Σ:

```
G−G† = G(D†−D)G†,
A=i(G−G†)=G Γtot G† ≽0.
```

Для full45 open finite-chain4 при local Γ repeated:
`A_c=P G45 (I9⊗Γ) G45† P†`. Никаких eta floor, clipping или bath
occupation для этого algebra не требуется. При singular undamped states
identity на sampled inverse неприменима; delta spectral measure не
восстанавливается умножением Γ=0. Это отдельный quadrature/atom вопрос.

`G Γ G†` избегает extraction of tiny spectral weight из O(1) coherent
complex entries. Но computed generic matrix congruence не является
**строго positivity-preserving floating primitive**. Для true PSD Γ=B B†
Gram form `(G B)(G B)†` и diagonal `sum abs2` сохраняют nonnegative trace
по арифметической конструкции. Factorization должна пройти без eigenvalue
clip/floor. Cholesky strict-PD case подходит; failure не означает право
добавить η. Near-semidefinite case требует qualified factorization/precision,
а не неявного spectral repair.

## 2. Native cavity spectral probe: что доказывает f=1

Actual accepted Core `e171caa60cce676c6c4238dd671e7c527f2bd977`,
`optimized/production_operators.jl::_cavity_side!` выполняет
local< + C g< C†, затем gR(...)gA и финальный C(...)C†.
Сохранить ту же **local retarded** family, подать технически

```
Σlocal< = i Γ,    Σlocal> = 0.
```

На внешнем конце cavity embedding0. Induction inward даёт

```
A_j = g_j [Γ+C A_(j−1) C†] g_j†,
g_j< = i A_j,  g_j> = 0,
Γemb,s = −i Σemb,s< = C_s A_s,4 C_s†.
```

C+=T+, C−=T+†; сумма двух embeddings плюс central Γ даёт
`Aprod=Gc(Γ+Γemb,+ +Γemb,−)Gc†`. Это exact full45 Schur-complement
identity. **f=1 — spectral test source**, не physical Fermi occupation,
не electronic density, не initialization at infinite μ. Не вызывать на
нём scattering/number normalization/paired constraint. Native cavity R
от этих L/G не зависит, поэтому same-input R comparison допустим.

Native covariance propagation сохраняет масштаб Γ в lesser products;
оно не извлекает Γemb из subtraction больших R entries. Но roundoff
congruences может дать tiny non-Hermitian/negative matrix eigenvalue.
Если нужен гарантированно nonnegative q, trace вычислять через qualified
Gram factors, без `max(q,0)`. Full A для FD/Φ и q должны иметь **один**
source. Старый `_seed_number_weights_production(GR)` снова извлечёт
−2ImdiagGR и не станет stable от того, что рядом создан другой A.

## 3. Fake broadening и где нужна Hermitian consistency

Если Λraw не Hermitian, local retarded broadening равна

```
i(ΣRraw−ΣRraw†)=Γ+i(Λraw−Λraw†).
```

FFT columns независимы; exact real PI operator сохраняет adjoints, но
computed complex columns могут нарушить это на масштабе eps||Λ||.
Такой absolute defect может быть намного больше physical1e−27 Γ.
Central raw h тоже надо проверить: native cavity plan symmetrizes hlocal,
а central Dyson использует переданный h. Coherent contribution
`i(h−h†)` столь же важен. T−/T+† должны быть согласованы.

Более того, даже при exact-Hermitian local Λ обычный floating
`C gR C†` в retarded cavity может получить anti-Hermitian noise, не
согласованный с Γemb из spectral probe. Поэтому **Herm(local Λ) alone
не является полным proof стабильности всей floating finite-chain map**.
На одном raw GR product with physical ΓK PSD может быть положителен,
хотя он не равен spectral function того же corrupted raw denominator.

Для будущего consistent arithmetic variant каждая Schur step может
разделять coherent Hermitian и spectral sources:

```
Λemb = Herm(C g_prev C†),   Γemb=C A_prev C† (qualified Gram),
g=[E I−Herm(h)−Herm(Λlocal)−Λemb+i(Γ+Γemb)/2]^−1,
A=g(Γ+Γemb)g†.
```

Это exact-equivalent algebra при exact input assumptions, но новый
floating implementation variant. Hermitian real-part representation
enforcement не clipping Γ и не new physical linewidth; raw correction
norms/blocks сохранить. В production нельзя молча объявить его прежним
bitwise native map. Сначала frozen probe; не обещать SCBA convergence.

## 4. Same-G decomposition отделяет два механизма

Для **зафиксированных computed arrays** Ghat,D (без смены GR):

```
ΓD=i(D†−D), ΓK=Γlocal+Γemb,probe,
Rright=I−Ghat D,
Asub=i(Ghat−Ghat†), AR=Ghat ΓD Ghat†, AK=Ghat ΓK Ghat†,
Asub−AK = Ghat(ΓD−ΓK)Ghat†
          + i(Ghat Rright†−Rright Ghat†).
```

Это exact algebraic identity; её вычисление имеет свой FP residual.
Первый term — retarded/physical broadening inconsistency (Λ/h/cavity).
Второй — inverse residual/extraction consistency, norm≤
`2||Ghat||2 ||Rright||2`, плюс rounding при evaluation. Она значительно
сильнее сравнения двух A без denominator audit. Сохранить signed traces,
eigenvalue witnesses, norms/scales и right/left inverse residuals.

Forward bound для стабильного product всё ещё нужен. Если ||δG||≤eG,
||δΓ||≤eΓ, g=||Ghat|| и Γscale — upper bound intended ||Γ||:

```
||δA|| <= Γscale*(2g eG+eG²)+g² eΓ+e_product.
```

Обычный inverse backward error может давать eG≈g² eD/(1−g eD),
если g eD<1. Вблизи pole это может быть большим. Gram positivity сама
по себе не доказывает accurate resonance position/weight/quadrature.
Для direct extraction error bound содержит O(eG), без малого Γscale;
это объясняет возможную огромную relative error tiny A.

## 5. Минимальный следующий frozen diagnostic — после resource-green

1. На exactly Γ-first/model/input attempt01 воспроизвести **одну R-only**
   evaluation, без μ fit/Φ/scattering iterations. Сохранить всю q vector,
   count negative/nonfinite, min/signed sum/weighted negative activity,
   e_minq и worst negative per-k trace на этих energies. Добавить несколько
   occupied/control points, максимум8 selected (E,k). Не выбирать indices
   по желаемому объяснению; current q-bad indices ещё неизвестны.
2. Один same-input native cavity spectral test Σ<=iΓ,Σ>=0. На selected
   points сохранять raw/h-Herm/Λ-Herm defects, Γinput, ΓD, ΓK, raw GR,
   Asub/AR/AK, predicted two terms, inverse residuals. Цель — определить,
   что доминирует, а не автоматически заменить A и продолжить solver.
3. Small selected full45 BigFloat256 oracle на **promoted same Float64
   data**: (a) raw Λ и actual native coherent blocks сохраняют fake
   broadening: central hraw, neighbours именно plan.hlocal (native уже
   symmetrizes их); (b) отдельно Herm(h) во всех cells, Herm(Λ) с тем же Γ.
   Full45 open H uses actual T/T†. В (b) сравнить
   direct high-precision i(G−G†) с `P G(I9⊗Γ)G†P†`/Cholesky Gram.
   Проверить Γ positive factorization без clipping и selected inverse/
   factor residuals. BigFloat384 repeat на худшем point — precision witness,
   не automatic interval certificate. Не строить full45 на всех30704×11.
4. Интерпретация: (a) reproduces negativity, (b) positive → corrupted
   coherent/retarded input contribution; same-data BigFloat positive,
   Double Asub negative → inverse/cavity arithmetic; Γ itself fails
   qualified positivity → input/cone issue. Возможна смесь механизмов.
   Если BigFloat residual/conditioning не позволяет решить — unresolved.

Scalar independent oracle: real ξ, Γ=γ>0, G=1/(ξ+iγ/2),
`A=γ/(ξ²+γ²/4)>0`, including γ=1e−27 and off-resonant ξ. Для complex
Hermitian manufactured chain с comparable coherent entries проверить
inverse-noise path отдельно: scalar case не воспроизводит matrix phase
cancellations. Никакого density renormalization, negative-q clipping,
Γ floor/η или tolerance relaxation ради прохождения fit.

Original raw first-Π STOP остаётся исходом baseline attempt. Даже если
stable A исправляет q, current μ/capacity/Φ/native strict audit должны
быть выполнены заново в явной variant identity. Это numerical arithmetic
change одной exact model, не уже принятый full solution.

## Первичная опора

[Anantram–Lundstrom–Nikonov, *Modeling of Nanoscale Devices*,
author paper](https://arxiv.org/pdf/cond-mat/0610247), Eqs156–158,
выводит spectral broadening identity и covariance correlations;
матрицы/порядок products сверены с нашим inverse proof, а не перенесены
по одной нотации. [Lake–Klimeck–Bowen–Jovanovic, JAP81,7845–7869(1997),
primary PDF](https://engineering.purdue.edu/gekcogrp/publications/pubs_src/DOC142_JApplPhys_81_7845.pdf)
задаёт recursive open-system/scattering NEGF treatment. Эти работы не
доказывают Float64 accuracy данного Core. Matrix residual/Gram arguments
и разделение двух numerical mechanisms выше выведены здесь независимо.
