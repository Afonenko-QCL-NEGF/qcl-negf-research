# Γ64: строгая локализация atoms/ghosts — feasibility, без запуска

**Да, полезный upgrade возможен без ещё одного Float eigen witness.** Вместо
повторения V4 dL/dR предлагается componentwise outward enclosure endpoint
Cauchy matrices и Gershgorin lower gaps для full45. Положительный enclosure
gap доказывает знак на exact stored-input model при явно оговорённой IEEE
arithmetic; неположительный означает **unavailable**, не atom и не |gap|.
Big192 eigensolve без residual+input bound сам по себе таким доказательством
не был бы. Actual Γ/H5 автор не читает; existing V4 mass≈1e−33 уже есть.

Пусть `Σ(z)=Σj Γj k((z−cj)/h)` — те же full positive triangular hats,
`Γj≻0` по принятому ALL exact-PD certificate, saved exact-Hermitian hH во
всех9 cells; H45 имеет upperT+, lowerT+†. Support `ℓ=c1−h,r=cN+h`.
Вне support `Σ(E)=Σ(E)†`, `Σ′(E)=−∫Γhat(t)/(2π(E−t)²)dt⪯0`, поэтому
`D′(E)=I45−I9⊗Σ′(E)⪰I45`. Если

`Cℓ=H45+I9⊗Σ(ℓ)−ℓI45≻0`, `Cr=rI45−H45−I9⊗Σ(r)≻0`,

то `D(E)≺0` для всех E≤ℓ и `D(E)≻0` для всех E≥r. Endpoint Σ конечна,
поскольку taper линейный; логарифмически бесконечная производная не мешает
интегрировать order inequality снаружи и брать one-sided limit. Ни exterior,
ни endpoint real poles/atoms тогда нет. ALL PD и перекрывающиеся adjacent
hats исключают interior singular mass; causal compressed resolvent имеет
total central massI5. Математическая опора — spectral dilation/вывод в
`../../pi-causality-window-certificate.md`,
[Gesztesy–Tsekanovskii, matrix Herglotz functions](https://arxiv.org/abs/funct-an/9712004),
[первичная публикация](https://doi.org/10.1002/1522-2616(200010)218:1%3C61::AID-MANA61%3E3.0.CO;2-D).

**FP geometry существенно:** V4 использует Float64 `cj=canonical[j]`, а не
точные real `a+(j−1)h`. Сертификат должен брать эти dyadic centres/h точно;
поддержка и distances строятся exact rational. Для t>1 положительный ряд
`k(t)=π⁻¹Σm≥0 t⁻(2m+1)/[(2m+1)(2m+2)]` имеет geometric upper tail.
Для t=1 `k=log2/π`. π и log2 ограничиваются exact rational Machin/atanh
series, затем проверяемыми dyadic enclosures. Real interval +,−,×,÷ с
prevfloat/nextfloat охватывает coefficient и accumulation errors; complex
entries — две real interval. Нет libm log, FFT, eigenvalue acceptance или
новой dense array. Требования: finite intermediates, IEEE basic operations,
gradual underflow, без fastmath/FMA reassociation. Это explicit numerical
domain, не end-to-end hardware/software formal verification.

Для native window[A,B] оценить uniform strip Hermitian gap через
`BL=Σj bjL Γj`, где `bjL=k((cj−A)/h)` при distance>1, иначе `bjL=1/π`;
аналогично BR при `(B−cj)/h`. Глобально `|k(t)|≤asinh(1)/π<1/π`, поэтому
`H45−AI−I9⊗BL⪰dL I`, `BI−H45−I9⊗BR⪰dR I` дают `||G45||≤1/d`.
При d>0 missing central mass ограничена

`βghost≤(Σj TrΓj [area_left_j/dL²+area_right_j/dR²])/(2π)`.

Areas — **exact** интегралы unit hats вне[A,B], включая возможное1ulp
перекрытие второй/предпоследней hat. Для идеальных endpoint strips это
прежний `h||Γend||/(4πd²)` (trace вместо norm только ослабляет bound),
никакого factor9. Sheet-number bound `gs Nb Σk wk βk`, optional f≤1.
Conservative positive interval gaps/β отделят atoms/ghosts от≈9.34e−5target
density difference; они **не** сертифицируют Mats tail/Float inverse/real-axis
quadrature errors. N512→1024 witness не заменяет эти bounds.

Предлагаемый diagnostic: один V4 source_frame (один Γarray, source hashes/
ALL-PD/checkpoint equality/hH), O(NE·11·25) interval pass,44 tiny gap bounds,
small TOML/H5; J1BLAS1/2GiB/120s. Никаких Green/Φ/fit/SCBA. При Gershgorin
failure — stop/unavailable; новый inverse/eigen framework не добавлять.
