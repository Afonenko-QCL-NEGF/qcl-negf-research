**Математика frozen boundary probe; write-only v2**

Для каждого исходного `(E,k)` фиксируются `h(k)`, `T=H[j,j+1]=Tplus`,
`Σlocᴿ(E,k)` и `D0=EI−h−Σlocᴿ`. Все клетки одинаковы, поле и Hartree равны
нулю. Γ получена суммой пяти LOCAL source retarded families, embedding
исключён; Σlocᴿ восстановлена accepted Core uniform linear-hat PI. Здесь
не вычисляется scattering map, не меняются μ, Γ, сетка или baths.

**Поверхности и центр.** Правая поверхность `gR=[D0−TgRT†]⁻¹`, левая
`gL=[D0−T†gLT]⁻¹`; центр бесконечной цепи
`G∞=[D0−TgRT†−T†gLT]⁻¹`. Для конечной цепи с `L` клетками по сторонам
применяются те же рекурсии от `gR,0=gL,0=0`, затем центральный inverse.
Таким образом конечные 4/8/16/32 computed fresh, без reuse bulk как surface.
В zero-field однородной задаче конечный `L` — одинаковый для всех E/k.

**Doubling / López Sancho.** Обозначим renormalized inverse bulk `B`, две
surface inverse `SR,SL`, эффективные hoppings `α,β`. Начальное состояние
`B=SR=SL=D0`, `α=T`, `β=T†`. В итерации один pivoted LU factorization B
и один in-place solve с `2Nb` RHS дают `X=B⁻¹α,Y=B⁻¹β`. Затем одновременно:

```
cR = αY; cL = βX
SR -= cR; SL -= cL; B -= cR+cL
αnew = αX; βnew = βY
```

После первого шага βnew **не** принимается равной αnew†: retarded D0
комплексен и диссипативен. Обновление обеих ориентаций вычисляется отдельно.
Стоп требует `max(||cR||F,||cL||F,||αnew||F,||βnew||F) <= atol+rtol*scale`,
где `scale=max(||D0||F,||SR||F,||SL||F)`. `atol=1e−13 eV/E0`,
`rtol=1e−12`; максимум 64 шага. Это численный stop энергии, не SCBA criterion.
No hidden η, clipping, symmetrization of physical source, or positivity repair.

В `algorithms.jl` плотные рабочие матрицы фиксированного размера выделяются
один раз на worker/k; ядра используют `mul!`, `lu!`, `ldiv!`. Штатный `lu!`
может выделять ограниченный pivot/metadata вектор; это не утверждение о zero
allocations. Floating diagnostics используют временные 5×5 матрицы.
Doubling стоит `O(iter Nb³)`, finite sweep до L32 — `O(32 Nb³)` на E/k,
память ядра — `O(Nb²)`; полные GR состояния не сохраняются. Native c4
alignment временно хранит один full A и штатные native buffers.

**Независимые свидетели.** Для каждой поверхности измеряется исходное
уравнение, например
`rR=||(D0−TgRT†)gR−I||F/(||D0−TgRT†||F ||gR||F+sqrt(Nb))`.
Этот residual безразмерен; делитель имеет ненулевой абсолютный член.
Тот же backward residual измеряется для центрального inverse. Диагностический
порог поверхности `1e−9`, отдельно κ центрального и surface denominator.
Малый backward error при большой κ не устанавливает малую forward error.
Причинность: `A_surface=i(g−g†)`; eigenvalue минимум и минимум, делённый на
`||g||F`, сохраняются. Negative relative witness ниже `−1e−10` помечается как
неразрешённый. Это порог floating диагностики, не физический допуск SCBA.
Все узлы с fail считают отдельно по stop/equation/causality; approximate
infinite integrals при них маркируются unresolved, не принимаются молча.

Scalar oracle использует causal analytic branch
`g=2/(z+sqrt(z−2|t|)*sqrt(z+2|t|))`, `Im z>0`, включая оба внеполосных
сектора. Для complex noncommuting 2×2 toy blocks finite recursion сравнивается
с independently assembled full block matrices L=1/2/4; infinite centre — с
Bloch integral. Это проверяет ориентации, а не только scalar commute case.

Дополнительно на actual source выбираются по одному occupied peak при
k-index 1,6,11: максимум `f(E;source μ,70K) Tr A_native_c4(E,k)`. На каждом
из этих трёх блоков сравниваются finite4 и independently inverted full45,
а `G∞` — с
`G_Bloch=(2π)⁻¹∫₀²π[D0−Te^{iθ}−T†e^{−iθ}]⁻¹dθ`.
Используется periodic trapezoid θ=2048→4096→8192, пока relative refinement
не больше `1e−7` и difference с decimation не больше `2e−7`. Оба числителя и
матрицы сохраняются; если refinement не достигнут, witness unresolved.
Это несколько actual E/k, не дорогая сетка NE×Nθ и не интерполяция Σ по E.

**Интегралы.** По исходным positive weights вычисляются per-k
`Mk=ΣE wE A/(2π)`, `Fk=ΣE wE f A/(2π)`,
`r_sum=maxk||Mk−I||F/sqrt(5)`, `N=gs Σk wk Tr Fk`. N не нормируется,
μ не подгоняется. Все M/F матрицы, числители weighted L2 differences для A,
fA, соответствующие reference norms и discrete triangle bounds сохраняются.
Для 7 выбранных пар есть finite4 reference и последовательные 4→8→16→32→∞.
Треугольные bounds относятся только к выбранной discrete measure; это не
continuum/grid/tail error estimates. Значения в новых boundary variants не
принимают старый finite_chain4 SCBA или бесконечную interacting задачу.

**Первичные опоры и provenance алгоритма.** Математический метод:
[López Sancho, López Sancho, Rubio, J. Phys. F 15, 851 (1985)](https://doi.org/10.1088/0305-4608/15/4/009).
Код здесь написан самостоятельно по приведённой recurrence, не скопирован из
GPL. Фактическое наличие схемы сверено с official SIESTA
`71c8602911fd632fa3331f2588ea55b4d6ed684e`,
[`SSR_sGreen_NoDOS`, m_ts_electrode.F90:428–637](https://gitlab.com/siesta-project/siesta/-/blob/71c8602911fd632fa3331f2588ea55b4d6ed684e/Src/m_ts_electrode.F90#L428):
один `zgesv` на 2Nb RHS (561), независимые α/β updates (571–574),
surface/bulk corrections (576–583). Его energy-absolute stop дополнен здесь
исходными surface equations. Foreign code не импортируется/не выполняется.
Core comparison: accepted
`e171caa60cce676c6c4238dd671e7c527f2bd977`,
[`finite_chain_green`, reference/periodicity.jl:172](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/reference/periodicity.jl#L172),
[`_cavity_side!`, production_operators.jl:235](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/e171caa60cce676c6c4238dd671e7c527f2bd977/src/numerics/optimized/production_operators.jl#L235).

Complex-contour density/Ozaki продолжение — отдельная работа Quantica agent;
данный probe его не выполняет и не смешивает с boundary sensitivity.

Дельта v2 по static review: дополнительно tiny scalar `Im(z)=1.3e−15`
с обеих внеполосных сторон и complex 2×2 `Γ=3e−15 vv†` (rank one, без
скрытого η) при E вне bound `||h||+2||T||`. Они сравниваются с analytic
branch / full finite matrix / Bloch2048→4096 и surface equations; все
свидетели должны быть finite. Это qualification weak/semidefinite damping
вне полюсов; оно не устанавливает точность неразрешённых in-band резонансов.
Итоговый `infinite_numerical_status` вычисляется только после selected
oracles: нужны все k, finite native/independent проверки, конечные интегралы,
ноль surface/finite/nonfinite failures, все selected stops и Bloch
refinement/agreement. Отдельный `surface_numerical_status` не означает
общий pass. Первое θ-значение refinement помечено `refinement_measured=false`
(NaN как отсутствие предыдущей сетки), а не принято за measured defect.
