# Статическая проверка M3 density postprocess: уравнение и граница вывода

WRITE-ONLY math review, без Julia/HDF5/новых inner maps/μ-fit. Предназначен для будущего отдельного source-bound postprocessor; это не review его ещё не полученного кода.

Пусть z=iω, h=hH−μI (central local block), bare reference — тот же полный Hermitian H45−μI. Для исходных equilibrium SCBA формул и normalized spectral moments:

`G_c(z)=I/z+h/z²+[(H45rel²)_cc+S0]/z³+[(H45rel³)_cc+h S0+S0 h+S1]/z⁴+…`,

`Σ(z)=S0/z+S1/z²+…`.

Bare reference имеет те же H45 terms, поэтому

`ΔG(z)=S0/z³+B(ρ)/z⁴+…`, `B=h S0+S0 h+S1`.

Ни9, ни дополнительные Tplus terms в B не возникают: bare/full chain terms отменяются; `(H45rel diag9S0)_cc=h S0`, его adjoint term аналогичен. Матрицы не коммутируют — оба ordered products обязательны.

Static direct C даёт `S0static=Kstatic[I]`, `S1static=Kstatic[h]`. LO Fock pole numerator даёт

`S0LO=(2NB+1)KLO[I]`,
`S1LO=KLO[(2NB+1)h+Ω(I−2ρ)]`.

Это следует из суммы residues `(NB+1−f)+(NB+f)=2NB+1` и first moment `(NB+1−f)(ε+Ω)+(NB+f)(ε−Ω)=(2NB+1)ε+Ω(1−2f)`. Следовательно предложенные

`S0=Kstatic[I]+(2NB+1)KLO[I]`,
`B(ρ)=hS0+S0h+Kstatic[h]+KLO[(2NB+1)h+Ω(I−2ρ)]`

**правильны** для infinite equilibrium spectral moments. Единицы: S0 — E0², B — E0³; h, Ω, β в тех же native energy units, μ вычтен один раз.

При Hermitian coefficients `Herm[S0/(iω)³]=0`, а `(iω)^−4=+ω^−4`, поэтому omitted density contribution leading order имеет **плюс** `t4 B`,

`t4=(2/β)Σ_(n=N)^∞ ωn^−4 = β³/(8π⁴) ζ(4,N+1/2)>0`.

Selfconsistent moment postprocess `ρM3=ρN+t4 B(ρM3)` эквивалентен линейной системе

`[I_super+2Ω t4 KLO]ρM3=ρN+t4 B0`,
`B0=hS0+S0h+Kstatic[h]+KLO[(2NB+1)h+ΩI]`.

Знак слева **плюс**. Решать на Hermitian fields/real coefficients либо явно измерять Hermitian defect; density не project/normalize. CP map как superoperator может иметь отрицательные eigenvalues (unitary conjugation — простой пример), поэтому «CP ⇒ I+aK положительно определён» неверно. Достаточная проверяемая well-posedness: `2Ω t4 ||KLO||<1`; затем Neumann inverse bound≤1/(1−2Ωt4||KLO||). Иначе нужна честная condition/residual measure линейной системы, не ссылка на CP. Сохранять linear residual abs/scale/weighted, raw+correctedρ, Pauli eigenvalues/number; не новые acceptance tolerances.

**Важная граница:** сохранённый finite-N LO operator использует outside Gref. Если его reference completion формально продолжить в high Mats n, correction в Σ1 использует **ρN**, а не ρM3: `−(2Ω/β)Σ_inside HermΔG =−2Ω(ρN−ρref)`. Поэтому selfconsistent M3 — новая asymptotic moment closure для оценки infinite tail, а не точная quadrature identity неизменного finite-N nonlinear map. Разница B(ρM3)−B(ρN) влияет на density на O(t4²), но это asymptotic statement без certified constant.

Даже после root finite-N raw residual0 M3 не доказывает infinite-N SCBA: higher moments/remainder, missing interacting LO tail feedback в низкие n и floating CP/operator error не ограничены этим членом. При Hermitian next moments ожидаемый следующий density term имеет порядок N^−5 (если достаточные моменты существуют); коэффициент не сертифицирован. Если 128/256 corrected densities сблизились, это отдельное measured asymptotic witness, не error bar. Fixed μ/dopant density/Poisson/native acceptance остаются отдельными вопросами. Предпочитать стабильную положительную t4 evaluation; total-minus-partial β³/48−partial может терять относительную точность маленького хвоста и требует numerical witness.
