# Thermal LO: точное mixed-axis продолжение Matsubara anchor

Research-only, 2026-11-16; формулы и будущий опыт, **не реализация/результат расчёта**. Источник модели: `adjacent/matsubara-all5/matsubara_all5.jl`, SHA256 `1434e86fe452f02c5831eaeb0fc83dc1dc12710dcc3d7e81f7f36cd985a31447` (`context`, `reference`, `lo_delta`, `greens`, `applyK`). Fine2x и исходные данные не изменены.

**Предложение:** продолжать Matsubara-решение функциональным уравнением с точными сдвигами ±Ω, без Hilbert-интеграла по действительной энергии. Метод Marsiglio–Schossmann–Carbotte формально точен; доступен abstract [PRB 37, 4965 (1988)](https://doi.org/10.1103/PhysRevB.37.4965), полный оригинал закрыт. Первичный полный текст [Margine–Giustino, PRB 87, 024505 (2013), §II.3, уравнения 26–31](https://arxiv.org/html/1211.3345) показывает mixed Matsubara/real-axis структуру и thermal kernel. Ниже — самостоятельный вывод для normal-state матриц, без их Fermi-surface/DOS/Nambu приближений и множителей iπ.

## Тождество, знаки и независимая проверка

Все энергии в одной native единице: ζ=z_phys−μ, g(ζ)=G(μ+ζ), β=1/kBT, ω_n=(2n+1)π/β, Ω=ħω_LO>0, f(x)=(e^{βx}+1)⁻¹, N_B=(e^{βΩ}−1)⁻¹. Линейная карта K_LO содержит исходные `qK`, `wk` ровно один раз; для одного уровня K[X]=g_c²X. Определим L(u)=2Ω/(Ω²−u²). Тогда

\[
\boxed{\Sigma_{LO}(\zeta)=C(\zeta)+[N_B+1-f(\zeta-\Omega)]K_{LO}[g(\zeta-\Omega)]
 +[N_B+f(\zeta+\Omega)]K_{LO}[g(\zeta+\Omega)]},\qquad
C(\zeta)=\frac1\beta\sum_{n\in\mathbb Z}L(\zeta-i\omega_n)K_{LO}[g(i\omega_n)].
\]

Доказательство: при ρ(x)=A(μ+x)/(2π) положим Q(w)=∫f(x)ρ(x)/(w−x)dx. Абсолютно сходящаяся сумма даёт
`Q(w)=f(w)g(w)+β⁻¹Σ_n g(iω_n)/(w−iω_n)`.
Подстановка в спектральный Fock
`K_LO[(NB+1)g(ζ−Ω)+NB g(ζ+Ω)+Q(ζ+Ω)−Q(ζ−Ω)]`
даёт рамочную формулу. Для g(w)=1/(w−ε) её anchor равен
\[
C/g_c^2=\frac{f(\zeta-\Omega)-f(\epsilon)}{\zeta-\Omega-\epsilon}
 -\frac{f(\zeta+\Omega)-f(\epsilon)}{\zeta+\Omega-\epsilon};\quad
\Sigma_{LO}=g_c^2\!\left[\frac{N_B+1-f(\epsilon)}{\zeta-\epsilon-\Omega}
 +\frac{N_B+f(\epsilon)}{\zeta-\epsilon+\Omega}\right].
\]
Это точное совпадение с заданным независимым oracle, включая знак и β. На ζ=iω_l оба shifted коэффициента равны нулю, а L(iν)=2Ω/(Ω²+ν²): восстанавливается Matsubara convolution. На действительной оси C=C† и коэффициенты ≥0, поэтому Γ_LO=a_-K[A(E−Ω)]+a_+K[A(E+Ω)] совпадает с native equilibrium FD/Bose алгеброй. Никакой коммутативности матриц не использовано.

## Вычислимый anchor и уравнение Dyson

Взять точный bare reference g₀(w)=Σ_aR_a/(w−ε_a) из **того же H45**, Δg_n=g(iω_n)−g₀(iω_n). Для d_f(w,ε)=[f(w)−f(ε)]/(w−ε), с пределом f′(ε)=−βf(ε)[1−f(ε)], считать
\[
C_0=K_{LO}\!\left[\sum_aR_a\{d_f(\zeta-\Omega,\epsilon_a)-d_f(\zeta+\Omega,\epsilon_a)\}\right],\quad
C=C_0+\beta^{-1}\sum_{n\ge0}\{L(\zeta-i\omega_n)B_n+L(\zeta+i\omega_n)B_n^\dagger\},\ B_n=K_{LO}[\Delta g_n].
\]
Это избегает вычитания больших искусственных bare-полюсов между `ΣLO_ref` и shifted Δg; divided differences требуют устойчивого `expm1`/предела. При Σ∞=0 имеем Δg=O(ω⁻³): абсолютный хвост O(Nω⁻⁴), после сопряжённого спаривания с Hermitian moments — O(Nω⁻⁵) при ограниченном ζ. Это асимптотики, **не конечный сертификат ошибки**; проверка tail/μ текущего Mats anchor обязательна.

Статические четыре канала: Σ_stat(ζ)=Σ_sK_s[g(ζ)]. Полная Σ_local=Σ_stat+Σ_LO повторяется во всех девяти 5×5 блоках; g — центральный блок `[(μ+ζ)I45−H45−I9⊗Σ_local]⁻¹`. Это тот же **finite45 / central-driven local5 / all5 direct-static** closure, не независимая cellwise45 SCBA и без множителя 9. Полный тепловой LO-Fock, фиксированный bath и hH; нового Hartree/tadpole нет.

На лестнице ζ_j=ξ+jΩ+iη нужны только соседи j±1; разные фазы ξ modulo Ω независимы. При комплексном ζ использовать **f(ζ±Ω)**, не f(Re ζ±Ω). Для первого опыта 0<η<π/β избегает искусственных Fermi-полюсов; η→0 проверяется отдельно. На конечных концах лестницы — bare/moment asymptote, затем удвоение длины, не нулевые ghost-значения.

```julia
# Только высокоуровневый алгоритм; не исполняемый новый solver.
ref = same_H45_reference(model, mu)
B = KLO(Gmats .- ref.Gmats)            # acts on each positive-n matrix/k slice
z = xi .+ (-M:M) .* Omega .+ im * eta
C = fermi_divdiff_anchor(ref, z) .+ paired_mats_anchor(B, z, beta)
g = bare_central45(ref, z)
for it in 1:cap
    lo = KLO(g)
    minus, plus = shift_one_node(lo; boundary=moment_asymptote)
    S = sum(Kstatic(s, g) for s in STATIC) .+ C
    S .+= (NB .+ 1 .- f.(z .- Omega)) .* minus
    S .+= (NB .+ f.(z .+ Omega)) .* plus
    phi = central45_resolvents(model, mu .+ z, S)
    record(raw_defect(phi, g), causal_witness(phi))
    g = mix(g, phi)                    # raw criterion before mixing; no PSD clipping
end
# Independently vary Nomega, M, eta; keep Mats mu fixed.
```

## Что различит следующий опыт

Это **иная identity**, чем finite-window PI-Hilbert Γ: direct-static K[G] сохраняет виртуальные хвосты, LO использует полную спектральную меру, а прежний путь обрезает/гасит её за окном. Γ-алгебра на границе совпадает, retarded closure целиком — ещё нет. Mats μ/number не заменять real-axis перенормировкой. Mats conditional q≈0.305 не доказывает contraction около действительных резонансов; при η>0 комплексные коэффициенты не дают простого CP-доказательства итерации. Ветвь должна быть causal и согласована с anchor/asymptotics; конечная ladder сама этого не доказывает. При T→0 LO-связь направлена от Fermi-области наружу; при 70 K малое NB позволяет outward sweep как preconditioner, **не разрешает отбросить absorption**.

После кампании: сначала exact одноуровневый и некоммутирующий 2×2 spectral oracle; затем одна Ω-лестница через residual anchor ~82.4 meV с фиксированным квалифицированным Mats input, сравнения Nω=128/256, M/2M и двух η. Измерить raw Dyson/Fock defects, causal eigenvalue witness, чувствительность к anchor-tail/границам; сходимость не обещана. Для Nω=256, Nk=11, Nb=5 один Mats/257-node field ≈1.1 MiB, несколько полей <20 MiB плюс kernels (~10 MiB), без оценки JIT/RSS. Один sweep 257×11=2827 полных 45×45 resolvent solves, 16 sweeps — 45232; anchor O(Nω·L·Nk·Nb²) однократно. Для полного fine спектра шаг Ω/1878 всё ещё требует многих фаз: экономия на отдельных запросах правдоподобна, экономия полного спектра/плотности **не установлена**.
