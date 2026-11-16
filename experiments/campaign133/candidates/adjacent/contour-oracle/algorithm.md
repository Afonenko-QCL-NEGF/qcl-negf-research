# Off-axis density той же frozen PI-модели

Алгоритм: **Γ → аналитическая Σ(z) → central resolvent полной цепи → Matsubara sum + моменты**. Нет Hilbert FFT, Core Dyson/cavity helper, real-axis energy integration, fitted μ или нового SCBA. Это independent observable integration одной неизменной frozen numerical model, а не новая self-consistency схема.

Для узлов E_j=a+(j−1)h и полноразмерных endpoint hats φ_j(E)=max(1−|E−E_j|/h,0), включая ghost strips [a−h,a], [b,b+h], задаётся Γ̃(E)=Σ_jΓ_jφ_j(E). Тогда

    Σ(z)=∫ Γ̃(E)/(2π(z−E)) dE=Σ_j c((z−E_j)/h) Γ_j,
    c(t)=[(t+1)log(t+1)−2t log t+(t−1)log(t−1)]/(2π).

Здесь **нет дополнительного h или trapezoid wE**. h входит в t. Principal logs при Im t>0 дают нужную аналитическую ветвь. Production coefficient вычисляется cancellation-free серией

    c(t)=(1/t+1/(6t³)+1/(15t⁵))/(2π) + R(t),
    |R(t)|≤1/[56π |t|⁷(1−|t|⁻²)], |t|>1.

Доказательство: разложить 1/(t−u) при |u|≤1; нечётные моменты triangle равны нулю, чётные μ₂m=2/[(2m+1)(2m+2)]. Первый отброшенный μ₆=1/28, последующие меньше; сумма модулей ограничена geometric series. Это bound **аналитической truncation коэффициента**, не всей Float64 программы. При 70 K min Im(t)=πkBT/h≈485. BigFloat256 exact F-second-difference отдельно проверяет коэффициенты на representative complex t; сохраняются mathematical remainder и фактическое Float64 отклонение, без заявления interval certification.

Σ(z) вычисляется chunked contraction `transpose(C)*reshape(Γ,NE,:)`, где C размера NE×32. Не используется adjoint(C): веса не должны conjugate. По каждому k строится H45 с девятью одинаковыми source hcell и source T+,T+†. D=zI45−H45−I9⊗Σ(z). LU решает D X=Bcentral (45×5), Gc=X[central,:]. Используется **Herm(Gc)=(Gc+Gc†)/2**, не `real.(Gc)`.

    ρc = ½ I5 + (2/β) Σ_(n≥0) Herm Gc(μ+iω_n),
    ω_n=(2n+1)π/β, β=1/kBT.

Эта spectral identity учитывает непрерывную и атомарную спектральную меру. Она не требует заранее исключать bound states. Для causal normalized G она следует из scalar Fermi Matsubara partial fractions и spectral representation. Σ — aggregate пяти локальных irreducible channels; g_s/wk отсутствуют внутри Σ.

Полные hats имеют local self-energy moments S0=hΣ_jΓ_j/(2π), S1=hΣ_jE_jΓ_j/(2π). Для central block resolvent:

    M1=H45_cc,
    M2=(H45²)_cc+S0,
    M3=(H45³)_cc+hcell*S0+S0*hcell+S1,
    a1=M1−μI, a3=M3−3μM2+3μ²M1−μ³I.

Отсюда Herm Gc(μ+iω)=−a1/ω²+a3/ω⁴+O(ω⁻⁶). После N положительных Matsubara frequencies возвращаются оба estimator:

    ρM1(N)=½I+(2/β)Σ_(n<N)HermG−(2/β)a1 T2(N),
    ρM3(N)=ρM1(N)+(2/β)a3 T4(N),
    T2=(β/π)²[π²/8−Σ_(n<N)(2n+1)⁻²],
    T4=(β/π)⁴[π⁴/96−Σ_(n<N)(2n+1)⁻⁴].

Tail scalars вычисляются BigFloat256, затем округляются в Float64, чтобы исключить cancellation малой разности. a3 содержит hcell*S0+S0*hcell: пропуск этих terms даёт неверный high-frequency tail при energy-dependent Σ. Для Γ=0 moments получаются из H45 и проверяются exact eigenvector f(H45). Стабилизация при N=128/256/512/1024 — numerical refinement witness; разность двух оценок **не формальный error bound**.

Сохраняются ρM1/ρM3/raw partial sum на каждом N для всех 11 k, moments, source/basis/scales hashes, kernel witnesses, selected Σ/G/D backward errors, сравнительные primary norms, sheet N=g_sΣwkTrρ. Реальная Π reference — `fresh_PI/F` из frozen-control; nodal и archive имеют самостоятельные labels. Γ PSD/absence of atoms не принимаются заново этим harness; см. `../../pi-ghost-mass-bound.md`. Разность contour−real может содержать ghost-strip weight, missed poles/atoms и finite energy quadrature error; без дополнительных условий она не локализует причину.

Дополнительный дешёвый witness без Green evaluations: вычисляются boundary Λ(a),Λ(b) прямой свёрткой Γ с теми же exact real-hat coefficients (small lags BigFloat256, large lags convergent moments). При M=asinh(1)/π, B_L=−Λ(a)+MΓ_first, B_R=Λ(b)+MΓ_last полные45×45 margins dL=λmin(H45−aI−I9⊗B_L), dR=λmin(bI−H45−I9⊗B_R). Только при dL,dR>0 записываются βL=h||Γ_first||F/(4πdL²), βR=h||Γ_last||F/(4πdR²). Frobenius norm — upper bound operator norm; нет factor9. Это conditional floating ghost-mass witness, не formal certificate; d≤0 даёт unavailable. Endpoint matrices сохраняются вместе с moments.
