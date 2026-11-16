# Что переносить из NESSi/libcntr в QCL-NEGF

Практический алгоритм: сначала решать и независимо проверять равновесную
задачу с общими μ/T и causal spectral representation, затем переносить
отклонение от неё в Keldysh map; точность преобразований проверять отдельно
от nonlinear residual. Для QCL полезнее этот порядок и обработка граничных
скачков/хвостов, чем перенос двухвременного Kadanoff–Baym solver. Ниже
предлагаются три ограниченные гипотезы: **KMS-aware interpolation residual,
Matsubara density oracle без analytic continuation, затем разделение
transform error floor и SCBA convergence**. Это предложения, не запуски и не
новые критерии принятия текущего baseline.

Работа read-only: исходники/статьи прочитаны; NESSi не устанавливалась,
не компилировалась и не исполнялась. Julia/Pkg/новые расчёты не запускались.
Core и frozen PI harness не менялись. Quantica/libNEGF/Kwant здесь не изучались.

## Что именно проверено

Официальный репозиторий [nessi-cntr/nessi](https://github.com/nessi-cntr/nessi),
shallow checkout в `adjacent/nessi-source`, **5.4 MiB**, clean HEAD
[`6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6`](https://github.com/nessi-cntr/nessi/commit/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6),
commit date 2026-05-01, merge `version2`. Лицензия в
[LICENSE.md](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/LICENSE.md) — MPL-2.0.
Настоящая записка предлагает самостоятельно реализуемые формулы/проверки;
копирования исходного C++ в Core нет.

Первичные публикации: Schüler et al., **Comput. Phys. Commun. 257, 107484
(2020)**, [DOI](https://doi.org/10.1016/j.cpc.2020.107484),
[полный preprint 1911.01211v1](https://arxiv.org/html/1911.01211v1);
Künzel et al., **NESSi 2.0**, [2605.02606v1, 4 May 2026](https://arxiv.org/html/2605.02606v1).
Версия 2 добавляет memory truncation и stationary `ness2`; это существеннее
для нашей задачи, чем старое описание только real-time. Официальный README
ведёт на [актуальную документацию](https://nessi.readthedocs.io/en/latest/);
старый [manual v1](https://nessi.tuxfamily.org/PMan.html) остаётся полезен для
контурных routines. При расхождении ниже приоритет exact source HEAD.

### Реальный algorithm map

| Часть | Реальный путь на указанном HEAD | Что это означает для QCL |
|---|---|---|
| Равновесие | `cntr::dyson_mat`: Fourier, steep/CG, fixed-point; затем согласованный initial contour slice | Отдельный equilibrium solve и independent oracle, не однократная FD-подмена несошедшегося G |
| Real time | Matsubara → simultaneous bootstrap первых k steps → polynomial predictor → несколько nonlinear correctors каждого нового time slice | Causal marching упрощает задачу в **physical time**, но не превращает stationary SCBA iteration в Volterra problem |
| Stationary v2 | `ness2::dyson`: transform Σ(t)→Σ(ω), invert GR, G<=GRΣ<GA, transform spectral/lesser обратно | Прямой родственный frequency-domain Dyson; временная сетка/окно становятся дополнительными controls |
| Электрон–фонон | uMig использует fixed thermal D0, sMig также решает D=D0+D0ΠD | Fixed bath и hot/dressed phonons — разные модели, а не режим оптимизации |

На stationary source `dyson_fft` выполняет
`GR=[(ω+μ)I−H−ΣR−Σbath]^-1`, `G<=GR(Σ<+Σbath<)GA`.
Это [ness2_dyson_impl.hpp:61–91](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_dyson_impl.hpp#L61-L91).
В QCL energy reference/μ входят иначе: нельзя добавить μ ещё раз к native h.
Также NESSi scalar DOS `a(E)=−Im GR/π`, тогда как наш matrix
`A=i(GR−GA)=2πa`. NESSi initialization с DOS даёт G<=2πif a, то есть
наш `ifA`; множитель 2π должен учитываться ровно один раз.
[Source equilibrium initialization:156–194](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_equilibrium_decl.hpp#L156-L194).

## Численные механизмы, которые стоит понять до переноса

### Matsubara/KMS: система согласованных компонентов

В Holstein example сначала итерационно обновляются Σ[G,D0] и Matsubara G;
первый Fourier solve служит стартом, затем программа переключается на
fixed-point refinement. Только после этого задаются contour initial values
и запускается bootstrap.
[Holstein_bethe_uMig.cpp:195–266](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/Holstein_bethe_uMig.cpp#L195-L266),
[initialization:295–326](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/Holstein_bethe_uMig.cpp#L295-L326).
`set_t0_from_mat` связывает mixed, lesser и retarded equal-time components,
а не заполняет их независимо;
[cntr_utilities_impl.hpp:467–477](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/cntr_utilities_impl.hpp#L467-L477).

Внутри Matsubara Fourier Dyson выделяется известный `1/iω` tail:
вычисляется остаток `[iω−H−Σ]^-1−1/iω`, отдельно возвращается аналитическая
часть в τ. Это реальная обработка fermionic jump, а не density normalization.
[dyson_mat_fourier_dispatch:389–431](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/cntr_dyson_impl.hpp#L389-L431),
[set_first_order_tail:32–94](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/cntr_matsubara_impl.hpp#L32-L94).
Важно различать внутренний **linear Dyson** fixed-point solve для заданной Σ
и внешний nonlinear update Σ[G]. Слово `fixpoint` в API не означает, что
внешняя SCBA уже принята.

### High-order Volterra и predictor/corrector

`integration::Integrator` содержит Gregory endpoint weights, BDF coefficients,
полиномиальные startup integrals и отдельную `rcorr` для Matsubara convolutions,
которые нельзя гладко продолжать через 0/β.
[integration.hpp:41–104](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/integration.hpp#L41-L104).
Например, bulk Gregory weights равны единице, correction сосредоточена у концов;
BDF в комментариях имеет порядок k+1. Реальная глобальная настройка solver
`MAX_SOLVE_ORDER=5`, Gregory tables допускают больше порядков. Нельзя обещать
«весь solver seventh-order» по README.
[global_settings:38](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/cntr_global_settings.hpp#L38),
[integration.hpp:137–163](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/integration.hpp#L137-L163).

Retarded timestep собирает history terms и решает bounded matrix linear system
с BDF diagonal и Gregory endpoint coefficient;
[cntr_dyson_impl.hpp:181–213](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/cntr_dyson_impl.hpp#L181-L213).
Predictor экстраполирует **physical time slices** полиномом по предыдущим
slices; corrector снова строит Σ и решает Dyson на новом slice.
[extrapolate_timestep_dispatch:19–73](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/cntr_utilities_impl.hpp#L19-L73),
[Holstein corrector:403–449](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/Holstein_bethe_uMig.cpp#L403-L449).
Для QCL переносимы согласованные quadrature/startup/endpoint rules и
extrapolated initial guess вдоль **разрешённой physical parameter continuation**.
Высоко-порядковая экстраполяция SCBA iteration history сама по себе не
получает Volterra stability и может уничтожить PSD.

### Causal transforms stationary v2

NESSi различает plain DFT и Fourier **integral**. В retarded transform есть
нижняя граница t=0: trapezoid branch вычитает половину G(0), затем умножает
на h. При обратном преобразовании берётся GR−GA и оставляется положительное
time support. Есть cubic interpolation correction routines, но public
`ness2::dyson` на этом HEAD отклоняет всё кроме `FFT_TRAPEZ`.
[transform_to_freq:230–284](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_herm_matrix_ness_impl.hpp#L230-L284),
[transform_to_time:339–395](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_herm_matrix_ness_impl.hpp#L339-L395),
[dyson dispatch:11–22](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_dyson_impl.hpp#L11-L22).

Periodic FFT arrays не доказывают физическую periodic energy continuation.
`dyson_fft` ещё и использует внутренний frequency band приблизительно в три
раза уже полного FFT band (`nfreq=(Nft/2−1)/3`). Это часть numerical padding/
cutoff scheme, которую нельзя молча применить к source QCL window.
[тот же source:67–88](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_dyson_impl.hpp#L67-L88).

## SCBA/electron–phonon и conservation

В uMig `Sigma_uMig` реализует `Σ(t,t′)=i g(t)D0(t,t′)G(t,t′)g(t′)` через
`Bubble2`, с фиксированным D0;
[source:174–187](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/Holstein_impurity_impl.cpp#L174-L187).
В sMig вычисляется polarization `Π=−2i gGg G(reverse times)`, обновляется
D через Dyson/VIE2, затем Σ использует dressed D;
[Sigma_Mig:24–50](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/Holstein_impurity_impl.cpp#L24-L50),
[step_D:258–272](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/Holstein_impurity_impl.cpp#L258-L272).

При fixed thermal phonons электронная энергия может уходить в bath: это не
ошибка conservation. При fully coupled closed electron+phonon dynamics
нужно проверять полную энергию. Это различие прямо пояснено в официальных
[Holstein examples, uMig/sMig](https://nessi.tuxfamily.org/P4.html).
Общий первичный критерий — self-consistent symmetry/functional structure
приближения, а не имя «Born/Migdal»;
[Baym, Phys. Rev. 127, 1391 (1962)](https://doi.org/10.1103/PhysRev.127.1391).
Для QCL fixed LO bath переносится **particle balance + energy flow to bath**,
не требование нулевого электронного energy collision integral. Local Holstein
vertex не воспроизводит наши momentum-dependent matrix kernels, LO form
factors, finite-chain embedding, impurity/IFR/alloy channels.

### Mixing: хороший pattern и граница

Stationary `ness2_bethe` измеряет `distance_norm2(Gcandidate,Gold)` **до**
linear mixing, затем смешивает. В Matsubara Holstein example error измеряется
**после** mixing0.5. Поэтому первое годится как pattern raw defect, второе
нельзя копировать как достаточную acceptance measure. Сам `distance_norm2`
есть unweighted Frobenius массива ret+les; при смене сетки это не наша
weighted/occupied scientific norm.
[ness2_bethe.cpp:129–140](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/ness2_bethe.cpp#L129-L140),
[Holstein example:228–258](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/Holstein_bethe_uMig.cpp#L228-L258),
[norm source:10–35](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_utils_impl.hpp#L10-L35).
Наш raw rK/rΣ уже решает часть этой проблемы; повторять инфраструктуру не нужно.

## Три ранжированные гипотезы и различающие опыты

### 1. KMS-aware interpolation: интерполировать nonequilibrium residual

**Гипотеза, наша:** fractional energy shifts портят equilibrium manifold
через некоммутативность `Sδ[fA] ≠ f(E+δ)Sδ[A]`. Одного высокого порядка Sδ
недостаточно. Мотив NESSi — согласованные thermal D0/G и компонентная
контурная algebra; собственный перенос в stationary QCL:

`Gn(E)=f(E)A(E)+δGn(E)`, `Gp(E)=(1−f(E))A(E)+δGp(E)`.

Интерполировать A и residuals, а f вычислять аналитически в **shifted physical
energy**. При δGn=δGp=0, общих μ/T и Bose NΩ точные Fermi/Bose identities
восстанавливают LO KMS algebra для выбранного interpolated A, даже при
fractional Ω/ΔE. Это отдельный discretization variant, не доказанная exact
оптимизация нынешнего nodal map. Вне равновесия residuals остаются свободными;
никакая FD projection исходного nonequilibrium G не разрешается.

**Минимальный опыт:** positive Hermitian 2×2 Lorentzian mixture с известным A,
NE129/257/513, Ω/ΔE=8 и8.3, общий μ/T, прежняя window/zero-outside policy.
Сравнить KMS numerator/occupied scale, eigenvalue minima Σin/out, signed
particle collision integral и weighted-adjoint defect shift pair. Если KMS
улучшена, но particle balance/PSD ухудшились, transplant не принят. После
этого один frozen actual source slice, без nonlinear solve. Окно опыта:
≤30s/≤512MiB/BLAS1, только после отдельного ресурсного решения.

**Граница:** positive interpolation A не гарантирует positivity суммы с
произвольным residual. Нужна проверка paired Σ</> и conservation, не clipping.
Текущий Γ-thermal/AAII pilot ведущего не повторяется: предлагаемая добавка —
**KMS-aware shift algebra для отклонения от equilibrium**, особенно перед
bias continuation.

```julia
function shifted_pair(A, Gn, Gp, E, δ, μ, θ)  # θ=kB*T, same energy units as E
    f0 = fermi.(E, μ, θ)
    δn, δp = Gn .- f0 .* A, Gp .- (1 .- f0) .* A
    Es = E .+ δ                         # explicit physical shift, no wrap
    As = positive_nodal_interpolate(A, E, Es)
    ns = fermi.(Es, μ, θ) .* As .+ interpolate(δn, E, Es)
    ps = (1 .- fermi.(Es, μ, θ)) .* As .+ interpolate(δp, E, Es)
    return ns, ps                      # no projection/normalization
end
```

### 2. Off-axis Matsubara density oracle без analytic continuation

**Гипотеза, наша:** equilibrium density лучше независимо проверить off-axis
causal resolvent и аналитическим tail reference, чем ещё раз сравнить два
близких real-axis quadratures. Это опора для уже существующего PI control,
не повтор PI-vs-nodal reconstruction и не полный Matsubara port.

Для **одного и того же** causal finite-support Σ(z), source h/T/μ,
`G(iω)=G(z=μ+iω)` непосредственно вычислим: analytic continuation из шумных
Matsubara data не требуется. Для ортонормированного полного propagator
`ρ = I/2 +(2/β) Σn≥0 Herm G(μ+iωn)`. Ускорить сходимость можно вычитанием
free reference Gref с теми же leading high-frequency moments и добавлением
его точной `ρref=f(Href−μ)`. Это наша численная конструкция по мотиву tail
subtraction NESSi, а не найденная там готовая функция QCL.

**Минимальный опыт:** сначала exactly solvable two-level Schur-complement
model; exact density — block of matrix Fermi function полного H. Nω32/64/128,
analytic rational Σ, сравнить full matrices, не только trace. Такой независимый
embedded oracle уже демонстрирует NESSi
[test_equilibrium.cpp:48–83](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/test_equilibrium.cpp#L48-L83).
Затем **один** наиболее occupied source k: analytic off-real Cauchy transform
той же PI hat Γ и тот же9-cell chain; сравнивать central5 block density,
без factor9, fitting tails и μ refit.
Предлагаемый cap≤120s/≤1GiB, никакого запуска сейчас.

**Что различит:** Matsubara tail convergence при фиксированной модели против
finite-window ∫fA. Несовпадение может означать missing ghost strips/poles/dark
atoms либо real-axis resolution; без дополнительных bounds оно не определяет
единственную причину. Важная защита: если считать только finite-window
spectral measure с zeroth moment M≠I, её identity имеет M/2, **не I/2**.
Нельзя подставить I и объявить исчезнувшим spectral deficit.

```julia
function density_off_axis(Σz, H, β, μ, Nω)
    E, U = eigen(Hermitian(H))
    ρ = U * Diagonal(fermi.(E, μ, inv(β))) * U'
    for n in 0:Nω-1
        z = μ + im * (2n + 1) * π / β
        G, Gref = inv(z*I - H - Σz(z)), inv(z*I - H)
        R = G - Gref                    # matched analytic asymptote
        ρ += (R + R') / β
    end
    return ρ                            # repeat Nω; bound remaining tail
end
```

### 3. Сначала calibrate transform floor, затем ускорять nonlinear map

**Гипотеза:** часть apparent SCBA stagnation определяется несовместимыми
convolution/boundary discretizations, а не выбором mixing. В stationary
NESSi accuracy зависит отдельно от time-step h и long-time cutoff Nft·h;
frequency resolution и window width не один параметр. Перенос — двухосевая
operator calibration и independent exact block oracle, а не перевод QCL в
O(Nt²) two-time storage.

Уже полученный frozen PI control ведущего переиспользуем:
`pi-frozen-control/attempt01/summary.toml`, source e6ff…, Core2a500740…,
μfixed; он измеряет operator sensitivity, не nonlinear convergence.
Новый минимальный опыт — на rational scalar/complex2×2 oracle отдельно менять
quadrature order, domain и endpoint/tail treatment; затем один frozen actual
operator application с двумя **обоснованными** accuracies. Сопоставить δmap
с raw nonlinear defect и δN/target. Если δmap больше нужного residual,
добавочные correctors/AAII не принимают исходное дискретное уравнение;
если δmap мала, исследовать conditioning/Jacobian map, сохраняя исходный
raw residual. Никакого автоматического ослабления tolerance.

Существующие PI/analytic tests и checkpoint checks использовать повторно.
Достаточно небольшого диагностического patch или Julia prototype; новая
abstraction/platform/third-party build не нужны. Не обещаем, что high-order
weights сохраняют PSD, KMS или current conservation — это разные проверки.

## Что нельзя принимать на доверии даже в соседнем solver

1. **Static sign discrepancy в `ness2::force_equilibrium`** на указанном HEAD.
   `FERMION=-1`; routine задаёт `ImG=(GR−GR†)/2`, затем
   `Gless=−(2*bosefermi)f ImG`. Для fermion GR=−i получается −2if, тогда как
   свой `green_equilibrium_ness_fermi` даёт +2if при a=1/π. Это точный
   scalar sign counterexample по тексту source, не результат исполнения
   NESSi и не доказательство ошибки всех её equilibrium paths.
   [force_equilibrium:448–480](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_herm_matrix_ness_impl.hpp#L448-L480),
   [FERMION:57](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/cntr_global_settings.hpp#L57),
   [independent initializer:184–191](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_equilibrium_decl.hpp#L184-L191).
   Следствие: заимствовать идею component consistency, проверять знаки своим
   single-level/KMS oracle. Не использовать эту routine как независимый
   oracle для QCL без отдельного исправления и runtime verification.
2. Cubic integral transforms не равны supported cubic Dyson; optional bath
   η — реальный дополнительный reservoir с lesser injection, а не невинный
   solver stabilizer. В QCL перенос η менял бы physical identity.
3. Upstream tests — полезные примеры контрактов, не полученные нами passes:
   [ness2_dyson.cpp:22–59](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/test/ness2_dyson.cpp#L22-L59)
   сравнивает stationary с contour solution; он не проверяет наш LO kernel,
   energy window или весь nonuniform solver.

## Решение для кампании

Сохранить native PI и текущий equilibrium pilot как самостоятельные варианты.
Первый следующий переносимый вопрос — KMS-aware residual interpolation с
одновременной conservation/PSD проверкой. Второй — independent off-axis
frozen density oracle. Полный real-time NESSi port, memory-truncated KBE,
phonon dressing и новый bath сейчас не обоснованы нашим stationary вопросом.
Ни source review, ни согласие агентов не принимают QCL SCBA, grid/window
adequacy или сравнение с экспериментом. Все предложенные numerical outcomes
пока `not_measured`.
