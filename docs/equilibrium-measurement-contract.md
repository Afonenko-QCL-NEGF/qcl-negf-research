Публикуемый source contract SCI03: **CP-DESIGN APPROVED только для bounded
measurement/pilot design** по независимому R2 от 2026-10-19. Review receipt SHA256
`351b58065dfd4dd45bbacf95a451a6daa22adaeacaee09149457f8e0ba295442`.
Numerical/scientific acceptance и CP-SCI runtime READY не установлены: finite
accuracy εQ, tail/truncation/interaction и nonlinear-conditioning bounds остаются
неизвестными; exact frozen plans, enforceable guards и fresh admission ещё нужны.

Технический корпус ниже — точные 40 812 bytes reviewed D2, SHA256
`788ab058f1c52520d79aa621dcecd0d525f2562b7fcf2ea10da2c664989c8018`.
Его исходный заголовок SCI03-D1, source snapshot и historical candidate markers
сохранены как provenance; актуальный статус CP-DESIGN указан выше. Четырнадцать
owning file hashes в §11 повторно сверены при переносе. Этот документ не меняет
formula, criteria, model, parameters или plans и не разрешает запуск.
Пути/aliases раскрыты после корпуса; exact body можно извлечь между markers,
без исходного чата или обязательных `/tmp` inputs.

<!-- BEGIN SCI03_REVIEWED_D2_BODY -->
# SCI03-D1 — научный бюджет ошибки и контракт равновесной проверки

2026-10-19. Автор: D / Sol ultra. **Кандидат для независимого CP-DESIGN**, не
принятый scientific contract и не CP-SCI READY. Этап — source/literature design;
SCBA/Poisson, Julia, продуктовые imports/tests/builds и новые scientific data не
выполнялись. Разрешение EXEC/B/G1/S1 без e–e передано координатором; PLAN-маркеры
исторического handoff сами по себе не отменяют это позднейшее поручение.

## 1. Вопрос, неизменяемые inputs и значение результата

Вопрос: допускает ли объявленная конечная модель stationary field-periodic
embedding самосогласованный thermal fixed point при **0 mV/period, lattice/LO
70 K**, с одной μ, правильным raw U/G, нейтральностью и нулевым интегральным
током; какие ограничения дискретизации остаются у этого вывода?

Сохраняются LO/acoustic/impurity/IFR/alloy, **electron–electron=none**, thermal
LO, fixed screening, parabolic dispersion, finite_chain embedding4. Точный
declared input и composition задают Research sources и gitlinks; SCI03 не меняет
структуру, dose, boundary, meshes, mixing, tolerances или scientific target.
Acoustic сейчас elastic; утверждение об уникальном возврате к FD из любого
initial не входит в проверку. Проверяется thermal fixed point и ограниченная
phase robustness, не глобальная эргодичность.

Snapshot чтения: root `7d179b28647f0b89347d2d6358f08cf2fe479653`;
Core `f5b73369a96d1225bd96fb34f5262c8dcf965f86`;
Research `b1a3be2b0f4b3adf423cc2a611bcca2b5e0a5ff3`;
Results `17574ca4c5495f74d86ed9bef51c60dea8265062`.
Это evidence locators, не второй maintained revision lock. Core/Research/Results
были чистыми; root имел concurrent dirty gitlinks/coordination artifacts,
перечисленные в report. Они не редактировались автором. Runner current dirty
source не принимается этим design за reviewed persistence implementation.

Прочитаны AGENTS, architecture, весь README/autonomy/operator-permissions/science
handoff, qcl-physics-audit и physics-checks, milestone, SCI02 owning passport.
В этих трёх components вложенных AGENTS.md не найдено. Старые S01 arrays,
защищённый checkout, run и archived attempts не читались и не менялись.

`run=not_created`; runtime/imported-producer, frozen plan bytes/hash/fingerprints,
execution/point/attempt/state IDs = `null`, причина `source-only design`.
Из source ни один численный итог не выводится как measured result.

## 2. Проверенные первичные формулы и границы attribution

| ID | Первичный источник, проверенный locator | Identity/hash и применимость |
| --- | --- | --- |
| P-FD | [Smirnov, Distribution of Fermions](https://www.iue.tuwien.ac.at/phd/smirnov/node25.html), (2.30)–(2.31): FD и определение общей μ числом частиц | HTML 12 970 bytes; SHA256 `8f47147bcdbe62fe340129d07e7638b4b52abf71141994872350efe49092ab76`, прочитан 2026-10-19. Никакой конечный tolerance из формулы не следует |
| P-FDR | [Fotso–Freericks, 2020](https://www.acsu.buffalo.edu/~hffotso/MyPapers/FrontiersPhys_FDT.pdf), DOI10.3389/fphy.2020.00324, printed p.8, (42): g< = −2if Im gR | PDF 1 195 786 bytes; SHA256 `6e97b1668f2d5a7e8b05a92c8fe51e8e7fbc3d0087d7dff9ddbbbf8dc75dfa8f`; p.8 извлечена непосредственно из этих bytes. Источник пишет scalar/local relation; matrix form ниже — явно выведенное применение equilibrium identity к единой матрице, не утверждение о конкретном QCL |
| P-DOS | [Lee–Wacker, cond-mat/0212059v1](https://arxiv.org/pdf/cond-mat/0212059v1), Appendix B, printed p.17, (B1)–(B2), и §II.A, printed p.4 | Уже существующий PDF, SHA256 `dfabd791bfb2aa4603233e9b545879e07a3bfc63e8906cbc810760726358ea0c`; адресно прочитаны. B2 даёт 2D DOS ρ₀=m/(πħ²), а не готовый broadened-QCL occupation oracle. §II.A обсуждает малые **inelastic** acoustic transfers, которых elastic closure не воспроизводит |
| P-J | [Wacker, 0706.4202v1](https://arxiv.org/pdf/0706.4202v1), printed p.2, (2)–(6), Fig.1 discussion | Уже существующий PDF SHA256 `091a39824c3fb4b5f1f3b0f0659d4a99663b48b6949806aa8586710b6d1d937b`; p.1–2 прочитаны. Полный coherence matrix нужен току; local-current variation зависит от конечного basis. Его signed electron charge convention не переносится молча на outward Results convention |
| P-MODEL | [Bosco et al., 1911.06582v1](https://arxiv.org/pdf/1911.06582v1), DOI10.1063/1.5110305, p.3 main text и p.4 Fig.3 | PDF SHA256 `b103b3da9fa977d1dc682d1e2caf96d74b5f081731e4db0da944af675a58dc70`, SCI02 primary receipt и опубликованный паспорт. Main-text cm⁻² и caption cm⁻³ конфликтуют; declared sheet4.5e14 m⁻² сохраняется с caveat. Расчёт200 K/56 mV и digitized NEGF curve не oracle для70 K/0 mV |

Последние два PDF повторно не загружались. Для FD/FDR bounded HTTP чтение шло
только в память с SHA256; новых PDF files не сохранено. Первое sandbox чтение
не разрешило DNS; исправленное permitted network чтение FD прошло, FDR упёрлось
в initial1 MiB reader cap, затем один адресный 4 MiB read получил полный PDF.
Это разные установленные причины чтения, не solver attempts и не установленные
окружения. Baym не используется как доказательство conserving implementation:
слово SCBA само по себе такого доказательства не даёт.

В ортонормированном общем basis, при общем μ/T и отсутствии drive:

```math
f(E;\mu,T)=[1+\exp((E-\mu)/(k_BT))]^{-1},\quad
A=i(G^R-G^A),\quad N=-iG^<,\quad P=iG^>,
\qquad N=fA,\quad P=(1-f)A.
```

Greater relation выводится вместе с `G>−G<=GR−GA`; matrix relation относится ко
всем a,b, а не только trace/diagonal. Для идеальных независимых parabolic 2D
subbands интегрирование DOS даёт
`n_i=g_i m_i kBT/(2πħ²) log(1+exp((μ−E_i)/(kBT)))` [m⁻²]. Это наш интеграл
P-DOS+P-FD, с объявленным degeneracy без двойного spin. Он годится для analytic
anchor, а не для принятия broadened/coherent localized populations full model.

## 3. Единая мера, наблюдаемые и source binding

Пути C и V ниже относительно owning Core и Results. Полные hashes в §11.
Core `grids.jl::build_grids` использует dimensionless E и k: energy trapezoid
`wE`, endpoint half weights; `wk_m=(κ_{m+1/2}²−κ_{m−1/2}²)/(4π)` уже включает
radial `1/(2π)`; z cell weights `wZ=Lp/Nz`. Определяем **однократно**

```math
W_{em}=\frac{g_s}{2\pi}wE_e\,wk_m\ge0,
\quad \|X\|_W^2=\sum_{em}W_{em}\|X_{em}\|_F^2,
\quad M_W=\sum_{em}W_{em}.
```

Дополнительный spin или `2π` не вводится. Все physical conversions берут
проверенные persisted scales, не guessed constants: number `/L0²` [m⁻²],
density `/L0³` [m⁻³], Hartree `×E0_eV` [eV], current `×J0` [A/m²].
`J0=e E0/(ħ_eV L0²)` в принятой положительной scale convention.

| Наблюдаемая Q | Определение/норма/единицы | Near-zero, domain и источник |
| --- | --- | --- |
| μ | Одна μ от `Σ W ReTr(A) f = NDsheet*L0²`; собственные fit μ для каждого уровня/raw/normalized запрещены | μ-fit — construction, не pass. Нужны finite positive T, target>0, неотрицательные spectral-number weights, capacity>target, bracket. Сохранять endpoint μ и uncertainty; `V::_fit_mu` |
| FDR raw/normalized | `η²=(‖N−fA‖W²+‖P−(1−f)A‖W²)/‖A‖W²`; full complex matrices, обе версии с той же μ | No floor. `_ratio(0,0)` — локальное helper convention measured0 с причиной, не полный FDR verdict. Весь reader сначала сохраняет applicability gate и требует defined μ; unavailable μ даёт inherited μ status/reason/value=None. При applicable context, неотрицательных spectral-number weights и capacity≤0 результат `not_measured`, reason=`empty_represented_spectrum`; helper zero/zero его не заменяет. Внутри helper denominator0/error>0 — undefined. V::_fit_mu(239–244),diagnose_stored_equilibrium(459–463,528–531),_Norm,_ratio; P-FDR |
| Sheet numbers и normalization | Raw stored N/P, normalized N/P, candidate raw N; occupied donor residual `∣N−ND∣/ND`; λ=Nraw/ND | ND>0; no floor. Empty number measured, donor residual `null/not_applicable`. Коррекции `‖Nnorm−Nraw‖W/‖Nraw‖W`, отдельно P; абсолютные numerators сохранять. V/C final_audit |
| n(z) | Все coherence terms: `L0^-3 ΣW Re[χ(z)^T N χ(z)*]`; `‖Δn‖L2(dz)/(NDsheet/√Lp)` | Не сравнивать только диагонали/localized labels. Один physical z domain, conservative transfer и отдельно transfer bound; no floor; C::observables._electron_density_bar |
| U(z) | Gauge-aligned `max_z∣ΔU∣` [eV], mean-zero gauge по одинаковому measure | Gauge offset не physical error; missing gauge/mapping — incomparable. Не выдавать small mixed step за raw map rU |
| J± | Отдельные raw/normalized outward `J0 ΣW Tr(Σ>G<−Σ<G>)`; `max(∣J+∣,∣J−∣)`, `∣J++J−∣`, мнимые части [A/m²] | Absolute errors; division/J=0 запрещено. Minus common-z=−outward; balanced nonzero outward fluxes не thermal J=0. V report и C::_boundary_flux_complex |
| Spectrum | На каждом k `M(k)=ΣwE A/(2π)`, `max_k‖M−I‖F/√Nb`, trace и boundary/occupied tail measures отдельно | No radial averaging этого worst-k gate. Tail loss и quadrature могут компенсироваться; neither edge ratio nor finite-window deficit is outside-tail bound |
| Collisions/powers | Подписанный integral channel до abs, каждый LO/acoustic/impurity/IFR/alloy отдельно; powers [W/m²] | Не abs каждого integrand для numerator. Поле0 даёт field work0; finite boundary энергетика должна соответствовать owning periodic model. Сумма powers не исключает ошибки компенсации channels |

Raw V diagnostics восстанавливают Keldysh с **stored** scattering+embedding_total
ровно один раз. Plus/minus используются для потоков; суммировать total+plus+minus
в raw нельзя. V не строит fresh map и не выдаёт acceptance. C fresh map включает
все пять channels и embedding total/plus/minus; source validity — SCI04/O08-03.

## 4. Existing strict policy — полностью сохраняется

Числа ниже — **project acceptance policy**, не выведенные литературные errors.
Они присутствуют в model YAML `solver.tolerances`, совпадают с Core typed defaults;
в прочитанных policy/study нет overrides этих tolerances. До submit это лишь
source-derived expected значения: frozen plan/runtime должны подтвердить каждый.
DiagnosticQuality/coarse thresholds никогда не заменяют strict final values.
Final evaluator принимает `value<=threshold`, inner/outer confirmations — `<`.

| Метрики | Actual norm / floor / units | Strict target и область |
| --- | --- | --- |
| rD stored | `max‖D GR−I‖F/(‖D‖F‖GR‖F+√Nb)`, block max, no quadrature | 1e−11; D с published U и storedΣ; backward, не forward error |
| rA stored | `‖A−GRΓGA‖F/(‖A‖F+1e−14)`, global unweighted | 1e−9; A scaled; no W |
| rD_candidate,rA_candidate | Те же нормы, freshΣ, тот же published U/G | **rK=1e−8**, не rD1e−11; нельзя смешивать budgets |
| rK | `‖G<−GRΣcandidate<GA‖F/(‖G<‖F+1e−14)` unweighted, candidate без fittedλ | 1e−8; fresh unmixed map |
| rΣ | Max family/component `‖Σnew−Σstored‖F/(‖Σnew‖F+1e−14)` | 1e−8; all enabled +total/plus/minus embedding; no α |
| rλ,r_lambda_candidate | `∣Nraw/ND−1∣`, ND positive, no floor, W | 1e−8; stored/raw и candidate отдельно |
| rP,rζ | Poisson `‖LUmap−b‖2/(‖b‖2+λP·1e−14)`; `∣ζ∣` scaled bordered system | 1e−11/1e−10; neutral periodic source, harmonic face ε, unweighted gauge |
| rU | `max∣Umap−Upublished∣`, scaled energy, no floor | 1e−7; same density/Green state, до mixing; **indicator×E0**, не solution bound |
| rn | Max last-outer/independent-final weighted-z density change divided ND/√Lp, no floor | 1e−7; outer history обязательно |
| rneutral | `∣ΣwZ(ND−n)∣/NDsheet`, no floor | 1e−10; stored density; normalized charge alone не raw certificate |
| rJchange(_scba) | `abs(J−Jprev)/max(absJ,absJprev,1e−12 A/m²)` | r_obs1e−4; first comparison unavailable; это temporal stability |
| rpopulation(_scba) | `‖p−pprev‖2/(‖p‖2+1e−14)`, p=diag(sheet matrix)/ND | r_obs1e−4; labels only within same basis, not FD |
| rPSD | `max(0,−λmin(HermX)−b)/max(‖HermX‖2,floatmin(Float64))`, `b=64 Nb eps max(‖HermX‖2,1)` | 1e−10; A,N/P raw+normalized,Γ; Herm only measures, never repairs |
| rcaus / identities / antihermiticity | Block spectral norm for ImΣ; identity numerator ΣR−ΣA−(Σ>−Σ<), denominator sum F norms+1e−14; lesser/greater antihermiticity F norm+1e−14 | 1e−10; scattering,embedding,total independently; no cancellation between families |
| rsum | `max_k‖ΣwE A/(2π)−I‖F/√Nb`, no floor | 1e−3; orthonormal finite basis, full-axis target |
| rJ | `∣Re(Φ++Φ−)∣/max(∣ReΦ+∣,∣ReΦ−∣,1e−12)` | 1e−6; floor в **scaled number flux**, physical floor=J0·1e−12 A/m²; не temporal floor и не absolute J gate |
| rC_channel | `c∣Re signedIntegral∣/(c activityIntegral+1e−14)`, c=gs/(2π), W without duplicate c | 1e−8; five named channels; near zero additive scaled rate floor |
| rpower | `∣ΣPch+JVp∣/(Σ∣Pch∣+∣JVp∣+1e−14P0)` | 1e−4; field-periodic energy balance; extra contact terms для иной model не придумываются |
| rImJ,rImC | Pointwise q=max∣Im trace∣, s=max∣Re trace∣; s>1e−12 ? q/s :100q | 1e−10; unweighted E,k maxima; branch switch scaled trace |
| rImP | `100 maxch∣Im Pch/P0∣`, no floor | 1e−10; absolute scaled power indicator |
| rroundJ/rroundC/rroundP | `∣F∣>1e−12 ? ∣F−Freverse∣/∣F∣ :100∣F−Freverse∣` | 1e−10; FJ=Φ±, prefactored scaled boundary flux с gs/(2π), max по sides; FC=signed из `_collision_sum`, **без** gs/(2π), отдельно по channel; FP=`_collision_power_bar`, prefactored dimensionless channel power с gs/(2π), max по internal channels. Switch1e−12 и absolute100·error применяются к этим разным operands. Empirical accumulation-order sensitivity, не rigorous error |
| rroundoff | FFT/direct `maxch‖Λ1−Λ2‖F/(‖Λ1‖F+1e−14)`; product integration uses max(norm,eps), отдельно | 1e−10; verify_fft_roundoff required; unmeasured не pass |
| rtail low/high/k | `∣Σselected p∣/Σall p`, p=wE wk ReTr N; energy fractions.05, momentum.1 | 1e−6; no floor; total≤0/nonfinite→Inf; empirical finite-domain indicator |
| redgeγ,redgeA | Endpoint maxΓsc /global maxΓsc; endpoint max∣ReTrA∣/global max∣ReTrA∣ | 1e−6; zero Γ denom evaluator returns0 (при enabled channels retain caveat), zero A denom→Inf; не tail bound |

Все нормы/branches — C `acceptance_metadata.jl`, `validation.jl`, `poisson.jl`,
`observables.jl`; registry — **описание единственного owning evaluator**, не
предложение повторить scientific gates в Results/UI. Metadata
`error_budget.estimated_bound=null, bound_status=not_measured` сохранены буквально.
e–e/hot-LO extra gates `not_applicable` по объявленному none/thermal с основанием;
отсутствие channel/mode evidence блокирует соответствующую applicability.

Непосредственные conversions, **не новые tolerances**: E0=.1 eV даёт rU strict
map indicator1e−8 eV; ND=4.5e14 m⁻² даёт raw-number indicator4.5e6 m⁻² при rλ1e−8
и stored-neutrality indicator4.5e4 m⁻² при1e−10. Nb5 PSD allowance coefficient
64Nb eps=7.105427357601002e−14, Nb7 —9.947598300641403e−14 (умножается на actual
max(norm,1)). Ни одно из этих чисел не объявляется accuracy всей модели.

## 5. Выводимые bounds, которые действительно можно использовать

### 5.1 μ algorithm uncertainty

При positive spectral weights `q_e=Σ_m W ReTr A`, positive T и bracket внутри
capacity, number strictly monotonic; одна μ разрешается этим scalar equation.
`dN/dμ=Σ_e q_e f_e(1−f_e)/(kBT)`; нужна lower derivative bound на **всём**
bracket, если из number error выводится μ error. Близость target к0/capacity
может сделать эту обратную чувствительность очень большой; guessed floor нельзя.

Из P-FD `|df/dμ|≤1/(4kBT)`. Для неподвижного A и bracket midpoint uncertainty
`uμ=(μhigh−μlow)/2` triangle inequality даёт

```math
|\eta(\widehat\mu)-\eta(\mu)|
\le \frac{\sqrt2}{4k_BT}\,u_\mu.
```

При70 K `kBT=(1.380649e−23/1.602176634e−19)·70≈.0060321332835 eV`,
коэффициент≈58.61166754 eV⁻¹. Это bound только ошибки μ **construction** при
фиксированном spectrum; physical μ error от изменившегося A/window/basis остаётся
отдельным. `8 eps max(Escale,kBT)` default у V — machine stopping rule, не scientific
accuracy. Nonbracketed/algorithm_limit/negative q не становятся FD pass.

### 5.2 Linear observables в фиксированном представленном basis

Для any finite ΔN и положительных W, Cauchy–Schwarz и trace inequality дают

```math
|\Delta N^{2D}|\le L_0^{-2}\sqrt{N_b M_W}\,\|\Delta N\|_W,
\quad
|\Delta n(z)|\le L_0^{-3}\|\chi(z)\|_2^2\sqrt{M_W}\,\|\Delta N\|_W.
```

Второй bound относится к actual complex χ с coherence terms; это не сравнение
диагональных populations. Он позволяет выразить correction/FDR defect в m⁻² или
m⁻³ для одной дискретизации. Ошибки basis/χ/weights не охвачены. Для перевода
unweighted residual numerator в W norm можно использовать `√max(W)`; надежнее
прямо измерять оба absolute numerator в одном frozen final audit. Относительное
rK само по себе не знает масштабов этих coefficients.

### 5.3 Абсолютный represented current без division/J

Для каждой boundary side отдельно, в dimensionless stored matrices определим
`Sn=−iΣ<`, `Sp=iΣ>`, `Γ=Sn+Sp`, `δS=Sn−fΓ`, `δN=N−fA`, `δP=P−(1−f)A`.
Algebraically, без commutation assumptions,

```math
\operatorname{Tr}(S_pN-S_nP)=
\operatorname{Tr}[(1-f)\Gamma\delta N-f\Gamma\delta P-
\delta S(A+\delta N+\delta P)].
```

По Frobenius/Cauchy:

```math
|J_s|\le B_{J,s}^{repr}=J_0\left[
\|(1-f)\Gamma_s\|_W\|\delta N\|_W+
\|f\Gamma_s\|_W\|\delta P\|_W+
\|\delta S_s\|_W(\|A\|_W+\|\delta N+\delta P\|_W)\right].
```

Это measurable upper bound **представленного exact-arithmetic trace** [A/m²];
применим отдельно raw stored, normalized, а при наличии SCI04 payload — fresh
candidate. Он не использует PSD/commutation и не требует cancellation потоков.
Σ FDR defect измеряется, а не предполагается0. Zero factors→absolute0; нет
защитного scientific denominator. Для complex J bound также ограничивает
мнимую часть. Неравенство не принимает произвольный большой B: конечный target
εJ_abs должен быть обоснован независимо; округление/неохваченные domain и
fixed-point errors добавляются. Existing reader ещё не вычисляет эти coefficient
norms; предлагаемое измерение требует owning Results reducer ticket с literal
2×2 oracle после CP-DESIGN, без переноса критериев из Core.

### 5.4 Алгебраическая и нелинейная ошибки различны

Для invertible block D, `G−D⁻¹=D⁻¹(DG−I)`. Backward rD даёт forward correction
только после upper bound `||D⁻¹||`; Hermitian eigenvalue gap нельзя использовать
вместо smallest singular value не-нормальной retarded matrix.

Для общего unmixed coupled map F, при доказанном contraction q<1 на области
с тем же fixed point: `||x−x*||≤||F(x)−x||/(1−q)`, затем для observable Lipschitz
LQ — `εiter,Q≤LQ||F(x)−x||/(1−q)`. Domain, weighted norm, map components и
local/nonlinear remainder должны быть объявлены. Историческое отношение двух
невязок или малая α не доказывает q. Сейчас q/LQ/general inverse-Jacobian
bound **не установлены**. Полные raw residuals остаются измерениями уравнений,
а не заявленной границей forward error.

## 6. Target accuracy и аддитивный бюджет

Ideal targets из equilibrium equations: ηraw=ηnormalized=0; each J±=0;
neutrality/raw donor defect=0; applicable algebraic/conservation identities=0;
spectral full-axis matrix=I. μ, n(z), U(z) имеют неизвестные exact solution
values. Поле0 не требует flat U или equal populations.

**Конечные numerical target accuracies** εη, εμ[eV], εn[m⁻³ или заданный weighted
relative norm], εsheet[m⁻²], εU[eV], εJ_abs[A/m²] в milestone/primary papers
**не заданы и не следуют из identities**. Поэтому они хранятся `null + reason`,
а не1e−4/seed linewidth/доля J0. Current strict r_obs1e−4 и rpower1e−4 сохраняются
в своих определениях; это не основание назначить η1e−4 или J/J0=1e−4.
Без требуемого различимого physical scale finite desired accuracy нельзя вывести
единственным образом. Design не скрывает это scientific input gap и не объявляет
физическую приёмку технически готовой. Вопрос не предлагается решить голосованием.

Для каждой Q и фиксированного target εQ контракт использует

```math
B_Q=B_{iter,Q}+B_{E,Q}+B_{window,Q}+B_{k,Q}+
B_{basis/z,Q}+B_{quadrature,Q}+B_{roundoff,Q}\le\epsilon_Q.
```

Coefficient каждого члена **1 по triangle inequality**, не1/7 и не fitted
allocation. Сохраняются absolute units/норма каждой Q, component upper bound,
source/assumptions, validity domain, missing reason. Такой учёт не требует
произвольно раздать fractional budgets. При нескольких источниках одного bias
нельзя одновременно списать один defect в разные независимые rows. Для coupled
cutoff/window interactions нужен отдельный explicit remainder; нет доказанного
remainder→B_Q unknown. Root-sum-squares запрещено без доказанной random independence.
Model/value-origin uncertainty (dose/ionization/material assumptions) и
experimental uncertainty сюда не включаются как численный шум.

| Член | Что уже определяется | Неизвестный coefficient / остаток | Измерение в принятом S1 |
| --- | --- | --- | --- |
| Biter | Raw strict residual norms, same final U/G; linear bounds §5 | Nonlinear q/inverse Jacobian/domain/LQ, actual rD conditioning; finite stopping defect | B0 fresh-map/history/Poisson identity; SCI04/O08-03 source contract. R только по доказанной cause, не новый algorithm sweep |
| BE | Actual h, LO fractional shift; E1−B0 and P1−B0 physical Q differences | Global quadrature smoothness constants/peak widths/remainder; Richardson order **не** из двух h | E1 NE61441 same window; P1 half-cell shift. Differences — sensitivities до valid bound |
| Bwindow | rsum/edge/occupied boundary fractions, common-domain difference | Outside-window spectrum/collision/Hilbert tail majorants и quadrature-tail separation | W window[−.45,1.05] NE38401 same h; E1/P1 помогают quadrature. Один wider window не даёт infinite-tail bound |
| Bk | K1−B0 on same cutoff, nonnegative annular weights | Radial quadrature derivative/remainder, unmeasured k>kmax contribution | K1 Nk21/kmax.60; C Nk16/kmax.9 +energy_max1.2 coupled, не pure cutoff |
| Bbasis/z | n/U common z transfer; physical projectors/current matrix comparisons | Missing-subspace/projector defect, basis-period3 truncation; χ/geometry voxelization change under Nz | Z Nz192; B Nb7. Period-count refinement/embedding4→8 не добавляются; при отсутствии bound этот contribution unknown |
| Bquadrature | Each angular/qz internal-channel difference | Remainders angular/qz; qz cutoff10 unchanged, no outside-qz proof | A angular64, Q qz161/max10; hypotheses freeze before execution |
| Broundoff | Reverse-order/Hilbert padding indicators; PSD construction allowance | Forward arithmetic error interval/condition constants; padding sensitivity не rigorous upper bound | B0 existing verify_fft_roundoff + independent source/oracle evidence, no new benchmark |

Для BE source arithmetic: B0 h39.0625 μeV, h/kBT≈.00647573556;
P1 shift19.53125 μeV; LO36.7 meV/h=939.52, E1 —1879.04. Интерполяционный KMS
defect поэтому проверяется, а не assumed0. h≪kBT не доказывает разрешение узких
spectral resonances; requested seed1.25 meV не конечная physical linewidth.

## 7. Сравнение variants: одно состояние, mapping, domain

Сравниваются только individually strict stationary states той же physics/model
с independent raw-final audit. Совпадение двух несошедшихся candidates ничего не
сертифицирует. Для μ выровнять common E_ref/gauge; для U вычесть mean с объявленным
z measure; n перенести conservatively на единый physical domain и измерить transfer
error. Basis changes: сравнивать physical density/projectors/coherence observables
с cross-basis map, не state label. Common energy/k integration отдельно от
добавленных/потерянных tails. Transfer/interpolation не прячет coarsening aliasing.

S1 C и W имеют разные окна: C[−.30,1.20], W[−.45,1.05] eV. Их общий domain
[−.30,1.05] eV и signed extra tails сохранять раздельно. Простое C−W не изолирует
k cutoff: остаются две различные boundary strips и self-consistent response.
Ни C−B0, ни вычитание W−B0 без interaction bound не дают чистого Bk.
Одновременный выбор «лучших» E1/K1/Z/B/A/Q не создаёт автоматически разрешённый
combined final config.

## 8. Pilot, freeze и заранее held-out H

**P0 measurement contract**: frozen exact source/model/plan identities,
все strict policy values §4, observable definitions §3, additive accounting §6,
raw/normalized/stored-vs-fresh distinctions, all known unknowns, exact output
dataset/weights/scales and H protocol. Freeze после independent CP-DESIGN и до
первого submit. Его content hash хранит coordinator; не source-selected число.

1. При missing εQ или remainder bounds разрешён только clearly labelled finite
   **pilot** из уже принятого S1. Сначала B0, max_attempts=1. Ни unknown field, ни
   normalization/seed/exit0 не дают pass. Fresh resource/admission/guards остаются
   обязательными; сегодня они missing, поэтому launch здесь не READY.
2. B0 должен закрыть relevant strict coupled convergence/final-state consistency,
   history confirmations и all applicable existing strict physical gates до
   routine accuracy refinements. Failure передаётся SCI13 cause review; это не
   разрешение слепо запускать E1/P1/K1 «для accuracy».
3. E1/P1/K1, далее условные W/Z/B/A/Q/C сохраняют измеренные differences и honest
   `bound_status=not_measured` при неизвестном remainder. Все13 slots,72h/attempt,
   936 execution node-hours, concurrency≤2 и P1 byte ceilings общие, без расширения.
4. **A1 acceptance version** может заполнить target/bounds только по независимому
   scientific scale/source/derivation. Если pilot помог выбрать coefficients или
   criteria, новые criteria не дают этому pilot retroactive acceptance. Его raw
   данные остаются exploratory evidence, а training-like dependencies объявляются.
5. Перед H freeze exact A1 bytes/hash, выбранный **уже разрешённый single variant**
   final candidate (или B0), source/plan и every required proof/control. Нет подходящего
   candidate/bound/target→H не запускается для подбора tolerances; gap сохраняется.

**Предложенный точный held-out protocol, объявленный до pilot:** у selected final
candidate оставить physics/basis/Nz/k/quadratures/ΔE/window width, метод и strict
policy; сдвинуть оба endpoint energies на **+ΔE/4**, где ΔE вычисляется из его
frozen pre-H plan. Создать новый H plan/identity, independent run с обычным FD seed
по новой сетке; seed не evidence. Это заранее выбранная quarter-cell phase,
отличная от P1 half-cell; знак/доля не выбираются по результату. Source assignment
проверяет SCI09/12: H является разрешённым held-out phase protocol, не скрытой
новой physics. Не добавлять одновременно initial-amplitude/model/method variants.
H raw не читается/не используется для выбора εQ/coefficient/candidate до A1 freeze.
После H только independent result verdict; failure сохраняется, criterion edit
не принимает тот же H. Если единственный H уже consumed, нового held-out solve
автоматически не появляется. H подтверждает заявленную phase robustness и final
contract на отдельном state, но сам не доказывает missing continuum tails/q.

Незатронутые analytic/operator controls остаются отдельными независимыми controls.
Общий production kernel в fresh map — useful audit, не независимый physics oracle.

## 9. Пять уровней verdict и правила ошибки

| Уровень | Необходимое evidence | Сейчас |
| --- | --- | --- |
| Execution | Frozen plan/source/runtime, run/attempt/state, process+durability+complete raw receipt | not_measured; run not_created |
| Equations | Unmixed fresh rDcandidate/rAcandidate/rK/rΣ/λ, aligned final U/G/reconstructed n, Poisson/gauge/history confirmations strict | not_measured |
| Physical checks | Applicable strict gates + full-matrix raw/normalized FD/FDR + absolute J± with justified ε/bounds; identities/PSD/balances/number corrections | not_measured; finite ε/bounds unresolved |
| Discretization | Individually converged comparable states + valid additive bounds in every claim-relevant axis, including tails/truncation/interactions | not_measured; finite catalog не гарантирует этот уровень |
| External validation | Independent comparable experimental data/model conditions/uncertainty | not_measured; не входит S1;200 K calculated curve не эксперимент |

`error`/`not_measured`/`not_applicable(reason)` различать; missing не0. Nonfinite,
bad receipt/identity/completeness, U/G/n mismatch, missing required fresh audit,
unknown applicable mandatory criterion, guard/resource/cap failure и ambiguous
submit останавливают affected operation. При budget exhaustion допустим negative
conclusion; pass не получают через ослабление threshold или physics. Two identical
failure fingerprints останавливают affected wave; source/design tasks продолжаются.

## 10. Конкретный остаток и следующие ready tickets

| Ticket / owner | Различающее evidence и условие завершения | Граница |
| --- | --- | --- |
| SCI03-REVIEW / independent R high | Проверить exact D1 hash, §5 inequalities, scope of primary formulas, all strict floors/units, H within S1, unknowns never pass | CP-DESIGN, без runtime; author не закрывает сам |
| SCI03-ACC / Core+Research D | Установить source-derived finite target accuracies/physical distinguishing scale либо оставить их explicitly null; показать necessity для заявленного уровня | Не invent1e−4 и не менять user target/priority; если доказательств нет, concrete scientific blocker сохраняется |
| SCI04/O08-03 / Core D→I | Согласовать final U/G/reconstructed density/fresh-map lineage; хранить raw numerator scales и unavailable flags; literal aligned/misaligned oracle | Не повторять solver kernel как independent physical oracle |
| SCI06 / Results D→I после review | Из same-state SL/SG+A/GR/GL/GG измерить §5.3 coefficient norms и absolute corrections; вывести represented-current bound с independent literal2×2 positive/negative oracle | Existing reader storedΣ only; no fresh map claim, no criteria registry в UI |
| SCI12 / Research+Core D | Freeze directional comparison mapping, window strips/C interactions, actual finite tail/conditioning proofs or explicit unknown; verify precise H phase encoding | Только existing S1 variants; no extra executions/embedding changes |
| SCI07/08 + SCI09 / owning lines | Fresh admission, all resident/cumulative I/O/log guards и exact plans/fingerprints; только затем pilot READY subset | Текущая native evidence missing; никаких probes/submit в SCI03 |

## 11. SHA256 проверенных owning files

| Locator (указанные immutable component commits) | SHA256 bytes |
| --- | --- |
| C/src/numerics/reference/acceptance_metadata.jl | `4d1880b7d6e475534d0bbfa229ced25dd40b062f674ec8a2ea86b45e0db1dd86` |
| C/src/numerics/reference/validation.jl | `29e1b92daaed51f35d676bab8f70d60836fd0b5dc65c25c38a58646e0395005a` |
| C/src/numerics/reference/final_audit.jl | `8384e2022de0c3c912384c15f6e5b3697879dac984d554f38f0d5bc2d916a673` |
| C/src/numerics/reference/poisson.jl | `8989ced0c4adc67e5d553b98863cd8ee3584705546ba445565cb5da0e7130645` |
| C/src/numerics/reference/grids.jl | `71c4cff98cbefb834eea818067549dc150cb0cd18141854f27dfbd549862d18e` |
| C/src/numerics/reference/positivity_diagnostics.jl | `796cdde49e736245bbd3cf53aab96b33368560e7c1518e27e96ba2d0d67d8039` |
| C/src/physics/observables.jl | `b147e41dcb04bb39792a30e0092f9298a64813d87cc5411037c1958bd56b4b9f` |
| C/src/domain/types.jl | `f428a6c7ce2c7d4796d6be2983decffcdbfc25f21ec400a600038863f673e71b` |
| C/src/domain/scaling.jl | `2049ea6638da871c406718dad9f1aa4eb5cee68f4ce8941a8270e050dd7eb9f0` |
| V/src/qcl_negf_results/equilibrium.py | `3f8ee56390741e15bc2c709de67d4b0c80260733c99752d8fbbb4f538e3042e6` |
| Research/docs/reference-2019-model-passport.md | `653dc08a55575043b3b95ee43c0ff5413914f39ae6b554c9993345514d5f8f3c` |
| Research/config/model/reference-2019-70k.yaml | `78444c7a358fabef315e07e028d3d65f65fec07efd119ac9ca238cc25c110829` |
| Research/config/policies/analytical-diagnostic.yaml | `fbbb2b4afe9eac694e71b29c77883930072416250f0f75e86d8ffaf3821c49b0` |
| Research/config/studies/equilibrium-0-70.yaml | `2fde52250d6ab85240d556f8c59b70ce1a065e04800d531dea95d9a2f3a8c0e9` |

## 12. Решения D1

| Утверждение | Свидетельство | Ограничение | Решение |
| --- | --- | --- | --- |
| Existing strict policy/norms сохраняются | Core commit+hash §4/§11, Research source values | Resolved/runtime plan ещё не создан | сохранить |
| FD/FDR и аналитические DOS/current anchors определяют exact equilibrium identities | Primary version/formula/hash §2, explicit derivations §5 | Не задают finite QCL forward-accuracy targets | сохранить |
| Имеются useful absolute bounds без произвольных coefficients | §5 μ derivative, Cauchy density/current identities | Fixed represented domain; measured coefficients/data отсутствуют | дополнительно измерить |
| Полная physical/discretization acceptance уже готова | Unknown εQ, nonlinear conditioning и domain/truncation bounds | Ни strict residual, ни two-grid sensitivity их не заменяют | данных недостаточно |
| Конечный pilot имеет reviewable measurement contract | S1 unchanged slots/caps, P0/transition/H §8 | Нужны CP-DESIGN, source/oracle, admission/guards/plan bindings; launch сегодня blocked | изменить |
| Геометрия/доза/200 K curve принимают70 K решение | SCI02 passport и primary caption conflict | Assumptions/value origins и иные условия сохраняют limits | данных недостаточно |
<!-- END SCI03_REVIEWED_D2_BODY -->

## Навигация по aliases и pinned источникам

| Alias | Канонический source locator |
| --- | --- |
| C | [Core f5b73369, дерево owning equations/evaluators](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/tree/f5b73369a96d1225bd96fb34f5262c8dcf965f86); C/src/... в корпусе — путь относительно этого дерева |
| R / Research | [Research b1a3be2b, входы и SCI02 passport](https://github.com/Afonenko-QCL-NEGF/qcl-negf-research/tree/b1a3be2b0f4b3adf423cc2a611bcca2b5e0a5ff3); Research/... — owning path |
| V | [Results17574ca4, independent stored-state reader](https://github.com/Afonenko-QCL-NEGF/qcl-negf-results/tree/17574ca4c5495f74d86ed9bef51c60dea8265062); V/src/... — путь относительно этого дерева |
| P1 | [Operator permissions и S1 caps, root9f8b08bb](https://github.com/Afonenko-QCL-NEGF/qcl-negf/blob/9f8b08bbf8b2f941eaa6426bc07c2b436a2a81d8/docs/implementation/review-closure-2026-10-19/operator-permissions.md) |
| milestone | [Первый научный рубеж, root9f8b08bb](https://github.com/Afonenko-QCL-NEGF/qcl-negf/blob/9f8b08bbf8b2f941eaa6426bc07c2b436a2a81d8/docs/system/scientific-milestone.md) |
| science / handoff | [SCI01–15, norms, finite waves и C01–C08, root9f8b08bb](https://github.com/Afonenko-QCL-NEGF/qcl-negf/blob/9f8b08bbf8b2f941eaa6426bc07c2b436a2a81d8/docs/implementation/review-closure-2026-10-19/science.md) |
| architecture | [Границы owning components, root9f8b08bb](https://github.com/Afonenko-QCL-NEGF/qcl-negf/blob/9f8b08bbf8b2f941eaa6426bc07c2b436a2a81d8/docs/architecture.md) |

Standalone Research checkout содержит собственные inputs/docs, но не включает
Core/Results/root corpus или primary PDF bytes. В таблице даны pinned внешние
references; сетевой доступ к ним при этом переносе не проверялся, offline
availability вне объявленного workspace неизвестна. Все source paths таблицы
проверены локально в integration workspace. Primary URL/version/hash из корпуса
сохраняются; эти PDF не перераспространяются и не подменяются локальными копиями
в Research. Временные `/tmp` paths в историческом корпусе — locators audit receipt,
а не runtime inputs или обязательный способ чтения опубликованных критериев.

| Утверждение | Свидетельство | Ограничение | Решение |
| --- | --- | --- | --- |
| Reviewed measurement design перенесён без изменения | Exact body40 812 bytes/SHA256 и R2 receipt | Independent owning CP-CODE ещё нужен | сохранить |
| Physical/discretization acceptance или runtime READY получены | Неизвестные accuracy/bounds/plans/guards/admission | Source contract не numerical result | данных недостаточно |
