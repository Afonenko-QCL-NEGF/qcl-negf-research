# Спектральная плотность при слабом уширении: стратегии открытых кодов

2026-11-16. Read-only: существующие clones, без сборок, пакетов, новых checkout,
реализации и runtime. Применён ранее прочитанный qcl-physics-audit. Исходный
superproject HEAD `b515c42ef14c1f61151626280da0bcf100d116be`; существующие чужие
изменения gitlinks сохранены. Сообщённый ведущим Γ-seed имеет минимальное
положительное значение 1.26e−27; **это ещё не измерение минимального собственного
значения Γ и не диагноз отрицательного q**.

**Наиболее полезен прецедент Kwant: LDOS из суммы квадратов incoming-channel
solutions. Для finite4 его аналог — Gram-сумма ВСЕХ локальных и cavity sources,
сохранённых независимо от малой anti-Hermitian части retarded resolvent.**
libNEGF/Quantica generic LDOS и NESSi FDR используют imaginary part Gr; готовой
η-free защиты от потери причинности в просмотренных путях не найдено.

## Что действительно делают pinned sources

| Код / ревизия | Алгоритм и точный источник | Применимость |
|---|---|---|
| Kwant `e1603f8817e1459f2544ae84831e78a3e4b79430` | [`common.py:552–561`](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/solvers/common.py#L552): один factorization, solve для incoming RHS, сумма `abs(psi)^2/(2π)`. [`leads.py:865–867`](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/physics/leads.py#L865): normalization на sqrt(abs(velocity)). | Положительная диагональ по конструкции; не вычитание Gr−Ga. Сам LDOS отвергает non-Hermitian H и selfenergy-only leads ([530–540](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/solvers/common.py#L530)); прямо импортировать API для SCBA нельзя. |
| libNEGF `3b4120f65b19107149a46034ecd268be00622042` | [`sparsekit_drv.F90::zspectral_dns:2854–2871`](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/sparsekit_drv.F90#L2854): буквально `i*(GreenR-GreenR†)`; CSR аналог :2812. [`iterative_cpu.F90:230–329`](https://github.com/libnegf/libnegf/blob/3b4120f65b19107149a46034ecd268be00622042/src/iterative_cpu.F90#L230): Gn через injection и Gr/Ga, с recursive Sigma_n. | Keldysh propagation — полезная форма; generic A не защищена Gram factorization. Внешние источники и local scattering суммируются, поэтому только lead Γ недостаточно. |
| Quantica `ed2bb1ce1598cf5e03f721e70fe1b547f8103473` | [`observables.jl:113–127`](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/observables.jl#L113): LDOS `-imag(G)/π`. [`spectrum.jl:123–129`](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/spectrum.jl#L123): density через Σ f_n ψ_nψ_n†. | Spectral sum положителен при f≥0, но `SpectralSum` ограничен bounded 0D Spectrum solver и отсутствием self-energy (:105–107); не замена energy-dependent interacting Σ. |
| NESSi `6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6` | [`ness2_dyson_impl.hpp:80–88`](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_dyson_impl.hpp#L80): direct inverse и `Gless=Gr*(Sigmaless+bathless)*Gr†`. [`ness2_herm_matrix_ness_impl.hpp:448–480`](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_herm_matrix_ness_impl.hpp#L448): `force_equilibrium` из M−M†. | Keldysh multiplication, но не factorized-positive spectral algorithm. Опциональный bath имеет согласованные ret/less components, а не невидимый retarded-only η. |

Kwant проверяет Hermiticity **входного H** относительно его масштаба
[`common.py:158–165`](https://gitlab.kwant-project.org/kwant/kwant/-/blob/e1603f8817e1459f2544ae84831e78a3e4b79430/src/kwant/solvers/common.py#L158).
Это проверка, а не symmetrization результата. Явного восстановления Hermitian Λ
после Hilbert FFT или общего PSD repair в перечисленных spectral paths не найдено.

Отдельно о регуляторе: Quantica для real frequency добавляет `5eps(T)`
([`greenfunction.jl:79–86`](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/greenfunction.jl#L79));
при неразрешимых Schur modes увеличивает Imω
([`schur.jl:390–415`](https://github.com/pablosanjose/Quantica.jl/blob/ed2bb1ce1598cf5e03f721e70fe1b547f8103473/src/solvers/green/schur.jl#L390)).
Это способ выбора retarded modes/регуляризации, **не доказательство η-independent
физики при Γ≪eps·energy_scale**. NESSi BATH_CONST явно задаёт Σr=−iη и
Σ<=2ifη ([`ness2_dyson_impl.hpp:25–35`](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/ness2_dyson_impl.hpp#L25)).
Вводить такую bath в S01 — отдельный physical variant.

## Практический кандидат для S01: структурированная Λ плюс Gram A

Полный optical-theorem finite-chain proof выполняет другой пакет; здесь только
precision strategy. Хранить Hermitian Λ и PSD Γ раздельно. При Hermitian input
Hilbert-оператор с действительными weights должен сохранять эту симметрию;
complex FFT может добавить roundoff anti-Hermitian leakage порядка
`eps*norm(Λ)`, потенциально больше малой Γ. Это **гипотеза**, пока не измерены
локаторы отрицательного q, величины и нормы. На диагонали i(Gr−Ga)=−2Im Gr:
переписать формулу через `imag` недостаточно — Im мог быть испорчен уже solve.

На фиксированных E,k диагностировать raw `Λ−Λ†`, Γ-local eigenspectrum,
source-derived Γ-cavity eigenspectrum, solve backward residual и разность
A_sub−A_Gram. Λ-Hermitian восстановление допустимо как реализация контрактной
симметрии, с сохранением исходного defect; это не PSD clipping. **Γeff нельзя
восстанавливать как i(Σeff−Σeff†), если малое уширение уже потеряно при сложении
с Λ**: переносить injection/broadening через cavities непосредственно. Для
одной source factor Γ=B B† перенос — X=g B, contribution X X†; собирать сумму
локального центра и обоих cavities. Γ, необходимая в центральном Dyson, должна
соответствовать тем же источникам. Веса density — исходные положительные wE,wk;
A имеет единицы 1/energy, q=Tr A/(2π), далее исходная spin/phase-space мера.

```julia
raw_herm_defect = norm(Λ - Λ')
Λh = Hermitian((Λ + Λ') / 2)                  # defect сохраняется
B = all_local_and_cavity_source_factors(state) # Γeff = B*B'; без clipping
Γeff = B * B'
M = E*I - H - Λh + 0.5im * Γeff              # η = 0, прежняя finite4
Y = M \ B                                    # solve, без inv
Agram = Y * Y'
q = sum(abs2, Y) / (2π)                       # trace, явно неотрицателен
save_audit(raw_herm_defect, M, Y, B, Agram, A_subtraction)
```

Псевдокод не реализован. При factorization локальной Γ нельзя молча обрезать
отрицательные eigenvalues или выбрасывать малые positive modes: сначала
различить input noncausality и roundoff. Gram сохраняет положительность
представления, но **не исправляет плохо решённый M** и не доказывает близость к
верному resolvent. Аудит нужен при backward-error масштабе исходного оператора,
не только относительно 1e−27. Настоящие bound/dark states с Γ=0 дают singular
resolvent и δ-вес: Gram continuum не восстанавливает такой вес автоматически.
Прежние spectral-sum/KMS/native gates остаются обязательны; raw subtraction
сохранять как независимую численную диагностику, не делать невыполнимый relative
sign gate для физического веса ниже backward-error floor.

Для подтверждения precision improvement достаточно сначала нескольких
**конкретных проблемных E,k blocks** с теми же H,Λ,Γ, finite4 и η=0: сравнить
обычную и structural/Gram формы с higher-precision block solve. Выбор точек
зависит от ещё ожидаемой диагностики; новых физических прогонов не предложено
как входное условие. Exact arithmetic equivalence — кандидат optimization;
введение η, clipping Γ или устранение физически значимой anti-Hermitian Λ —
смена модели/closure, требующая отдельной identity.

Первичный методический контекст: [Groth et al., Kwant, NJP 16,063065 (2014),
arXiv v2](https://arxiv.org/abs/1309.2926v2). Конкретный Gram-LDOS алгоритм выше
установлен по pinned code, а не приписан статье по одному abstract. Никаких
утверждений о баге upstream или причине текущего S01 rejection здесь нет.
