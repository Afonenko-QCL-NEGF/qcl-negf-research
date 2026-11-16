# Independent limited value audit: refined μ-root M3

**Вердикт: PASS финальных charge/M3 value checks нового finite-N варианта.**
На одном новом μ=0,08458841693811386 native (8,458841693811387meV)
независимый пересчёт подтверждает заданный charge target и number-cutoff
witness N512/N1024. Принят именно ограниченный fixed-charge linear-M3
asymptotic variant с frozen physical/K/hH frame. Original real-axis
baseline, infinite-N SCBA, Poisson и эксперимент этим не принимаются.

Primary: `adjacent/matsubara-mu-root-refined/attempt02/summary.toml`
SHAe041cca82a7e9c08aa1a823fe24c09f906f87249d634c87cd4d11ab86f9e8bf3;
финальные call-6-N512 и call-7-N1024 `density-M3.h5`, `density.toml`,
latestρ, frozen module/input/source packet. Wrapper cb1d60a4…,
Refined module ecf7b1cc… и original1434e86f…/M3 module9182ae80… проверены
по archived/current SHA. Qualification/review source packet полностью
согласован; все указанные text, definition и M3 actual-source SHA проверены.
Source model700ae176… и actual saved hH42fe0a38… не заменялись;
Γ source arrays audit не читает и старый Γ64 PD proof не переносит на
новый nonlinear Mats result.

**Область независимой проверки.**
`reviews/matsubara_mu_refined_value_audit.py` — только stdlib+H5dump,
AS≤128MiB/alarm30s; фактически1,05s/41656KiB/exit0.
Сохранены `matsubara-mu-refined-value-audit.json` и
`matsubara-mu-refined-value-h5commands.json`. Новых Julia/Φ/Dyson/μ-fit/N
запусков нет. Primary rawΣ/G/Φ/112-map convergence проверяет отдельный
review; данный audit не подменяет его.

Малые density arrays прочитаны полностью. Rawρ в `density-M3.h5` точно
совпадает с final latestρ. Проверены source logical SHA всех5 kernels/qK,
wk и hH; выполнена **независимая буквальная шестии́ндексная свёртка**,
без native contraction helper. При новом μ hrel=hH−μI вычитается один раз;
восстановлены все4 static+LO S0/B0 и оба порядка hS0+S0h.
Relative differences: S0=2,28774e−16, B0=2,70550e−16. Это подтверждает,
что correction действительно соответствует **новому**, не source μ.

Все charge numbers независимо получены как gs·Σk wk·Re trρ(k), gs2,
target0,045000000000000005 и source native measure. Никакой повторной
нормировки density или изменения wk/spin/E0/L0 не сделано.

| Финальный показатель | N512 | N1024 |
|---|---:|---:|
| raw finite-N number |0,044999978868865705|0,04499999734866421|
| one-shot M3 number |0,04499999999896239|0,044999999989930074|
| linear-M3 number |0,04499999999896239|0,044999999989930074|
| relative error к target |−2,3058098640e−11|−2,2377624445e−10|
| M3 correction Frobenius |5,685849611e−9|7,107322166e−10|
| M3−one-shot Frobenius |3,65116861e−17|3,26636644e−17|
| independent linear relative residual |7,9393231115e−17|1,0265408151e−16|
| independent weighted relative residual |9,4812985370e−17|1,1247555586e−16|

M3 number1024−512=−9,03231656135e−12, relative target
−2,00718145808e−10; matrix Frobenius difference3,12941947712e−12.
По **заранее заданным** charge thresholds1e−9/1e−8 и number-cutoff1e−8
все три проверки проходят. Это measured two-cutoff witness, не строгий
error bar infinite-N. Matrix Frobenius и radial-weighted charge — разные
нормы; малость первой не использовалась вместо charge criterion.

**Linear closure/Pauli.** Direct residual пересчитан буквальным KLO из
savedρ в уравнении `(I+2Ωt4KLO)ρ=ρraw+t4B0`. Сохранённый one-shot также
воспроизводится. η≈1,07309e−10/1,34136e−11 и reported cond≈1 поддерживают
устойчивость этого small linear solve; его conditional CP/Neumann bound
не является finite-N/infinite-N physics error certificate. One-shot и
linear charge здесь совпадают до Double digits; дополнительный feedback
на уровне roundoff не доказывает преимущества новой feedback closure.

Pauli spectra проверены отдельным real-symmetric Jacobi алгоритмом для
10×10 representation Herm(ρ), без LAPACK/Julia. Контроль — аналитическая
complex2×2 matrix; для всех11 k residual≤8,84e−17, duplicated-eigenvalue
defect0. Получены minima7,151255230137108e−8/7,151170845213298e−8 и maxima
0,2174988717669255/0,21749887176689767, согласованные с исходными diagnostics.
Hermiticity defects1,80523e−17/2,46213e−17. Это независимая **численная**
Pauli проверка, не interval/eigenvalue certificate; Herm здесь только
диагностическое представление, сохранённая density не ремонтировалась.

**Предел вывода/выполнение.** Metadata сообщает6 outer calls плюс1 N1024
witness,112 maps из лимита144,BLAS1 и scientific/native/infinite_N flags
false; эти границы согласуются с квалификацией именно нового finite-N
fixed-charge M3 варианта. Глобальные false flags не означают провал
выполненных локальных charge checks. По metadata/logs runtime≈85s,
peak952148KiB/swap0/data635703192B. Самостоятельной новой solver проверки
или воспроизведения112 maps этот ограниченный value audit не выполнял.

M3 остаётся asymptotic closure, а конечный LO map использует reference tail.
Higher moments/remainder, interacting omitted-frequency feedback и
infinite-cutoff error не сертифицированы; совпадение512/1024 не доказывает
их отсутствия. U=0 fixed-Hartree, saved hH-all9, finite45, full-line
reference-tail ALL5 и новая μ/closure identity должны сохраняться в любом
дальнейшем утверждении. Raw-h native baseline/Poisson/эксперимент — отдельные,
здесь не измеренные задачи.
