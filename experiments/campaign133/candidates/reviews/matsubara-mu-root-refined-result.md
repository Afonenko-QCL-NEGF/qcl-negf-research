# Независимый результат: fixed-charge Matsubara ALL5 + linear M3, N512/1024

**Вердикт: подтверждены численная стационарность конечных Matsubara-отображений и заранее заданные проверки заряда нового linear-M3 варианта.** При одном найденном μ=0,08458841693811386 native = 8,458841693811387 meV выполнены charge target на N512, независимый N1024 witness и двухотсечное сравнение числа частиц. Это не приёмка исходного real-energy baseline и не доказательство решения при N→∞. Сохранённые `native_pass`, `scientific_acceptance`, `infinite_N_certified` остаются false.

Первичный артефакт: `adjacent/matsubara-mu-root-refined/attempt02/summary.toml`, SHA256 `e041cca82a7e9c08aa1a823fe24c09f906f87249d634c87cd4d11ab86f9e8bf3`. Независимый audit: `matsubara-mu-root-refined-result-audit.py/.json/.log`; ограничение AS128MiB, фактически 0,29s/31384KiB. Прочитаны малые TOML всех 7 calls/112 итераций и binary Σ/Φ только двух финальных `latest.h5` через h5dump. Новых Julia, Φ, Dyson, μ-fit и чтения real-energy Γ не выполнялось.

## Что решалось

Для каждого μ заново от Σ₀=0 выполнены 16 несмешанных шагов α=1 конечного отображения Σ_N=Φ_N(Σ_N;μ). Модель сохраняет all5: четыре direct-static канала и LO с full-line bare reference и конечной Matsubara-суммой разности. Используется finite45/open9-cell Hermitian hH модель, Nb5/Nk11, local self-energy повторяется по всем девяти ячейкам; это не независимая cellwise45 SCBA. T=70K, U=V=0, E₀=0,1eV, gₛ=2; физические матрицы/вершины/радиальные веса сохранены. Изменение μ сдвигает относительные полюса на −Δμ, а не меняет физический hH.

Raw density и отдельная asymptotic linear-M3 closure имеют вид

\[
\rho_N=\rho_0(\mu)+\frac{2}{\beta}\sum_{n<N}\operatorname{Herm}(G_n-G^0_n),\qquad
(I+2\Omega t_4 K_{LO})\rho_{M3}=\rho_N+t_4 B_0(\mu),
\]

\[
\mathcal N=g_s\sum_k w_k\operatorname{Re}\operatorname{tr}\rho_{M3}(k).
\]

Таким образом подтверждается стационарность конечного Φ_N и fixed-charge условие именно для отдельно определённой M3-плотности. Нет перенормировки ρ, PSD-clipping или подстановки старого real-energy Γ. M3-коррекция не является утверждением о точном interacting tail всех пропущенных частот.

## Счётчики, μ и невязка

Сверены все сохранённые requests, histories и per-iteration TOML: 6 outer N512 + 1 witness N1024, ровно 16 шагов каждый; 112 callbacks и 112 reserved maps при структурном пределе 144. Нет скрытого подбора μ внутри этих calls. Начальная скобка μ_source±0,001 native имеет разные знаки charge defect; все четыре последующих secant proposals воспроизведены независимо. Call6 выбирает μ; call7 запускается с Σ₀=0 при **том же** μ. Оба финальных состояния относятся к итерации16, а не к смешанному следующему шагу.

Квалификация inner call требует одновременно raw Frobenius relative≤1e−9 и absolute≤1e−12 native в последних трёх итерациях. По всем семи calls максимумы соответственно 1,596342709820291e−18 и 1,5831138727864743e−19. Relative использует норму raw Φ; wk-weighted норма сохранена как отдельная диагностика, не подменяет критерий.

В двух финальных HDF5 независимо прочитаны полные Σ_real/imag и Φ_real/imag: все значения finite, обе пары **побитово равны**, поэтому их raw absolute/relative defect точно 0 в сохранённой Float64 арифметике. Независимые Frobenius scales Φ равны 0,09917130345227167 и 0,09931957240899796. Проверена ориентация HDF5 physical `(5,5,11,N)` против Julia logical `(N,11,5,5)`.

| Проверка числа частиц | N512, call6 | N1024, call7 |
|---|---:|---:|
| Raw finite-N number | 0,0449999788688657 | 0,044999997348664215 |
| Linear-M3 number | 0,04499999999896239 | 0,044999999989930074 |
| (M3−target)/target | −2,3058098640e−11 | −2,2377624445e−10 |
| Заданный charge tolerance | 1e−9 | 1e−8 |

Target=0,045000000000000005 native. (N1024−N512)/target=−2,0071814580783587e−10 при заданном cutoff tolerance1e−8. Raw-number ошибки −4,6958076e−7 и −5,8918573e−8 сами эти charge tolerances не проходят: M3 closure существенна для этого результата и не должна замалчиваться.

## Независимые значения M3 и Pauli

Отдельный буквальный шестии́ндексный audit всех пяти K выполнен в `matsubara-mu-root-refined-value.md` (SHA256 `ee0f9940a46640581aadc346ec35cfd3d4dcc6f6eedfcd20a7e6d76ec0128aa4`). Он подтверждает source kernels, wk, hH, однократный сдвиг новым μ, S₀/B₀ relative differences 2,29e−16/2,71e−16, исходные и исправленные числа частиц. Невязки linear-M3 equations независимо получены 7,94e−17/1,03e−16; conditions≈1. Значения η≈1,07e−10/1,34e−11 характеризуют именно малый linear solve, не ошибку infinite-frequency физики.

Полные correctedρ проверены без native eigensolver, real10×10 Jacobi для Herm(ρ): minima 7,151255230137108e−8/7,151170845213298e−8, maxima 0,2174988717669255/0,21749887176689767; eigen residual≤8,84e−17. Hermiticity defects≈1,81e−17/2,46e−17. Это численное Pauli свидетельство конечной correctedρ, без ремонта массива и без interval certificate. Matrix difference двух correctedρ Frobenius=3,12941947712e−12; эта норма не заменяла charge criterion. M3−one-shot Frobenius≈3,65e−17/3,27e−17; преимущество feedback над one-shot здесь не измерено.

## Provenance и ресурсы

Сверены copied/source/input SHA и исходные bindings. Рабочий wrapper V2 SHA256 `cb1d60a42c0be18761d5fa8e69846c1e5c8a81a08eb8b5c95ebc22483c31f144`; actual input `0057daf02b9c5fd80f0384350d6c18de35ac07795ed0369dc9045175459a8191`. Refined module `ecf7b1cce70ecebe054905cb0e53be7c5fedd91c3061ed9b2634af4bbbdd1998`: обратная замена имени, doc domain и allowlist512/1024 восстанавливает parent1434 byte-for-byte. Original source adapter94c6, parent1434 и M3 module9182 не менялись. Accepted Core e171caa6/tree742ed5a6 и native environment bindings сохранены.

Полные финальные latest SHA256: N512 `c5d83670abd7fce35d6474ffcee548d81feae50ce3789ed2a9eeb076c24f755f`; N1024 `6877b34851dcead3300e133b203249c88bcf2375a5b62c8f201ab65304d6d816`.

Root actual log: exit0, script84,953s/process85,33s, peak952148KiB, swap0, данные635703192B; пределы2GiB/1200s/1GiB data/file32MiB соблюдены. Julia1/BLAS getter1. До actual пройдены source-free13 tests и отдельный loader4 controls. Полная новая генерация 112 maps независимым review не повторялась.

История сохранена: старый128/256 controller UNRUN; его source cutoff недостаточен для прежнего target, но это не failed run нового μ. Отдельный ошибочный command и attempt01 world-age failure завершились до model/maps; attempt01 пуст. Единственная runtime V2 поправка — `Base.invokelatest(() -> Main.load_modules(sconf))` после dynamic include; математика/лимиты неизменны. Предварительный tiny с неверным qualified module name также сохранён отдельно от успешного corrected command.

## Пределы вывода

Результат даёт содержательное **численное равновесие нового all5 direct-static/full-line-LO/reference-tail M3 finite-N fixed-charge варианта**. Точное сохранённое Float fixed point и согласие512/1024 не доказывают математическую единственность, N→∞ SCBA, строгий tail remainder, отсутствие общих ошибок алгоритмов или точность исходной real-energy дискретизации. Pauli correctedρ не является сертификатом causal real-frequency self-energy. Нет проверки basis/k/cell convergence, Poisson, finite-bias transport, экспериментальной валидации или приёмки native real-energy baseline. Ток не измерялся: предложенный whole-cell Matsubara bond-current control алгебраически нулевой для произвольной локальной Σ и не дал бы независимой transport/conservation приёмки.
