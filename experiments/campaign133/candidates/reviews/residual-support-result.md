# Независимый разбор residual-support cont128

**Raw residual локализован преимущественно в почти незанятой части окна
и в нескольких узких группах E-узлов. Это локализация, а не установленная
причина stalled convergence. Занятая часть решения тоже не принята.**

Проверены только сохранённые CSV/TOML/log и малые source identities.
Reviewer не читал HDF5 и не запускал Julia, новый Φ, inverse или fit.
Воспроизводимый bounded stdlib audit:
`reviews/residual-support-audit.py`; первичные числа и hashes:
`reviews/residual-support-audit.json`. Audit all checks passed,0.18s,
peak39060KiB, self-enforced address-space cap256MiB.

## Что независимо сверено

Actual `residual-support/attempt01` имеет copied script SHA b77dac503c3bf35e9bc7bc927242fbcbdeeafe7dec608e91ab5e0c05fb81edee
и input46a74dd52a61b81de7d97c994f6ecb76ae2916d07a07208d3a8d31096c5fcc18,
равные frozen requested artifacts и summary identities. Малые map/input/
identity/writer hashes совпадают с summary. Большие model/terminal hashes
сверял actual probe перед streaming; reviewer их повторно не читал.
Source — accepted Core e171/tree742ed5a6, Hermitian-Λ PI finite-chain4,
NE30704,Nk11,Nb5, map128; source не сошёлся.

Audit проверил30704 perE rows,11 perK rows и20 selected joint rows:
finite values; суммы raw/weighted numerator и scale; обе маргинальные
нормировки долей; извлечение всех шести norms из squared sums; replay
source map; независимый расчёт трёх зон и stable-Fermi f² columns;
полное f² numerator/scale/relative; selected E/k coordinates, f,
fraction/width-over-spacing и descending order. Максимальное относительное
отличие replay от source norm меньше4.6e−12 при допуске1e−10.

| Raw norm | Replayed |
|---|---:|
| absolute |1.245095667631454e−4|
| scale |17.31149651993275|
| relative |7.192305218661078e−6|
| weighted absolute |9.69867133467482e−7|
| weighted scale |0.16279633580611688|
| weighted relative |5.957548913278676e−6|

Measure wE·wk и reference raw Φ сохранены; нет2π/spin, working-map
Hermitization или residual floor. Arithmetic replay не является новым
scientific tolerance. Сумма squared weighted numerator9.406422565804306e−13.
Joint top20 accounting согласован с маргиналами и frozen algorithm;
без HDF5 reviewer не пересортировал независимо все joint(E,k) blocks.

## Энергетическая и k локализация

Stored μ=8.4596136688meV, kBT=6.0321332835meV. Зональные границы
μ±5kBT:−21.7010527487 и38.6202800863meV относительно source Eref.

| Зона | Доля squared weighted raw numerator |
|---|---:|
| E<μ−5kBT |0.0000311868%|
| |E−μ|≤5kBT |0.4643910041%|
| E>μ+5kBT |99.5355778091%|

Один E-index9785 (82.399148meV) несёт36.19420535%; его сосед9786
(82.438232meV) —18.62291816%. Вместе54.81712351%. Top10 E-rows
несут90.73771414%, top20 —96.48296198%, top100 —99.92218747%.
Это fractions squared norm, не linear norm или charge.

Largest joint node(e9785,k4,κ1.8) несёт5.77091079% полного weighted
numerator; все joint top20 вместе56.52427743%. Они распределены по
нескольким k: максимум perK приκ2.4 —14.53242572%; диапазонκ1.2…3.6
содержит66.43490538%. Следовательно, support узок по E, но не сводится
к одному k. k=0 сохраняет0.90831446%, а не искусственный ноль.

В source `build_context` используется reference_parameters; accepted
Core `src/domain/presets.jl:72` задаёт LO36.7meV. Здесь939 energy steps
равны36.7meV (defect≈−6.9e−18eV). На той же grid phase действительно есть
группы, совпадающие с целыми LO shifts. В таблице каждая доля суммирует
**три** E-rows(anchor−1,anchor,anchor+1), не всю спектральную линию:

| Shift от e9785 | Anchor E,meV | Доля squared weighted numerator |
|---|---:|---:|
| −LO |45.699148 |0.47847330%|
| 0 |82.399148 |54.81712351%|
| +LO |119.099148 |15.72680439%|
| +2LO |155.799148 |9.25088337%|
| +3LO |192.499148 |0.27132251%|

Другой видимый anchor e9182=58.831416meV имеет9.15196226% в ±1row.
Его ближайшие LO-aligned anchors22.131416 и95.531416meV содержат
0.23856550% и0.92893227%; −2LO попадает на−14.568584meV. Совпадение
с LO step согласуется с energy-shift coupling, но snapshot не отделяет
действие LO operator от спектральной структуры, вклада других channels
или исторической динамики iteration. Это не causal attribution LO.

## f² diagnostic и selected spectral annotations

При f² reweighting absolute=3.4924644934e−8, scale=0.00185866004757,
relative=1.87902273898e−5. Relative больше global weighted relative:
нельзя из уменьшения absolute заключать, что occupied solution принят.
Эта величина вообще не является error bound для N/current или native
SCBA acceptance.99.97377690% f² numerator находится в средней зоне.
Пять соседних rows7303…7307 (−14.607668…−14.451331meV) несут90.11502356%
f² numerator; крупнейший row7304=−14.568584meV —31.70466952%.
Они почти не видны в unoccupied-dominated global ranking, но не исчезают.

У largest joint node82.399148meV,κ1.8 source f=4.74882058e−6,
traceΓ/Nb=30.8390463μeV и ΔE=39.0841321μeV, ratio=.78904263.
На соседнем E-index9786 при том же k proxy уже995.438092μeV,
ratio25.46911086. Это сильная локальная вариация sampled Γ;
traceΓ/Nb не является linewidth certificate или eigenchannel linewidth.
Selected HermΓ min eigenvalue в первом node1.82196409μeV, anti-Hermitian
norm0. Сохранённый condition там32.54148; диапазон selected top20
conditions≈3.6143…95.1148. Это не доказательство, что stall вызван
плохой обусловленностью; source map global maximum≈1280.616 относится
к другому, не локализованному этим CSV node. Atrace и condition доступны
и finite в top20, но не составляют полный spectral-resolution audit.

## Цена и границы вывода

Root log: actual exit0, wall8.99s, user8.07s/system.59s,
peak623868KiB≈.595GiB, swaps0. Summary: full-file SHA2.47195s отдельно,
последующий diagnostic4.34646s; остаток общего wall включает startup/JIT
и orchestration. Admission фиксирует MemoryMax2GiB,High1.75GiB,Swap0,
Julia1/BLAS1,timeout180s. Output attempt01=5,415,453bytes≈5.16MiB,
меньше32MiB. Измеренный короткий read-only snapshot diagnostic не
оценивает цену новой SCBA/Φ iteration или full-state solve.

Поддержан вывод о concentration raw residual в высокоэнергетических
узлах и нескольких LO-aligned группах, при отдельной оставшейся
f²-weighted структуре в термически занятой области. Не установлены
Jacobian eigenvalues/nullmode, причина или способ лечения stall,
linewidths, достаточность ΔE/window/Nk, native SCBA convergence,
Poisson solution или density acceptance. Для причинного объяснения
нужно отдельное различающее свидетельство, не переименование snapshot.
