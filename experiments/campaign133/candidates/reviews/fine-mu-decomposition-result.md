# Independent result review: fine μ attribution

**Удержание source μ устраняет old-node collocation difference до
наблюдаемого floating уровня, но практически не меняет dominant defect
на новых midpoints. Следовательно, μ refit не объясняет полный измеренный
fine-grid defect этого frozen Γ.** Это finite-operator counterfactual,
не доля заряда, не SCBA convergence и не доказанная причина stagnation.

Проверены малые TOML/CSV/log/source hashes; Julia/HDF5 arrays reviewer
не использовал. Saved audit `fine-mu-result-audit.py/.json/.log` PASS,
self-enforced128MiB,0.20s/26768KiB; дополнительная проверка admission/
Manifest также PASS. Source equation отдельно прослежена в static review.

## Source relation и выполненная работа

Runtime input SHA39a67a1ee3c3eca4f948d02b6781772703ce7d603169cd9ccd3b75b6694924e0,
main68ccd9d03dc6e4861c31798bad714ca6ed0b74a129882d44b2fa5c4d938cb8f6
и copied tiny/protocol hashes independently совпали с summary/current
freeze. Seven source hashes заполнены, совпадают input/summary; четыре
HDF5 hashes совпадают с previously checked recovery metadata, три малых
source metadata hashes пересчитаны. Refinement definitions5084/7d71/565a,
Manifestda0780fb и original summary3135bba1 согласованы. Большие files
reviewer повторно не хешировал; actual driver проверил их до/после work
и записал sources_unchanged_after=true.

Actual counters:1 retarded,2 LO,0Φ/Pi/native. Fresh N при обеих stored μ
и LO(μf) Frobenius norm reproduce saved values **exact0 error**. Это
контролирует fresh source binding/operator, не научную приёмку. Root
correct `--tiny-mu`16/16 PASS,16.23s,815940KiB; прежний root `--tiny`
завершился usage error до tests/source calls. Этот command failure
сохранён отдельно; он не был изменением алгоритма или failed LO oracle.

Actual log exit0:46.71s process,41.97005s script,4456660KiB≈4.25GiB,
Swap0. Root admission OS/numerical12.125GiB,High11.5GiB,J1/BLAS1,
CPU100%,wall300s,output64MiB согласованы. No extra calls/retries.

## Уравнение и meaning

При фиксированных fine Γ/R/A/cache thermal-LO broadening affine по f:

    ΓLO(f)=K[(NLO+1)S−A + NLO S+A + S+(fA)−S−(fA)].

Поэтому M=LO_f(μf)−LO_f(μc)=K[S+(ΔfA)−S−(ΔfA)], а elastic
contractions не зависят от μ. Для saved D=Φf(μf)−PΦc(μc) верно
F=D−M=Φf(μc)−PΦc(μc), без третьей LO/новой Φ/fit. This identity
относится к тому же finite native operator и actual fine A. Все complex
entries и weights wE_f·wk включены; old-node partition тоже использует
fine weights (coarse weights вдвое больше). Denominator table —PΦc,
не Φf и не target N.

μc=.08459613668783181, μf=.08459119027608916; Δμ=−.4946411743μeV
при E0=.1eV. Удержанный μc на fine spectrum даёт N=.045003199994375624,
то есть N/target−1=7.11109861249e−5 (.0071111%). Это сохранённый
counterfactual number, не number repair или новая accepted state.

## Независимо выведенные ratios

Проверены18 scalar moments/full=old+mid partitions, square-root/relative
definitions, complex inner-product identities/Cauchy bounds и agreement
D с recovery. На5445 selected rows independently проверены matching
coordinates/D values, nonnegative squares, D²=M²+F²+2Re〈M,F〉,
Re〈D,M〉=M²+Re〈M,F〉 и trace bound. Все checks passed.

| Fine measure partition | ‖D‖w/‖PΦc‖w | ‖F‖w/‖PΦc‖w | ‖F‖w/‖D‖w |
|---|---:|---:|---:|
| Old nodes |1.12262133577e−6 |8.14611604956e−15 |7.256334607e−9|
| New midpoints |3.75366264163e−3 |3.75365977600e−3 |.999999236578|
| All fine |2.65398310487e−3 |2.65398096003e−3 |.999999191838|

Таким образом old-node D почти целиком совпадает с M; остаток F имеет
weighted absolute9.37735e−16. На midpoints остаётся99.9999236578%
исходной нормы D, глобально99.9999191838%. Weighted ‖M‖/‖D‖ всего
.000299114 на midpoints и.000423023 globally. Это нормы разных
векторов; их отношение не является аддитивной attribution fraction.

Signed inner products подтверждают вывод без такого допущения:

- full Re〈D,M〉w/‖D‖²w=8.97635876e−7, cosine≈.00212195;
- midpoint projection coefficient8.08156257e−7, cosine≈.00270183;
- old-node projection≈.999999999984 и cosine≈1.

При вычитании M squared norm изменяется на
(2Re〈D,M〉−‖M‖²)/‖D‖²:1.61632291e−6 globally и1.52684335e−6
на midpoints. Это algebraic signed norm reduction, **не fraction of
charge/physical error**. Для old nodes вычитание почти равных quantities
теряет относительную точность: rounded moments дают cosine чуть выше1
(≈1+1.3e−14). Это floating summation scale, не physical violation и
не повод клиппировать исходные moments. Direct F norm надёжнее
catastrophically cancelled «100% removed» estimate около нуля.

Важно: устранение old-node **D** не устраняет raw fixed-point residual.
Old-node R−M relative=5.95754891323e−6 кPΦc, то есть восстанавливает
прежний coarse raw residual, а не ноль. Full ‖R−M‖w/‖R‖w=
.999999193209: global raw residual тоже практически не меняется.
Unweighted norms дают тот же qualitative result (full F/D=.999998630883).

## Границы grid-vs-μ вывода

Recovery установил99.999991% weighted D² на midpoints. Этот run
показывает, что этот dominant component сохраняется при held μ; значит
его нельзя объяснить scalar μ refit. Результат поддерживает существенную
**frozen collocation sensitivity при добавлении nodes**, не обещание,
что fine SCBA сойдётся или что residual stall вызван только energy grid.
Endpoint ghost-width variant, interpolation/retarded/native coupling и
floating policies остаются частью controlled comparison; отдельная
continuous-causality qualification не появилась.

Selected μ CSV по-прежнему содержит495E×11k, old9785 LO ladder.
Он покрывает5.18917214% weighted D² и5.18917795% weighted F²,
лишь3.19525903e−6 полного weighted M². Global maximum recovery
e39394≈469.820607meV отсутствует. Поэтому selected CSV не устанавливает
главную energy localization ни D/F, ни M. Global/partition выводы
выше опираются на all-node streaming moments; узкие spectral shapes,
linewidths, Jacobian/nullmode и causal mechanism этим не определены.
