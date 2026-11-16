# Final independent result: one-time FDseed + full native128

**Выполнение и initializer-control подтверждены; обе ветви не сошлись.** Primary `native-bias1mV-FDseed128-max26/attempt01/summary.toml` SHA256 `f7b2f4fe194744d004b2f71925e8b70c6235e6651b911792e37ef1f10f62e21a`: exit0,128completed full native map evaluations/case, status=max_iterations, quality=unresolved, converged=false/scientific_acceptance=false. Все256progress records имеют fixed_point_passed=false/streak0. Исходный coarse32 использован как seed, не как принятое решение. По поручению пользователя дальнейших опытов нет; этот review завершает анализ текущего расчёта.

Научный вариант: matched0/1mV **на период**, all5/no correlation e–e, fixedU0/noPoisson, N24001[−0.3,0.9]eV/50μeV/Nk11/Nb5, direct-static C/LO PIFFT/finite-chain4/70K/αΣ0.2. Один global-FD reset **только исходных L/G**, все R/GR/A сохранены; затем fullnativeR/L/G maps без дальнейшей FD-проекции. Execution policy J26/BLAS1/residualworkers4, others retainrequested≤26. Finite-bias globalFD — искусственное начальное приближение, не equilibrium condition или новый reservoir.

## Audit и provenance

`native-bias1mV-FDseed128-max26-result-audit.py/.json`: stdlib Python AS≤128MiB, все256CSVrows,34marker/balance files (1,8,…,128/case), обаseed/boundary/initialization reports, identities/copiedsource/Manifest иphase log. **Никаких новых Julia/Pkg/Φ/μ-fit или чтения больших HDF arrays.** Eigenvalue/full-field/hash assertions seed-ledger — свидетельство pinned исполнившегося source; reviewer проверяет ledger/source/арифметику, не повторяет array eigensolves.

Actual source hashes сверены: driver `3dc4d1d2c6b81e6667eea6424a321fd54395879a8e28b73bf4fd036ba5b746a6`, loader `47a5d267274d8ff11be80ce90197b0e60c1c090d73c52ed1b97e0148b947f0d6`, seed `16be9aee01bb9191526eeb7c4603fd3d6762515397760276876d7ea10491b40a`, trueinput `a73fb8a4938c9a490b6ef8dd0284edd5d235336476e5c3ada0e159d899709f04`. Sourceparent bindings совпадают с owning coarse32cases и runtimeinitialization; inheritedsameV/model/frame/units/privateownership checks не менялись. Повторяемый первый boundary уFD сохраняет GR/A, но намеренно меняет LG: одинаковые128map budgets сunmodifiedcontrol не означают одинаковые novel-state counts.

## Источник и однократная FD-пара

Ledgers: A+8ΓR, каждый264011blocks, итого2376099/case; все9 source scans прошли finite/Hermiticity/nativePSDbudget,failedblocks0. Ни один L/G не меняется до завершения этих проверок. **Это не strictPSD certificate:** например retained0V ΓR_totalembedding minλ=−1.33645e−13, admitted своим unchanged numericalbudget; smallnegative tails не clipping/repair. Exact causal model или window accuracy этим не доказаны.

Общий μ дляk/channels каждойcase:0.08482370427731846/0.07804747686020258 native =8.48237043/7.80474769meV. Число независимо пересчитано из сохранённыхq·f, derivative изq·f·hole/kBT: обаfit удовлетворяют target0.045 сrelativeerror−1.54198e−13/−1.15340e−13≤2e−13. Свереныcapacity,stablef/hole, все8R-unchanged flags, GR/A-unchanged иpairchange hashes. Γ использована изi(R−R†), не изΓK; sign conventions L=ifΓ/G=−iholeΓ верны. Embeddingtotal и± не удваиваются вnative total; additive defects до/после записаны.

First fresh GR/A replay reports bitwise zero differences кparent. FreshLG сравнивается сновым ожидаемымfA диагностически, не со старымparentLG. Firstλ=0.9999999999997962/0.999999999999824, firstJ≈8.18e−5/6.80e−5A/m² — **конструируемые seed identities**, неphysical acceptance. Последующие128полныхmaps заново обновляют allmechanisms иshiftedcavityembedding, безFDafterinit.

## Финальные raw witnesses

| Последняя метрика |0mV/period|1mV/period|
|---|---:|---:|
| rK |0.009354297632|0.002129794568|
| rΣ |0.01534610113|0.006508286376|
| rλ |0.000473433815|0.000175533276|
| rD |4.67e−16|2.46e−16|
| rA |7.98e−15|5.32e−15|
| Raw number |0.04502130452|0.04500789900|
| Paired number |0.04500000000|0.04500000000|
| J, A/cm², **UNCONVERGED** |−0.168491259|2.605534432|

Сверены всеλ=rawN/target, rλ=|λ−1|, pairedN,perchannelrΣ иparticle/power scalar balances. Малые rD/rA — current linear-Dyson/spectral algebra residuals, а не сходимость свежейSCBA-карты или forward-error certificate. Pairedtarget создаётся iterationconstraint, не доказывает исходное число. Finalparticle diagnostics0.00026750/0.00009920, power ratios0.57855/0.25317; physical conservation/current acceptance не установлена.

## Что именно сейчас ограничивает rΣ

Coree171 `production_iteration.jl::_production_fixed_point_measurement` берёт **максимум** поlocal5,totalembedding,±embedding иR/</>:

\[
r_\Sigma=\max_{c,s}\frac{\|\Sigma^{fresh}_{c,s}-\Sigma^{current}_{c,s}\|_F}
 {\|\Sigma^{fresh}_{c,s}\|_F+10^{-14}}.
\]

Это **не** pooledglobalratio и не fractionofcharge. F-нормы здесь unweighted поE/k/matrix. Native max независимо воспроизведён изmarkers накаждойmeasuredrow (разница≤6.8e−13). Вfinalbothcases максимум даёт **embedding_plus, lesser**:

| Raw component |0mV|1mV|
|---|---:|---:|
| Difference Frobenius |0.00243838911116|0.00108025541865|
| Fresh component scale |0.158893069334|0.165981543572|
| Их ratio сnativefloor |0.0153461011319|0.00650828637570|

Floor составляет лишь6.29e−14/6.02e−14denominator: это не near-zero denominator artifact. Следующие компоненты — totalembedding</minuslesser, ratio≈0.01504/0.01497 и0.006249/0.006243. Вместе сrK, который проверяет `GL−GR·ΣL_fresh·GA`, это локализует **текущий незамкнутый кинетический/embedding-lesser блок**. Это не причинаего незамыкания: cavitycandidate зависит иотfreshlocalΣ; изmax нельзя заключить, чтоошибочнаименноboundaryphysics, выявленJacobiannullmode илиисключены quadrature/mixing/conditioningeffects. RetardedSCBA не объявляетсярешённой изsmallrD.

## Ресурсы и границы сравнения

Root actual1228.409s/process,1227.938s/script, case659.39/561.24s; cgroupCPU8927.381s,среднее7.2674cores. Peak16765247488B приhard20.75GiB/high20.25GiB/swap0; high/max/OOM events0. Output6850929200B≤8GiB; prospectiveglobal72GiB recordedдоlaunch. Phase log:2cases×128progress, balancednonoverlappingbegin/end, **residualwidth4 во всех256**, другиеheavyphases26;counter_scope=process. Phasewall1087.035s/phaseCPU8605.318s ограниченыwholeprocess/cgroupvalues. Это проверкаисполненияpolicy, не A/Bend-to-end speedup: FDinitialization **и** resourcepolicy отличаютсяотunmodified128.

Приравном128map budget финальные rawдефекты меньшеunmodifiedcontrol, но обеветвипродолжаютнепроходитьnativegates. Не установлены стационарная I–V, математический fixed point, единственность, сходимость по сетке/окну/базису/k/cavity, Poisson, causal real-axis accuracy или экспериментальная валидация. DeclaredfinalstateSHA:0V `b175f676c28f80b330a43c162ee18acdcb68ec5c30c45e92b618b961c47530a9`,1mV `a34db685669bdde89351e7854f1eadd74e8f93a545462a5711cc7384a5922279`; здесьонинеперехешировались. Результатотрицательныйвотношенииприёмкинаданномbudget, невотношениисуществованияфизическогорешения. Новыхопытовилипродолжениянеразрешено.
