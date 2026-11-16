# Matsubara ALL5: независимая проверка результата

**Решение: приняты два finite-N fixed-μ fixed point нового direct-static+LO/reference-tail варианта.** Raw residual на N128/N256 достигает точного floating fixed point; это подтверждено binary HDF5, а не округлением scalar log. **Плотность по cutoff не квалифицирована:** N128→N256 меняет number на1.17407123e−6≈2.60905e−5 dopant target. Infinite-cutoff/native real-energy/density-fit/Poisson/experimental acceptance не установлены.

Actual `adjacent/matsubara-all5/attempt01`, summaryb9397bf71b7e483bfc3062c0c002a78612c9e896ad1ee5a0b1c443e55113131e; module1434e86f…/adapter94c6f2d4…/frozeninputd920768e…/actualinput748d6093…. Superseded adapter9b18 архивирован **без actual numerical исхода**. Runtime input меняет толькоauthorized bool. Accepted Coree171/tree742/Julia1.13.3; sourceframe verified model700ae…/spectral42fe…/hH9748…/Tplus466b…/nativewk b6de…. Фиксированная μ=.08459614262982018native=.008459614262982019eV; E0=.1eV, β16.57788303078583native−1, Ω=.367native, V0/70K/Nb5/Nk11/finite4 all9hH/thermal bath. Source Γ/Φ arrays не импортируются, CP/PD старой Γ не принимает новую closure.

## Что решено и что проверено

Code/math просмотрены независимо и сопоставлены с `reviews/matsubara-all5-static.md` (включая final94c6 review). StaticΣs=Ks[G] использует те же raw5 kernels/qK/wk и **direct complex G**, поэтому это новая full-line identity относительно finite-energy HermPI исходного пилота. LO имеет правильный положительный1/β bosonic kernel2Ω/(ν²+Ω²), negative fermion packing G_(−n−1)=G_n†, analytic bare45 pole reference и finite convolution толькоΔG. В full45 одна central-driven Σ повторяется all9cells; это не cellwise SCBA и не infinite lead. No factor9/spin/2π/μ double shift. Densityρ=ρref+(2/β)ΣHermΔG и number=gsΣwkTrρ, дополнительного energy1/(2π) здесь нет.

Два независимых стартаΣ0, N128/N256 по16 α1 maps; requested=completed32, statuses iteration_budget сохранены, native/scientific_pass=false, sourceΓused=false, density_normalized=false. Causal guard/checkpoint/callback до update и немутирующий S обеспечивают aligned inputΣ/G/Φ/ρ. Request ledgers ограничивают16/32, нет скрытых дополнительных inner steps. μ-fit/Poisson/normalization не вызываются.

Independent binary H5dump подтвердил **все сохранённые real/imag bytes latestΣ==latestΦ** на обоихN, ||Φ−Σ||=0. По32 small records нулевые raw absolute+relative начинаются сi11 и остаются доi16. Дополнительно дляi1..15 каждогоN восстановлены raw norm/relative из successive savedΣ: приα1 Φ_i=Σ_(i+1); все30 replay совпали (без нового Φ/Julia). Начальные N128relative[1,.00330423,3.48158e−5,3.38617e−7,3.96039e−9,4.38822e−11,…], N256 аналогичны. Exact0 означает floating fixed point **этого finite map**, не mathematically exact infinite solution.

All5 actual17 non-Hermitian-slice native/six-index errors≤1.43252e−15. Independent onepole raw/exact residual norms по complex records и refinement перепроверены; actual module two-pole wiring controls Grelative≈1.54e−16/referenceLO0, ΔLO abs3.061e−9→4.005e−10/density2.871e−9→3.589e−10. Отдельный Root tiny22/22PASS сохраняется primary qualification native primitive wiring; reviewer его не повторял.

## Плотность и causality

Все32 measured Mats nodes diagnostics дают ImΣmax<0/ImGmax<0 и Dyson backward<4e−17; это finite-node floating causality/linear residual, не real-axis analytic continuation certificate. CP ratio L/ω0²=.30529308144788564 независимо пересчитан из L/.189505056…², но **условен на complete positivity kernels**, не strict exact CP certificate и не finite-truncation/infinite-SCBA uniqueness proof. Малый conditional bound не равен измеренной скорости nonlinear iteration.

Из rawρ и nativewk независимо получены N128=.045003687922660796, N256=.04500486199389256; difference−1.174071231763163e−6 совпал с report. rho Frobenius difference3.156267158674883e−7. Finalρ Hermitian defects≤6.94e−18, floating complex LDL дляρ иI−ρ положителен во всех11k обоихN (не exact certificate). Reported конечный Pauli spectrum: N128min8.25463e−9/max.2175201674, N256min6.36194e−8/max.2175201827.

**Промежуточная Pauli оговорка:** N128-i02ρmin=−3.191154307277137e−9. Независимый rawρ floating LDL/Rayleigh witness приk11 даёт normalized v†ρv=−2.9378378645433084e−9 (imag≈6.6e−34); в audit сохранён сам vector. Начальныеi01 обоихN имеют ещё≈−9.6e−15 roundoff-scale minima. Поэтому нельзя утверждать Pauli positivity всех iterations. Никакойρ repair не было; finite reference-tail reconstruction при intermediateΣ не обязан быть physical full-frequency density. Это не опровергает положительные конечныеρ, но ограничивает общий CP/causal-density вывод.

Shared low32 raw Grelative1.5166459534983262e−13/Srelative8.053484047706052e−11 полностью перепроверены из сохранённых HDF; nonlinear residual обоихN равен0, поэтому здесь cutoff witness не спутан с незавершённым inner iteration. Однако его малость не контролирует integrated density tail: number difference выше независимо подтверждён. N128/N256 отклоняются от.045 соответственно+8.19538e−5/+1.08044e−4 relative; **это fixedμ outputs**, μ root не решался. Infinity extrapolation/tail completion не выполнялись reviewer и остаются отдельным этапом.

## Стоимость, воспроизводимость, границы

Actual logexit0,25.12s process/24.940049s script, MaxRSS888616KiB=909,942,784B,Swap0; bytes60,325,650<64MiB. Root scope2GiB/High1.75GiB/J1/BLAS1/CPU≤200%/600s/file32MiB. Adapter проверяет memory.max/swap/CPU/file, High задаётся externally (прежняя caveat сохранена). Source build_context со старым native memory estimate — model reconstruction, не native workload admission под2GiB.

Независимый bounded audit `matsubara_all5_result_audit.py`/`matsubara-all5-result-audit.json`:26checks PASS, stdlibPython≤128MiB+H5dump только bounded result states и11wk scalars; новая Φ/Julia/Pkg/solver не вызывались. Source model/Γ HDF5 большие поля повторно не загружались. Подтверждены finiteN raw equation fixed points, source/frame и tail/ρ limitation; nativereal-energy root/dopant fit/∞cutoff не приняты. Следующее различающее свидетельство — отдельный tail analysis или другой заявленный cutoff test, а не добавление итераций уже floating-converged finite maps.
