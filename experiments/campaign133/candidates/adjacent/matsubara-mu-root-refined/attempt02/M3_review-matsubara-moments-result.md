# Independent actual result: M3 moment postprocess V2

**Вердикт:** completed V2 корректно выполняет заявленный frozen-source
moment postprocess. Нормы/числа подтверждены независимым арифметическим
аудитом saved matrices и исходных kernels. Number-cutoff разность128→256
уменьшилась в128,0055раз после M3, но это **асимптотическое свидетельство**,
не infinite-cutoff/SCBA/native acceptance. Charge control также не принят:
cutoff difference2,03823e−7 target превышает уровень1e−8 в20,38раз,
а при сохранённом μ плотность ещё отличается от target на1,11801e−4.

Primary: `adjacent/matsubara-moments/attempt02/summary.toml`
SHAcb2dcc61903627ecd822f1152ba3d04c1f9679ed6c94344a28e6c95e0774ecc5,
small HDF N128/N256, scalar oracle, archived input/module/adapter и log
`actual-v2-root.log`. Module9182ae80… и adapter7edebeca… совпали с frozen
источниками; actual input отличается от false input только authorization.
Все source packet SHA, включая оба latest HDF, model и small spectral
HDF, проверены. BLAS actual/start/final=1, Julia=1; record inner_maps=0,
mu_fits=0, native_pass=false; функции postprocess не вызывают физический Φ,
Dyson или μ fit.

**Независимое воспроизведение:**
`python3 .build/campaign133/reviews/matsubara_m3_result_audit.py`;
JSON `matsubara-m3-result-audit.json`, точные H5dump commands в
`matsubara-m3-result-h5commands.json`. Один bounded stdlib процесс,
AS≤128MiB/alarm30s; фактически1,04s/41204KiB/exit0. Julia и новые solver
maps не запускались. Читались только small Mats states/density matrices
и пять saved kernel tensors (не Γ/Green real-axis arrays).

Аудит проверяет logical SHA kernels/qK/wk/hH и буквально сворачивает
`Σmp,cc,d wk[mp]·K[m,mp,a,cc,d,b]·X[mp,cc,d]` с qK — без native BLAS
contraction helper. Восстановлены все4 static+LO moments и ordered hS0+S0h:
relative discrepancy S0=2,28774e−16, B0=2,69009e−16.
One-shot из savedρN воспроизведён; direct linear residual вычислен
повторным независимым KLO, не прочитан готовым из TOML:

| Число/норма | N128 | N256 |
|---|---:|---:|
| saved exact-bare-subtracted number |0,04500368792266079|0,045004861993892566|
| one-shot M3 number |0,045005040205109544|0,045005031033066724|
| selfconsistent M3 number |0,045005040205103|0,045005031033066606|
| M3 correction Frobenius |3,63882707020e−7|4,54863794166e−8|
| M3−one-shot Frobenius |1,91330035695e−15|5,94731441423e−17|
| independent linear relative residual |9,14059486e−17|1,26400983e−16|
| independent weighted relative residual |9,81488248e−17|1,43466950e−16|

Числа повторно получены из всех11×5×5 density entries с gs2 и source wk;
нормы полей — unweighted Frobenius, weighted residual использует только wk,
без spin/E weights. Stored feedback точно равен M3−one-shot. Bare ρ в output
совпадает с source latestρ; module не нормализует его и не clip/project.
Числа Python/Julia различаются максимум в последних Double digits из-за
суммирования; JSON сохраняет собственные значения аудита.

**Cutoff interpretation.** Recorded number256−128:
bare+1,174071231777e−6; one-shot−9,17204282e−9;
selfconsistent−9,17203639e−9. Matrix Frobenius difference уменьшается
3,15626716e−7→3,17362966e−9, то есть примерно99,45раз — это иной показатель,
его не следует обозначать number-фактором128. t4/correction при удвоении N
уменьшаются примерно в8раз, согласуясь с удалением leading N⁻³ density tail.
Основной measured эффект уже есть в **one-shot**; дополнительный linear
selfconsistent feedback очень мал, при N256 близок к Double roundoff.
Это не свидетельство необходимости именно новой feedback closure.

η=6,86757e−9/8,58466e−10, cond(L)=1,00000000584/1,00000000073;
linear solve хорошо обусловлен в данном вычислении. Reported conditional
Neumann error bound относится только к уравнению M3 closure при CP/numerical
предпосылках, не к неизвестному finite-N/infinite-N remainder. Pauli minima
M3≈7,24e−8/7,15e−8, maxima≈0,21752; Hermitian defects1,94e−17/2,77e−17
численно благоприятны и не являются exact PSD/total physics certificate.
Scalar auxiliary-pole oracle подтверждает sign/moment/tail в конкретном
аналитическом случае; он не ограничивает remainder реальных пяти каналов.

Finite-N raw Σ=Φ source root остаётся отдельным принятым утверждением о
reference-tail map. Selfconsistent M3 меняет moment feedback и не становится
его точной finite-N density identity. Higher moments, omitted interacting LO
feedback в low n и cutoff nonlinear error всё ещё не ограничены. Fixed μ
8,459614meV не подбирался; ни малая linear residual, ни cutoff agreement не
принимают допантную плотность, Poisson или original raw-h native solver.

**Execution/history.** V2 completed16,1669s script/16,56s process,
722136KiB/swap0/exit0. В log остаются world-age warnings; они не сорвали
этот exact Julia1.13.3 run и не доказывают будущую совместимость wrapper.
V1 exit124/120,12s сохраняется как timeout с partial N128, без финального
summary и без записанного BLAS getter. V1 N128 four reported numbers равны
V2; rounding-level matrix/residual differences после смены BLAS reductions
не выдаются за физическое изменение. Полный speedup из censored V1 не
вычисляется. Explicit BLAS setter/getter устранил наблюдавшийся execution
gap, а не научную ошибку M3.

Новые N512/1024/domain/wrapper и любой μ-charge control требуют собственного
source/parameter/budget qualification. В этом review они не запускались и
не принимаются результатом старой128/256 пары.
