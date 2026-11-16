# Refined μ-root N512/N1024: отдельная identity, WRITE-ONLY/HOLD

В actual M3 V2 при sourceμ число изменилось N256−N128 на
−9.172036394278216e−9 native =−2.03823e−7 target: M3 улучшил bare
разность1.1740712317770408e−6, но не достиг1e−8target. Поэтому прежний
μ128/256 вариант остаётся UNRUN, его source/input сохранены отдельно.
Здесь N512outer/N1024witness — заранее объявленная более строгая проверка
при прежнем cutoff tolerance, не автоматическое продление старого опыта.

Identity: fixed charge, V0/U0, all5 direct-static/thermal-LO, Hermitian
finite45 с hH во всех9cells, full-line bare reference/finite-N LO
reference completion и отдельная **linear M3 asymptotic density closure**.
Никаких Γ-loads/Poisson/native-pass/ρ normalization или clipping.
Это не infinite-N SCBA certificate и не certified tail error.

Исходный module1434/adapter94c6/source frame неизменны. Их source/Core/
project/units guards проверяются wrapper; original gate workload не вызывается.
Original load_modules/source загружают context при sourceμ
0.08459614262982018 native. Новая копия MatsubaraAll5Refined отличается
ровно module name, pilot allowed N=(512,1024), domain comment; все numerical
function bodies byte-identical после обратной замены. Pilot/controls берутся
из refined module, MM.correct использует original applyK — тот же body.
При μ меняются только μ/atoms.eps−Δμ; fullH,K,wk,beta,Omega,V shared.

1. До runtime — frozen M3 API+independent result review+completed actual
   qualification, новый wrapper/module static review и Root source-free tiny.
   Semantic gate связывает completed M3 report с module/adapter/source/input,
   BLAS1 и finite positive linear qualification; иначе HOLD без fallback.
2. N512 endpoints sourceμ±0.001 native, затем safeguarded secant/bisection
   строго внутри sign bracket, максимум8 calls включая endpoints. Bracket
   не расширять, monotonicity не предполагать. Каждый call:Σ0,α1,16 maps.
3. Inner qualification: status=iteration_budget, last3 raw REL≤1e−9 **и**
   ABS≤1e−12 nativeΣ. Норма unweighted F по uniform n,k,a,b; wk raw norm
   отдельно diagnostic, без wE/spin/floor. До qualification нет N(μ).
4. MM.correct сохраняет raw/one-shot/linearρ, S0/B0, t4, Pauli/Herm,
   residual/scale/η/condition/feedback. Candidate data записываются до gates.
   Target0.045000000000000005; outer linear-M3 charge REL≤1e−9.
5. Только qualified root даёт один независимый N1024 при том жеμ/Σ0/α1/16.
   Требовать те же inner guards, witness charge REL≤1e−8 и
   |N1024−N512|/target≤1e−8. Это **number** cutoff witness, не matrixρ bound.
   Нет дополнительных N/calls/tolerance relaxation; failure сохраняется.

Envelope:≤9×16=144maps;2GiB/swap0/1200s/Julia1/BLAS1/OMP1/CPU≤200%,
file32MiB, output1GiB. Старый512MiB cap повышен **до** этого нового опыта:
8×16 raw N512 states≈576MB плюс N1024≈144MB, latest/maps/metadata дают
ориентир<~850MB. Одно N1024 complex field4,505,600bytes, LO D1024×2048
Float64≈16.8MB. Source actualRSS~0.7–0.9GB+small fields даёт ориентир<1.1GB;
это memory estimate, не measurement будущего процесса. Root hard caps
обязательны, wrapper проверяет wall/output; превышение не разрешает рост.

Каждый call-start/reserved16 пишется до pilot. Каждый map сохраняет aligned
Σ/G/rawρ и scalar witnesses, latest добавляетΦ. maps_completed_callbacks —
честный callback counter; structural upper bound144 отдельно. Состояния между
μ не перезаписываются. M3 exceptions остаются not_measured; transient rawρ
Pauli violations не проектируются и не объявляются final pass.
