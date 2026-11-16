**Независимый static review stable-A-only16**

Решение: **кандидат одобряется для root tiny integration check и затем
одного bounded pilot16 после qualifier/result/resource admission**.
Блокирующих math/sign/API/graph ошибок в frozen коде не найдено.
Это experimental reduced arithmetic identity; согласие reviewer не
заменяет actual run, selected Big-oracle accuracy или untouched native
acceptance. Автор review Julia/Pkg/HDF5/solver не запускал, исходники
кандидата не менял. Выполнены чтение кода и механическое сравнение текстов
и малых input TOML.

Frozen package `equilibrium-gamma-stableA16/`:

| Файл | SHA256 |
|---|---|
| gamma_functions.jl | ffde0384a99887669ccaa729428e0eb1f885155925c6db465dcc36f8ca5b6a74 |
| stable_spectral.jl | eeca1f29710d9bb60bf6ce70ce0ebd1620fb86d4c77314b4bd8ffa9a5bcd5aa5 |
| run_gamma.jl | 28e505f046ca8c97d8a9a447421d218247b4761e492f238016fa96d3244f577d |
| stable_checks.jl | 7d7d33bbd9f4c934336322d6b87c12d8099aeafbe2fb25e518dce522d246c008 |
| pilot16.toml | 3e58c1c3df1b4afd09656b6c6b6b60d41a34b6be0c2cd66693a113f4098fe463 |

Сопоставлен actual parent cf923bfb… из affine4096/attempt01, primary
preflight a34b1da… и его result06d6fa46…; последний независимо рассмотрен
в `positive-source-preflight-result.md`. Source remains e171/tree742ed5a.
Runtime false / qualifier SHA UNSET в frozen input — ожидаемая граница:
root заполняет exact approved summary digest и authorization в отдельной
runtime копии с новым SHA, не меняя seed/physics/algorithm controls.

**1. Минимальная дельта и входы**

Механическое сравнение подтверждает byte-identical AA `aa_trial` body
(9330chars) и `native_stage` (3453chars). `audit_native` меняет только
scope annotation; public native calls/residual gate arithmetic не изменены.
`covariance_A` и `number_weights_A` совпадают byte-for-byte с квалифицированным
preflight. Их generic F64 congruence не назван strict Gram/PSD certificate.

Actual cold target-input comparison даёт **ровно** четыре inherited
различия: identity, runtime authorization, affine L1 10→4096 и data cap
6→8GiB. Warm loader whitelist разрешает именно эти controls; все прочие
cold inherited keys, включая physics/grids/budgets, проверяются. Extra
spectral-operator/native-scope keys объявляют новую identity, не исключают
старые physical keys из проверки. Это не попытка replay позднего Γlast128:
инициализация original warm Γ4142e06d… с empty history через manifest184397…,
source FD128 a690… и model6b41…. SHA/schema/axes/scales/operator verification
в `load_warm_seed(c,...)` выполняются **на original raw-h Context**.

После loader дополнительно требуется положительный **computed** seed
eigenvalue. Это floating witness, а не строгая рациональная PSD-квалификация.
Все последующие trials остаются под прежним cone backward-error allowance;
они могут иметь tiny negatives. Поэтому новое название stable-A не обещает
положительность каждого последующего входа/Φ и не создаёт скрытый Γ repair.

**2. Единый reduced оператор**

`stable_context(c)` копирует exact `c.cavity.hlocal`, строит новый plan
с этим hH и проверяет `isequal(newplan.hlocal,hH)`. Original c.h/c.cavity
не изменяются. Новая central Dyson и все8 side cells имеют одну hH политику,
а native plan stale check получает именно passed hH.

`_stable_spectral_state` выполняет native PI Λraw→explicit HermΛH с
прежним correction1e−10 и сохранением preassembly buffer. ΣR=ΛH−iΓ/2,
Σ<=iΓ,Σ>=0 — spectral f=1 source, не physical occupancy. Native retarded
recursion неизменна; physical total source Γ−iΣemb< имеет правильный sign.
GR содержит no extra η. A_cov=GR physical GR†, Asub — тот же computed GR;
effective source `i[(hH+ΣR)−(hH+ΣR)†]` позволяет отдельно увидеть mismatch.
Preassembly ΛH не восстанавливается назад из assembled ΣR.

q вычисляется только из A_cov с gs·wE·wk/(2π); каждый entered Π заново
fit current μ, затем **тот же A_cov** идёт в FD lesser/greater, static4
и LO raw Φ. Старый GR-diagonal number-weight extractor в reduced evaluate
больше не вызывается. Native FD pair identity/condition/scale соответствуют
этому GR; нет density/A/Φ rescaling, floor, clipping или PSD projector.

Exact source-congruence proof из положительного operator plan применим
к Hermitian all9 model. Floating Γeff/source mismatch и forward accuracy
не доказаны одним положительным q: selected Big split-input oracle остаётся
отдельным свидетельством. В spectral ledger это честно сформулировано.

**3. Невязки, failures и конечные вызовы**

`_evaluate` сначала input cone, затем `actual_call!(pi=true)`. Ошибка cone
до входа не расходует Π; inverse/A/μ/rawΦ/control failure после входа
расходует его. Новый wrapper сохраняет separate `domain-failure-NNN.toml`
с before/after consumption и reason, затем rethrows. Programming/budget
exceptions не преобразуются в Infeasible/fallback.

Raw Φ проверяется и при failure сохраняется до representation conversion.
`gamma_residuals` отдельно возвращает исходный raw Φ−Γ и Herm working
residual. Trigger, streak×3, best-score и reported residual используют
**raw** относительную/weighted норму с denominator ||raw Φ||, без additive
convergence floor. Малый смешанный/AA step не заменяет этот критерий.
Herm conversion для AA явно ledgered; PSD projection отсутствует.

AAII m3/L14096/cond1e8/α.2/β1, affine admission и один convex fallback
не менялись. После entered failure fallback расходует дополнительный Π;
верхний предел16 TOTAL проверяет `actual_call!`, включая failures. Trial
cone rejection и AA inadmissibility сами не выданы за вычисленную карту.
Reduced-stage finally сохраняет last valid/best Γ, не заменяет их failed
trial. Счётчики и domain records позволяют восстановить actual work.

Spectral selected ledger/H5 пишутся до μ и последующих domain failures:
q_cov/sub, minima/masses, h/rawh/T, Γ, Λraw/preassembly ΛH, assembled local/
totalR, GR, both A, physical/effective sources и source weights. Это не
добавляет maps и не пишет full Green на каждом шаге. Ранний inverse failure
без Green честно имеет только failure record. Output/data accounting
переиспользует штатный padded-HDF admission.

**4. Native выход не подменён reduced hH map**

На initial same-G consumer candidate per-channel KMS/retarded checks
сохранены; initial GR control теперь **matched hH all9** и именно так
помечен. Его pass не назван original raw-h graph equivalence.

После trigger reduced history/QR/старые Green освобождаются. `native_warm`
получает stable-A FD только для материализации scattering, затем вызывает
**original c.h/c.cavity raw-h** embedding и currentGreen. Public native
SCBA и fresh audit работают на original raw-h Context, старых физических
convergence/collision/FDR/cone gates. Новый hH Context наружу не возвращён
и не заменяет native state; matched precheck не принимает task.

Native structural accounting1+30+1≤32 сохранено. Initial stable response
в materialization является retarded construction для одного scattering
map, не лишним hidden Π pilot. Scientific_acceptance не устанавливается
автоматически reduced trigger; исходные native результаты подлежат
независимому result review. Даже успех нового reduced map не переименует
old raw-h failure в pass.

**5. Содержательность tiny integration test и ресурс**

Новый `stable_checks.jl` строит actual configured toy65×3×2 finite4, вносит
малый anti-Hermitian perturbation rawh, проверяет неизменность original
Context/plan и ненулевую h-correction. Один 2-band centre проверяется
независимой assembled18×18 inversion и explicit double full-source sum,
не native cavity recursion. Same-A q→μ→FD сверяется с native number
functional, затем raw Φ сравнивается со всеми5 native kinetic families.
Ledger проходит TOML roundtrip. Это integration behavior, а не assertions
на случайную внутреннюю версию. Сильный/weak complex noncommuting primitive
Big192 tiny8 был квалифицирован ранее exact copied function bytes;
новый tiny не выдаётся за повтор всего physical104 suite.

Reviewer tiny/API execution не запускал; после static pass root выполнил
новый integration tiny, и его `equilibrium-gamma-stableA16/tiny-root.log`
прочитан: **11/11 pass**, exit0, wall13.02s, peak771648KiB, Swap0,
J4/BLAS1. Frozen ffde/eeca/28e source hashes при этом не изменились.
Это закрывает tiny integration gate; actual S01/time/RSS остаются
not_measured до pilot. **Итоговый static+tiny gate: APPROVE actual16**
после root source/qualifier/resource admission без дополнительной кодовой
дельты или повторения tiny. Big-oracle qualification может идти отдельно;
результат pilot остаётся provisional arithmetic experiment, не принятый
SCBA. Большие source/model HDF5 reviewer не читал.
Capsule freezes finite proposal: hard12.125GiB/high11.5GiB,
Swap0/J4BLAS1/CPU400%, pre180/pilot600/native900/overall1800s,
16 Π/native32, data8GiB/file5GiB. Source gate сравнивает active cgroup
с input, outer timeout/admission остаются у root. Его legacy maximum
memory bound12.375GiB шире текущего12.125GiB input; это не расширяет
полномочия: reviewed runtime copy может менять только false/qualifier SHA,
а фактический hard cap должен остаться13019119616B.

Против старого warm Herm16 новая траектория объединяет hH/A_cov и AA
cap10→4096; изоляция эффекта только A была лишь в matched preflight.
Будущий отдельный PSD-trial experiment обязан сохранить этот operator,
seed/controls/empty history и raw Φ acceptance. В текущем candidate он
не реализован. Одобрение ограничено этим заявленным pilot и не утверждает
нелинейное улучшение, discretization sufficiency или scientific acceptance.
