# Независимый review Γ warm initialization: protocol

Проверен `gamma-warm-plan.md` SHA256 `b6a5b6eca6634ec51bc4fc5ca7087cf8da8b41f7fcaa4f830c9f9dfc343cfc3c` до наличия реализации. **Замысел scientific control согласован:** меняется только начальное поле Γ и его provenance, target m939 Π/AA/criteria/caps прежние. Это не source-state exact resume и не перенос native convergence. Frozen code ещё не просмотрен; runtime/import/admissibility не измерены рецензентом. Julia/Pkg/full source H5 не запускались/не читались, проверяемый код не менялся.

Source FD64-continuation128 stateSHA `a690ba56918d4c9803d70d759f2c228f4e9133fdc7f4bcede132f19229c0042f` остаётся unresolved. Main primary seedΓK=Σlocal i(Σ>−Σ<), crosscheckΓR=Σlocal i(ΣR−ΣR†), оба **matrix adjoint** с комплексной ориентацией(a,b); elementwise−2imagR неверен. Embedding/± исключены. Source ΓR/K comparison требует numerator/denominator/zero policy, weights и maxblock; старый near-converged native residual не даёт bound target residual после regrid и смены retarded closure.

Actual-node interpolation с общими nonnegative weights всех matrix entries/k сохраняет convex cone в exact arithmetic, включая zero ghost matrix. Правый targetendpoint выходит заsource на.106496µeV, поэтому declared linear taper к source ghost, а не clamp/constant extrapolation, корректно задаёт рецепт. Targetendpoint и собственные PI ghost hats не менять. Этот initialization не обязан сохранять integralΓ, и нормировкаΓ не разрешена. Допущенная Herm representation correction измеряется отдельно; no eigenvalue clipping/floor. Raw/working cone passage — floating admission, не rigorous PSD proof.

**Проверяемые условия frozen implementation:**

1. Streaming full source SHA и exact source-input/frame guards до использования большихΣ arrays; sourceΓK/ΓR читаются из пяти stored local R/L/G hyperslabs. HDF disk(b,a,k,E) против Julia(E,k,a,b) подтвердить noncommuting complex2×2 fixture и actual axes, не догадкой оtranspose. Kernel/basis/frame comparison источника и actual target ограничить соответствующими собственными measure/scales, не использовать sourceenergy как target grid.
2. Builder выделяет честный отдельный Julia resource window3–4GiB/1thread/300s и bounded chunk; не новая solver map. Old ReΣR/GR/A/μ/embedding/AA history не загружаются в solver state. Target всегда пересчитывает PI localR/cavity/GR/currentμ на warmΓ.
3. Actual source/target nodes и row indices/weights, zero-ghost flags, constant-interior и row-sum/nonnegative witnesses сохраняются. TargetΓraw и Herm correction/seed cone проверены; ΓR−ΓK больше frozen arithmetic target1e−10 означает stop, не выбор удобного поля.
4. Seed SHA+builder/source/target identity+norms/scales сохранены. Solver загружает новый independent dense field, не mutating shared source; AA starts empty. Seed/mapping output учтены в raw data budget. Originalcold attempt02 и source immutable.
5. Current m939 shift/cache не импортируется изsource; initialΠ(warmseed) входит в общий16 и заново проходит source/native-five/KMS/same-input cavity control. Нет claim KMS всех будущихΓ поold input илиinitialpass. Native1+30+1 exit только по прежнему rawtrigger, затемunconstrainednoFDprojection иfreshnativeaudits.

Критерий сравнения cold/warm — raw numerator/scale и weighted rΓ по одинаковым calls1/2/4/8/16, feasibility/fallback/AA counts, currentμ/conditioning/cost, не normalized density и не близостьμ. Сравнение старогоsource иtarget имеет отличияgrid/PI restriction; scientificacceptance изwarmmetadata не следует. При отсутствииtrigger nativeaudit остаётсяnot_measured; отдельное решение о64calls не меняетзаднимчислом cap16.

Verdict protocol: **no blocker**, дождаться exact freeze и проверить initialization delta перед runtime. Source ΓK vs ΓR и matrix orientation — главные functional risks; следующие — boundary taper и нечаянный import old R/history. Реализацию/runtime этим protocol review не принимаем.

## Frozen code delta и bounded builder verdict

Frozen `equilibrium-gamma-warm/SHA256SUMS`: все12entries проверены OK. New `warm_seed.jl` SHA256 `d786df94a73aaeebb1bda1a2633806f3a5fef578e389bd1524d464fd9a048644`, builder `82fece0d77c0b4edd84bb7c64bd46d39625593f970e33f0094d3f0e6ffa59393`, warm checks `7262d8144e9ff9f89d8474a17d795ce2828141f8096d004bf141373b0d96a71b`, functions `99d70fc1e68cfacaa7975d385c12dd1ca1ad375e2ea5874cb224cd9562affe4b`, main `0ff25f63672e98543ce939c89f07deaf14069ae18a02e23d60b5d93836e82108`. Cold79checkscc0e143b неизменны; old cold sources/raw не редактированы. Protocol/header честно обозначают новую warm identity, не прежний baseline.

**Static no blocker. Explicit green только для предложенного bounded builder:** Julia1/BLAS1/heap2G,hard4GiB/high3/Swap0/CPU100/300s,file/data512MiB. Warm16 ещё не принимается: сначала нужен actual builder result/manifest/seed/admissibility review и freeze manifest SHA в authorized input. Runtime/source arrays рецензент не запускал/не читал.

Root `small-root-v1.log`:79cold+19warm=98pass,15.94s process wall,peak694004KiB,exit0,swap0,NoActualS01Maps. Tiny HDF external h5dump fixture проверяет disk(b,a,k,E) against Julia(E,k,a,b), не shared write/read alone. CP interpolation, zero-ghost taper, constant interior, no outside-ghost extrapolation, complex matrix Im vs elementwise Im, model mismatch и frozen option guards покрыты. Полный actual source consumer ими ещё не проверен.

### Прочитанный owning путь

`build_seed.jl` сначала проверяет source/import Core e171/tree742/clean, Julia/resource cgroup/FSIZE; потом streaming fullsourceSHAa690, source-input SHA6784, targetmodelSHA6b41 и targetinputc084. Hashes стоят **до sourceΣ hyperslabs**. Model input guards требуют oldf50/Runner18f45, sourceV0/alloy/finite-chain4/NE30721 и exact targetNE30704/endpoint; source schema/zeroU/E0/L0/spin, k/x/weights, seven basis matrices и actualnormalized kernels+qK сравниваются с frozen targetmodel. Идентификационные TOML copies сохраняются; их substantive contents проверены, fileSHA fullsource/model/input дают основную array provenance.

`build_source_gamma` читает≤1024E R/L/G slabs в logical(E,k,a,b) и суммирует только пять local families в фиксированном порядке. ΓK=i(G−L),ΓR=i(R−R†) используют правильные signs/adjoint indices; source ΓR real part никуда не переносится. Per-channel/total raw norms denominatorΓK, weightedsourcewEwk и maxblock locator сохраняются; превышение1e−10 unweighted либоweighted останавливает builder. Source rawΓK и cones сохранены без clipping. Exact sourceSHA защищает все family shapes; явный logical shape30721×11×5×5 проверен.

`warm_interpolation_plan` требует конечные строго возрастающие actualsource/target energies; row searchsortedlast даёт общие scalar nonnegative weights, включая ghostleft0/ghostright0. На targetexactsourceN весghost0, на targetendpoint заN taper кzero; выход заadjacent ghosts запрещён. Chunk regrid читает только необходимый диапазон sourceΓK, затем создаёт **раздельно** rawY иZ=HermY. Raw/working cones и correction norms сохраняются. Нормировка integralΓ, eigen-clipping и additive floor отсутствуют. Ошибка в данных в пределах native allowance остаётся floating admission, не formal cone certificate.

Seed manifest связывает source/model/input/builder-code/Core identities и actualseedSHA. Loader проверяет manifestSHA до parse, executioncomplete/Core/provenance, coldinputSHA и все прежние cold options кромеidentity/runtime. SeedSHA и modelSHA проверяются streaming. Потом `verify_problem_model` сверяет **фактически заново построенный** target problem: всеgrids/profiles, basisoperators, projectedh, kernels/qK; seed axes/scales/spin/shape проверяются перед чтением одного fresh dense Γ field и полнымcone. Packed targettrustedmask сопоставляется denseInt8. Это защищает regrid кстарой/не той grid и использование sourceR вместо текущегоPI. Дополнительная exact equality может честно остановить run при недетерминированном rebuild; обхода такого отказа нет.

`reduced_stage` заменяет только constant initializer наload_warm_seed. После этого unchanged initialΠ строит новый PI R/cavity/GR/A/currentμ и nativefive/KMS/sameinput control; states/residual history создаются пустыми после initialevaluation. OldG/μ/embedding/history не импортируются. Native exit/strict criteria,onefallback/caps16/32 прежние; stdout-only fix теперь сообщает acceptedconvexfallback честно. Default source fields/input scheme сравниваются с cold даже для numerically irrelevant seed_halfwidth metadata.

### Что ещё требуется до warm16

Actual builder должен завершиться вcap, предъявить fullsource hash, per-channel+aggregate ΓR/K ratios/scales, source/target cones, row/endpoint weights, rawHerm correction, SHA seed/manifest и закрытые HDF outputs. При недостигнутом guard это negative builder outcome, не permission ослабить1e−10 либо выбратьΓR. Seed должен войти в globaldata ledger; builder output хранится один раз, pilot manifest references егоSHA. Внешний cap обязателен: internal checks междуslabs не останавливают один HDF/SHA call и не гарантируютfinally послеOOM/SIGKILL. Wholebuilder peak/time покаnot_measured. Никаких actualwarmΠ/KMS/nativeconvergence claims до следующего stage нет.

## Actual builder result: seed admissible, warm16 green

**Explicit green для ровно warm16 с прежними Φ/AA/criteria/caps и новым initial current-μ/native gate.** Результат builder закрыт и пригоден как initializer; это не warm Π/native/science acceptance. Authorized pilot input обязан зафиксировать manifest SHA ниже; обычный loader ещё раз проверит seed/actual rebuilt model до первой Π. Выход за16, ослабление guards или импорт старого G этим review не разрешаются.

`seed/builder-summary.toml` executioncomplete/smallpass; `builder-attempt01.log` exit0,27.78s process wall,23.013s harness,peak766824KiB,swap0, data412207360B<512MiB. Builder source e171/tree742/import path записан; build guard доΣ import прошёл streaming source a690/input6784/model6b41/targetinputc084. Рецензент исходный3.8GB state повторно не читал; fixed builder code и closed completion дают свидетельство исполнения этих guards. Source/target k/x/weights, basis7 и kernels5/qK имеют сохранённые successful exact fingerprints; сам source остаётся old f50 unresolved.

Независимый streaming hash малых builder artifacts совпал:

- manifest `18439792c37a85685831ee5037d624760b764b3ec9babf42b3c5811a0155a0be`;
- gamma-seed `4142e06db48b20772d99332f277baf3ecb71bb3d79f8bdfc5951c98f807e42e8`;
- source-gamma `aa0b317629ea0947042de2604d369d745cd81216718cf168eb41a6ad21013d80`.

Source ΓK−ΓR global relatives3.768e−16–4.435e−16, weighted3.768e−16–4.257e−16; aggregate3.83917e−16/weighted3.80951e−16. Numerators/scales/maxblock/locators сохранены;1e−10 target пройден без изменения формулы. Пятёрка channel IDs и embedding exclusion соответствуют frozen recipe. Source cone IFR имеет tiny negative eigenvalue−1.0967e−18 при allowance7.1054e−14; это ожидаемый **admitted backward-error** case, не strict exact PSD theorem. Aggregate source min≈1.262e−27, targetraw/working min≈1.25853e−27 также намного меньше allowance: rigorous PD/dark-state exclusion отсюда не следуют.

Raw→Herm target conversion absolute4.08356e−15/scale17.28219, relative2.36287e−16, weighted2.34664e−16. Raw anti-Hermitian max≈3.78e−16, working0; no eigen-clipping. SourceΓ и raw target import сохранены отдельно, поэтому correction не скрывает исходный state.

Tiny reproducible read-only audit сохранён `reviews/warm_seed_result_audit.py` + `warm-seed-result-audit.json` с точными h5dump paths/selections. Проверены только axis/row vectors и три selected5×5 slabs, не full Γ/source state:0.514s,Python23432KiB. Первый запуск самого review-script имел неверный порядок h5dump−s после global options и остановился без научного вывода; исправлен только review-script, последующий audit завершён. Все30704 interpolation rows независимо пересчитаны Pythonbisect по actualnodes: indices/ghost flags совпали, maxweightdifference0,0≤w≤1, ровноoneghostrow. Endpointsource[−2.9999999999999996,9],targetlast9.000001064962726; rightghost9.000390625. Последний sourceweight=.9972736954214245,ghost=.002726304578575515.

Selected rows e1/k1,e15352/k5,e30704/k11 независимо воспроизведены scalar convex combination builder sourceΓK, затем matrix Hermitian conversion. Raw interpolation и working Herm differences0; interior indices15360/15361,rightweight.50108626198, last30721/30722zeroGhost. Это actual orientation/boundary witness, не all-entry independent ΓR/K reconstruction. Whole-field finite/cone/identity checks обеспечиваются прошедшим immutable builder и повторным loadercone, не подменяются тремя выборками.

Сама warm Π ещё **not_measured**. Следующий outcome должен повторно предъявить actual rebuilt-target fingerprints, source/manifest/seed binding, seedcone/currentμ/initial KMS и trajectoryraw residuals. Density-fit и resemblance со старымsource ни сходимость, ни physical acceptance заранее не устанавливают.
