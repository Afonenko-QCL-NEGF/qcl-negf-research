# Независимый static review affine4096 pilot32

**Static green для одного bounded exploratory pilot32 после new tiny pass и штатного resource admission.** Единственное изменение численного выбора AA — предел affine L1:10→4096. Φ, PI/HermΛ, μ-fit, physical/domain/native tolerances не изменены. Новый предел не является доказанным bound ошибки Π или гарантией сходимости. Рецензент Julia/Pkg/maps и больших parent arrays не запускал/не читал.

Freeze `equilibrium-gamma-affine4096/SHA256SUMS`:15/15OK. Functions SHA256 `cf923bfb7e5dc56d5aa296c6863e8615a92984f82e125433ea1002bac9c9ee63`; restart `b020d599a5eb87c1dc658d703b30acac93fd4524da86764fef4c3a032534229e`; main `876a8c3e4fa1c2c76a7962c5f4b16ecc47bfe852946c67aaf1502b4531228b49`; input `b5d09e3e8f2cf91016cd0ec9c8eea5b9083249be730d77c9590ec09bdbe3236a`; patch `d43f75989a48ab9b7a75683c8390219b708bbcad68f9b01dae3507dca1393397`.

## Source и initializer

Pinned source — **завершённый** cont128, не live/current selector. Малые parent input fddde931, identity c63682ba, map1288a34e3ac и operator57f2 hashes independently совпали input. Parent summarycomplete/128calls/pilot_call_cap,finalrawrΓ7.19230522e−6,nativenot_measured. Большие Γlastb30faef1/terminalfb7cea34/model47682acc hashes будет streaming проверять loader; рецензент их повторно не читал.

Restart contract разрешает только новые lineage metadata/run budget/affine cap; independent scalar diff не нашёл других inherited option changes. Explicitparent_call128 проверяется вgate/loader/HDFcall. Initializer читает **Γ input map128**, без repair/interpolation/rescaling, старыхG/μ/history. Current rebuilt frame comparison и raw cone прежние. First Π call1 повторяет rawΦ128/scalars по256eps и nativefive/KMS/cavity controls до continued AA. Это fresh empty-history run с максимум31 дальнейшей Π, неexactAAresume.

## AA delta и telemetry

Просмотрен полный AA diff. Порядок existing weighted QR projections/R solve/affine coefficients/trial updates сохранён. Histm3,QRcond≤1e8,β1,linearα.2 прежние. Новая finite-coefficient проверка только **report flag**; scientific/numerical guard по-прежнему existing finiteL1+limit. Для того жеcap10 proposed candidate и все old info fields должны быть bitwise unchanged.

Telemetry сохраняет каждый attempted suffix в Vector{Dict{String,Any}}: retained firstpair/differences/normalized rank/Rdiag/initial scales, condition с measured flag, coefficient-finite flag, affine coefficients/L1/reason. Rank/condition/nonfinite/L1 причины не смешиваются. Никакой скрытой дополнительной Π или new Field reduction нет: m1d использует уже вычисленный QR initialnorm, m1a/b — cached existing working weighted_absolute map scalars. Cache push/pop/clear синхронизирован с states/residuals, включая fallback reset. Эти cached norms могут отличаться от inner-based norms arithmetic rounding; они не управляют candidate.

Научный смысл L1≤4096: условная чувствительность affine combination к **уже ограниченным** ошибкам прошлых pairs. Вычисленный Cauchy m1 bound≈3722.76 мотивирует4096 для измеренной old parent pair; он не гарантирует следующий fresh-history pair и не ограничивает nonlinear Π noise. FFTroundoff≈3e−16 не охватывает μ/cavity/contractions/condition effects, поэтому4096×FFTroundoff не является целой error budget. Более крупный extrapolation может не пройти cone/q/domain; прежний отказ/один convex fallback/caps должны остаться решающими.

## Tests и границы

`aa_reference.jl` действительно exact parent AA function body после переименования: это независимо проверено текстовым сравнением. Tiny reference защищает механическую telemetry-equivalence cap10, а не независимую correctness AA math. Новые tests также различают cap10 reject и512/4096 admit на scalar contractionλ=.9 (L1≈99), проверяют same tiny candidate bits для512/4096, nested TOML producer/consumer, parent128 contract и budget32. Toy convergence не доказывает QCL convergence.

Root `small-root-v1.log` затем завершён и проверен:24/24pass,5.06s process wall,peak502824KiB,swap0,exit0,NoS01Maps. Tiny prerequisite выполнен. Physical104/restart15 evidence можно использовать по unchanged bodies, без нового научного окружения или повторного большого расчёта.

Fixed run envelope:32TOTAL including boundary/rejected evaluated maps,precheck180/pilot600/native900/overall1800s; native1+30+1≤32 только прежнимtrigger. Hard12.125GiB/high11.5/heap8G/J4/BLAS1/Swap0,file5GiB,data8GiB; padded-HDF admission unchanged. Global32GiB admission и сохранность priorraw остаются ответственностью Root wrapper; массивы не удаляются. Main binds32/4096/parent128 и same source/resourcegate. No automatic extension.

Decision: no static blocker, tiny pass подтверждён; **explicit green для одного bounded affine4096 pilot32** при штатном resource admission. Replay/pass/actualsuffixguards и resulting raw trajectory нужно измерить до любых выводов. Ни large coefficients сами по себе, ни consumer/replay pass, ни normalized density не дают restricted/native/Poisson/discretization/experiment acceptance.
