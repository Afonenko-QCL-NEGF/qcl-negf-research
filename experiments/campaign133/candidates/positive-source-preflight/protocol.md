# Однократная квалификация положительного спектрального источника

Вопрос: устраняет ли вычисление `A_cov=GR Γ_source GR†` отрицательные спектральные веса при исходном положительном warm Γ, сохраняя равновесный пятиканальный map? Это новая арифметическая identity; запуск программы не принимает native SCBA. Core e171caa/tree742ed5a неизменен. Root запускает после независимого review; автор Julia/Pkg не запускал.

Основание выбора: affine4096/attempt01 завершил 32 Π без native handoff: r32=5.81008885e−6, best r30=5.80634820e−6. Семь AA trial отвергнуты cone guard, хотя некоторые имели L1<4096; поэтому увеличение affine cap уже не изолирует препятствие. Сам последний checkpoint допускает tiny negative Γ в исходном backward-error allowance и не служит строгим положительным seed. Берётся только исходный warm seed SHA4142e06d…; источники/полные хеши заданы input.toml.

Минимальное изменение и контроль:

1. Собрать исходные S01/m939/NE30704/11k/5b/finite_chain4/70K/fixedU=0 и прежние операторы из замороженного helper cf923bfb…. Проверить фактический Core, модель SHA6b41…, seed SHA4142…, schema, все axes/scales. Γ seed должен иметь положительное измеренное минимальное собственное значение; ремонта нет.
2. Взять **точный Float64 `c0.cavity.hlocal`** в качестве hH и для центральной ячейки, и для восьми соседних. Пересобрать native plan; сохранить raw h, hH, абсолютную/относительную коррекцию. Λraw=PI_FFT[Γ], ΛH=Herm(Λraw); исходный arithmetic target1e−10 не меняется. Сохранить ΛH **до** сложения с −iΓ/2.
3. Один native spectral probe: ΣR=ΛH−iΓ/2, Σ<=iΓ, Σ>=0. Эти L/G — технический спектральный источник, не физические occupations. Cavity R не Herm-проектируется. `Γ_source=Γ−iΣembedding<`; `A_cov=GR Γ_source GR†`. Это обычная double матричная конгруэнция, а не строгий Gram-факторный сертификат PSD.
4. qE=gs/(2π) wE Σk wk ReTr(A_cov); текущий μ заново по **этому** q; G<=if(μ)A_cov, G>=−i(1−f)A_cov; raw Φ=sum static4 K[A_cov]+LO[GFD]. q, μ, FD и Φ используют одно A. Никаких η, clip, normalization или ослабления исходных feasibility guards.
5. Контроль `A_sub=i(GR−GR†)` использует **те же Γ/hH/GR**; собственный текущий μ и тот же Φ только если qsub допустим. Иначе честное not_measured с причиной; qsub не исправляется. Поэтому hH и способ A внутри этого сравнения не смешаны. Старый HermΛ trajectory имел central raw h и не является этим однопараметрическим контролем.
6. Один полный native пятисемейный candidate на переданном GFD даёт независимый consumer control K/KMS каждого канала, сравнение totalΣR с Π[rawΦ]. Scope: supplied stable-A Green; это не fresh unconstrained native Green и не exit audit. Численные нормы/масштабы и pass исходного KMS1e−10 сохраняются, не заменяются успешным exit.

Выход: summary.toml, whole-energy q_cov/q_sub, selected all-k Γ/Λraw/**ΛH_preassembly**/localR/totalR/GR/A_cov/A_sub/Γsource/Γeff, rawΦcov/control, hraw/hH/Tplus и исходные axes/weights. Узлы: minima q/A/Φ, последний E, μ-nearest; данные до известной ошибки сохраняются finally. Выбранные матрицы позволяют существующему full45 Big oracle различить primitive assembly, cavity и inverse; Big qualification не заявляется автоматически. Строгая nonnegativity не следует из generic covariance; если она/точность ledger не проходит, следующий отдельный вариант — structured Gram source и causal Schur recurrence, не repair этого результата.

Конечный бюджет: 1 retarded probe, максимум 2 reduced rawΦ и 1 native candidate, started/completed считаются до/после; 300s, MemoryMax12.125GiB, High11.5GiB, Swap0, CPU400%, J4/BLAS1, heap8G, LimitFSIZE/data32MiB. Сохраняются только selected matrices, полный Green не пишется. runtime_authorized=false в заморозке. Root обеспечивает внешний wall stop и cgroup.

## Следующий 16-map pilot: только протокол

После квалификации оператора первый различающий pilot — **stable-A only**, 16 TOTAL Π including failures/boundary, прежние AA m3/affine4096/cond1e8/α.2/β1, пустая history, тот же строго положительный seed. Затем отдельная identity при необходимости добавляет **только PSD projection AA trial**, с тем же stable-A map/seed и контролем первых одинаковых карт. Не объединять два изменения и не приписывать результат одному.

Для будущего projector: Herm trial→eigen positive part разрешается только как алгоритмическое предложение; raw Φ никогда не проектировать. Сохранять correction abs/weighted, negative part, число изменённых блоков/ranks, raw trial и принятую Γ, residual исходного Φ−Γ. Проекция может создать rank deficiency и стационарность projected step без корня; feasibility и свежий μ обязательны, floor не вводится. Реализация projector здесь отсутствует.

Сходимость reduced Γ по прежнему1e−9 ещё не принимает задачу: прежние materialization/native32/три свежих unprojected passes, rK/rΣ/λ/kinetic spectral/collision/cone/native FDR gates остаются обязательны. Новый A не подменяет native acceptance. Pilot memory12.125/high11.5/file5GiB/data8GiB; precheck180/pilot600/native900/overall1800s, без автоматического увеличения бюджета. Только новый разрешённый finite attempt после qualified operator; текущий preflight ничего из этого не запускает.
