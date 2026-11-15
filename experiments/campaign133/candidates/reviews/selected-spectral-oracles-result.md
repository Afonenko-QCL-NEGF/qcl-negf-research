# Независимый review selected spectral oracles V2

**192-bit oracles подтверждают отрицательный summed spectral trace уже для frozen exact-double primitives raw retarded operator на двух проверенных энергиях.** HermΛ устраняет sign obstruction этих q, Hermh-only — нет. Однако HermΛ-only сохраняет matrix-indefinite A; отдельная Float64 inverse/subtraction ошибка также существует. Это локализация двух разных numerical effects, не universal fix, не Core acceptance и не SCBA result.

## Primary и исполнение

`warm-spectral-obstruction/attempt01`: script SHA256 `08e29ff1fc3be2c31345cb5795bea80c0602bd6afab883eb025010d842fbe602`, summary `20ae721f53e1cf1d3942cccb902f7e3d602c5512ed19a6b1b9972ed1a0d3fc39`, witnessesH5 `aaf51d926245cd34af10d385003ce73b35e1d983c0081347f80acc9a79684509`, tracesCSV `55b13b820902ae85cd5213446650d9b516ff8c573be516909ac2878e6cee8d61`. Inputs hash selected45fa3f… и coldmodel6b41… проверяются до вычисления. NoCore import: JuliaLinearAlgebra/HDF5 и сохранённые primitive matrices; no fresh full real-axis pass,μfit,Φ/SCBA.

V2 tiny9/9pass иmain9/9pass, exit0,12.58s processwall,6.572s harness,peak642976KiB,swap0. V1 repeatedrawh не применим кnative attribution, сохранён вinvalid-repeated-h-v1; его результат не смешивается сV2.

Рецензент не запускалJulia/новыеoracle solves. `reviews/selected_oracles_audit.py` +`selected-oracles-audit.json`: read-only stdlibDecimal90/малыеHDF,2.885s script/3.19s process,peak24236KiB. Пересчитаны96CSVtrace rows, highprecision weighted aggregates иselected downcast central-row physical products.

## Какой operator действительно проверен

NativeCoree171 `production_operators.jl:309–315` строит sidehlocal как **Float64(h+h†)/2**; centralDyson используетrawh. Sidecouplings370/384 — adjoint(Tplus) слева иTplus справа. V2 full45 chain содержит8sideHermcells иcentralrawh, uniformupperTplus/lowerTplus†. ActualTminus−Tplus†=0 сохранён отдельным diagnostic, не использован дляусредненияcouplings.

ВBig192 sideh сначала формируется вFloat64, затем конвертируется какточныйdouble, так чтоprecision variant не подменяет native sideprojection более точным h. Energy/rawh/localΣR/Γ/T — exact savedFloat64 primitives. **D затем собирается вBigFloat**, ане копируется уже собранныйFloat64D. Поэтому Float-vs-Bigdifference включаетassembly rounding+LU+spectral subtraction; он не изолируетLU единственным фактором. Нельзя назвать rawBig idealcontinuumPI oracle: upstream Float64 Hilbert/source matrices сохранены.

FiveRHS решения D†Y=P, row=Y† дают именноcentral **row** полногоGR, правильную дляΣ_j row_j Γ_j row_j†. D\P дал быcolumn иошибочный physicalproduct дляnon-Hermitianoperator. Per-cellΓeff_j=i[(onsite_j+ΣRlocal)−h.c.] включает central/sideh различие. Γphysical repeated на9cells; nofactor9 передcentralcompression. Recomputed downcastrow5×45 наworstk7 даётphysicalmatrixproduct match≈6.65e−39absolute. Это structural audit маленьких продуктов, не повторныйBigLU.

Четыре counterfactual:raw; HermΛ; Hermcentralh; обеprojection. Λ реконструируется **вworkingprecision** какsavedlocalR+iΓ/2, затемHermΛ−iΓ/2. Snapshotне содержит original pre-subtraction FFTΛ, поэтому эта операция не является прямым all-node FFTerror certificate. Разныеvariants — разныеоператоры, не новое увеличениеprecision одногооператора.

## Precision и агрегирование

Primary192-bit traces/LUresiduals сохранены decimalstrings; HDFmatrices — Float64downcasts. q192 считается изstrings192-arithmetic, **не** изdowncastmatrixdiagonal. CSVstrings совпали summaryexactly; independentDecimal90 aggregation вgs*wE*wk/(2π) совпала до≈1e−58relative. Вкодеprefactor `gs/(2pi)` сначалаFloat64, затем promoted; это намеренно тотже discreteprefactor sourceq, не192-bitπ. Остальные quadratureweights также exactsavedDoubles.

Selectedscope:24(E,k)nodes×4variants. Energies104(firstnegative) и4145(worstnegative) имеют **все11k**; e7894/k1 иe30704/k11 — толькоpartialweightedmass. Именаnative_whole_q вpartialrecords — reference, не результатsinglek aggregate. Ни1429negativeq, ниall30704nodes highprecision этимoracle не перепроверены.

192-bit maxLUresidual≈2.24e−57, resolventA−effectiveidentity residualstrings малы. Precision384 повторена только наe4145/k7 — strongestnegative weightedk исходногоsnapshot — для4variants: traceabsolutechanges≈.8e−74–2.05e−74,relative≈2.8e−57–2.26e−52. Это убедительная precision устойчивость этой точки, **не per-node/global rigorous interval bound** и не384-bit wholeq sum. Остальныеnodes имеют192LU/identityevidence, но не384repeat.

## Wholeq на двух энергиях

| variant | q192(e104), all11k | q192(e4145), all11k |
|---|---:|---:|
|raw|−1.69173955043e−23|−6.69198874492e−22|
|Hermh|−1.69173955043e−23|−6.69198874492e−22|
|HermΛ|+8.12137069418e−30|+3.17661647417e−26|
|Hermboth|+8.12137042347e−30|+3.17661647471e−26|

Raw negative q совпадает с native Float trace до floating differences. Он не исчезает при Big192 и в проверенной 384-bit точке, поэтому **только transient Float64 subtraction/inverse cancellation** не объясняет знак данного summed trace. Вместе с Root full-node controlled HermΛ probe и independent same-G source decomposition это поддерживает numerical anti-Hermitian Λ contamination как главный factor initial q obstruction. Это не доказательство отрицательной физической плотности или отсутствия решения; source Γ в tail крайне мала.

## qtrace не является matrixPSD

Поsaved192-computedA, downcastтолько дляeigensolve, raw/HermΛ/Hermh имеютnegativeHermitianeigenvalues на23/24nodes. Минимумы соответственно≈−1.405e−16,−1.397e−16,−1.808e−18; Hermbothpositive24/24,minimum≈2.91e−29. Этиcomputedsigns/norm-relativebackwardscales — numericalwitnesses, неintervalcertificate.

СледовательноHermΛ-only решает **этотqtrace sign**, ноanti-Hermitianrawh продолжает формироватьindefiniteoffdiagonalA вexact-doublecounterfactual. SelectednormhdefectбольшеΛdefect, хотяегоtrace contribution почтисокращается; такоеcancellationнеустраняетmatrixproblem. ДляPSDidentityнадо отдельно соблюдатьHermitianH иconsistentΓphysical, либо обосновывать budgets численно.

Вторая numericalпроблема: Float64matrixA vsBigA может значительно различаться даже послеHermboth — maxrelative≈7.19e10, потомучтоphysicalAоченьмала. BigGR vsnativeGRприэтомrelative≈4.3e−16. МалыйGRrelativeerror/LUresidualне гарантируетдостоверноеtinyi(G−G†). Наe4145/k7 subtraction **downcastBigGR** отличается отsavedBigA на≈1.43e−17(raw/HermΛ),≈1.19e−19(Hermh),≈4.97e−33(Hermboth). Downcastматрицы нельзя использовать дляпереоценки192-bit signs. ВFloatHermboth oracle поsummaryестьmatrixerror, отсутствующийвBigphysicalidentity: высокаяprecision/stablephysicalrepresentationостаютсяотдельнойзадачей, анеследствиемоднойΛsymmetry.

## Decision

Принять narrow selected localization и precision witness. Не принимать Core fix/general causality из этого одного набора. Exploratory HermΛ warm16 остаётся ограниченным опытом с same guards/native audit; его nonnegative q само не принимает matrix PSD/native SCBA. Отдельный future implementation должен защитить Hermitian operator representation и stable spectral identity с явным error budget; clipping/floor/tolerance weakening не следуют из oracle. Φ/native/Poisson/discretization/experiment outcomes этот selected diagnostic не измеряет.
