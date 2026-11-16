# Intended split-model selected Big qualifier

Статус: write-only freeze V1. Автор не запускал Julia/Pkg/HDF5/source arrays.
Selected SHA a74dd7... предоставлен Root/независимым reviewer; малые input,
summary и preflight source SHA пересчитаны автором. Никаких Core/Runner edits,
solver maps, density fitting или tolerance changes. Root запускает только
после independent static review и tiny pass, под2GiB/120s.

Цель — точность Double covariance A относительно **intended split primitives**,
а не к ΣR, в которой tiny offdiagonal Γ могла исчезнуть при Float assembly.
Source preflight a34... завершён с3started/completed operations; input de3...
и summary06d6... pinned. Новый qualifier читает только selected1.3MB HDF.
Все4 energies[7893,8038,30696,30704]×11k =44blocks, ни одного полного real-grid
retarded pass. Source grid30704/11k/5basis и Julia field(E_selected,k,a,b)
проверены. По запросу Root это фиксированный selected set, не auto-current.

## Уравнение intended model

В BigFloat192 **до операций** переносятся saved Float64 h_H, original
Lambda_H_preassembly, Gamma, Tplus и actual E. Собирается45×45 open chain
с9 одинаковыми cells:

    Dcc = E I - h_H - Lambda_H_preassembly + i Gamma/2,
    D[c,c+1] = -Tplus,  D[c+1,c] = -Tplus†.

Никакого повторного Herm projection h/Lambda/Gamma нет. Все9cells используют
тот же h_H, включаяcentral; это отличается от прежнего raw-h obstruction
oracle с8sideHerm/centralraw и соответствует данному preflight variant.
No eta. Γ/G/R из старыхsolverstates не загружаются. Original localR/totalR
не используются для Big assembly и не подменяют preassembly Lambda buffer.

Adjoint solve D†Y=P с5RHS даёт centralROW=(Y†), не centralCOLUMN D\P.
Gc=row[:,central]. Main Big observable:

    Aprod = Σ_(9cells) row_c Gamma row_c†.
    Asub_big = i(Gc-Gc†).

Сохраняются Big product/subtraction consistency, LU backward residual/scales,
Double covariance payloadA/Asub/GR absolute+relative errors относительноBig,
Float64 SVD condition estimate **downcast intendedD**. Последнее не является
verified condition bound exact Big model. Double-vs-Big discrepancies включают
split-vs-assembledΣR, cavity, solve и A arithmetic; одной LU они не приписываются.
No arbitrary norm floor; 0/0=0, nonzero/0=Inf.

## Exact selected Γ PD

Base Rational{BigInt}(Float64) использует exact tol=0 conversion (Julia Base
rational.jl:175). Для каждого5×5Γ проверяется exact Hermiticity и no-pivot
LDL† с rational arithmetic. Пять положительных real pivots **и exact
reconstruction Q=L D L†** дают strictPD proof данной exact binary matrix.
Сохраняются exact numerators/denominators/signs всехpivots. Ни eigenvalues,
ни pivots не clip'ятся; nonpositive pivot означает no strictPD certificate,
не автоматически «matrix indefinite» при zero pivot. Non-Hermitian input
останавливает declared Hermitian model contract, а не чинится.

Вычисленные Float Γ eigensigns и64n eps*normΓ — только numerical comparison.
ExactPD proof относится к44selectedΓ, не всему grid и не uncertainty исходной
QCL model. При exactHerm h/Lambda иPDΓ intendedD invertible: ImD=I9⊗Γ/2 PD,
а exact Aprod PD как compression repeated-source congruence. Это algebraic
selected-model premise; approximate Big inverse всё ещё требует residual/
precision evidence, и full-grid causality/SCBA из него не следует.

## Precision, aggregation и output

Все44Big evaluations192bit. Primary traces/LUresiduals остаются decimal
strings192, а q для каждойselectedE агрегируется поall11k в192bit arithmetic
изэтихtraces. Prefactor gs/(2π) и sourcewE/wk — promoted exact Float64 values,
какsource numerical measure; это не пересчёт continuum π/weights. qselected
не является sheet densityintegral по всем E. Stored Doubleq_cov/q_sub имеют
явные labels и не выдаются за Big aggregates.

Один384-bit repeat выбирается по **максимальному relative DoubleAcov vs BigAprod
error среди44nodes**, правило фиксировано доrun. Сравниваются обеBig traces и
consistency; это не384proof всех44 или wholegrid. Никакого выбора «красивого»
результата послеизмерения. Critical HDF содержит downcasts192матриц дляworstCov
иworstSub nodes; original192/384 traces/residuals вTOML — primaryprecisiondata.
traces.csv сохраняет44per-k строки для independentaggregation.

## Tiny и ресурсы

Tiny controls: exact0.1 binaryrational; complex2×2ΓPD с известным вторымLDL
pivot19/8; negative/nonHermitian controls; independent complex/non-normal
three-cell Schur formula дляGR иsourceA; weakΓ1e−25case сBigproduct/subtraction;
nestedTOML proof roundtrip. Written, **not pass**. Tiny не читает actual HDF.

Root внешний scope: Julia1/BLAS1/heap1G, hard2GiB/Swap0/CPU100%, timeout120s,
file8MiB. Internal110s checks междуnodes не заменяютOSlimit; процесс/JIT
можетначатьсяраньшеinternaltimer. Expectedoutput<2MiB, onlysmallcriticalmatrices.
Environment existingaccepted используется дляHDF5/stdlib, безCore imports/Pkg.

    julia --startup-file=no --threads=1 --heap-size-hint=1G \
      --project=.build/campaign133/environment/accepted \
      .build/campaign133/positive-source-selected-big/oracle.jl --tiny

Actual: Root делает authorized inputcopy рядом сdefault, меняя толькоruntime
flag и сохраняя SHA, затемтотже entrypoint INPUT_RUN NEW_OUTPUT. ДоsourceHDF
read streaming проверяютсяselected/input/summary/preflightscriptSHAs.
Executioncomplete, selectedmodelPD иaccuratecovariance — разныеутверждения.
Scientificacceptancefalse всегда; никакого wholeΠ/fixedpoint/native/Poisson/
discretization/experiment pass этот qualifier не объявляет.
