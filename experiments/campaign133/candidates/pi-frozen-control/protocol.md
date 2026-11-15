# Frozen local Γ: uniform PI versus nodal FFT, campaign133

Статус: протокол/код для отдельного разрешаемого запуска; сейчас write-only,
без Julia/Pkg/SCBA. Core не меняется. Требуются review и resource-green ведущего.

## Вопрос и фиксированные входы

Изменит ли замена nodal omitted-cell Hilbert quadrature на существующий
piecewise-linear product integration существенно frozen retarded spectrum и
его конечные интегралы? Это sensitivity одного retarded operator, не новый
SCBA шаг, не доказательство baseline bug и не приёмка nonlinear решения.

Источник `.build/campaign133/raw/S01-warm-FD64/state.h5`, v3,
SHA256 `e6fffb997f36397f877aee1e2c98b1dbad7d97417d45458a042b5f074ce57a6d`.
Input/identity/scales и семь saved basis operators связаны с этим H5.
Core B v2 `2a500740a82dfd5b44f6be966ef63dd54010d097` (или явно записанный
после review descendant с той же scientific реализацией). Проверить actual
import path/HEAD/dirty status. Source input SHA из identity
`0fa0d0bc2fc27602fa99b8fe0797071eae99c66ab2661dfde05bd00ebe19b86e`.
Source unresolved FD64: SCBA не сошлась; Hartree=0, T=70K, E=[−.3,.9]eV,
NE30721, Nk11, Nb5; finite_chain4, source frame из v3. Никаких kernel builds.

Frozen local ΣR = source IFR+acoustic+impurity+alloy+LO **ровно один раз**.
Все embedding families исключены. Γlocal=i(ΣR−ΣR†), один и тот же array для
обоих reconstructors. No clipping, positivity repair или tail modification.
μ=0.008495393439818052eV из source FD64 marker: только fixed FD weight,
не новая fitted μ, не восстановление actual nonequilibrium Gn.

## Два свежих retarded arms и archived control

1. Native nodal `ProductionFFTHilbertPlan` / weighted omitted-cell convolution.
2. Native `_ProductionPIFFTHilbertPlan` / identical frozen Γ, full ghost hats,
   без дополнительного wE factor. Не вызывать dense H при NE30721.
3. Для каждого ΣR=Λ−iΓ/2 построить **свежий** native finite_chain4 embedding,
   тот же h(k), T+, basis/energy/k/weights, нулевой bias. Затем native Dyson
   GR=[EI−h−ΣRlocal−ΣRembedding]⁻¹ без η, A=i(GR−GR†).

Отдельного native retarded-only cavity entrypoint в frozen B нет. Harness
использует `_cavity_embedding_production` с общим нулевым lesser/greater
buffer и потребляет только retarded output; dummy Keldysh outputs сразу
освобождаются. Это техническое выполнение native retarded recursion, не
физическая Keldysh задача. Не запускать scattering/candidate/normalization/
SCBA/Poisson; не сохранять новые full states.

Три labels: `archived_source_A`, `fresh_nodal_frozen_local_gamma`,
`fresh_PI_frozen_local_gamma`. Source embedding может быть смешанным, не
fresh cavity текущего local Σ. Поэтому PI−fresh nodal различает выбранный
operator; fresh nodal−archive показывает alignment/history differences и
не приписывается одной причине.

## Наблюдаемые и меры

До embedding: nodal rebuild−archived local ΣR (history/linearity sanity),
PI−nodal Λ absolute и relative Frobenius, weighted global/occupied norms,
Hermiticity Λ, reconstructed Γ differences, min Γ eigenvalue/scale и witness.
Сохранить primary numerators/denominators, Γ hash, plan buffers, PI geometry,
arithmetic witnesses L/2L; witness не bound Hilbert/window error.

После embedding: M_k=∫A(E,k)dE/(2π), F_k=∫f(E;source μ,70K)A(E,k)dE/(2π)
по **исходной positive trapezoid measure**. N[fA]=gs Σwk TrF, без нормировки.
r_sum=maxk||M_k−I5||F/√5. Сохранять full small matrices, per-k и global
integral differences, N/target и δN/target. Nodewise spectral differences
измерять в weighted L2 для A и fA; occupied-relative denominator norm(fA)
и absolute numerator сохраняются отдельно. Bounds на integral differences:
||δM_k||F≤∫||δA||F dE/(2π), |δN|/target≤gsΣwk√5∫f||δA||F dE/(2π)/target.
Это discrete triangle bounds выбранной меры, не continuum error bounds.
PSD A/Γ и conditioning — floating witnesses; no sign-only acceptance.

## Дополнительные matrix edge witnesses

Если бюджет позволяет, ideal PI ghost edges ℓ=a−h и r=b+h. Вычислить прямо
Λ(ℓ)=Σj k(−j)Γj, Λ(r)=Σj k(N+1−j)Γj, stable native coefficients up to lagN,
без quadrature/extrapolation/wE. Для каждого11 k сохранить D45=I9⊗[EI5−h−Λ]
с offdiagonal −T+/−T+† и его eigenvalues/eigenvectors residual/orthogonality,
anti-Hermitian defect, λmax(Dleft), λmin(Dright). Source T− согласовать с T+†.
Не заменять full45 sign знаком d5. Numerical eigensigns/residuals не interval
certificate без input/arithmetic error bounds. Ideal uniform geometry имеет
Float64 representation caveat. Certificate expert рассматривает monotonic
D′≥I outside ghost support; original [a,b] всё ещё отсекает ghost strips;
rank-deficient Γ допускает interior dark atoms. Failed witness не означает,
что pole виден в central5 block; positive witness не принимает SCBA/window.

## Ресурс и выходы

Проспективно: один процесс Julia4, BLAS1, OS cap6GiB, timeout300s с TERM/
KILL grace10s, outputs≤32MiB. Никакого install/Pkg. До запуска убедиться в
отсутствии другого тяжёлого опыта. Сохранить `/usr/bin/time -v` и фактический
cap из запуска ведущего. Script проверяет elapsed/output на фазовых границах;
OS cap/timeout задаются снаружи штатными средствами.

Outputs: summary.toml, small integrals/witness matrices H5, selected block
witnesses и source/input/code/project/Manifest hashes. Full Γ/Σ/GR/A не
сохраняются. Исходный H5 не изменяется. При failure сохранять stage/status,
никакой нулевой replacement недоступных результатов. Без root launch и
независимого review runtime остаётся not_measured.

Проспективная команда (только после root resource-green и static review;
OS cap6GiB/CPU400% задаёт ведущий снаружи, не эта команда):

```sh
OPENBLAS_NUM_THREADS=1 /usr/bin/time -v timeout --signal=TERM --kill-after=10s 300s julia --startup-file=no --check-bounds=yes --threads=4 --project=.build/campaign133/environment/pi-B .build/campaign133/pi-frozen-control/frozen_control.jl .build/campaign133/pi-frozen-control/run-v1
```
