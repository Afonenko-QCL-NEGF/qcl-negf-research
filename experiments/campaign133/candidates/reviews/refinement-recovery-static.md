# Static review: completed two-map streaming recovery

**Статических блокеров для root tiny→bounded recovery нет.** Область —
postprocessing уже завершённых двух Φ; новая Φ/retarded/context не
вызываются. Reviewer не запускал Julia/Pkg и не читал HDF5/source arrays;
только code/малые metadata, stdlib128MiB cap. Original attempt неизменён.

Frozen main a1b2d1c4312cda8b4559a5b55f16314ada019252826d3eea9fada85d2885f9e7,
input510cea294689feb34706d6b6512e6bfba28c334135cfe880990aa6c86d9123f0,
protocol89e91b754889d74e83ad17ba36245b189cea1909d9b9373d880448b980feef18
independently совпали с RECOVERY_SHA256SUMS. Source summary3135bba16e5867d2d89cb3f8112d7c2191fdc767f5665d1aa64377e3949c1a4f
проверен независимо:2 started/2 completed, stopped/Pi_complete,
exact InexactError Int64(1.2), coarse/fine map records существуют.
Это известный postprocessing failure после сохранения fine Φ, а не
failed scientific map. История отказа не переименовывается в успех.

## Индексы, меры, знаменатели

Общий fine_index(e)=2*e−1 используется actual и tiny, scientific-literal
ловушка устранена. Selected ladder shift bounds покрывают всё окно;
integer locator fields исключают прежнюю implicit Float64→Int conversion.
source_rows корректно включает right neighbour even последней строки
chunk. Actual-coordinate prolongation копирует все odd source values,
even θ∈[0,1] одинаков для всех complex entries/k.

Streaming covers все61407 fine rows один раз, включая30703 even и30704
odd. Coarse reads могут пересекаться ради neighbour interpolation, но
coarse accumulator учитывает только odd fine rows, поэтому двойного
учёта нет. Signature accumulate(ref,value) соблюдена:

- all/even D=Φf−PΦc, denominator PΦc, fine wE·wk;
- old-node D, denominatorΦc, coarse wE·wk;
- raw coarse/fine numerator ‖Γ−Φ‖, denominatorΦ своей сетки;
- Γ prolongation numerator Γf−PΓc, denominator PΓc;
- identity reference r_f=Φf−Γf, value Pr_c+D_f=Φf−PΓc.

Последний identity difference есть Γf−PΓc, поэтому отдельно измеренный
prolongation replay нужен и корректен. Это не независимый oracle Φ.
Нормы используют complex abs² всех entries и положительные source
weights, без spin/2π, clipping/Herm projection или нового μ fit. Zero
ratios следуют source convention; replay намеренно finite-only для
данного exact summary с finite positive scales. Метрики/dict locator
indices типизированы Float64/Int; нет прежнего неявного приведения.

## Input integrity/ресурсы

Exact source failure guard и summary SHA, copied attempt code/input SHA,
Core metadata проверяются до numerical reads. Все4 whole-file SHA
вычисляются до datasets. Model fixed expected47682acc…; три output
sentinels COMPUTE_AT_RUNTIME явно отличаются от historical expected SHA.
Они дают новое content binding, не доказательство отсутствия прошлой
мутации. Atomic writer history, full dataset shapes/sizes, все finite
streamed Γ/Φ chunks и replay6 raw norm fields по обеим сеткам дают
дополнительную фактическую проверку. Lambda shape проверена, но Lambda
не используется и не объявляется повторно validated array.

Nested energy/endpoints, exact source trapezoid wE, fine-input/output
axes и model wk byte hash связаны с immutable summary. Sourceμ/N
diagnostics только копируются с явным not_recomputed label. Full D norms,
subset/coarse norms и selected CSV сохраняются отдельно от original
summary. Source files открываются r; новый output обязателен.

Chunk≤4096 (actual1024) даёт O(chunk) поля: CΓ/CΦ, PΓ/PΦ, Γf/Φf и
небольшие arithmetic temporaries; полные arrays не загружаются. Source
hash/stream work входит во внутренние170s, external180s покрывает также
imports/I/O. Actual cgroup3GiB/Swap0/CPU≤100%,Julia1/BLAS1 проверяются,
output64MiB контролируется checkpoint. Hard timeout/cap должен задать
root; script сам не прерывает зависший primitive. Existing output
отвергается, checkpoint/failure сохраняются, original files не меняются.

Tiny защищает общий индекс/prolongation/accumulator, dense complex norm
oracle и known-failure admission. Это ещё не actual recovery evidence.
После completion new D metrics требуют отдельного result audit; новый
SCBA/Π run не нужен и не разрешается. Main cont64 должен закончиться
до root runtime по текущему общему ресурсному окну.
