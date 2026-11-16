# Static review: fine μ / thermal-LO attribution

**Математический diagnostic и source-free tiny допустимы. Unmodified
actual пока HOLD по F1 performance; source binding также ещё UNSET.**
Рекомендован конечный accumulator-only fix, затем новая frozen delta
и root tiny, до отдельного actual окна после active main cont64.
Reviewer Julia/Pkg/HDF5 не запускал; читал code/малые metadata со
stdlib memory cap128MiB, source не редактировал.

Main f3bfd44533439f141e563de2117692ea1265670787d0fdb5bd886de692788a2d,
tiny aa843f234a9b30e2d739ae8b3c0490d377f2467294a61720c0ec966e5e47843e,
input3ecab429ce80c2a3c4a878368251b179ea85a389f9bcdebd6ee542507f8403be
и protocola8a3d7a9eac37eb1ed10a5e5184b733d65da3aa85530b762a896e93388dba266
independently совпали с SHA256SUMS. Imported refinement5084/operators7d71/
tiny565a и Gamma/warm/restart literal hashes имеют корректные bindings
к уже reviewed source. Include main guards не запускают original maps.

## Алгебра и observables

Для одного и того же frozen fine A, thermal Bose population N и
kernel K native LO даёт

    ΓLO(f)=K[(N+1)S−A + N S+A + S+(fA) − S−(fA)].

Следовательно M=ΓLO(μf)−ΓLO(μc)=K[S+(Δf A)−S−(Δf A)], без Bose
factor в delta. Это согласуется с independent tiny oracle, включая
sign greater−lesser и shift orientation. Four static Γ contractions
зависят от A, поэтому не меняются при μ на том же A. Тогда

    D=Φf(μf)−PΦc(μc), F=D−M=Φf(μc)−PΦc(μc).

Знак subtraction правильный. R=Φf−Γf и R−M тоже относятся к одному
fine input. Это algebraic μ attribution **данного конечного operator**,
не новая Φ, nonlinear convergence или принятие fixed-density solution.

Saved metadata: μf=.08459119027608916, μc=.08459613668783181;
fine fixed-μ witness использовал именно μc. Δμ=−4.946411742645651e−6
dimensionless. Fresh q/N replay на обеих stored μ до LO и LO(μf) norm
replay используют original1e−10 arithmetic guard. μ fit отсутствует.
q включает g_s·wE·wk/(2π) ровно один раз; LO содержит native qK/wk.
LO norm replay — полезная sanity check, не full field equivalence proof.

Accumulator использует complex abs² и real(conj(X)Y), поэтому
‖D‖²=‖M‖²+‖F‖²+2Re〈M,F〉 и weighted analogue корректны без
предположений об Hermiticity. Denominator всех reported relative norms
есть PΦc; это явно отличается от raw residual denominatorΦf. Norm
ratios/projection coefficients нельзя выдавать за неотрицательную
«долю причины»: M/F могут усиливать или компенсировать друг друга.
Old-node/midpoint partitions используют fine weights, не coarse norm.
Все61407 rows/all complex entries покрыты; selected CSV только supplement.
Explicit2*old−1 индекс не содержит прежней exponent-literal ловушки.

## F1: inner-loop allocation cost, actual HOLD

`add_mu_moments!` создаёт18 динамических строковых keys на каждый вызов
через name*"_squared"/etc и обновляет Dict entries. На fine grid
61407·11·25=16,886,925 scalar entries, helper вызывается дважды (full и
partition): около607,929,300 string-key updates/concatenations, плюс
Dict lookups. Это статический operation count, не runtime measurement.
Живой payload может оставаться малым, но cumulative pool allocation/GC
не защищён cap12.125GiB и может поглотить300s. Tiny этих размеров не
проверяет. Это повторяет класс риска прежнего boxed/GC-heavy audit.

Minimal closure: typed numeric accumulator/fields или fixed-index
Float64 vector без string concatenation/Dict updates в orbital inner
loop; conversion к named summary один раз в finish. Hoist full/partition
references из scalar loop. Preserve formulas/weights/identity/tiny
oracle; новый backend/physics/process framework не нужен. Root/автор
уведомлены; reviewer source не правил и не запускает benchmark.

## F2: finite stream guards

В текущем stream отсутствует finite check PΦc/fineΦ/Γ chunks. NaN
только в PΦc даёт NaN D/reference, но raw fine residual replay от PΦc
не зависит; theoretical corrupted input мог бы завершиться с NaN metrics.
Для actual immutable source прежний GP.cone и completed recovery с
exact matching whole-file SHAs могут отдельно установить finite fields;
это делает F2 nonblocking при таком input binding. Cheap explicit
all(isfinite) chunks/delta и finite accumulator outputs лучше замыкают
failure behavior непосредственно в этом diagnostic и рекомендованы
вместе с F1. Нельзя принять NaN как измеренный ноль.

## Source/resource lifecycle

Перед первым source HDF5 обязательно7 whole-file hashes; input с UNSET
не готов к actual. Model/h/basis/K/frame сверяются exact, fine Γ digest
с saved summary, μ/N и LO norm с сохранёнными scalars; Manifest совпадает.
Whole-file hashes повторяются после diagnostic. Original sources только r.
Core e171/clean и native guard сохранены. Imported parent workers4 при
Julia1 разрешается Core min(requested,nthreads) в одного actual worker;
явная нормализация parent workers кcfg1 улучшила бы provenance, но это
не API failure. Actual numerical и OS budgets оба12.125GiB,High11.5GiB,
Swap0; скрытого8GiB scope/12GiB options mismatch нет.

Один retarded call site, два LO calls с counters до работы, zero full
Φ/fit/native solve/retries. First LO returned field удерживается как
delta, second LO вычитается in place; Green/LO intermediates не выходят
из barrier. После второго LO full GR/A/etc освобождаются; остаются один
delta field и O(chunk) saved Φ/Γ chunks. Это ниже reviewed24F+3GiB
envelope по payload, но не измеренный RSS. Main issue F1 — allocation
cost, а не необходимость большего live-memory budget.

Внешние J1/BLAS1/CPU≤100%/memory/high/swap/file limits читаются gate;
soft290s/external300s,data64MiB/file32MiB. Root задаёт actual scope и
нет-overlap window после cont64. Без binding closure и accumulator
delta этот review не принимает actual execution; source-free tiny
может проверять native thermal-LO independent oracle и norms отдельно.

## V2 delta closure:68ccd freeze

**F1/F2 закрыты. Статических блокеров для bounded1R+2LO actual после
заполнения7 whole SHA/root resource admission нет.** Старое actual HOLD
по accumulator относится к V1, не к этой исправленной frozen delta.
Новая физическая Φ или μ fit по-прежнему не разрешены.

Проверены все5 SHA256SUMS entries: main68ccd9d03dc6e4861c31798bad714ca6ed0b74a129882d44b2fa5c4d938cb8f6,
tiny a9ac669891f72f613b581df80c0738a6ce6a76629c578dc1f50fcadf441e25db,
protocol1c6bf7beb80db7b1fb23846bd8c7a40f629939b640f053f2438f1a5b26310053,
unchanged input3ecab429…, patch8d4e2327895a1b87c20b009b28e606d86fb1b5d6a52feee90b44948d20f69d61.
Archived V1 четыре hashes совпали с предыдущим freeze.

Полный patch прочитан:18 slots correctly map D/M/F/reference/R/R−M
squared+weighted(1:12), затем D·M/D·F/M·F inner+weighted(13:18).
Columns1/2/3 —full/old/mid. Каждое scalar expression и per-slot addition
order сохранены; finish_once mapping labels/denominators/identities
согласован. Hot loop не создаёт строки и не обновляет Dict. Это
устранение выявленного structural allocation risk, не заявленный benchmark.
Chunk fine/PΦc/Γ finite guard стоит до arithmetic accumulation.
New tiny включает bitwise сравнение всех18 original dictionary values
и отдельные independent norm/dot identities, NaN rejection для всех3
source chunks. LO/retarded/source/resource code не менялся.

Прочитан actual root `tiny-v2-correct-root.log`:3+13=16/16 PASS,
exit0,16.23s,peak815940KiB,Swap0. Это source-free fixture; source HDF5
и actual1R2LO ещё не приняты этим review. Seven UNSET bindings должны
быть заполнены отдельно и сохраняться в runtime input hash. После
recovery обнаружен global D maximum e39394 вне old9785 ladder; selected
CSV coverage ≈5.19% weighted D². Full moments всё равно покрывают все
nodes; добавление runtime anchor19697 может показать этот maximum в
supplement CSV без изменения operator/counts, но должно быть явно
записано в новой runtime identity. Нельзя выдавать sparse CSV за полную
локализацию даже при успешном μ attribution.
