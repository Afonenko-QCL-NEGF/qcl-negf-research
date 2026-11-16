# Independent delta review: M3 postprocess BLAS V2

**GREEN для одного нового bounded attempt02**, после Root getter preflight,
в прежних120s/2GiB/J1/file32MiB/output32MiB/Swap0. Это исправление выполнения,
не новый научный результат; reviewer Julia/HDF не запускал.

Frozen adapter7edebecad9df0327894dab5e374b78be79fc0d687c5f328092340ea2a651bdf8;
false inpute197b6afee61250ffc3750abb22d01f464d3feb8d08f4118f96dc7a83a1a0ef7.
Все7 SHA256SUMS и4 text guards проверены; parsed input против archived V1
меняет **только** `text_sha256.adapter`. Module9182ae80… и tiny4c76270d…
сохранились,9/9 прежних tiny проверок относятся к неизменённой математике.
Old adapter94c6f2… не менялся. Transient input f945… имел ошибочный old-adapter
SHA, reviewer отметил её; в final e197… она устранена до выполнения V2.

**Фактический исход V1.** External wrapper exit124/120,12s,peak755012KiB,
swap0,CPU187%; N128 сохранён, N256 прерван в KLO basis construction.
Whole-pair summary отсутствует: не completed и не scientific failure.
Raw/log/source/частичный N128 остаются исходными, archived V1 adapter3fe642…,
inputc9b451… и manifest сохранены отдельно. Новых Φ или μ-fit не было.
N128 нельзя выдавать за принятый complete128/256 postprocess.

В V1 фактическое число BLAS threads **не записывалось**. Stack
zgemm_thread_nn/GOMP_parallel и отдельный Root getter diagnostic совместимы
с oversubscription: OPENBLAS env1, фактический getter128; после explicit
set_num_threads(1) getter1. Диагностический getter не является прямым
измерением уже завершённого V1 процесса. Прежний static review не обнаружил
отсутствие явного setter в reused source-only path: old `run_actual`
делал set1, а old gate/source этого не делали. Одного env1 на этой сборке
недостаточно. Тайминг V1 не представляет требуемый BLAS1 режим.

**Проверенная дельта.** После импорта LinearAlgebra, до moment include и
любой source/model/matrix работы, `blas1` вызывает setter1 и требует getter1.
Затем повторяется после dynamic imports, до reused gate/source и после
source, перед каждым M3 calculation. При несовпадении getter выполняется
error; нет fallback, расширения времени или тихого изменения модели.
Startup log, report до вычислений и final report содержат фактические
BLAS/Julia getters, не только ожидаемые input options. J1 и прежние OS
пределы проверяются неизменённым gate; внешние caps остаются обязательными.

Diff не меняет moments/tails/K/source binding/linear system/нормы/μ/N128/256
или критерии. Точность будущего результата и floating differences при другом
порядке BLAS reduction требуют actual данных; побитовая эквивалентность
V1 partial здесь не заявляется. Новая пара должна быть сохранена в новом
output с runtime-authorized copy; scientific scope остаётся прежним:
one-shot и новая selfconsistent asymptotic M3 closure, не infinite-tail
certificate, infinite-N SCBA acceptance или proof ускорения алгоритма.
