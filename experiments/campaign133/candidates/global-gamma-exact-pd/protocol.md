# Exact PD Γ-last64: ограниченный протокол

Владелец этого диагностического кода — b0_review; выполнение только Root.
Никаких Julia/Pkg/HDF запусков автор не выполнял. Код не меняет Core/Runner,
Γ, h, Λ, μ, окно, допуски или baseline; HDF5 открывается только для чтения.
Frozen input имеет runtime_authorized=false и незаполненные actual SHA.

## Вопрос и доказательство

Проверяются **точные бинарные числа Float64** каждого из 337744 блоков
5×5 `Gamma-last.h5` stableA-cont64, т.е. Γ-вход принятой карты64.
Это НЕ raw Φ64, не Γ-best, не прежний Γ128, не исходная математическая
SCBA Γ с оценённой ошибкой. Номер64 и identity проверяются в самом HDF.

1. Все компоненты finite; Hermiticity проверяется точно, без усреднения.
   Единственное числовое отождествление IEEE — +0 и −0: они оба точный ноль.
2. IEEE bits дают x=m·2^e (normal: m=±(2^52+fraction), e=expbits−1075;
   subnormal: m=±fraction, e=−1074). Выбирается минимум e среди ненулевых
   компонент блока. M=2^(−e_min)Γ — матрица гауссовых целых BigInt.
   Все операции после декодирования целочисленные; Γ не ремонтируется.
3. Одношаговая fraction-free схема Bareiss без перестановок:
   B_next[i,j]=(p·B[i,j]−B[i,k]·B[k,j])/p_previous, p_previous(initial)=1.
   На шаге k, p=B[k,k]=det(M[1:k,1:k]). Деление обеих целых компонент
   проверяется `divrem`, остаток обязан быть точно0. Нереальный p или
   ненулевой остаток — INTERNAL FAILURE; научный pass невозможен.
4. Каждый p должен быть вещественным и строго положительным. По критерию
   Сильвестра это эквивалентно PD для exact Hermitian M. Положительный
   scale сохраняет PD; Δ_k(Γ)=p_k·2^(k·e_min). Pivoting не требуется:
   при первом p≤0 строгая PD уже исключена, продолжать деление не нужно.
   p=0 не доказывает отрицательную собственную величину; это отказ strictPD.

Первичный источник fraction-free тождества:
E. H. Bareiss, *Sylvester's Identity and Multistep Integer-Preserving Gaussian
Elimination*, Mathematics of Computation22(103),565–578(1968),
[DOI10.1090/S0025-5718-1968-0226829-0](https://doi.org/10.1090/S0025-5718-1968-0226829-0),
[publisher PDF](https://www.ams.org/mcom/1968-22-103/S0025-5718-1968-0226829-0/S0025-5718-1968-0226829-0.pdf).
Web-поиск подтвердил первичную запись/аннотацию AMS и индексированный текст
самой статьи (разделыIV–VI: integer-preserving operations, exact binary
input, commutative elements, условия pivoting). Непосредственная загрузка
PDF через web и urllib получила403, поэтому не заявляется чтение всего PDF.
Ранняя авторская версия: [ANL-7213, Bareiss1966, UNT/OSTI](https://digital.library.unt.edu/ark:/67531/metadc1035277/).
Доказательство применимости здесь не зависит от недоступного текста: формула
выше — тождество окаймлённых миноров, проверяемая точная делимость и
пять знаков явно заданы. Tiny oracle — другой алгоритм, rational LDL†
с точным восстановлением матрицы и сравнением ВСЕХ5 principal determinants.

## Источники и исполнение

Root заполняет все9 whole-file SHA в отдельном input-run. Все хеши проверяются
потоково **до первого открытия HDF**: Γ/model/input/identity/map64/operator/
spectral/ledger/summary. Дополнительно source identity связывает input,
оператор и stable_spectral; map64 и standalone ledger должны совпадать.
Input обязан иметь project_AA_trial=false и intended all9 hH policy.
Γ остаётся unprojected; Hermitian hH и Λ policy относятся к источнику.

Рекомендуемый порядок (команды запускает Root в существующей среде):

1. `julia --project=... -t1 run_certificate.jl --tiny`.
2. `julia --project=... -t1 run_certificate.jl --fixtures input-run.toml`:
   те же manufactured/negative/zero/nonHerm/nonfinite/subnormal/scaling тесты
   плюс44 frozen5×5 Γ fixtures из прежнего selected qualifier. Все fixture
   determinants сравниваются с exact RationalLDL; никакого полного source pass.
3. `mode="benchmark5000"`: только первые5000 блоков (E-major, k-minor).
   Whole source SHA/HDF checks всё равно обязательны. Измерить wall/RSS/GC.
   Этот результат может принять только данный subset. Не экстраполировать
   уверенность PD или формальную стоимость на оставшиеся332744 блоков.
4. Если фактическая стоимость допускает: отдельный input/hash `mode="all"`,
   новый пустой output. ALL проходит все30704×11=337744 блоков.

Внешний предел:300s/4GiB/Julia4/CPU400%/Swap0/file32MiB, без другого главного
расчёта. Внутренний wall≤280s проверяется между chunks≤1024E, по умолчанию270s.
Он не ограничивает время внутри одного chunk/whole SHA: внешний deadline
обязателен. Чтение HDF только main thread, chunks512E (~2×1.1MiB);
workers независимы и владеют5×5 Float/BigInt scratch. Нет per-block
Rational/Dict/strings; Dict/strings только first/failure≤8perworker и
chunk summary. Обновление минимального exact determinant witness копирует
только fixed-size5×5 Float и5BigInt pivots в owning worker-local scratch;
строки/Dict при обновлении worst не создаются. Tiny cache-control проверяет
сохранность witness после перезаписи working scratch и TOML roundtrip.
Физический minΔ5
сравнивается целочисленно с учётом scale, не по Float-приближению и не
выдаётся за min eigenvalue. Mutable BigInt unsafe GMP API не используется.

На chunk boundary сохраняется summary с execution=running. Только complete
ALL+337744 PD+без исключений даёт all_blocks_strict_PD_certified=true.
Wall/interrupt/partial/не-Herm/nonfinite/p≤0 не принимаются. Вывод включает
counts, divrem-count, integer bit bound, первый/минимальныйdet/первые8failure
witnesses с exact integer principal minors и exponent; максимум32MiB.
ALL положительный ожидает ровно337744×60=20264640 component divrem checks.
V1 сохранён полностью в unexecuted-v1; V2 меняет только worker worst cache
и момент его сериализации, exact_pd.jl и source/input gates неизменны.

## Научная граница

Успех устанавливает strict PD всех сохранённых бинарных Γ-узлов именно
этого Γ-last64. Положительная scalar scale единиц не меняет этот вывод.
Convex hat interpolation переносит PSD между узлами при соответствующем
интерполянте; support/ghosts/outside-window atoms и intended Hermitian
all9 hH должны отдельно соответствовать будущему Mats-source варианту.
PD — необходимая конкретная опора этого будущего анализа, а не полный
сертификат причинности, отсутствия атомов или плотности.
Нет вывода о точности fullΠ, глобальной forward error, SCBA convergence,
native acceptance, достаточности energy grid, эксперименте или Γ128.
