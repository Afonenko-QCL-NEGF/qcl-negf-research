# Independent static review: exact PD of unprojected ΓLast64

**Блокеров для root tiny→44 fixture controls→bounded benchmark5000 не
найдено.** ALL337744 требует отдельного root cost/resource decision
после controls/benchmark, нового input mode/hash и полного coverage.
Reviewer Julia/Pkg/HDF5/scientific compute не запускал; только code,
малые hashes/TOML, stdlib≤128MiB. Γ/Core/Runner не редактировались.

## Frozen source

Все5 SHA256SUMS entries independently verified:

| Artifact | SHA256 |
|---|---|
| exact_pd.jl |c6b14339d05720f4d55152b35edd9d43dca910507bfd16c8e750a41cdf02b6a9|
| run_certificate.jl V2 |d2395572bd724f9493be4d5c5cfb92fe9859e609719cafa0bced68f7de02934f|
| pending input |78305e280bbf590bd7ef233248a527cddf5e727c8e168b70b563a9fddc67c8d8|
| protocol |a56f0a2cfb167f8a9ef296b8983528e57da22defca1067f6715755e31a93a30d|
| V1→V2 patch |e30a48452c8f8bd31b0b057d15bede604f69952fbc47c2b68b3ac9bdec28fc17|

V1 exact_pd/input byte-identical; V2 changes only owning worst-witness
cache/deferred serialization and adds cache/TOML tests. All9 pending
binding paths resolve to the intended stableA-cont64 attempt. Их hashes
ещё ROOT_FILL: actual не готов до root runtime input freeze. Root-owned
stableA-cont64-source-digests.json даёт exact Γ-last cab81029…, model
700ae176… и7 small bindings. Малые7 hashes independently confirmed;
большие HDF5 reviewer не открывал/не хешировал.

source_gate проверяет все9 whole SHA до первого HDF open, exact source
identity/call64/dimensions, project_AA_trial=false и intended spectral
policy. Identity links operator/input/stable_spectral; standalone ledger
равен map64 ledger, complete summary last_map равен map64. Archived
writer gamma_checkpoint:409 сохраняет именно input Gamma, schema
checkpoint-v1/identity/pi_calls_total; finally:762 сохраняет последний
accepted gamma. Last accepted map64 guard исключает failed-final-trial
подмену. Не читаются raw_Phi, Γ-best или projected64 files. Actual source
code не делает PSD clipping: Γ=accepted AA trial; Hermitian working-Φ
и all9hH policy не переименовываются в PSD projection Γ.

## Exact arithmetic proof

dyadic декодирует IEEE64 bits без Float arithmetic: normal signed
mantissa53bits и exponent ef−1075; subnormal fraction/exponent−1074.
Int64 mantissa помещает ±(2^53−1); Inf/NaN reject. Signed zeros обоих
знаков математически0. emin берётся только по nonzero components;
BigInt shifts неотрицательны, B=2^(−emin)Γ — Gaussian-integer matrix,
положительный real scale сохраняет PD.

Hermiticity проверяется exact binary equality real transpose/imag
negative transpose, с единственным отождествлением signed zeros. На
diagonal это требует imag=0. Нет усреднения/rounded eigensolver tolerance
или repair. Matrix size5×5/eltypeFloat64 фиксируются actual HDF guards.

Bareiss trailing update

    B_next[i,j]=(p·B[i,j]−B[i,k]·B[k,j])/previous

использует неизменённые row/column k: update indices оба>k, поэтому
in-place порядок не портит operands. По тождеству окаймлённых миноров
p на stage k есть det original leading k×k; previous=leading(k−1),
initial1. На каждом positive stage предыдущий divisor ненулевой real.
Exact complex numerator делится на real BigInt по двум divrem; оба
остатка проверяются, нарушение — INTERNAL error, не научный pass.

Для exact Hermitian матрицы strict positivity всех5 leading minors
эквивалентна PD (Sylvester criterion). Первый p≤0 достаточно отвергает
strictPD; p=0 не доказывает отрицательную eigenvalue. Pivoting не нужен
для positive route, а после отказа certificate уже невозможен. Physical
minor rule Δk(Γ)=p_k·2^(k·emin) верен **в исходных stored native units**;
это не eV^k без дополнительного E0^k. Positive unit conversion не меняет
sign; minimum Δ5 не является minimum eigenvalue/condition estimate.

На PD block выполняются2(4²+3²+2²+1²)=60 component exact divisions;
ALL337744 требует20,264,640 таких checks. Tiny сравнивает все5 minors
с independently reconstructed exact rational LDL†, включая complex Gram,
negative/zero/nonHerm/nonfinite/subnormal/scaling/signed-zero controls.
Rational Float conversion ранее independently checked exact (tol0).
Это другой алгоритм проверки, не Float eig sign.44 старых selected
fixtures полезны как realistic exact-input algorithm controls; они не
являются44 проверенными блоками **нового Γ64** и не входят в ALL count.

Primary algorithm provenance: Bareiss1966/ANL-7213
[UNT primary-report record](https://digital.library.unt.edu/ark:/67531/metadc1035277/).
Publisher1968 AMS page снова403, OSTI full-report endpoint502; reviewer
не заявляет чтение полного PDF. Применимость формулы выше проверена
через exact minor identity и код; библиография не заменяет это доказательство.

## Ownership, counters и fail-safe outcome

Каждый Worker владеет Scratch; @threads static назначает одному worker
отдельный strided subset indices. Для каждого chunk indices1:count
partition без пропусков/повторов, E-major/k-minor mapping корректен.
HDF read main thread до parallel work; aggregate после barrier.
Only ordinary allocating BigInt API: aliases pivots/previous/worst_integer
не мутируются unsafe GMP, а vector cells заменяются новыми values.
V2 worst cache копирует X/pivots/exponent в owning Scratch; перезапись
working scratch не меняет cached witness. Строки/Dict создаются только
first/failure bounded records и chunk summary; per-block rational/Dict/
string construction отсутствует, кроме неизбежных BigInt арифметических
allocations. Cache tiny проверяет worst после overwrite и TOML roundtrip.

benchmark5000 читает455E rows, но учитывает ровно первые5000 blocks;
last5 прочитанных k не становятся checked. ALL покрывает30704·11=
337744. Statuses nonfinite/nonHerm/nonpositive/pd входят в totalcounts.
ALL flag требует complete +modeall +executioncomplete +pdcount337744;
subset flag не подменяет ALL. Ошибка/timeout/incomplete не принимаются.
Early summaries сохраняют all flag=false; fail/finally не превращает
partial в pass. Наличие failed block может coexist с complete scan,
но PD flag тогда false.

First8 failures и first/worst positive records bounded; witness row-major
serialization согласована, BigInt minors становятся strings только
на output boundary. Worst comparison dyadic_less учитывает exponent
5emin exact, не сравнивает несоизмеримые integerized determinants.

## Resource / scientific boundary

Root должен enforce4GiB/J≤4/CPU≤400%/Swap0/300s/file32MiB и no-overlap.
Код ограничивает Threads≤4, chunks≤1024E, internal≤280s и output summary
≤32MiB, но сам не задаёт cgroup/timeout. Hashes входят во время;
проверки wall между chunks не прерывают отдельный hash/chunk. Chunk512E
две Float slabs суммарно≈2.15MiB, fixed5×5 worker scratch; BigInt temporary
cost нельзя оценить только payload. First5000 не являются формальным
worst-case time/memory bound для остальных332744; их exponent/condition
distribution может отличаться. Сохранять partial evidence, не повышать
caps автоматически. Actual runtime/RSS/ALL пока not_measured.

Успех ALL установит exact strictPD **сохранённых binary Γ64 узлов**.
Positive hats/convex interpolation могут дать PSD continuum конкретного
interpolant, но future intended all9 hH/source/ghost definitions требуют
отдельной matched identity. Точная Hermitian hH и causal Cauchy model
не подменяются native rounded Lambda/source ledger одной PD меткой.
Это prerequisite для будущей causal-Mats постановки, не её density/tail/
atom/window/error certificate, не forward accuracy полногоΠ и не SCBA
convergence. Старый Γ128 Mats HOLD и qualification прошлых raw-h
Mats functionals не снимаются этим новым Γ64 certificate.
