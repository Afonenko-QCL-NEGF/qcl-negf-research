# Independent result review: ALL exact PD ΓLast64

**Узкий ALL certificate принят: все337744 сохранённых binary Float64
Γ-блока unprojected stableA-cont64 точно Hermitian и strict positive
definite.** Это сертификат stored inputs, не forward accuracy полногоΠ,
SCBA/native acceptance или Matsubara density. Старый Γ128 HOLD не снят.

Reviewer использовал только code/малые TOML/log и две witness-матрицы
из summary; source HDF5/Julia не запускал. Small stdlib audit
`global-gamma-exact-pd-result-audit.py/.json/.log` PASS, RLIMIT128MiB,
0.02s/17624KiB. Это независимая проверка execution/proof structure и
двух exact witnesses, **не второй полный HDF5 scan всех337744 blocks**.

## Fixed execution и источник

Summary SHA a703ec4c3e2feda1ae9571ab3d1d20b69dd9de48898f3d620ada008da2583335,
ALL input27da5a97ab61e9f290f5f8c7a2979d5873ec55e81f45e6b336fba99c495711d4,
runner d2395572bd724f9493be4d5c5cfb92fe9859e609719cafa0bced68f7de02934f,
exact math c6b14339d05720f4d55152b35edd9d43dca910507bfd16c8e750a41cdf02b6a9
independently совпали с actual files. Все9 source bindings равны
root source-digests,7 малых source SHA independently пересчитаны.
Γ file cab81029eff5cfccf8a71eea41d1571c16902b6a22a7bc81fb08a0e700e8db55,
model700ae17684d9ae52623ca6cfb7e07432ab1f43eb5151748061da53681b686540
принадлежат unprojected accepted-input map64; это не rawΦ64/Γbest,
не projected branch и не Γ128. Whole-file checks перед HDF и producer
checkpoint/input/ledger guards рассмотрены в static review.

Actual mode=all,execution=complete, requested=checked=30704×11=337744;
pd337744, nonpositive/nonHermitian/nonfinite0, failure witnesses empty,
exception field отсутствует. Exact component divrem count20,264,640=
60×337744. Это согласуется с полным traversal reviewed runner и
проверенным Bareiss/Sylvester positive route. ALL flag=true не основан
на selected fixture/benchmark flag. Benchmark имеет ровно5000 checked
и ALL=false; его counts не прибавлялись к ALL.

## Независимые exact witnesses

Julia TOML Float decimals roundtrip переведены в Python Fraction без
complex-Float промежуточных операций. Для first(e1,k1) и сохранённого
minimum-Δ5(e382,k11) отдельно выполнены exact rational LDL†, проверены
все10 positive pivots, все10 published integer minors с2^(k·emin)
scale и exact reconstruction всех50 matrix entries. Exact Hermiticity
обоих blocks также подтверждена. Это другой алгоритм относительно
production Gaussian-integer Bareiss.

Minimum recorded determinant имеет log10≈−196.6164 в stored native
units; точные integer minors/exponent−247 — первичный witness.
Он строго положителен, а не округлённый ноль, но **не** minimum
eigenvalue/condition bound или quantified distance всего source до
PSD boundary. Минимальность по всем nodes опирается на reviewed exact
dyadic comparison/owning worst cache и full execution, не на повторный
reviewer scan. Maximum integer determinant length641bits.

## Цена и научная область

Прочитан root ALL log:exit0,17.45s process/16.14106s script,
1180640KiB≈1.126GiB,Swap0,J4,CPU average208%. Admission4GiB/High3.5GiB,
CPU400%,300s/internal270s/output32MiB совпадает с отдельным ALL input.
Предыдущий fixture log:51+4+309=364 assertions PASS, включая44 старых
realistic matrices; они проверяют алгоритм, не новые Γ64 nodes. Benchmark
2.178s и его rough147.1s extrapolation были cost screen, а не guarantee;
actual17.45s не следует переименовывать в warmed speedup.

Exact PD всех frozen Γ64 nodes с положительными hats позволяет отдельно
построить PSD continuum данного интерполянта. Для intended exact
Hermitian all9 hH+Cauchy self-energy это существенная causal prerequisite;
matched h/Γ/T/ghost/source model и source real-axis convention всё ещё
нужно явно связать с future Mats run. PD одного input не доказывает
floating resolvent accuracy, absence/size real atoms, tail/window/grid
adequacy, correct nonlinear solution или experimental validation.
Γ64 pilot остаётся несошедшимся; ALL certificate не принимает его SCBA.
