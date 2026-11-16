**Независимый narrow result review stable-A cont64**

Факт: continuation завершилась по cap64, не по convergence. Raw Φ−Γ
на last/best call64: relative4.125583337831235e−6,
weighted3.412502855963494e−6; target1e−9 не достигнут, native0/
scientific_acceptance=false корректны. First boundary whole rawΦ replay
имеет saved absolute/weighted difference0; девять scalar values
независимо сверены bitwise по малым TOML, все совпадают с parent map16.

Accounting проверен раздельно: **64 entered Π =62 completed maps+2
q-failures (calls5,23)**. Всего65 attempts: accepted seed1/linear4/
AAII54/convex fallback3 и3 rejected AAII. Третий reject — до входа в Π
при counter28: Γ minimum−2.5592666e−13 (e8832,k1), хуже budget7.1054e−14;
Π не расходовал. После каждого отказа выполнен один принятый fallback.
Map gaps ровно5/23, spectral ledgers64, domain records3; счётчики не
путают trial и фактическую карту.

Оба entered failures — **finite отрицательные q**, не NaN: call5
qmin−1.11246e−37 (e394), call23−2.20263e−38 (e533), по одному entry.
Relative negative masses2.4721e−36 и4.8947e−37 target. Guard сработал
без floor/clipping; generic covariance не гарантирует допустимость всех
AA trials. Это не доказательство причины отрицательности каждого trial
без дополнительных selected input/precision checks.

На62 принятых картах computed inputΓ/rawΦ minima положительны
(минимумы6.2020e−41/7.8287e−41), но это floating witnesses. Last qcov
positive minimum1.3284e−43; same-G qsub содержит64 negative entries,
relative negative mass1.7601e−35. Область надёжной матричной арифметики
улучшена, strict global PSD/future-map accuracy не доказаны.

Прежний HermΛ/raw-central-h cont128 с AA cap10 имел final/best
7.192305218634285e−6 (weighted5.95755e−6); новый final составляет0.57361
от него при другом call budget. **Нельзя приписать разницу одной A**:
h-политика/AA cap и restarted trajectory различаются. В этом cont64 уже
приняты6 AA proposals с L1>10 (calls8,26,44,62,63,64), максимум29.4372;
поэтому прежний cap10 действительно ограничивал бы часть предложений.
Исходный stableA16, напротив, повторял Herm16 trajectory на machine level.

Execution519.6968s script/528.69s process, peak7616820KiB, Swap0,
data2512015759B: в объявленных hard12.125GiB/high11.5GiB/overall1800s/
data8GiB/file5GiB. Native не вызывался из-за отсутствия reduced trigger.
Новый checkpoint остаётся unresolved reduced state, а не принятым SCBA.

Reviewer проверил только small TOML/log, counters, ratios и parent
comparison под128MiB Python; Julia/HDF5/full fields не открывал.
`reviews/stableA-cont64-result-audit.json` сохраняет exact accounting,
failure records, accepted trajectory extrema, provenance и AA-cap evidence.
Ни elapsed improvement, ни two unresolved outcomes не доказывают target
root; следующий trial-admission вариант должен сохранить raw Φ/gates и
явно отличать изменение AA proposal от исправления исходных уравнений.
