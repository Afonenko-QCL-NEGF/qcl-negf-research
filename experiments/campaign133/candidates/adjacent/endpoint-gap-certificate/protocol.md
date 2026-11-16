# Frozen Γ64 endpoint/ghost enclosure: write-only protocol

Цель — upgrade существующих V4 floating witnesses до консервативных lower
gaps/upper ghost mass. Ни actual Γ, ни HDF5 автор не читал; Julia/Pkg/tiny/
postprocessing не запускались. Не относится к новому Matsubara all5 SCBA
root без этого Γ source. Отдельная numerical identity, не SCBA acceptance.

Источники зафиксированы: V4 script814dff50…, actual V4 inputf6560732…,
actual completed summary323cc177…, canonical Γ64 checkpointcab81029…,
terminal6eb3ca27…, exact saved hH9748cff3…. Используется **существующий V4**
`source_frame`/`require_all_PD_certificate`: whole337744 ALL exact-PD binding,
terminal=checkpoint entrywise, logical array hashes, saved exact Hermitian hH
all9, Tminus=Tplus†, native A/μ comparator только для прежнего input binding.
Одна Γ array остаётся в памяти; нет новых GR/A/Φ/occupation/μ-fit/SCBA maps.
Source V4 включается только как definitions; его `main` не вызывается.

`input.toml` false. Будущий authorized copy может отличаться только
`execution_authorized=true`; paths привязаны к frozen каталогу. До actual:
independent source/math review, затем root tiny (`tiny_checks.jl --tiny`)
в отдельном малом процессе. Сам tiny не читает actual H5: basic interval
primitives/exact rational constants/independent Big256 coefficient formula,
complex Gershgorin sign, triangle sliver area, mass factor и negative-gap
rejection. Big256 здесь тестовый witness, не runtime certificate/eigensolve.

Actual один процесс: J1/BLAS1, MemoryMax2GiB2147483648B/Swap0/CPU<=100%,
external120s (internal START до imports), небольшой output<4MiB по fixed
shape; root может поставить LimitFSIZE64MiB. Peak/time не измерены. Γ loader
унаследован от actual V4 (~917MiB process peak там), но новый peak не обещан.
На memory/time/input/certificate/overflow/domain fail исходы сохраняются
как failed/unavailable, не pass. Не добавлять новый heavy fallback.

Сертифицируемая точная модель — V4 full equal-width hats с **stored Float64
canonical centres** как exact dyadic cj и dyadic h. Это учитывает1ulp centre
rounding; не заменяет их молча ideal exact `a+(j−1)h`. Native window[A,B] —
реальные saved grid endpoints. Adjacent support overlap проверяется outward;
точные supportℓ=c1−h,r=cN+h и areas вне[A,B] считаются rational. Tiny overlap
второй/предпоследней hat не опускается. Возможный nodewise mismatch rounded
hat superposition с native ideal-Toeplitz nodal closure этим не ограничивается;
его нельзя скрыто называть quadrature error.

Численный proof domain: exact integer/rational операции GMP; basic Float64
IEEE +,−,×,÷ в RoundNearest, gradual underflow (`get_zero_subnormals=false`),
без fastmath/BLAS/FMA reassociation в enclosure pass, finite intervals.
Каждая scalar operation расширена prevfloat/nextfloat. Поддержка IEEE
subnormals/nextfloat описана в [официальных Julia numbers docs](https://docs.julialang.org/en/v1/base/numbers/)
и [floating-point manual](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/).
Это auditable conditional numerical certificate, не formal verification
CPU/compiler/HDF5 implementation. Любые nonfinite intermediates отвергаются.

π bounding не доверяет libm constant: `tan(4atan(1/5)−atan(1/239))=1` и
диапазон угла дают Machin identity;32 alternating rational terms ограничивают
оба atan. `log2=2atanh(1/3)` ограничен40 положительными rational terms плюс
geometric tail. Dyadic bounds проверяются **exact rational comparisons**.
Для t>1 коэффициент ряда positive и tail после24 terms ограничен
`t^−49/[π·49·50·(1−t^−2)]`; t=1 использует log2/π. Нет log/cancellation FFT
на endpoint. M берётся более консервативный1/π>asinh1/π.

Выход `lower_gaps` имеет порядок exterior_left/right,strip_left/right.
Positive full45 Gershgorin lower bounds доказывают нужные Loewner signs в
данной exact-input модели; никакого Float eigenvalue acceptance нет. При
nonpositive bound — unavailable, без |gap|/floor. Bound для ghost mass
использует traceΓ>=||Γ||2 и exact area, поэтому может быть слабее прежнего
floating Frobenius witness, оставаясь upper bound. В per-k mass factor9 нет;
sheet bound включает gs*Nb*wk, без второго2π. Missing FD number≤mass bound
поскольку0≤f≤1; μ фиксирована.

H5 сохраняет interval lower/upper endpoint matrices, q-order, support bounds,
native window. TOML сохраняет все44 signed lower gaps, conditional statuses,
per-k upperβ, optional sheet bound/target, π/k1 bounds, source hashes/units/PD
binding. Архивируются новый script/input и использованные V4 source/input/
summary. Ошибка после output allocation попадает в summary; ранняя source
guard ошибка — в root captured log. No data overwrite.

Если bounds положительны, exterior/endpoint atoms и ghost mass данного
frozen model исключаются/ограничиваются. Чтобы **строго** приписать остаток
≈9.3377e−5target real-axis quadrature, всё ещё нужны Mats tail+arithmetic и
native spectral evaluation/geometry bounds: `N512−N1024≈4.4e−8target` — лишь
witness. Этот diagnostic отделяет missing spectral mass, не доказывает
SCBA/grid/экспериментальную валидацию и не заменяет новый root computation.
