# Независимое static scientific review Matsubara ALL5 (V1)

**Решение: математика допускает bounded prototype; actual HOLD до минимального raw-preservation fix, actual-module scalar tiny и review дельты.** Источники не менял; Julia/Pkg/HDF5/compute не запускал. Root запускает qualification и source-bound runtime отдельно.

Reviewed V1: `adjacent/matsubara-all5/matsubara_all5.jl`, SHA **7e847d0f8dfd048b1eee18472bd43f1629127cf68edb5ab5024c2335da6c84b0** (191LOC; автор сохраняет unexecuted-v1). feasibility SHA9b998591…, usage SHAcd8cb504…. Native bindings просмотрены на accepted Core **e171caa60cce676c6c4238dd671e7c527f2bd977**. Это новая fixed-μ equilibrium direct-static C + LO Fock closure, не исходный finite-energy PI/root density-fit solve.

## Findings перед actual

**P1 — raw map теряется при штатном soft stop.** `pilot:142–143` возвращает только `S/history`; `S` уже изменён mixing в конце предыдущего шага и не имеет возвращённых согласованных `G/Phi/rho/w`. `callback(rec):150` передаёт только scalars. Root не может сохранять full raw candidate через API во время расчёта, а OS stop также может оставить только history. Исправление узкое: callback получает coherent snapshot `(S,G,Phi,rho,w,record)` **до guard/mixing**, update S без mutation сохранённого snapshot, soft stop возвращает last coherent evaluated state. Root adapter сохраняет каждый малый raw map/atomic latest (N≤256 arrays малы). Это не критерий сходимости и не смена физики. Автор уведомлён и согласился исправить.

**P2 — скрытый eps floor в нормировках.** `qualify_native_ops:65`, `pilot:148`, `onelevel_oracle:177`, `compare_cutoffs:186–187` используют `max(norm,eps)`. Primary raw relative при околонулевом Phi тогда зависит от произвольной абсолютной единицы. Сейчас ложной acceptance нет (`status` всегда budget/noncausal), но доказательное сравнение требует `rel(a,b)=b==0 ? (a==0 ?0:Inf) : a/b`, absolute+scale отдельно. Для roundoff diagnostic causality observed-scale floor — отдельный явно объявленный backward allowance, не tolerance raw SCBA. Автор принял исправление relative policy.

**Qualification gap — standalone oracle не упражняет реальный module wiring.** `onelevel_oracle` корректно сравнивает независимую raw Mats sum с рациональным ответом, но не вызывает `reference`, `lo_delta`, `greens`, `density`. До ALL5 нужен минимальный actual-module scalar context: nb=nk=1, Tplus=0, девять одинаковых uncoupled cells, static0/LOg²; referenceLO и density против независимого one-level pole/Fermi, plus **nonzero ΔG** и direct symmetric sum для отрицательных частот/β packing. Это защитит связи, не дублирует actual native six-index17-slice check. Автор сообщил подготовку такого tiny.

P2 scope: `context:12` требует positive wk; в actual Core annular weights положительны даже при κ0=0 (`reference/grids.jl:31–39`), поэтому здесь не найдено incompatibility. Root source gate должен отдельно фиксировать finite beta/Omega/mu, positive qK и original all5 kernels/hH/Tplus/geometry; в модуле нет provenance/physical-source loader и это явно вынесено Root.

## Проверенная алгебра и code bindings

**LO знак и β.** При D+(τ)=(NB+1)e^−Ωτ+NB e^Ωτ и Gε(τ)=−(1−fε)e^−ετ конвенция Στ=D+ KLO[Gτ]. Fourier D+=2Ω/(ν²+Ω²), поэтому положительный prefactor1/β в `lo_delta:91–97` корректен. Bare pole даёт

`ΣLO(z)=g²[(NB+1−fε)/(z−ε−Ω)+(NB+fε)/(z−ε+Ω)]`.

Независимая связь с native lesser/greater: из `reference/kernels.jl:439–459` ΓLO at E=ε+Ω имеет коэффициент `(NB+1)(1−f)+NB f=NB+1−f`; at E=ε−Ω коэффициент `(NB+1)f+NB(1−f)=NB+f`. Следовательно `reference:79` и standalone oracle:171 не перепутали emission/absorption или знак ImΣ. Thermal bath undressed; Hartree/tadpole отсутствует в обоих сравниваемых Fock вариантах.

Первичный внешний source действительно использует Mats Bubble2 `−D_M G_M`; при D_M=−D+ это тот же знак. [NESSi Bubble2 exact commit,83–95](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/libcntr/cntr/cntr_bubble_impl.hpp#L83). Fixed-bath uMig берёт D0 в Fock bubble, без phonon Dyson. [NESSi uMig,174–190](https://github.com/nessi-cntr/nessi/blob/6bb4c832a209d2c079fdbbdbc5c05ab8d5503aa6/examples/programs/Holstein_impurity_impl.cpp#L174). Оба primary source просмотрены; их conventions не подменяют приведённый native sign derivation.

**Negative fermion packing.** В `lo_delta` j∈[−N,N−1], G_(−n−1)=G_n†, negative j index−j; difference ωn−ωj=2π(n−j)/β. В reshape energy-frequency axis первый в Julia, обе стороны совпадают. Ошибки transpose/conj/one-index не найдены статически; nonzero-ΔG tiny всё равно необходим.

**Units/κ/static contraction.** h/μ/Ω/ω/Σ должны использовать единицу E0; β_native=E0/(kBT), μ_native=0.0084596eV/E0, Ω_native=0.0367eV/E0; δμ второй раз в h не добавлять. `_production_kernel` и `_production_static_contraction!` получают нормированные K и тот же qK/wk. Independent `literalK:42–50` совпадает с Core `reference/kernels.jl:418–436`: K[m,m′,a,c,d,b]X[m′,c,d] **без conjugate X**, qK один раз. Static Σ=K[G] — direct Green C, не finite-window Hilbert baseline. qK сохраняет native scaled coupling; дополнительного Bose factor в static channels нет.

**9-cell/45×45 geometry.** `context:16–28` H Hermitian, hH повторён all9, upper coupling Tplus/lower Tplus†, centre rows4nb+1:5nb. Это соответствует V=0 native finite4 recursion: left coupling adjoint(Tplus), right Tplus (`production_operators.jl:361–385`). `greens` повторяет **одну central-driven Σ** во всех9cells, не cellwise SCBA. Central residues V atom V† имеют сумму I_nb, exact bare reference включает все45 atoms, не только5 central eigenvalues. Никаких extra factor9.

**Density.** ρ=ρref+(2/β)Σn≥0 HermΔG_n; заполняется Mats reference f(H45−μ) centre, где common μ сдвинут один раз. Number=gsΣk wk Trρ. Дополнительный1/(2π) здесь был бы ошибкой: он уже вошёл через energy spectral integral →ρ, native wk содержит annular k measure1/(2π). Core seed N=gs/(2π)∫dEΣwk Af даёт то же. Pauli eigenvalues лишь диагностика; clipping/normalization нет.

## Causality, tail и предел научного вывода

Согласованный Hermitian full45 и ImΣ_n≼0 дают ImG_n≼0 на положительных Mats nodes. Exact bare Gref — Stieltjes/Herglotz resolvent, но **reference + arbitrary finite ΔG не является автоматически global causal analytic completion**. При CP kernels, conjugate frequency packing и causal samples есть более узкий nodewise argument: paired LO weights D(ωn−ωj)−D(ωn+ωj)≥0 для n,j≥0. Симметричная inside correction эквивалентна inside G + outside bare Gref, каждый pair сохраняет знак. Это conditional Mats-node property, не analytic continuation/real-axis A certification. Actual pointwise `causal_max` должен сохраняться с measured scales; global512eps*maxabs allowance может скрывать ошибки очень слабых blocks, поэтому не является строгим PSD/global theorem.

При CP L=||ΣstaticK[I]+(2NB+1)KLO[I]||sup-k, ||R45||≤1/ω и ||diag9 δΣ||=||δΣ||, поэтому conditional L/ω0² без9 корректен. Если ratio≥1, contraction не доказан; если CP/native operator source не подтверждён, scalar bound вообще условен. Не применять общую MDE theorem к central-driven closure автоматически.

При Σ∞=0 bare/full moments I/z,H/z² совпадают и ΔG=O(z^-3); HermΔG=O(ω^-4) при Hermitian moments. Density residual tailO(N^-3) — asymptotic **без конечной certified constant**. LO reference correction оставляет omitted ΔG tail, не exact infinite SCBA. N128/256 различия на первых32ω иρ смешивают cutoff и nonlinear error, пока rawΣ residual заметен; 16 updates каждой сетки не гарантируют inner convergence. Fixedμ не решает dopant Ntarget и не свидетельствует native public pass, Poisson, density constraint или experimental validation. Малый Dyson backward — качество линейного inverse, не SCBA residual.

Итог: signs/units/local45/κ-reference construction научно согласованы; **V1 HOLD по raw preservation**. После узкой V2 дельты, actual-module tiny, source provenance и Root resource admission разумен только bounded fixedμ exploration128/256≤16 each. При budget status сохранять raw и not_converged; acceptance/posterior cutoff conclusion не выводить автоматически.

## V2 delta review — module APPROVE, conditional tiny/source/runtime

Module **1434e86fe452f02c5831eaeb0fc83dc1dc12710dcc3d7e81f7f36cd985a31447**, 218LOC, просмотрен diff к pinned V1. Математический reference/lo_delta/greens/density/κ-код не изменён.

P1 закрыт в API: `latest=(S,G,Phi,rho,w,record)` передаётся callback до causal guard/mixing, новый S создаётся без mutation, soft return содержит последний согласованный evaluated snapshot; initial soft stop до первого map честно обозначен unaligned. Root adapter обязан действительно записывать callback snapshots и не вызывать `compare_cutoffs` для unaligned/missing arrays. OS total600s/32 updates всё ещё внешние; сам модуль их суммарно не enforce.

P2 relative floor закрыт явным zero policy. Roundoff scale allowance остался отдельным диагностическим механизмом. Raw norm scale желательно записать scalar field явно, но он восстанавливается из сохраняемого Phi и не блокирует qualification.

`wired_onelevel` задаёт independent auxiliary pole self-energy `a/(z−eb)` при Tplus0; actual full9 solver должен давать `1/(z−e−a/(z−eb))`, с двумя poles `(e+eb±sqrt((e−eb)²+4a))/2` и residues `(lp−eb)/(lp−lm), 1−wp`. Эти формулы и знак/residues в fixture корректны. Actual `reference` проверяется bare onepole expression, actual `lo_delta` имеет nonzero delta и сверяется с суммой analytic Fock poles; density сверяется с weighted Fermi poles. Это закрывает wiring gap в source; **результат tiny ещё not_measured автором/рецензентом** и принадлежит Root execution.

Решение V2: **APPROVE module source для Root tiny и затем bounded fixedμ pilot при успешном actual tiny, independent source-bound adapter review/admission**. Новых mathematical blockers нет. Это не acceptance inner/cutoff/native problem, не causal real-axis continuation и не разрешение неограниченного μ-loop. Отдельный новый runtime adapter здесь не проверен, пока он under construction. Source/arrays не редактировал, Julia/HDF5/solver не запускал.

## Runtime adapter/source review — APPROVE bounded actual on final freeze

Reviewed adapter **94c6f2d48851e9b894ce25f6a353d66accbc34e2ebb269ee05b040f4cd717c1e** (187LOC), input **d920768e48aef32a9af3bba5c2778b7f7cd6659fc15981a72b39d3bd65b4b85e**, module1434 unchanged. Superseded9b18ac had an actual source-key blocker: parent8d5c… contains no `voltage_mV`. Автор/рецензент независимо нашли это до запуска; final checks actual `iszero(_scaled_physics(...).Eᵖ)` после pinned builder вместо несуществующего input key. Missing-key blocker закрыт; physics default по отсутствующему key не выдуман. Final Tminus/Tplus adjoint сравнивается elementwise `==` (signed-zero differences не физическая non-Hermiticity).

Полный static source path просмотрен: exact accepted project/manifest/Core HEAD/tree/clean; GP50bc и три included definitions; source model/input/identity/summary/map64/spectral ledger/selected; V4 helper192e/digestb121. Byte SHA каждой перечисленной definition/text/source file и текущих accepted Project/Manifest дополнительно проверены **streaming bytes only**. Собственные HDF dataset reads/Julia отсутствуют. Фактическое source input совпадает с pinned8d5c. Runtime HDF/model-array binding остаётся Root execution, не объявлен выполненным этим review.

Frame корректен: raw native model `verify_problem_model` перед saved hH; exact logical hH digest, hraw/Tplus/axes/E0 checks; E0=.1eV/L0=1e−8m/spin2; physical V0 и thermal LO; β=1/ctx.kBT, Ω=scaled ħωLO. Fixed μ — точное map64 0.08459614262982018 native =0.008459614262982018eV. Не загружаются Γ/Φ source arrays или PD certificate; saved hH source identity не объявляет новый Mats point принятым. Core CP contraction ещё проверяется17 non-Hermitian frequencies independent loop на actual kernels до iterations.

Limits/сохранение: top-level PROCESS_STARTED до imports плюс внешняя OS wall600, no >32 structural requested updates (два N по16), каждый независимо Σ0, α1. Callback до guards сохраняет per-map inputΣ/G/ρ и scalar, atomic latest additionally rawΦ/w; только named latest own-output replace разрешён. При α1 rawΦ_j совпадает с последующим stored inputΣ_(j+1) по значениям, а lastΦ каждогоN сохранён latest; signed-zero bit pattern не обещан. Admission учитывает old latest+temp. Payload≈56.4MiB before headers плюс biggest latest temp≈3.2MiB; existing64MiB gate проверяется до/после каждой записи. Global kill может оставить tmp; request ledger честно ограничивает started/completed, не выдаёт unknown за0. Native_public_calls=0, μ-root не запускается. Полное raw/source сохранение даёт воспроизводимый bounded negative outcome.

P2 resource caveat: adapter gate читает memory.max/swap.max/CPU/file ceiling, **не читает memory.high**, хотя input задаёт1.75GiB; Root должен явно установить High1879048192 внешним scope. Это не blocker при указанном Root enforcement, но script сам High не подтверждает. Source GP build_context memory estimate сохраняет старый native envelope; это model reconstruction/cache build, не admission native solve под2GiB. Фактический RSS этой части ещё not_measured; bounded refusal остаётся допустимым отрицательным исходом без cap increase.

Root сообщил/лог tiny-root.log просмотрен: **22/22 PASS**,4.80s/maxRSS491948KiB/Swap0/exit0. Actual wired Grel≈1.26e−16, refLO0; ΔG LO abs128→256:3.061e−9→4.005e−10, density2.87e−9→3.59e−10; snapshot alignment/immutability и alpha1 inheritance проверены. Это actual-module qualification, не ALL5 result.

**Окончательное решение: APPROVE final adapter94c6/inputd920/module1434 для одного bounded actual128/256≤16each при Root OS2GiB/High1.75GiB/600s/J1BLAS1/Swap0/file32/data64MiB.** Source arrays/source gate должны реально пройти. При convergence/tail/causality сомнениях сохранять unresolved; новые μ/N/iterations не разрешаются этим review. Реализацию не менял, расчёты не запускал.
