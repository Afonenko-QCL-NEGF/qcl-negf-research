# Independent static review: Matsubara M3 moments postprocess

**Вердикт:** frozen source пригоден для Root tiny qualification и затем
одного ограниченного source-bound postprocess при успешных tiny controls.
Научный результат/исполнение actual здесь **not_measured**; автор review
не запускал Julia/HDF5, не менял код или данные. Новая selfconsistent M3
closure корректно выделена отдельно от unchanged finite-N map.

Freeze: module `moment_tail.jl`
9182ae802fd26c41e6105bdcf2979bb9d2406d0174bd35b31b5ccfee793b80ca;
adapter3fe6428e700d5916637f98431be85ca5350264d7383ff65918299aab02656d41;
tiny4c76270df6fed02a57cf53ebefc865b6c91765a9d6ecf926c934b02a4d268c4e;
false inputc9b45119888c902d8fc4f9e2d539f790e961d8599a5e4f297daaffecfec8b82f.
Все6 SHA256SUMS entries проверены stdlib hashing; SHA text dependencies и
малых actual TOML совпали;7 ROOT source-digest entries совпали с input.
Большие HDF/source arrays этот review не читал и не перехешировал.

**Формулы/единицы.** `moments` включает ровно4 static+LO maps. `h` —
central block saved Hermitian H45 минус μI **один раз**. Оба упорядоченных
произведения hS0+S0h сохранены. Bare45 subtraction отменяет H45²/H45³
terms; дополнительного factor9 или Tplus term в interaction B нет.
Получаются S0=Σstatic K[I]+(2NB+1)KLO[I] и
B(ρ)=hS0+S0h+Σstatic K[h]+KLO[(2NB+1)h+Ω(I−2ρ)], согласованные с
`matsubara-all5/M3-density-derivation-review.md` и прямым разложением
thermal pole weights. Размерности: S0:E0², B:E0³, t4:E0⁻³, ρ dimensionless.

`tails` считает (2/β)Σtailω⁻²/⁻⁴, не забывая2/β; t4>0. BigFloat256
total-minus-prefix для128/256 дополнительно проверяется независимыми
положительными integral lower/upper bounds. Это numerical scalar witness,
не interval enclosure; огромные N этим finite prototype не принимаются.
При K≡0 S0=B0=0, L=I и correction0: сохранённая exact-bare-subtracted
ρN не получает второй раз bare tail. Classic M1 от saved G вычисляется
и сохраняется **отдельно**, с −t2h, без подмены исходной ρN.

`correct` сохраняет both one-shot ρN+t4B(ρN) и solution
(I+2Ωt4KLO)ρ=ρN+t4B0; знак слева **плюс**. Код строит general complex
275×275 superoperator на всех elementary complex matrices, решает LU,
вычисляет residual повторным K application. Hermitian/Pauli diagnostics
не проектируют результат. Не используется ложное CP⇒SPD и не решаются
normal equations. All4 ordered static/LO components не потеряны.

η=2Ωt4 max_k||KLO[I]_k||₂<1 — условная CP C*-norm/Neumann опора,
не измеренная строгая bound rounded operator. При η≥1 prototype останавливается;
cond(L), residual absolute/scale/relative и weighted equivalents сохраняются.
Weighted norm — sqrt(Σk wk||r_k||F²), без E-weight/spin. Zero/zero→0,
nonzero/zero→Inf, без floor. Bound ||r||block/(1−η) относится только к
linear moment equation при этих предпосылках/точной арифметике, не к
infinite-N density error или всему floating pipeline.

**Oracle/source.** Tiny scalar oracle использует независимые два pole/residue
для G=[z−ε−a/(z−b)]⁻¹: interaction third moment a(2ε+b), direct thermal
residue density, finite-N bare subtraction и sign/refinement M3. Tiny также
проверяет scalar native K fixture, analytic scalar linear/one-shot solution,
positive t4 bounds/zero policy/weighted residual normalization. Он ещё
не выполнен в этом review и не является whole physical model oracle.

Input связывает N128/N256 latest HDF SHA61cebaa7…/fd4752c1… с actual ALL5
32 completed maps, raw Σ=Φ0, μ/frame/iteration16/ω/grid/ρ-number. ROOT
digest JSON согласован; whole source SHA runtime проверяются до HDF reads.
Adapter reused old94c6 gate/source проверяет Coree171/tree/accepted env,
saved hH/K/qK/units и source-binding. Читаются только малые saved states;
source model reconstruction разрешена и не названа новым Φ. Вызовов
old run_actual/pilot, Dyson/full fixed-point map, μ-fit нет. Model K
используется лишь как linear moment map; имеющиеся G/ρ не перезаписываются.

**Scope/бюджет.** One-shot соответствует leading high-frequency coefficient
с savedρN. Selfconsistentρ подставляет corrected density в этот coefficient
и тем самым вводит **новую asymptotic moment closure**, не exact identity
старого finite-N reference-tail map. Их O(t4²) различие — асимптотическое
утверждение без certified constant. Ни small linear residual, ни близость
corrected128/256 не доказывают infinite-N SCBA, omitted LO feedback в low n,
higher-moment remainder, правильную fixed density/μ или Poisson/native pass.

План Root: tiny сначала; actual после tiny ≤120s/2GiB/Swap0/J1-BLAS1/
CPU≤200%, file и суммарный output≤32MiB, отдельный output. Frozen input
runtime_authorized=false; допускается отдельная actual copy с изменением
только этого флага. Newinput equality/source gates и external OS caps
нужны совместно: внутренний tick не прерывает отдельный LU/import сегмент.
Предварительно найденный пробел weighted residual evidence устранён до
этого freeze. Дополнительных blockers в зафиксированной версии не найдено.
