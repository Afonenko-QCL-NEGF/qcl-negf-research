# Positive-source operator: отдельный qualification plan

Этот delta-документ уточняет frozen `../stable-spectral-plan.md`; опубликованная исходная записка сохраняется неизменной. Только proof/design, без implementation/runtime.


Resource identity следующего proposed pilot: actual MemoryMax13019119616B
(12.125GiB), Swap0, Julia4/BLAS1/CPU400%, pilot600s/overall1800s.
Прежние12.375GiB/180s из раннего proposal не являются текущим ограничением.

До изменения accelerator выбрать **один** original positive FD128-derived
regridded Γ seed, без project/clipping стартового Γ. Минимальная arithmetic
ветвь для operator-only проверки:

1. hH=(hraw+hraw†)/2 вычислить в Float64, как native plan.hlocal, и
   использовать **одинаково в central и side cells**. T+/T+† unchanged.
2. Local uniform PI Λraw → Herm(Λraw), correction norm/locator сохранить.
   При mathematical Hermitian Γ exact PI real part Hermitian; операция
   восстанавливает это representation свойство, а не добавляет scattering.
3. Native all-component cavity один раз с localR=ΛH−iΓ/2,
   local<=iΓ, local>=0. Central native GR и physical-source Γemb дать
   Acov=GR(Γ+Γemb)GR†; q/FD/scattering потребляют именно это A.
   Пока **не** переписывать retarded Schur recursion layer-by-layer:
   сначала измерить её ΓR/ΓK discrepancy и error относительно small Big oracle.
4. Не project Φ. После retarded qualification допустим один fresh raw Π
   на том же input: current μ/capacity, per-channel Φ/PSD/KMS/raw norms.
   AA/projection отсутствуют. Если Phi negative, различить positive-source
   kernel arithmetic/CP deficiency, а не скрыть negative output.

Exact model этой ветви: all9 Hermitian hH, same positive-hat analytic Σ,
same Γ and native finite4, следовательно causal A=source congruence по
доказательству выше. В Double GR может иметь retarded/covariance mismatch;
positive Gram q не доказывает accurate GR. Фактический mismatch должен
быть qualified, прежде чем объявить ветвь исправленным operator.

Herm(hraw) сохраняет Hermitian part и является identity в точной underlying
Hamiltonian theory. Но **literal floating raw central h изменяется**.
Relative correction eps не доказывает harmless effect возле Γ1e−27:
для delta h bound на GR содержит g²||delta h||/(1−g||delta h||).
Raw h/hH, Λraw/ΛH, ΣR assembled, source/cavity ΓR/ΓK, GR/Acov/Asub,
number weights и precision/conditioning witnesses должны остаться в ledger.

В selected snapshot особенно полезно сохранить **preassembly ΛH** до
`−iΓ/2`, рядом со stored Float64 ΣR. Tiny Γ off-diagonal может потеряться
уже в representation сложения с большой coherent Λ. Big intended-model
oracle строить из promoted split hH/ΛH/Γ/T, не только promoting ΣR, в которой
эта часть могла исчезнуть. Реконструкция Λ=ΣR+iΓ/2 остаётся отдельным
witness и не заменяет original preassembly buffer. Достаточно нескольких
заранее объявленных5×5 blocks, не нового full Λ output.

Acceptance зависит от заявленного результата: если утверждается прежний
unchanged native fixed point, unprojected final native audit должен
использовать его raw-h equations. Для явного нового hH/positive-source
arithmetic variant возможна scientific acceptance согласованных
variant-equations с прежними physical tolerances; old raw-h residual тогда
является comparability diagnostic. Нельзя выдавать native API, которая
внутри rebuild hraw, за hH audit или переименовывать её fail. Operator-only
qualification сама не принимает SCBA ни по одному из этих контрактов.
