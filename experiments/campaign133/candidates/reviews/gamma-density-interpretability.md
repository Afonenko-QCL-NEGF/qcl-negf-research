# Γ128 Mats: prerequisite и область интерпретации

**Рекомендация: отложить actual Mats run, если его вопрос — certified
continuous number той же source state.** Freeze-IO полезен и может
продолжаться. Нынешний V3 разрешает определить и посчитать только
tail-corrected Matsubara functional. Его отождествление с real-contour
spectral density пока условно, а strict same-operator relation нарушена
raw-H policy. Это не диагноз noncausality фактического Γ128.

Нужная предпосылка для spectral identity: H=H† и Hermitian Γ(ε)≽0
в точной определённой модели. Для положительных PI hats
Σ(z)=∫Γ(ε)/(2π(z−ε))dε при y=Imz>0 имеем ImΣ(z)≼0 и
Im[zI−H−Σ(z)]≽yI. Значит upper-half-plane zeros отсутствуют,
‖G(z)‖≤1/y и действует Schwarz reflection. Именно эти свойства
позволяют заменить real contour Matsubara sum, с отдельным учётом
конечных windows/ghosts/atoms и normalization.

Γ128 cone report допускает λmin≈−1e−18 в floating backward allowance.
Такой результат не доказывает отрицательность **точного** spectrum,
но и не сертифицирует PSD точного source matrix. Small global norm
не исключает zeros произвольно близко к real axis. Простой аналитический
контрпример: H=μ+iη,Γ=0,0<η<πkBT. Тогда G(z)=1/(z−μ−iη)
имеет upper pole при любом η>0. На всех sampled Mats points он хорошо
определён, а V3 functional1/2+(2/β)ΣHermG(μ+iωn) равен1/2.
Real-axis A=−2ImG — отрицательная Lorentzian, и её полный Fermi integral
равен−1/2. Разность1 не исчезает при η→0. Это пример недостаточности
малой anti-Hermitian нормы, не утверждение о таком pole в Γ128.

Итак, backward error LU, bounded kernel/tail remainder, Σ-sign tests на
нескольких frequencies и512→1024 agreement контролируют sampled series,
но не пустоту strip0<Imz<πkBT. Upper-pole residues нельзя исключить
или оценить их вклад в density по одним этим witnesses. Операция Herm(G)
в summand также не восстанавливает missing Schwarz identity raw model.

Есть отдельная policy mismatch. V3 `chain_hamiltonian:104` копирует
raw h во все9 cells и лишь допускает norm(H−H†)≤1e−11. Native e171
`production_operators.jl:300` формирует hlocal=(h+h†)/2 для cavity sides;
central Dyson использует raw h (source writer `_retarded_state:254`).
Поэтому V3 all-raw9 и source raw-central/Herm8 не тождественны даже
до вопроса PSD. Их малая разность может быть корректным соседним
operator control, но не exact same-source density identity.

Для заявленного continuous-number вывода нужны либо qualified causal
source/model и matched operator, либо отдельное достаточное доказательство
analyticity/no upper poles и допустимости contour deformation для exact
используемой модели. Hermitian H и PSD-defined Γ с positive hats —
прямой достаточный путь. Если для этого вводятся Herm projection/PSD
construction, это **новый frozen operator с новой identity**; нельзя
молча сравнить его с прежним saved A как exact same-input oracle.
Matched reference A при прежнем source μ должен относиться к тому же
новому operator; старый A остаётся отдельным comparison control.

Даже exact causal source не делает finite nodal fA integral равным
full continuous number без quadrature/window/atom/ghost qualification.
Эти ограничения прежнего backend review сохраняются. При отдельном
разрешении raw V3 run допустим с label «Matsubara functional of raw
frozen operator; spectral-density identity not certified», без density
acceptance и attribution mismatch quadrature. Сейчас такой run не
устранит обозначенную prerequisite, поэтому разумно удержать actual.

Это уточнение относится и к прежнему FD64 Mats backend5249: его
`chain_hamiltonian` тоже использовал raw h во всех9 cells (тот же
побайтово сохранённый numerical prefix). Прежний N1024 — действительный
**вычисленный functional** с сохранённым tail/refinement witness;
continuous-density interpretation остаётся условной при exact
causality/Hermiticity assumptions. Даже computed PSD FD64 Γ сам по
себе не устраняет all-raw-H оговорку. Прежний результат не следует
публиковать как rigorous independent real-axis error certificate.
Это уточнение области вывода, а не обнаружение actual upper poles
или отмена корректно сохранённого вычисления. Исторические frozen
scripts/inputs/results не изменяются; publication может явно сослаться
на эту позднейшую qualification.

Только static source/math assessment. Reviewer не выполнял Julia/HDF5,
pole search, PSD certification или новые scientific вычисления.
