# Cauchy tree prototype: independent static math/code review

**Решение:** в authoritative module `a74ba641…` blocking formula/index/
sign error не найден. Prototype согласован с ограниченным transform-only
вопросом; он не предъявляет formal whole-program error/causality certificate.
Math review не разрешает actual source probe/новые вычисления сам по себе.
Julia/Pkg/HDF5 source/Green/Γ runtime рецензент не запускал, код не менял.

## Frozen source и независимое runtime evidence

Все четыре записи `SHA256SUMS` проверены, OK:

| Файл | SHA256 |
| --- | --- |
| CauchyPanelTree.jl | `a74ba641e403ddc60fe68a2c9a481b878c7a4d32ba8e3de30ca696b7c805257d` |
| tiny_tests.jl | `3be8f0d8b1237d025b7976b732cc6d60eededbb6008be86c052a12562c58d62f` |
| source_selected_k.jl | `48df6f1224875278aa4c6b66820a327e171337a236681484c47eb74a6cacb321` |
| protocol.md | `a24298cc8b9cfb992ced2c5f3a191634c3bd8b41502ec49a97fad2ac70b8b865` |

Старая preexecution `5ce…` не использовалась вместо a74. Parent сообщил
и я прочитал `tiny-root-v1.log`, SHA256
`698ed19a4ff8a78167e977e76e080148f348d56537116f46c146665a65d03ffc`:
67/67, exit0, wall6.41s, peak405864KiB, swap0. Это actual **Root** tiny
evidence, не мой запуск и не actual Γ/source/density/SCBA result.
Native reference source accepted HEAD `e171caa60cce676c6c4238dd671e7c527f2bd977`,
tree `742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8` прочитан для geometry scope.

## Проверенные математические части

1. **Leaf moments, lines24–43/84.** Gauss13 имеет degree25 exactness при
   exact nodes/weights. Для pmax24 highest moment q23 умножается на
   linear Γ, integrand degree24: rule достаточна. Фактор width*w/(4π)
   правильно включает affine map [−1,1]→panel и 1/(2π). `a*=t` строит
   normalized powers t^q. Float64 tridiagonal eigensolve не объявлен exact;
   его shared node/weight errors не оцениваются reverse-order witness.
2. **Translations, lines45–59.** Коэффициент
   `binomial(q,j)*((c_child−c_parent)/r_parent)^(q−j)
   *(r_child/r_parent)^j` правильный. Ненормированный лишний radius/Jacobian
   не добавлен. При p≤24 binomial integers не близки overflow. Reverse
   child/order paths действительно меняют порядок accumulation, но не
   coefficients — sensitivity, не независимый forward bound.
3. **Matrix mass, lines95/101.**
   `sum width*(normΓleft_F+normΓright_F)/(4π)` majorizes
   integral||Γ̃||F/(2π) по convex norm inequality, даже без PSD. Поэтому
   `mass*q^p/(d−r)` ограничивает analytic omitted field **во Frobenius**,
   а следовательно и operator norm. Это корректная, более консервативная
   замена trace mass. Zero-mass shortcut корректен в exact arithmetic;
   floating underflow/outward rounding общей программой не certified.
4. **Phase/Horner, lines162–174.** `ratio=r/(z−c)`, division на z−c и
   p terms от q0 до q(p−1) правильны. Abs используется лишь в geometry/
   bound, не подменяет complex phase. Horner и explicit powers согласованы.
5. **Near PV, lines143–160.** Общий constant term
   `Γ(E)[log(E−a)−log(b−E)]` имеет правильный sign. Residual panel integral
   `−B*width+[Γleft+B(E−left)−Γ(E)]*log|ratio|` правильный. На touching
   panel coefficient исчезает аналитически, его пропуск реализует finite
   limit вместо Inf−Inf. На outer ghost Γ(E)=0 constant log не вычисляется.
   Current differences-of-logs устраняют прежний overflow ratio при
   nextfloat knot. Остаётся заявленный near floating error, не exact log.
6. **Direct panels, lines124–141.** При L=log(z−l)−log(z−u),
   B=(z−l)L/width−1, coefficients `(L−B)Γl+BΓu` правильны. Upper-half-plane
   principal log difference даёт правильную Cauchy branch. Remote scalar
   BigFloat192 защищает cancellation coefficients, но не превращает
   ComplexF64 matrix accumulation в rigorous arithmetic bound.
7. **Real/complex contracts, lines175–240.** Real path вычисляет PV и
   ровно один `−iΓ(E)/2`; complex path требует Imz>0 и не добавляет его
   второй раз. Literal nonuniform diagonal не принудительно занулён:
   asymmetric hat даёт Γ log(hL/hR)/(2π), покрыто tiny oracle.
8. **Partition/call controls, lines179–231.** Near interval включает
   containing panel и touching neighbours; wholly-near tree nodes skipped,
   partially overlapping nodes split/direct-leaf path пропускает только
   near panels. Count-balanced disjoint children гарантируют partition
   structurally; covered-panel count подтверждает его. Far acceptance
   требует no-overlap, d>r, q≤.5 и p≤pmax. При нехватке p открываются children/
   direct leaves, tolerance не повышается. Allocation proportional node
   mass даёт суммарный exact-arithmetic truncation≤tau_far.

## Layout, input preservation, source support

`Gamma(E,a,b)` reshaped в `(E,nb²)` и result обратно в `(nb,nb)` с одним
column-major a/b ordering. При Nb5 обрабатываются **25 entries**, включая
все complex off-diagonal, а не20. Tiny five-basis test проверяет этот путь.
Raw data не Hermitize/clip; Hermitian wrapper служит только witness.
При nonzero declared hermitian/psd tolerance exact PSD не следует из gate,
и protocol корректно не называет его formal certificate.

Loader не выполняет source reads на inclusion. Он связывает slice с
полным source-file SHA и фиксированным extraction order five **local**
retarded families, исключая embedding/±. Formula
`Gamma[e,a,b]=i*(R[e,a,b]−conj(R[e,b,a]))` имеет правильные transpose/sign;
elementwise imaginary part вместо matrix adjoint не используется.
Single-k SHA не выдан за prior aggregate all-k SHA. Full-array source
сохранён; only selected-k hyperslabs читаются при будущем явном вызове.

**Qualification actual FFT comparison:** loader задаёт literal stored E
и ghosts из первых/последних adjacent gaps. Это correct nonuniform
full-hat model и совпадает с dense-reference source convention. Accepted
uniform PI FFT использует canonical range/lattice spacing, сохраняя
reported Float64 geometry deviation. Поэтому их tiny coordinate/ghost
rounding difference не следует называть bitwise equality continuous
interpolants. В actual probe сохранить обе geometry/ghost sets и разность;
независимый direct panel oracle относится именно к tree literal model.
Diagnostic allowance может включать такой effect лишь с явным описанием;
нельзя автоматически объявлять discrepancy formula bug или менять Γ/grid.

## Смысл ledger и предел вывода

Returned ledger честно разделяет far truncation estimate, forward/reverse
moment sensitivity, Horner/explicit sensitivity, near/final sum-order
sensitivity и uncertified geometry/log/outward-rounding fields.
`formal_whole_certificate=false` соответствует source. Разность порядков
может пропустить common Gauss/log/translation errors; её нельзя сложить
с far estimate и подставить как rigorous δΣ в GF resolvent theorem.
Actual frozen protocol прямо задаёт far+numerical allowance как diagnostic
gate, запрещает его молчаливое ослабление. Это appropriate exploratory scope.

Positive-hat exact integral causal; target-dependent truncated evaluator
не обязан быть globally analytic/Herglotz. Для exact Hermitian input real
Λ имеет Hermitian algebra и exact Γ half-term, но finite accumulation
errors/declared input allowances остаются. Нулевой sampled causality metric
не примет KK/spectral mass. Никаких G, density или Matsubara-integral
conclusions этот transform-only module не производит.

Tiny reference использует independent BigFloat256 nonuniform-hat divided
difference F(t)=t log t. Он не зовёт panel/tree/moment helpers; relevant
phase/endpoint/asymmetric/zero/layout checks — содержательные. Root-moment
comparison использует отдельный no-translation integral; источник Gauss
rule там общий, что ограничивает independence именно этого subtest.
Direct actual probes и independent hat reference по frozen protocol
остаются нужными; tiny pass не принимает30721-node source accuracy/RSS/time.

**Итог static gate:** current formulas/shape/domain/ledger подходят для
первого bounded actual **transform** comparison после решения ведущего.
Не найдено основания задержать его из-за math defect. Source Γ witnesses,
actual HDF5 consumer, full-source performance и all error budgets пока
not_measured этим reviewer; warm-primary priority не меняется.
