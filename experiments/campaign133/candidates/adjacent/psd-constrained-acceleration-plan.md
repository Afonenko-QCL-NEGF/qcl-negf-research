# PSD-constrained Γ acceleration: один следующий structural variant

**Выбор — A: projection algorithmic AA trial на PSD cone**, с исходной
unprojected Φ и native acceptance. Factor coordinates B пока не выбирать.
Это следующий возможный опыт **после** результата текущего32, не причина
его остановить. Сообщённые Root cone rejection при L1≈2389 и r≈7.18e−6
показывают новый barrier после снятия coefficient lockout, но ещё не
доказывают ни physical impossibility, ни пользу projection. Здесь только
proof/design; implementation/Julia/Pkg/Core/runtime не выполнялись.

## 1. Arithmetic foundation перед новым accelerator

Оба варианта требуют positive-source spectral A вместо извлечения tiny
weight из complex GR: при consistent Hermitian coherent parts
`A=Gc(Γ+Γemb,source)Gc†`, где native test source Σ<=iΓ, Σ>=0 даёт
Γemb,source=−iΣemb<. Это technical spectral probe, **не occupation**.
Gram representation/trace sumabs² сохраняет positive number weights;
из этого же A строятся current μ, GFD и physical Φ. Нельзя исправить
только q, оставив другое A в scattering. Proof и same-G decomposition —
`../stable-spectral-plan.md`.

Herm(h)/Herm(Λ) и consistent side/central policy должны быть явно
зафиксированы; raw coherent corrections и ΓR/ΓK witnesses сохраняются.
True physical kernel map PSD требует соответствующей CP/positive-source
опоры; numerical contractions всё ещё могут дать tiny negative Φ.
**Φ не project/clip.** Если её negative part не объясняется controlled
arithmetic, это отдельный source/operator obstruction. Из одного stable A
не следует strict floating PSD каждого channel.

## 2. Projection trial сохраняет исходные roots, но не гарантирует convergence

Для Hermitian trial X=U diag(λ)U† определить

```
P+(X)=U diag(max(λ,0)) U†,
δproj=X−P+(X),  ||δproj||F²=sum_(λ<0)λ².
```

Это orthogonal projection closed convex PSD cone в real Hermitian
Frobenius space: по unitary invariance nearest PSD diagonal — max(λ,0).
Она nonexpansive и для любого PSD root Γ*:
`||P+(X)−Γ*||F≤||X−Γ*||F`. То же blockwise для positive weighted norm
(для zero weights — seminorm statement). **Positive floor не вводится.**
Численный eigen/reconstruction backward error измеряется отдельно;
малое computed negative eigenvalue output не даёт formal cone certificate.

Если Γ*=Φ(Γ*)≽0, нулевой raw residual делает AA/linear trial равным Γ*,
и projection его не меняет: исходный feasible fixed point сохранён.
Но projection может создать неподвижность **алгоритма** без root:
scalar Φ(x)=x−1, x≥0 имеет P+(Φ(0))=0, тогда как raw Φ(0)−0=−1.
Поэтому ни ||Γnext−Γ||, ни projected residual, ни повторение projected
boundary state не являются acceptance. Проверять неизменные raw
`||Φraw−Γ||/||Φraw||`, absolute/weighted scales и full native gates.

Более сильное ограничение: для Hermitian raw Φ и любого Γ≽0
`||Φ−Γ||F≥||Φ−P+(Φ)||F`. Negative-part norm — measurable lower bound
unprojected residual; projection trial не способна скрыть его в этом metric.
Анти-Hermitian raw defect также входит в raw residual отдельно.
Это lower bound в данном candidate, не доказательство отсутствия roots
в других областях Γ-domain.

## 3. Почему B сейчас менее удачный различающий шаг

При Φ(Γ)≽0 map `T(L)=principal_sqrt(Φ(L L†))` имеет equivalent roots:
fixed L principal PSD ⇒ Γ=L²=Φ(Γ); любой PSD Γ root даёт L=√Γ.
True negative Φ делает T неопределённой; исправлять её eigenclip нельзя.
Таким образом B не решает current output-arithmetic obstruction сама.

На PD root derivative sqrt удовлетворяет `L δL+δL L=δΓ`, то есть
`δLij=δΓij/(sqrt(γi)+sqrt(γj))`. При γ→0 координаты плохо обусловлены/
не smooth. Interior change of variables даёт **similar** Jacobians
`J_L=Dsqrt J_Γ (Dsqrt)^−1`, поэтому intrinsic spectral radius не меняется.
Гарантированного излечения instability отсюда нет. Principal Hermitian
factors + real AA убирают arbitrary unitary gauge; complex unconstrained
L содержит redundant L→LU directions. Даже Hermitian AA trial L может
сменить sign branches, хотя LL† PSD. При любом B нужна Γraw residual,
а не L-step. Для нынешнего барьера A меняет меньше и сохраняет Γ history.

## 4. Один bounded следующий protocol, если Root выберет после32

Новая identity, например `Gamma-positive-source-AAtrialPSD-v1`;
старые no-clipping baselines сохраняются. Их local protocol не является
запретом пользователя на явно объявленный новый numerical method.
Старт — original archived **positive FD128-derived regridded Γ seed**,
не исправление отрицательных entries текущего Γ checkpoint.

1. Один fresh positive-source Π preflight на этом seed; текущие μ/root,
   per-channel ΓR/ΓK/KMS/raw Φ и arithmetic witnesses. До attribution
   accelerator эффекта этот operator должен быть квалифицирован отдельно:
   сравнение с прежним i(G−G†) run включает arithmetic change.
2. Same real AA-II/history3/affine4096/QR guards строит Xtrial. Сохранить
   raw Xtrial, min eigenvalue и negative-part activity; **project только
   этот trial**, затем одна fresh Π(P+Xtrial). History хранит actual
   projected Γ и fresh unprojected residual, не rejected raw Xtrial.
3. Inverse/capacity/χ/domain failure — existing one convex α=.2 fallback
   и history restart. Convex fallback допустим PSD только при PSD Φ;
   если numerical Φ делает его inadmissible, stop/arithmetic diagnosis,
   без hidden projection Φ, η или relaxation physics gates.
4. Log ||δproj|| absolute/weighted, changed blocks/rank/locator, trace
   change, raw residual **после** projection, μ/number conditioning,
   counted rejected calls. Large coefficients/projection могут ухудшить
   nonlinear response; projection-distance inequality этого не исключает.

Prospective cap: **16 total Π**, один процесс, actual MemoryMax
**13019119616B (12.125GiB)**/Swap0/Julia4/BLAS1/CPU400%, pilot≤600s,
overall≤1800s включая build/JIT и предусмотренные stages, small scalar evidence
и bounded raw checkpoints. Это proposal, не запущенный опыт/обещание цены.
При raw trigger сохраняется обязательный **unprojected** all-family native
materialization и strict3-pass audit в отдельно заранее предусмотренном
budget; если его нет/он fails — Γ-only unresolved. Никакой persistent η,
post-fit density scaling или ослабления tolerances. Π output negativity,
domain failure, budget либо неизменный raw residual дают честный stop.

## Первичные работы и предел переноса

[Higham, LAA103,103–118(1988), author PDF](https://www.maths.manchester.ac.uk/~higham/narep/narep126.pdf)
даёт nearest-PSD construction; complex-Hermitian extension здесь следует
из приведённого unitary proof. [Higham–Strabić, Numerical Algorithms72,
1021–1042(2016)](https://link.springer.com/article/10.1007/s11075-015-0078-3)
демонстрирует AA в alternating PSD/correlation projections, **не** QCL SCBA.
[Zhang–O’Donoghue–Boyd, SIOPT30,3170–3197(2020), primary PDF](https://web.stanford.edu/~boyd/papers/pdf/scs_2.0_v_global.pdf)
и [Garstka–Cannon–Goulart, ECC2022,435–440, author manuscript](https://ora.ox.ac.uk/objects/uuid%3A46531826-6e0a-4cd4-92d7-1f9d49eb8e95/files/sht24wk58m)
исследуют safeguarded AA при nonexpansive hypotheses. Эти hypotheses
для Γ-map не доказаны; papers не дают этому projected-AA variant theorem
о QCL contraction/convergence. Источники проверены primary browsing.
