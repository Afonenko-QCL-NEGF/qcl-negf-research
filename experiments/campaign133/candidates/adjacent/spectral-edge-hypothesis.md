# Spectral edge: правдоподобная, пока не установленная гипотеза

Сообщённый Root рост TrΓ/Nb≈30.8→995µeV между82.399/82.438meV
**совместим** с self-consistent spectral edge и замедлением; он также
совместим с coarse collocation/iteration error. Γtrace average не является
linewidth отдельной occupied mode или density of states. Здесь proof/design
only: solver/Julia/H5/arrays не запускались, archived plans не менялись.

Для scalar CP-SCBA с s>0 (energy²)

```
G=(z−e−sG)^−1,
G=[z−e−sqrt((z−e)²−4s)]/(2s), Im G<0,
Γ=−2Im(sG)=sqrt(4s−(E−e)²)  внутри band.
```

Γ имеет square-root onset. Для retarded fixed-point map
`J_G[v]=sG²v`, stability `B=1−sG²` стремится к0 на edge;
response растёт как inverse square root distance. **Внутри scalar band
|sG²|=1 тоже**, поэтому один modulus≈1 не локализует edge: нужны
distance eigenvalue к+1/near-null B и mode shape. У damped map eigenvalue
равен1−α+αsG². Это slowing, не доказательство отсутствия root.

Если буквально fit две Γ-точки zero-background scalar onset,
`Γ²≈C(E−Ec)`, то edge лежал бы всего≈0.037µeV ниже первой точки
(δ/h=30.8²/(995²−30.8²)≈.000959). Такая grid phase возможна и одновременно
плохо разрешена. Background channels/matrix averages делают этот fit
неидентифицируемым; скачок не является evidence настоящей discontinuity.
LO copies естественно переносят feature;82.4−2×36.7≈9meV близко source μ,
поэтому unoccupied high-energy feature нельзя автоматически отбросить.

Matrix аналог для **линейной CP self-energy** S[G] имеет stability
`B=I−(X→G S[X] G)`. [Ajanki–Erdős–Krüger, PTRF173(2019),293–373](https://doi.org/10.1007/s00440-018-0835-z)
исследует MDE stability; [Alt–Erdős–Krüger, Doc.Math.25(2020),1421–1539,
§§4–7, primary paper](https://arxiv.org/pdf/1804.07752) выводит square-root
regular edges при дополнительных positivity/flatness/boundedness hypotheses
(их upper-Im m соответствует m=−G). Для QCL эти hypotheses не предъявлены:
actual static **PI K[A]**, dynamic LO shifts, finite window/chain и current-μ
constraint не являются просто local linear S[G]. Theorem нельзя перенести
по слову Dyson. Primaries проверены browsing.

**Один различающий следующий test — frozen patch consistency/response audit.**
После positive-source qualification freeze тот же continuous positive-hat
Γ, hH/T, ghosts и source μ. На patch вокруг82.4meV и его ±Ω/−2Ω copies
вычислить actual raw Φ в old nodes **и midpoints**, включая необходимые
all-k contractions и shifted A на тех же physical energies. Self-energy
должна быть Cauchy transform **того же source Γ**, не linear interpolation
storedΣR и не новый ghost taper. Без SCBA updates сравнить unweighted/
occupied norms `Φ(E)−Γtilde(E)` между nodes и в nodes. В рамках того же
frozen audit один smooth patch tangent из предыдущего Jv plan проверяет
near-null actual restricted B, с reported leakage за patch; fixed-μ
response не выдавать за full current-μ Jacobian.

Большой midpoint defect над qualified arithmetic/nodal residual означает
coarse-collocation contribution. Малый midpoint defect вместе с verified
near-null B и устойчивым mode поддержит edge-criticality hypothesis,
но не докажет physical edge/global convergence. Если off-grid evaluator/
tangent ещё не qualified, результат `not_measured`; не запускать новый
fine SCBA ради подтверждения заранее выбранного диагноза.
