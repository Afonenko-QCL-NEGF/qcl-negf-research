# Почему full-chain bond zero не выбран как новый physical check

**Новый helper/run не реализован. Current остаётся not_measured;
дополнительной transport/conservation acceptance нет.** Предложенный
контроль net current через8 cell boundaries зануляется алгебраически
для любого cell-local Σ на каждом Matsubara node, даже до SCBA/μ/root
или cutoff convergence. Он мог бы проверить LU/index/transpose wiring,
но не даёт независимой физической опоры сверх существующих controls.
Root исключил новый runtime после этого вывода. Reviewer Julia/HDF5
не запускал; исходные Mats files не менял.

## Тождество для actual finite45 construction

Mats module1434e86f… `context:12–31` строит9-cell open chain:
H[j,j+1]=Tplus,H[j+1,j]=Tplus†, onsite hH. `greens:98–109`
использует D(z)=zI−H−diagΣ(z) и full inverse, сохраняя central block.
Σ повторяется локально по cells, межклеточной Σ нет.

Пусть P_j выбирает все orbitals первых j cells, j=1…8. Для exact
G=D⁻¹ имеем DG−GD=0. Trace каждого onsite commutator равен0;
внутренние связи prefix сокращаются, через cut проходит одна связь:

    a_j=tr(Tplus G[j+1,j]) = tr(Tplus† G[j,j+1])=b_j.

Для HermG=(G+G†)/2 отсюда

    2 Im tr(Tplus HermG[j+1,j])=Im(a_j−b_j)=0.

Не нужны TR/reality hH, PSD/causality Σ, thermal FDR, SCBA stationarity
или правильная μ. Даже ошибочная локальная Σ даст этот ноль при верном
inverse и pairing. Complex physical basis не отменяет trace identity;
внутренние orbital loop currents не равны net flow через весь cell cut.
Open chain существенна: ring имеет дополнительную связь через prefix
и не подчиняется выводу об отдельной bond из одного cut.

Bare fullρ0=fβ(H−μI) коммутирует с Hermitian H и даёт тот же ноль.
В reference-subtracted estimator

    ρ_N=ρ0+(2/β)Σ_(n<N) Herm[G_n−G0_n]

каждое слагаемое имеет zero net bond current. Любой finite N уже даёт
ноль в exact arithmetic. Поэтому N128→256 zero-current agreement не
ограничивает density/current tail, не проверяет μ fit или continuum
spectral identity. Ненулевой floating witness отражал бы solve/assembly/
indices/rounding; малый witness не доказывал бы equilibrium transport.

## Знак и нормировка

При Hupper=Tplus и ρ_ab=⟨c_b†c_a⟩ source native `bond_current`
в e171 `physics/observables.jl:266` использует2Re(T G<right,left).
Так как ∫G</(2π)=iρ, native electron-flow convention в source units:

    Jbar_j=−2 g_s Σ_k wk Im tr(Tplus ρ_k[j+1,j]).

Предложенный +2Im имеет противоположную ориентацию; ноль от знака не
зависит. Дополнительного2π в ρ-based expression нет: интеграл уже
включён вρ. Physical conversion использует native positive scale
J0=e E0/(ħ L0²), затем при необходимости A/m²→A/cm²; hidden η не нужна
и не вводилась. Знак/units здесь выведены для сопоставимости, а не
для заявления выполненного current measurement.

Вывод ограничен данным proposed observable/construction. Он не
исключает другие независимые physical checks, но не разрешает их
новый runtime и не принимает current/fixed-point/charge conservation.
