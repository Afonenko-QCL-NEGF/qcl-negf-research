# Physical LO width V2: integrability — отдельное необходимое условие

Future research-only delta к V1 SHA
f13c0d4b3b3f6212dc1baa2d49d0c955a3c55dd8c65b5036721d242a52a7f403;
V1 не изменён. Основание — независимая review
`reviews/noEE-frozen-response-and-linewidth-review.md`. Реализации/Julia/Pkg/
H5/новых maps нет, actual model change не одобрен.

Odd causal phonon spectral function и positive-half b(ε) **недостаточны** для
finite Born bath. В дополнение к общей causal R/L/G specification и
energy-dependent Bose KMS требуется normalization/zero-width limit и,
например, sufficient thermal-weighted condition

`∫₀∞ b(ε)[2N_B(ε)+1] ||Kε|| dε < ∞`.

Здесь bKε — coupling spectral measure в одной energy convention, а ||Kε||
ограничивает соответствующий operator на retained matrix/k fields. Это
sufficient condition для finite thermal-weighted operator; для specific G
на границе нужны также существующие resolvent/spectral limits. При ε→0
`N_B(ε)~kBT/ε`: произвольная positive-half Lorentzian с b(0)>0 и bounded
ненулевым Kε даёт log infrared divergence. Такая interpolation не является
готовым physical bath. При ω→∞ также нужны tail/moment conditions под
заявленные quantities; одно слово causal их не задаёт.

Future variant обязан задать low-energy behavior (часто b(ε)=O(ε^α),α>0
при bounded vertex достаточно около0), high-energy tail, units и oscillator
strength/sum-rule normalization, а не подбирать width под mesh/convergence.
Emission и absorption используют одинаковые b,Kε и N_B(ε); greater functional
есть `−i∫b Kε[(Nε+1)p(E−ε)+Nεp(E+ε)]dε`. Это сохраняет bath KMS termwise.
ΣR дисперсия и identity `ΣR−ΣA=Σ>−Σ<` должны соответствовать тем же sources;
electronη/independent smoothing этого не заменяют.

V1 structural conclusion сохраняется: разрешённые off-sector transfers могут
разрушить ideal energy-phase graph, но mixing/uniqueness/convergence не следуют
из finite linewidth автоматически. T1/T2/HWHM calibration текущего70K device
по-прежнему отсутствует. e–e не предлагается. Перед physical variant выбран
no-model-change control с **fresh-base replay** из соседней V2-delta.
