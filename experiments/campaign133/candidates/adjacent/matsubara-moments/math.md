# M3 density tail after completed Mats ALL5, write-only

Actual source: fixedμ/C-static/full-line thermal LO with reference-tail closure;
N128/256,32/32 inner calls completed. This module performs **zero new inner Φ,
Dyson solves or μ fits**. Existing main/module/actual artifacts remain immutable.

Let z=iω, Hrel=H45−μI, Σloc(z)=S0/z+S1tilde/z²+O(z⁻³), repeated in all9
cells. Resolvent series gives central moments I,hrel,(Hrel²)cc+S0 and
(Hrel³)cc+hrelS0+S0hrel+S1tilde. Exact bare45 subtraction removes bare moments
at every order. Herm[S0/(iω)³]=0 for Hermitian S0; (iω)⁻⁴=+ω⁻⁴. Thus the
leading omitted density term is **+t4 B**, t4=(2/β)Σ_(n≥N)ω_n⁻⁴,
B=hrelS0+S0hrel+S1tilde. No factor9: local self-energy insertion is diagonal;
only hrel from the same central block enters this interaction moment term.

Static covariance gives S0_el=ΣKs[I], S1_el=ΣKs[hrel]. Independently expand the
known thermal LO pole functional: weights NB+1−f(x), NB+f(x), poles x±Ω.
Their sum is2NB+1, their first weighted moment is(2NB+1)x+Ω(1−2f(x)). Therefore
S0_LO=(2NB+1)KLO[I], S1_LO=KLO[(2NB+1)hrel+Ω(I−2ρ)]. This uses normalized
central spectral mass I and first moment hrel, **not** an imposed density.
Consequently B=B0−2ΩKLO[ρ]; solve

    (I+2Ω t4 KLO) ρM3 = ρbareN + t4 B0.

Sign is plus on the left. KLO is the same exact native qK/wk operator, applied
only to small moment matrices. A dense275×275 general complex LU system suffices (CP is not SPD);
its residual/condition are reported. Conditional norm η=2Ωt4||KLO[I]||∞<1
(for the declared CP map) implies inverse bound1/(1−η). The reported residual
bound concerns this **linear moment equation**, not floating total error,
nonlinear SCBA or an asymptotic finite-N remainder certificate.

`tails` reuses V4 odd-series identity at BigFloat256 precision:
Σodd⁻⁴=π⁴/96−prefix. Independent positive integral lower/upper bounds check
the resulting scalar t4. BigFloat→Float64 is numerical arithmetic, not interval
certification. `scalar_oracle` uses G=[z−ε−a/(z−b)]⁻¹ with independently known
two poles/residues: B=a(2ε+b)=Σweights λ³−ε³. Exact Fermi residues independently
check M3 sign and128→256 improvement without a self-consistent kernel.

Stored ρbareN already includes exact bare45 reference density (all bare orders).
It is **not** merely classic M1. Report classic M1 separately from saved G;
apply interaction-only M3 correction to the stored bare-subtracted density.
Report BOTH one-shot `ρN+t4 B(ρN)` and self-consistent moment closure above,
their difference, Pauli eigenvalues, Hermiticity and physical numbers;
no rescale, projection or forced DOPE fit.

The finite-N LO map substitutes Gref outside ±N. Its extended high-frequency Σ1 would involve the saved ρN. One-shot uses
that coefficient; the self-consistent moment completion is a NEW asymptotic
tail-closure identity, changing feedback by O(t4²), not an exact identity of
the saved finite-N map or the infinite self-consistent Mats problem.
Measured low32 Σ difference8.05e−11/G1.52e−13 is useful existing evidence, not
an exact bound on all omitted frequencies. Agreement of corrected128/256 is
a diagnostic; remaining cutoff/nonlinear/reference-tail error remains explicit.
N512/1024 would be a separately authorized new experiment if needed.
