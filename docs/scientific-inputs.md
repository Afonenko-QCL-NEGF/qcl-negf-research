# Scientific inputs

## Models

`config/model/reference-2019-70k.yaml` describes a two-well GaAs/AlGaAs period
with 25% Al barriers. Its layer sequence is 3.26 / 7.99 / 1.90 / 8.40 / 2.90 /
5.16 nm, with the 2.90 nm segment doped and a sheet density of
4.5 × 10¹⁴ m⁻². Geometry and sheet doping trace to DOI
[10.1063/1.5110305](https://doi.org/10.1063/1.5110305).

The present model temperature is 70 K. The included digitized publication curve
is at 200 K. Those are distinct conditions: the curve is provenance and a resource
for a separately designed comparison, not a pass/fail threshold for these studies.

The synthetic one-well geometry combines
`config/model/common-comparison-70k.yaml` with
`config/model/synthetic-one-well-70k.yaml`. It uses 15% Al barriers, a 120 meV
conduction-band offset and parabolic masses. It is a numerical control without a
claim of reproducing an experimental device.

## Configuration composition

Each study lists `configuration.sources` explicitly. Sources are merged in order;
mapping values merge recursively and sequences replace earlier sequences. Study
overrides apply next, then the selected variant's overrides. The solver resolves
temperatures, voltages and output policy into the frozen execution plan.

`config/policies/analytical-diagnostic.yaml` selects the diagnostic discretization
and algorithm policy. Its name does not imply that the stationary transport
problem has an analytical solution. The reference and grid studies explicitly
raise the iteration budgets for stationary calculations.

The operator studies isolate different operator identities and numerical
approximations. The spectrum studies vary basis size, basis period window and
spatial resolution; a failed variant is diagnostic evidence, not a reason to
silently change the model or weaken a tolerance.

## Reference candidate and grid comparisons

`reference-48` computes an independent candidate at 48 mV per period and 70 K.
`grid-refinement-48` runs nine declared variants covering energy spacing and
window, energy-grid phase, transverse momentum density/cutoff, spatial resolution,
basis size and embedding. Their `controlled_paths` identify what changes.

These studies are independent inputs. Their inclusion in one bundle does not
declare the reference accepted and does not enforce a scientific acceptance gate
in the scheduler. Assess comparisons using the recorded physical observables,
convergence status and diagnostics from both calculations. In particular,
agreement of two unconverged calculations does not establish discretization
accuracy.

Record the research repository commit, solver release, frozen plan fingerprint,
execution environment and result artifact identities with any scientific claim.
AiiDA records execution provenance; the result contracts retain the solver's
scientific classifications. CPU allocation and technical exit status are separate
from those classifications.

## Supported scope

The catalog uses ordinary scientific plans and independently scheduled execution
nodes. It has no evidence-driven resource broker, reserve task admission or
automatic experimental acceptance policy. Cluster concurrency, memory and wall
time belong to the AiiDA/Slurm deployment.
