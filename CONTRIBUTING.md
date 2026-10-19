# Contributing scientific inputs

Create a study under `config/studies/` with explicit configuration sources,
units, scientific purpose and controlled variants. Include it in the appropriate
top-level plan. Keep every configuration fragment reachable and avoid copying
common model parameters into multiple studies.

For documentation-only text or link corrections, review the diff and check changed
links, anchors and formatting. No Julia, package installation, full suite, plan
execution or transport calculation is required. Record checks that were not run;
changes to executable examples, study inputs, schemas or scientific data need the
relevant checks below.

For changed study inputs or validation code, use the relevant schema validator and
tests described in the [README](README.md#validate-inputs) in a prepared matching
environment. Run `qcl-negf plan` against each changed entry point with the matching
`QCLNEGFRunner.jl` release as an assigned semantic check. It requires Julia; it
does not run the transport solve. The README's environment preparation and full
suite are separate assigned stages with a finite shared budget. Missing runtime
limits the completed checks and does not block an independent text correction.
Expensive transport calculations belong on allocated Slurm nodes with a declared
scientific stage, frozen model/plan/criteria, fresh admission and finite budget;
structural checks do not replace scientific convergence and accuracy studies.

Changes to physical assumptions, layer dimensions, scattering channels or
acceptance interpretations need a scientific justification. Cite publications by
title and DOI. Preserve data units, digitization uncertainty and conditions, and
do not turn an unverified reference candidate into an accepted result by changing
its label.

Submit new results to your research data/provenance system rather than committing
machine-specific paths, large numerical artifacts, credentials or execution logs
to this repository. Integrated CI runs from the
[qcl-negf superproject](https://github.com/Afonenko-QCL-NEGF/qcl-negf) on trusted pushes or
manual dispatch; run untrusted contributions only on disposable isolated runners.
