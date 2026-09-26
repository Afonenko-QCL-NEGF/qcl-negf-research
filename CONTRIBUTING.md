# Contributing scientific inputs

Create a study under `config/studies/` with explicit configuration sources,
units, scientific purpose and controlled variants. Include it in the appropriate
top-level plan. Keep every configuration fragment reachable and avoid copying
common model parameters into multiple studies.

Use the schema validator and tests described in the [README](README.md#validate-inputs).
Run `qcl-negf plan` against each changed entry point with the matching
`QCLNEGFRunner.jl` release. Expensive
transport calculations belong on allocated Slurm nodes; structural checks do not
replace scientific convergence and accuracy studies.

Changes to physical assumptions, layer dimensions, scattering channels or
acceptance interpretations need a scientific justification. Cite publications by
title and DOI. Preserve data units, digitization uncertainty and conditions, and
do not turn an unverified reference candidate into an accepted result by changing
its label.

Submit new results to your research data/provenance system rather than committing
machine-specific paths, large numerical artifacts, credentials or execution logs
to this repository. Integrated CI runs from the
[qcl-negf superproject](https://github.com/AfonenkoA/qcl-negf) on trusted pushes or
manual dispatch; run untrusted contributions only on disposable isolated runners.
