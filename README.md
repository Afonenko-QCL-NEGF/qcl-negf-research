# QCL-NEGF Research

Versioned scientific inputs and literature provenance for quantum cascade laser
NEGF calculations. The numerical implementation lives in
[QCLNEGF.jl](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl), with scientific plans and CLI in
[QCLNEGFRunner.jl](https://github.com/Afonenko-QCL-NEGF/QCLNEGFRunner.jl); execution and
provenance on Slurm are provided by
[qcl-negf-aiida](https://github.com/Afonenko-QCL-NEGF/qcl-negf-aiida).

This repository contains input definitions, not precomputed scientific results.
Every included study is reachable from one of the three entry points below.

| Entry point | Scientific purpose | Execution |
| --- | --- | --- |
| `config/operators.yaml` | Seven algebraic, metric, representation, continuum, cavity, LO equilibrium and spectral resolution checks | Diagnostic calculations |
| `config/spectra.yaml` | Structure spectra and basis comparisons for a literature-derived structure and a synthetic well | Diagnostic calculations |
| `config/convergence.yaml` | Independent 48 mV candidate and nine controlled numerical comparisons at 70 K | Expensive stationary calculations |

## Run a study

Install `qcl-negf` from a compatible `QCLNEGFRunner.jl` release. From this repository:

```sh
mkdir -p plans
qcl-negf plan config/operators.yaml > plans/operators.json
qcl-negf run-plan plans/operators.json results/operators
```

The frozen JSON plan contains resolved inputs and a fingerprint. Keep it with the
results. To submit the same plan to the department cluster, use the
[AiiDA submission interface](https://github.com/Afonenko-QCL-NEGF/qcl-negf-aiida)
with a registered `qcl-negf` Code and explicit Slurm resources. Configuration files
contain scientific choices; they do not allocate cluster resources.

Start with the diagnostic entry points. The convergence bundle uses 30,721 energy
nodes for its reference candidate and up to 46,081 for a comparison. It is not a
laptop smoke test. A technically completed process does not certify convergence,
operator checks, numerical accuracy or agreement with experiment; inspect the
solver's recorded scientific status and diagnostics.

## Validate inputs

Use the Python 3.14 environment prepared by the
[qcl-negf superproject](https://github.com/Afonenko-QCL-NEGF/qcl-negf). From its root:

```sh
deno task prepare
uv run --locked python components/qcl-negf-research/tools/validate.py
uv run --locked pytest components/qcl-negf-research/tests
```

Validation checks the study graph, every merged configuration variant against
the shared JSON schemas, duplicate YAML keys, missing or escaping source paths,
unreachable configuration files and literature data checksums. This is structural
validation; `qcl-negf plan` additionally applies the solver's semantic checks.
Neither procedure runs the physical convergence experiment.

The superproject Git tree fixes the contracts revision and its generated
`uv.lock` fixes external dependencies. Integrated catalog checks run in this same
environment on its self-hosted GitHub Actions runner.

## Scientific interpretation and contributions

- [Model assumptions and controlled comparisons](docs/scientific-inputs.md)
- [Literature data, provenance and comparison limits](literature/README.md)
- [Contribution checks](CONTRIBUTING.md)

Code and project documentation use the [MIT license](LICENSE). Publication
metadata and digitized numerical data retain their source attribution in
`literature/reference-2019/source.yaml`; the source article is not redistributed.
