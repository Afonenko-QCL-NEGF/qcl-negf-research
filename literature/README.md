# Literature data

`reference-2019` identifies the publication **Thermoelectrically cooled THz
quantum cascade laser operating up to 210 K**, Applied Physics Letters 115,
010601 (2019), [DOI 10.1063/1.5110305](https://doi.org/10.1063/1.5110305).
The label names a reference dataset, not an original experiment performed by
this software project.

| File | Content |
| --- | --- |
| `reference-2019/source.yaml` | Publication identity, digitization method, source artifact hashes, comparison limits and local data hashes |
| `reference-2019/structure.csv` | Published dimensions, doping and device characteristics with units and source descriptions |
| `reference-2019/negf_iv_digitized.csv` | Six digitized samples of the publication's calculated NEGF current curve, with digitization uncertainties |

The CSV current curve is a digitization of a **calculated NEGF curve**, not raw
experimental observations. Its active-region temperature is 200 K. The horizontal
calibration uses the publication's 56 mV-per-period peak statement. Quoted
uncertainties describe digitization only; they are not estimates of model error.
Do not extrapolate this curve to score points outside 40–60 mV per period.

Compare current density against voltage **per period**. Terminal voltage/current,
threshold current, optical output power, maximum operating temperature and
absolute gain require additional device physics and cannot be inferred from
agreement with these data. The 405-period geometry value in the example model is
not directly established by the publication and must not become a terminal
voltage validation target.

The source PDFs are not bundled. Their recorded hashes identify the source
artifacts used for digitization; local CSV hashes allow accidental data changes
to be detected without downloading those PDFs. When replacing or extending a
digitization, document the source curve, units, method and uncertainty, then
update its hash deliberately.
