# ras-ss

A **summary-statistic extension of the Regional Association Score (RAS)**
method, so that regional association testing can run on GWAS summary statistics
instead of individual-level data.

RAS (Jiang & Zhang, 2025) finds trait-associated *regions* by converting
regional association scores into a time series and applying change-point
detection, outperforming single-SNP and set-based tests on power in challenging
(i.e. sparse-causal) cases while maintaining low false-positive rate. Its one
limitation is that it needs individual-level genotypes and phenotypes, which
most published GWAS don't share. `ras-ss` replaces RAS's per-window regression
with a 1-df weighted burden statistic (`T_burden = w'Z / sqrt(w'Rw)`) derived
from marginal Z-scores and an external LD reference, then inputs the resulting
profile into RAS's **own, unmodified** change-point detection.

**Reproducible write-up:** <https://anson-li8.github.io/ras-ss/>

## Status

Simulation-validated; a real summary-statistic run is pending. On a 300-SNP
toy chromosome (500 null / 200 signal replicates, exact binomial CIs):

| Run | Type I error | Power |
| --- | --- | --- |
| ras-ss | 0.038 (0.023–0.059) | 0.935 (0.891–0.965) |
| individual-level, same data (paired) | 0.036 (0.021–0.056) | 0.940 (0.898–0.969) |
| individual-level, original `ras()` (`num_rep = 5`) | 0.010 (0.003–0.023) | 1.000 (0.982–1.000) |

On identical simulated data, ras-ss is consistent with the individual-level method on
**199 of 200** signal replicates and **499 of 500** null replicates, showing the
measured cost of the summary-statistic substitution. Type I error is within the
original paper's stated range. The ~6-point power gap to the 5-averaged
original demonstrates the gap without resampling-average, which summary statistics cannot
reproduce (the per-individual split assignments are not in marginal Z-scores).
It is the cost of the no individual-level data constraint the extension was built for, not
a defect in the statistic. ras-ss is meant for the summary-stat-only studies the original cannot run on,
not as a replacement when individual-level data exists.

## Read it in order

1. [Problem Statement](https://anson-li8.github.io/ras-ss/00_problem_statement.html) — original RAS, the data-access problem, and the summary-stat `T_burden` derivation.
2. [Simulation](https://anson-li8.github.io/ras-ss/01_simulation.html) — simulation with a window-size sweep and the `slope_check_window_size` detection fix.
3. [Validation](https://anson-li8.github.io/ras-ss/02_simulation.html) — three-run comparison and validation diagnostics to confirm consistency with original.

## Repository layout

This is a [workflowr] project

- `analysis/` — the R Markdown pages (problem statement, two simulations,
  about, index, license).
- `code/ras_ss.R` — the shared scan / detection functions, sourced by the
  simulation pages.
- `docs/` — the built website (served by GitHub Pages).
- `analysis/_site.yml` — the site navigation.

## Citation

Jiang & Zhang (2025), *Empowering genome-wide association studies via a
visualizable test based on the regional association score*, PNAS 122(9)
e2419721122 for the original RAS method and the
[RAS package](https://github.com/hepingzhangyale/RAS).

## License

Code: MIT (see `LICENSE`). Written content: CC BY 4.0.

[workflowr]: https://github.com/workflowr/workflowr
