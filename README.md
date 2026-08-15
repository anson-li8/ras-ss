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

## Status

Simulation-validated at two scales; a real summary-statistic run is pending.

**Paper Replication (1,000-SNP, paper design):** With 500 null reps per trait and 100 signal reps per scenario, ras-ss controls type I error at exactly 0.05 for both traits (exact 95% CIs ≈ 0.03–0.08), is close to the individual-level method on power across all scenarios, and calls ~0 false regions (FPR-A ≈ 0).
Its core profile-generation step is ~10× faster than the paired individual-level scan (0.045 vs 0.48 s/rep, a closed-form burden statistic vs per-window regressions), and total per-rep runtimes are similar (6.8 vs 7.7 s) only because both go through the same shared change-point-detection cost, and the 5-averaged original runs for ~30 s/rep. ras-ss's main value is that it runs on summary-stat-only studies the original cannot.

**Validation (300-SNP toy):** On identical simulated data, ras-ss is consistent 
with the individual-level method on **199 of 200** signal replicates and 
**499 of 500** null replicates. Type I error (0.038) is within the original 
paper's stated range. The ~6-point power gap to the 5-averaged original 
demonstrates the gap without the resampling-average, which summary statistics 
cannot reproduce. ras-ss is meant for the summary-stat-only studies the 
original cannot run on, not as a replacement when individual-level data exists.

## Read it in order

1. [Problem Statement](00_problem_statement.html) — original RAS, the data-access problem, and the summary-statistic burden statistic (T_burden) derivation.
2. [Pilot Study](01_simulation.html) — the 300-SNP tuning run: window-size sweep and the slope_check_window_size detection fix.
3. [Validation](02_simulation.html) — ras-ss vs the individual-level method on matched 300-SNP data vs the original ras() with num_rep = 5.
4. [Paper Replication](03_simulation.html) — the main study: the paper's simulation design (continuous/dichotomous × M = 1/3 × q) on a 1,000-SNP chromosome, with type I, power, FPR-A, and runtime.

## Repository layout

This is a [workflowr] project

- `analysis/` — the R Markdown pages (problem statement, pilot study, validation,
  paper replication, about, index, license).
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
