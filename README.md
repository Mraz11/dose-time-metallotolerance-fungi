# Dose–time metallotolerance screening in filamentous fungi (As, Pb, Hg) — reproducible R pipeline

![R](https://img.shields.io/badge/R-%E2%89%A54.2-blue)
![Figures](https://img.shields.io/badge/Figures-1%E2%80%938-success)
![Outputs](https://img.shields.io/badge/Exports-PDF%20%2B%20600dpi%20TIFF-informational)

This repository contains the **full reproducible R workflow** used to generate **publication-ready Figures 1–8** and associated **CSV tables** for the manuscript on dose–time metallotolerance profiling of *Aspergillus niger*, *Trichoderma asperellum*, and *Rhizopus stolonifer* under **As(III), Pb(II), and Hg(II)** stress (0–1000 mg L⁻¹) across **24–72 h**.

The pipeline converts a simple agar radial-growth assay into quantitative endpoints:
- colony diameter (mean ± SD)
- tolerance index (TI) and growth inhibition (%)
- plate-level growth rate (slope, mm h⁻¹)
- endpoint contrasts with Tukey lettering (48 h)
- TI heatmap (48 h)
- LL.4 dose–response curves for TI at 48 h
- **IC50 (TI50) with bootstrap 95% CI + censoring** and a summary plot (Fig 8)

---

## Reproduce all figures (one command)

From the repository root:

```bash
Rscript scripts/run_analysis_make_figures.R
# dose-time-metallotolerance-fungi
