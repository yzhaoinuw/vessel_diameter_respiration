# Project Overview

This repository contains the MATLAB workflow and reporting artifacts for heart-rate and vessel-diameter correlation analyses in dexmedetomidine mouse recordings. It is an analysis repository rather than an application: scripts read local subject MAT files, generate ignored intermediate/results artifacts, and feed the curated figures and conclusions in `correlation_analysis_writeup.md`.

## Authoritative Scope

`correlation_analysis_writeup.md` defines the retained scientific scope. The current repository supports:

- quality-filtered one-minute ECG HR;
- overlapping ten-minute vessel windows advanced by one minute;
- curated and all-recording Spearman sensitivity analyses;
- filtered ECG/Yue/vessel HR quality comparisons; and
- the older non-sliding Yue-HR Spearman result retained as a supplement.

An analysis not represented in that document belongs in `archive/` unless the maintainer explicitly reopens it.

## Active Runtime Path

### 1. Legacy comparator generation

[`PlotResults_continuous_v5.m`](PlotResults_continuous_v5.m) generates `results.mat`. The retained workflow uses its `avg_HR_Yue_10min` ECG estimate and `hr` vessel-spectrum peak frequency for quality control and the older supplementary analysis.

### 2. ECG HR extraction

[`extract_hr_amp_data_all_subjects.m`](extract_hr_amp_data_all_subjects.m) scans `data/*/signals.mat` and `data/*/im.mat` and writes `hr_amp_data.mat` per recording. The primary overlap analysis consumes `hr_filter_timepts_min` and `hr_minute_bpm_filtered` from these files.

### 3. Overlapping vessel features

[`extract_vessel_features_10min_overlap.m`](extract_vessel_features_10min_overlap.m) writes `vessel_features_10min_overlap.mat` with ten-minute windows advanced every minute. Retained features are slow vasomotion amplitude (`amp_svm`), slow vasomotion power (`bp_svm`), and cardiac-band pulsation amplitude (`amp_c`).

### 4. HR quality control

[`hr_hr_correlation/`](hr_hr_correlation/) retains the filtered ECG-versus-Yue, filtered ECG-versus-vessel, and Yue-versus-vessel comparisons. [`plot_three_hr_comparison.m`](plot_three_hr_comparison.m) produces the representative three-series overlay used by the write-up.

### 5. Primary overlap Spearman analysis

[`spearman_correlation_overlap/`](spearman_correlation_overlap/) contains the correlation, Fisher-z group test, diagnostic plot, and group-summary plot scripts. The same code supports the curated cohort and the broader all-recording sensitivity analysis through `overlap_analysis_mode`.

### 6. Older non-sliding supplement

[`spearman_correlation_yue/`](spearman_correlation_yue/) reproduces the Yue-HR Spearman analysis used for Supplementary Figure 2.

## Repository Structure

```text
project_root/
|- AGENTS.md
|- README.md
|- correlation_analysis_writeup.md
|- project_overview.md
|- next_steps.md
|- work_log.md
|- work_log_archive/
|- PlotResults_continuous_v5.m
|- extract_hr_amp_data_all_subjects.m
|- extract_vessel_features_10min_overlap.m
|- plot_three_hr_comparison.m
|- hr_hr_correlation/             # retained HR-quality scripts
|- spearman_correlation_overlap/  # primary + sensitivity analysis
|- spearman_correlation_yue/      # older supplementary analysis
|- writeup_figures/               # four tracked curated figures
|- archive/                       # historical scripts and documents
|- data/                          # ignored local subject data
`- results/output folders         # ignored generated artifacts
```

## Active Versus Archived

### Active and relevant

- The four root MATLAB scripts shown above.
- The three retained script directories shown above.
- `correlation_analysis_writeup.md` and its four tracked figures.

### Archived

`archive/scripts/` contains respiration coupling, Pearson variants, earlier five-minute sliding and cross-correlation workflows, F153-only pilots, superseded HR-quality scripts, and other exploratory coupling code. `archive/docs/` contains superseded planning, pipeline, methods, and coupling documents.

Archived paths may depend on old directory layouts, variable schemas, or generated files. They are retained for provenance, not as supported entry points.

## Data Contract

Active raw extraction expects one folder per recording beneath `data/`.

Required for current HR/vessel extraction:

- `signals.mat`: `ECG`, `t_ECG`, `detected_r_peaks`, `r_peak_confidence`, and `good_r_peaks` as required by the extractor.
- `im.mat`: vessel diameter, frame timing, ECG timing offset, and `umperpix` scaling. The scripts contain fallbacks for known timing/diameter field variants.

Required only for legacy comparators and the non-sliding supplement:

- `results.mat`: `timepts`, `avg_HR_Yue_10min`, `hr`, and relevant legacy vessel features.

All raw data and regenerated output tables/figures remain ignored. Do not force-add them without a deliberate reproducibility or release decision.

## Verification

There is no automated test suite. Proportionate verification consists of:

1. `git diff --check`;
2. `treaty validate .` for collaboration docs;
3. MATLAB execution of each affected analysis stage;
4. inspection of generated CSV/MAT schemas; and
5. visual inspection of regenerated figures when plot code changes.

## Practical Reading Order

1. [`README.md`](README.md)
2. [`correlation_analysis_writeup.md`](correlation_analysis_writeup.md)
3. [`AGENTS.md`](AGENTS.md)
4. The applicable retained extraction/correlation scripts
5. [`archive/README.md`](archive/README.md) only when historical provenance is needed

## Known Boundary

The tracked write-up figures are curated copies of generated outputs. If analysis code or cohort rules change, regenerate the source outputs, copy the accepted figures into `writeup_figures/`, and update the write-up and work log together.
