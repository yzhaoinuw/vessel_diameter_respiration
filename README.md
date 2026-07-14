# Vessel Diameter and Heart-Rate Correlation Analysis

[![Agent Collab Treaty](https://raw.githubusercontent.com/yzhaoinuw/agent_collab_treaty/main/assets/treaty-adopted.svg)](https://github.com/yzhaoinuw/agent_collab_treaty)

This MATLAB repository evaluates whether ECG-derived heart rate covaries with vessel-diameter dynamics in dexmedetomidine recordings. The authoritative methods, results, caveats, and curated figures are in [`correlation_analysis_writeup.md`](correlation_analysis_writeup.md).

## Retained Analyses

The active code supports four reported components:

1. A primary Spearman analysis pairing quality-filtered one-minute ECG HR with overlapping ten-minute vessel windows advanced in one-minute steps.
2. An all-recording sensitivity analysis that removes the HR-agreement exclusion while retaining the filtered-ECG coverage requirement.
3. Quality-control comparison of filtered ECG HR, legacy Yue ECG HR, and vessel-derived cardiac frequency.
4. An older non-overlapping Yue-HR Spearman analysis retained as a supplementary comparison.

Pearson variants, respiration coupling, five-minute sliding analyses, cross-correlation experiments, and F153-only pilots are preserved under [`archive/`](archive/) but are not part of the reported workflow.

## Active Pipeline

```text
signals.mat + im.mat
        |
        +-- extract_hr_amp_data_all_subjects.m
        |       -> hr_amp_data.mat (filtered one-minute ECG HR)
        |
        +-- extract_vessel_features_10min_overlap.m
        |       -> vessel_features_10min_overlap.mat
        |
        +-- spearman_correlation_overlap/
        |       -> primary and all-recording Spearman summaries
        |
        +-- PlotResults_continuous_v5.m / results.mat
                -> Yue HR and vessel-derived cardiac-frequency comparators
                -> hr_hr_correlation/, spearman_correlation_yue/,
                   and plot_three_hr_comparison.m
```

## Running The Analysis

Activate the project environment and run MATLAB from the repository root:

```powershell
conda activate ecg
matlab -batch "extract_hr_amp_data_all_subjects"
matlab -batch "extract_vessel_features_10min_overlap"
```

The overlap scripts live under `spearman_correlation_overlap/`. They support `overlap_analysis_mode = "curated"` and `"all_subjects"`; see [`AGENTS.md`](AGENTS.md) for the exact sequence and repository conventions.

## Data And Outputs

Raw subject data and generated outputs are intentionally ignored. A subject folder used by the active extraction requires:

- `signals.mat` with ECG timing, detected R peaks, R-peak confidence, and accepted peaks.
- `im.mat` with vessel diameter, frame timing, ECG offset, and micron-per-pixel scaling.
- `results.mat` only for the legacy Yue HR and vessel-derived cardiac-frequency comparisons.

The curated figures embedded by the write-up are tracked under `writeup_figures/`.

## Repository Map

- [`project_overview.md`](project_overview.md) — detailed active-path and data-contract map.
- [`next_steps.md`](next_steps.md) — current unfinished work.
- [`work_log.md`](work_log.md) — recent implementation and verification history.
- [`archive/README.md`](archive/README.md) — historical scripts and superseded documents.
