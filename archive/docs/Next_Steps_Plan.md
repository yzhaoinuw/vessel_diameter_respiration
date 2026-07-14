# Next Steps Plan

## Current State

The project now has an all-subject 5-minute-bin HR/vessel feature dataset.

Completed pieces:

- `extract_hr_amp_data_all_subjects.m`
  - creates `hr_amp_data.mat` in every subject folder.
- `compute_hr_diam_corr.m`
  - computes early, late, and overall Spearman correlations between `hr_median_bpm` and all vessel features.
- `plot_hr_diam_amp_coupling.m`
  - plots HR and vessel features over time for included subjects.
- `plot_hr_amp_corr_boxplots.m`
  - plots across-subject correlation distributions for early, late, and overall periods.

Primary variables now used:

- HR feature: `hr_median_bpm`
- Vessel features: `amp_svm`, `bp_svm`, `amp_c`, `diam_mean`
- HR QC: `dubious_hr_bin_fraction`, `is_dubious_hr_bin`

Current exclusion rule:

- Exclude subjects with `dubious_hr_bin_fraction > 0.50`.

Currently excluded by that rule:

- `M167-dex40-04032025`

## Immediate QC Tasks

### Inspect M167 HR

`M167-dex40-04032025` has all HR bins marked dubious.

Next checks:

1. Plot short ECG subsections with accepted R-peaks.
2. Inspect the distribution of RR intervals and `hr_inst_bpm`.
3. Check whether the issue is:
   - bad R-peak detection,
   - wrong ECG/image time alignment,
   - overly strict HR sanity bounds,
   - a real physiological or recording anomaly.

Decision after QC:

- exclude M167 from all HR analyses,
- fix/relax HR QC for this subject,
- or regenerate R-peak variables if the detections are bad.

### Inspect Other Outliers

Use `results/hr_amp_corr_summary.csv` and the time-series plots under `figures/hr_amp_coupling/` to identify subjects or variables with unusually strong correlations.

For any strong effect:

1. Check whether it is driven by one or two time bins.
2. Inspect the corresponding HR and vessel traces.
3. Check whether large `amp_svm` or `bp_svm` peaks coincide with visible vessel artifacts.

## Cross-Correlation Script

The main remaining analysis script is the all-subject cross-correlation study.

Create or update:

- `compute_hr_diam_cross_corr.m`

Goal:

Compute lagged cross-correlations between `hr_median_bpm` and each vessel feature from `hr_amp_data.mat`.

Inputs:

- `data/*/hr_amp_data.mat`

Use the same exclusion rule:

- exclude `dubious_hr_bin_fraction > 0.50`
- display excluded subjects in the command window
- save excluded subjects to CSV

Vessel variables:

- `amp_svm`
- `bp_svm`
- `amp_c`
- `diam_mean`

Recommended lag setup:

- 5-minute bins.
- Start with `max_lag_bins = 6`, equivalent to plus/minus 30 minutes.
- Consider `max_lag_bins = 12`, equivalent to plus/minus 60 minutes, as sensitivity analysis.

Output table should include:

- subject
- HR variable
- vessel variable
- lag values in minutes
- normalized cross-correlation at each lag
- peak positive correlation and lag
- peak absolute correlation and lag
- `dubious_hr_bin_fraction`

Suggested CSV outputs:

- `results/hr_amp_xcorr_summary.csv`
- `results/hr_amp_xcorr_excluded_subjects.csv`

Suggested figures:

- one cross-correlation plot per subject and vessel variable under `figures/hr_amp_xcorr/`
- optional group-level boxplots or mean curves by vessel variable

Important convention to document:

- Define clearly whether positive lag means HR leads vessel or vessel leads HR.
- Put that convention in plot labels and table metadata.

## Correlation Refinements

Potential refinements for `compute_hr_diam_corr.m`:

1. Add Pearson correlations as optional secondary outputs.
2. Add sensitivity analysis excluding dubious HR bins entirely, instead of using interpolated HR bins.
3. Add robust scatter plots for each HR-vessel pair.
4. Add group-level summary statistics for each vessel variable and time window.

## Plotting Refinements

Potential refinements:

1. Add individual subject points overlaid on the correlation boxplots.
2. Label excluded subjects directly in figure captions or sidecar CSVs.
3. Add separate folders for boxplots and summary figures if `figures/` becomes crowded.
4. Generate time-series plots only for selected variables if the full batch becomes too noisy to browse.

## Documentation Updates

Update the pipeline guide after the cross-correlation script is finalized:

- `HR_Diameter_Pipeline_Guide.md`

Add sections for:

- all-subject `hr_amp_data.mat` extraction,
- final variable names, especially `diam_mean`,
- HR QC/exclusion rule,
- correlation script outputs,
- cross-correlation script outputs and lag convention.

## Preferred Next Coding Order

1. QC M167 R-peaks and HR bins.
2. Update `compute_hr_diam_cross_corr.m` to use `hr_amp_data.mat`.
3. Run cross-correlation for all included subjects.
4. Generate cross-correlation figures and summaries.
5. Add group-level visualization for cross-correlation results.
6. Update `HR_Diameter_Pipeline_Guide.md` with the finalized analysis workflow.
