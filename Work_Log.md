# Work Log

## 2026-05-06 HR Agreement and Legacy-Yue Follow-Up

### Goal

We extended the recent HR correlation work in three directions:

- compare our rebinned beat-to-beat ECG HR against legacy `avg_HR_Yue_10min`,
- compare both ECG-HR summaries against legacy vessel-derived `results.hr`,
- rerun the HR-vessel feature correlation workflow using `avg_HR_Yue_10min` instead of `hr_median_bpm`.

### HR-HR Comparison with Rebinned ECG HR

A new folder was created:

- `hr_hr_correlation/`

New script:

- `hr_hr_correlation/compute_hr_hr_corr.m`

This script:

- loads `hr_amp_data.mat` and `results.mat` for each mouse,
- rebins beat-to-beat `hr_inst_bpm` onto the 10-minute `results.timepts` windows,
- compares that rebinned ECG HR against `results.hr`,
- computes per-mouse Pearson and Spearman correlations for early, late, and overall windows,
- runs Fisher `z` one-sample t-tests across mice,
- excludes subjects with `dubious_hr_bin_fraction > 0.50`.

Outputs:

- `hr_hr_correlation/hr_hr_corr_by_mouse.csv`
- `hr_hr_correlation/hr_hr_corr_ttest_summary.csv`
- `hr_hr_correlation/hr_hr_corr_results.xlsx`
- `hr_hr_correlation/hr_hr_corr_excluded_subjects.csv`
- `hr_hr_correlation/figures/hr_hr_corr_ttest_summary.png`

Result summary:

- this ECG-HR vs vessel-HR comparison produced significant positive group-level Fisher `z` t-test results for both Pearson and Spearman in the late and overall windows,
- `M167-dex40-04032025` was excluded again for 100% dubious HR bins.

### HR Agreement with `avg_HR_Yue_10min`

New script:

- `hr_hr_correlation/compute_hr_agreement_yue.m`

This script:

- compares rebinned beat-to-beat ECG HR against `results.avg_HR_Yue_10min`,
- computes per-mouse Pearson and Spearman correlations,
- computes Fisher `z` group summaries,
- writes agreement metrics including overlap count, mean difference, mean absolute difference, and RMSE,
- saves one time-series comparison figure per mouse.

Outputs:

- `hr_hr_correlation/hr_agreement_yue_by_mouse.csv`
- `hr_hr_correlation/hr_agreement_yue_ttest_summary.csv`
- `hr_hr_correlation/hr_agreement_yue_results.xlsx`
- `hr_hr_correlation/hr_agreement_yue_excluded_subjects.csv`
- `hr_hr_correlation/figures/hr_agreement_yue_ttest_summary.png`
- `hr_hr_correlation/figures/timeseries_by_mouse/*_hr_agreement.png`

Interpretation:

- across the 10 included mice, `avg_HR_Yue_10min` generally agreed strongly with the rebinned beat-to-beat ECG HR,
- some individual mice still showed sizable offsets or weaker agreement, which became easier to inspect in the overlay plots.

### Legacy Yue HR vs Vessel HR

New script:

- `hr_hr_correlation/compute_hr_hr_corr_yue_vs_vessel.m`

This script:

- compares `results.avg_HR_Yue_10min` against `results.hr`,
- computes early, late, and overall Pearson and Spearman correlations by mouse,
- runs Fisher `z` one-sample t-tests across mice,
- writes the same table and figure structure as the first HR-HR workflow.

Outputs:

- `hr_hr_correlation/hr_hr_corr_yue_vs_vessel_by_mouse.csv`
- `hr_hr_correlation/hr_hr_corr_yue_vs_vessel_ttest_summary.csv`
- `hr_hr_correlation/hr_hr_corr_yue_vs_vessel_results.xlsx`
- `hr_hr_correlation/hr_hr_corr_yue_vs_vessel_excluded_subjects.csv`
- `hr_hr_correlation/figures/hr_hr_corr_yue_vs_vessel_ttest_summary.png`

Result summary:

- this Yue-HR vs vessel-HR analysis also showed positive group-level correlations across mice,
- the Yue-based comparison was somewhat cleaner than the first rebinned-HR vs vessel-HR run.

### Yue-Based HR-Vessel Correlation Workflow

Two new parallel folders were created:

- `pearson_corr_yue/`
- `spearman_correlation_yue/`

New compute scripts:

- `pearson_corr_yue/compute_hr_diam_corr_pearson_yue.m`
- `spearman_correlation_yue/compute_hr_diam_corr_spearman_yue.m`

These workflows:

- use `results.avg_HR_Yue_10min` as the heart-rate variable,
- use legacy 10-minute `results.mat` vessel features: `amp_svm`, `bp_svm`, `amp_c`, and `diam_tp`,
- keep the same mouse-level unit-of-inference convention,
- reuse the same dubious-HR exclusion rule from `hr_amp_data.mat`,
- write by-mouse tables, Wilcoxon summaries, excluded-subject lists, and Excel workbooks.

Outputs:

- `pearson_corr_yue/hr_amp_corr_yue_pearson_results.xlsx`
- `pearson_corr_yue/hr_amp_corr_yue_pearson_by_mouse.csv`
- `pearson_corr_yue/hr_amp_corr_yue_pearson_wilcoxon_summary.csv`
- `pearson_corr_yue/hr_amp_corr_yue_pearson_excluded_subjects.csv`
- `spearman_correlation_yue/hr_amp_corr_yue_spearman_results.xlsx`
- `spearman_correlation_yue/hr_amp_corr_yue_spearman_by_mouse.csv`
- `spearman_correlation_yue/hr_amp_corr_yue_spearman_wilcoxon_summary.csv`
- `spearman_correlation_yue/hr_amp_corr_yue_spearman_excluded_subjects.csv`

### Yue-Based Diagnostic and T-Test Plots

New plotting and t-test scripts:

- `pearson_corr_yue/compute_hr_diam_corr_pearson_yue_ttest.m`
- `spearman_correlation_yue/compute_hr_diam_corr_spearman_yue_ttest.m`
- `pearson_corr_yue/plot_hr_diam_corr_pearson_yue_diagnostics.m`
- `spearman_correlation_yue/plot_hr_diam_corr_spearman_yue_diagnostics.m`
- `pearson_corr_yue/plot_hr_diam_corr_pearson_yue_ttest_summary.m`
- `spearman_correlation_yue/plot_hr_diam_corr_spearman_yue_ttest_summary.m`

These scripts produce:

- per-mouse scatter diagnostics under each folder’s `figures/scatter_by_mouse/`,
- Wilcoxon summary figures under each folder’s `figures/wilcoxon_summary/`,
- Fisher `z` t-test summary figures under each folder’s `figures/ttest_summary/`.

Main result:

- the Yue-based HR-vessel feature analysis gave broadly similar conclusions to the earlier `hr_median_bpm` workflow,
- no group-level Pearson or Spearman HR-vessel feature test crossed `p < 0.05`,
- the closest tendency was a negative late-window `amp_svm` effect in the Spearman t-test, but it remained above threshold.

## 2026-05-01 Correlation Inference Session

### Goal

We extended the HR-vessel correlation analysis from per-mouse correlation tables to group-level inference across mice.

The main questions were:

- whether Pearson correlation should be compared against Spearman,
- how to test whether correlations are consistently nonzero across mice,
- how to visualize and interpret the resulting group-level uncertainty.

### Pearson Correlation Workflow

A dedicated Pearson analysis folder was created:

- `pearson_corr/`

New script:

- `pearson_corr/compute_hr_diam_corr_pearson.m`

This script:

- loads `data/*/hr_amp_data.mat`,
- uses `hr_median_bpm` as the HR feature,
- computes per-mouse Pearson correlations for `amp_svm`, `bp_svm`, `amp_c`, and `diam_mean`,
- reports early, late, and overall windows,
- excludes subjects with `dubious_hr_bin_fraction > 0.50`,
- writes both CSV and Excel outputs.

Pearson outputs:

- `pearson_corr/hr_amp_corr_pearson_results.xlsx`
- `pearson_corr/hr_amp_corr_pearson_by_mouse.csv`
- `pearson_corr/hr_amp_corr_pearson_wilcoxon_summary.csv`
- `pearson_corr/hr_amp_corr_pearson_excluded_subjects.csv`

### Spearman Comparison Workflow

A matching comparison folder was created:

- `spearman_correlation/`

New script:

- `spearman_correlation/compute_hr_diam_corr_spearman.m`

This mirrors the Pearson workflow but uses Spearman correlation and writes parallel outputs:

- `spearman_correlation/hr_amp_corr_spearman_results.xlsx`
- `spearman_correlation/hr_amp_corr_spearman_by_mouse.csv`
- `spearman_correlation/hr_amp_corr_spearman_wilcoxon_summary.csv`
- `spearman_correlation/hr_amp_corr_spearman_excluded_subjects.csv`

### Group-Level Wilcoxon Tests

For both Pearson and Spearman folders, the `Wilcoxon` summary sheet and CSV test the mouse-level correlations across animals.

Important inference convention:

- the Wilcoxon tests do **not** pool 5-minute bins across mice,
- the unit of inference is mouse,
- each test uses one correlation value per mouse for a given vessel feature and time window.

Interpretation:

- per-mouse `p` values answer whether a within-mouse correlation is nonzero,
- Wilcoxon `p` values answer whether the correlation is consistently above or below zero across mice.

### Scatter Diagnostics and Wilcoxon Plots

New plotting scripts:

- `pearson_corr/plot_hr_diam_corr_pearson_diagnostics.m`
- `spearman_correlation/plot_hr_diam_corr_spearman_diagnostics.m`

These create:

- per-mouse HR-vs-feature scatter plots under `figures/scatter_by_mouse/`,
- group-level Wilcoxon summary figures under `figures/wilcoxon_summary/`.

The scatter plots:

- separate early and late bins by color,
- show a linear fit guide,
- annotate each panel with the per-mouse correlation summary.

The Wilcoxon summary plots:

- show one point per mouse,
- show a zero reference line,
- show the median correlation line,
- annotate each panel with the Wilcoxon `p` value.

### Fisher Z One-Sample T-Tests

After reviewing the rationale for Fisher `z`, we added one-sample t-tests on the mouse-level correlations without rebuilding the underlying correlation pipeline.

New scripts:

- `pearson_corr/compute_hr_diam_corr_pearson_ttest.m`
- `spearman_correlation/compute_hr_diam_corr_spearman_ttest.m`

These scripts:

- read the existing per-mouse correlation tables,
- apply `z = atanh(r)` after clipping `r` away from `-1` and `1`,
- run one-sample t-tests of the Fisher-transformed correlations against `0`,
- back-transform the mean and confidence interval to correlation units,
- append a `TTest` sheet to the existing Excel workbooks,
- write CSV summaries.

T-test outputs:

- `pearson_corr/hr_amp_corr_pearson_ttest_summary.csv`
- `spearman_correlation/hr_amp_corr_spearman_ttest_summary.csv`

### T-Test Summary Plots

New plotting scripts:

- `pearson_corr/plot_hr_diam_corr_pearson_ttest_summary.m`
- `spearman_correlation/plot_hr_diam_corr_spearman_ttest_summary.m`

These produce:

- `pearson_corr/figures/ttest_summary/pearson_ttest_summary.png`
- `spearman_correlation/figures/ttest_summary/spearman_ttest_summary.png`

Interpretation of the plotted red lines:

- solid red line: back-transformed mean correlation,
- red dotted lines: 95% confidence interval bounds computed in Fisher `z` space and transformed back to correlation units,
- black dashed line: zero-correlation reference.

The t-test summary plots were later refined to:

- reduce axis-label size,
- shorten x-axis labels to compact subject IDs like `F168`,
- force subplot titles to black,
- omit the `diam_mean` row from the t-test figures while leaving `diam_mean` in the tables.

### Main Statistical Result

Across the 10 included mice, no group-level test crossed the conventional `p < 0.05` threshold.

This held for:

- Pearson + Wilcoxon,
- Pearson + Fisher `z` t-test,
- Spearman + Wilcoxon,
- Spearman + Fisher `z` t-test.

Closest trends were in the early window, especially:

- Spearman `bp_svm`,
- Spearman `diam_mean`,
- Spearman `amp_c`.

These remained suggestive rather than statistically significant.

### Recommended Interpretation

A reasonable short conclusion from the current analysis is:

- individual mice can show significant within-mouse HR-vessel correlations,
- but the direction and magnitude are not yet consistent enough across animals to support a significant group-level effect.

## 2026-05-13 Filtered ECG HR Update

Added a confidence-filtered ECG HR branch to the all-subject extractors:

- `extract_hr_amp_data_all_subjects.m`
- `extract_hr_amp_data_all_subjects_sliding.m`

The new branch:

- builds 1-minute ECG HR bins from accepted beat-to-beat HR,
- marks a 1-minute bin bad when at least 50% of candidate R-peaks in that minute have `r_peak_confidence < 0.5`,
- leaves bad 1-minute bins uncalculated initially,
- linearly interpolates or extrapolates missing 1-minute HR values from valid minute bins,
- saves the filtered minute series as `t_hr_filtered`, `hr_inst_bpm_filtered`, and `valid_inst_hr_filtered`,
- also saves 5-minute or sliding-window summaries: `hr_mean_bpm_filtered`, `hr_median_bpm_filtered`, and `hr_iqr_bpm_filtered`.

The original HR variables are still saved unchanged for comparison:

- `t_hr`
- `hr_inst_bpm`
- `valid_inst_hr`
- `hr_median_bpm`

Added:

- `hr_hr_correlation/compute_hr_hr_corr_filtered.m`

This script rebins `hr_inst_bpm_filtered` onto the legacy 10-minute `results.mat` windows and correlates it against `results.hr` converted from Hz to BPM.

Generated outputs:

- `hr_hr_correlation/hr_hr_corr_filtered_by_mouse.csv`
- `hr_hr_correlation/hr_hr_corr_filtered_ttest_summary.csv`
- `hr_hr_correlation/hr_hr_corr_filtered_excluded_subjects.csv`
- `hr_hr_correlation/hr_hr_corr_filtered_results.xlsx`
- `hr_hr_correlation/figures/hr_hr_corr_filtered_ttest_summary.png`
- per-subject overlays in `hr_hr_correlation/figures/timeseries_filtered_by_mouse/`

Run results:

- Included subjects: 9
- Excluded by filtered HR QC: `F168-dex40-10032025`, `M167-dex40-04032025`
- Overall group-level filtered ECG HR vs `results.hr` correlations:
  - Spearman Fisher-z t-test: back-transformed mean `r = 0.755`, `p = 0.0146`
  - Pearson Fisher-z t-test: back-transformed mean `r = 0.854`, `p = 0.0135`

Added `hr_hr_correlation/compute_hr_agreement_yue_filtered.m` to compare filtered ECG HR against `avg_HR_Yue_10min`.
It saves Fisher-z t-test and Wilcoxon summaries plus per-subject diagnostic plots with:

- filtered ECG HR vs Yue scatter,
- filtered ECG HR trace,
- `avg_HR_Yue_10min` trace,
- `results.hr * 60` vessel-HR trace.

Run results:

- Included subjects: 9
- Excluded by filtered HR QC: `F168-dex40-10032025`, `M167-dex40-04032025`
- Overall filtered ECG HR vs `avg_HR_Yue_10min`:
  - Spearman Fisher-z t-test: back-transformed mean `r = 0.902`, `p = 0.00382`; Wilcoxon `p = 0.00781`
  - Pearson Fisher-z t-test: back-transformed mean `r = 0.946`, `p = 0.00167`; Wilcoxon `p = 0.00781`

## 2026-04-22 Follow-Up Session

### All-Subject Extraction

All mouse subject folders were updated with raw `signals.mat` and `im.mat`, then the F153-specific compact extractor was generalized into:

- `extract_hr_amp_data_all_subjects.m`

The all-subject extractor:

- loops over every subject under `data/`,
- aligns ECG time using `t_ecgoffset`,
- accepts R-peaks using `r_peak_confidence >= 0.5` and `good_r_peaks`,
- computes beat-to-beat `hr_inst_bpm`,
- bins HR and vessel features into 5-minute windows starting at 14 minutes,
- computes `amp_svm`, `bp_svm`, `amp_c`, and `diam_mean`,
- tracks `is_dubious_hr_bin` and `dubious_hr_bin_fraction`,
- saves `hr_amp_data.mat` inside each subject folder,
- writes `results/hr_amp_extraction_summary.csv`.

Compatibility fallbacks were added for variation across subjects:

- use `t_frames` when available, otherwise fall back to `t_vect`;
- use `diam` when available, otherwise fall back to `diam_dlc`.

All 11 subjects now have `hr_amp_data.mat`.

Important QC result:

- `M167-dex40-04032025` saved successfully using `diam_dlc`, but has `33 / 33` dubious HR bins, so it should be excluded from HR analyses until its R-peaks are inspected or the HR QC rule is revised.
- `M163-dex40-24022025` has `1 / 36` dubious HR bins.
- All other subjects had `0` dubious HR bins.

### Variable Rename

The old variable name `diam_tp` was renamed to:

- `diam_mean`

This is clearer because the variable is the mean raw vessel diameter in each bin, in microns. All regenerated `hr_amp_data.mat` files now contain `diam_mean` and no longer contain `diam_tp`.

### Correlation and Time-Series Plotting

The existing scripts were updated to consume `hr_amp_data.mat` instead of legacy `results.mat`:

- `compute_hr_diam_corr.m`
- `plot_hr_diam_amp_coupling.m`

`compute_hr_diam_corr.m` now:

- uses `hr_median_bpm` as the default heart-rate feature,
- computes Spearman correlations against `amp_svm`, `bp_svm`, `amp_c`, and `diam_mean`,
- reports early, late, and overall correlations in one table,
- excludes subjects with `dubious_hr_bin_fraction > 0.50`,
- writes `results/hr_amp_corr_summary.csv`,
- writes `results/hr_amp_corr_excluded_subjects.csv`.

`plot_hr_diam_amp_coupling.m` now:

- plots HR and each vessel feature from `hr_amp_data.mat`,
- excludes subjects with `dubious_hr_bin_fraction > 0.50`,
- saves figures under `figures/hr_amp_coupling/`,
- writes an excluded-subject list.

Both scripts were sanity-run successfully. `M167-dex40-04032025` was excluded due to 100% dubious HR bins.

### Correlation Boxplots

A new script was added:

- `plot_hr_amp_corr_boxplots.m`

It reads `results/hr_amp_corr_summary.csv` and saves one boxplot per vessel variable directly under `figures/`. Each figure shows three boxes across subjects:

- Early
- Late
- Overall

Generated figures:

- `figures/hr_corr_boxplot_amp_svm.png`
- `figures/hr_corr_boxplot_bp_svm.png`
- `figures/hr_corr_boxplot_amp_c.png`
- `figures/hr_corr_boxplot_diam_mean.png`

## 2026-04-22 Session

### Project Context

We worked on the mouse ECG and vessel-diameter analysis project. The main scientific goal is to study relationships between ECG-derived heart rate and vessel diameter dynamics, especially slow vasomotion amplitude/power and cardiac-band vessel pulsation.

### Pipeline Understanding

We reviewed the existing analysis flow, primarily from `PlotResults_continuous_v5.m`, and documented it in:

- `HR_Diameter_Pipeline_Guide.md`

Key points captured:

- `signals.mat` contains ECG, `t_ECG`, candidate R-peaks, confidence values, and prefiltered good R-peak variables.
- `im.mat` contains vessel diameter, frame times, pixel-to-micron scale, and ECG/image time offset.
- ECG time is aligned to imaging time with `t_ECG = t_ECG - t_ecgoffset`.
- Vessel diameter is detrended with a 10-minute moving average.
- Legacy slow vasomotion feature extraction uses `Fc_svm = [0.01 0.3]`, implemented as a 0.3 Hz lowpass after detrending.
- Legacy cardiac-band diameter extraction for dex subjects uses `Fc_c = [3.5 15]`.
- `results.hr` is not ECG heart rate. It is the vessel-diameter PSD peak frequency in the cardiac band.
- ECG-derived binned HR in the legacy output is stored as `avg_HR_Yue_10min`.

### Fine-Resolution Pilot

We decided to use 5-minute bins as the current fine-resolution choice.

Reasoning:

- The slowest vasomotion frequency of interest is 0.01 Hz.
- A 0.01 Hz oscillation has a 100-second period.
- A 5-minute bin contains about 3 cycles at 0.01 Hz.
- This is a reasonable lower practical limit for IQR and band-power summaries while improving time resolution for cross-correlation.

We created:

- `pilot_fine_hr_diam_features_F153.m`

This script runs on `data/F153-Dex40-06022025/` and:

- Applies `r_peak_confidence >= 0.5`.
- Requires candidate peaks to also be in `good_r_peaks`.
- Computes beat-to-beat instantaneous HR from RR intervals.
- Uses broad mouse anesthesia sanity bounds of 250-750 BPM.
- Computes 5-minute binned vessel features starting at 14 minutes.
- Saves a pilot `results_fine.mat`.
- Plots a short ECG subsection with accepted R-peaks as red dashed vertical lines and instantaneous HR on the same time axis.

The pilot plot was saved to:

- `figures/F153-Dex40-06022025_pilot_ecg_rpeaks_inst_hr.png`

### Compact HR/Amplitude Dataset

We created the non-plotting extraction script:

- `extract_hr_amp_data_F153.m`

It writes:

- `data/F153-Dex40-06022025/hr_amp_data.mat`

The compact MAT file includes:

- `t_rpeaks`
- `r_peak_inds`
- `t_hr`
- `hr_inst_bpm`
- `valid_inst_hr`
- `hr_mean_bpm`
- `hr_median_bpm`
- `hr_iqr_bpm`
- raw versions of the HR bin summaries before interpolation
- `is_dubious_hr_bin`
- `dubious_hr_bin_fraction`
- `n_beats_bin`
- `amp_svm`
- `bp_svm`
- `amp_c`
- `diam_tp`
- bin metadata and frequency-band metadata

We intentionally excluded:

- `amp_noise`
- `bp_svm2`
- `amp_svm2`

For F153, the run reported:

```text
Dubious HR bins interpolated: 0 / 33 (0.000)
```

### Plotting

We restored:

- `plot_hr_diam_amp_coupling.m`

from `archive/`.

We created:

- `plot_hr_amp_svm_fine_F153.m`

This script loads `hr_amp_data.mat` and plots:

- `hr_median_bpm`
- `amp_svm`

on the same time axis with dual y-axes.

It saves:

- `figures/F153-Dex40-06022025_fine_hr_amp_svm.png`

We compared this 5-minute-bin plot against the old 10-minute-bin `amp_svm` plot. Both showed a vessel `amp_svm` spike around the same time window, supporting the idea that the spike is in the data rather than caused by the new 5-minute binning.

### Project Cleanup

We moved older exploratory MATLAB scripts and the older coupling report into:

- `archive/`

The project root now keeps the active working files:

- `compute_hr_diam_corr.m`
- `compute_hr_diam_cross_corr.m`
- `extract_hr_amp_data_F153.m`
- `HR_Diameter_Pipeline_Guide.md`
- `pilot_fine_hr_diam_features_F153.m`
- `PlotResults_continuous_v5.m`
- `plot_hr_amp_svm_fine_F153.m`
- `plot_hr_diam_amp_coupling.m`

### Notes and Caveats

- The compact extraction script is currently subject-specific for F153.
- Several subject folders currently have only `results.mat`; all-subject fine extraction needs `signals.mat` and `im.mat` in every subject folder.
- The next correlation scripts should use ECG-derived HR variables from `hr_amp_data.mat`, not legacy `results.hr`.
- The 5-minute bins are a good current compromise for slow vasomotion and cross-correlation.

## 2026-05-20 Session

### Confidence-Filtered Minute HR

Updated the all-subject HR extraction logic so the confidence-filtered 1-minute ECG HR no longer fills rejected bins. Bins where at least 50% of R-peaks have confidence below 0.5 remain `NaN` in `hr_minute_bpm_filtered`, and the filtered binned HR summaries are no longer interpolated or extrapolated.

Updated scripts:

- `extract_hr_amp_data_all_subjects.m`
- `extract_hr_amp_data_all_subjects_sliding.m`

Generated a three-way HR comparison plot workflow:

- `plot_three_hr_comparison.m`

The comparison plots overlay:

- filtered 1-minute ECG HR with `NaN` gaps left unconnected
- `avg_HR_Yue_10min`
- vessel-derived `results.hr` converted from Hz to BPM

The generated plots were saved locally under `three_hr_comparison/`.

### Overlapping Vessel Features

Added:

- `extract_vessel_features_10min_overlap.m`

This computes `amp_svm`, `bp_svm`, and `amp_c` using the legacy `PlotResults_continuous_v5.m` definitions, but with 10-minute windows advanced every 1 minute. Each subject gets a local `vessel_features_10min_overlap.mat` file.

Validation showed that values at the old legacy 10-minute start times match `results.mat` to numerical precision:

```text
amp_svm max error: 4.44e-16
bp_svm max error: 1.78e-15
amp_c max error: 6.66e-16
```

### Overlap Correlation Analyses

Created overlap correlation analysis folders:

- `pearson_corr_overlap/`
- `spearman_correlation_overlap/`

Each folder contains scripts to compute by-mouse correlations, Wilcoxon summaries, Fisher-z t-test summaries, and diagnostic figures. The overlap correlations align filtered 1-minute ECG HR to the 10-minute vessel-feature window start times, using complete finite pairs only.

Tracked scripts include:

- `pearson_corr_overlap/compute_hr_diam_corr_pearson_overlap.m`
- `pearson_corr_overlap/compute_hr_diam_corr_pearson_overlap_ttest.m`
- `pearson_corr_overlap/plot_hr_diam_corr_pearson_overlap_diagnostics.m`
- `pearson_corr_overlap/plot_hr_diam_corr_pearson_overlap_ttest_summary.m`
- `spearman_correlation_overlap/compute_hr_diam_corr_spearman_overlap.m`
- `spearman_correlation_overlap/compute_hr_diam_corr_spearman_overlap_ttest.m`
- `spearman_correlation_overlap/plot_hr_diam_corr_spearman_overlap_diagnostics.m`
- `spearman_correlation_overlap/plot_hr_diam_corr_spearman_overlap_ttest_summary.m`

Initial overlap correlation runs excluded subjects with more than 50% aligned filtered-HR `NaN` bins. Later, the overlap Spearman analysis was rerun with manual HR-disagreement exclusions requested by the user:

- `F168-dex40-10032025`
- `F169-dex40-10032025`
- `M163-Dex40-26022025`
- `M163-dex40-24022025`
- `M166-dex40-04032025`
- `M167-dex40-04032025`

That final Spearman-overlap run retained:

- `F153-Dex40-06022025`
- `F154-dex40-19022025`
- `F155-dex40-16022025`
- `F156-Dex40-04022025`
- `F157-dex40-16022025`

The strongest recurring pattern remained a negative relationship between filtered ECG HR and `amp_svm`, especially in the early window.
