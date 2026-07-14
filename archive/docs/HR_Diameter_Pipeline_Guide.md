# Heart Rate - Vessel Diameter Pipeline Guide

This guide summarizes how the project currently extracts ECG-derived heart-rate features and vessel-diameter features for downstream correlation and cross-correlation analyses. It is based mainly on `PlotResults_continuous_v5.m`, with `compute_hr_diam_corr.m`, `compute_hr_diam_cross_corr.m`, `analyze_hr_diam_coupling.m`, `preprocess_vessel_diam.m`, and one inspected subject folder used as context.

The main study goal is to quantify how instantaneous or binned heart rate relates to vessel diameter dynamics across frequency bands.

## Data Layout

Each mouse subject lives in a subfolder under `data/`.

Important files:

- `signals.mat`: ECG, respiration, ECG time base, R-peak detections, and R-peak quality variables.
- `im.mat`: vessel diameter, frame time base, spatial scaling, and ECG-to-image time offset.
- `results.mat`: derived variables produced by the processing/plotting pipeline.

In the inspected subject `data/F154-dex40-19022025`, the key raw variables are:

- From `signals.mat`
  - `t_ECG`: ECG sample times.
  - `ECG`: raw ECG signal.
  - `resp`: respiration signal.
  - `detected_r_peaks`: candidate R-peak sample indices into `t_ECG`.
  - `r_peak_confidence`: confidence value for each candidate R-peak.
  - `good_r_peaks`: prefiltered R-peak indices, apparently confidence-thresholded.
  - `pk_inds_keep`: two-column list of accepted adjacent peak pairs used for RR intervals.
  - `pk_inds_isvalid`: logical validity mask for accepted peaks.
- From `im.mat`
  - `diam`: vessel diameter in pixels.
  - `t_frames`: imaging frame times.
  - `t_vect`: older/alternate frame time vector.
  - `umperpix`: micron-per-pixel scale factor.
  - `t_ecgoffset`: offset used to align ECG time to imaging time.
  - `diam_dlc`: optional alternate diameter trace when DLC-derived diameter is selected.

Some subject folders currently contain only `results.mat` and not raw `signals.mat` or `im.mat`. Raw extraction or reprocessing therefore requires checking file availability subject by subject.

## Time Alignment

The ECG and imaging clocks are aligned by subtracting the stored offset:

```matlab
t_ECG = t_ECG - t_ecgoffset;
```

After this correction, `t_ECG` and `t_frames` should be in the same time coordinate. Most feature extraction then bins both ECG and diameter signals using absolute time in seconds or minutes.

The legacy script prefers `t_frames` when available:

```matlab
if exist('t_frames','var')
    t_vect = t_frames';
end
```

Downstream analyses should use `t_frames` consistently when it exists.

## Vessel Diameter Preprocessing

### Legacy Processing in `PlotResults_continuous_v5.m`

The core legacy steps are:

1. Load `diam`, `t_frames` or `t_vect`, and `umperpix`.
2. Optionally replace `diam` with `diam_dlc` if the data list flags DLC use.
3. Estimate the imaging sampling rate:

```matlab
fs_median = 1 / median(diff(t_vect));
```

4. Compute a slow moving-average baseline using a 10-minute window:

```matlab
diam_avg = smoothdata(diam, 'movmean', 10 * 60 * fs_median);
```

5. Detrend by subtracting `diam_avg`:

```matlab
diam_det = diam - diam_avg;
```

6. Filter the detrended signal into frequency bands.

The script keeps the filtered traces in pixel units and multiplies summary amplitudes by `umperpix` to report microns.

### Cleaner Normalized Processing in `preprocess_vessel_diam.m`

`preprocess_vessel_diam.m` is a cleaner helper that:

1. Converts diameter from pixels to microns.
2. Flags large frame-to-frame jumps, NaNs, and nonpositive diameters.
3. Interpolates over flagged frames.
4. Removes a moving-average baseline.
5. Returns normalized percent diameter change:

```matlab
diam_norm = (diam_clean - baseline) ./ baseline * 100;
```

This normalized representation is useful when comparing subjects with different vessel diameters. The legacy `results.mat` variables, however, mostly use absolute micron units for summary amplitudes.

## Vessel Diameter Frequency Bands

`PlotResults_continuous_v5.m` defines these primary bands:

- Slow vasomotion, broad lowpass:

```matlab
Fc_svm = [0.01 0.3];
diam_filt_svm = lowpass(diam - diam_avg, Fc_svm(2), fs_median);
```

Note: even though `Fc_svm` stores `[0.01 0.3]`, the actual legacy filter is a lowpass at `0.3 Hz`, so it includes components below `0.01 Hz` after detrending.

- Slow vasomotion, narrower bandpass:

```matlab
Fc_svm2 = [0.05 0.3];
diam_filt_svm2 = bandpass(diam - diam_avg, Fc_svm2, fs_median);
```

- Cardiac diameter pulsation:

```matlab
if KX
    Fc_c = [1.5 10];
elseif dex
    Fc_c = [3.5 15];
elseif PPA or PINP
    Fc_c = [1.5 15];
else
    Fc_c = [6 15];
end
diam_filt_c = bandpass(diam - diam_avg, Fc_c, fs_median);
```

- High-frequency noise:

```matlab
Fc_noise = 15;
diam_filt_noise = highpass(diam - diam_avg, Fc_noise, fs_median);
```

For the current dex40 data, the cardiac vessel band is expected to be `[3.5 15] Hz` in the legacy logic.

## Binned Diameter Features

The legacy script summarizes data in 10-minute bins:

```matlab
segtime = 10; % minutes
timepts = [14 25:segtime:t_vect(end)/60-segtime];
```

This skips the early recording and starts bins at 14 minutes, then 25, 35, 45, etc. For each bin, it computes:

- `diam_tp`: mean vessel diameter in microns.
- `amp_svm`: IQR amplitude of `diam_filt_svm`, in microns.
- `amp_svm2`: IQR amplitude of `diam_filt_svm2`, in microns.
- `amp_c`: IQR amplitude of cardiac-band diameter, in microns.
- `amp_noise`: IQR amplitude of high-frequency noise, in microns.
- `bp_svm`: Welch PSD band power in the slow vasomotion band.
- `bp_svm2`: Welch PSD band power in the narrower slow vasomotion band.
- `bpp_svm`: slow-band power divided by total power above `Fc_svm(1)`.
- `svm_tails`: difference between the 99.9th and 95th percentiles of the slow diameter trace.
- `hr`: peak frequency of the vessel diameter PSD inside `Fc_c`. Despite its name, this is not ECG-derived heart rate.
- `w`: a rough width/variability measure of the cardiac-band PSD peak.

The most useful diameter targets for HR-vessel coupling are likely:

- `amp_svm`: slow vasomotion amplitude.
- `amp_svm2`: slow vasomotion amplitude after excluding ultra-slow residuals.
- `bp_svm` or `bp_svm2`: slow-band power.
- `amp_c`: cardiac pulsation amplitude in the vessel wall.
- `diam_tp`: mean vessel diameter.

Be careful with `hr`: in `results.mat`, `hr` is the vessel-spectrum peak frequency, not the ECG heart rate.

## ECG R-Peak Quality and Heart Rate

### Candidate Peaks and Confidence

The raw R-peak information is stored as:

- `detected_r_peaks`: candidate ECG sample indices.
- `r_peak_confidence`: confidence score for each candidate.

The inspected data also include `good_r_peaks`, which appears to already apply a confidence threshold, and `pk_inds_keep`, which stores accepted adjacent R-peak pairs.

For a confidence-based reconstruction from raw candidates, the natural starting point is:

```matlab
peak_mask = r_peak_confidence >= threshold;
good_peaks = detected_r_peaks(peak_mask);
t_rpeaks = t_ECG(good_peaks);
```

The current scripts often use the already-saved `good_r_peaks` or `pk_inds_keep` instead of rebuilding from `r_peak_confidence`.

### Beat-to-Beat Instantaneous Heart Rate

The legacy script computes beat-to-beat HR from accepted adjacent peak pairs:

```matlab
hr_inst = 1 ./ (t_ECG(pk_inds_keep(:,2)) - t_ECG(pk_inds_keep(:,1)));
```

This is in Hz. Multiply by 60 for BPM.

This is a true instantaneous RR-based heart-rate estimate, but in the current legacy script it is mainly used to compute alternative binned summaries:

```matlab
HR_avg(i) = mean(hr_inst(ind_HRpeaks));
HR_med(i) = median(hr_inst(ind_HRpeaks));
HR_var(i) = iqr(hr_inst(ind_HRpeaks));
```

Those variables are not saved in the inspected `results.mat`, while `avg_HR_Yue_10min`, `med_HR`, and `var_HR` are saved.

### Yue/Andy 10-Minute HR Summary

The preferred binned HR in `PlotResults_continuous_v5.m` is:

```matlab
[avg_HR_Yue_10min(i), med_HR(i), var_HR(i)] = ...
    Cal_avg_HeartRate_Yue(ind_ECG, detected_r_peaks, fs_ECG, peak_mask);
```

Here `peak_mask` is made by marking candidate peaks as valid only if they are included in `pk_inds_keep(:,1)`. This method is described in comments as better than the simple `hr_inst` summary, because it uses continuous valid-peak sections.

The exact implementation of `Cal_avg_HeartRate_Yue` is not in this repo, so its details should be treated as an external dependency. For reproducibility inside this repo, we may eventually want to replace or wrap it with an explicit local implementation.

### Good-Section Approach in `find_good_r_peak_sections.m`

`find_good_r_peak_sections.m` provides a local way to identify usable R-peak runs:

1. Compute RR intervals from `t_rpeaks`.
2. Reject intervals longer than 200 ms, corresponding to very low mouse HR.
3. Mark both peaks around a bad interval as bad.
4. Keep consecutive runs of surviving peaks.
5. Drop sections shorter than `min_duration`, default 5 seconds.

This helper is useful for building a cleaner instantaneous HR time series:

```matlab
t_rpeaks = t_ECG(good_r_peaks);
t_rpeaks = t_rpeaks(t_rpeaks > 0 & t_rpeaks < t_frames(end));
good_sections = find_good_r_peak_sections(t_rpeaks);
```

Then pool `good_sections(k).peak_time` and compute RR intervals or binned HR.

## Recommended Heart Rate Representations

For downstream correlation, keep these distinct:

- Beat-level HR: one value per RR interval, at time `t_rpeaks(2:end)` or RR midpoints.
- Binned HR: mean or median HR per analysis bin, aligned to diameter features.
- Continuous HR: beat-level HR interpolated onto `t_frames` or another common grid.

Recommended construction:

1. Align ECG time: `t_ECG = t_ECG - t_ecgoffset`.
2. Select good R-peaks using either `good_r_peaks` or `r_peak_confidence >= threshold`.
3. Clean into valid sections with `find_good_r_peak_sections`.
4. Compute RR intervals:

```matlab
rr = diff(t_good);
t_hr = t_good(2:end);
hr_hz = 1 ./ rr;
```

5. Reject physiologically implausible HR values, for example outside 150-1200 BPM.
6. For bin-level analyses, summarize `hr_hz` within the same bins as diameter features.
7. For continuous analyses, interpolate with care:

```matlab
hr_on_frames = interp1(t_hr, hr_hz, t_frames, 'pchip', NaN);
```

Avoid extrapolating over long gaps or across rejected R-peak sections.

## Downstream Correlation Examples

### `compute_hr_diam_corr.m`

This script loops over `data/*/results.mat`, loads binned results, and computes Spearman correlations:

```matlab
hr = results.avg_HR_Yue_10min;
amp = results.(band);
r = corr(hr', amp', 'Type', 'Spearman');
```

It splits bins into:

- early: `timepts < 60`
- late: `timepts >= 60`
- overall: all bins

Current caution: the script sets `band = 'hr'`, which means it correlates ECG HR against `results.hr`. In `results.mat`, `hr` is the vessel-diameter cardiac peak frequency, not vessel amplitude. That is meaningful only if the question is whether the vessel cardiac spectral peak tracks ECG HR. For vessel amplitude coupling, use one of:

```matlab
band = 'amp_svm';
band = 'amp_svm2';
band = 'bp_svm';
band = 'bp_svm2';
band = 'amp_c';
band = 'diam_tp';
```

Also consider using `corr(..., 'Rows', 'complete')` to avoid NaN problems.

### `compute_hr_diam_cross_corr.m`

This script computes normalized cross-correlation across binned time series:

## Current Group-Level Correlation Workflow

The project now has two parallel correlation-analysis folders:

- `pearson_corr/`
- `spearman_correlation/`

Both workflows use the same mouse-level inputs:

- `data/*/hr_amp_data.mat`

Both use the same HR feature:

- `hr_median_bpm`

Both use the same vessel features:

- `amp_svm`
- `bp_svm`
- `amp_c`
- `diam_mean`

Both use the same time windows:

- Early: `timepts < 60`
- Late: `timepts >= 60`
- Overall: all bins

Both use the same exclusion rule:

- exclude subjects with `dubious_hr_bin_fraction > 0.50`

At present, this excludes:

- `M167-dex40-04032025`

### Pearson Workflow

Primary script:

- `pearson_corr/compute_hr_diam_corr_pearson.m`

Outputs:

- `hr_amp_corr_pearson_results.xlsx`
- `hr_amp_corr_pearson_by_mouse.csv`
- `hr_amp_corr_pearson_wilcoxon_summary.csv`
- `hr_amp_corr_pearson_excluded_subjects.csv`

### Spearman Workflow

Primary script:

- `spearman_correlation/compute_hr_diam_corr_spearman.m`

Outputs:

- `hr_amp_corr_spearman_results.xlsx`
- `hr_amp_corr_spearman_by_mouse.csv`
- `hr_amp_corr_spearman_wilcoxon_summary.csv`
- `hr_amp_corr_spearman_excluded_subjects.csv`

### Important Unit-of-Inference Convention

The project now distinguishes clearly between:

- per-mouse correlation tests across 5-minute bins,
- group-level tests across mice.

Per-mouse correlation outputs:

- one correlation coefficient and one `p` value per mouse, vessel feature, and time window.

Group-level tests:

- use one mouse-level correlation per mouse,
- do **not** pool all 5-minute bins across subjects,
- treat mouse as the unit of replication.

This avoids pseudoreplication from combining repeated bins across animals.

## Group-Level Wilcoxon Testing

Both correlation folders include a Wilcoxon summary table.

This summary tests whether the distribution of mouse-level correlations is centered away from `0`.

Interpretation:

- a per-mouse correlation `p` value answers whether a single mouse shows a nonzero HR-vessel association across time bins;
- a Wilcoxon `p` value answers whether the correlation is consistently positive or negative across mice.

The Wilcoxon output is currently paired with:

- the median mouse-level correlation,
- a per-mouse point plot,
- a zero reference line.

## Fisher Z One-Sample T-Tests

For group-level parametric inference, both folders now include a Fisher `z` t-test layer.

Scripts:

- `pearson_corr/compute_hr_diam_corr_pearson_ttest.m`
- `spearman_correlation/compute_hr_diam_corr_spearman_ttest.m`

Method:

1. Read the existing per-mouse correlation coefficients.
2. Clip `r` values slightly away from `-1` and `1`.
3. Apply the Fisher transform:

```matlab
z = atanh(r);
```

4. Run a one-sample t-test of `z` against `0`.
5. Convert the mean and confidence interval back to correlation units:

```matlab
r_mean = tanh(mean_z);
```

The t-test outputs are written to:

- `pearson_corr/hr_amp_corr_pearson_ttest_summary.csv`
- `spearman_correlation/hr_amp_corr_spearman_ttest_summary.csv`

and added to the Excel workbooks as a `TTest` sheet.

### Confidence Interval Interpretation

The 95% confidence intervals in the t-test summaries:

- are computed in Fisher `z` space,
- then back-transformed to raw correlation units for display and reporting.

In the summary plots:

- solid red line: back-transformed mean group-level correlation,
- red dotted lines: back-transformed 95% CI bounds,
- black dashed line: zero-correlation reference.

If the CI includes `0`, the group-level effect is not significant at roughly the `0.05` level.

## Diagnostic Plots

The current workflow generates two main plot types in both folders:

- per-mouse scatter diagnostics,
- group-level summary plots.

Scatter diagnostics:

- show HR vs vessel feature for each mouse,
- color early and late bins separately,
- overlay a linear fit guide,
- annotate the per-mouse correlation result.

Group-level summary plots:

- show one point per mouse for each vessel feature and time window,
- summarize either Wilcoxon or Fisher `z` t-test results,
- support quick comparison of Pearson and Spearman inference.

## Current Interpretation

With the present dataset of 10 included mice, the current group-level analyses did not identify a statistically significant consistent HR-vessel association in any tested vessel feature or time window.

This was true for:

- Pearson plus Wilcoxon,
- Pearson plus Fisher `z` t-test,
- Spearman plus Wilcoxon,
- Spearman plus Fisher `z` t-test.

Some individual mice showed significant within-mouse correlations, but those effects were not consistent enough across mice to produce a significant group-level result.

```matlab
hr = results.hr;
amp = results.(band);
[xc, lags] = xcorr(hr - mean(hr), amp - mean(amp), max_lag, 'normalized');
```

Current caution: this uses `results.hr`, the vessel-spectrum peak frequency, not ECG-derived HR. For HR-vessel coupling, it should usually use:

```matlab
hr = results.avg_HR_Yue_10min;
```

The lag units are bins. With 10-minute bins, `max_lag = 3` means plus/minus 30 minutes.

## Higher-Resolution Analysis Direction

The current `results.mat` pipeline is bin-level and useful for coarse trends. For "instantaneous HR versus vessel amplitude change in different frequency bands", the stronger analysis path is:

1. Build a continuous, cleaned HR time series from R-peaks.
2. Filter diameter into the frequency band of interest.
3. Convert filtered diameter to an amplitude envelope if the question is about band amplitude:

```matlab
diam_band = bandpass(diam_norm, Fc, fs_cam);
diam_amp = abs(hilbert(diam_band));
```

4. Put HR and diameter amplitude on a common time grid.
5. Optionally smooth both signals at the timescale of interest.
6. Compute:
   - Spearman or Pearson correlation.
   - Lagged cross-correlation.
   - Magnitude-squared coherence.
   - Windowed correlations to test time-varying coupling.

Good common grids:

- `t_frames` for frame-level vessel analysis.
- A lower-rate grid, such as 1-10 Hz, for HR variability and slow-band analysis.
- Bin centers for coarse 10-minute summaries.

Avoid interpreting frame-level samples as independent observations after heavy filtering or interpolation. Use effective degrees of freedom, windowed statistics, permutation tests, or subject-level summaries.

## Frequency Band Interpretation

Useful band meanings in this project:

- Below about 0.3 Hz: slow vasomotion and slow vessel tone changes.
- 0.05-0.3 Hz: cleaner slow vasomotion band with less ultra-slow drift.
- About 3.5-15 Hz in dex mice: cardiac-band vessel pulsation.
- Above 15 Hz: high-frequency noise proxy.

Mouse HR is often around 5-15 Hz, depending on anesthesia and condition. For each subject, the cardiac band should be checked against the ECG-derived HR distribution. The fixed `Fc_c` values are reasonable defaults but may be too broad or too low/high for some subjects.

## Quality-Control Checks

Before trusting correlations:

1. Confirm `t_ECG - t_ecgoffset` aligns to `t_frames`.
2. Plot candidate and accepted R-peaks over ECG for a few windows.
3. Plot RR interval distributions and rejected sections.
4. Confirm the ECG-derived median HR lies inside the chosen cardiac diameter band.
5. Inspect `diam`, `diam_avg`, and filtered diameter traces for motion or focus artifacts.
6. Compare `amp_c` to `amp_noise`; cardiac amplitude near noise level should be treated cautiously.
7. Use `diam_norm` or subject-normalized features when comparing across animals.
8. Track NaNs and missing raw files subject by subject.

## Important Variable Name Cautions

- `results.hr`: vessel diameter PSD peak frequency in the cardiac band, in Hz. It is not ECG HR.
- `avg_HR_Yue_10min`: ECG-derived average HR per 10-minute bin.
- `med_HR`, `var_HR`: ECG-derived binned HR summaries from the external Yue/Andy function.
- `hr_inst`: beat-to-beat ECG HR in Hz inside `PlotResults_continuous_v5.m`, but not saved in the inspected `results.mat`.
- `amp_svm` and `amp_svm2`: IQR amplitude of slow diameter fluctuations, in microns.
- `bp_svm` and `bp_svm2`: PSD band power of slow diameter fluctuations.
- `amp_c`: IQR amplitude of cardiac-band vessel pulsation, in microns.

## Suggested Standard Pipeline for This Study

For a reproducible HR-diameter coupling analysis:

1. For each subject, load `signals.mat` and `im.mat`.
2. Align `t_ECG` to imaging time using `t_ecgoffset`.
3. Clean vessel diameter with spike interpolation and moving-average detrending.
4. Save both absolute micron change and normalized percent change.
5. Build cleaned R-peak sections from `good_r_peaks` or `r_peak_confidence`.
6. Compute beat-level HR and a continuous HR trace.
7. Define analysis bands:
   - slow: `[0.01 0.3]` or `[0.05 0.3]` Hz,
   - cardiac: subject-specific around ECG HR, or legacy `Fc_c`,
   - optional noise: `>15 Hz`.
8. Extract diameter-band traces and amplitude envelopes.
9. Aggregate HR and diameter features on matched time supports.
10. Run correlation, cross-correlation, coherence, and windowed analyses.
11. Store subject-level outputs in `results/` with unambiguous names, such as `ecg_hr_bpm`, `diam_peak_freq_hz`, `diam_svm_amp_um`, and `diam_cardiac_amp_um`.

## Near-Term Cleanup Recommendations

- Rename or document `results.hr` as `diam_cardiac_peak_freq_hz`.
- Update `compute_hr_diam_cross_corr.m` to use `avg_HR_Yue_10min` when the goal is ECG HR coupling.
- Change the default `band` in `compute_hr_diam_corr.m` from `'hr'` to a true vessel feature when studying diameter amplitude.
- Add `Rows`, `complete` to correlation calls.
- Save beat-level or continuous HR outputs in `results.mat` so downstream scripts do not depend only on 10-minute summaries.
- Bring the external `Cal_avg_HeartRate_Yue` dependency into the repo or replace it with a local function based on `find_good_r_peak_sections.m`.
- Prefer normalized diameter features for cross-subject analyses, while keeping micron units for within-subject physiology.
