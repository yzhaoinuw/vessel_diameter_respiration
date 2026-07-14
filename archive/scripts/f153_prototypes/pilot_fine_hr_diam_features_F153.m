%% Pilot fine-resolution HR and vessel diameter feature extraction for F153
% This script mirrors the core feature extraction in PlotResults_continuous_v5.m
% but uses finer bins and saves only the variables needed for HR-diameter
% correlation / cross-correlation studies.

clear; close all;

data_folder = fullfile(pwd, "data");
subject_name = "F153-Dex40-06022025";
subject_folder = fullfile(data_folder, subject_name);

out_file = fullfile(subject_folder, "results_fine.mat");
fig_folder = fullfile(pwd, "figures");
if ~exist(fig_folder, "dir")
    mkdir(fig_folder);
end

%% Parameters
conf_thresh = 0.5;
start_min = 14;
bin_minutes = 5;

% A 5-min window contains 1.5 cycles at 0.005 Hz, 3 cycles at 0.01 Hz,
% and 90 cycles at 0.3 Hz. This is a reasonable first pilot resolution for
% [0.01 0.3] Hz vasomotion summaries; shorter bins become unstable near
% the lower edge of the band.
segtime = bin_minutes;

Fc_svm = [0.01 0.3];
Fc_svm2 = [0.05 0.3];
Fc_c = [3.5 15];       % dex subject, following PlotResults_continuous_v5.m
Fc_noise = 15;

hr_valid_bpm = [250 750];       % typical broad anesthesia sanity range
min_beats_per_bin = 20;
ecg_plot_start_min = start_min;
ecg_plot_duration_sec = 3;

%% Load data
signal_data = load(fullfile(subject_folder, "signals.mat"), ...
    "ECG", "t_ECG", "detected_r_peaks", "r_peak_confidence", "good_r_peaks");
img_data = load(fullfile(subject_folder, "im.mat"), ...
    "diam", "t_frames", "t_ecgoffset", "umperpix");

ECG = double(signal_data.ECG(:));
t_ECG = double(signal_data.t_ECG(:)) - double(img_data.t_ecgoffset);
t_frames = double(img_data.t_frames(:));
diam = double(img_data.diam(:));
umperpix = double(img_data.umperpix);

fs_ECG = 1 / median(diff(t_ECG));
fs_median = 1 / median(diff(t_frames));

%% Confidence-thresholded R-peaks and instantaneous HR
detected_r_peaks = double(signal_data.detected_r_peaks(:));
r_peak_confidence = double(signal_data.r_peak_confidence(:));
good_r_peaks = double(signal_data.good_r_peaks(:));

is_confident = r_peak_confidence >= conf_thresh;
is_good_peak = ismember(detected_r_peaks, good_r_peaks);
r_peak_inds = detected_r_peaks(is_confident & is_good_peak);
r_peak_inds = unique(r_peak_inds(:));
r_peak_inds = r_peak_inds(r_peak_inds >= 1 & r_peak_inds <= numel(t_ECG));

t_rpeaks = t_ECG(r_peak_inds);
in_recording = t_rpeaks >= 0 & t_rpeaks <= t_frames(end);
t_rpeaks = t_rpeaks(in_recording);
r_peak_inds = r_peak_inds(in_recording);

rr = diff(t_rpeaks);
t_hr = t_rpeaks(2:end);
hr_inst_hz = 1 ./ rr;
hr_inst_bpm = 60 .* hr_inst_hz;

valid_inst_hr = isfinite(hr_inst_bpm) & ...
    hr_inst_bpm >= hr_valid_bpm(1) & hr_inst_bpm <= hr_valid_bpm(2);

%% Vessel diameter processing, following PlotResults_continuous_v5.m
diam_avg = smoothdata(diam, "movmean", 10 * 60 * fs_median);
diam_det = diam - diam_avg;

diam_filt_svm = lowpass(diam_det, Fc_svm(2), fs_median);
diam_filt_svm2 = bandpass(diam_det, Fc_svm2, fs_median);
diam_filt_c = bandpass(diam_det, Fc_c, fs_median);
diam_filt_noise = highpass(diam_det, Fc_noise, fs_median);

%% Fine bins
last_start_min = floor(t_frames(end) / 60 - segtime);
timepts = start_min:segtime:last_start_min;
n_bins = numel(timepts);

amp_svm = nan(1, n_bins);
amp_svm2 = nan(1, n_bins);
bp_svm = nan(1, n_bins);
bp_svm2 = nan(1, n_bins);
amp_c = nan(1, n_bins);
amp_noise = nan(1, n_bins);
diam_tp = nan(1, n_bins);

hr_mean_bpm_raw = nan(1, n_bins);
hr_median_bpm_raw = nan(1, n_bins);
hr_iqr_bpm_raw = nan(1, n_bins);
n_beats_bin = zeros(1, n_bins);
is_dubious_hr_bin = false(1, n_bins);

for i = 1:n_bins
    t0 = timepts(i) * 60;
    t1 = (timepts(i) + segtime) * 60;

    ind = find(t_frames > t0, 1):find(t_frames > t1, 1);
    if isempty(ind) || any(isnan(ind))
        continue
    end

    diam_tp(i) = mean(diam(ind), "omitnan") * umperpix;
    amp_svm(i) = iqr(diam_filt_svm(ind)) * umperpix;
    amp_svm2(i) = iqr(diam_filt_svm2(ind)) * umperpix;
    amp_c(i) = iqr(diam_filt_c(ind)) * umperpix;
    amp_noise(i) = iqr(diam_filt_noise(ind)) * umperpix;

    [pxx, f] = pwelch(diam_det(ind) * umperpix, [], [], [], fs_median);
    svm_mask = f > Fc_svm(1) & f < Fc_svm(2);
    svm2_mask = f > Fc_svm2(1) & f < Fc_svm2(2);
    bp_svm(i) = trapz(f(svm_mask), pxx(svm_mask));
    bp_svm2(i) = trapz(f(svm2_mask), pxx(svm2_mask));

    hr_in_bin = hr_inst_bpm(t_hr >= t0 & t_hr < t1 & valid_inst_hr);
    n_beats_bin(i) = numel(hr_in_bin);

    if n_beats_bin(i) >= min_beats_per_bin
        hr_mean_bpm_raw(i) = mean(hr_in_bin, "omitnan");
        hr_median_bpm_raw(i) = median(hr_in_bin, "omitnan");
        hr_iqr_bpm_raw(i) = iqr(hr_in_bin);
    end
end

is_dubious_hr_bin = n_beats_bin < min_beats_per_bin | ...
    ~isfinite(hr_median_bpm_raw) | ...
    hr_median_bpm_raw < hr_valid_bpm(1) | ...
    hr_median_bpm_raw > hr_valid_bpm(2);

hr_mean_bpm = hr_mean_bpm_raw;
hr_median_bpm = hr_median_bpm_raw;
hr_iqr_bpm = hr_iqr_bpm_raw;

reasonable_bins = ~is_dubious_hr_bin;
if nnz(reasonable_bins) >= 2
    hr_mean_bpm(is_dubious_hr_bin) = interp1(timepts(reasonable_bins), ...
        hr_mean_bpm_raw(reasonable_bins), timepts(is_dubious_hr_bin), ...
        "linear", "extrap");
    hr_median_bpm(is_dubious_hr_bin) = interp1(timepts(reasonable_bins), ...
        hr_median_bpm_raw(reasonable_bins), timepts(is_dubious_hr_bin), ...
        "linear", "extrap");
    hr_iqr_bpm(is_dubious_hr_bin) = interp1(timepts(reasonable_bins), ...
        hr_iqr_bpm_raw(reasonable_bins), timepts(is_dubious_hr_bin), ...
        "linear", "extrap");
else
    warning("Fewer than two reasonable HR bins; leaving dubious HR bins as NaN.");
end

%% Save compact fine-resolution result
save(out_file, ...
    "subject_name", "conf_thresh", "start_min", "bin_minutes", "timepts", ...
    "Fc_svm", "Fc_svm2", "Fc_c", "Fc_noise", ...
    "t_hr", "hr_inst_bpm", "valid_inst_hr", ...
    "hr_mean_bpm", "hr_median_bpm", "hr_iqr_bpm", ...
    "hr_mean_bpm_raw", "hr_median_bpm_raw", "hr_iqr_bpm_raw", ...
    "is_dubious_hr_bin", "n_beats_bin", "hr_valid_bpm", ...
    "amp_svm", "amp_svm2", "bp_svm", "bp_svm2", "amp_c", ...
    "amp_noise", "diam_tp", "fs_median", "fs_ECG");

%% Plot ECG subsection with accepted R-peaks and corresponding instantaneous HR
plot_t0 = ecg_plot_start_min * 60;
plot_t1 = plot_t0 + ecg_plot_duration_sec;
ecg_mask = t_ECG >= plot_t0 & t_ECG <= plot_t1;
peak_mask_plot = t_rpeaks >= plot_t0 & t_rpeaks <= plot_t1;
hr_mask_plot = t_hr >= plot_t0 & t_hr <= plot_t1 & valid_inst_hr;

fig = figure("Color", "w", "Position", [100 100 1100 430]);
yyaxis left
plot(t_ECG(ecg_mask) / 60, ECG(ecg_mask), "k-", "LineWidth", 0.9);
hold on;
peak_times_min = t_rpeaks(peak_mask_plot) / 60;
for k = 1:numel(peak_times_min)
    xline(peak_times_min(k), "r--", "LineWidth", 0.8);
end
ylabel("ECG");

yyaxis right
plot(t_hr(hr_mask_plot) / 60, hr_inst_bpm(hr_mask_plot), ...
    "Color", [0 0.4470 0.7410], "LineWidth", 1.2, "Marker", ".");
ylabel("Instantaneous HR (BPM)");
ylim(hr_valid_bpm);

xlabel("Time (min)");
title(sprintf("%s ECG R-peaks and instantaneous HR, %.3f-%.3f min", ...
    subject_name, plot_t0 / 60, plot_t1 / 60));
grid on;

ax = gca;
ax.Color = "w";
ax.XColor = "k";
ax.YAxis(1).Color = "k";
ax.YAxis(2).Color = [0 0.4470 0.7410];
ax.Toolbar.Visible = "off";
xlim([plot_t0 plot_t1] / 60);
set(findall(fig, "Type", "text"), "Color", "k");

fig_file = fullfile(fig_folder, subject_name + "_pilot_ecg_rpeaks_inst_hr.png");
exportgraphics(fig, fig_file, "Resolution", 200);

fprintf("Saved fine-resolution result: %s\n", out_file);
fprintf("Saved pilot plot: %s\n", fig_file);
fprintf("Reasonable HR bins: %d / %d\n", nnz(reasonable_bins), n_bins);
