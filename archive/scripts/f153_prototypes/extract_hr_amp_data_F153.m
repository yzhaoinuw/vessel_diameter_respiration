%% Extract compact HR and vessel amplitude data for F153
% Saves a fine-resolution MAT file for HR-diameter correlation studies.

clear;

data_folder = fullfile(pwd, "data");
subject_name = "F153-Dex40-06022025";
subject_folder = fullfile(data_folder, subject_name);
out_file = fullfile(subject_folder, "hr_amp_data.mat");

%% Parameters
conf_thresh = 0.5;
start_min = 14;
bin_minutes = 5;

Fc_svm = [0.01 0.3];
Fc_c = [3.5 15];       % dex subject, following PlotResults_continuous_v5.m

hr_valid_bpm = [250 750];
min_beats_per_bin = 20;

%% Load data
signal_data = load(fullfile(subject_folder, "signals.mat"), ...
    "t_ECG", "detected_r_peaks", "r_peak_confidence", "good_r_peaks");
img_data = load(fullfile(subject_folder, "im.mat"), ...
    "diam", "t_frames", "t_ecgoffset", "umperpix");

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
hr_inst_bpm = 60 ./ rr;
valid_inst_hr = isfinite(hr_inst_bpm) & ...
    hr_inst_bpm >= hr_valid_bpm(1) & hr_inst_bpm <= hr_valid_bpm(2);

%% Vessel diameter processing, matching the legacy choices where possible
diam_avg = smoothdata(diam, "movmean", 10 * 60 * fs_median);
diam_det = diam - diam_avg;

diam_filt_svm = lowpass(diam_det, Fc_svm(2), fs_median);
diam_filt_c = bandpass(diam_det, Fc_c, fs_median);

%% Fine bins
last_start_min = floor(t_frames(end) / 60 - bin_minutes);
timepts = start_min:bin_minutes:last_start_min;
n_bins = numel(timepts);
bin_centers_min = timepts + bin_minutes / 2;

amp_svm = nan(1, n_bins);
bp_svm = nan(1, n_bins);
amp_c = nan(1, n_bins);
diam_mean = nan(1, n_bins);

hr_mean_bpm_raw = nan(1, n_bins);
hr_median_bpm_raw = nan(1, n_bins);
hr_iqr_bpm_raw = nan(1, n_bins);
n_beats_bin = zeros(1, n_bins);

for i = 1:n_bins
    t0 = timepts(i) * 60;
    t1 = (timepts(i) + bin_minutes) * 60;

    idx0 = find(t_frames > t0, 1);
    idx1 = find(t_frames > t1, 1);
    if isempty(idx0) || isempty(idx1) || idx1 <= idx0
        continue
    end
    ind = idx0:idx1;

    diam_mean(i) = mean(diam(ind), "omitnan") * umperpix;
    amp_svm(i) = iqr(diam_filt_svm(ind)) * umperpix;
    amp_c(i) = iqr(diam_filt_c(ind)) * umperpix;

    [pxx, f] = pwelch(diam_det(ind) * umperpix, [], [], [], fs_median);
    svm_mask = f > Fc_svm(1) & f < Fc_svm(2);
    bp_svm(i) = trapz(f(svm_mask), pxx(svm_mask));

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
dubious_hr_bin_fraction = mean(is_dubious_hr_bin);

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

save(out_file, ...
    "subject_name", "conf_thresh", "start_min", "bin_minutes", ...
    "timepts", "bin_centers_min", ...
    "Fc_svm", "Fc_c", ...
    "t_rpeaks", "r_peak_inds", "t_hr", "hr_inst_bpm", "valid_inst_hr", ...
    "hr_mean_bpm", "hr_median_bpm", "hr_iqr_bpm", ...
    "hr_mean_bpm_raw", "hr_median_bpm_raw", "hr_iqr_bpm_raw", ...
    "is_dubious_hr_bin", "dubious_hr_bin_fraction", "n_beats_bin", ...
    "hr_valid_bpm", "min_beats_per_bin", ...
    "amp_svm", "bp_svm", "amp_c", "diam_mean", ...
    "fs_ECG", "fs_median");

fprintf("Saved %s\n", out_file);
fprintf("Dubious HR bins interpolated: %d / %d (%.3f)\n", ...
    nnz(is_dubious_hr_bin), n_bins, dubious_hr_bin_fraction);
