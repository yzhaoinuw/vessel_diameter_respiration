%% Extract compact HR and vessel amplitude data for all subjects with sliding bins
% Saves one hr_amp_data_sliding.mat file in each subject folder.

clear;

data_folder = fullfile(pwd, "data");
summary_file = fullfile(pwd, "results", "hr_amp_extraction_summary_sliding.csv");
if ~exist(fullfile(pwd, "results"), "dir")
    mkdir(fullfile(pwd, "results"));
end

%% Parameters
conf_thresh = 0.5;
start_min = 14;
bin_minutes = 5;
bin_step_minutes = 1;
hr_filter_bin_minutes = 1;
low_conf_rpeak_fraction_threshold = 0.50;
min_filtered_hr_samples_per_bin = 1;

Fc_svm = [0.01 0.3];
hr_valid_bpm = [250 750];
min_beats_per_bin = 20;

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, "."));

summary = table();

for s = 1:numel(subjects)
    subject_name = string(subjects(s).name);
    subject_folder = fullfile(data_folder, subject_name);
    signals_file = fullfile(subject_folder, "signals.mat");
    im_file = fullfile(subject_folder, "im.mat");
    out_file = fullfile(subject_folder, "hr_amp_data_sliding.mat");

    has_signals = exist(signals_file, "file") == 2;
    has_im = exist(im_file, "file") == 2;

    if ~has_signals || ~has_im
        fprintf("Skipping %s: signals.mat=%d, im.mat=%d\n", ...
            subject_name, has_signals, has_im);
        summary = [summary; make_summary_row(subject_name, false, ...
            has_signals, has_im, NaN, NaN, NaN, "missing raw file")];
        continue
    end

    try
        result = extract_one_subject(subject_name, subject_folder, ...
            conf_thresh, start_min, bin_minutes, bin_step_minutes, ...
            Fc_svm, hr_valid_bpm, min_beats_per_bin, ...
            hr_filter_bin_minutes, low_conf_rpeak_fraction_threshold, ...
            min_filtered_hr_samples_per_bin);

        save(out_file, "-struct", "result");

        fprintf("Saved %s | bins=%d | dubious=%d/%d (%.3f)\n", ...
            out_file, result.n_bins, nnz(result.is_dubious_hr_bin), ...
            result.n_bins, result.dubious_hr_bin_fraction);

        summary = [summary; make_summary_row(subject_name, true, ...
            has_signals, has_im, result.n_bins, ...
            nnz(result.is_dubious_hr_bin), ...
            result.dubious_hr_bin_fraction, "")];
    catch ME
        warning("Failed %s: %s", subject_name, ME.message);
        summary = [summary; make_summary_row(subject_name, false, ...
            has_signals, has_im, NaN, NaN, NaN, string(ME.message))];
    end
end

writetable(summary, summary_file);
fprintf("Wrote extraction summary: %s\n", summary_file);

%% Local functions
function result = extract_one_subject(subject_name, subject_folder, ...
    conf_thresh, start_min, bin_minutes, bin_step_minutes, Fc_svm, ...
    hr_valid_bpm, min_beats_per_bin, hr_filter_bin_minutes, ...
    low_conf_rpeak_fraction_threshold, min_filtered_hr_samples_per_bin)

    signals_file = fullfile(subject_folder, "signals.mat");
    im_file = fullfile(subject_folder, "im.mat");

    signal_data = load(signals_file, ...
        "t_ECG", "detected_r_peaks", "r_peak_confidence", "good_r_peaks");

    im_vars = string({whos("-file", im_file).name});
    requested_im_vars = ["diam", "diam_dlc", "t_frames", "t_vect", ...
        "t_ecgoffset", "umperpix"];
    vars_to_load = cellstr(requested_im_vars(ismember(requested_im_vars, im_vars)));
    img_data = load(im_file, vars_to_load{:});

    t_ECG = double(signal_data.t_ECG(:)) - double(img_data.t_ecgoffset);
    if isfield(img_data, "t_frames")
        t_frames = double(img_data.t_frames(:));
        frame_time_source = "t_frames";
    elseif isfield(img_data, "t_vect")
        t_frames = double(img_data.t_vect(:));
        frame_time_source = "t_vect";
    else
        error("Neither t_frames nor t_vect was found in im.mat.");
    end
    if isfield(img_data, "diam")
        diam = double(img_data.diam(:));
        diam_source = "diam";
    elseif isfield(img_data, "diam_dlc")
        diam = double(img_data.diam_dlc(:));
        diam_source = "diam_dlc";
    else
        error("Neither diam nor diam_dlc was found in im.mat.");
    end
    umperpix = double(img_data.umperpix);

    fs_ECG = 1 / median(diff(t_ECG));
    fs_median = 1 / median(diff(t_frames));
    Fc_c = cardiac_band_for_subject(subject_name);

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

    [hr_filter_timepts_min, hr_filter_bin_centers_min, ...
        hr_filter_r_peak_count, hr_filter_low_conf_fraction, ...
        is_bad_conf_hr_filter_bin, hr_minute_bpm_raw, ...
        hr_minute_bpm_filtered, t_hr_filtered, hr_inst_bpm_filtered, ...
        valid_inst_hr_filtered] = filter_hr_by_minute_confidence( ...
        detected_r_peaks, r_peak_confidence, t_ECG, t_frames(end), ...
        t_hr, hr_inst_bpm, valid_inst_hr, conf_thresh, ...
        low_conf_rpeak_fraction_threshold, hr_filter_bin_minutes, ...
        hr_valid_bpm);

    diam_avg = smoothdata(diam, "movmean", 10 * 60 * fs_median);
    diam_det = diam - diam_avg;

    diam_filt_svm = lowpass(diam_det, Fc_svm(2), fs_median);
    diam_filt_c = bandpass(diam_det, Fc_c, fs_median);

    last_start_min = floor(t_frames(end) / 60 - bin_minutes);
    timepts = start_min:bin_step_minutes:last_start_min;
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
    hr_mean_bpm_filtered_raw = nan(1, n_bins);
    hr_median_bpm_filtered_raw = nan(1, n_bins);
    hr_iqr_bpm_filtered_raw = nan(1, n_bins);
    n_filtered_hr_bin = zeros(1, n_bins);

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

        hr_filtered_in_bin = hr_inst_bpm_filtered( ...
            t_hr_filtered >= t0 & t_hr_filtered < t1 & valid_inst_hr_filtered);
        n_filtered_hr_bin(i) = numel(hr_filtered_in_bin);

        if n_filtered_hr_bin(i) >= min_filtered_hr_samples_per_bin
            hr_mean_bpm_filtered_raw(i) = mean(hr_filtered_in_bin, "omitnan");
            hr_median_bpm_filtered_raw(i) = median(hr_filtered_in_bin, "omitnan");
            hr_iqr_bpm_filtered_raw(i) = iqr(hr_filtered_in_bin);
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
        warning("%s: fewer than two reasonable HR bins; leaving dubious HR bins as NaN.", ...
            subject_name);
    end

    is_dubious_hr_filtered_bin = n_filtered_hr_bin < min_filtered_hr_samples_per_bin | ...
        ~isfinite(hr_median_bpm_filtered_raw) | ...
        hr_median_bpm_filtered_raw < hr_valid_bpm(1) | ...
        hr_median_bpm_filtered_raw > hr_valid_bpm(2);
    dubious_hr_filtered_bin_fraction = mean(is_dubious_hr_filtered_bin);

    hr_mean_bpm_filtered = hr_mean_bpm_filtered_raw;
    hr_median_bpm_filtered = hr_median_bpm_filtered_raw;
    hr_iqr_bpm_filtered = hr_iqr_bpm_filtered_raw;

    reasonable_filtered_bins = ~is_dubious_hr_filtered_bin;
    if nnz(reasonable_filtered_bins) >= 2
        hr_mean_bpm_filtered(is_dubious_hr_filtered_bin) = interp1( ...
            timepts(reasonable_filtered_bins), ...
            hr_mean_bpm_filtered_raw(reasonable_filtered_bins), ...
            timepts(is_dubious_hr_filtered_bin), "linear", "extrap");
        hr_median_bpm_filtered(is_dubious_hr_filtered_bin) = interp1( ...
            timepts(reasonable_filtered_bins), ...
            hr_median_bpm_filtered_raw(reasonable_filtered_bins), ...
            timepts(is_dubious_hr_filtered_bin), "linear", "extrap");
        hr_iqr_bpm_filtered(is_dubious_hr_filtered_bin) = interp1( ...
            timepts(reasonable_filtered_bins), ...
            hr_iqr_bpm_filtered_raw(reasonable_filtered_bins), ...
            timepts(is_dubious_hr_filtered_bin), "linear", "extrap");
    elseif nnz(reasonable_filtered_bins) == 1
        hr_mean_bpm_filtered(is_dubious_hr_filtered_bin) = ...
            hr_mean_bpm_filtered_raw(reasonable_filtered_bins);
        hr_median_bpm_filtered(is_dubious_hr_filtered_bin) = ...
            hr_median_bpm_filtered_raw(reasonable_filtered_bins);
        hr_iqr_bpm_filtered(is_dubious_hr_filtered_bin) = ...
            hr_iqr_bpm_filtered_raw(reasonable_filtered_bins);
    else
        warning("%s: no reasonable filtered HR bins; leaving filtered HR bins as NaN.", ...
            subject_name);
    end

    result = struct();
    result.subject_name = subject_name;
    result.conf_thresh = conf_thresh;
    result.start_min = start_min;
    result.bin_minutes = bin_minutes;
    result.bin_step_minutes = bin_step_minutes;
    result.is_sliding_bin = true;
    result.timepts = timepts;
    result.bin_centers_min = bin_centers_min;
    result.n_bins = n_bins;
    result.Fc_svm = Fc_svm;
    result.Fc_c = Fc_c;
    result.frame_time_source = frame_time_source;
    result.diam_source = diam_source;
    result.t_rpeaks = t_rpeaks;
    result.r_peak_inds = r_peak_inds;
    result.t_hr = t_hr;
    result.hr_inst_bpm = hr_inst_bpm;
    result.valid_inst_hr = valid_inst_hr;
    result.hr_filter_bin_minutes = hr_filter_bin_minutes;
    result.low_conf_rpeak_fraction_threshold = low_conf_rpeak_fraction_threshold;
    result.hr_filter_timepts_min = hr_filter_timepts_min;
    result.hr_filter_bin_centers_min = hr_filter_bin_centers_min;
    result.hr_filter_r_peak_count = hr_filter_r_peak_count;
    result.hr_filter_low_conf_fraction = hr_filter_low_conf_fraction;
    result.is_bad_conf_hr_filter_bin = is_bad_conf_hr_filter_bin;
    result.hr_minute_bpm_raw = hr_minute_bpm_raw;
    result.hr_minute_bpm_filtered = hr_minute_bpm_filtered;
    result.t_hr_filtered = t_hr_filtered;
    result.hr_inst_bpm_filtered = hr_inst_bpm_filtered;
    result.valid_inst_hr_filtered = valid_inst_hr_filtered;
    result.hr_mean_bpm = hr_mean_bpm;
    result.hr_median_bpm = hr_median_bpm;
    result.hr_iqr_bpm = hr_iqr_bpm;
    result.hr_mean_bpm_filtered = hr_mean_bpm_filtered;
    result.hr_median_bpm_filtered = hr_median_bpm_filtered;
    result.hr_iqr_bpm_filtered = hr_iqr_bpm_filtered;
    result.hr_mean_bpm_raw = hr_mean_bpm_raw;
    result.hr_median_bpm_raw = hr_median_bpm_raw;
    result.hr_iqr_bpm_raw = hr_iqr_bpm_raw;
    result.hr_mean_bpm_filtered_raw = hr_mean_bpm_filtered_raw;
    result.hr_median_bpm_filtered_raw = hr_median_bpm_filtered_raw;
    result.hr_iqr_bpm_filtered_raw = hr_iqr_bpm_filtered_raw;
    result.is_dubious_hr_bin = is_dubious_hr_bin;
    result.dubious_hr_bin_fraction = dubious_hr_bin_fraction;
    result.is_dubious_hr_filtered_bin = is_dubious_hr_filtered_bin;
    result.dubious_hr_filtered_bin_fraction = dubious_hr_filtered_bin_fraction;
    result.n_beats_bin = n_beats_bin;
    result.n_filtered_hr_bin = n_filtered_hr_bin;
    result.hr_valid_bpm = hr_valid_bpm;
    result.min_beats_per_bin = min_beats_per_bin;
    result.min_filtered_hr_samples_per_bin = min_filtered_hr_samples_per_bin;
    result.amp_svm = amp_svm;
    result.bp_svm = bp_svm;
    result.amp_c = amp_c;
    result.diam_mean = diam_mean;
    result.fs_ECG = fs_ECG;
    result.fs_median = fs_median;
end

function [timepts_min, bin_centers_min, r_peak_count, low_conf_fraction, ...
    is_bad_conf_bin, hr_minute_raw, hr_minute_filtered, t_hr_filtered, ...
    hr_inst_bpm_filtered, valid_inst_hr_filtered] = ...
    filter_hr_by_minute_confidence(detected_r_peaks, r_peak_confidence, ...
    t_ECG, recording_end_sec, t_hr, hr_inst_bpm, valid_inst_hr, ...
    conf_thresh, low_conf_rpeak_fraction_threshold, bin_minutes, hr_valid_bpm)

    last_start_min = floor(recording_end_sec / 60 - bin_minutes);
    if last_start_min < 0
        timepts_min = [];
    else
        timepts_min = 0:bin_minutes:last_start_min;
    end

    n_bins = numel(timepts_min);
    bin_centers_min = timepts_min + bin_minutes / 2;
    r_peak_count = zeros(1, n_bins);
    low_conf_fraction = nan(1, n_bins);
    is_bad_conf_bin = true(1, n_bins);
    hr_minute_raw = nan(1, n_bins);

    peak_inds = round(double(detected_r_peaks(:)));
    peak_conf = double(r_peak_confidence(:));
    valid_peak = isfinite(peak_inds) & peak_inds >= 1 & peak_inds <= numel(t_ECG);
    peak_inds = peak_inds(valid_peak);
    peak_conf = peak_conf(valid_peak);
    peak_times = t_ECG(peak_inds);
    in_recording = peak_times >= 0 & peak_times <= recording_end_sec;
    peak_times = peak_times(in_recording);
    peak_conf = peak_conf(in_recording);

    for i = 1:n_bins
        t0 = timepts_min(i) * 60;
        t1 = (timepts_min(i) + bin_minutes) * 60;

        peaks_in_bin = peak_times >= t0 & peak_times < t1;
        r_peak_count(i) = nnz(peaks_in_bin);
        if r_peak_count(i) > 0
            low_conf_fraction(i) = mean(peak_conf(peaks_in_bin) < conf_thresh);
            is_bad_conf_bin(i) = ...
                low_conf_fraction(i) >= low_conf_rpeak_fraction_threshold;
        end

        if ~is_bad_conf_bin(i)
            hr_in_bin = hr_inst_bpm(t_hr >= t0 & t_hr < t1 & valid_inst_hr);
            if ~isempty(hr_in_bin)
                hr_minute_raw(i) = median(hr_in_bin, "omitnan");
            end
        end
    end

    hr_minute_filtered = hr_minute_raw;
    finite_bins = isfinite(hr_minute_raw);
    missing_bins = ~finite_bins;
    if nnz(finite_bins) >= 2
        hr_minute_filtered(missing_bins) = interp1( ...
            bin_centers_min(finite_bins), hr_minute_raw(finite_bins), ...
            bin_centers_min(missing_bins), "linear", "extrap");
    elseif nnz(finite_bins) == 1
        hr_minute_filtered(missing_bins) = hr_minute_raw(finite_bins);
    end

    t_hr_filtered = bin_centers_min * 60;
    hr_inst_bpm_filtered = hr_minute_filtered;
    valid_inst_hr_filtered = isfinite(hr_inst_bpm_filtered) & ...
        hr_inst_bpm_filtered >= hr_valid_bpm(1) & ...
        hr_inst_bpm_filtered <= hr_valid_bpm(2);
end

function Fc_c = cardiac_band_for_subject(subject_name)
    if contains(subject_name, "KX", "IgnoreCase", true)
        Fc_c = [1.5 10];
    elseif contains(subject_name, "dex", "IgnoreCase", true)
        Fc_c = [3.5 15];
    elseif contains(subject_name, "PPA", "IgnoreCase", true) || ...
            contains(subject_name, "PINP", "IgnoreCase", true)
        Fc_c = [1.5 15];
    else
        Fc_c = [6 15];
    end
end

function row = make_summary_row(subject_name, success, has_signals, has_im, ...
    n_bins, n_dubious_bins, dubious_fraction, message)

    row = table(subject_name, success, has_signals, has_im, n_bins, ...
        n_dubious_bins, dubious_fraction, message, ...
        'VariableNames', {'Subject', 'Success', 'HasSignals', 'HasIm', ...
        'NBins', 'NDubiousHrBins', 'DubiousHrBinFraction', 'Message'});
end
