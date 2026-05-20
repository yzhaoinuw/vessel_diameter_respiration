%% Extract overlapping 10-minute vessel features at 1-minute resolution
% Saves one vessel_features_10min_overlap.mat file in each subject folder.

clear;

data_folder = fullfile(pwd, "data");
summary_file = fullfile(pwd, "results", "vessel_features_10min_overlap_summary.csv");
if ~exist(fullfile(pwd, "results"), "dir")
    mkdir(fullfile(pwd, "results"));
end

%% Parameters
bin_minutes = 10;
step_minutes = 1;
Fc_svm = [0.01 0.3];

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, "."));

summary = table();

for s = 1:numel(subjects)
    subject_name = string(subjects(s).name);
    subject_folder = fullfile(data_folder, subject_name);
    im_file = fullfile(subject_folder, "im.mat");
    out_file = fullfile(subject_folder, "vessel_features_10min_overlap.mat");

    has_im = exist(im_file, "file") == 2;
    if ~has_im
        fprintf("Skipping %s: missing im.mat\n", subject_name);
        summary = [summary; make_summary_row(subject_name, false, ...
            has_im, NaN, NaN, "missing im.mat")];
        continue
    end

    try
        result = extract_one_subject(subject_name, subject_folder, ...
            bin_minutes, step_minutes, Fc_svm);

        save(out_file, "-struct", "result");

        fprintf("Saved %s | bins=%d | %.1f-%.1f min centers\n", ...
            out_file, result.n_bins, result.bin_centers_min(1), ...
            result.bin_centers_min(end));

        summary = [summary; make_summary_row(subject_name, true, ...
            has_im, result.n_bins, result.n_valid_bins, "")];
    catch ME
        warning("Failed %s: %s", subject_name, ME.message);
        summary = [summary; make_summary_row(subject_name, false, ...
            has_im, NaN, NaN, string(ME.message))];
    end
end

writetable(summary, summary_file);
fprintf("Wrote extraction summary: %s\n", summary_file);

%% Local functions
function result = extract_one_subject(subject_name, subject_folder, ...
    bin_minutes, step_minutes, Fc_svm)

    im_file = fullfile(subject_folder, "im.mat");
    im_vars = string({whos("-file", im_file).name});
    requested_im_vars = ["diam", "diam_dlc", "t_frames", "t_vect", "umperpix"];
    vars_to_load = cellstr(requested_im_vars(ismember(requested_im_vars, im_vars)));
    img_data = load(im_file, vars_to_load{:});

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
    fs_median = 1 / median(diff(t_frames), "omitnan");
    Fc_c = cardiac_band_for_subject(subject_name);

    diam_avg = smoothdata(diam, "movmean", 10 * 60 * fs_median);
    diam_det = diam - diam_avg;
    diam_filt_svm = lowpass(diam_det, Fc_svm(2), fs_median);
    diam_filt_c = bandpass(diam_det, Fc_c, fs_median);

    last_start_min = floor(t_frames(end) / 60 - bin_minutes);
    if last_start_min < 0
        timepts = [];
    else
        timepts = 0:step_minutes:last_start_min;
    end

    n_bins = numel(timepts);
    bin_centers_min = timepts + bin_minutes / 2;

    amp_svm = nan(1, n_bins);
    bp_svm = nan(1, n_bins);
    amp_c = nan(1, n_bins);
    n_frames_bin = zeros(1, n_bins);

    for i = 1:n_bins
        t0 = timepts(i) * 60;
        t1 = (timepts(i) + bin_minutes) * 60;

        idx0 = find(t_frames > t0, 1);
        idx1 = find(t_frames > t1, 1);
        if isempty(idx0) || isempty(idx1) || idx1 <= idx0
            continue
        end

        ind = idx0:idx1;
        n_frames_bin(i) = numel(ind);

        amp_svm(i) = iqr(diam_filt_svm(ind)) * umperpix;
        amp_c(i) = iqr(diam_filt_c(ind)) * umperpix;

        [pxx, f] = pwelch(diam_det(ind) * umperpix, [], [], [], fs_median);
        svm_mask = f > Fc_svm(1) & f < Fc_svm(2);
        bp_svm(i) = trapz(f(svm_mask), pxx(svm_mask));
    end

    result = struct();
    result.subject_name = subject_name;
    result.bin_minutes = bin_minutes;
    result.step_minutes = step_minutes;
    result.timepts = timepts;
    result.bin_centers_min = bin_centers_min;
    result.n_bins = n_bins;
    result.n_valid_bins = nnz(isfinite(amp_svm) & isfinite(bp_svm) & isfinite(amp_c));
    result.Fc_svm = Fc_svm;
    result.Fc_c = Fc_c;
    result.frame_time_source = frame_time_source;
    result.diam_source = diam_source;
    result.fs_median = fs_median;
    result.amp_svm = amp_svm;
    result.bp_svm = bp_svm;
    result.amp_c = amp_c;
    result.n_frames_bin = n_frames_bin;
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

function row = make_summary_row(subject_name, success, has_im, ...
    n_bins, n_valid_bins, message)

    row = table(subject_name, success, has_im, n_bins, n_valid_bins, message, ...
        'VariableNames', {'Subject', 'Success', 'HasIm', ...
        'NBins', 'NValidBins', 'Message'});
end
