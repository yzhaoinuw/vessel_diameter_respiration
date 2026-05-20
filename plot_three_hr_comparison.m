%% Plot three heart-rate estimates for each subject
% Compares confidence-filtered 1-minute ECG HR, Yue 10-minute ECG HR, and
% vessel diameter-derived HR from results.mat.

clear;

data_folder = fullfile(pwd, "data");
fig_folder = fullfile(pwd, "figures", "three_hr_comparison");
if ~exist(fig_folder, "dir")
    mkdir(fig_folder);
end

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, "."));

summary = table();

for s = 1:numel(subjects)
    subject_name = string(subjects(s).name);
    subject_folder = fullfile(data_folder, subject_name);
    hr_file = fullfile(subject_folder, "hr_amp_data.mat");
    results_file = fullfile(subject_folder, "results.mat");

    if exist(hr_file, "file") ~= 2 || exist(results_file, "file") ~= 2
        fprintf("Skipping %s: hr_amp_data.mat=%d, results.mat=%d\n", ...
            subject_name, exist(hr_file, "file") == 2, ...
            exist(results_file, "file") == 2);
        continue
    end

    hr_data = load(hr_file, "hr_filter_bin_centers_min", ...
        "hr_minute_bpm_filtered", "is_bad_conf_hr_filter_bin", ...
        "hr_filter_low_conf_fraction");
    result_data = load(results_file, "timepts", "avg_HR_Yue_10min", "hr");

    t_ecg_min = double(hr_data.hr_filter_bin_centers_min(:));
    filtered_ecg_bpm = double(hr_data.hr_minute_bpm_filtered(:));

    t_results_min = double(result_data.timepts(:));
    yue_bpm = double(result_data.avg_HR_Yue_10min(:));
    vessel_hr_bpm = double(result_data.hr(:)) * 60;

    fig = figure("Visible", "off", "Color", "w", ...
        "Position", [100 100 1200 650]);
    ax = axes(fig);
    set(ax, "Color", "w", "XColor", "k", "YColor", "k", ...
        "GridColor", [0.75 0.75 0.75]);
    hold(ax, "on");

    plot(ax, t_ecg_min, filtered_ecg_bpm, "o-", ...
        "Color", [0.00 0.45 0.74], ...
        "MarkerFaceColor", [0.00 0.45 0.74], ...
        "MarkerSize", 4, "LineWidth", 1.2, ...
        "DisplayName", "Filtered ECG HR, 1-min bins");
    plot(ax, t_results_min, yue_bpm, "^-", ...
        "Color", [0.85 0.33 0.10], ...
        "MarkerFaceColor", [0.85 0.33 0.10], ...
        "MarkerSize", 7, "LineWidth", 1.5, ...
        "DisplayName", "avg_HR_Yue_10min");
    plot(ax, t_results_min, vessel_hr_bpm, "s-", ...
        "Color", [0.47 0.67 0.19], ...
        "MarkerFaceColor", [0.47 0.67 0.19], ...
        "MarkerSize", 7, "LineWidth", 1.5, ...
        "DisplayName", "results.hr from vessel diam x 60");

    grid(ax, "on");
    box(ax, "off");
    xlabel(ax, "Time (min)", "Color", "k");
    ylabel(ax, "Heart rate (BPM)", "Color", "k");
    title(ax, sprintf("%s: three HR estimates", subject_name), ...
        "Interpreter", "none", "Color", "k");
    lgd = legend(ax, "Location", "best", "Interpreter", "none");
    set(lgd, "Color", "w", "TextColor", "k", "EdgeColor", [0.4 0.4 0.4]);

    all_t = [t_ecg_min; t_results_min];
    all_y = [filtered_ecg_bpm; yue_bpm; vessel_hr_bpm];
    finite_t = isfinite(all_t);
    finite_y = isfinite(all_y);
    if any(finite_t)
        xlim(ax, [min(all_t(finite_t)) max(all_t(finite_t))]);
    end
    if any(finite_y)
        y_range = [min(all_y(finite_y)) max(all_y(finite_y))];
        pad = max(10, 0.05 * diff(y_range));
        ylim(ax, y_range + [-pad pad]);
    end

    out_png = fullfile(fig_folder, subject_name + "_three_hr_comparison.png");
    exportgraphics(fig, out_png, "Resolution", 300);
    close(fig);

    n_filtered_na = nnz(~isfinite(filtered_ecg_bpm));
    n_filtered_total = numel(filtered_ecg_bpm);
    n_bad_conf = nnz(hr_data.is_bad_conf_hr_filter_bin);
    summary = [summary; table(subject_name, n_filtered_na, ...
        n_filtered_total, n_bad_conf, string(out_png), ...
        'VariableNames', {'Subject', 'NFilteredMinuteNa', ...
        'NFilteredMinuteBins', 'NBadConfidenceMinuteBins', 'FigureFile'})];

    fprintf("Saved %s | filtered ECG NaN bins=%d/%d\n", ...
        out_png, n_filtered_na, n_filtered_total);
end

summary_file = fullfile(fig_folder, "three_hr_comparison_summary.csv");
writetable(summary, summary_file);
fprintf("Wrote summary: %s\n", summary_file);
