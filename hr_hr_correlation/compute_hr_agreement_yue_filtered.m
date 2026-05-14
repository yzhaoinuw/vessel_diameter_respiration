%% Compare confidence-filtered ECG HR against legacy avg_HR_Yue_10min
% Rebins filtered 1-minute ECG HR from hr_amp_data.mat onto the legacy
% 10-minute windows used in results.mat, then computes per-mouse Pearson
% and Spearman correlations against results.avg_HR_Yue_10min. Per-mouse
% figures show the filtered-vs-Yue scatter plus filtered ECG HR,
% avg_HR_Yue_10min, and results.hr converted to BPM.

clear; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));
data_folder = fullfile(project_root, "data");
output_folder = fileparts(mfilename('fullpath'));
fig_folder = fullfile(output_folder, "figures");
mouse_fig_folder = fullfile(fig_folder, "filtered_yue_by_mouse");

if ~exist(output_folder, "dir")
    mkdir(output_folder);
end
if ~exist(fig_folder, "dir")
    mkdir(fig_folder);
end
if ~exist(mouse_fig_folder, "dir")
    mkdir(mouse_fig_folder);
end

cutoff_min = 60;
dubious_exclusion_threshold = 0.50;
legacy_bin_minutes = 10;
min_points_for_corr = 3;

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, "."));

excluded_table = table();
by_mouse_table = table();

for i = 1:numel(subjects)
    subject_name = string(subjects(i).name);
    hr_amp_file = fullfile(data_folder, subject_name, "hr_amp_data.mat");
    results_file = fullfile(data_folder, subject_name, "results.mat");

    if ~exist(hr_amp_file, "file") || ~exist(results_file, "file")
        missing_bits = strings(0, 1);
        if ~exist(hr_amp_file, "file")
            missing_bits(end + 1) = "hr_amp_data.mat"; %#ok<AGROW>
        end
        if ~exist(results_file, "file")
            missing_bits(end + 1) = "results.mat"; %#ok<AGROW>
        end
        fprintf("Skipping %s - missing %s\n", subject_name, strjoin(missing_bits, ", "));
        continue
    end

    hr_vars = string({whos("-file", hr_amp_file).name});
    required_hr_vars = ["t_hr_filtered", "hr_inst_bpm_filtered", ...
        "valid_inst_hr_filtered", "dubious_hr_filtered_bin_fraction"];
    missing_hr_vars = required_hr_vars(~ismember(required_hr_vars, hr_vars));
    if ~isempty(missing_hr_vars)
        fprintf("Skipping %s - missing filtered HR variables: %s\n", ...
            subject_name, strjoin(missing_hr_vars, ", "));
        excluded_table = [excluded_table; table(subject_name, NaN, ...
            "missing filtered HR variables", ...
            'VariableNames', {'Subject', 'DubiousHrBinFraction', 'Reason'})]; %#ok<AGROW>
        continue
    end

    hr_data = load(hr_amp_file, required_hr_vars{:});
    results_data = load(results_file, "timepts", "avg_HR_Yue_10min", "hr");

    dubious_fraction = hr_data.dubious_hr_filtered_bin_fraction;
    if dubious_fraction > dubious_exclusion_threshold
        fprintf("Excluding %s - dubious filtered HR bins %.1f%% > %.1f%%\n", ...
            subject_name, 100 * dubious_fraction, ...
            100 * dubious_exclusion_threshold);
        excluded_table = [excluded_table; table(subject_name, dubious_fraction, ...
            "dubious filtered bins exceed threshold", ...
            'VariableNames', {'Subject', 'DubiousHrBinFraction', 'Reason'})]; %#ok<AGROW>
        continue
    end

    timepts_min = double(results_data.timepts(:));
    yue_hr_bpm = double(results_data.avg_HR_Yue_10min(:));
    vessel_hr_bpm = 60 * double(results_data.hr(:));
    filtered_hr_bpm = rebin_ecg_hr_to_legacy_bins( ...
        double(hr_data.t_hr_filtered(:)), ...
        double(hr_data.hr_inst_bpm_filtered(:)), ...
        logical(hr_data.valid_inst_hr_filtered(:)), ...
        timepts_min, legacy_bin_minutes);

    n_common = min([numel(timepts_min), numel(yue_hr_bpm), ...
        numel(vessel_hr_bpm), numel(filtered_hr_bpm)]);
    timepts_min = timepts_min(1:n_common);
    yue_hr_bpm = yue_hr_bpm(1:n_common);
    vessel_hr_bpm = vessel_hr_bpm(1:n_common);
    filtered_hr_bpm = filtered_hr_bpm(1:n_common);

    row = table(subject_name, dubious_fraction, ...
        'VariableNames', {'Subject', 'DubiousHrBinFraction'});

    valid_overall = isfinite(filtered_hr_bpm) & isfinite(yue_hr_bpm);
    row.N_Overlap = nnz(valid_overall);
    row.MeanDiffBPM = mean(filtered_hr_bpm(valid_overall) - ...
        yue_hr_bpm(valid_overall), "omitnan");
    row.MeanAbsDiffBPM = mean(abs(filtered_hr_bpm(valid_overall) - ...
        yue_hr_bpm(valid_overall)), "omitnan");
    row.RMSE_BPM = sqrt(mean((filtered_hr_bpm(valid_overall) - ...
        yue_hr_bpm(valid_overall)).^2, "omitnan"));

    window_defs = {
        "Early",   timepts_min < cutoff_min;
        "Late",    timepts_min >= cutoff_min;
        "Overall", true(size(timepts_min))
    };

    for w = 1:size(window_defs, 1)
        window_name = window_defs{w, 1};
        mask = window_defs{w, 2};

        [r_spear, p_spear, n_spear] = corr_complete( ...
            filtered_hr_bpm(mask), yue_hr_bpm(mask), ...
            "Spearman", min_points_for_corr);
        [r_pear, p_pear, n_pear] = corr_complete( ...
            filtered_hr_bpm(mask), yue_hr_bpm(mask), ...
            "Pearson", min_points_for_corr);

        row.("R_Spearman_" + window_name) = r_spear;
        row.("P_Spearman_" + window_name) = p_spear;
        row.("N_Spearman_" + window_name) = n_spear;
        row.("R_Pearson_" + window_name) = r_pear;
        row.("P_Pearson_" + window_name) = p_pear;
        row.("N_Pearson_" + window_name) = n_pear;
    end

    by_mouse_table = [by_mouse_table; row]; %#ok<AGROW>

    plot_mouse_diagnostics(subject_name, timepts_min, filtered_hr_bpm, ...
        yue_hr_bpm, vessel_hr_bpm, mouse_fig_folder, row, cutoff_min);
end

ttest_table = build_ttest_table(by_mouse_table);
wilcoxon_table = build_wilcoxon_table(by_mouse_table);

by_mouse_csv = fullfile(output_folder, "hr_agreement_yue_filtered_by_mouse.csv");
ttest_csv = fullfile(output_folder, "hr_agreement_yue_filtered_ttest_summary.csv");
wilcoxon_csv = fullfile(output_folder, "hr_agreement_yue_filtered_wilcoxon_summary.csv");
excluded_csv = fullfile(output_folder, "hr_agreement_yue_filtered_excluded_subjects.csv");
excel_file = fullfile(output_folder, "hr_agreement_yue_filtered_results.xlsx");

writetable(by_mouse_table, by_mouse_csv);
writetable(ttest_table, ttest_csv);
writetable(wilcoxon_table, wilcoxon_csv);
writetable(excluded_table, excluded_csv);
writetable(by_mouse_table, excel_file, "Sheet", "ByMouse");
writetable(ttest_table, excel_file, "Sheet", "TTest");
writetable(wilcoxon_table, excel_file, "Sheet", "Wilcoxon");
writetable(excluded_table, excel_file, "Sheet", "Excluded");

plot_group_summary(by_mouse_table, ttest_table, wilcoxon_table, fig_folder);

fprintf("Saved %s\n", by_mouse_csv);
fprintf("Saved %s\n", ttest_csv);
fprintf("Saved %s\n", wilcoxon_csv);
fprintf("Saved %s\n", excluded_csv);
fprintf("Saved %s\n", excel_file);

%% Local functions
function hr_binned = rebin_ecg_hr_to_legacy_bins(t_hr, hr_inst_bpm, ...
    valid_inst_hr, timepts_min, bin_minutes)

    hr_binned = nan(size(timepts_min));
    keep = isfinite(t_hr) & isfinite(hr_inst_bpm) & valid_inst_hr;
    t_hr = t_hr(keep);
    hr_inst_bpm = hr_inst_bpm(keep);

    for k = 1:numel(timepts_min)
        t0 = timepts_min(k) * 60;
        t1 = (timepts_min(k) + bin_minutes) * 60;
        in_bin = t_hr >= t0 & t_hr < t1;
        if any(in_bin)
            hr_binned(k) = median(hr_inst_bpm(in_bin), "omitnan");
        end
    end
end

function [r, p, n] = corr_complete(x, y, corr_type, min_points)
    x = x(:);
    y = y(:);
    valid = isfinite(x) & isfinite(y);
    n = nnz(valid);

    if n < min_points
        r = NaN;
        p = NaN;
        return
    end

    [r, p] = corr(x(valid), y(valid), "Type", corr_type);
end

function ttest_table = build_ttest_table(by_mouse_table)
    corr_types = ["Spearman", "Pearson"];
    windows = ["Early", "Late", "Overall"];
    ttest_table = table();

    for c = 1:numel(corr_types)
        corr_type = corr_types(c);
        for w = 1:numel(windows)
            window_name = windows(w);
            r_col = "R_" + corr_type + "_" + window_name;
            r_values = by_mouse_table.(r_col);
            r_values = r_values(isfinite(r_values));

            [p_ttest, t_stat, df, mean_z, ci_low_z, ci_high_z, n_subjects] = ...
                fisher_ttest_complete(r_values);

            ttest_table = [ttest_table; table(corr_type, window_name, ...
                n_subjects, mean(r_values, "omitnan"), ...
                median(r_values, "omitnan"), mean_z, tanh(mean_z), ...
                tanh(ci_low_z), tanh(ci_high_z), p_ttest, t_stat, df, ...
                'VariableNames', {'CorrelationType', 'Window', ...
                'N_Subjects', 'MeanR', 'MedianR', 'MeanFisherZ', ...
                'BackTransformedMeanR', 'CI_Low_R', 'CI_High_R', ...
                'P_TTest', 'TStatistic', 'DegreesOfFreedom'})]; %#ok<AGROW>
        end
    end
end

function wilcoxon_table = build_wilcoxon_table(by_mouse_table)
    corr_types = ["Spearman", "Pearson"];
    windows = ["Early", "Late", "Overall"];
    wilcoxon_table = table();

    for c = 1:numel(corr_types)
        corr_type = corr_types(c);
        for w = 1:numel(windows)
            window_name = windows(w);
            r_col = "R_" + corr_type + "_" + window_name;
            p_col = "P_" + corr_type + "_" + window_name;
            r_values = by_mouse_table.(r_col);
            mouse_p_values = by_mouse_table.(p_col);
            valid = isfinite(r_values);
            r_values = r_values(valid);
            mouse_p_values = mouse_p_values(valid);

            [p_wilcoxon, signed_rank, n_subjects] = signrank_complete(r_values);

            wilcoxon_table = [wilcoxon_table; table(corr_type, window_name, ...
                n_subjects, mean(r_values, "omitnan"), ...
                median(r_values, "omitnan"), ...
                min(r_values, [], "omitnan"), max(r_values, [], "omitnan"), ...
                p_wilcoxon, signed_rank, mean(mouse_p_values, "omitnan"), ...
                median(mouse_p_values, "omitnan"), ...
                'VariableNames', {'CorrelationType', 'Window', ...
                'N_Subjects', 'MeanR', 'MedianR', 'MinR', 'MaxR', ...
                'P_Wilcoxon', 'SignedRankStatistic', ...
                'MeanMousePValue', 'MedianMousePValue'})]; %#ok<AGROW>
        end
    end
end

function [p, t_stat, df, mean_z, ci_low_z, ci_high_z, n] = ...
    fisher_ttest_complete(r_values)

    r_values = r_values(:);
    r_values = r_values(isfinite(r_values));
    n = numel(r_values);

    if n < 2
        p = NaN;
        t_stat = NaN;
        df = NaN;
        mean_z = NaN;
        ci_low_z = NaN;
        ci_high_z = NaN;
        return
    end

    r_values = min(max(r_values, -0.999999), 0.999999);
    z_values = atanh(r_values);
    mean_z = mean(z_values, "omitnan");

    [~, p, ci, stats] = ttest(z_values, 0);
    t_stat = stats.tstat;
    df = stats.df;
    ci_low_z = ci(1);
    ci_high_z = ci(2);
end

function [p, stats_value, n] = signrank_complete(x)
    x = x(:);
    x = x(isfinite(x));
    n = numel(x);

    if n == 0
        p = NaN;
        stats_value = NaN;
        return
    end

    if all(x == 0)
        p = 1;
        stats_value = 0;
        return
    end

    [p, ~, stats] = signrank(x, 0);
    if isstruct(stats) && isfield(stats, "signedrank")
        stats_value = stats.signedrank;
    else
        stats_value = NaN;
    end
end

function plot_group_summary(by_mouse_table, ttest_table, wilcoxon_table, fig_folder)
    corr_types = ["Spearman", "Pearson"];
    windows = ["Early", "Late", "Overall"];

    fig = figure("Color", "w", "Position", [100 100 1260 720], ...
        "Visible", "off");
    set(fig, "DefaultAxesFontSize", 9);
    set(fig, "DefaultTextFontSize", 9);

    for c = 1:numel(corr_types)
        corr_type = corr_types(c);
        for w = 1:numel(windows)
            window_name = windows(w);
            subplot_idx = (c - 1) * numel(windows) + w;
            subplot(numel(corr_types), numel(windows), subplot_idx);
            hold on;

            r_col = "R_" + corr_type + "_" + window_name;
            r_values = by_mouse_table.(r_col);
            subjects_here = by_mouse_table.Subject;
            valid = isfinite(r_values);
            r_values = r_values(valid);
            subjects_here = subjects_here(valid);

            x = 1:numel(r_values);
            scatter(x, r_values, 58, [0.20 0.45 0.85], "filled");
            plot([0.5, max(numel(r_values), 1) + 0.5], [0 0], "k--", ...
                "LineWidth", 1);

            t_row = ttest_table.CorrelationType == corr_type & ...
                ttest_table.Window == window_name;
            w_row = wilcoxon_table.CorrelationType == corr_type & ...
                wilcoxon_table.Window == window_name;
            if any(t_row)
                mean_r = ttest_table.BackTransformedMeanR(t_row);
                ci_low = ttest_table.CI_Low_R(t_row);
                ci_high = ttest_table.CI_High_R(t_row);
                p_t = ttest_table.P_TTest(t_row);
                plot([0.5, max(numel(r_values), 1) + 0.5], [mean_r mean_r], ...
                    "-", "Color", [0.85 0.20 0.20], "LineWidth", 1.5);
                plot([0.5, max(numel(r_values), 1) + 0.5], [ci_low ci_low], ...
                    ":", "Color", [0.85 0.20 0.20], "LineWidth", 1.2);
                plot([0.5, max(numel(r_values), 1) + 0.5], [ci_high ci_high], ...
                    ":", "Color", [0.85 0.20 0.20], "LineWidth", 1.2);
            else
                mean_r = NaN;
                p_t = NaN;
            end

            if any(w_row)
                p_w = wilcoxon_table.P_Wilcoxon(w_row);
            else
                p_w = NaN;
            end

            title_handle = title(sprintf("%s | %s\nz t p=%.3g, W p=%.3g, mean r=%.3f", ...
                corr_type, window_name, p_t, p_w, mean_r), ...
                "Interpreter", "none", "FontSize", 10);
            title_handle.Color = "k";

            set(gca, "XTick", x, ...
                "XTickLabel", cellstr(compact_subject_ids(subjects_here)), ...
                "XTickLabelRotation", 45, "FontSize", 8);
            ylim([-1 1]);
            xlim([0.5, max(numel(r_values), 1) + 0.5]);
            ylabel("Mouse-level HR-HR correlation");
            grid on;

            ax = gca;
            ax.Color = "w";
            ax.XColor = "k";
            ax.YColor = "k";
            ax.Toolbar.Visible = "off";
        end
    end

    sgtitle("Filtered ECG HR vs avg_HR_Yue_10min group summary", ...
        "Interpreter", "none", "FontSize", 12);

    out_file = fullfile(fig_folder, "hr_agreement_yue_filtered_group_summary.png");
    exportgraphics(fig, out_file, "Resolution", 220);
    close(fig);
    fprintf("Saved %s\n", out_file);
end

function plot_mouse_diagnostics(subject_name, timepts_min, filtered_hr_bpm, ...
    yue_hr_bpm, vessel_hr_bpm, out_folder, row, cutoff_min)

    fig = figure("Color", "w", "Position", [100 100 980 780], ...
        "Visible", "off");
    tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact");

    nexttile;
    early = timepts_min < cutoff_min;
    late = timepts_min >= cutoff_min;
    valid_early = early & isfinite(filtered_hr_bpm) & isfinite(yue_hr_bpm);
    valid_late = late & isfinite(filtered_hr_bpm) & isfinite(yue_hr_bpm);
    scatter(filtered_hr_bpm(valid_early), yue_hr_bpm(valid_early), 58, ...
        [0.10 0.45 0.80], "filled", "DisplayName", "Early");
    hold on;
    scatter(filtered_hr_bpm(valid_late), yue_hr_bpm(valid_late), 58, ...
        [0.85 0.35 0.10], "filled", "DisplayName", "Late");
    valid_all = isfinite(filtered_hr_bpm) & isfinite(yue_hr_bpm);
    if nnz(valid_all) >= 2
        x_fit = linspace(min(filtered_hr_bpm(valid_all)), ...
            max(filtered_hr_bpm(valid_all)), 100);
        coeffs = polyfit(filtered_hr_bpm(valid_all), yue_hr_bpm(valid_all), 1);
        plot(x_fit, polyval(coeffs, x_fit), "k-", ...
            "LineWidth", 1.4, "DisplayName", "Fit");
    end
    xlabel("Filtered ECG HR (BPM)");
    ylabel("avg\_HR\_Yue\_10min (BPM)");
    title_handle = title(sprintf("%s | Overall Pearson r=%.3f | Spearman r=%.3f", ...
        subject_name, row.R_Pearson_Overall, row.R_Spearman_Overall), ...
        "Interpreter", "none");
    title_handle.Color = "k";
    legend("Location", "best");
    grid on;
    ax1 = gca;
    ax1.Color = "w";
    ax1.XColor = "k";
    ax1.YColor = "k";
    ax1.Toolbar.Visible = "off";

    nexttile;
    plot(timepts_min, filtered_hr_bpm, "-o", "Color", [0.10 0.45 0.80], ...
        "LineWidth", 1.4, "MarkerFaceColor", [0.10 0.45 0.80], ...
        "DisplayName", "Filtered ECG HR");
    hold on;
    plot(timepts_min, yue_hr_bpm, "-s", "Color", [0.85 0.35 0.10], ...
        "LineWidth", 1.4, "MarkerFaceColor", [0.85 0.35 0.10], ...
        "DisplayName", "avg\_HR\_Yue\_10min");
    plot(timepts_min, vessel_hr_bpm, "-^", "Color", [0.20 0.60 0.35], ...
        "LineWidth", 1.4, "MarkerFaceColor", [0.20 0.60 0.35], ...
        "DisplayName", "results.hr x 60");
    xline(cutoff_min, "k--", "LineWidth", 1, "DisplayName", "60 min cutoff");
    xlabel("Time (min)");
    ylabel("HR (BPM)");
    subtitle(sprintf("Filtered - Yue: MAD = %.2f BPM | RMSE = %.2f BPM", ...
        row.MeanAbsDiffBPM, row.RMSE_BPM), "Interpreter", "none");
    legend("Location", "best");
    grid on;
    ax2 = gca;
    ax2.Color = "w";
    ax2.XColor = "k";
    ax2.YColor = "k";
    ax2.Toolbar.Visible = "off";

    out_file = fullfile(out_folder, subject_name + "_filtered_yue_diagnostics.png");
    exportgraphics(fig, out_file, "Resolution", 220);
    close(fig);
    fprintf("Saved %s\n", out_file);
end

function ids = compact_subject_ids(subjects_here)
    ids = strings(size(subjects_here));
    for i = 1:numel(subjects_here)
        token = regexp(char(subjects_here(i)), "^[A-Za-z]\d+", "match", "once");
        if isempty(token)
            ids(i) = subjects_here(i);
        else
            ids(i) = string(token);
        end
    end
end
