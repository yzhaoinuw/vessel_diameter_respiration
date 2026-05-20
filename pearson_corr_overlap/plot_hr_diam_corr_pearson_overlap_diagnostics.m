%% Plot overlap Pearson HR-vessel scatter diagnostics and Wilcoxon summaries

project_root = fileparts(fileparts(mfilename('fullpath')));
data_folder = fullfile(project_root, 'data');
output_folder = fullfile(project_root, 'pearson_corr_overlap');
fig_folder = fullfile(output_folder, 'figures');
scatter_folder = fullfile(fig_folder, 'scatter_by_mouse');
summary_folder = fullfile(fig_folder, 'wilcoxon_summary');
if ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end
if ~exist(scatter_folder, 'dir')
    mkdir(scatter_folder);
end
if ~exist(summary_folder, 'dir')
    mkdir(summary_folder);
end

by_mouse_csv = fullfile(output_folder, 'hr_amp_corr_overlap_pearson_by_mouse.csv');
wilcoxon_csv = fullfile(output_folder, 'hr_amp_corr_overlap_pearson_wilcoxon_summary.csv');

if ~exist(by_mouse_csv, 'file') || ~exist(wilcoxon_csv, 'file')
    error('Missing overlap Pearson outputs. Run compute_hr_diam_corr_pearson_overlap.m first.');
end

vessel_vars = {'amp_svm', 'bp_svm', 'amp_c'};
cutoff_min = 60;

by_mouse_table = readtable(by_mouse_csv, 'TextType', 'string');
wilcoxon_table = readtable(wilcoxon_csv, 'TextType', 'string');

subjects = unique(by_mouse_table.Subject, 'stable');

for i = 1:numel(subjects)
    subject_name = string(subjects(i));
    subject_folder = fullfile(data_folder, subject_name);
    hr_file = fullfile(subject_folder, 'hr_amp_data.mat');
    vessel_file = fullfile(subject_folder, 'vessel_features_10min_overlap.mat');
    if ~exist(hr_file, 'file') || ~exist(vessel_file, 'file')
        continue
    end

    hr_data = load(hr_file, 'hr_filter_timepts_min', 'hr_minute_bpm_filtered');
    vessel_data = load(vessel_file, 'timepts', 'amp_svm', 'bp_svm', 'amp_c');

    [timepts, hr] = align_hr_to_vessel_starts( ...
        vessel_data.timepts, hr_data.hr_filter_timepts_min, ...
        hr_data.hr_minute_bpm_filtered);
    early = timepts < cutoff_min;
    late = timepts >= cutoff_min;

    fig = figure('Color', 'w', 'Position', [100 100 1280 520], 'Visible', 'off');

    for v = 1:numel(vessel_vars)
        vessel_var = vessel_vars{v};
        vessel = vessel_data.(vessel_var)(:);

        subplot(1, 3, v);
        hold on;

        valid_early = early & isfinite(hr) & isfinite(vessel);
        valid_late = late & isfinite(hr) & isfinite(vessel);
        scatter(hr(valid_early), vessel(valid_early), 32, [0.10 0.45 0.85], ...
            'filled', 'MarkerFaceAlpha', 0.65);
        scatter(hr(valid_late), vessel(valid_late), 32, [0.90 0.45 0.10], ...
            'filled', 'MarkerFaceAlpha', 0.65);

        valid_all = isfinite(hr) & isfinite(vessel);
        if nnz(valid_all) >= 2
            coeffs = polyfit(hr(valid_all), vessel(valid_all), 1);
            xfit = linspace(min(hr(valid_all)), max(hr(valid_all)), 100);
            yfit = polyval(coeffs, xfit);
            plot(xfit, yfit, 'k-', 'LineWidth', 1.5);
        end

        row = by_mouse_table.Subject == subject_name & ...
            by_mouse_table.VesselVariable == string(vessel_var);

        if any(row)
            title(sprintf('%s\nOverall r=%.3f, p=%.3g | Early r=%.3f | Late r=%.3f', ...
                vessel_label(vessel_var), ...
                by_mouse_table.R_Overall(row), by_mouse_table.P_Overall(row), ...
                by_mouse_table.R_Early(row), by_mouse_table.R_Late(row)), ...
                'Interpreter', 'none');
        else
            title(vessel_label(vessel_var), 'Interpreter', 'none');
        end

        xlabel('Filtered ECG HR (BPM)');
        ylabel(vessel_label(vessel_var));
        grid on;
        legend({'Early bins', 'Late bins', 'Pearson fit'}, 'Location', 'best');

        ax = gca;
        ax.Color = 'w';
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.Toolbar.Visible = 'off';
    end

    title_handle = sgtitle(sprintf('%s: overlap Pearson HR-vessel scatter diagnostics', subject_name), ...
        'Interpreter', 'none');
    title_handle.Color = 'k';

    out_file = fullfile(scatter_folder, subject_name + "_pearson_overlap_scatter.png");
    exportgraphics(fig, out_file, 'Resolution', 220);
    close(fig);
    fprintf('Saved %s\n', out_file);
end

windows = {'Early', 'Late', 'Overall'};
fig = figure('Color', 'w', 'Position', [100 80 1260 980], 'Visible', 'off');
set(fig, 'DefaultAxesFontSize', 9);
set(fig, 'DefaultTextFontSize', 9);

for v = 1:numel(vessel_vars)
    vessel_var = vessel_vars{v};
    vessel_rows = by_mouse_table.VesselVariable == string(vessel_var);

    for w = 1:numel(windows)
        window_name = windows{w};
        subplot_idx = (v - 1) * numel(windows) + w;
        subplot(numel(vessel_vars), numel(windows), subplot_idx);
        hold on;

        r_col = "R_" + window_name;
        r_values = by_mouse_table.(r_col)(vessel_rows);
        subjects_here = by_mouse_table.Subject(vessel_rows);
        valid = isfinite(r_values);
        r_values = r_values(valid);
        subjects_here = subjects_here(valid);
        subject_labels = compact_subject_ids(subjects_here);

        x = 1:numel(r_values);
        scatter(x, r_values, 58, [0.20 0.45 0.85], 'filled');
        plot([0.5, numel(r_values) + 0.5], [0 0], 'k--', 'LineWidth', 1);

        if ~isempty(r_values)
            med_r = median(r_values, 'omitnan');
            plot([0.5, numel(r_values) + 0.5], [med_r med_r], '-', ...
                'Color', [0.85 0.20 0.20], 'LineWidth', 1.5);
        end

        ylim([-1 1]);
        xlim([0.5, max(numel(r_values), 1) + 0.5]);
        grid on;

        summary_row = wilcoxon_table.VesselVariable == string(vessel_var) & ...
            wilcoxon_table.Window == string(window_name);
        if any(summary_row)
            p_w = wilcoxon_table.P_Wilcoxon(summary_row);
            med_r = wilcoxon_table.MedianR(summary_row);
            title_handle = title(sprintf('%s | %s\nWilcoxon p=%.3g, median r=%.3f', ...
                vessel_label(vessel_var), window_name, p_w, med_r), ...
                'Interpreter', 'none', 'FontSize', 10);
        else
            title_handle = title(sprintf('%s | %s', vessel_label(vessel_var), window_name), ...
                'Interpreter', 'none', 'FontSize', 10);
        end
        title_handle.Color = 'k';

        set(gca, 'XTick', x, 'XTickLabel', cellstr(subject_labels), ...
            'XTickLabelRotation', 45, 'FontSize', 7);
        ylabel('Mouse-level Pearson r', 'FontSize', 9);

        ax = gca;
        ax.Color = 'w';
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.Toolbar.Visible = 'off';
    end
end

title_handle = sgtitle('Across-mouse overlap Pearson correlations with Wilcoxon summary', ...
    'Interpreter', 'none', 'FontSize', 12);
title_handle.Color = 'k';

out_file = fullfile(summary_folder, 'pearson_overlap_wilcoxon_summary.png');
exportgraphics(fig, out_file, 'Resolution', 220);
close(fig);
fprintf('Saved %s\n', out_file);

%% Local functions
function [timepts, hr] = align_hr_to_vessel_starts(vessel_timepts, hr_timepts, hr_values)
    timepts = double(vessel_timepts(:));
    hr_timepts = double(hr_timepts(:));
    hr_values = double(hr_values(:));

    [matched, loc] = ismember(timepts, hr_timepts);
    hr = nan(size(timepts));
    hr(matched) = hr_values(loc(matched));
end

function label = vessel_label(vessel_var)
    switch string(vessel_var)
        case "amp_svm"
            label = "amp\_svm";
        case "bp_svm"
            label = "bp\_svm";
        case "amp_c"
            label = "amp\_c";
        otherwise
            label = char(string(vessel_var));
    end
end

function ids = compact_subject_ids(subjects_here)
    ids = strings(size(subjects_here));
    for i = 1:numel(subjects_here)
        token = regexp(char(subjects_here(i)), '^[A-Za-z]\d+', 'match', 'once');
        if isempty(token)
            ids(i) = subjects_here(i);
        else
            ids(i) = string(token);
        end
    end
end
