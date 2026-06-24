%% Plot Fisher z t-test summaries for overlap Spearman correlations

project_root = fileparts(fileparts(mfilename('fullpath')));
if ~exist('overlap_analysis_mode', 'var')
    overlap_analysis_mode = "curated";
else
    overlap_analysis_mode = string(overlap_analysis_mode);
end

switch overlap_analysis_mode
    case "curated"
        output_folder = fullfile(project_root, 'spearman_correlation_overlap');
        analysis_label = "quality-curated recordings";
    case "all_subjects"
        output_folder = fullfile(project_root, 'spearman_correlation_overlap_all_subjects');
        analysis_label = "all recordings";
    otherwise
        error('Unknown overlap_analysis_mode: %s', overlap_analysis_mode);
end

fig_folder = fullfile(output_folder, 'figures', 'ttest_summary');
if ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end

by_mouse_csv = fullfile(output_folder, 'hr_amp_corr_overlap_spearman_by_mouse.csv');
ttest_csv = fullfile(output_folder, 'hr_amp_corr_overlap_spearman_ttest_summary.csv');

if ~exist(by_mouse_csv, 'file') || ~exist(ttest_csv, 'file')
    error('Missing overlap Spearman t-test inputs. Run compute_hr_diam_corr_spearman_overlap_ttest.m first.');
end

by_mouse_table = readtable(by_mouse_csv, 'TextType', 'string');
ttest_table = readtable(ttest_csv, 'TextType', 'string');
vessel_vars = ["amp_svm", "bp_svm", "amp_c"];
window_names = ["Early", "Late", "Overall"];

fig = figure('Color', 'w', 'Position', [100 80 1260 980], 'Visible', 'off');
set(fig, 'DefaultAxesFontSize', 9);
set(fig, 'DefaultTextFontSize', 9);

for v = 1:numel(vessel_vars)
    vessel_var = vessel_vars(v);
    vessel_rows = by_mouse_table.VesselVariable == vessel_var;

    for w = 1:numel(window_names)
        window_name = window_names(w);
        subplot_idx = (v - 1) * numel(window_names) + w;
        subplot(numel(vessel_vars), numel(window_names), subplot_idx);
        hold on;

        r_col = char("R_" + window_name);
        r_values = by_mouse_table.(r_col)(vessel_rows);
        subjects_here = by_mouse_table.Subject(vessel_rows);
        valid = isfinite(r_values);
        r_values = r_values(valid);
        subjects_here = subjects_here(valid);
        subject_labels = compact_subject_ids(subjects_here);

        x = 1:numel(r_values);
        scatter(x, r_values, 58, [0.20 0.45 0.85], 'filled');
        plot([0.5, numel(r_values) + 0.5], [0 0], 'k--', 'LineWidth', 1);

        summary_row = ttest_table.VesselVariable == vessel_var & ...
            ttest_table.Window == window_name;
        if any(summary_row)
            mean_r = ttest_table.BackTransformedMeanR(summary_row);
            ci_low = ttest_table.CI_Low_R(summary_row);
            ci_high = ttest_table.CI_High_R(summary_row);
            p_t = ttest_table.P_TTest(summary_row);

            plot([0.5, numel(r_values) + 0.5], [mean_r mean_r], '-', ...
                'Color', [0.85 0.20 0.20], 'LineWidth', 1.5);
            plot([0.5, numel(r_values) + 0.5], [ci_low ci_low], ':', ...
                'Color', [0.85 0.20 0.20], 'LineWidth', 1.2);
            plot([0.5, numel(r_values) + 0.5], [ci_high ci_high], ':', ...
                'Color', [0.85 0.20 0.20], 'LineWidth', 1.2);

            title_handle = title(sprintf('%s | %s\nFisher p=%.3g, mean rho=%.3f', ...
                vessel_label(vessel_var), window_name, p_t, mean_r), ...
                'Interpreter', 'none', 'FontSize', 10);
        else
            title_handle = title(sprintf('%s | %s', vessel_label(vessel_var), window_name), ...
                'Interpreter', 'none', 'FontSize', 10);
        end
        title_handle.Color = 'k';

        set(gca, 'XTick', x, 'XTickLabel', cellstr(subject_labels), ...
            'XTickLabelRotation', 45, 'FontSize', 7);
        ylim([-1 1]);
        xlim([0.5, max(numel(r_values), 1) + 0.5]);
        ylabel('Mouse-level Spearman rho', 'FontSize', 9);
        grid on;

        ax = gca;
        ax.Color = 'w';
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.Toolbar.Visible = 'off';
    end
end

title_handle = sgtitle(sprintf('Sliding-window Spearman correlations (%s)', analysis_label), ...
    'Interpreter', 'none', 'FontSize', 12);
title_handle.Color = 'k';

out_file = fullfile(fig_folder, 'spearman_overlap_ttest_summary.png');
exportgraphics(fig, out_file, 'Resolution', 220);
close(fig);
fprintf('Saved %s\n', out_file);

%% Local functions
function label = vessel_label(vessel_var)
    switch string(vessel_var)
        case "amp_svm"
            label = "Slow vasomotion amplitude";
        case "bp_svm"
            label = "Slow vasomotion power";
        case "amp_c"
            label = "Cardiac pulsation amplitude";
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
