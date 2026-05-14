%% Plot Fisher z t-test summaries for mouse-level Pearson correlations

project_root = fileparts(fileparts(mfilename('fullpath')));
output_folder = fullfile(project_root, 'pearson_corr');
fig_folder = fullfile(output_folder, 'figures', 'ttest_summary');
if ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end

by_mouse_csv = fullfile(output_folder, 'hr_amp_corr_pearson_by_mouse.csv');
ttest_csv = fullfile(output_folder, 'hr_amp_corr_pearson_ttest_summary.csv');

if ~exist(by_mouse_csv, 'file') || ~exist(ttest_csv, 'file')
    error('Missing Pearson t-test inputs. Run compute_hr_diam_corr_pearson_ttest.m first.');
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

            title_handle = title(sprintf('%s | %s\nT-test p=%.3g, mean r=%.3f', ...
                vessel_label(vessel_var), window_name, p_t, mean_r), ...
                'Interpreter', 'none', 'FontSize', 10);
        else
            title_handle = title(sprintf('%s | %s', vessel_label(vessel_var), window_name), ...
                'Interpreter', 'none', 'FontSize', 10);
        end
        title_handle.Color = 'k';

        set(gca, 'XTick', x, 'XTickLabel', cellstr(subject_labels), 'XTickLabelRotation', 45, ...
            'FontSize', 7);
        ylim([-1 1]);
        xlim([0.5, max(numel(r_values), 1) + 0.5]);
        ylabel('Mouse-level Pearson r', 'FontSize', 9);
        grid on;

        ax = gca;
        ax.Color = 'w';
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.Toolbar.Visible = 'off';
    end
end

sgtitle('Across-mouse Pearson correlations with Fisher z t-test summary', ...
    'Interpreter', 'none', 'FontSize', 12);

out_file = fullfile(fig_folder, 'pearson_ttest_summary.png');
exportgraphics(fig, out_file, 'Resolution', 220);
close(fig);
fprintf('Saved %s\n', out_file);

%% Local functions
function label = vessel_label(vessel_var)
    switch string(vessel_var)
        case "amp_svm"
            label = "amp\_svm";
        case "bp_svm"
            label = "bp\_svm";
        case "amp_c"
            label = "amp\_c";
        case "diam_mean"
            label = "diam\_mean";
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
