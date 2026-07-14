%% Plot Spearman HR-vessel scatter diagnostics and Wilcoxon summaries

project_root = fileparts(fileparts(mfilename('fullpath')));
data_folder = fullfile(project_root, 'data');
output_folder = fullfile(project_root, 'spearman_correlation');
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

by_mouse_csv = fullfile(output_folder, 'hr_amp_corr_spearman_by_mouse.csv');
wilcoxon_csv = fullfile(output_folder, 'hr_amp_corr_spearman_wilcoxon_summary.csv');

if ~exist(by_mouse_csv, 'file') || ~exist(wilcoxon_csv, 'file')
    error('Missing Spearman correlation outputs. Run compute_hr_diam_corr_spearman.m first.');
end

hr_var = 'hr_median_bpm';
vessel_vars = {'amp_svm', 'bp_svm', 'amp_c', 'diam_mean'};
cutoff_min = 60;
dubious_exclusion_threshold = 0.50;

by_mouse_table = readtable(by_mouse_csv, 'TextType', 'string');
wilcoxon_table = readtable(wilcoxon_csv, 'TextType', 'string');

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, '.'));

for i = 1:length(subjects)
    subject_name = string(subjects(i).name);
    data_file = fullfile(data_folder, subject_name, 'hr_amp_data.mat');
    if ~exist(data_file, 'file')
        continue
    end

    data = load(data_file);
    dubious_fraction = data.dubious_hr_bin_fraction;
    if dubious_fraction > dubious_exclusion_threshold
        continue
    end

    timepts = data.timepts(:);
    hr = data.(hr_var)(:);
    early = timepts < cutoff_min;
    late = timepts >= cutoff_min;

    fig = figure('Color', 'w', 'Position', [100 100 1200 860], 'Visible', 'off');

    for v = 1:numel(vessel_vars)
        vessel_var = vessel_vars{v};
        vessel = data.(vessel_var)(:);

        subplot(2, 2, v);
        hold on;

        valid_early = early & isfinite(hr) & isfinite(vessel);
        valid_late = late & isfinite(hr) & isfinite(vessel);
        scatter(hr(valid_early), vessel(valid_early), 48, [0.10 0.45 0.85], ...
            'filled', 'MarkerFaceAlpha', 0.75);
        scatter(hr(valid_late), vessel(valid_late), 48, [0.90 0.45 0.10], ...
            'filled', 'MarkerFaceAlpha', 0.75);

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
            title(sprintf('%s\nOverall rho=%.3f, p=%.3g | Early rho=%.3f | Late rho=%.3f', ...
                vessel_label(vessel_var), ...
                by_mouse_table.R_Overall(row), by_mouse_table.P_Overall(row), ...
                by_mouse_table.R_Early(row), by_mouse_table.R_Late(row)), ...
                'Interpreter', 'none');
        else
            title(vessel_label(vessel_var), 'Interpreter', 'none');
        end

        xlabel('Heart rate (BPM)');
        ylabel(vessel_label(vessel_var));
        grid on;
        legend({'Early bins', 'Late bins', 'Linear fit guide'}, 'Location', 'best');

        ax = gca;
        ax.Color = 'w';
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.Toolbar.Visible = 'off';
    end

    sgtitle(sprintf('%s: Spearman HR-vessel scatter diagnostics', subject_name), ...
        'Interpreter', 'none');

    out_file = fullfile(scatter_folder, subject_name + "_spearman_scatter.png");
    exportgraphics(fig, out_file, 'Resolution', 220);
    close(fig);
    fprintf('Saved %s\n', out_file);
end

windows = {'Early', 'Late', 'Overall'};
fig = figure('Color', 'w', 'Position', [100 80 1260 980], 'Visible', 'off');
set(fig, 'DefaultAxesFontSize', 9);
set(fig, 'DefaultTextFontSize', 9);

summary_vessel_vars = {'amp_svm', 'bp_svm', 'amp_c'};
for v = 1:numel(summary_vessel_vars)
    vessel_var = summary_vessel_vars{v};
    vessel_rows = by_mouse_table.VesselVariable == string(vessel_var);

    for w = 1:numel(windows)
        window_name = windows{w};
        subplot_idx = (v - 1) * numel(windows) + w;
        subplot(numel(summary_vessel_vars), numel(windows), subplot_idx);
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

        set(gca, 'XTick', x, 'XTickLabel', subjects_here, 'XTickLabelRotation', 45);
        ylim([-1 1]);
        xlim([0.5, max(numel(r_values), 1) + 0.5]);
        grid on;

        summary_row = wilcoxon_table.VesselVariable == string(vessel_var) & ...
            wilcoxon_table.Window == string(window_name);
        if any(summary_row)
            p_w = wilcoxon_table.P_Wilcoxon(summary_row);
            med_r = wilcoxon_table.MedianR(summary_row);
            title_handle = title(sprintf('%s | %s\nWilcoxon p=%.3g, median rho=%.3f', ...
                vessel_label(vessel_var), window_name, p_w, med_r), ...
                'Interpreter', 'none', 'FontSize', 10);
        else
            title_handle = title(sprintf('%s | %s', vessel_label(vessel_var), window_name), ...
                'Interpreter', 'none', 'FontSize', 10);
        end
        title_handle.Color = 'k';

        set(gca, 'XTick', x, 'XTickLabel', cellstr(subject_labels), ...
            'XTickLabelRotation', 45, 'FontSize', 7);
        ylabel('Mouse-level Spearman rho', 'FontSize', 9);

        ax = gca;
        ax.Color = 'w';
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.Toolbar.Visible = 'off';
    end
end

sgtitle('Across-mouse Spearman correlations with Wilcoxon summary', ...
    'Interpreter', 'none', 'FontSize', 12);

out_file = fullfile(summary_folder, 'spearman_wilcoxon_summary.png');
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
