%% Plot subject-level HR-vessel correlation distributions
% Reads results/hr_amp_corr_summary.csv and saves one boxplot per vessel
% feature directly under figures/.

clear; close all;

results_file = fullfile(pwd, "results", "hr_amp_corr_summary.csv");
fig_folder = fullfile(pwd, "figures");
if ~exist(fig_folder, "dir")
    mkdir(fig_folder);
end

if ~exist(results_file, "file")
    error("Missing %s. Run compute_hr_diam_corr.m first.", results_file);
end

corr_table = readtable(results_file, "TextType", "string");
vessel_vars = unique(corr_table.VesselVariable, "stable");

for i = 1:numel(vessel_vars)
    vessel_var = vessel_vars(i);
    rows = corr_table.VesselVariable == vessel_var;

    values = [
        corr_table.R_Early(rows);
        corr_table.R_Late(rows);
        corr_table.R_Overall(rows)
    ];
    groups = [
        repmat("Early", nnz(rows), 1);
        repmat("Late", nnz(rows), 1);
        repmat("Overall", nnz(rows), 1)
    ];

    valid = isfinite(values);
    values = values(valid);
    groups = categorical(groups(valid), ["Early", "Late", "Overall"], ...
        "Ordinal", true);

    fig = figure("Color", "w", "Position", [100 100 560 460], ...
        "Visible", "off");

    boxchart(groups, values, "BoxFaceColor", [0.2 0.45 0.85], ...
        "MarkerStyle", "o", "JitterOutliers", "on");
    hold on;
    yline(0, "k--", "LineWidth", 1);
    ylim([-1 1]);
    ylabel("Spearman correlation with HR");
    xlabel("Time window");
    title_handle = title(sprintf("HR vs %s correlation across subjects", vessel_var), ...
        "Interpreter", "none");
    title_handle.Color = "k";
    grid on;

    ax = gca;
    ax.Color = "w";
    ax.XColor = "k";
    ax.YColor = "k";
    ax.Toolbar.Visible = "off";
    set(findall(fig, "Type", "text"), "Color", "k");

    out_file = fullfile(fig_folder, "hr_corr_boxplot_" + vessel_var + ".png");
    exportgraphics(fig, out_file, "Resolution", 200);
    close(fig);
    fprintf("Saved %s\n", out_file);
end
