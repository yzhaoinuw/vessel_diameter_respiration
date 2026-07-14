%% Plot HR and vessel variables from fine-resolution hr_amp_data.mat files

data_folder = ".\data\";
fig_folder = ".\figures\hr_amp_coupling\";
if ~exist(fig_folder, "dir")
    mkdir(fig_folder);
end

hr_var = "hr_median_bpm";
vessel_vars = ["amp_svm", "bp_svm", "amp_c", "diam_mean"];
dubious_exclusion_threshold = 0.50;

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, "."));

excluded_table = table();

for i = 1:length(subjects)
    subject_name = string(subjects(i).name);
    data_file = fullfile(data_folder, subject_name, "hr_amp_data.mat");
    if ~exist(data_file, "file")
        fprintf("Skipping %s - no hr_amp_data.mat found\n", subject_name);
        continue
    end

    data = load(data_file);
    dubious_fraction = data.dubious_hr_bin_fraction;

    if dubious_fraction > dubious_exclusion_threshold
        fprintf("Excluding %s - dubious HR bins %.1f%% > %.1f%%\n", ...
            subject_name, 100 * dubious_fraction, ...
            100 * dubious_exclusion_threshold);
        excluded_table = [excluded_table; table(subject_name, dubious_fraction, ...
            'VariableNames', {'Subject', 'DubiousHrBinFraction'})];
        continue
    end

    if isfield(data, "bin_centers_min")
        t_plot = data.bin_centers_min;
    else
        t_plot = data.timepts + data.bin_minutes / 2;
    end
    hr = data.(hr_var);

    for v = 1:numel(vessel_vars)
        vessel_var = vessel_vars(v);
        vessel = data.(vessel_var);

        fig = figure("Position", [100 100 900 500], "Visible", "off", ...
            "Color", "w");
        set(fig, "Renderer", "painters");

        yyaxis left
        plot(t_plot, hr, "ro-", "LineWidth", 1.8, "MarkerFaceColor", "r");
        ylabel("Heart rate (BPM)");

        yyaxis right
        plot(t_plot, vessel, "b*-", "LineWidth", 1.8, "MarkerSize", 7);
        ylabel(vessel_label(vessel_var));

        xlabel("Time (min)");
        title(sprintf("%s: %s vs %s", subject_name, hr_var, vessel_var), ...
            "Interpreter", "none");
        legend("Heart rate", vessel_var, "Location", "best", ...
            "Interpreter", "none");
        grid on;

        ax = gca;
        ax.Color = "w";
        ax.XColor = "k";
        ax.YAxis(1).Color = "r";
        ax.YAxis(2).Color = "b";
        ax.Toolbar.Visible = "off";

        text(0.02, 0.95, ...
            sprintf("Dubious HR bins: %.1f%%", 100 * dubious_fraction), ...
            "Units", "normalized", "Color", "k", ...
            "FontSize", 9, "VerticalAlignment", "top");

        out_name = subject_name + "_" + vessel_var + ".png";
        exportgraphics(fig, fullfile(fig_folder, out_name), "Resolution", 200);
        close(fig);
        fprintf("Saved %s\n", fullfile(fig_folder, out_name));
    end
end

fprintf("\nExcluded subjects due to > %.1f%% dubious HR bins:\n", ...
    100 * dubious_exclusion_threshold);
disp(excluded_table);

writetable(excluded_table, fullfile(fig_folder, "excluded_subjects.csv"));

%% Local functions
function label = vessel_label(vessel_var)
    switch string(vessel_var)
        case "amp_svm"
            label = "Slow vasomotion IQR (um)";
        case "bp_svm"
            label = "Slow vasomotion band power (um^2)";
        case "amp_c"
            label = "Cardiac-band IQR (um)";
        case "diam_mean"
            label = "Mean vessel diameter (um)";
        otherwise
            label = string(vessel_var);
    end
end
