%% Plot fine-resolution ECG heart rate and slow vasomotion amplitude for F153
% Uses the compact data file produced by extract_hr_amp_data_F153.m.

clear; close all;

data_folder = fullfile(pwd, "data");
fig_folder = fullfile(pwd, "figures");
subject_name = "F153-Dex40-06022025";
subject_folder = fullfile(data_folder, subject_name);

if ~exist(fig_folder, "dir")
    mkdir(fig_folder);
end

data_file = fullfile(subject_folder, "hr_amp_data.mat");
if ~exist(data_file, "file")
    error("Missing %s. Run extract_hr_amp_data_F153.m first.", data_file);
end

data = load(data_file, ...
    "bin_centers_min", "timepts", "bin_minutes", ...
    "hr_median_bpm", "hr_median_bpm_raw", ...
    "amp_svm", "is_dubious_hr_bin", "dubious_hr_bin_fraction", ...
    "Fc_svm", "subject_name");

if isfield(data, "bin_centers_min")
    t_plot = data.bin_centers_min;
else
    t_plot = data.timepts + data.bin_minutes / 2;
end

fig = figure("Color", "w", "Position", [100 100 900 460]);

yyaxis left
plot(t_plot, data.hr_median_bpm, "ro-", "LineWidth", 1.8, ...
    "MarkerFaceColor", "r", "DisplayName", "Median HR");
hold on;
if any(data.is_dubious_hr_bin)
    plot(t_plot(data.is_dubious_hr_bin), ...
        data.hr_median_bpm(data.is_dubious_hr_bin), "ks", ...
        "MarkerFaceColor", "y", "MarkerSize", 7, ...
        "DisplayName", "Interpolated HR bin");
end
ylabel("Heart rate (BPM)");

yyaxis right
plot(t_plot, data.amp_svm, "b*-", "LineWidth", 1.8, ...
    "MarkerSize", 7, "DisplayName", "amp\_svm");
ylabel(sprintf("Slow vasomotion IQR (um), %.2f-%.1f Hz", ...
    data.Fc_svm(1), data.Fc_svm(2)));

xlabel("Time (min)");
title(sprintf("%s: 5-min HR and slow vasomotion amplitude", data.subject_name));
grid on;

ax = gca;
ax.Color = "w";
ax.XColor = "k";
ax.YAxis(1).Color = "r";
ax.YAxis(2).Color = "b";
ax.Toolbar.Visible = "off";

legend("Location", "best");

annotation_text = sprintf("Dubious HR bins interpolated: %.1f%%", ...
    100 * data.dubious_hr_bin_fraction);
text(0.02, 0.95, annotation_text, "Units", "normalized", ...
    "Color", "k", "FontSize", 9, "VerticalAlignment", "top");

fig_file = fullfile(fig_folder, subject_name + "_fine_hr_amp_svm.png");
exportgraphics(fig, fig_file, "Resolution", 200);
fprintf("Saved %s\n", fig_file);
