%% Sliding-bin HR-vessel cross-correlation from hr_amp_data_sliding.mat files
%
% Lag convention:
%   positive lag = HR leads vessel by that many minutes.
%   negative lag = vessel leads HR by that many minutes.

clear;

data_folder = ".\data\";
fig_folder = ".\figures\xcorr_sliding\";
per_subject_folder = fullfile(fig_folder, "per_subject");
res_folder = ".\results\";
if ~exist(fig_folder, "dir"); mkdir(fig_folder); end
if ~exist(per_subject_folder, "dir"); mkdir(per_subject_folder); end
if ~exist(res_folder, "dir"); mkdir(res_folder); end

hr_var = "hr_median_bpm";
vessel_vars = ["amp_svm", "bp_svm", "amp_c", "diam_mean"];
dubious_exclusion_threshold = 0.50;
max_lag_min = 30;
data_file_name = "hr_amp_data_sliding.mat";

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, "."));

xcorr_table = table();
peak_table = table();
excluded_table = table();

for v = 1:numel(vessel_vars)
    vessel_var = vessel_vars(v);
    subject_curves = [];
    subject_names = strings(0, 1);
    lag_minutes = [];

    for i = 1:length(subjects)
        subject_name = string(subjects(i).name);
        data_file = fullfile(data_folder, subject_name, data_file_name);
        if ~exist(data_file, "file")
            fprintf("Skipping %s - no %s found\n", subject_name, data_file_name);
            continue
        end

        data = load(data_file);
        dubious_fraction = data.dubious_hr_bin_fraction;

        if dubious_fraction > dubious_exclusion_threshold
            fprintf("Excluding %s - dubious HR bins %.1f%% > %.1f%%\n", ...
                subject_name, 100 * dubious_fraction, ...
                100 * dubious_exclusion_threshold);
            excluded_table = add_excluded_once(excluded_table, subject_name, dubious_fraction);
            continue
        end

        hr = data.(hr_var);
        vessel = data.(vessel_var);
        bin_step_minutes = data.bin_step_minutes;
        max_lag_bins = round(max_lag_min / bin_step_minutes);
        lags = -max_lag_bins:max_lag_bins;
        lag_minutes = lags * bin_step_minutes;

        [xc, n_pairs] = lagged_corr(hr, vessel, lags);

        lag_labels = "Lag_" + signed_label(lag_minutes) + "min";
        row = array2table(xc, "VariableNames", cellstr(lag_labels));
        row = addvars(row, subject_name, hr_var, vessel_var, ...
            dubious_fraction, bin_step_minutes, max_lag_min, ...
            "Before", 1, "NewVariableNames", ...
            ["Subject", "HrVariable", "VesselVariable", ...
            "DubiousHrBinFraction", "BinStepMinutes", "MaxLagMinutes"]);
        xcorr_table = [xcorr_table; row];

        [peak_pos_corr, idx_pos] = max(xc);
        [peak_abs_corr, idx_abs] = max(abs(xc));
        peak_table = [peak_table; table(subject_name, hr_var, vessel_var, ...
            peak_pos_corr, lag_minutes(idx_pos), xc(idx_abs), ...
            lag_minutes(idx_abs), max(n_pairs), dubious_fraction, ...
            'VariableNames', {'Subject', 'HrVariable', 'VesselVariable', ...
            'PeakPositiveCorr', 'PeakPositiveLagMinutes', ...
            'PeakAbsCorr', 'PeakAbsLagMinutes', 'MaxNPairs', ...
            'DubiousHrBinFraction'})];

        subject_curves = [subject_curves; xc];
        subject_names(end + 1, 1) = subject_name;

        plot_subject_xcorr(subject_name, hr_var, vessel_var, ...
            lag_minutes, xc, per_subject_folder);

        fprintf("Done: %s %s vs %s\n", subject_name, hr_var, vessel_var);
    end

    if ~isempty(subject_curves)
        plot_group_xcorr(hr_var, vessel_var, lag_minutes, ...
            subject_curves, subject_names, fig_folder);
    end
end

fprintf("\nExcluded subjects due to > %.1f%% dubious HR bins:\n", ...
    100 * dubious_exclusion_threshold);
disp(excluded_table);

fprintf("\nSliding-bin cross-correlation peaks:\n");
disp(peak_table);

writetable(xcorr_table, fullfile(res_folder, "hr_amp_xcorr_summary_sliding.csv"));
writetable(peak_table, fullfile(res_folder, "hr_amp_xcorr_peaks_sliding.csv"));
writetable(excluded_table, fullfile(res_folder, "hr_amp_xcorr_excluded_subjects_sliding.csv"));

%% Local functions
function [xc, n_pairs] = lagged_corr(hr, vessel, lags)
    hr = hr(:);
    vessel = vessel(:);
    xc = nan(1, numel(lags));
    n_pairs = zeros(1, numel(lags));

    for j = 1:numel(lags)
        lag = lags(j);
        if lag > 0
            x = hr(1:end-lag);
            y = vessel(1+lag:end);
        elseif lag < 0
            shift = abs(lag);
            x = hr(1+shift:end);
            y = vessel(1:end-shift);
        else
            x = hr;
            y = vessel;
        end

        valid = isfinite(x) & isfinite(y);
        n_pairs(j) = nnz(valid);
        if n_pairs(j) >= 3 && std(x(valid)) > 0 && std(y(valid)) > 0
            xc(j) = corr(x(valid), y(valid), "Type", "Pearson");
        end
    end
end

function labels = signed_label(values)
    labels = strings(size(values));
    for i = 1:numel(values)
        if values(i) > 0
            labels(i) = "plus" + string(values(i));
        elseif values(i) < 0
            labels(i) = "minus" + string(abs(values(i)));
        else
            labels(i) = "0";
        end
    end
end

function excluded_table = add_excluded_once(excluded_table, subject_name, dubious_fraction)
    if isempty(excluded_table) || ~any(excluded_table.Subject == subject_name)
        excluded_table = [excluded_table; table(subject_name, dubious_fraction, ...
            'VariableNames', {'Subject', 'DubiousHrBinFraction'})];
    end
end

function plot_subject_xcorr(subject_name, hr_var, vessel_var, lag_minutes, xc, out_folder)
    fig = figure("Position", [100 100 650 430], "Visible", "off", "Color", "w");
    set(fig, "Renderer", "painters");
    stem(lag_minutes, xc, "filled", "LineWidth", 1.2);
    xline(0, "k--");
    yline(0, "k:");
    xlabel("Lag (min): negative = vessel leads HR; positive = HR leads vessel");
    ylabel("Pearson cross-correlation");
    title(sprintf("%s: sliding %s vs %s", subject_name, hr_var, vessel_var), ...
        "Interpreter", "none");
    text(0.02, 0.94, "Vessel leads HR", "Units", "normalized", ...
        "HorizontalAlignment", "left", "FontSize", 9);
    text(0.98, 0.94, "HR leads vessel", "Units", "normalized", ...
        "HorizontalAlignment", "right", "FontSize", 9);
    grid on;
    ax = gca;
    ax.Toolbar.Visible = "off";
    ylim([-1 1]);

    out_name = subject_name + "_" + hr_var + "_vs_" + vessel_var + ...
        "_xcorr_sliding.png";
    exportgraphics(fig, fullfile(out_folder, out_name), "Resolution", 200);
    close(fig);
end

function plot_group_xcorr(hr_var, vessel_var, lag_minutes, subject_curves, ...
    subject_names, out_folder)

    mean_curve = mean(subject_curves, 1, "omitnan");
    sem_curve = std(subject_curves, 0, 1, "omitnan") ./ ...
        sqrt(sum(isfinite(subject_curves), 1));

    fig = figure("Position", [100 100 760 500], "Visible", "off", "Color", "w");
    set(fig, "Renderer", "painters");
    hold on;

    h_sem = fill([lag_minutes, fliplr(lag_minutes)], ...
        [mean_curve + sem_curve, fliplr(mean_curve - sem_curve)], ...
        [0.75 0.82 0.92], "EdgeColor", "none", "FaceAlpha", 0.6);
    h_subjects = plot(lag_minutes, subject_curves', "Color", [0.70 0.70 0.70], ...
        "LineWidth", 0.8);
    h_mean = plot(lag_minutes, mean_curve, "b-", "LineWidth", 2.2);
    xline(0, "k--");
    yline(0, "k:");

    xlabel("Lag (min): negative = vessel leads HR; positive = HR leads vessel");
    ylabel("Pearson cross-correlation");
    title(sprintf("Across subjects: sliding %s vs %s (n=%d)", ...
        hr_var, vessel_var, numel(subject_names)), "Interpreter", "none");
    text(0.02, 0.94, "Vessel leads HR", "Units", "normalized", ...
        "HorizontalAlignment", "left", "FontSize", 9);
    text(0.98, 0.94, "HR leads vessel", "Units", "normalized", ...
        "HorizontalAlignment", "right", "FontSize", 9);
    legend([h_sem, h_subjects(1), h_mean], ...
        ["Mean +/- SEM", "Subjects", "Mean"], "Location", "best");
    grid on;
    ylim([-1 1]);
    ax = gca;
    ax.Toolbar.Visible = "off";

    out_name = hr_var + "_vs_" + vessel_var + "_xcorr_sliding_across_subjects.png";
    exportgraphics(fig, fullfile(out_folder, out_name), "Resolution", 200);
    close(fig);
end
