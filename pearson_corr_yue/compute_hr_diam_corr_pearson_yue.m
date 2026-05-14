%% Compute Pearson HR-vessel correlations using results.avg_HR_Yue_10min

project_root = fileparts(fileparts(mfilename('fullpath')));
data_folder = fullfile(project_root, 'data');
output_folder = fullfile(project_root, 'pearson_corr_yue');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

hr_var = "avg_HR_Yue_10min";
vessel_vars = ["amp_svm", "bp_svm", "amp_c", "diam_tp"];
window_names = ["Early", "Late", "Overall"];
cutoff_min = 60;
dubious_exclusion_threshold = 0.50;

excel_file = fullfile(output_folder, 'hr_amp_corr_yue_pearson_results.xlsx');
by_mouse_csv = fullfile(output_folder, 'hr_amp_corr_yue_pearson_by_mouse.csv');
wilcoxon_csv = fullfile(output_folder, 'hr_amp_corr_yue_pearson_wilcoxon_summary.csv');
excluded_csv = fullfile(output_folder, 'hr_amp_corr_yue_pearson_excluded_subjects.csv');

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, '.'));

excluded_table = table();
by_mouse_table = table();

for i = 1:length(subjects)
    subject_name = string(subjects(i).name);
    results_file = fullfile(data_folder, subject_name, 'results.mat');
    hr_amp_file = fullfile(data_folder, subject_name, 'hr_amp_data.mat');
    if ~exist(results_file, 'file')
        fprintf("Skipping %s - no results.mat found\n", subject_name);
        continue
    end
    if ~exist(hr_amp_file, 'file')
        fprintf("Skipping %s - no hr_amp_data.mat found for QC\n", subject_name);
        continue
    end

    qc_data = load(hr_amp_file, 'dubious_hr_bin_fraction');
    dubious_fraction = qc_data.dubious_hr_bin_fraction;

    if dubious_fraction > dubious_exclusion_threshold
        fprintf("Excluding %s - dubious HR bins %.1f%% > %.1f%%\n", ...
            subject_name, 100 * dubious_fraction, ...
            100 * dubious_exclusion_threshold);
        excluded_table = [excluded_table; table(subject_name, dubious_fraction, ...
            'VariableNames', {'Subject', 'DubiousHrBinFraction'})]; %#ok<AGROW>
        continue
    end

    results_data = load(results_file, 'timepts', char(hr_var), 'amp_svm', 'bp_svm', 'amp_c', 'diam_tp');
    hr = results_data.(char(hr_var));
    timepts = results_data.timepts;
    early = timepts < cutoff_min;
    late = timepts >= cutoff_min;

    for v = 1:numel(vessel_vars)
        vessel_var = vessel_vars(v);
        vessel = results_data.(char(vessel_var));

        [r_early, p_early, n_early] = corr_complete(hr(early), vessel(early));
        [r_late, p_late, n_late] = corr_complete(hr(late), vessel(late));
        [r_overall, p_overall, n_overall] = corr_complete(hr, vessel);

        by_mouse_table = [by_mouse_table; table(subject_name, hr_var, vessel_var, ...
            r_early, p_early, n_early, r_late, p_late, n_late, ...
            r_overall, p_overall, n_overall, dubious_fraction, ...
            'VariableNames', {'Subject', 'HrVariable', 'VesselVariable', ...
            'R_Early', 'P_Early', 'N_Early', ...
            'R_Late', 'P_Late', 'N_Late', ...
            'R_Overall', 'P_Overall', 'N_Overall', ...
            'DubiousHrBinFraction'})]; %#ok<AGROW>
    end
end

wilcoxon_table = table();
for v = 1:numel(vessel_vars)
    vessel_var = vessel_vars(v);
    rows = by_mouse_table.VesselVariable == vessel_var;

    for w = 1:numel(window_names)
        window_name = window_names(w);
        r_col = char("R_" + window_name);
        p_col = char("P_" + window_name);

        r_values = by_mouse_table.(r_col)(rows);
        mouse_p_values = by_mouse_table.(p_col)(rows);
        valid = isfinite(r_values);
        r_values = r_values(valid);
        mouse_p_values = mouse_p_values(valid);

        [p_wilcoxon, signed_rank, n_subjects] = signrank_complete(r_values);

        wilcoxon_table = [wilcoxon_table; table(vessel_var, window_name, ...
            n_subjects, mean(r_values, 'omitnan'), median(r_values, 'omitnan'), ...
            min(r_values, [], 'omitnan'), max(r_values, [], 'omitnan'), ...
            p_wilcoxon, signed_rank, mean(mouse_p_values, 'omitnan'), ...
            median(mouse_p_values, 'omitnan'), ...
            'VariableNames', {'VesselVariable', 'Window', ...
            'N_Subjects', 'MeanR', 'MedianR', 'MinR', 'MaxR', ...
            'P_Wilcoxon', 'SignedRankStatistic', ...
            'MeanMousePValue', 'MedianMousePValue'})]; %#ok<AGROW>
    end
end

writetable(by_mouse_table, by_mouse_csv);
writetable(wilcoxon_table, wilcoxon_csv);
writetable(excluded_table, excluded_csv);

writetable(by_mouse_table, excel_file, 'Sheet', 'ByMouse');
writetable(wilcoxon_table, excel_file, 'Sheet', 'Wilcoxon');
writetable(excluded_table, excel_file, 'Sheet', 'ExcludedSubjects');

fprintf("\nSaved Yue Pearson outputs to:\n%s\n", output_folder);
fprintf("- %s\n", excel_file);
fprintf("- %s\n", by_mouse_csv);
fprintf("- %s\n", wilcoxon_csv);
fprintf("- %s\n", excluded_csv);

%% Local functions
function [r, p, n] = corr_complete(x, y)
    x = x(:);
    y = y(:);
    valid = isfinite(x) & isfinite(y);
    n = nnz(valid);
    if n < 3
        r = NaN;
        p = NaN;
        return
    end
    [r, p] = corr(x(valid), y(valid), 'Type', 'Pearson');
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
    if isstruct(stats) && isfield(stats, 'signedrank')
        stats_value = stats.signedrank;
    else
        stats_value = NaN;
    end
end
