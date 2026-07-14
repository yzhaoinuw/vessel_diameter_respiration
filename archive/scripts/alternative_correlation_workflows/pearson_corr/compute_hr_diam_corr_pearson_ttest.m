%% Run Fisher z one-sample t-tests on mouse-level Pearson correlations

project_root = fileparts(fileparts(mfilename('fullpath')));
output_folder = fullfile(project_root, 'pearson_corr');
by_mouse_csv = fullfile(output_folder, 'hr_amp_corr_pearson_by_mouse.csv');
excel_file = fullfile(output_folder, 'hr_amp_corr_pearson_results.xlsx');
ttest_csv = fullfile(output_folder, 'hr_amp_corr_pearson_ttest_summary.csv');

if ~exist(by_mouse_csv, 'file')
    error('Missing %s. Run compute_hr_diam_corr_pearson.m first.', by_mouse_csv);
end

by_mouse_table = readtable(by_mouse_csv, 'TextType', 'string');
vessel_vars = unique(by_mouse_table.VesselVariable, 'stable');
window_names = ["Early", "Late", "Overall"];

ttest_table = table();

for v = 1:numel(vessel_vars)
    vessel_var = vessel_vars(v);
    rows = by_mouse_table.VesselVariable == vessel_var;

    for w = 1:numel(window_names)
        window_name = window_names(w);
        r_col = char("R_" + window_name);

        r_values = by_mouse_table.(r_col)(rows);
        valid = isfinite(r_values);
        r_values = r_values(valid);
        subjects_here = by_mouse_table.Subject(rows);
        subjects_here = subjects_here(valid); %#ok<NASGU>

        [p_ttest, t_stat, df, mean_z, ci_low_z, ci_high_z, n_subjects] = ...
            fisher_ttest_complete(r_values);

        mean_r_back = tanh(mean_z);
        ci_low_r_back = tanh(ci_low_z);
        ci_high_r_back = tanh(ci_high_z);
        median_r = median(r_values, 'omitnan');

        ttest_table = [ttest_table; table(vessel_var, window_name, n_subjects, ...
            mean(r_values, 'omitnan'), median_r, mean_z, ...
            mean_r_back, ci_low_r_back, ci_high_r_back, ...
            p_ttest, t_stat, df, ...
            'VariableNames', {'VesselVariable', 'Window', 'N_Subjects', ...
            'MeanR', 'MedianR', 'MeanFisherZ', ...
            'BackTransformedMeanR', 'CI_Low_R', 'CI_High_R', ...
            'P_TTest', 'TStatistic', 'DegreesOfFreedom'})];
    end
end

writetable(ttest_table, ttest_csv);
writetable(ttest_table, excel_file, 'Sheet', 'TTest');

fprintf('Saved %s\n', ttest_csv);
fprintf('Updated %s with TTest sheet\n', excel_file);

%% Local functions
function [p, t_stat, df, mean_z, ci_low_z, ci_high_z, n] = fisher_ttest_complete(r_values)
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
    mean_z = mean(z_values, 'omitnan');

    [~, p, ci, stats] = ttest(z_values, 0);
    t_stat = stats.tstat;
    df = stats.df;
    ci_low_z = ci(1);
    ci_high_z = ci(2);
end
