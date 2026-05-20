%% Compute Spearman correlations for filtered 1-minute HR vs overlapping vessel features

project_root = fileparts(fileparts(mfilename('fullpath')));
data_folder = fullfile(project_root, 'data');
output_folder = fullfile(project_root, 'spearman_correlation_overlap');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

hr_var = "hr_minute_bpm_filtered";
vessel_vars = ["amp_svm", "bp_svm", "amp_c"];
window_names = ["Early", "Late", "Overall"];
cutoff_min = 60;
hr_nan_exclusion_threshold = 0.50;
manual_excluded_mouse_ids = ["F168", "F169", "M163", "M166", "M167"];

excel_file = fullfile(output_folder, 'hr_amp_corr_overlap_spearman_results.xlsx');
by_mouse_csv = fullfile(output_folder, 'hr_amp_corr_overlap_spearman_by_mouse.csv');
wilcoxon_csv = fullfile(output_folder, 'hr_amp_corr_overlap_spearman_wilcoxon_summary.csv');
excluded_csv = fullfile(output_folder, 'hr_amp_corr_overlap_spearman_excluded_subjects.csv');

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, '.'));

excluded_table = table();
by_mouse_table = table();

for i = 1:length(subjects)
    subject_name = string(subjects(i).name);
    subject_folder = fullfile(data_folder, subject_name);
    hr_file = fullfile(subject_folder, 'hr_amp_data.mat');
    vessel_file = fullfile(subject_folder, 'vessel_features_10min_overlap.mat');

    if ~exist(hr_file, 'file') || ~exist(vessel_file, 'file')
        fprintf("Skipping %s - missing hr or vessel overlap file\n", subject_name);
        continue
    end

    excluded_mouse_id = matching_mouse_id(subject_name, manual_excluded_mouse_ids);
    if strlength(excluded_mouse_id) > 0
        fprintf("Excluding %s - manual exclusion for HR disagreement (%s)\n", ...
            subject_name, excluded_mouse_id);
        excluded_table = [excluded_table; table(subject_name, NaN, ...
            NaN, NaN, "manual HR disagreement exclusion", ...
            'VariableNames', {'Subject', 'FilteredHrNaFraction', ...
            'NFiniteAlignedHrBins', 'NAlignedBins', 'Reason'})];
        continue
    end

    hr_data = load(hr_file, 'hr_filter_timepts_min', 'hr_minute_bpm_filtered');
    vessel_data = load(vessel_file, 'timepts', 'bin_centers_min', ...
        'amp_svm', 'bp_svm', 'amp_c');

    [timepts, hr, matched] = align_hr_to_vessel_starts( ...
        vessel_data.timepts, hr_data.hr_filter_timepts_min, ...
        hr_data.hr_minute_bpm_filtered);
    hr_nan_fraction = mean(~isfinite(hr(matched)));

    if hr_nan_fraction > hr_nan_exclusion_threshold
        fprintf("Excluding %s - filtered HR NaN bins %.1f%% > %.1f%%\n", ...
            subject_name, 100 * hr_nan_fraction, ...
            100 * hr_nan_exclusion_threshold);
        excluded_table = [excluded_table; table(subject_name, hr_nan_fraction, ...
            nnz(isfinite(hr(matched))), nnz(matched), ...
            "filtered HR NaN fraction > threshold", ...
            'VariableNames', {'Subject', 'FilteredHrNaFraction', ...
            'NFiniteAlignedHrBins', 'NAlignedBins', 'Reason'})];
        continue
    end

    early = timepts < cutoff_min;
    late = timepts >= cutoff_min;

    for v = 1:numel(vessel_vars)
        vessel_var = vessel_vars(v);
        vessel = vessel_data.(char(vessel_var));

        [r_early, p_early, n_early] = corr_complete(hr(early), vessel(early));
        [r_late, p_late, n_late] = corr_complete(hr(late), vessel(late));
        [r_overall, p_overall, n_overall] = corr_complete(hr, vessel);

        by_mouse_table = [by_mouse_table; table(subject_name, hr_var, vessel_var, ...
            r_early, p_early, n_early, r_late, p_late, n_late, ...
            r_overall, p_overall, n_overall, hr_nan_fraction, ...
            'VariableNames', {'Subject', 'HrVariable', 'VesselVariable', ...
            'R_Early', 'P_Early', 'N_Early', ...
            'R_Late', 'P_Late', 'N_Late', ...
            'R_Overall', 'P_Overall', 'N_Overall', ...
            'FilteredHrNaFraction'})];
    end
end

wilcoxon_table = summarize_mouse_level_tests(by_mouse_table, vessel_vars, window_names);

writetable(by_mouse_table, by_mouse_csv);
writetable(wilcoxon_table, wilcoxon_csv);
writetable(excluded_table, excluded_csv);

writetable(by_mouse_table, excel_file, 'Sheet', 'ByMouse');
writetable(wilcoxon_table, excel_file, 'Sheet', 'Wilcoxon');
writetable(excluded_table, excel_file, 'Sheet', 'ExcludedSubjects');

fprintf("\nSaved overlap Spearman outputs to:\n%s\n", output_folder);
fprintf("- %s\n", excel_file);
fprintf("- %s\n", by_mouse_csv);
fprintf("- %s\n", wilcoxon_csv);
fprintf("- %s\n", excluded_csv);

%% Local functions
function mouse_id = matching_mouse_id(subject_name, mouse_ids)
    subject_name = string(subject_name);
    mouse_id = "";
    for i = 1:numel(mouse_ids)
        candidate = string(mouse_ids(i));
        if startsWith(subject_name, candidate, "IgnoreCase", true)
            mouse_id = candidate;
            return
        end
    end
end

function [timepts, hr, matched] = align_hr_to_vessel_starts(vessel_timepts, hr_timepts, hr_values)
    timepts = double(vessel_timepts(:));
    hr_timepts = double(hr_timepts(:));
    hr_values = double(hr_values(:));

    [matched, loc] = ismember(timepts, hr_timepts);
    hr = nan(size(timepts));
    hr(matched) = hr_values(loc(matched));
end

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
    [r, p] = corr(x(valid), y(valid), 'Type', 'Spearman');
end

function summary_table = summarize_mouse_level_tests(by_mouse_table, vessel_vars, window_names)
    summary_table = table();
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

            summary_table = [summary_table; table(vessel_var, window_name, ...
                n_subjects, mean(r_values, 'omitnan'), median(r_values, 'omitnan'), ...
                min(r_values, [], 'omitnan'), max(r_values, [], 'omitnan'), ...
                p_wilcoxon, signed_rank, mean(mouse_p_values, 'omitnan'), ...
                median(mouse_p_values, 'omitnan'), ...
                'VariableNames', {'VesselVariable', 'Window', ...
                'N_Subjects', 'MeanR', 'MedianR', 'MinR', 'MaxR', ...
                'P_Wilcoxon', 'SignedRankStatistic', ...
                'MeanMousePValue', 'MedianMousePValue'})];
        end
    end
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
