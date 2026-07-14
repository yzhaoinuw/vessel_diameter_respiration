%% Compute HR-vessel correlations from fine-resolution hr_amp_data.mat files

data_folder = ".\data\";
results_folder = ".\results\";
if ~exist(results_folder, "dir")
    mkdir(results_folder);
end

hr_var = "hr_median_bpm";
vessel_vars = ["amp_svm", "bp_svm", "amp_c", "diam_mean"];
cutoff_min = 60;
dubious_exclusion_threshold = 0.50;

subjects = dir(data_folder);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, "."));

excluded_table = table();
results_table = table();

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

    hr = data.(hr_var);
    timepts = data.timepts;

    for v = 1:numel(vessel_vars)
        vessel_var = vessel_vars(v);
        vessel = data.(vessel_var);

        early = timepts < cutoff_min;
        late = timepts >= cutoff_min;

        [r_early, p_early, n_early] = corr_complete(hr(early), vessel(early));
        [r_late, p_late, n_late] = corr_complete(hr(late), vessel(late));
        [r_overall, p_overall, n_overall] = corr_complete(hr, vessel);

        results_table = [results_table; table(subject_name, hr_var, vessel_var, ...
            r_early, p_early, n_early, r_late, p_late, n_late, ...
            r_overall, p_overall, n_overall, dubious_fraction, ...
            'VariableNames', {'Subject', 'HrVariable', 'VesselVariable', ...
            'R_Early', 'P_Early', 'N_Early', ...
            'R_Late', 'P_Late', 'N_Late', ...
            'R_Overall', 'P_Overall', 'N_Overall', ...
            'DubiousHrBinFraction'})];
    end
end

fprintf("\nExcluded subjects due to > %.1f%% dubious HR bins:\n", ...
    100 * dubious_exclusion_threshold);
disp(excluded_table);

fprintf("\nCorrelation results:\n");
disp(results_table);

writetable(results_table, fullfile(results_folder, "hr_amp_corr_summary.csv"));
writetable(excluded_table, fullfile(results_folder, "hr_amp_corr_excluded_subjects.csv"));

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
    [r, p] = corr(x(valid), y(valid), "Type", "Spearman");
end
