%% ========================================================================
%  TEST 1 - Respiration-Triggered Vessel Diameter Averaging
%  Signals needed: respiration + vessel diameter (no ECG required)
%  Question: does vessel diameter show a stereotyped response to each breath?
%  ========================================================================

clear; clc; close all;

%% 0. Data loading
signal_data = load(".\data\signals.mat");
img_data    = load(".\data\im.mat");

t_ecg    = double(signal_data.t_ECG(:));
resp_raw = double(signal_data.resp(:));
t_frames = double(img_data.t_frames(:));

fs_cam = 1 / median(diff(t_frames));
fs_ecg = 1 / median(diff(t_ecg));

t_ecgoffset = double(img_data.t_ecgoffset);
assert(t_ecgoffset > 0, 't_ecgoffset must be positive.')
t_phys = t_ecg - t_ecgoffset;

%% 1. Vessel diameter preprocessing
diam_norm = preprocess_vessel_diam(img_data);

%% 2. Respiration preprocessing
% 2a. Lowpass at 5 Hz to remove cardiac contamination
%     Mouse breathing: ~1-4 Hz.  Mouse cardiac: ~8-12 Hz.
Fc_resp = 5;
[b_lp, a_lp] = butter(4, Fc_resp / (fs_ecg/2), 'low');
resp_filt = filtfilt(b_lp, a_lp, resp_raw);

% 2b. Sliding-window std to find good-amplitude breathing sections
%     2-second window -> captures 4-8 breath cycles at 2-4 Hz breathing,
%     giving a stable amplitude estimate.
win_std_sec  = 2;
win_std_samp = round(win_std_sec * fs_ecg);
resp_std     = movstd(resp_filt, win_std_samp);

% Visualise - use this to tune the threshold
figure('Name', 'Resp std distribution');
histogram(resp_std, 100, 'FaceColor', [0.27 0.51 0.71], 'EdgeColor', 'none');
xlabel('Sliding-window std (a.u.)'); ylabel('Count');
title('Respiration amplitude distribution - choose threshold');

% Start with 25th percentile; adjust after inspecting the histogram
std_thresh = prctile(resp_std, 30);
xline(std_thresh, '--r', 'LineWidth', 1.5, ...
      'Label', sprintf('Threshold = %.3f', std_thresh));
fprintf('Resp std threshold (25th pctile): %.4f\n', std_thresh);

resp_good_mask = resp_std > std_thresh;

% 2c. Detect breath peaks in filtered signal
%     MinPeakProminence rejects the tiny wiggles between real breaths.
%     Start with median prominence / 2 as threshold; tune after visual check.
min_peak_dist = round(0.20 * fs_ecg);          % 200 ms (~5 Hz ceiling)
resp_prom     = median(abs(resp_filt)) * 1;   % data-driven prominence threshold

[~, pk_locs] = findpeaks( resp_filt, 'MinPeakDistance', min_peak_dist, ...
                          'MinPeakProminence', resp_prom);
[~, tr_locs] = findpeaks(-resp_filt, 'MinPeakDistance', min_peak_dist, ...
                          'MinPeakProminence', resp_prom);

% Keep only peaks/troughs inside good-std regions
pk_good = pk_locs(resp_good_mask(pk_locs));
tr_good = tr_locs(resp_good_mask(tr_locs));

fprintf('Breath peaks: %d total, %d in good regions\n', ...
        numel(pk_locs), numel(pk_good));

% 2d. Sanity check - breath rate distribution
breath_rate_Hz = 1 ./ diff(t_phys(pk_good));

figure('Name', 'Breath rate distribution');
histogram(breath_rate_Hz, 50, 'FaceColor', [0.47 0.67 0.19], 'EdgeColor', 'none');
xlabel('Breath rate (Hz)'); ylabel('Count');
title('Instantaneous breathing rate (peak-to-peak)');
fprintf('Breath rate: median = %.1f Hz, IQR = [%.1f, %.1f] Hz\n', ...
        median(breath_rate_Hz), prctile(breath_rate_Hz, 25), prctile(breath_rate_Hz, 75));

% 2e. Visual check - plot from camera start (where diam data exists)
figure('Name', 'Resp peak detection check');
t_plot_start = t_frames(1);        % start of vessel data
t_plot_end   = t_plot_start + 10;  % 10 seconds
idx_show = t_phys >= t_plot_start & t_phys <= t_plot_end;
plot(t_phys(idx_show), resp_filt(idx_show), 'k', 'LineWidth', 0.8); hold on;
pk_show = pk_good(t_phys(pk_good) >= t_plot_start & t_phys(pk_good) <= t_plot_end);
tr_show = tr_good(t_phys(tr_good) >= t_plot_start & t_phys(tr_good) <= t_plot_end);
plot(t_phys(pk_show), resp_filt(pk_show), 'rv', 'MarkerFaceColor', 'r');
plot(t_phys(tr_show), resp_filt(tr_show), 'b^', 'MarkerFaceColor', 'b');
xlabel('Time (s)'); ylabel('Resp (a.u.)');
title('Filtered respiration - 10 s from camera start');
legend('Filtered resp', 'Insp. peak', 'Exp. trough');
hold off;

%% 3. Identify good resp sections and filter to >=5 s
%     We only need this to ensure breath peaks come from sustained clean
%     breathing, not isolated good samples surrounded by noise.
min_section_sec = 5;
resp_good_runs  = bwconncomp(resp_good_mask);

% Build a mask of samples that belong to sufficiently long good runs
resp_long_mask = false(size(resp_good_mask));
for r = 1:resp_good_runs.NumObjects
    idx = resp_good_runs.PixelIdxList{r};
    dur = (idx(end) - idx(1)) / fs_ecg;
    if dur >= min_section_sec
        resp_long_mask(idx) = true;
    end
end

% Final set of breath peaks: in good-std AND in a >=5 s run
pk_final = pk_good(resp_long_mask(pk_good));
fprintf('Breath peaks in good sections (>=%d s): %d\n', ...
        min_section_sec, numel(pk_final));

%% 4. Vectorised resp-triggered ETA (same pattern as R-peak -> diam)

% Convert breath peak times to nearest camera frame indices
t_pk_final = t_phys(pk_final);
pk_frame_idx = interp1(t_frames, 1:numel(t_frames), t_pk_final, 'nearest');
pk_frame_idx = pk_frame_idx(~isnan(pk_frame_idx));

% ETA window: +/-500 ms (covers one full breath cycle at ~2 Hz)
win_sec  = 0.500;
win_samp = round(win_sec * fs_cam);
t_win    = (-win_samp:win_samp) / fs_cam * 1000;  % ms

% Remove edge cases
valid = pk_frame_idx - win_samp >= 1 & pk_frame_idx + win_samp <= numel(diam_norm);
pk_frame_idx = pk_frame_idx(valid);

% Build index matrix and extract all snippets in one shot
offsets      = -win_samp:win_samp;
idx_matrix   = round(pk_frame_idx(:)) + offsets;
all_snippets = diam_norm(idx_matrix);

n_triggers = size(all_snippets, 1);
eta_mean   = mean(all_snippets, 1);
eta_sem    = std(all_snippets, 0, 1) / sqrt(n_triggers);

fprintf('Resp-triggered ETA: %d breath triggers used\n', n_triggers);

%% 5. Plot
figure('Position', [100 100 900 400]);

% Left: individual traces + mean
subplot(1,2,1);
hold on;
n_show = min(50, n_triggers);
idx_show = randperm(n_triggers, n_show);
for k = 1:n_show
    plot(t_win, all_snippets(idx_show(k),:), 'Color', [0.7 0.7 0.7 0.3]);
end
plot(t_win, eta_mean, 'b', 'LineWidth', 2);
xline(0, 'r--', 'Insp. peak');
xlabel('Time from inspiration peak (ms)');
ylabel('Vessel diameter (% from baseline)');
title(sprintf('Individual traces (n=%d shown) + mean', n_show));
hold off;

% Right: mean +/- SEM
subplot(1,2,2);
hold on;
fill([t_win fliplr(t_win)], ...
     [eta_mean+eta_sem fliplr(eta_mean-eta_sem)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_win, eta_mean, 'b', 'LineWidth', 2);
xline(0, 'r--', 'Insp. peak');
xlabel('Time from inspiration peak (ms)');
ylabel('Vessel diameter (% from baseline)');
title(sprintf('Mean \\pm SEM (n=%d breaths)', n_triggers));
hold off;

%% Bootstrap CI on the ETA
n_boot = 1000;
boot_means = zeros(n_boot, size(all_snippets, 2));

for i = 1:n_boot
    % Resample with replacement from actual triggered snippets
    idx = randi(n_triggers, n_triggers, 1);
    boot_means(i,:) = mean(all_snippets(idx,:), 1);
end

boot_upper = prctile(boot_means, 97.5, 1);
boot_lower = prctile(boot_means, 2.5, 1);

figure('Name', 'Bootstrap CI');
fill([t_win fliplr(t_win)], [boot_upper fliplr(boot_lower)], ...
     [0.7 0.85 1], 'EdgeColor', 'none'); hold on;
plot(t_win, eta_mean, 'b', 'LineWidth', 2);
xline(0, 'r--', 'Insp. peak');
xlabel('Time from inspiration peak (ms)');
ylabel('Vessel diameter (% from baseline)');
title(sprintf('ETA with bootstrap 95%% CI  (n = %d breaths)', n_triggers));
hold off;