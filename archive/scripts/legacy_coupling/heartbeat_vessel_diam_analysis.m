signal_data = load(".\data\signals.mat");
img_data = load(".\data\im.mat");

t_ecg = double(signal_data.t_ECG(:));
resp  = double(signal_data.resp(:));

t_frames = double(img_data.t_frames(:));
fs_cam = 1 / median(diff(t_frames));
fs_ecg = 1 / median(diff(t_ecg));
t_ecgoffset = double(img_data.t_ecgoffset);
assert(t_ecgoffset > 0, 't_ecgoffset must be positive.')
% Align phys timestamps to camera clock
t_phys = t_ecg - t_ecgoffset;

% ecg
good_r_peaks = signal_data.good_r_peaks; % good R peaks: conf >= 0.5
t_rpeaks = t_phys(good_r_peaks);
t_rpeaks = t_rpeaks(t_rpeaks>0 & t_rpeaks<t_frames(end));
good_sections = find_good_r_peak_sections(t_rpeaks);

diam_norm = preprocess_vessel_diam(img_data);
%% Extract matched segments of ECG, respiration, and vessel diameter
n_sections = numel(good_sections);

%% Resample respiration to camera rate
% t_frames is the authoritative camera time axis
resp_resampled = interp1(t_phys, resp, t_frames, 'linear', 'extrap');

%% Extract matched segments
for s = 1:n_sections

    t_start = good_sections(s).t_start;
    t_end   = good_sections(s).t_end;

    % Single mask — both signals now share t_frames
    mask = t_frames >= t_start & t_frames <= t_end;

    good_sections(s).resp = resp_resampled(mask);
    good_sections(s).diam = diam_norm(mask);

    % Sanity check — both should give identical sample counts
    n_resp = numel(good_sections(s).resp);
    n_diam = numel(good_sections(s).diam);
    dur    = n_resp / fs_cam;
    fprintf('Section %d:  %.2f s  |  resp %d samples  |  diam %d samples\n', ...
            s, dur, n_resp, n_diam);

    if n_resp ~= n_diam
        fprintf('  WARNING: sample count mismatch — should never happen\n');
    end
end

%% ============================================================
%% Phase 3a: Event-Triggered Averaging (ETA)
%% ============================================================

%% Vectorized ETA: R-peak → Vessel Diameter

win_sec  = 0.150;
win_samp = round(win_sec * fs_cam);
t_win    = (-win_samp:win_samp) / fs_cam * 1000;  % ms

% Convert all R-peak times to nearest camera frame indices
rp_idx = interp1(t_frames, 1:numel(t_frames), t_rpeaks, 'nearest');
rp_idx = rp_idx(~isnan(rp_idx));  % drop any outside t_frames range

% Remove edge cases where window would exceed signal bounds
valid = rp_idx - win_samp >= 1 & rp_idx + win_samp <= numel(diam_norm);
rp_idx = rp_idx(valid);

% Build index matrix: each row is one snippet
offsets = -win_samp:win_samp;                      % 1 x (2*win_samp+1)
idx_matrix = round(rp_idx(:)) + offsets;            % n_beats x (2*win_samp+1)

% Extract all snippets in one shot
all_snippets = diam_norm(idx_matrix);               % n_beats x (2*win_samp+1)

n_beats  = size(all_snippets, 1);
eta_mean = mean(all_snippets, 1);
eta_sem  = std(all_snippets, 0, 1) / sqrt(n_beats);

fprintf('ETA: %d beats collected\n', n_beats);
%% Plot ETA + overlay a handful of individual snippets for intuition
figure('Position', [100 100 900 400]);

% Left: individual traces + mean
subplot(1,2,1);
hold on;
n_show = min(50, n_beats);
idx_show = randperm(n_beats, n_show);
for k = 1:n_show
    plot(t_win, all_snippets(idx_show(k),:), 'Color', [0.7 0.7 0.7 0.3]);
end
plot(t_win, eta_mean, 'b', 'LineWidth', 2);
xline(0, 'r--', 'R-peak');
xlabel('Time from R-peak (ms)');
ylabel('Vessel diameter (norm)');
title(sprintf('Individual traces (n=%d shown) + mean', n_show));
hold off;

% Right: mean ± SEM
subplot(1,2,2);
hold on;
fill([t_win fliplr(t_win)], ...
     [eta_mean+eta_sem fliplr(eta_mean-eta_sem)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_win, eta_mean, 'b', 'LineWidth', 2);
xline(0, 'r--', 'R-peak');
xlabel('Time from R-peak (ms)');
ylabel('Vessel diameter (% from baseline)');
title(sprintf('Mean \\pm SEM (n=%d beats)', n_beats));
hold off;

%% Bootstrap CI on R-peak → Vessel Diameter ETA
n_boot = 1000;
boot_means = zeros(n_boot, size(all_snippets, 2));
for i = 1:n_boot
    idx = randi(n_beats, n_beats, 1);
    boot_means(i,:) = mean(all_snippets(idx,:), 1);
end
boot_upper = prctile(boot_means, 97.5, 1);
boot_lower = prctile(boot_means, 2.5, 1);

figure('Position', [100 100 600 400]);
fill([t_win fliplr(t_win)], [boot_upper fliplr(boot_lower)], ...
     [0.7 0.85 1], 'EdgeColor', 'none', 'DisplayName', '95% bootstrap CI');
hold on;
plot(t_win, eta_mean, 'b', 'LineWidth', 2, 'DisplayName', 'Mean');
xline(0, 'r--', 'R-peak');
xlabel('Time from R-peak (ms)');
ylabel('Vessel diameter (% from baseline)');
title(sprintf('R-peak \\rightarrow Vessel ETA with bootstrap 95%% CI (n=%d beats)', n_beats));
legend('Location', 'best');
hold off;

%% Phase 3b: Cross-correlation — RR intervals ↔ Vessel Diameter
%% ============================================================
% Step 1: Compute RR interval series
rr_intervals = diff(t_rpeaks) * 1000;
t_rr = t_rpeaks(1:end-1) + diff(t_rpeaks)/2;

% Step 2: Lowpass vessel diameter to remove cardiac pulsation
% Keep only fluctuations below ~2 Hz (vasomotion + respiratory range)

fc_slow = 2;  % Hz — change this to explore
[b_lp, a_lp] = butter(4, fc_slow / (fs_cam/2), 'low');
diam_slow = filtfilt(b_lp, a_lp, diam_norm);
filter_label = sprintf('Vessel LP < %.2f Hz', fc_slow);

% Step 3: Cross-correlate within each good section
max_lag_sec = 5;
min_section_dur = 3 * max_lag_sec;
max_lag_samp = round(max_lag_sec * fs_cam);
lags_sec = (-max_lag_samp:max_lag_samp) / fs_cam;

all_xcorr      = [];
section_weights = [];
rr_z_store     = {};   % store for shuffle test
diam_z_store   = {};

for s = 1:n_sections
    t_start = good_sections(s).t_start;
    t_end   = good_sections(s).t_end;
    if (t_end - t_start) < min_section_dur, continue; end

    rr_mask = t_rr >= t_start & t_rr <= t_end;
    if sum(rr_mask) < 10, continue; end

    seg_mask = t_frames >= t_start & t_frames <= t_end;
    t_seg    = t_frames(seg_mask);

    rr_interp = interp1(t_rr(rr_mask), rr_intervals(rr_mask), t_seg, 'linear', 'extrap');
    diam_seg  = diam_slow(seg_mask);

    rr_z   = (rr_interp - mean(rr_interp)) / std(rr_interp);
    diam_z = (diam_seg  - mean(diam_seg))  / std(diam_seg);

    [xc, ~] = xcorr(rr_z, diam_z, max_lag_samp, 'coeff');
    all_xcorr       = [all_xcorr; xc(:)'];
    section_weights = [section_weights; numel(diam_seg)];
    rr_z_store{end+1}   = rr_z;
    diam_z_store{end+1} = diam_z;
end

w          = section_weights / sum(section_weights);
xcorr_mean = w' * all_xcorr;
[~, peak_idx] = max(abs(xcorr_mean));
peak_lag  = lags_sec(peak_idx);
peak_corr = xcorr_mean(peak_idx);
fprintf('Cross-correlation computed across %d sections\n', size(all_xcorr, 1));

% Step 4: Phase randomization null distribution (200 shuffles)
N_SHUF     = 200;
n_sec      = numel(rr_z_store);
null_peaks = zeros(N_SHUF, 1);

for sh = 1:N_SHUF
    shuf_xcorr = zeros(n_sec, numel(lags_sec));
    for s = 1:n_sec
        rr_z   = rr_z_store{s};
        diam_z = diam_z_store{s};
        % Circular shift by random offset — at least 1 s worth of samples
        min_shift = round(fs_cam);
        shift = randi([min_shift, numel(rr_z) - min_shift]);
        rr_z_shuf = circshift(rr_z, shift);
        [xc, ~] = xcorr(rr_z_shuf, diam_z, max_lag_samp, 'coeff');
        shuf_xcorr(s, :) = xc(:)';
    end
    shuf_mean = w' * shuf_xcorr;
    null_peaks(sh) = max(abs(shuf_mean));
end

null_p95  = prctile(null_peaks, 95);
null_p99  = prctile(null_peaks, 99);
is_sig    = abs(peak_corr) > null_p95;
fprintf('Observed peak: %.3f at %.2f s\n', peak_corr, peak_lag);
fprintf('Null 95th pct: %.3f | 99th pct: %.3f | Significant: %s\n', ...
        null_p95, null_p99, string(is_sig));

%%
% Step 5: Plot
figure('Position', [100 100 900 400]);

subplot(1,2,1);
plot(lags_sec, xcorr_mean, 'b', 'LineWidth', 1.5);
hold on;
yline(null_p95,  'r--', '95%', 'LineWidth', 1);
yline(-null_p95, 'r--',        'LineWidth', 1);
yline(0, 'k:');
yline( null_p95, 'r--', sprintf(' %.3f', null_p95),  'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
yline(-null_p95, 'r--', sprintf('%.3f', -null_p95), 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
plot(peak_lag, peak_corr, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
text(peak_lag, peak_corr, sprintf('  (%.2f s, %.3f)', peak_lag, peak_corr), ...
    'Color', 'r', 'FontSize', 9, 'VerticalAlignment', 'bottom');
ylim([-1 1]);
xlabel('Lag (s)  [positive = RR leads vessel]');
ylabel('Cross-correlation');
title(sprintf('RR intervals \\leftrightarrow %s', filter_label));
hold off;

subplot(1,2,2);
hold on;
for k = 1:size(all_xcorr, 1)
    plot(lags_sec, all_xcorr(k,:), 'Color', [0.6 0.6 0.6 0.4]);
end
plot(lags_sec, xcorr_mean, 'b', 'LineWidth', 2);
yline(0, 'k:');

ylim([-1 1]);
xlabel('Lag (s)  [positive = RR leads vessel]');
ylabel('Cross-correlation');
title(sprintf('Individual sections + weighted mean | %s', filter_label));
hold off;