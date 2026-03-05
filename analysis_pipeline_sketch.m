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
rr_intervals = diff(t_rpeaks) * 1000;       % ms
t_rr = t_rpeaks(1:end-1) + diff(t_rpeaks)/2; % midpoint timestamps

% Step 2: Lowpass vessel diameter to remove cardiac pulsation
% Keep only fluctuations below ~2 Hz (vasomotion + respiratory range)

fc_slow = 0.3;  % Hz — change this to explore
[b_lp, a_lp] = butter(4, fc_slow / (fs_cam/2), 'low');
diam_slow = filtfilt(b_lp, a_lp, diam_norm);
filter_label = sprintf('Vessel LP < %.2f Hz', fc_slow);


% Step 3: Cross-correlate within each good section
max_lag_sec = 5;                             % ±5 seconds
min_section_dur = 3 * max_lag_sec;           % need at least 15 seconds

max_lag_samp = round(max_lag_sec * fs_cam);
lags_sec = (-max_lag_samp:max_lag_samp) / fs_cam;

all_xcorr = [];
section_weights = [];

for s = 1:n_sections
    t_start = good_sections(s).t_start;
    t_end   = good_sections(s).t_end;
    dur = t_end - t_start;

    % Skip short sections — need enough data for meaningful xcorr
    if dur < min_section_dur, continue; end

    % RR intervals in this section — interpolate to camera grid
    rr_mask = t_rr >= t_start & t_rr <= t_end;
    if sum(rr_mask) < 10, continue; end

    seg_mask = t_frames >= t_start & t_frames <= t_end;
    t_seg = t_frames(seg_mask);

    rr_interp = interp1(t_rr(rr_mask), rr_intervals(rr_mask), ...
                         t_seg, 'linear', 'extrap');

    % Vessel diameter (slow) for this section
    diam_seg = diam_slow(seg_mask);

    % Normalize both to zero mean, unit variance
    rr_z   = (rr_interp - mean(rr_interp)) / std(rr_interp);
    diam_z = (diam_seg - mean(diam_seg)) / std(diam_seg);

    % Cross-correlation (normalized)
    [xc, lags] = xcorr(rr_z, diam_z, max_lag_samp, 'coeff');

    all_xcorr = [all_xcorr; xc(:)'];
    section_weights = [section_weights; numel(diam_seg)];
end

figure('Position', [100 100 900 400]);
durs = [good_sections.t_end] - [good_sections.t_start];
fprintf('Section durations: min %.1f, median %.1f, max %.1f s\n', ...
         min(durs), median(durs), max(durs));
histogram(durs, 30);
xlabel('Duration (s)');
ylabel('Count');

% Weighted average across sections (longer sections contribute more)
w = section_weights / sum(section_weights);
xcorr_mean = w' * all_xcorr;

fprintf('Cross-correlation computed across %d sections\n', size(all_xcorr, 1));

%% Step 4: Plot
figure('Position', [100 100 900 400]);

subplot(1,2,1);
plot(lags_sec, xcorr_mean, 'b', 'LineWidth', 1.5);
hold on;
xline(0, 'r--');
yline(0, 'k:');
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
xline(0, 'r--');
yline(0, 'k:');
ylim([-1 1]);
xlabel('Lag (s)  [positive = RR leads vessel]');
ylabel('Cross-correlation');
title(sprintf('Individual sections + weighted mean | %s', filter_label));
hold off;

