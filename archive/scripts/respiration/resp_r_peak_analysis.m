%% TEST 02 — Respiration → Heart Rate ETA
%  Trigger on inspiration peaks, average instantaneous heart rate around
%  each trigger.  Tests for respiratory sinus arrhythmia: does the mouse's
%  heart speed up during inspiration and slow down during expiration?
%
%  Signals used : respiration + ECG (no vessel diameter)
%  Good-section catalog : resp ∩ ECG
%  Established conventions:
%    - All times in camera clock (t_phys = t_ECG − t_ecgoffset)
%    - Beat-indexed ETA on native R-peak times, plotted vs mean time offset
%    - 2-second sliding-window std for resp amplitude
%    - MinPeakProminence = 0.5 × median |resp_filt| alongside MinPeakDistance
%    - Diagnostic time axis anchored to t_frames(1)

%% ========================================================================
%  1. DATA LOADING
%  ========================================================================
signal_data = load('.\data\signals.mat');
img_data    = load('.\data\im.mat');          % needed for t_ecgoffset, t_frames

% --- physiology signals (1000 Hz) ---
t_ECG  = double(signal_data.t_ECG(:));
resp   = double(signal_data.resp(:));
fs_phys = 1 / median(diff(t_ECG));            % ≈ 1000 Hz

% --- R-peak indices (from CNN detector, conf >= 0.5) ---
good_r_peaks = signal_data.good_r_peaks;

% --- camera timing (for clock reference only) ---
t_frames    = double(img_data.t_frames(:));
t_ecgoffset = double(img_data.t_ecgoffset);

%% ========================================================================
%  2. TIME ALIGNMENT  (everything → camera clock)
%  ========================================================================
t_phys = t_ECG - t_ecgoffset;                 % physiology times in camera clock

min_dur_sec = 5;                               % minimum section length (used throughout)

%% ========================================================================
%  3. RESPIRATION PREPROCESSING  (via preprocess_resp)
%  ========================================================================
rr = preprocess_resp(resp, fs_phys, min_dur_sec, [0.5 5], 2, false, true);

resp_filt = rr.resp_filt;

% Convert index-based good_sections to time-based (camera clock)
resp_sections = table( ...
    t_phys(rr.good_sections.idx_start), ...
    t_phys(rr.good_sections.idx_end), ...
    rr.good_sections.duration, ...
    'VariableNames', {'t_start','t_end','duration'});

% Inspiration peak times & amplitudes (camera clock)
insp_peak_times = t_phys(rr.insp_peaks.locs);
insp_peak_amps  = rr.insp_peaks.amps;

%% ========================================================================
%  4. ECG GOOD SECTIONS  (via find_good_r_peak_sections)
%  ========================================================================
%  Derive R-peak times from good_r_peaks indices, clip to valid camera range
t_rpeaks = t_phys(good_r_peaks);
t_rpeaks = t_rpeaks(t_rpeaks > 0 & t_rpeaks < t_frames(end));

ecg_sec_struct = find_good_r_peak_sections(t_rpeaks, min_dur_sec);

% Convert struct array → table to match resp_sections format
ecg_sections = table( ...
    [ecg_sec_struct.t_start]', ...
    [ecg_sec_struct.t_end]', ...
    [ecg_sec_struct.duration]', ...
    'VariableNames', {'t_start','t_end','duration'});

fprintf('ECG good sections:   %d  (%.1f s total)\n', ...
    height(ecg_sections), sum(ecg_sections.duration));

%% ========================================================================
%  5. INTERSECT  resp ∩ ECG  good sections
%  ========================================================================
good_sections = intersect_sections(resp_sections, ecg_sections, min_dur_sec);
fprintf('Resp∩ECG sections:   %d  (%.1f s total)\n', ...
    height(good_sections), sum(good_sections.duration));

%% ========================================================================
%  6–8. PER-SECTION: grab ±N beats around each inspiration peak, collect ETA
%  ========================================================================
%  For each intersection section:
%    1) Get validated R-peak times that fall inside it
%    2) Filter precomputed inspiration peaks to this section
%    3) For each trigger, find nearest R-peak, grab ±N beats
%    4) Compute HR = 60 / diff(peak_times)
%  Everything stays within one section — no boundary issues.

n_beats_half = 5;
beat_offsets = -n_beats_half : n_beats_half;
n_beat_win   = numel(beat_offsets);

% Collect validated R-peak times from ECG sections (consecutive intervals
%  already guaranteed plausible by find_good_r_peak_sections)
t_rpeaks_valid = [];
for s = 1:numel(ecg_sec_struct)
    t_rpeaks_valid = [t_rpeaks_valid; ecg_sec_struct(s).peak_time(:)]; %#ok<AGROW>
end
t_rpeaks_valid = sort(t_rpeaks_valid);

hr_snippets     = [];                          % n_events × n_beat_win
time_off_snippets = [];                        % same size, actual time offsets
resp_amp_all    = [];                          % breath amplitude per trigger
insp_times_all  = [];                          % trigger times (for diagnostics)

for s = 1:height(good_sections)
    ts = good_sections.t_start(s);
    te = good_sections.t_end(s);

    % --- R-peaks in this section (from validated ECG sections only) ---
    rp = t_rpeaks_valid(t_rpeaks_valid >= ts & t_rpeaks_valid <= te);
    if numel(rp) < 2*n_beats_half + 2, continue; end  % need enough beats

    % --- Inspiration peaks in this section (precomputed) ---
    in_sec = insp_peak_times >= ts & insp_peak_times <= te;
    trig_times = insp_peak_times(in_sec);
    trig_amps  = insp_peak_amps(in_sec);
    if isempty(trig_times), continue; end

    % --- For each trigger, find nearest R-peak, grab ±N beats ---
    n_rp_sec = numel(rp);
    trig_rp_idx = interp1(rp, 1:n_rp_sec, trig_times, 'nearest');
    trig_rp_idx = round(trig_rp_idx);

    % Keep only triggers with ±N+1 beats available (need +1 for last HR)
    ok = trig_rp_idx - n_beats_half >= 1 & ...
         trig_rp_idx + n_beats_half + 1 <= n_rp_sec & ...
         ~isnan(trig_rp_idx);
    trig_rp_idx = trig_rp_idx(ok);
    trig_times  = trig_times(ok);
    trig_amps   = trig_amps(ok);
    if isempty(trig_rp_idx), continue; end

    % Index matrix into this section's R-peaks
    idx_mat = trig_rp_idx(:) + beat_offsets;   % n_trig × n_beat_win

    % HR snippets: 60 / (rp(i+1) − rp(i))
    snip = 60 ./ (rp(idx_mat + 1) - rp(idx_mat));

    % Time offsets relative to trigger's nearest R-peak
    t_off = rp(idx_mat) - rp(trig_rp_idx(:));

    hr_snippets     = [hr_snippets;     snip];      %#ok<AGROW>
    time_off_snippets = [time_off_snippets; t_off];  %#ok<AGROW>
    resp_amp_all    = [resp_amp_all;    trig_amps];  %#ok<AGROW>
    insp_times_all  = [insp_times_all;  trig_times]; %#ok<AGROW>
end

n_events = size(hr_snippets, 1);
mean_time_offsets = mean(time_off_snippets, 1);

eta_mean = mean(hr_snippets, 1);
eta_sem  = std(hr_snippets, 0, 1) / sqrt(n_events);

fprintf('ETA: %d triggers, %d beats per snippet\n', n_events, n_beat_win);
fprintf('Mean beat spacing ≈ %.1f ms\n', 1000*mean(diff(mean_time_offsets)));

%% ========================================================================
%  9. DIAGNOSTIC PLOTS
%  ========================================================================

% ---- Figure 1: example good section (resp + HR + peak markers) ----------
figure('Position', [100 100 900 350], 'Color', 'w');

% Pick the longest section for the example
[~, best] = max(good_sections.duration);
ex_start = good_sections.t_start(best);
ex_end   = good_sections.t_end(best);
% Show at most 10 seconds for readability
ex_end   = min(ex_end, ex_start + 10);

% Respiration trace (filtered, physiology grid)
mask_resp = t_phys >= ex_start & t_phys <= ex_end;
yyaxis left
plot(t_phys(mask_resp), resp_filt(mask_resp), 'Color', [0.2 0.5 0.8], ...
    'LineWidth', 0.8);
ylabel('Resp (filtered, a.u.)');

% HR trace — beat-by-beat values from R-peaks in this window
rp_in_win = t_rpeaks_valid(t_rpeaks_valid >= ex_start & t_rpeaks_valid <= ex_end);
rr_win = diff(rp_in_win);
hr_win = 60 ./ rr_win;
yyaxis right
plot(rp_in_win(1:end-1), hr_win, '.-', 'Color', [0.85 0.33 0.1], ...
    'MarkerSize', 8, 'LineWidth', 0.8);
ylabel('Heart rate (bpm)');

% Mark inspiration peaks in this window
peaks_in_win = insp_times_all(insp_times_all >= ex_start & insp_times_all <= ex_end);
yyaxis left; hold on;
for k = 1:numel(peaks_in_win)
    xline(peaks_in_win(k), '--', 'Color', [0.4 0.4 0.4], 'Alpha', 0.5);
end
hold off;
xlabel('Time (s)');
title(sprintf('Test 02: Example good section (%.1f s)', ex_end - ex_start));
xlim([ex_start ex_end]);

% ---- Figure 2: ETA — individual traces (left) and mean ± SEM (right) ---
figure('Position', [100 500 1100 350], 'Color', 'w');

t_eta = mean_time_offsets * 1000;

% Left: 50 random individual traces + mean
subplot(1,2,1);
trace_idx = randperm(n_events, min(50, n_events));
plot(t_eta, hr_snippets(trace_idx, :)', 'Color', [0.6 0.75 0.9 0.3], 'LineWidth', 0.5);
hold on;
plot(t_eta, eta_mean, '-', 'Color', [0.15 0.3 0.7], 'LineWidth', 2);
xline(0, '--k', 'LineWidth', 1);
hold off;
xlabel('Time from inspiration peak (ms)');
ylabel('Heart rate (bpm)');
title(sprintf('Individual traces (n=50 of %d)', n_events));
xlim([t_eta(1)*1.05  t_eta(end)*1.05]);

% Right: mean ± SEM
subplot(1,2,2);
fill([t_eta fliplr(t_eta)], ...
     [eta_mean + eta_sem, fliplr(eta_mean - eta_sem)], ...
     [0.8 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
plot(t_eta, eta_mean, '.-', 'Color', [0.15 0.3 0.7], ...
    'LineWidth', 1.8, 'MarkerSize', 14);
xline(0, '--k', 'LineWidth', 1);
hold off;
xlabel('Time from inspiration peak (ms)');
ylabel('Heart rate (bpm)');
title(sprintf('Mean ± SEM  (n = %d triggers)', n_events));
xlim([t_eta(1)*1.05  t_eta(end)*1.05]);

sgtitle('Test 02: Resp → HR ETA');

%%
% ---- Figure 3: Bootstrap 95% CI on ETA shape ----------------------------
figure('Position', [100 900 900 350], 'Color', 'w');

n_boot = 2000;
rng(42);
boot_means = zeros(n_boot, n_beat_win);
for b = 1:n_boot
    idx = randi(n_events, n_events, 1);
    boot_means(b, :) = mean(hr_snippets(idx, :), 1);
end
ci_lo = prctile(boot_means, 2.5, 1);
ci_hi = prctile(boot_means, 97.5, 1);

fill([t_eta fliplr(t_eta)], [ci_hi fliplr(ci_lo)], ...
     [0.8 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
plot(t_eta, eta_mean, '.-', 'Color', [0.15 0.3 0.7], ...
    'LineWidth', 1.8, 'MarkerSize', 14);
xline(0, '--k', 'LineWidth', 1);
hold off;
xlabel('Time from inspiration peak (ms)');
ylabel('Heart rate (bpm)');
title(sprintf('Test 02: Resp → HR ETA — bootstrap 95%% CI  (n = %d triggers, %d resamples)', ...
    n_events, n_boot));
xlim([t_eta(1)*1.05  t_eta(end)*1.05]);

%% ---- Figure 2: amplitude-stratified ETA --------------------------------
figure('Position', [100 100 700 400], 'Color', 'w');

% Breath amplitude (resp_amp_all) was collected in the main loop,
%  already aligned row-for-row with hr_snippets.
resp_amp = resp_amp_all;

% Split into tertiles
edges = prctile(resp_amp, [0 33.3 66.7 100]);
tertile_labels = {'Low', 'Medium', 'High'};
colors = [0.6 0.75 0.9;  0.2 0.5 0.8;  0.05 0.2 0.5];

hold on;
for g = 1:3
    if g == 3
        grp = resp_amp >= edges(g) & resp_amp <= edges(g+1);  % include upper edge
    else
        grp = resp_amp >= edges(g) & resp_amp < edges(g+1);
    end
    grp_mean = mean(hr_snippets(grp, :), 1);
    grp_sem  = std(hr_snippets(grp, :), 0, 1) / sqrt(sum(grp));

    fill([t_eta fliplr(t_eta)], ...
         [grp_mean + grp_sem, fliplr(grp_mean - grp_sem)], ...
         colors(g,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15, ...
         'HandleVisibility', 'off');
    plot(t_eta, grp_mean, '.-', 'Color', colors(g,:), ...
        'LineWidth', 1.5, 'MarkerSize', 10, ...
        'DisplayName', sprintf('%s (n=%d)', tertile_labels{g}, sum(grp)));
end
xline(0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
hold off;

xlabel('Time from inspiration peak (ms)');
ylabel('Heart rate (bpm)');
title('RSA by breath amplitude (tertiles)');
legend('Location', 'best');
xlim([t_eta(1)*1.05  t_eta(end)*1.05]);

%% ---- Figure 3: resp amplitude vs HR (mean & swing) ---------------------
figure('Position', [100 100 1000 400], 'Color', 'w');

% Per-trigger metrics
hr_mean_per_trigger = mean(hr_snippets, 2);    % mean HR across the ±5 beats
hr_swing_per_trigger = max(hr_snippets, [], 2) - min(hr_snippets, [], 2);  % max−min

% ---- Left: mean HR vs resp amplitude (shared autonomic state) -----------
subplot(1,2,1);
scatter(resp_amp, hr_mean_per_trigger, 4, [0.3 0.5 0.8], 'filled', 'MarkerFaceAlpha', 0.05);
hold on;
% Bin-averaged trend line for clarity
n_trend_bins = 20;
[~, ~, bin_id] = histcounts(resp_amp, n_trend_bins);
bin_counts   = accumarray(bin_id(bin_id>0), 1, [n_trend_bins 1]);
bin_hr_mean  = accumarray(bin_id(bin_id>0), hr_mean_per_trigger(bin_id>0), [n_trend_bins 1], @mean, NaN);
bin_amp_mean = accumarray(bin_id(bin_id>0), resp_amp(bin_id>0), [n_trend_bins 1], @mean, NaN);
ok_bins = bin_counts >= 5;                     % need ≥5 points per bin
plot(bin_amp_mean(ok_bins), bin_hr_mean(ok_bins), '.-k', 'LineWidth', 1.5, 'MarkerSize', 12);
hold off;

% Correlation
[r1, p1] = corr(resp_amp(:), hr_mean_per_trigger(:), 'type', 'Spearman');
xlabel('Resp amplitude at trigger (a.u.)');
ylabel('Mean HR around trigger (bpm)');
title(sprintf('Shared state:  r_s = %.3f,  p = %.1e', r1, p1));

% ---- Right: HR swing vs resp amplitude (RSA dose-response) --------------
subplot(1,2,2);
scatter(resp_amp, hr_swing_per_trigger, 4, [0.3 0.5 0.8], 'filled', 'MarkerFaceAlpha', 0.05);
hold on;
bin_swing_mean = accumarray(bin_id(bin_id>0), hr_swing_per_trigger(bin_id>0), [n_trend_bins 1], @mean, NaN);
plot(bin_amp_mean(ok_bins), bin_swing_mean(ok_bins), '.-k', 'LineWidth', 1.5, 'MarkerSize', 12);
hold off;

[r2, p2] = corr(resp_amp(:), hr_swing_per_trigger(:), 'type', 'Spearman');
xlabel('Resp amplitude at trigger (a.u.)');
ylabel('HR swing per trigger (bpm)');
title(sprintf('RSA dose-response:  r_s = %.3f,  p = %.1e', r2, p2));

sgtitle('Resp amplitude vs Heart rate');

%% ========================================================================
%  HELPER: intersect_sections
%  ========================================================================
function out = intersect_sections(A, B, min_dur)
%INTERSECT_SECTIONS  Overlap of two section catalogs.
%  Each input is a table with columns t_start, t_end, duration.
%  Returns only overlapping intervals ≥ min_dur seconds.

    rows = [];
    for ia = 1:height(A)
        for ib = 1:height(B)
            s = max(A.t_start(ia), B.t_start(ib));
            e = min(A.t_end(ia),   B.t_end(ib));
            d = e - s;
            if d >= min_dur
                rows = [rows; s, e, d]; %#ok<AGROW>
            end
        end
    end

    if isempty(rows)
        out = table([], [], [], 'VariableNames', {'t_start','t_end','duration'});
    else
        out = array2table(rows, 'VariableNames', {'t_start','t_end','duration'});
    end
end