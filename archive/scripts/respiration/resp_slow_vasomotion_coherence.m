%% Test 5: Resp ↔ Slow Vessel Coherence
%
% Computes magnitude-squared coherence and cross-spectrum phase between
% filtered respiration and the slow (<0.3 Hz) vessel diameter component.
% Uses the resp∩vessel section catalog — no ECG required.
%
% Method: Welch's averaged periodogram computed independently per section;
% coherence estimates pooled across sections weighted by the number of
% Welch windows each section contributes. Bootstrap CI by resampling
% sections with replacement. Significance threshold via the standard
% formula for K independent Welch windows.
%
% Phase is shown only where coherence exceeds the significance threshold
% (phase is unreliable where there is no detectable coupling).

%% ── Parameters ──────────────────────────────────────────────────────────
show_plots  = true;
verbose     = true;

FC_SLOW      = 0.3;    % Hz — lowpass cutoff for slow vessel component
WIN_SEC      = 30;     % Welch window length (s) → freq resolution ≈ 1/30 Hz
OVERLAP_FRAC = 0.5;    % Welch window overlap fraction
MIN_SEC      = 30;     % minimum section duration (s); must fit ≥1 Welch window
F_MAX_PLOT   = 5.0;    % Hz — upper frequency limit for plots
ALPHA        = 0.05;   % significance level for coherence threshold
N_BOOT       = 1000;   % bootstrap resamples for CI
BOOT_CI      = 95;     % CI level (%)

%% ── Load data ────────────────────────────────────────────────────────────
signal_data = load(".\data\signals.mat");
img_data    = load(".\data\im.mat");

% Physiology time axis aligned to camera clock
t_ECG       = double(signal_data.t_ECG(:));
t_ecgoffset = double(img_data.t_ecgoffset);
t_phys      = t_ECG - t_ecgoffset;
fs_phys     = round(1 / median(diff(t_phys)));   % Hz (should be 1000)

% Raw respiration signal
resp_raw = double(signal_data.resp(:));

% Vessel
t_frames = double(img_data.t_frames(:));
fs_cam   = 1 / median(diff(t_frames));

%% ── Preprocessing ───────────────────────────────────────────────────────

% Respiration — returns struct; good_sections uses sample indices into resp_raw
resp_out = preprocess_resp(resp_raw, fs_phys, 5, [0.5 5], 2, false, verbose);

resp_filt     = resp_out.resp_filt;
good_resp_sec = resp_out.good_sections;   % table: idx_start, idx_end, duration

% Vessel diameter (percent change, detrended)
diam_pct = preprocess_vessel_diam(img_data);

%% ── Build resp ∩ vessel good sections ───────────────────────────────────
% Convert resp sample indices → times using t_phys, then clip to camera window.

n_resp_sec = height(good_resp_sec);
good_sections = struct('t_start', {}, 't_end', {}, 'duration', {});
n_sec = 0;

for rs = 1:n_resp_sec
    % Convert sample indices to times
    t0_resp = t_phys(good_resp_sec.idx_start(rs));
    t1_resp = t_phys(good_resp_sec.idx_end(rs));

    % Clip to camera recording window
    t0  = max(t0_resp, t_frames(1));
    t1  = min(t1_resp, t_frames(end));
    dur = t1 - t0;
    if dur < MIN_SEC
        continue
    end
    n_sec = n_sec + 1;
    good_sections(n_sec).t_start  = t0;
    good_sections(n_sec).t_end    = t1;
    good_sections(n_sec).duration = dur;
end

if verbose
    total_dur = sum([good_sections.duration]);
    fprintf('Test 5: %d resp∩vessel sections ≥ %g s  (total %.1f s)\n', ...
        n_sec, MIN_SEC, total_dur);
end

if n_sec == 0
    error('Test 5: no qualifying sections. Lower MIN_SEC or check signal quality.');
end

%% ── Fixed Welch parameters ───────────────────────────────────────────────
nwin     = round(WIN_SEC * fs_cam);
noverlap = round(OVERLAP_FRAC * nwin);
nfft     = nwin;

% Frequency axis (identical for every section)
f_coh = (0 : nfft/2)' * (fs_cam / nfft);   % [nf × 1]
nf    = numel(f_coh);

%% ── Per-section coherence ────────────────────────────────────────────────
Coh_all = zeros(nf, n_sec);
Pxy_all = zeros(nf, n_sec) + 0i;
K_all   = zeros(1,  n_sec);

for is = 1:n_sec
    t0 = good_sections(is).t_start;
    t1 = good_sections(is).t_end;

    % Vessel frames within section
    frame_mask    = t_frames >= t0 & t_frames <= t1;
    t_sec         = t_frames(frame_mask);
    diam_pct_sec  = diam_pct(frame_mask);

    % Slow vessel: lowpass at FC_SLOW
    diam_slow_sec = lowpass(diam_pct_sec, FC_SLOW, fs_cam, ...
        'ImpulseResponse', 'iir', 'Steepness', 0.85);

    % Resp resampled onto camera frame times
    resp_sec = interp1(t_phys, resp_filt, t_sec, 'linear', NaN);

    % Drop any frames outside physiology time axis
    valid         = ~isnan(resp_sec) & ~isnan(diam_slow_sec);
    resp_sec      = resp_sec(valid);
    diam_slow_sec = diam_slow_sec(valid);

    K_s = floor((numel(resp_sec) - noverlap) / (nwin - noverlap));

    if K_s < 1
        if verbose
            fprintf('  Section %d: only %d samples (need %d), skipping.\n', ...
                is, numel(resp_sec), nwin);
        end
        continue
    end

    Coh_all(:, is) = mscohere(resp_sec, diam_slow_sec, ...
        hamming(nwin), noverlap, nfft, fs_cam);

    Pxy_all(:, is) = cpsd(resp_sec, diam_slow_sec, ...
        hamming(nwin), noverlap, nfft, fs_cam);

    K_all(is) = K_s;
end

% Remove skipped sections
valid_sec = K_all > 0;
Coh_all   = Coh_all(:, valid_sec);
Pxy_all   = Pxy_all(:, valid_sec);
K_all     = K_all(valid_sec);
n_valid   = sum(valid_sec);

if verbose
    fprintf('  Using %d of %d sections for coherence estimation.\n', n_valid, n_sec);
end

if n_valid == 0
    error('Test 5: no sections had sufficient data. Lower MIN_SEC or WIN_SEC.');
end

%% ── Pool across sections (weighted by K) ────────────────────────────────
w          = K_all / sum(K_all);
coh_mean   = Coh_all * w';
phase_mean = angle(Pxy_all * w');

%% ── Bootstrap CI ─────────────────────────────────────────────────────────
boot_coh = zeros(nf, N_BOOT);
rng(0, 'twister');
for b = 1:N_BOOT
    idx  = randi(n_valid, 1, n_valid);
    w_b  = K_all(idx) / sum(K_all(idx));
    boot_coh(:, b) = Coh_all(:, idx) * w_b';
end

alpha_tail = (100 - BOOT_CI) / 2;
coh_ci     = prctile(boot_coh, [alpha_tail, 100 - alpha_tail], 2)';  % [2 × nf]

%% ── Significance threshold ───────────────────────────────────────────────
K_total    = sum(K_all);
coh_thresh = 1 - ALPHA^(1 / (K_total - 1));

if verbose
    fprintf('  K_total = %d Welch windows → significance threshold = %.4f\n', ...
        K_total, coh_thresh);
end

%% ── Plots ────────────────────────────────────────────────────────────────
if ~show_plots; return; end

fmask     = f_coh <= F_MAX_PLOT;
f_plt     = f_coh(fmask);
coh_plt   = coh_mean(fmask);
ci_lo_plt = coh_ci(1, fmask);
ci_hi_plt = coh_ci(2, fmask);
phase_plt = phase_mean(fmask);

resp_band = [2.0, 4.0];  % Hz — typical mouse respiratory band (adjust if known)

col_main  = [0.18 0.39 0.69];
col_shade = [0.63 0.79 0.95];
col_resp  = [0.85 0.85 0.85];

figure('Color', 'w', 'Position', [100 100 800 580]);

% ── Panel 1: Magnitude-squared coherence ─────────────────────────────────
ax1 = subplot(2, 1, 1);
hold on;

fill([resp_band(1) resp_band(2) resp_band(2) resp_band(1)], ...
     [0 0 1 1], col_resp, 'EdgeColor', 'none', 'FaceAlpha', 0.6);

fill([f_plt; flipud(f_plt)], [ci_lo_plt'; flipud(ci_hi_plt')], ...
     col_shade, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

plot(f_plt, coh_plt, 'Color', col_main, 'LineWidth', 1.6);

yline(coh_thresh, '--r', 'LineWidth', 1.2, ...
    'Label', sprintf('p=%.2f  (K=%d)', ALPHA, K_total), ...
    'LabelHorizontalAlignment', 'right', ...
    'LabelVerticalAlignment',   'bottom');

hold off;
xlim([0, F_MAX_PLOT]);
ylim([0, max(1, max(ci_hi_plt) * 1.1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude-squared coherence');
title(sprintf('Resp ↔ Slow vessel coherence  (n=%d sections, %.0f s total)', ...
    n_valid, sum([good_sections(valid_sec).duration])));

text(mean(resp_band), ax1.YLim(2) * 0.97, 'Resp band', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Color', [0.4 0.4 0.4]);

% ── Panel 2: Cross-spectrum phase ────────────────────────────────────────
ax2 = subplot(2, 1, 2);
hold on;

fill([resp_band(1) resp_band(2) resp_band(2) resp_band(1)], ...
     [-180 -180 180 180], col_resp, 'EdgeColor', 'none', 'FaceAlpha', 0.6);

yline(0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);

% Mask phase where coherence is below threshold (unreliable)
phase_masked = rad2deg(phase_plt);
phase_masked(coh_plt < coh_thresh) = NaN;

plot(f_plt, phase_masked, 'Color', col_main, 'LineWidth', 1.6);

hold off;
xlim([0, F_MAX_PLOT]);
ylim([-180, 180]);
yticks(-180:90:180);
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Cross-spectrum phase  (shown only where coherence > threshold)');

text(mean(resp_band), 175, 'Resp band', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Color', [0.4 0.4 0.4]);

linkaxes([ax1 ax2], 'x');