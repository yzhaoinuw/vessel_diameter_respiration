%% Test 4: Resp ↔ HR Cross-Correlation (Beat-Indexed)
%
% Strategy: stay in beat-indexed domain — no interpolation onto a uniform
% HR time grid. At each RR midpoint, compute instantaneous HR (60/RR) and
% sample resp_filt at that beat time. Cross-correlate the two beat-indexed
% sequences with lags expressed in beats (and converted to seconds via
% mean RR). Significance via bootstrap CI (resample sections with
% replacement).

show_plots = true;
verbose    = true;

%% ── 1. Load raw data ────────────────────────────────────────────────────
signal_data = load(".\data\signals.mat");
img_data    = load(".\data\im.mat");

% Physiology time axis (aligned to camera)
t_ECG       = double(signal_data.t_ECG(:));
t_ecgoffset = double(img_data.t_ecgoffset);
t_phys      = t_ECG - t_ecgoffset;

% Raw respiration signal (same length as t_phys, 1000 Hz)
resp_raw = double(signal_data.resp(:));

FS_PHYS = round(1 / median(diff(t_phys)));   % should be 1000 Hz

%% ── 2. Preprocess respiration ───────────────────────────────────────────
resp_out          = preprocess_resp(resp_raw, FS_PHYS, 5, [0.5 5], 2, false, verbose);
resp_filt         = resp_out.resp_filt;
good_resp_sections = resp_out.good_sections;   % table: idx_start, idx_end, duration

%% ── 3. Detect good ECG sections ─────────────────────────────────────────
% Extract R-peak times from signal_data and align to camera clock
t_rpeaks        = t_phys(double(signal_data.good_r_peaks));
ecg_sec_struct  = find_good_r_peak_sections(t_rpeaks, 5, false, verbose);

%% ── 4. Parameters ───────────────────────────────────────────────────────
MAX_LAG_BEATS = 25;      % ± lags to compute
N_BOOT        = 1000;    % bootstrap resamples
BOOT_CI       = 95;      % CI level (%)
MIN_BEATS     = 2 * MAX_LAG_BEATS + 10;   % minimum beats per section

%% ── 5. Pool all good R-peak times ───────────────────────────────────────
peak_times_all = [];
for ie = 1:numel(ecg_sec_struct)
    peak_times_all = [peak_times_all; ecg_sec_struct(ie).peak_time(:)]; %#ok<AGROW>
end
peak_times_all = sort(unique(peak_times_all));

%% ── 6. Build resp ∩ ECG section catalog ─────────────────────────────────
% Convert resp section sample indices → times using t_phys
resp_sec_t_start = t_phys(good_resp_sections.idx_start);
resp_sec_t_end   = t_phys(good_resp_sections.idx_end);

n_resp = height(good_resp_sections);
n_ecg  = numel(ecg_sec_struct);

overlap_t_start = [];
overlap_t_end   = [];

for ir = 1:n_resp
    rs = resp_sec_t_start(ir);
    re = resp_sec_t_end(ir);
    for ie = 1:n_ecg
        es = ecg_sec_struct(ie).t_start;
        ee = ecg_sec_struct(ie).t_end;
        ov_s = max(rs, es);
        ov_e = min(re, ee);
        if ov_e > ov_s
            overlap_t_start(end+1) = ov_s;
            overlap_t_end(end+1)   = ov_e;
        end
    end
end

n_sections = numel(overlap_t_start);
if verbose
    fprintf('resp ∩ ECG sections: %d\n', n_sections);
end

%% ── 7. Per-section cross-correlations ───────────────────────────────────
all_xcorr      = [];
section_nbeats = [];
all_rr_used    = [];   % RR intervals from sections that passed min-beats check

for is = 1:n_sections
    t_s = overlap_t_start(is);
    t_e = overlap_t_end(is);

    beats_in = peak_times_all(peak_times_all >= t_s & peak_times_all <= t_e);

    if numel(beats_in) < MIN_BEATS
        if verbose
            fprintf('  Section %d: only %d beats, skipping.\n', is, numel(beats_in));
        end
        continue
    end

    % Instantaneous HR at each inter-beat midpoint
    rr      = diff(beats_in);
    all_rr_used = [all_rr_used; rr]; %#ok<AGROW>
    hr_beat = 60 ./ rr;                       % bpm
    t_hr    = beats_in(1:end-1) + rr / 2;     % midpoint times

    % Sample resp_filt at midpoint times (nearest sample — no interpolation)
    idx_resp  = round((t_hr - t_phys(1)) * FS_PHYS) + 1;
    idx_resp  = max(1, min(numel(resp_filt), idx_resp));
    resp_beat = resp_filt(idx_resp);

    % Detrend within section (mean-centres both signals)
    hr_beat   = detrend(hr_beat,   'linear');
    resp_beat = detrend(resp_beat, 'linear');

    % Normalised cross-correlation
    [xc, lags] = xcorr(hr_beat, resp_beat, MAX_LAG_BEATS, 'normalized');

    all_xcorr(end+1, :)   = xc;
    section_nbeats(end+1) = numel(hr_beat);
end

n_valid = size(all_xcorr, 1);
if n_valid == 0
    error('No sections passed the minimum-beats threshold (%d beats).', MIN_BEATS);
end
if verbose
    fprintf('Sections used: %d / %d\n', n_valid, n_sections);
end

%% ── 8. Weighted mean cross-correlation ──────────────────────────────────
w          = section_nbeats(:) / sum(section_nbeats);
xcorr_mean = sum(all_xcorr .* w, 1);

%% ── 9. Bootstrap CI (resample sections with replacement) ─────────────────
boot_xcorr = zeros(N_BOOT, size(all_xcorr, 2));
rng(42);
for ib = 1:N_BOOT
    idx_boot = randi(n_valid, 1, n_valid);
    w_b      = section_nbeats(idx_boot) / sum(section_nbeats(idx_boot));
    boot_xcorr(ib, :) = sum(all_xcorr(idx_boot, :) .* w_b(:), 1);
end

alpha    = (100 - BOOT_CI) / 2;
xcorr_ci = prctile(boot_xcorr, [alpha, 100 - alpha], 1);   % [2 × nlags]

%% ── 10. Lag axes ─────────────────────────────────────────────────────────
mean_rr   = mean(all_rr_used);
lag_beats = lags;
lag_sec   = lags * mean_rr;

if verbose
    fprintf('Mean RR = %.1f ms  (mean HR ≈ %.0f bpm)\n', mean_rr * 1000, 60 / mean_rr);
end

%% ── 11. Summary ──────────────────────────────────────────────────────────
if verbose
    [pk_val,  pk_idx] = max(xcorr_mean);
    [tr_val,  tr_idx] = min(xcorr_mean);
    ci_zero = xcorr_ci(:, lag_beats == 0);
    fprintf('\n── Cross-correlation summary ──\n');
    fprintf('  Peak:   r = %+.3f at lag %+.0f ms (%+.1f beats)\n', ...
        pk_val, lag_sec(pk_idx) * 1000, lag_beats(pk_idx));
    fprintf('  Trough: r = %+.3f at lag %+.0f ms (%+.1f beats)\n', ...
        tr_val, lag_sec(tr_idx) * 1000, lag_beats(tr_idx));
    fprintf('  Bootstrap CI at lag 0: [%.3f, %.3f]\n', ci_zero(1), ci_zero(2));
end

%% ── 12. Plots ────────────────────────────────────────────────────────────
if show_plots

    figure('Name', 'Resp–HR Cross-Correlation', 'Color', 'w', ...
        'Position', [100 100 700 480]);

    ax = axes; hold(ax, 'on');

    fill(ax, [lag_sec, fliplr(lag_sec)], ...
         [xcorr_ci(1,:), fliplr(xcorr_ci(2,:))], ...
         [0.7 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

    plot(ax, lag_sec([1 end]), [0 0], 'k--', 'LineWidth', 0.8);
    plot(ax, lag_sec, xcorr_mean, 'b-', 'LineWidth', 2);

    [~, pk_idx] = max(abs(xcorr_mean));
    plot(ax, lag_sec(pk_idx), xcorr_mean(pk_idx), 'ro', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(ax, lag_sec(pk_idx), xcorr_mean(pk_idx) - 0.03 * sign(xcorr_mean(pk_idx)), ...
        sprintf('  r=%.3f\n  %+.0f ms', xcorr_mean(pk_idx), lag_sec(pk_idx)*1000), ...
        'FontSize', 9, 'Color', 'r', 'VerticalAlignment', 'top');

    xlabel(ax, 'Lag (s)  [positive = HR leads resp]');
    ylabel(ax, 'Normalised cross-correlation');
    title(ax, sprintf('Resp ↔ HR cross-correlation  (n=%d sections, %d%% CI)', n_valid, BOOT_CI));
    xlim(ax, lag_sec([1 end]));
    box(ax, 'on');

    % Beat-lag axis on top
    ax2 = axes('Position', ax.Position, 'XAxisLocation', 'top', ...
               'Color', 'none', 'XLim', lag_beats([1 end]), ...
               'YLim', ax.YLim, 'YTick', [], 'YColor', 'none');
    xlabel(ax2, 'Lag (beats)');

end