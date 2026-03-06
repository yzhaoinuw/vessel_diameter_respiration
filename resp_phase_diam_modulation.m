%% test_03_resp_phase_cardiac_PAC.m
%
% Test 3: Phase-amplitude coupling (PAC)
%   Phase signal  : instantaneous respiratory phase (from preprocess_resp)
%   Amplitude signal: Hilbert envelope of cardiac-band vessel diameter
%
% Question: does the amplitude of cardiac pulsations in vessel diameter
% vary systematically with respiratory phase? (i.e., does breathing gate
% how strongly heartbeats move the vessel?)
%
% Sections catalog used: resp ∩ vessel (no ECG required)
% Metric: Modulation Index (Tort et al. 2010, J Neurophysiol)
% Significance: permutation test via circular phase shifts
%
% Follows the same data loading / preprocessing pattern as
% test_01_resp_vessel_ETA.m

clear; clc; close all;

%% ── Parameters ────────────────────────────────────────────────────────────
arguments_block = struct( ...
    'band_cardiac',   [5 15],  ... % Hz — cardiac band in vessel signal
    'detrend_win_s',  100,     ... % s  — moving-mean detrend window (10× slowest vasomotion)
    'min_sec_dur',    5,       ... % s  — minimum section length
    'n_bins',         18,      ... % phase bins (20° each, Tort et al. standard)
    'n_perm',         500,     ... % circular-shift permutations for null MI distribution
    'show_plots',     true,    ...
    'verbose',        true     ...
);
p = arguments_block;

%% ── 1. Load data ──────────────────────────────────────────────────────────
signal_data = load('.\data\signals.mat');
img_data    = load('.\data\im.mat');

t_frames  = double(img_data.t_frames(:));
fs_cam    = 1 / median(diff(t_frames));

t_ecgoffset = double(img_data.t_ecgoffset);
t_phys      = double(signal_data.t_ECG(:)) - t_ecgoffset;  % physiology on camera clock
resp_raw    = double(signal_data.resp(:));
fs_phys     = 1 / median(diff(t_phys));

%% ── 2. Preprocess vessel diameter ────────────────────────────────────────
% preprocess_vessel_diam(img_data, baseline_window, show_plot)
% baseline_window: 100 s = 10× slowest vasomotion period (~10 s)
diam_norm = preprocess_vessel_diam(img_data, p.detrend_win_s, false);

%% ── 3. Preprocess respiration ────────────────────────────────────────────
resp_struct = preprocess_resp(resp_raw, fs_phys, ...
    p.min_sec_dur, [0.5 5], 2, false, p.verbose);

% resp_struct fields used here:
%   .resp_filt      — bandpass-filtered resp at fs_phys
%   .phase          — instantaneous phase (angle of Hilbert), fs_phys
%   .good_sections  — table: idx_start, idx_end (indices into t_phys)

%% ── 4. Resample resp phase to camera grid ─────────────────────────────────
% preprocess_resp computes phase at 1000 Hz (t_phys).
% We need it at camera rate to pair with vessel data.
resp_phase_cam = interp1(t_phys, resp_struct.phase, t_frames, 'linear');

%% ── 5. Compute cardiac-band Hilbert envelope of vessel diameter ───────────
% Bandpass vessel at camera rate, then take abs(hilbert) to get envelope.
% Note: fs_cam ≈ 55.56 Hz → Nyquist ≈ 27.8 Hz, so 5–15 Hz band is fine.
diam_card     = bandpass(diam_norm, p.band_cardiac, fs_cam);
diam_card_env = abs(hilbert(diam_card));   % instantaneous amplitude envelope

%% ── 6. Find resp ∩ vessel good sections ──────────────────────────────────
% resp good sections are indexed into t_phys (1000 Hz).
% Convert to wall-clock times, then find overlapping camera frames.

n_resp_sec = height(resp_struct.good_sections);
good_sections = struct('t_start', {}, 't_end', {}, 'cam_idx', {}, 'duration', {});
sec_count = 0;

for s = 1:n_resp_sec
    t_start = t_phys(resp_struct.good_sections.idx_start(s));
    t_end   = t_phys(resp_struct.good_sections.idx_end(s));

    % Camera frames within this resp section
    cam_mask = t_frames >= t_start & t_frames <= t_end;
    if sum(cam_mask) < round(p.min_sec_dur * fs_cam)
        continue
    end

    % Quick vessel quality check: reject if >20% of frames are NaN
    env_seg = diam_card_env(cam_mask);
    if mean(isnan(env_seg)) > 0.20
        continue
    end

    sec_count = sec_count + 1;
    good_sections(sec_count).t_start  = t_start;
    good_sections(sec_count).t_end    = t_end;
    good_sections(sec_count).cam_idx  = find(cam_mask);
    good_sections(sec_count).duration = t_end - t_start;
end

n_sections   = numel(good_sections);
total_data_s = sum([good_sections.duration]);

if p.verbose
    fprintf('\nTest 3 — Resp phase → vessel cardiac envelope PAC\n');
    fprintf('  Good sections (resp ∩ vessel): %d\n', n_sections);
    fprintf('  Total data: %.1f s\n', total_data_s);
end

if n_sections == 0
    error('No valid resp ∩ vessel sections found. Check preprocessing outputs.');
end

%% ── 7. Pool (phase, amplitude) pairs across all sections ─────────────────
all_phase = [];
all_amp   = [];

for s = 1:n_sections
    idx = good_sections(s).cam_idx;
    ph  = resp_phase_cam(idx);
    am  = diam_card_env(idx);

    % Drop any NaN frames within the section
    valid = ~isnan(ph) & ~isnan(am);
    all_phase = [all_phase; ph(valid)];  
    all_amp   = [all_amp;   am(valid)];  
end

%% ── 8. Compute observed MI (Tort et al. 2010) ────────────────────────────
% MI = (log(N) + sum(p .* log(p))) / log(N)
% where p(j) = mean amplitude in bin j / sum of mean amplitudes across bins
% MI = 0 → amplitude uniform across phase (no coupling)
% MI = 1 → all amplitude in one bin (perfect coupling)

n_bins  = p.n_bins;
edges   = linspace(-pi, pi, n_bins + 1);
bin_amp = zeros(n_bins, 1);

for b = 1:n_bins
    mask = all_phase >= edges(b) & all_phase < edges(b+1);
    bin_amp(b) = mean(all_amp(mask), 'omitnan');
end

% Normalize to probability distribution
p_dist    = bin_amp / sum(bin_amp);
MI_obs    = (log(n_bins) + sum(p_dist .* log(p_dist + eps))) / log(n_bins);

bin_centers = (edges(1:end-1) + edges(2:end)) / 2;  % radians, for plotting

if p.verbose
    fprintf('  Observed MI: %.5f\n', MI_obs);
end

%% ── 9. Permutation test (circular phase shifts) ──────────────────────────
% Circular-shift the phase vector by a random offset (at least 1 full breath
% cycle = ~fs_cam / resp_rate ≈ 14–55 frames) to destroy phase-amplitude
% pairing while preserving each signal's autocorrelation structure.
n_perm     = p.n_perm;
MI_null    = zeros(n_perm, 1);
N_samples  = numel(all_phase);

% Minimum shift: 1 breath cycle at slowest expected rate (0.5 Hz)
min_shift  = round(fs_cam / 0.5);

rng(42);  % reproducibility
for k = 1:n_perm
    shift  = min_shift + randi(N_samples - 2*min_shift);
    ph_shf = circshift(all_phase, shift);

    bin_amp_null = zeros(n_bins, 1);
    for b = 1:n_bins
        mask = ph_shf >= edges(b) & ph_shf < edges(b+1);
        bin_amp_null(b) = mean(all_amp(mask), 'omitnan');
    end
    p_null         = bin_amp_null / sum(bin_amp_null);
    MI_null(k)     = (log(n_bins) + sum(p_null .* log(p_null + eps))) / log(n_bins);
end

MI_null_mean = mean(MI_null);
MI_null_std  = std(MI_null);
MI_zscore    = (MI_obs - MI_null_mean) / MI_null_std;
MI_pval      = mean(MI_null >= MI_obs);  % one-sided: fraction of null >= observed

if p.verbose
    fprintf('  Null MI: mean=%.5f, std=%.5f\n', MI_null_mean, MI_null_std);
    fprintf('  MI z-score: %.2f\n', MI_zscore);
    fprintf('  Permutation p-value: %.4f (n_perm=%d)\n', MI_pval, n_perm);
end

%% ── 9b. Equalized-bin verification ───────────────────────────────────────
% Uneven phase sampling (bottom-left panel) can inflate MI even when no
% true coupling exists. To verify the z-score is not a sampling artifact,
% subsample each bin to the minimum bin count and recompute MI + permutation.

% Collect per-bin sample indices
bin_indices = cell(n_bins, 1);
for b = 1:n_bins
    bin_indices{b} = find(all_phase >= edges(b) & all_phase < edges(b+1));
end
bin_counts  = cellfun(@numel, bin_indices);
n_eq        = min(bin_counts);   % equalized count = smallest bin

rng(42);
bin_amp_eq = zeros(n_bins, 1);
for b = 1:n_bins
    sub_idx        = bin_indices{b}(randperm(bin_counts(b), n_eq));
    bin_amp_eq(b)  = mean(all_amp(sub_idx));
end
p_dist_eq  = bin_amp_eq / sum(bin_amp_eq);
MI_eq      = (log(n_bins) + sum(p_dist_eq .* log(p_dist_eq + eps))) / log(n_bins);

% Permutation test on equalized data
MI_null_eq = zeros(n_perm, 1);
% Pool the equalized samples into flat phase/amp vectors for shifting
all_phase_eq = [];
all_amp_eq   = [];
for b = 1:n_bins
    sub_idx      = bin_indices{b}(randperm(bin_counts(b), n_eq));
    all_phase_eq = [all_phase_eq; all_phase(sub_idx)]; %#ok<AGROW>
    all_amp_eq   = [all_amp_eq;   all_amp(sub_idx)];   %#ok<AGROW>
end
N_eq = numel(all_phase_eq);

for k = 1:n_perm
    shift     = min_shift + randi(N_eq - 2*min_shift);
    ph_shf_eq = circshift(all_phase_eq, shift);

    bin_amp_null_eq = zeros(n_bins, 1);
    for b = 1:n_bins
        mask = ph_shf_eq >= edges(b) & ph_shf_eq < edges(b+1);
        bin_amp_null_eq(b) = mean(all_amp_eq(mask), 'omitnan');
    end
    p_null_eq        = bin_amp_null_eq / sum(bin_amp_null_eq);
    MI_null_eq(k)    = (log(n_bins) + sum(p_null_eq .* log(p_null_eq + eps))) / log(n_bins);
end

MI_eq_null_mean = mean(MI_null_eq);
MI_eq_null_std  = std(MI_null_eq);
MI_eq_zscore    = (MI_eq - MI_eq_null_mean) / MI_eq_null_std;
MI_eq_pval      = mean(MI_null_eq >= MI_eq);

if p.verbose
    fprintf('\n  — Equalized-bin verification (n_eq=%d per bin) —\n', n_eq);
    fprintf('  Equalized MI        : %.5f\n', MI_eq);
    fprintf('  Null MI (eq)        : mean=%.5f ± %.5f\n', MI_eq_null_mean, MI_eq_null_std);
    fprintf('  z-score (eq)        : %.2f\n', MI_eq_zscore);
    fprintf('  p (eq, permutation) : %.4f\n', MI_eq_pval);
    if MI_eq_pval >= 0.05
        fprintf('  → z-score collapses after equalization: original result was sampling artifact\n');
    else
        fprintf('  → z-score survives equalization: effect is real\n');
    end
end

%% ── 10. Plots ─────────────────────────────────────────────────────────────
if ~p.show_plots; return; end

figure('Name', 'Test 3: Resp Phase → Vessel Cardiac Envelope PAC', ...
       'Position', [100 80 1100 820]);

% ── 10a. Mean cardiac envelope by respiratory phase ──
subplot(2,2,1);
bar(rad2deg(bin_centers), bin_amp, 'FaceColor', [0.27 0.60 0.85], 'EdgeColor', 'none');
hold on;
plot(rad2deg(bin_centers), bin_amp, 'k-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'MarkerFaceColor', 'k');
xlabel('Respiratory phase (°)');
ylabel('Mean cardiac envelope (μm)');
title(sprintf('PAC: resp phase → vessel cardiac amplitude\nMI = %.4f,  z = %.1f,  p = %.3f', ...
    MI_obs, MI_zscore, MI_pval));
xticks(-180:60:180);
xline(0, '--', 'color', [0.5 0.5 0.5]);
xlim([-180 180]);

% ── 10b. Null distribution of MI ──
subplot(2,2,2);
histogram(MI_null, 40, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
hold on;
xline(MI_obs, 'r-', 'LineWidth', 2);
xline(MI_eq,  'g--', 'LineWidth', 2);
xlabel('MI (null)');
ylabel('Count');
title(sprintf('Permutation null distribution\nObserved MI (red): %.5f  |  Equalized MI (green): %.5f', MI_obs, MI_eq));
legend('Null MI', 'Observed', 'Equalized', 'Location', 'northeast');

% ── 10c. Phase histogram (data coverage) ──
subplot(2,2,3);
histogram(rad2deg(all_phase), 36, 'FaceColor', [0.4 0.75 0.4], 'EdgeColor', 'none');
xlabel('Respiratory phase (°)');
ylabel('Sample count');
title('Phase sampling coverage');
xticks(-180:60:180);
xlim([-180 180]);

% ── 10d. Polar plot of amplitude by phase ──
subplot(2,2,4);
theta_polar = [bin_centers(:); bin_centers(1)];
rho_polar   = [bin_amp(:);     bin_amp(1)];
polarplot(theta_polar, rho_polar, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor','b');
hold on;
% Mark inspiration peak (phase = 0, convention from preprocess_resp)
polarplot([0 0], [0 max(bin_amp)*1.05], 'r--', 'LineWidth', 1.2);
title('Polar: amplitude by resp phase');
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'counterclockwise';

sgtitle(sprintf('Test 3 — Resp phase → vessel cardiac envelope PAC  |  %d sections, %.0f s total', ...
    n_sections, total_data_s), 'FontWeight', 'bold');

%% ── 11. Summary to console ────────────────────────────────────────────────
fprintf('\n── Test 3 Summary ──\n');
fprintf('  Sections used   : %d\n',    n_sections);
fprintf('  Total data      : %.1f s\n', total_data_s);
fprintf('  Phase bins      : %d (%.0f° each)\n', n_bins, 360/n_bins);
fprintf('  Observed MI     : %.5f\n',  MI_obs);
fprintf('  Null MI (mean)  : %.5f ± %.5f\n', MI_null_mean, MI_null_std);
fprintf('  z-score         : %.2f\n',  MI_zscore);
fprintf('  p (permutation) : %.4f\n',  MI_pval);
fprintf('  — Equalized-bin (n=%d per bin) —\n', n_eq);
fprintf('  Equalized MI    : %.5f\n',  MI_eq);
fprintf('  z-score (eq)    : %.2f\n',  MI_eq_zscore);
fprintf('  p (eq)          : %.4f\n',  MI_eq_pval);
if MI_eq_pval >= 0.05
    fprintf('  → Null confirmed: z-score collapses after equalization (sampling artifact)\n');
else
    fprintf('  → PAC is REAL: z-score survives equalization\n');
end