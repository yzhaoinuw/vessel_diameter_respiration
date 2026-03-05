%% TEST 02 — Respiration → Heart Rate ETA
%  Trigger on inspiration peaks, average instantaneous heart rate around
%  each trigger.  Tests for respiratory sinus arrhythmia: does the mouse's
%  heart speed up during inspiration and slow down during expiration?
%
%  Signals used : respiration + ECG (no vessel diameter)
%  Good-section catalog : resp ∩ ECG
%  Established conventions:
%    - All times in camera clock (t_phys = t_ECG − t_ecgoffset)
%    - Beat-indexed ETA on native R-peak times (no interpolation)
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

%% ========================================================================
%  3. RESPIRATION PREPROCESSING  (native 1000 Hz)
%  ========================================================================
band_resp = [0.5 5];                           % Hz — mouse breathing range
resp_filt = bandpass(resp, band_resp, fs_phys);

% --- sliding-window std for amplitude estimation (2-s window) ---
amp_win_sec  = 2;                              % captures 4–8 breath cycles
amp_win_samp = round(amp_win_sec * fs_phys);
resp_std = movstd(resp, amp_win_samp);

% --- identify quiet-breathing sections ---
%  High resp_std → active/clipped breathing; low → quiet usable epochs.
%  Threshold: use median of resp_std as a reasonable split.
%  (User may adjust based on diagnostic plot.)
resp_std_thresh = median(resp_std);
quiet_mask = resp_std < resp_std_thresh;       % logical, 1000 Hz

% Convert to contiguous sections (start/end indices in t_phys)
d_quiet = diff([0; quiet_mask(:); 0]);
resp_starts_idx = find(d_quiet ==  1);
resp_ends_idx   = find(d_quiet == -1) - 1;

min_dur_sec = 5;
resp_durs   = (resp_ends_idx - resp_starts_idx + 1) / fs_phys;
keep_resp   = resp_durs >= min_dur_sec;

resp_sections = table( ...
    t_phys(resp_starts_idx(keep_resp)), ...
    t_phys(resp_ends_idx(keep_resp)),   ...
    resp_durs(keep_resp), ...
    'VariableNames', {'t_start','t_end','duration'});

fprintf('Resp good sections:  %d  (%.1f s total)\n', ...
    height(resp_sections), sum(resp_sections.duration));

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
%  6. COLLECT R-PEAK TIMES FROM GOOD ECG SECTIONS
%  ========================================================================
%  HR will be computed beat-by-beat in the ETA step. Here we just need
%  all valid R-peak times in one sorted vector for nearest-peak lookup.

t_rpeaks_good = [];
for s = 1:numel(ecg_sec_struct)
    t_rpeaks_good = [t_rpeaks_good; ecg_sec_struct(s).peak_time(:)]; %#ok<AGROW>
end
t_rpeaks_good = sort(t_rpeaks_good);

fprintf('Good R-peaks:  %d  spanning %.1f – %.1f s\n', ...
    numel(t_rpeaks_good), t_rpeaks_good(1), t_rpeaks_good(end));

%% ========================================================================
%  7. FIND INSPIRATION PEAKS  (within good sections only)
%  ========================================================================
%  Peak detection on bandpass-filtered respiration (1000 Hz).
%  MinPeakDistance from typical mouse breathing: min ~0.15 s (= 1/~6.5 Hz)
%  MinPeakProminence: 0.5 × median |resp_filt| to avoid catching wiggles.

min_peak_dist_sec  = 0.15;                     % seconds
min_peak_dist_samp = round(min_peak_dist_sec * fs_phys);
min_prom = 0.5 * median(abs(resp_filt));

insp_times_all = [];                           % collect across sections

for s = 1:height(good_sections)
    idx_start = find(t_phys >= good_sections.t_start(s), 1, 'first');
    idx_end   = find(t_phys <= good_sections.t_end(s),   1, 'last');
    seg       = resp_filt(idx_start:idx_end);

    [~, locs] = findpeaks(seg, ...
        'MinPeakDistance',   min_peak_dist_samp, ...
        'MinPeakProminence', min_prom);

    insp_times_all = [insp_times_all; t_phys(idx_start - 1 + locs)]; %#ok<AGROW>
end

fprintf('Inspiration peaks found: %d\n', numel(insp_times_all));

%% ========================================================================
%  8. ETA — BEAT-INDEXED
%  ========================================================================
%  For each inspiration peak, find the nearest R-peak, then grab ±N beats.
%  HR for beat i = 60 / (t_rpeak(i+1) − t_rpeak(i)).
%  Beat index 0 = the RR interval starting at the nearest R-peak.
%  All triggers align by beat number — no interpolation, no binning.

n_beats_half = 5;                              % ±5 beats around trigger
beat_offsets = -n_beats_half : n_beats_half;   % -5 ... 0 ... +5
n_beat_win   = numel(beat_offsets);

% For each trigger, find nearest R-peak index
n_rp = numel(t_rpeaks_good);
insp_rp_idx = interp1(t_rpeaks_good, 1:n_rp, insp_times_all, 'nearest');
insp_rp_idx = round(insp_rp_idx);

% Boundary guard: need n_beats_half extra peaks on each side for HR calc
%  HR at beat i requires peak i and peak i+1, so rightmost beat offset
%  needs peak at (idx + n_beats_half + 1)
valid = insp_rp_idx - n_beats_half >= 1 & ...
        insp_rp_idx + n_beats_half + 1 <= n_rp & ...
        ~isnan(insp_rp_idx);
insp_rp_idx = insp_rp_idx(valid);
n_events = numel(insp_rp_idx);
fprintf('Valid ETA events: %d / %d\n', n_events, numel(insp_times_all));

% Build index matrix: each row is one trigger, columns are beat offsets
%  idx_matrix(k, j) → index into t_rpeaks_good for beat offset j
idx_matrix = insp_rp_idx(:) + beat_offsets;    % n_events × n_beat_win

% HR for each beat = 60 / (next_peak − this_peak)
%  So for beat offset j, HR = 60 / (t_rpeaks_good(idx+1) − t_rpeaks_good(idx))
hr_snippets = 60 ./ (t_rpeaks_good(idx_matrix + 1) - t_rpeaks_good(idx_matrix));

% --- mean ± SEM ---
eta_mean = mean(hr_snippets, 1);
eta_sem  = std(hr_snippets, 0, 1) / sqrt(n_events);

fprintf('ETA: %d triggers, %d beats per snippet\n', n_events, n_beat_win);

%% ========================================================================
%  9. DIAGNOSTIC PLOTS
%  ========================================================================

figure('Position', [100 100 900 700], 'Color', 'w');

% ---- Panel 1: example good section (resp + HR + peak markers) -----------
ax1 = subplot(2,1,1);

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

% HR trace — beat-by-beat values plotted at R-peak times
rr_good = diff(t_rpeaks_good);
hr_good = 60 ./ rr_good;
t_hr_plot = t_rpeaks_good(1:end-1);            % HR lives at each beat
mask_hr = t_hr_plot >= ex_start & t_hr_plot <= ex_end;
yyaxis right
plot(t_hr_plot(mask_hr), hr_good(mask_hr), '.-', 'Color', [0.85 0.33 0.1], ...
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
title(sprintf('Example good section (%.1f s)', ex_end - ex_start));
xlim([ex_start ex_end]);

% ---- Panel 2: ETA  (mean ± SEM, beat-indexed) ---------------------------
ax2 = subplot(2,1,2);

fill([beat_offsets fliplr(beat_offsets)], ...
     [eta_mean + eta_sem, fliplr(eta_mean - eta_sem)], ...
     [0.8 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
plot(beat_offsets, eta_mean, '.-', 'Color', [0.15 0.3 0.7], ...
    'LineWidth', 1.8, 'MarkerSize', 14);
xline(0, '--k', 'LineWidth', 1);
hold off;

xlabel('Beat relative to inspiration');
ylabel('Heart rate (bpm)');
title(sprintf('Resp → HR ETA  (n = %d triggers, mean ± SEM)', n_events));
xlim([beat_offsets(1)-0.5  beat_offsets(end)+0.5]);
xticks(beat_offsets);

sgtitle('Test 02: Respiratory Sinus Arrhythmia');

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