%% TEST 02 — Respiration → Heart Rate ETA
%  Trigger on inspiration peaks, average instantaneous heart rate around
%  each trigger.  Tests for respiratory sinus arrhythmia: does the mouse's
%  heart speed up during inspiration and slow down during expiration?
%
%  Signals used : respiration + ECG (no vessel diameter)
%  Good-section catalog : resp ∩ ECG
%  Established conventions:
%    - All times in camera clock (t_phys = t_ECG − t_ecgoffset)
%    - Vectorized snippet extraction via index matrix
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
%  6. INSTANTANEOUS HEART RATE  (from R-peak times within good sections)
%  ========================================================================
%  RR intervals come from consecutive peaks inside ecg_sec_struct, which
%  already guarantees plausible intervals (no extra artifact rejection needed).

tHR_all = [];
HR_all  = [];

for s = 1:numel(ecg_sec_struct)
    pt = ecg_sec_struct(s).peak_time(:);       % R-peak times in this section
    rr = diff(pt);                              % seconds
    hr = 60 ./ rr;                              % bpm
    t_mid = pt(1:end-1) + rr/2;                % midpoint of each interval
    tHR_all = [tHR_all; t_mid];                %#ok<AGROW>
    HR_all  = [HR_all;  hr];                   %#ok<AGROW>
end

% --- interpolate onto uniform 100-Hz grid ---
%  100 Hz is plenty for HR fluctuations (changes at most every ~100 ms)
fs_hr = 100;
t_hr  = (t_rpeaks(1) : 1/fs_hr : t_rpeaks(end))';
HR_interp = interp1(tHR_all, HR_all, t_hr, 'linear', NaN);

% Light smoothing (0.2-s movmean ≈ 20 samples) to reduce interpolation jitter
HR_interp = smoothdata(HR_interp, 'movmean', round(0.2 * fs_hr));

fprintf('HR grid:  %.1f – %.1f s,  %d samples at %d Hz\n', ...
    t_hr(1), t_hr(end), numel(t_hr), fs_hr);

%% ========================================================================
%  7. FIND INSPIRATION PEAKS  (within good sections only)
%  ========================================================================
%  Peak detection on bandpass-filtered respiration (1000 Hz).
%  MinPeakDistance from typical mouse breathing: min ~0.15 s (= 1/~6.5 Hz)
%  MinPeakProminence: 0.5 × median |resp_filt| to avoid catching wiggles.

min_peak_dist_sec  = 0.15;                     % seconds
min_peak_dist_samp = round(min_peak_dist_sec * fs_phys);
min_prom = 1 * median(abs(resp_filt));

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
%  8. ETA — VECTORIZED SNIPPET EXTRACTION
%  ========================================================================
%  Window: symmetric ±0.5 s around each inspiration peak (one full breath)
pre_sec  = 0.5;
post_sec = 0.5;

n_pre  = round(pre_sec  * fs_hr);
n_post = round(post_sec * fs_hr);
n_win  = n_pre + n_post + 1;

% Convert inspiration times → nearest indices on the HR grid
insp_idx = round((insp_times_all - t_hr(1)) * fs_hr) + 1;

% Build offset vector:  -n_pre ... 0 ... +n_post
offsets = (-n_pre : n_post);                   % 1 × n_win

% Index matrix:  n_events × n_win
idx_matrix = insp_idx(:) + offsets;            % broadcast

% Boundary guard: discard events whose window exceeds the HR grid
valid_events = all(idx_matrix >= 1, 2) & all(idx_matrix <= numel(HR_interp), 2);
idx_matrix   = idx_matrix(valid_events, :);
n_events     = size(idx_matrix, 1);
fprintf('Valid ETA events:  %d / %d\n', n_events, numel(insp_times_all));

% Extract snippets (vectorised)
hr_snippets = HR_interp(idx_matrix);           % n_events × n_win

% Also discard any snippet with NaN (HR gap)
nan_rows    = any(isnan(hr_snippets), 2);
hr_snippets = hr_snippets(~nan_rows, :);
n_clean     = size(hr_snippets, 1);
fprintf('Clean ETA snippets (no NaN): %d\n', n_clean);

% --- compute mean ± SEM ---
eta_mean = mean(hr_snippets, 1);
eta_sem  = std(hr_snippets, 0, 1) / sqrt(n_clean);
t_eta    = offsets / fs_hr;                    % seconds relative to trigger

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

% HR trace (100-Hz grid)
mask_hr = t_hr >= ex_start & t_hr <= ex_end;
yyaxis right
plot(t_hr(mask_hr), HR_interp(mask_hr), 'Color', [0.85 0.33 0.1], ...
    'LineWidth', 1.2);
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

% ---- Panel 2: ETA  (mean ± SEM shaded) ----------------------------------
ax2 = subplot(2,1,2);

fill([t_eta fliplr(t_eta)], ...
     [eta_mean + eta_sem, fliplr(eta_mean - eta_sem)], ...
     [0.8 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
plot(t_eta, eta_mean, 'Color', [0.15 0.3 0.7], 'LineWidth', 1.8);
xline(0, '--k', 'LineWidth', 1);
hold off;

xlabel('Time from inspiration peak (s)');
ylabel('Heart rate (bpm)');
title(sprintf('Resp → HR ETA  (n = %d triggers)', n_clean));
xlim([t_eta(1) t_eta(end)]);

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