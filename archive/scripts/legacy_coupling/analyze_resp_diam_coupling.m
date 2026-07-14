%% analyze_resp_diam_coupling.m

% =========================================================================
%% 0. Load data
% =========================================================================
data_folder    = ".\data\";
subject_name   = "F154-dex40-19022025";
subject_folder = fullfile(data_folder, subject_name);

signal_data = load(fullfile(subject_folder, "signals.mat"));
img_data    = load(fullfile(subject_folder, "im.mat"));

t_ecg = double(signal_data.t_ECG(:));
resp  = double(signal_data.resp(:));

% Sampling rates
fs_ecg = 1 / median(diff(t_ecg));
fs_cam = 1 / median(diff(double(img_data.t_frames(:))));

% Align ECG/resp time axis to camera
t_ecgoffset = double(img_data.t_ecgoffset);
assert(t_ecgoffset > 0, 't_ecgoffset must be positive.')
t_ecg = t_ecg - t_ecgoffset;

% =========================================================================
%% 1. Vessel diameter preprocessing
% =========================================================================
baseline_window = 10 * 60;  % 10-min moving-average baseline
diam_norm = preprocess_vessel_diam(img_data, baseline_window, false);

% =========================================================================
%% 2. Respiration preprocessing
%    Low-pass at 5 Hz to remove cardiac contamination.
%    Mouse breathing: ~1-4 Hz.  Mouse cardiac: ~8-12 Hz.
% =========================================================================
Fc_resp       = 5;
[b_lp, a_lp] = butter(4, Fc_resp / (fs_ecg/2), 'low');
resp_filt     = filtfilt(b_lp, a_lp, detrend(resp));

% =========================================================================
%% 3. Minute bin edges and centers
% =========================================================================
skip_sec    = 14 * 60;      % skip first 14 minutes for dominant freq
t_start     = t_ecg(1);
t_end       = t_ecg(end);
bin_edges   = t_start : 60 : t_end;
bin_centers = bin_edges(1:end-1) + 30;
n_bins      = length(bin_centers);

% =========================================================================
%% 4. Spectrograms
%    Display band  : [0.1, 5] Hz — full range keeps slow vasomotion visible
%    Dominant freq : [1.5, 5] Hz — respiratory band only (Option A)
% =========================================================================
win_sec     = 60;           % 60 s window -> 1/60 Hz frequency resolution
overlap_pct = 0.75;         % 75% overlap for smooth time axis
f_disp      = [0.1, 5];     % Hz — full display range for both spectrograms
f_resp_band = [1.5, 5];     % Hz — respiratory band for dominant freq extraction

% --- Respiration spectrogram (runs at fs_ecg) ---
win_r  = round(win_sec * fs_ecg);
hop_r  = round(win_r * (1 - overlap_pct));
nfft_r = 2^nextpow2(win_r);
[S_resp, f_resp, t_resp] = spectrogram(resp_filt, hann(win_r), ...
                                        win_r - hop_r, nfft_r, fs_ecg);
t_resp   = t_resp + t_ecg(1);
f_disp_r = f_resp >= f_disp(1)      & f_resp <= f_disp(2);
f_dom_r  = f_resp >= f_resp_band(1) & f_resp <= f_resp_band(2);

% --- Vessel diameter spectrogram (runs at fs_cam) ---
win_d  = round(win_sec * fs_cam);
hop_d  = round(win_d * (1 - overlap_pct));
nfft_d = 2^nextpow2(win_d);
[S_diam, f_diam, t_diam] = spectrogram(diam_norm, hann(win_d), ...
                                        win_d - hop_d, nfft_d, fs_cam);
t_diam   = t_diam + t_ecg(1);
f_disp_d = f_diam >= f_disp(1)      & f_diam <= f_disp(2);
f_dom_d  = f_diam >= f_resp_band(1) & f_diam <= f_resp_band(2);

% =========================================================================
%% 5. Dominant frequency per minute bin (skip first 14 min)
%    Both signals searched within the respiratory band [1.5, 5] Hz
%    Vessel diam: bins below the 25th-percentile peak power are rejected as noise
% =========================================================================
dom_freq_resp = NaN(n_bins, 1);
dom_freq_diam = NaN(n_bins, 1);

% Pre-pass: collect peak power per bin for vessel diam to set noise threshold
peak_power_diam = NaN(n_bins, 1);
for i = 1:n_bins
    cols_d = t_diam >= bin_edges(i) & t_diam < bin_edges(i+1);
    if any(cols_d)
        P = mean(abs(S_diam(f_dom_d, cols_d)).^2, 2);
        peak_power_diam(i) = max(P);
    end
end
power_thresh = prctile(peak_power_diam, 50);

% Main loop
for i = 1:n_bins
    if bin_centers(i) - t_start < skip_sec
        continue
    end

    % Respiration
    cols_r = t_resp >= bin_edges(i) & t_resp < bin_edges(i+1);
    if any(cols_r)
        P                = mean(abs(S_resp(f_dom_r, cols_r)).^2, 2);
        [~, fi]          = max(P);
        f_sub            = f_resp(f_dom_r);
        dom_freq_resp(i) = f_sub(fi);
    end

    % Vessel diameter — skip bins that don't clear the noise floor
    cols_d = t_diam >= bin_edges(i) & t_diam < bin_edges(i+1);
    if any(cols_d)
        P = mean(abs(S_diam(f_dom_d, cols_d)).^2, 2);
        if max(P) >= power_thresh
            [~, fi]          = max(P);
            f_sub            = f_diam(f_dom_d);
            dom_freq_diam(i) = f_sub(fi);
        end
    end
end

% =========================================================================
%% 6. Plot
% =========================================================================
figure('Position', [100 100 1300 1000]);
t_min = bin_centers / 60;   % shared x-axis in minutes

% --- Subplot 1: Respiration spectrogram ---
ax1 = subplot(3,1,1);
imagesc(t_resp/60, f_resp(f_disp_r), ...
        10*log10(abs(S_resp(f_disp_r,:)).^2 + eps));
axis xy;
colormap(ax1, 'turbo');
cb = colorbar; cb.Label.String = 'Power (dB)';
hold on;
plot(t_min, dom_freq_resp, 'w-o', ...
     'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 4);
xlabel('Time (min)');
ylabel('Frequency (Hz)');
title('Respiration spectrogram  |  dominant freq searched in 1.5–5 Hz');
ylim(f_disp);
clim_percentile(ax1, S_resp(f_disp_r,:), [2 98]);

% --- Subplot 2: Vessel diameter spectrogram ---
ax2 = subplot(3,1,2);
imagesc(t_diam/60, f_diam(f_disp_d), ...
        10*log10(abs(S_diam(f_disp_d,:)).^2 + eps));
axis xy;
colormap(ax2, 'turbo');
cb = colorbar; cb.Label.String = 'Power (dB)';
hold on;
plot(t_min, dom_freq_diam, 'w-o', ...
     'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 4);
xlabel('Time (min)');
ylabel('Frequency (Hz)');
title('Vessel diameter spectrogram  |  dominant freq searched in 1.5–5 Hz');
ylim(f_disp);
clim_percentile(ax2, S_diam(f_disp_d,:), [2 98]);

% --- Subplot 3: Dominant frequency comparison ---
ax3 = subplot(3,1,3);
plot(t_min, dom_freq_resp, 'b-o', ...
     'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 4);
hold on;
plot(t_min, dom_freq_diam, 'r-o', ...
     'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 4);
xline(skip_sec/60, '--k', '14 min', 'LabelVerticalAlignment', 'bottom');
xlabel('Time (min)');
ylabel('Dominant frequency (Hz)');
title('Dominant frequency comparison  |  1.5–5 Hz band');
legend('Respiration', 'Vessel diameter', 'Location', 'best');
ylim(f_resp_band);
grid on;

linkaxes([ax1, ax2, ax3], 'x');

% =========================================================================
%% Local functions
% =========================================================================
function clim_percentile(ax, S, pct)
    P_db = 10*log10(abs(S(:)).^2 + eps);
    clim(ax, [prctile(P_db, pct(1)), prctile(P_db, pct(2))]);
end