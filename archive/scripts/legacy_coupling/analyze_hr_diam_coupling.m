data_folder = ".\data\";
subject_name = "F154-dex40-19022025";
subject_folder = fullfile(data_folder, subject_name);

signal_data = load(fullfile(subject_folder, "signals.mat"));
img_data    = load(fullfile(subject_folder, "im.mat"));

t_ecg    = double(signal_data.t_ECG(:));
t_frames = double(img_data.t_frames(:));

fs_cam = 1 / median(diff(t_frames));
fs_ecg = 1 / median(diff(t_ecg));

t_ecgoffset = double(img_data.t_ecgoffset);
assert(t_ecgoffset > 0, 't_ecgoffset must be positive.')
t_ecg = t_ecg - t_ecgoffset;

%% 1. Vessel diameter preprocessing
t_skip   = 14 * 60;          % [s] skip first 14 minutes
baseline_window = 10 * 60; % 10-minute window for moving average
diam_norm = preprocess_vessel_diam(img_data, baseline_window, false);


% 2. ECG preprocessing
good_r_peaks = signal_data.good_r_peaks; % good R peaks: conf >= 0.5
t_rpeaks = t_ecg(good_r_peaks);
t_rpeaks = t_rpeaks(t_rpeaks>0 & t_rpeaks<t_frames(end));
good_sections = find_good_r_peak_sections(t_rpeaks);

n_sections = numel(good_sections);

%% Spectrogram of Normalized Vessel Diameter – Mouse
% Parameters
win_dur  = 60 * 10;
overlap_frac = 0.5;

nWin     = round(win_dur * fs_cam);
nOverlap = round(nWin * overlap_frac);
nFFT     = 2^nextpow2(nWin);

% Frequency axis limits
f_low  = 4; % can be as low as 3
f_high = 15;

% Spectrogram
[S, F, T] = spectrogram(diam_norm, hann(nWin), nOverlap, nFFT, fs_cam);

% Crop to frequency band of interest
fMask  = F >= f_low & F <= f_high;
F_crop = F(fMask);
S_crop = S(fMask, :);

% Power spectral density (dB re 1 %²/Hz)
PSD_dB = 10 * log10(abs(S_crop).^2 + eps);

%% Peak frequency per window (skip first 14 minutes)
tMask     = T >= t_skip;           % logical mask over time windows

[~, peakIdx]  = max(PSD_dB, [], 1);     % index of max PSD row for every window
f_peak        = F_crop(peakIdx);        % map indices → Hz
f_peak(~tMask) = NaN;                   % blank out the first 14 min


%% Heart Rate (Hz) in 1-minute bins, skipping the first 14 minutes

% --- Parameters ---
bin_dur  = win_dur * (1 - overlap_frac);               % [s] bin width

% Pool all good peak times into a single vector for easy binning
t_good = [];
for k = 1:numel(good_sections)
    t_good = [t_good, good_sections(k).peak_time(:)'];
end
t_good = sort(t_good);

% --- Define bin edges starting at t_skip ---
t_start_all = t_skip;
t_end_all   = floor((max(t_good) - t_skip) / bin_dur) * bin_dur + t_skip;
bin_edges   = t_start_all : bin_dur : t_end_all + bin_dur;
bin_centers = bin_edges(1:end-1) + bin_dur/2;   % centre of each bin [s]
n_bins      = numel(bin_centers);

% --- Compute HR in Hz per bin ---
hr_hz = nan(1, n_bins);

for b = 1:n_bins
    % Good peaks that fall inside this bin
    in_bin = t_good >= bin_edges(b) & t_good < bin_edges(b+1);
    t_bin  = t_good(in_bin);

    if numel(t_bin) < 2          % need at least 2 peaks for an RR interval
        continue
    end

    rr_mean   = mean(diff(t_bin));   % mean RR interval [s]
    hr_bin = 1 / rr_mean;
    hr_hz(b)  = 1 / rr_mean;        % HR in Hz (beats/s)
end

%% Plot
% Spectrogram
figure('Color','w','Position',[100 100 1100 450]);

imagesc(T/60, F_crop, PSD_dB);
axis xy;
colormap(turbo);
cb = colorbar;
cb.Label.String = 'PSD (dB re 1 \%^{2}/Hz)';

hold on;

% Peak frequency line (thick dashed white)
plot(T/60, f_peak, '--ok', ...
     'LineWidth', 2, 'MarkerSize', 4, ...
     'MarkerFaceColor', 'k');

% Heart rate overlay
plot(bin_centers/60, hr_hz, '--ow', ...
     'LineWidth', 2, 'MarkerSize', 4, ...
     'MarkerFaceColor', 'w');

% Axes
xlabel('Time (min)');
ylabel('Frequency (Hz)');
title("Vessel Diam Spectrogram – Mouse " + subject_name + "(\DeltaD / \langle D \rangle, %)");
ylim([f_low, f_high]);
clim([prctile(PSD_dB(:), 5), prctile(PSD_dB(:), 99)]);

set(findall(gcf, 'Type','text'), 'Color','k');
set(gca, 'XColor','k', 'YColor','k');
cb.Label.Color = 'k';
cb.Color = 'k';              % colorbar tick labels
lg = legend({'Peak freq.', 'Heart rate (Hz)'}, ...
       'Location','northeast', 'Color','none', 'EdgeColor','none');
lg.TextColor = 'k';

%%
% Interpolate HR onto spectrogram time axis (post-skip only)
tMask_T      = T >= t_skip;
T_valid      = T(tMask_T);
hr_interp    = interp1(bin_centers, hr_hz, T_valid, 'linear');
f_peak_valid = f_peak(tMask_T);

% Remove NaNs from either vector
hr_interp    = hr_interp(:);
f_peak_valid = f_peak_valid(:);
valid = ~isnan(hr_interp) & ~isnan(f_peak_valid);

[r, pval] = corr(f_peak_valid(valid), hr_interp(valid), 'Type', 'Spearman');
fprintf('Spearman r = %.3f,  p = %.4f\n', r, pval);

%% Figure 2 — Frequency Ratio
ratio = f_peak_valid ./ hr_interp;

figure('Color','w', 'Position', [100 100 900 350]);
plot(T_valid/60, ratio, '-k', 'LineWidth', 1.5);
yline(1, '--r', 'LineWidth', 1.5, 'Label', '1:1 lock', ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
xlabel('Time (min)');
ylabel('f_{peak} / HR');
title('Frequency Ratio — Vessel Peak / Heart Rate');
ylim([0.5 1.5]);
grid on; box off;
set(gca, 'Color','w', 'XColor','k', 'YColor','k');
set(get(gca,'XLabel'), 'Color','k');
set(get(gca,'YLabel'), 'Color','k');
set(get(gca,'Title'),  'Color','k');

%% Figure 3 — Coherence

% Build evenly-sampled RR interval series
rr_intervals = diff(t_good);
t_rr         = t_good(2:end);

% Resample both signals to a common grid at fs_cam
t_common  = t_rr(1) : 1/fs_cam : t_rr(end);
rr_interp = interp1(t_rr, rr_intervals, t_common, 'pchip');
diam_trim = interp1(t_frames, diam_norm,  t_common, 'linear');

% Crop to post-skip period
post_skip = t_common >= t_skip;
rr_coh    = rr_interp(post_skip);
diam_coh  = diam_trim(post_skip);

% Coherence
win_coh   = round(60 * fs_cam);
nOver_coh = round(win_coh * 0.5);
nFFT_coh  = 2^nextpow2(win_coh);
[Cxy, F_coh] = mscohere(diam_coh, rr_coh, hann(win_coh), nOver_coh, nFFT_coh, fs_cam);

% Significance threshold
n_windows  = 2 * floor(numel(diam_coh) / win_coh) - 1;
alpha      = 0.05;
coh_thresh = 1 - alpha^(1 / (n_windows - 1));

% Median HR (within band only, for the xline)
med_hr = median(hr_hz(hr_hz >= f_low & hr_hz <= f_high), 'omitnan');

figure('Color','w', 'Position', [100 100 900 350]);
plot(F_coh, Cxy, 'k', 'LineWidth', 1.5);
xline(med_hr, '--r', 'LineWidth', 1.5, ...
      'Label', sprintf('Median HR = %.1f Hz', med_hr), ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
yline(coh_thresh, ':k', 'LineWidth', 1.2, ...
      'Label', sprintf('p<%.2f (n=%d windows)', alpha, n_windows), ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
xlabel('Frequency (Hz)');
ylabel('Magnitude-Squared Coherence');
title('Vessel Diameter – RR Interval Coherence');
xlim([f_low, f_high]);
ylim([0 1]);
grid on; box off;
set(gca, 'Color','w', 'XColor','k', 'YColor','k');
set(get(gca,'XLabel'), 'Color','k');
set(get(gca,'YLabel'), 'Color','k');
set(get(gca,'Title'),  'Color','k');
set(findall(gcf, 'Type','text'), 'Color','k');