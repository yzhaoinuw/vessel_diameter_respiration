%% Explore respiration across different time windows
% Pick windows from different parts of the recording
signal_data = load(".\data\signals.mat");
img_data = load(".\data\im.mat");

t_ecg = double(signal_data.t_ECG(:));
resp  = double(signal_data.resp(:));

t_frames = double(img_data.t_frames(:));
fs_cam = 1 / median(diff(t_frames));
fs_ecg = 1 / median(diff(t_ecg));
t_ecgoffset = double(img_data.t_ecgoffset);
t_phys = t_ecg - t_ecgoffset;

%% Explore respiration signal
figure('Position', [100 100 1200 600]);

% Raw signal overview
subplot(3,1,1);
plot(t_phys, resp, 'k', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Resp (raw units)');
title('Full respiration signal');

% Zoom in to ~5 seconds to see individual breaths
subplot(3,1,2);
t_zoom = [10 15];  % adjust if this isn't a good window
mask = t_phys >= t_zoom(1) & t_phys <= t_zoom(2);
plot(t_phys(mask), resp(mask), 'k', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Resp (raw units)');
title('Zoomed: 5 seconds');

% Power spectrum
subplot(3,1,3);
[pxx, f] = pwelch(resp, round(fs_ecg*5), [], [], fs_ecg);
plot(f, 10*log10(pxx), 'k', 'LineWidth', 1);
xlim([0 20]);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Respiration power spectrum');

% Basic stats
fprintf('Resp range: [%.2f, %.2f]\n', min(resp), max(resp));
fprintf('Resp mean: %.4f, std: %.4f\n', mean(resp), std(resp));
fprintf('Resp sampling rate: %.1f Hz\n', fs_ecg);
fprintf('Duration: %.1f s\n', (numel(resp)-1)/fs_ecg);
%%
windows = [
    10   15;    % early (what we already saw)
    500  505;   % mid-early
    2000 2005;  % transition zone
    3500 3505;  % settled period
    6000 6005;  % mid-recording
    9000 9005;  % late recording
];


n_win = size(windows, 1);
figure('Position', [100 100 1400 800]);

for w = 1:n_win
    t_zoom = windows(w, :);
    mask = t_phys >= t_zoom(1) & t_phys <= t_zoom(2);

    % Waveform
    subplot(n_win, 2, (w-1)*2 + 1);
    plot(t_phys(mask), resp(mask), 'k', 'LineWidth', 0.8);
    ylabel('Resp');
    title(sprintf('%.0f – %.0f s', t_zoom(1), t_zoom(2)));
    if w == n_win, xlabel('Time (s)'); end

    % Local power spectrum
    subplot(n_win, 2, (w-1)*2 + 2);
    seg = resp(mask);
    [pxx, f] = pwelch(seg, round(fs_ecg*1), [], [], fs_ecg);
    plot(f, 10*log10(pxx), 'k', 'LineWidth', 1);
    xlim([0 20]);
    ylabel('dB');
    title(sprintf('Spectrum: %.0f–%.0f s', t_zoom(1), t_zoom(2)));
    if w == n_win, xlabel('Frequency (Hz)'); end
end

sgtitle('Respiration signal across recording');