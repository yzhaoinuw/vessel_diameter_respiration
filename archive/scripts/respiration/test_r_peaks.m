signal_data = load(".\data\signals.mat");
img_data = load(".\data\im.mat");

t_ecgoffset = double(img_data.t_ecgoffset);
t_ecg = double(signal_data.t_ECG(:));
% Align phys timestamps to camera clock
t_phys = t_ecg - t_ecgoffset;
good_r_peaks = signal_data.good_r_peaks; % good R peaks: conf >= 0.5
t_rpeaks = t_phys(good_r_peaks);

% Define time window: start of ECG up to the 10th R-peak
t_start = t_phys(1);
t_end   = t_rpeaks(10);

% Mask ECG signal to that window
mask = t_phys >= t_start & t_phys <= t_end;

figure;
plot(t_phys(mask), signal_data.ECG(mask), 'Color', [0.27 0.51 0.71], 'LineWidth', 0.8);
hold on;

% Vertical dashed lines at each of the first 10 R-peaks
for i = 1:10
    xline(t_rpeaks(i), '--r', 'LineWidth', 1.2, ...
          'Label', sprintf('R%d', i), 'LabelVerticalAlignment', 'bottom');
end

xlabel('Time (s)');
ylabel('Amplitude');
title('ECG — First 10 R-peaks');
xlim([t_start, t_end]);
hold off;