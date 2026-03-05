signal_data = load(".\data\signals.mat");
img_data = load(".\data\im.mat");

t_ecg = double(signal_data.t_ECG(:));
resp  = double(signal_data.resp(:));

t_frames = double(img_data.t_frames(:));
fs_cam = 1 / median(diff(t_frames));

t_ecgoffset = double(img_data.t_ecgoffset);
assert(t_ecgoffset > 0, 't_ecgoffset must be positive.')
% Align phys timestamps to camera clock
t_phys = t_ecg - t_ecgoffset;
fs_phys = 1 / median(diff(t_phys));

% ecg
good_r_peaks = signal_data.good_r_peaks; % good R peaks: conf >= 0.5
t_rpeaks = t_phys(good_r_peaks);
t_rpeaks = t_rpeaks(t_rpeaks>0 & t_rpeaks<t_frames(end));
good_sections = find_good_r_peak_sections(t_rpeaks);

diam_norm = preprocess_vessel_diam(img_data);
%% Extract matched segments of ECG, respiration, and vessel diameter
n_sections = numel(good_sections);

%% Resample respiration to camera rate
% t_frames is the authoritative camera time axis
resp_resampled = interp1(t_phys, resp, t_frames, 'linear', 'extrap');

%% Compare resp waveforms from low-std vs high-std sections
% Compute std for each section
resp_stds = zeros(n_sections, 1);
for s = 1:n_sections
    mask = t_phys >= good_sections(s).t_start & ...
           t_phys <= good_sections(s).t_end;
    resp_stds(s) = std(resp(mask));
end

% Pick examples from low, mid, high std
[~, idx_sorted] = sort(resp_stds);
picks = [idx_sorted(5), ...                              % low std
         idx_sorted(round(n_sections/2)), ...             % mid std
         idx_sorted(end-5)];                              % high std

figure('Position', [100 100 1200 600]);
for p = 1:3
    s = picks(p);
    t_start = good_sections(s).t_start;
    t_end   = good_sections(s).t_end;
    
    % Show up to 5 seconds
    t_show_end = min(t_end, t_start + 5);
    mask = t_phys >= t_start & t_phys <= t_show_end;
    
    subplot(3,1,p);
    plot(t_phys(mask), resp(mask), 'k', 'LineWidth', 0.8);
    ylabel('Resp');
    title(sprintf('Section %d | std = %.3f | t = %.0f–%.0f s', ...
           s, resp_stds(s), t_start, t_end));
end
xlabel('Time (s)');
sgtitle('Low / Mid / High resp variability sections');
