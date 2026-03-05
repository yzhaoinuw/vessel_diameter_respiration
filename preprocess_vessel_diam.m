function diam_norm = preprocess_vessel_diam(img_data, baseline_window, show_plot)
    arguments
        img_data
        baseline_window = 10 * 60 % 10-minute window for moving average
        show_plot = false
    end
    
    %% pull out key variables (keep them as column vectors)
    t_frames = double(img_data.t_frames(:));
    umperpix = double(img_data.umperpix);
    % camera sampling rate estimate
    fs_cam = 1 / median(diff(t_frames));
    
    diam = double(img_data.diam(:));
    % convenience: diameter in microns
    diam_um = diam * umperpix;
    
    %% 1. Spike removal
    % Flag frames where diameter changes by more than 20% of baseline in one step
    baseline_diam = median(diam_um);  % robust baseline estimate
    spike_thresh  = 0.20 * baseline_diam;  % 20% of baseline — adjust if needed
    
    frame_deltas  = abs(diff(diam_um));
    is_spike      = [false; frame_deltas > spike_thresh];  % flag the arriving frame
    
    % Also flag NaNs or zeros (failed edge detection frames)
    is_spike = is_spike | ~isfinite(diam_um) | diam_um <= 0;
    
    fprintf('Spike frames detected: %d / %d (%.1f%%)\n', ...
            sum(is_spike), numel(diam_um), 100*mean(is_spike));
    
    % Interpolate over spikes
    t_good     = t_frames(~is_spike);
    d_good     = diam_um(~is_spike);
    diam_clean = interp1(t_good, d_good, t_frames, 'linear', 'extrap');
    
    %% 2. Detrending via moving average baseline
    % Window should be long enough to capture slow drift but not vasomotion
    % Rule of thumb: at least 10x the slowest physiological period of interest
    % Slowest vasomotion ~0.1 Hz → period ~10 s → window ~100 s
    win_sec   = baseline_window;
    win_frames = round(win_sec * fs_cam);
    if mod(win_frames, 2) == 0
        win_frames = win_frames + 1;  % make odd for symmetric window
    end
    
    baseline  = movmean(diam_clean, win_frames, 'omitnan');
    diam_det  = diam_clean - baseline;  % detrended (absolute change from baseline, microns)
    
    % 3. Normalise to percentage change from baseline
    diam_norm = (diam_det ./ baseline) * 100;  % units: % change from baseline
    
    % 4. Quick summary plot
    if show_plot    
        figure;
        subplot(3,1,1);
        plot(t_frames, diam_um, 'Color', [0.7 0.7 0.7]); hold on;
        plot(t_frames(is_spike), diam_um(is_spike), 'rx', 'MarkerSize', 4);
        plot(t_frames, baseline, 'k', 'LineWidth', 1.2);
        ylabel('Diameter (\mum)'); title('Raw + spikes + baseline'); 
        legend('raw', 'spikes', 'baseline');
        
        subplot(3,1,2);
        plot(t_frames, diam_det, 'Color', [0.27 0.51 0.71]);
        ylabel('\DeltaDiameter (\mum)'); title('Detrended (absolute)');
        
        subplot(3,1,3);
        plot(t_frames, diam_norm, 'Color', [0.18 0.63 0.45]);
        ylabel('\DeltaD/D_0 (%)'); title('Normalised (% change from baseline)');
        xlabel('Time (s)');
    end
end