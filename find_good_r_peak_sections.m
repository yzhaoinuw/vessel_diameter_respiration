function good_sections = find_good_r_peak_sections(t_rpeaks, min_duration, show_plots, verbose)
    arguments
        t_rpeaks
        min_duration = 5
        show_plots = false
        verbose = false
    end
    
    rr_intervals = diff(t_rpeaks);  % in seconds
    
    %% 1. Filter R-peaks based on physiologically valid mouse RR intervals
    % Mouse HR: ~400–700 bpm → RR: ~85–150 ms; allow generous upper bound of 200 ms
    rr_max_ms = 200;   % adjust as needed
    rr_max    = rr_max_ms / 1000;  % convert to seconds
    
    % An R-peak is "bad" if the interval BEFORE or AFTER it is too long.
    % rr_intervals(i) = t_rpeaks(i+1) - t_rpeaks(i)
    bad_intervals = rr_intervals > rr_max;
    
    % Mark peaks bounding each bad interval as bad
    bad_peaks = false(size(t_rpeaks));
    bad_peaks(find(bad_intervals))     = true;  % peak i   (left bound)
    bad_peaks(find(bad_intervals) + 1) = true;  % peak i+1 (right bound)
    
    good_peaks = ~bad_peaks;
    t_rpeaks_good = t_rpeaks(good_peaks);
    
    if verbose
        fprintf('Total R-peaks:     %d\n', numel(t_rpeaks));
        fprintf('Good R-peaks:      %d\n', numel(t_rpeaks_good));
        fprintf('Removed R-peaks:   %d\n', numel(t_rpeaks) - numel(t_rpeaks_good));
    end
    
    %% 2. Identify consecutive good sections
    % A "good section" = run of consecutive surviving peaks in the original index list
    good_idx = find(good_peaks);  % indices into t_rpeaks of surviving peaks
    
    % Find breaks: wherever two surviving peaks are not adjacent in the original array
    breaks = find(diff(good_idx) > 1);  % positions in good_idx where a gap occurs
    
    % Build section start/end indices into good_idx
    sec_starts = [1;          breaks + 1];
    sec_ends   = [breaks;     numel(good_idx)];
    
    num_sections = numel(sec_starts);
    %fprintf('\nFound %d good consecutive section(s):\n', num_sections);
    k = 1;
    good_sections = struct();
    for s = 1:num_sections
        idx_range = good_idx(sec_starts(s)) : good_idx(sec_ends(s));  % original peak indices
        t_start = t_rpeaks(idx_range(1));
        t_end = t_rpeaks(idx_range(end));
        duration = t_end - t_start;
        if duration < min_duration
            continue
        end
    
        good_sections(k).peak_indices = idx_range;
        good_sections(k).peak_time = t_rpeaks(idx_range);
        good_sections(k).t_start      = t_start;
        good_sections(k).t_end        = t_end;
        good_sections(k).n_peaks      = numel(idx_range);
        good_sections(k).duration = duration;
        %fprintf('  Section %d: t = [%.2f, %.2f] s,  %d peaks, duration = %.2f s.\n', ...
        %        k, good_sections(k).t_start, good_sections(k).t_end, good_sections(k).n_peaks, good_sections(k).duration);
        k = k + 1; 
    end

    if show_plots
        plot_rr_interval_distribution(t_rpeaks)
        plot_good_r_peak_sections(t_rpeaks, good_sections, rr_max_ms)
    end
end

%% helper functions
function plot_rr_interval_distribution(t_rpeaks)
    % Compute and plot RR intervals
    rr_intervals = diff(t_rpeaks);  % in seconds
    rr_ms = rr_intervals * 1000;
    
    figure;
    histogram(rr_ms, 'BinWidth', 1, ...          % 2 ms bins — adjust to taste
              'FaceColor', [0.27 0.51 0.71], 'EdgeColor', 'white', 'LineWidth', 0.5);
    xlabel('RR Interval (ms)');
    ylabel('Count');
    title('Distribution of RR Intervals');
    %xlim([prctile(rr_ms, 0.05) * 0.5, prctile(rr_ms, 0.95) * 1.5]);    % tight x-axis with small padding
    xlim([50, 200]);    % tight x-axis with small padding
    
    xline(median(rr_ms), '--r', 'LineWidth', 1.5, 'Label', 'Median');
end

function plot_good_r_peak_sections(t_rpeaks, good_sections, rr_max_ms)
    % Visualise good sections on the RR interval trace
    arguments
        t_rpeaks
        good_sections
        rr_max_ms = 200
    end

    rr_intervals = diff(t_rpeaks);

    figure;
    plot(t_rpeaks(2:end), rr_intervals * 1000, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    hold on;
    
    num_good_sections = numel(good_sections);
    colors = lines(num_good_sections);
    for k = 1:num_good_sections
        idx = good_sections(k).peak_indices;
        % RR intervals within this section (interval i lives between peak i and i+1)
        rr_idx = idx(1:end-1);  % indices into rr_intervals = diff(t_rpeaks)
        plot(t_rpeaks(rr_idx + 1), rr_intervals(rr_idx) * 1000, ...
             'Color', colors(k,:), 'LineWidth', 1.5);
    end
    
    yline(rr_max_ms, '--r', 'LineWidth', 1.2, 'Label', sprintf('Max = %d ms', rr_max_ms));
    xlabel('Time (s)');
    ylabel('RR Interval (ms)');
    title('RR Intervals — Good Sections Highlighted');
    hold off;
end

function results = compute_HRV_metrics(good_sections)
    % Compute HRV metrics for each good section
    pnn_thresh_ms = 6;  % pNN threshold for mice (~6 ms vs 50 ms in humans)
    
    n_sections = numel(good_sections);
    results = struct();
    
    for s = 1:n_sections
    
        % Extract RR intervals for this section (in ms)
        idx  = good_sections(s).peak_indices;
        rr   = diff(t_rpeaks(idx)) * 1000;  % ms
    
        if numel(rr) < 2
            fprintf('Section %d: too few intervals, skipping.\n', s);
            continue
        end
    
        % --- Time-domain metrics ---
        meanRR  = mean(rr);
        meanHR  = 60000 / meanRR;
        SDNN    = std(rr);
        RMSSD   = sqrt(mean(diff(rr).^2));
        pNN     = 100 * mean(abs(diff(rr)) > pnn_thresh_ms);
        CV      = SDNN / meanRR * 100;
    
        % --- Store results ---
        results(s).section  = s;
        results(s).t_start  = good_sections(s).t_start;
        results(s).t_end    = good_sections(s).t_end;
        results(s).duration = good_sections(s).duration;
        results(s).n_beats  = numel(rr);
        results(s).meanRR   = meanRR;
        results(s).meanHR   = meanHR;
        results(s).SDNN     = SDNN;
        results(s).RMSSD    = RMSSD;
        results(s).pNN      = pNN;
        results(s).CV       = CV;
    
        %fprintf('Section %d (%.1f s, %d beats):\n', s, good_sections(s).duration, numel(rr));
        %fprintf('  meanRR=%.1f ms  HR=%.1f bpm  SDNN=%.2f ms  RMSSD=%.2f ms  pNN%d=%.1f%%  CV=%.2f%%\n\n', ...
        %        meanRR, meanHR, SDNN, RMSSD, pnn_thresh_ms, pNN, CV);
    end
end