function results = preprocess_resp(resp, fs, min_duration, band_resp, amp_win_sec, show_plots, verbose)
%PREPROCESS_RESP  Bandpass-filter respiration, find quiet-breathing sections,
%  detect inspiration peaks, and compute instantaneous phase/amplitude.
%
%  results = preprocess_resp(resp, fs)
%  results = preprocess_resp(resp, fs, min_duration, band_resp, amp_win_sec, show_plots, verbose)
%
%  INPUTS
%    resp           – raw respiration signal (column vector, e.g. thermocouple)
%    fs             – sampling rate (Hz), e.g. 1000
%    min_duration   – minimum quiet-section length in seconds  [default 5]
%    band_resp      – [lo hi] bandpass edges in Hz              [default [0.5 5]]
%    amp_win_sec    – sliding-window std window in seconds      [default 2]
%    show_plots     – show diagnostic plots                     [default false]
%    verbose        – print summary to command window           [default false]
%
%  OUTPUT  (struct with fields)
%    resp_filt      – bandpass-filtered respiration (same length as resp)
%    phase          – instantaneous phase from Hilbert transform (rad, −π to π)
%                     0 at inspiration peak, ±π at expiration trough
%    amplitude      – instantaneous amplitude envelope (Hilbert)
%    resp_std       – sliding-window std of raw signal (for diagnostics)
%    resp_std_thresh – threshold used to split quiet vs active breathing
%    good_sections  – table (idx_start, idx_end, duration) of quiet-breathing epochs
%                     idx_start/idx_end are sample indices into resp; caller
%                     converts to seconds via their own time vector
%    insp_peaks     – struct with fields:
%                       .locs  – indices into resp/resp_filt at each peak
%                       .amps  – resp_filt value at each peak
%
%  NOTES
%    - Quiet-breathing threshold = median(resp_std). Adjust resp_std_thresh
%      after calling if the default split is too aggressive/lenient.
%    - Inspiration peaks detected with MinPeakDistance = 0.15 s and
%      MinPeakProminence = 0.5 × median(|resp_filt|).
%    - Phase is computed via Hilbert transform on resp_filt. Phase = 0 at
%      positive peaks (inspiration), phase = ±π at troughs (expiration).

    arguments
        resp          (:,1) double
        fs            (1,1) double
        min_duration  (1,1) double = 5
        band_resp     (1,2) double = [0.5 5]
        amp_win_sec   (1,1) double = 2
        show_plots    (1,1) logical = false
        verbose       (1,1) logical = false
    end

    %% 1. Bandpass filter
    resp_filt = bandpass(resp, band_resp, fs);

    %% 2. Instantaneous phase & amplitude (Hilbert transform)
    analytic  = hilbert(resp_filt);
    phase     = angle(analytic);               % −π to π; 0 near positive peaks
    amplitude = abs(analytic);                 % envelope

    %% 3. Sliding-window std for amplitude estimation
    amp_win_samp = round(amp_win_sec * fs);
    resp_std = movstd(resp, amp_win_samp);

    %% 4. Identify quiet-breathing sections
    resp_std_thresh = median(resp_std);
    quiet_mask = resp_std < resp_std_thresh;

    % Contiguous quiet runs
    d_quiet    = diff([0; quiet_mask(:); 0]);
    starts_idx = find(d_quiet ==  1);
    ends_idx   = find(d_quiet == -1) - 1;
    durs       = (ends_idx - starts_idx + 1) / fs;

    keep = durs >= min_duration;
    good_sections = table( ...
        starts_idx(keep), ends_idx(keep), durs(keep), ...
        'VariableNames', {'idx_start', 'idx_end', 'duration'});

    %% 5. Detect inspiration peaks (within quiet sections only)
    min_peak_dist_samp = round(0.15 * fs);     % ~6.5 Hz max breathing rate
    min_prom = 0.5 * median(abs(resp_filt));

    all_locs = [];
    for s = 1:height(good_sections)
        i1 = good_sections.idx_start(s);
        i2 = good_sections.idx_end(s);
        seg = resp_filt(i1:i2);

        [~, locs] = findpeaks(seg, ...
            'MinPeakDistance',   min_peak_dist_samp, ...
            'MinPeakProminence', min_prom);

        all_locs = [all_locs; i1 - 1 + locs]; %#ok<AGROW>
    end

    insp_peaks.locs  = all_locs;               % indices into resp / resp_filt
    insp_peaks.amps  = resp_filt(all_locs);    % filtered amplitude at each peak

    %% 6. Pack output
    results.resp_filt       = resp_filt;
    results.phase           = phase;
    results.amplitude       = amplitude;
    results.resp_std        = resp_std;
    results.resp_std_thresh = resp_std_thresh;
    results.good_sections   = good_sections;
    results.insp_peaks      = insp_peaks;

    %% 7. Verbose summary
    if verbose
        fprintf('Resp preprocessing (fs = %g Hz, band = [%.1f %.1f] Hz)\n', ...
            fs, band_resp(1), band_resp(2));
        fprintf('  Quiet sections:      %d  (%.1f s total)\n', ...
            height(good_sections), sum(good_sections.duration));
        fprintf('  Inspiration peaks:   %d\n', numel(all_locs));
        fprintf('  Std threshold:       %.4f\n', resp_std_thresh);
    end

    %% 8. Diagnostic plots
    if show_plots
        plot_resp_overview(resp, resp_filt, resp_std, resp_std_thresh, ...
            good_sections, all_locs, fs);
    end
end

%% ========================================================================
%  HELPER: diagnostic plot
%  ========================================================================
function plot_resp_overview(resp, resp_filt, resp_std, thresh, sections, peak_locs, fs)
    t = (0:numel(resp)-1)' / fs;

    figure('Position', [100 100 1000 700], 'Color', 'w');

    % ---- Panel 1: raw + filtered resp with quiet sections shaded ---------
    ax1 = subplot(3,1,1);
    plot(t, resp, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4); hold on;
    plot(t, resp_filt, 'Color', [0.2 0.5 0.8], 'LineWidth', 0.8);
    for s = 1:height(sections)
        xs = t(sections.idx_start(s));
        xe = t(sections.idx_end(s));
        yl = ylim;
        patch([xs xe xe xs], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.8 0.9 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    ylabel('Resp (a.u.)');
    title('Raw (grey) + filtered (blue), quiet sections shaded');
    legend('Raw', 'Filtered', 'Quiet', 'Location', 'best');

    % ---- Panel 2: sliding-window std with threshold ----------------------
    ax2 = subplot(3,1,2);
    plot(t, resp_std, 'Color', [0.85 0.33 0.1], 'LineWidth', 0.8);
    yline(thresh, '--k', 'LineWidth', 1.2, 'Label', 'Threshold');
    ylabel('Resp std');
    title(sprintf('Sliding-window std (threshold = %.4f)', thresh));

    % ---- Panel 3: filtered resp with inspiration peaks marked ------------
    ax3 = subplot(3,1,3);
    plot(t, resp_filt, 'Color', [0.2 0.5 0.8], 'LineWidth', 0.8); hold on;
    plot(t(peak_locs), resp_filt(peak_locs), 'rv', 'MarkerSize', 4, ...
        'MarkerFaceColor', [0.85 0.2 0.2]);
    ylabel('Resp filtered (a.u.)');
    xlabel('Time (s)');
    title(sprintf('Inspiration peaks (%d detected)', numel(peak_locs)));

    linkaxes([ax1 ax2 ax3], 'x');
    sgtitle('Respiration Preprocessing Overview');
end