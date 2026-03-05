signal_data = load(".\data\signals.mat");
img_data    = load(".\data\im.mat");

%% pull out key variables (keep them as column vectors)
t_frames = double(img_data.t_frames(:));
%t_vect = double(img_data.t_vect(:));   % exists but t_frames is usually the one to trust
umperpix = double(img_data.umperpix);
% camera sampling rate estimate
fs_cam = 1 / median(diff(t_frames));

diam = double(img_data.diam(:));
% convenience: diameter in microns
diam_um = diam * umperpix;

%% =========================
%  SECTION 2: Vessel preprocessing (robust)
%  Assumes Section 1 already created:
%    diam_um, t_frames, fs_cam
% ==========================

% ---- create analysis directory ----
analysis_dir = ".\analysis";
if ~exist(analysis_dir, "dir")
    mkdir(analysis_dir);
end

% ---- detrend baseline (10-min moving mean) ----
baseline_minutes = 10;
win = round(baseline_minutes * 60 * fs_cam);
win = max(win, 3);

diam_avg  = movmean(diam_um, win, "Endpoints", "shrink");
diam_detr = diam_um - diam_avg;

% ---- fill missing values BEFORE filtering ----
diam_detr = fillmissing(diam_detr, "linear");
diam_detr = fillmissing(diam_detr, "nearest"); % edges

% ---- define bands ----
Fc_slow_hi = 0.30;      % Hz (vasomotion upper cutoff)
band_card  = [5 15];    % Hz (cardiac band)
order = 4;

% ---- Slow component: lowpass (robust, matches legacy intent) ----
% This avoids unstable near-DC bandpass. The 10-min detrend already
% removes ultra-slow drift, so lowpass captures ~0.01–0.3-ish behavior.
diam_slow = lowpass(diam_detr, Fc_slow_hi, fs_cam, ...
                    "ImpulseResponse","iir", "Steepness", 0.85);

% ---- Cardiac component: bandpass (stable) ----
[b_card, a_card] = butter(order, band_card/(fs_cam/2), "bandpass");
diam_card = filtfilt(b_card, a_card, diam_detr);

% ---- Cardiac envelope ----
diam_card_env = abs(hilbert(diam_card));

% ---- QC plot: 10 minutes ----
%Tshow = 10*60;
%t0 = t_frames(1);

Tshow = 1*60;
t0 = 20*60;

mask = t_frames >= t0 & t_frames <= (t0 + Tshow);

figure("Name","QC: vessel preprocessing (1 min sample)","Color","w");
tiledlayout(4,1,"Padding","compact","TileSpacing","compact");

nexttile;
plot(t_frames(mask), diam_um(mask));
ylabel("diam (\mum)");
title("Raw diameter (microns)");

nexttile;
plot(t_frames(mask), diam_detr(mask));
ylabel("detr (\mum)");
title(sprintf("Detrended (baseline = %d min movmean)", baseline_minutes));

nexttile;
plot(t_frames(mask), diam_slow(mask));
ylabel("slow (\mum)");
title(sprintf("Slow = lowpass(detr, %.2f Hz)", Fc_slow_hi));

nexttile;
plot(t_frames(mask), diam_card(mask)); hold on;
plot(t_frames(mask), diam_card_env(mask));
ylabel("card (\mum)");
xlabel("time (s)");
title(sprintf("Cardiac bandpass [%.1f %.1f] Hz + envelope", band_card(1), band_card(2)));
legend(["diam\_card","env"], "Location","best");

% ---- numeric sanity ----
fprintf("diam_slow: NaNs=%d, std=%.6f um\n", sum(isnan(diam_slow)), std(diam_slow,"omitnan"));
fprintf("diam_card: NaNs=%d, std=%.6f um\n", sum(isnan(diam_card)), std(diam_card,"omitnan"));

% ---- save stage 2 outputs ----
save(fullfile(analysis_dir,"stage02_vessel_preproc.mat"), ...
    "baseline_minutes","Fc_slow_hi","band_card", ...
    "diam_avg","diam_detr","diam_slow","diam_card","diam_card_env", ...
    "-v7.3");

fprintf("Saved %s\n", fullfile(analysis_dir,"stage02_vessel_preproc.mat"));

%% =========================
%  SECTION 3: Heart + Resp preprocessing (aligned to camera time)
%  Assumes Section 1 already created:
%    signal_data, img_data, t_frames, fs_cam
%  Outputs:
%    t_phys, resp_cam, resp_filt_cam, resp_phase_cam
%    tR, tHR, HR_inst_bpm, HR_cam
% ==========================


% ---- Pull physio signals ----
t_ECG = double(signal_data.t_ECG(:));
resp  = double(signal_data.resp(:));

% ---- Alignment offset (camera start relative to phys start) ----
% In your example: t_ecgoffset is stored in im.mat (img_data)
t_ecgoffset = double(img_data.t_ecgoffset);

% Align phys timestamps to camera clock
t_phys = t_ECG - t_ecgoffset;

% ---- R-peak indices (MATLAB 1-based) ----
pk_inds_keep = double(signal_data.pk_inds_keep);  % Nx2

i1 = pk_inds_keep(:,1);
i2 = pk_inds_keep(:,2);

% sanity bounds
assert(all(i1>=1 & i1<=numel(t_phys)), "pk_inds_keep(:,1) out of bounds");
assert(all(i2>=1 & i2<=numel(t_phys)), "pk_inds_keep(:,2) out of bounds");

% Beat times and instantaneous HR
tR  = t_phys(i1);                  % R-peak times (sec, aligned)
tHR = 0.5*(t_phys(i1)+t_phys(i2)); % timestamps for HR_inst (sec)
RR  = t_phys(i2) - t_phys(i1);     % seconds
RR(RR<=0) = NaN;
HR_inst_bpm = 60 ./ RR;

% Optional cleanup: remove implausible HR (awake mouse usually ~250–900 bpm)
bad = HR_inst_bpm < 150 | HR_inst_bpm > 1200;
HR_inst_bpm(bad) = NaN;

% ---- Interpolate HR onto camera timeline ----
good_hr = isfinite(tHR) & isfinite(HR_inst_bpm);
HR_cam = interp1(tHR(good_hr), HR_inst_bpm(good_hr), t_frames, "linear", "extrap");

% Mild smoothing helps for later Granger / slow coupling
HR_cam = smoothdata(HR_cam, "movmean", round(2*fs_cam));  % 2-sec smoothing

% ---- Interpolate respiration onto camera timeline ----
resp_cam = interp1(t_phys, resp, t_frames, "linear", "extrap");

% ---- Resp filtering + phase on camera grid ----
% Start band for awake mouse; we can tune if needed
band_resp = [0.5 5.0]; % Hz

% Use bandpass (stable) then Hilbert phase
resp_filt_cam = bandpass(resp_cam, band_resp, fs_cam);
resp_phase_cam = angle(hilbert(resp_filt_cam));  % [-pi, pi]

% ---- QC plot (first 10 minutes) ----
%Tshow = 10*60;
%t0 = t_frames(1);

Tshow = 1*60;
t0 = 20*60;
mask = t_frames >= t0 & t_frames <= (t0 + Tshow);

figure("Name","QC: phys preprocessing (1 min sample)","Color","w");
tiledlayout(4,1,"Padding","compact","TileSpacing","compact");

nexttile;
plot(t_frames(mask), HR_cam(mask));
ylabel("HR (bpm)");
title("HR interpolated to camera timeline (smoothed)");

nexttile;
plot(t_frames(mask), resp_cam(mask));
ylabel("Resp (a.u.)");
title("Resp interpolated to camera timeline");

nexttile;
plot(t_frames(mask), resp_filt_cam(mask));
ylabel("Resp filt");
xlabel("time (s)");
title(sprintf("Resp bandpass [%.1f %.1f] Hz + phase computed", band_resp(1), band_resp(2)));

nexttile;
plot(t_frames(mask), resp_phase_cam(mask));
ylabel("Card phase");
xlabel("time (s)");
title("bandpass + Hilbert phase");

% ---- numeric sanity ----
fprintf("HR_cam: median=%.1f bpm, IQR=%.1f bpm\n", median(HR_cam,"omitnan"), iqr(HR_cam));
fprintf("resp_filt_cam: std=%.6f\n", std(resp_filt_cam,"omitnan"));

% ---- save stage 3 outputs ----
save(fullfile(analysis_dir,"stage03_phys_preproc.mat"), ...
    "t_ecgoffset","t_phys", ...
    "tR","tHR","HR_inst_bpm","HR_cam", ...
    "resp","t_ECG","resp_cam","resp_filt_cam","resp_phase_cam","band_resp", ...
    "-v7.3");

fprintf("Saved %s\n", fullfile(analysis_dir,"stage03_phys_preproc.mat"));

%% =========================
%  SECTION 4A: Uniform analysis grid
% ==========================

%analysis_dir = ".\analysis";

% analysis sampling rate (match camera)
fsA = fs_cam;

tA = (t_frames(1) : 1/fsA : t_frames(end))';
N = numel(tA);

% Interpolate signals onto uniform grid
diam_card_A     = interp1(t_frames, diam_card,     tA, "linear", "extrap");
diam_card_env_A = interp1(t_frames, diam_card_env, tA, "linear", "extrap");
diam_slow_A     = interp1(t_frames, diam_slow,     tA, "linear", "extrap");

resp_filt_A  = interp1(t_frames, resp_filt_cam,  tA, "linear", "extrap");
HR_A         = interp1(t_frames, HR_cam,         tA, "linear", "extrap");

fprintf("Uniform grid: %d samples @ %.1f Hz\n", N, fsA);

%% =========================
%  SECTION 4B-1: Heartbeat-triggered vessel response
% ==========================

pre_sec  = 0.25;
post_sec = 0.50;

pre_samp  = round(pre_sec  * fs_cam);
post_samp = round(post_sec * fs_cam);

% Map R-peaks to camera indices
idxR = interp1(t_frames, 1:numel(t_frames), tR, "nearest");
idxR = idxR(idxR-pre_samp > 0 & idxR+post_samp <= numel(diam_card));

tau = (-pre_samp:post_samp) / fs_cam;

ETA = zeros(numel(idxR), numel(tau));
for k = 1:numel(idxR)
    i = idxR(k);
    ETA(k,:) = diam_card(i-pre_samp:i+post_samp);
end

ETA_mean = mean(ETA,1,'omitnan');
ETA_sem  = std(ETA,[],1,'omitnan') / sqrt(size(ETA,1));

% amplitude + latency
[~,pk_idx] = max(ETA_mean(tau>=0));
pk_latency = tau(tau>=0);
pk_latency = pk_latency(pk_idx);

ETA_amp = prctile(ETA_mean,95) - prctile(ETA_mean,5);

fprintf("ETA: %d beats, amp=%.3f um, latency=%.3f s\n", ...
        size(ETA,1), ETA_amp, pk_latency);

% Plot
figure("Name","Heartbeat-triggered vessel response","Color","w");
plot(tau, ETA_mean, "LineWidth",1.5); hold on;
fill([tau fliplr(tau)], ...
     [ETA_mean+ETA_sem fliplr(ETA_mean-ETA_sem)], ...
     [0.8 0.8 0.8], "EdgeColor","none");
xline(0,'--');
xlabel("Time from R-peak (s)");
ylabel("Vessel pulsation (\mum)");
title("Heartbeat-triggered vessel diameter");

%% =========================
%  SECTION 4C: Resp phase modulation of cardiac envelope
% ==========================

nBins = 18;
edges = linspace(-pi, pi, nBins+1);
centers = (edges(1:end-1)+edges(2:end))/2;

mean_env = nan(1,nBins);

for b = 1:nBins
    mask = resp_phase_cam >= edges(b) & resp_phase_cam < edges(b+1);
    mean_env(b) = mean(diam_card_env(mask),'omitnan');
end

mod_index = (max(mean_env)-min(mean_env)) / mean(mean_env,'omitnan');

fprintf("Resp→cardiac modulation index = %.3f\n", mod_index);

figure("Name","Resp phase modulation","Color","w");
plot(centers, mean_env, '-o');
xlabel("Respiration phase (rad)");
ylabel("Cardiac envelope (\mum)");
title(sprintf("Resp-phase modulation (MI=%.3f)", mod_index));

%% =========================
%  SECTION 4D-1: Coherence (slow vessel ↔ HR)
% ==========================

win_sec = 60;
nwin = round(win_sec * fsA);

[Cohr,f] = mscohere(diam_slow_A, HR_A, hamming(nwin), [], [], fsA);

figure("Name","Coherence: slow vessel ↔ HR","Color","w");
plot(f, Cohr);
xlim([0 1]);
xlabel("Frequency (Hz)");
ylabel("Coherence");
title("Coherence between slow vessel motion and HR");

%% =========================
%  SECTION 4D-2: Coherence (slow vessel ↔ respiration)
% ==========================

[Cresp,f] = mscohere(diam_slow_A, resp_filt_A, hamming(nwin), [], [], fsA);

figure("Name","Coherence: slow vessel ↔ respiration","Color","w");
plot(f, Cresp);
xlim([0 1]);
xlabel("Frequency (Hz)");
ylabel("Coherence");
title("Coherence between slow vessel motion and respiration");

%% =========================
%  SECTION 5: Directionality (Windowed Conditional Granger)
% ==========================

analysis_dir = ".\analysis";
if ~exist(analysis_dir, "dir")
    mkdir(analysis_dir);
end

% ---- Build low-rate uniform grid for VAR/Granger ----
fsG = 10;  % Hz (good trade-off)
tG = (t_frames(1) : 1/fsG : t_frames(end))';

% Resample onto tG
diamG = interp1(t_frames, diam_slow,     tG, "linear", "extrap");
hrG   = interp1(t_frames, HR_cam,        tG, "linear", "extrap");
respG = interp1(t_frames, resp_filt_cam, tG, "linear", "extrap");

% Clean NaNs (VAR hates NaNs)
diamG = fillmissing(diamG, "linear"); diamG = fillmissing(diamG, "nearest");
hrG   = fillmissing(hrG,   "linear"); hrG   = fillmissing(hrG,   "nearest");
respG = fillmissing(respG, "linear"); respG = fillmissing(respG, "nearest");

% Standardize (z-score) globally to keep scale comparable
diamGz = zscore(diamG);
hrGz   = zscore(hrG);
respGz = zscore(respG);

% Stack: columns are variables
% 1 = vessel slow, 2 = HR, 3 = respiration
X = [diamGz, hrGz, respGz];

% ---- Windowing parameters ----
win_sec  = 120;   % window length (s)
step_sec = 30;    % step (s)

win  = round(win_sec  * fsG);
step = round(step_sec * fsG);

% VAR lag order (in samples at fsG)
% Rule of thumb: keep modest; too high overfits.
pLag = 10;  % = 1 second history at 10 Hz

nW = floor((size(X,1) - win)/step) + 1;
t_mid = nan(nW,1);

% Directions tested (src -> dst)
pairs = [
    2 1  % HR  -> vessel
    3 1  % Resp-> vessel
    1 2  % vessel -> HR
    1 3  % vessel -> Resp
    2 3  % HR -> Resp
    3 2  % Resp -> HR
];

nPairs = size(pairs,1);
Fstat = nan(nW, nPairs);
pval  = nan(nW, nPairs);

% ---- Granger per window using gctest ----
% gctest(x,y,'NumLags',p) tests whether y Granger-causes x
% in a bivariate setting. For conditional Granger with >2 vars,
% gctest supports multivariate by passing a matrix and specifying
% response/predictor columns.
%
% We'll use gctest with matrix form:
%   [h,p,stat] = gctest(Y, X, 'NumLags', pLag)
% where Y is response series and X includes predictors.
% Conditional tests: include all variables in X, then test exclusion.

for w = 1:nW
    a = (w-1)*step + 1;
    b = a + win - 1;

    seg = X(a:b, :);
    t_mid(w) = mean(tG([a b]));

    % Demean within window (helps stationarity)
    seg = seg - mean(seg, 1);

    % For each direction src->dst, test whether src adds predictive power
    for k = 1:nPairs
        src = pairs(k,1);
        dst = pairs(k,2);

        % Response is dst column
        Y = seg(:, dst);

        % Full predictor set: all variables
        Xfull = seg;

        % gctest syntax: gctest(Y, X, ...) tests if X Granger-causes Y
        % BUT we want to test src -> dst conditionally on the others.
        % We do that by comparing:
        %   full model: all predictors
        %   reduced model: predictors excluding src
        %
        % MATLAB's gctest returns results for whether X causes Y; for conditional,
        % we compute F-test via nested regressions manually.

        % Build lagged design matrices
        [Yt, Xlag_full] = build_var_design(Y, Xfull, pLag);

        % Reduced predictors: remove src column from Xfull before lagging
        Xred = Xfull;
        Xred(:, src) = [];  % drop src variable
        [~,  Xlag_red]  = build_var_design(Y, Xred,  pLag);

        % Fit OLS: Yt = [1 Xlag]*beta
        % Full
        beta_full = Xlag_full \ Yt;
        resid_full = Yt - Xlag_full * beta_full;
        RSS_full = sum(resid_full.^2);

        % Reduced
        beta_red = Xlag_red \ Yt;
        resid_red = Yt - Xlag_red * beta_red;
        RSS_red = sum(resid_red.^2);

        % Compute F-stat
        df1 = size(Xlag_full,2) - size(Xlag_red,2);
        df2 = numel(Yt) - size(Xlag_full,2);
        F = ((RSS_red - RSS_full)/df1) / (RSS_full/df2);

        % p-value
        p = 1 - fcdf(F, df1, df2);

        Fstat(w,k) = F;
        pval(w,k) = p;
    end
end

% ---- Plot -log10(p) ----
labels = {
    "HR \rightarrow Vessel"
    "Resp \rightarrow Vessel"
    "Vessel \rightarrow HR"
    "Vessel \rightarrow Resp"
    "HR \rightarrow Resp"
    "Resp \rightarrow HR"
};

figure("Name","Windowed conditional Granger (-log10 p)","Color","w");
hold on;
for k = 1:nPairs
    plot(t_mid, -log10(pval(:,k) + 1e-300), "LineWidth", 1.2);
end
yline(-log10(0.05), "--", "p=0.05");
xlabel("Time (s)");
ylabel("-log10(p)");
title(sprintf("Windowed conditional Granger (win=%ds, step=%ds, fs=%gHz, lag=%d)", ...
      win_sec, step_sec, fsG, pLag));
legend(labels, "Location","best", "FontSize", 9);
grid on;

% ---- Save outputs ----
save(fullfile(analysis_dir,"stage05_granger.mat"), ...
    "fsG","tG","t_mid","pairs","labels","Fstat","pval", ...
    "win_sec","step_sec","pLag", ...
    "-v7.3");

fprintf("Saved %s\n", fullfile(analysis_dir,"stage05_granger.mat"));

%% -------- helper function (keep at bottom of script) --------
function [Yt, Xlag] = build_var_design(Y, X, pLag)
% Build lagged regression design matrix for predicting Y from lags of X.
% Returns:
%   Yt   = Y(pLag+1:end)
%   Xlag = [ones, X(t-1..t-pLag for all columns)]

    Y = Y(:);
    T = size(X,1);
    K = size(X,2);

    Yt = Y(pLag+1:end);

    % Build lag matrix: [X(t-1) ... X(t-pLag)]
    Xlag_no_const = zeros(T-pLag, K*pLag);

    col = 1;
    for lag = 1:pLag
        Xlag_no_const(:, col:col+K-1) = X(pLag+1-lag : T-lag, :);
        col = col + K;
    end

    Xlag = [ones(size(Xlag_no_const,1),1), Xlag_no_const];
end

%%
RR_all = diff(tR);
mean_RR = mean(RR_all);
frames_per_beat = mean_RR * fs_cam;
Nphase = round(frames_per_beat);

fprintf("Mean frames per beat: %.1f\n", frames_per_beat);

% Inputs you already have:
% t_frames   : time stamps of camera frames (s), size [Nframes x 1]
% diam_card  : vessel diameter trace (same length as t_frames), size [Nframes x 1]
% tR         : R-peak times (s), size [Nbeats x 1]
% fs_cam     : camera sampling rate (Hz)  (optional for plotting, not required)

% --------- settings ----------
%Nphase = 200;          % number of samples per normalized cycle (0..1)
phase_grid = linspace(0,1,Nphase);

% Ensure column vectors
t_frames  = t_frames(:);
diam_card = diam_card(:);
tR        = tR(:);

% Optional: only keep R-peaks within frame time range
tR = tR(tR >= t_frames(1) & tR <= t_frames(end));

% We'll need beat k and beat k+1, so usable beats are 1..end-1
nBeats = numel(tR) - 1;

% Preallocate
ETA_phase = nan(nBeats, Nphase);
RR = nan(nBeats,1);

% --------- main loop ----------
for k = 1:nBeats
    t0 = tR(k);
    t1 = tR(k+1);
    RR(k) = t1 - t0;

    % Skip crazy intervals (artifact / missed detection)
    if RR(k) <= 0 || RR(k) > 1
        continue
    end

    % Extract the portion of the diameter trace within this beat
    inBeat = (t_frames >= t0) & (t_frames < t1);
    if nnz(inBeat) < 5
        continue
    end

    tb = t_frames(inBeat);
    xb = diam_card(inBeat);

    % Convert to phase in [0,1)
    ph = (tb - t0) ./ RR(k);

    % Resample onto common phase grid
    % Use 'linear' interpolation; 'pchip' is also fine if you prefer smoother
    ETA_phase(k,:) = interp1(ph, xb, phase_grid, "linear", nan);

    % Optional baseline: subtract mean of early phase (e.g., first 10%)
    % ETA_phase(k,:) = ETA_phase(k,:) - mean(ETA_phase(k, phase_grid <= 0.1), "omitnan");
end

% Keep only valid rows
valid = all(isfinite(ETA_phase), 2);
ETA_phase = ETA_phase(valid,:);
RR = RR(valid);

% --------- summary stats ----------
ETA_mean = mean(ETA_phase, 1, "omitnan");
ETA_sem  = std(ETA_phase, 0, 1, "omitnan") ./ sqrt(size(ETA_phase,1));

% Optional amplitude metrics on the mean waveform
amp_p2p = max(ETA_mean) - min(ETA_mean);
fprintf("Phase-warped ETA: %d beats used, mean RR=%.3f s, p2p amp=%.4f (units of diam_card)\n", ...
    size(ETA_phase,1), mean(RR,"omitnan"), amp_p2p);

% --------- plot mean ± SEM ----------
figure("Name","Heartbeat phase-warped vessel waveform","Color","w");
hold on;

% Shaded SEM band
fill([phase_grid fliplr(phase_grid)], ...
     [ETA_mean + ETA_sem, fliplr(ETA_mean - ETA_sem)], ...
     [0.85 0.85 0.85], "EdgeColor","none");

% Mean curve
plot(phase_grid, ETA_mean, "LineWidth", 1.8);

xlabel("Cardiac phase (0 \rightarrow 1)");
ylabel("Vessel diameter (same units as diam\_card)");
title("Phase-warped heartbeat-triggered vessel waveform");
xlim([0 1]);
box off;
