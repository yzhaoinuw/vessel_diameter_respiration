# Cardio-Respiratory-Vascular Coupling Analysis
## Directionality and Coupling Study in Mice

*Technical Analysis Report*

---

## Executive Summary

This report presents the analysis of coupling relationships between blood vessel pulsatility, cardiac activity (ECG/heart rate), and respiration in mice. The analysis pipeline processes synchronized physiological recordings to quantify both the strength and directionality of interactions between these cardiovascular and respiratory systems.

**Key findings** indicate weak linear coupling between slow vessel dynamics and cardio-respiratory signals, modest but detectable respiratory gating of cardiac pulsatility (2.7% modulation), and Granger causality results that require methodological refinement due to apparent artifacts in portions of the recording.

---

## 1. Introduction and Scientific Background

### 1.1 Physiological Context

Blood vessel diameter in living tissue exhibits fluctuations across multiple timescales, driven by distinct physiological mechanisms:

**Cardiac pulsation (5–15 Hz):** Passive mechanical distension of vessel walls with each heartbeat as the pulse pressure wave arrives. In mice with heart rates of 400–700 bpm (6.7–11.7 Hz), this creates high-frequency oscillations in vessel diameter.

**Vasomotion (<0.3 Hz):** Intrinsic rhythmic contractions of vascular smooth muscle cells, independent of cardiac or neural input. This represents active vessel wall dynamics.

**Respiratory modulation (0.5–5 Hz):** Breathing affects vessel diameter through multiple pathways including intrathoracic pressure changes, autonomic modulation, and direct mechanical effects.

### 1.2 Research Questions

This analysis addresses three fundamental questions about cardio-respiratory-vascular interactions:

1. How strongly are slow vessel fluctuations coupled to heart rate and respiration in the frequency domain?
2. Does respiration modulate the amplitude of cardiac pulsatility in vessels?
3. What are the directional relationships between these signals—which physiological process predicts or drives changes in others?

---

## 2. Analysis Pipeline Overview

### 2.1 Data Acquisition and Preprocessing

The analysis begins with synchronized recordings from two data sources:
- **Imaging data:** Vessel diameter measurements with frame timestamps and spatial calibration
- **Physiological signals:** ECG timestamps, respiration waveforms, and detected R-peak indices

**Vessel diameter preprocessing** involves baseline detrending using a 10-minute moving average to remove ultra-slow drift (photobleaching, focus changes), followed by frequency band separation. A lowpass filter (cutoff 0.3 Hz) isolates the slow vasomotion component, while a 4th-order Butterworth bandpass filter (5–15 Hz) extracts cardiac pulsatility. The Hilbert transform then computes the instantaneous envelope of the cardiac component.

**Physiological signal preprocessing** includes temporal alignment between camera and ECG clocks using a stored offset, instantaneous heart rate computation from R-R intervals with artifact rejection (excluding values outside 150–1200 bpm), and respiration bandpass filtering (0.5–5 Hz) with Hilbert phase extraction.

### 2.2 Analytical Methods

**Heartbeat-triggered averaging (ETA):** Extracts epochs of cardiac-filtered vessel diameter locked to each R-peak (−250 ms to +500 ms), then averages across all beats. This reveals the stereotyped vessel response to each heartbeat, including pulse transit time and waveform morphology.

**Respiratory phase-amplitude coupling:** Bins all timepoints by respiratory phase (18 bins from −π to +π), computes mean cardiac envelope in each bin, and calculates a modulation index: MI = (max − min) / mean. Values near zero indicate no respiratory modulation; higher values indicate stronger phase-dependent gating.

**Coherence analysis:** Computes magnitude-squared coherence between slow vessel fluctuations and both HR and respiration using Welch's method (60-second windows). Coherence ranges from 0 (no linear relationship) to 1 (perfect linear coupling) at each frequency.

**Conditional Granger causality:** Tests whether one signal predicts another beyond what the third signal already predicts. Implemented via windowed VAR models (120 s windows, 30 s steps, 10 Hz sampling, 10 lags = 1 s history). For each direction, an F-test compares full vs. reduced regression models.

---

## 3. Results and Interpretation

### 3.1 Coherence: Slow Vessel Motion and Respiration

![Coherence between slow vessel motion and respiration](slow_vessel___respiration.png)

*Figure 1. Coherence between slow vessel diameter fluctuations and filtered respiration signal.*

Coherence values are extremely low throughout the 0–1 Hz frequency range, with a maximum of approximately 0.015 near 0.57 Hz. The sparse peaks around 0.45–0.6 Hz likely correspond to the respiratory frequency band, but the coupling strength is negligible. For reference, meaningful physiological coherence typically exceeds 0.1–0.2, with strong coupling indicated by values above 0.3.

**Interpretation:** Slow vessel diameter changes are not linearly coupled to respiration at matching frequencies. Any respiratory influence on vessel tone operates through nonlinear mechanisms or at timescales not captured by spectral coherence.

### 3.2 Coherence: Slow Vessel Motion and Heart Rate

![Coherence between slow vessel motion and HR](slow_vessel___HR.png)

*Figure 2. Coherence between slow vessel diameter fluctuations and interpolated heart rate.*

Similar to the respiration coherence, values remain at noise floor (0–0.014) across all frequencies. Peaks near 0.2–0.3 Hz and 0.45–0.5 Hz may reflect shared respiratory modulation of both signals (respiratory sinus arrhythmia affects HR while respiration mechanically affects vessels), but coupling strength is negligible.

**Interpretation:** Slow vessel fluctuations and heart rate variability are not strongly linearly coupled. This suggests that vasomotion operates largely independently of beat-to-beat cardiac dynamics in this preparation.

### 3.3 Heartbeat-Triggered Vessel Response

![Heartbeat-triggered vessel response](Heartbeat-triggered_vessel_response.png)

*Figure 3. Event-triggered average of cardiac-filtered vessel diameter locked to R-peaks.*

The waveform shows a regular oscillatory pattern at approximately 15–20 Hz with peak-to-peak amplitude of roughly 0.5 μm. A sharp negative deflection occurs at t=0 (the R-peak trigger), followed by continued oscillation that gradually dampens after 0.3 seconds.

**Methodological concern:** This result likely reflects an artifact of the analysis approach rather than a true physiological response. The oscillation frequency (15–20 Hz) exceeds the mouse heart rate (approximately 10 Hz for 600 bpm), and the pre-trigger oscillations mirror those post-trigger. This pattern suggests the ETA is displaying the autocorrelation structure of the cardiac bandpass signal rather than a genuine pulse-triggered response.

A physiologically meaningful ETA would show: a relatively flat baseline before t=0, a single peak or biphasic deflection at some characteristic latency (typically 10–50 ms representing pulse transit time), and return to baseline. The current result requires methodological revision.

### 3.4 Respiratory Phase Modulation of Cardiac Envelope

![Respiratory phase modulation](Resp_phase_modulation.png)

*Figure 4. Mean cardiac envelope as a function of respiratory phase. MI = modulation index.*

The cardiac envelope shows a clear sinusoidal relationship with respiratory phase. Minimum amplitude (approximately 0.631 μm) occurs near phase = −0.5 radians, while maximum amplitude (approximately 0.648 μm) occurs near phase = ±π. The modulation index of 0.027 indicates a 2.7% variation in cardiac pulsation amplitude across the respiratory cycle.

**Interpretation:** This represents a genuine physiological finding. Respiratory-cardiac interactions are well-established: intrathoracic pressure changes during breathing modulate venous return and stroke volume, affecting pulse amplitude. The smooth sinusoidal pattern and consistent relationship across phase bins support the biological validity of this result.

The 2.7% modulation is modest but meaningful in physiological terms. Assuming phase = 0 corresponds to peak inspiration, the data suggest cardiac pulsation is weakest during mid-expiration and strongest at end-expiration/early inspiration.

### 3.5 Windowed Conditional Granger Causality

![Windowed conditional Granger causality](Windowed_conditional_Granger__-log10_p_.png)

*Figure 5. Time-resolved conditional Granger causality for all six directional relationships. Values shown as −log₁₀(p-value).*

The results show extreme spikes in significance (reaching −log₁₀(p) > 300, corresponding to p-values below 10⁻³⁰⁰) clustered in specific time windows: approximately 1000–1500 s, 2000 s, and 3500–4500 s. After approximately 5500 s, values drop dramatically and remain near the p=0.05 significance threshold for the remainder of the recording.

**Critical concerns:** These results are highly suspicious and should not be interpreted as reflecting genuine directional causality. Several features indicate methodological problems:

1. **Implausible p-values:** P-values of 10⁻³⁰⁰ are not physiologically plausible and suggest the model is fitting artifacts or experiencing numerical issues.

2. **Simultaneous spikes:** All six directional relationships spike simultaneously at the same time points, which indicates a shared confound (such as motion artifact or state change) rather than true directional relationships.

3. **Temporal clustering:** The extreme values cluster in the first half of the recording, suggesting the recording quality or animal state differs substantially between early and late portions of the session.

**Interpretation:** The Granger causality analysis is dominated by artifacts or non-stationarity in the first portion of the recording. The later portion (>6000 s) shows more reasonable values hovering near the significance threshold, which is actually consistent with the weak coherence results—if linear coupling is weak, Granger causality should also be weak or absent.

---

## 4. Summary of Findings

| Analysis | Finding | Confidence / Notes |
|----------|---------|-------------------|
| Vessel–Resp Coherence | Very weak / negligible | High; values at noise floor |
| Vessel–HR Coherence | Very weak / negligible | High; values at noise floor |
| Resp Phase Modulation | Modest effect (MI = 2.7%) | Medium-high; smooth sinusoidal pattern |
| Heartbeat-Triggered Avg | Likely artifact | Low; methodological revision needed |
| Granger Causality | Unreliable | Low; artifacts dominate early recording |

---

## 5. Conclusions and Speculations

### 5.1 Biological Interpretation

The weak coherence results suggest that slow vessel diameter fluctuations (vasomotion) operate largely independently of beat-to-beat cardiac and respiratory dynamics. This is consistent with the understanding that vasomotion represents intrinsic smooth muscle activity rather than a passive response to cardio-respiratory forcing.

The modest but real respiratory modulation of cardiac pulsation amplitude (MI = 2.7%) likely reflects mechanical coupling through intrathoracic pressure: inspiration increases venous return and preload, potentially affecting stroke volume and pulse pressure. The 2.7% effect size is physiologically reasonable for a peripheral vessel where multiple factors modulate the arriving pulse wave.

The apparent absence of strong directional relationships (in the trustworthy portions of the Granger analysis) is consistent with the coherence findings and suggests these physiological systems may interact through higher-order or state-dependent mechanisms not captured by linear time-series methods.

### 5.2 Alternative Hypotheses

Several factors could explain the weak observed coupling:

- **Local autoregulation:** The vessel being imaged may be in a vascular bed with strong local autoregulation that buffers against systemic cardio-respiratory influences.

- **Physiological state:** The animal's state during recording (anesthesia depth, arousal level) may affect coupling strength.

- **Frequency band mismatch:** The frequency bands chosen for analysis may not optimally capture the timescales of interaction.

- **Nonlinear relationships:** The relationships may be fundamentally nonlinear, requiring different analytical approaches such as information-theoretic measures or phase-phase coupling.

---

## 6. Recommendations for Next Steps

### 6.1 Immediate Methodological Fixes

**Heartbeat-triggered average revision:** Use the raw detrended signal rather than the cardiac-filtered signal for ETA computation. This avoids the autocorrelation artifact. Alternatively, if using the filtered signal, ensure the filter bandwidth is narrow enough that beats are isolated, or use a different triggering approach.

**Granger causality troubleshooting:** Visually inspect raw signals during the extreme-spike time periods (1000–1500 s, 3500–4500 s) to identify artifacts or state changes. Consider excluding the first 5000 s of recording or segmenting analysis by behavioral state. Add stationarity tests (ADF test) per window and exclude non-stationary segments.

**Statistical rigor:** Add coherence significance thresholds based on the number of segments. Implement permutation testing for the phase-amplitude modulation to confirm the MI = 2.7% exceeds chance. Apply FDR correction for multiple comparisons in windowed Granger analysis.

### 6.2 Extended Analyses

**Nonlinear coupling measures:** Consider mutual information or transfer entropy as alternatives to Granger causality. These capture nonlinear relationships that coherence and linear VAR models miss.

**State-dependent analysis:** If behavioral state annotations are available, analyze coupling separately during different states (rest, movement, different arousal levels). Coupling may be strong in some states and absent in others.

**Cross-frequency coupling:** The current phase-amplitude analysis examines respiration phase vs. cardiac envelope amplitude. Consider also testing cardiac phase vs. slow vessel amplitude, or phase-phase coupling between respiration and cardiac rhythms.

**Adaptive frequency bands:** Compute the actual heart rate and respiratory rate from the data and adjust filter bands accordingly, rather than using fixed values that may not match the specific animal.

### 6.3 Data Quality Verification

Before proceeding with refined analyses, verify:

- The temporal alignment between camera and ECG is accurate throughout the recording (not just at the start)
- R-peak detection quality is consistent across the recording
- The vessel diameter measurement is not contaminated by motion or focus artifacts during the periods showing Granger anomalies

---

## Appendix: Technical Parameters

| Parameter | Value |
|-----------|-------|
| Baseline detrend window | 10 minutes (moving mean) |
| Slow component cutoff | 0.30 Hz (lowpass) |
| Cardiac band | 5–15 Hz (4th order Butterworth) |
| Respiration band | 0.5–5.0 Hz (bandpass) |
| HR artifact rejection | 150–1200 bpm accepted |
| ETA window | −250 ms to +500 ms from R-peak |
| Phase modulation bins | 18 bins from −π to +π |
| Coherence window | 60 seconds (Hamming) |
| Granger sampling rate | 10 Hz (downsampled) |
| Granger window / step | 120 s / 30 s |
| VAR lag order | 10 samples (1 second history) |

---

*Report generated from analysis pipeline: sketch.m*
