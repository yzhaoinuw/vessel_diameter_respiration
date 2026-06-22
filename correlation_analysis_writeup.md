# Correlation Analysis of Heart Rate and Vessel Diameter Dynamics

## Quick Summary

### Aim

The aim of the analysis is to determine whether changes in heart rate were associated with changes in vessel diameter dynamics in the dex40 mice. We examined mean vessel caliber and three features of diameter fluctuation: slow vasomotion amplitude, cardiac-frequency pulsation amplitude, and slow vasomotion power. Since the ECG recordings can be noisy for estimating the heart rate, We also compared three representations of heart rate to evaluate whether the ECG-derived measurements were internally consistent and whether they agreed with the cardiac frequency visible in the vessel signal.

### Analysis

Heart rate and vessel features were summarized across repeated time windows within each mouse. We used quality-filtered one-minute ECG heart rate values paired with overlapping ten-minute vessel windows, also advanced in one-minute steps. Within each recording, heart rate was correlated separately with each vessel feature during the early period (< 60 minutes), late period (>=60 minutes), and full recording. The resulting correlation coefficients were then summarized across mice, with the mouse retained as the unit of replication.

The heart rate was computed from a ECG processing pipeline in which a machine learning model was trained to identify R-peaks with a confidence value between 0 and 1. The R-peaks were then used to compute the heart rate. We experimented two different methods for post-processing the R-peaks and their assocaited confidence value. They differ in time-bin duration and in how low-quality periods or recordings were identified and handled. Since these two methods only differ in the prost-processing and thus not independent from each other, we also had a third estimate that was obtained from the cardiac-frequency peak in the vessel diameter spectrum and provided a cross-modal comparison.

### Results

The overlapping-window analyses produced a recurring exploratory negative association between heart rate and slow vasomotion amplitude: mice tended to show smaller slow vessel oscillations when heart rate was higher. This association reached nominal significance in the early and late periods using Pearson correlation in the broader quality-qualified cohort and in the early period using Spearman correlation in the more stringently curated cohort. However, it was not consistently significant across all windows, correlation measures, or cohort definitions. Slow vasomotion power, cardiac-frequency pulsation amplitude, and mean vessel diameter did not show a stable relationship with heart rate.

The three heart rate representations generally followed similar temporal patterns. In particular, the two ECG summaries agreed strongly despite their different binning and quality-handling rules, and both showed broad agreement with the vessel-derived cardiac frequency. These comparisons support the technical reliability of the heart rate measurement, but they do not by themselves demonstrate biological coupling between heart rate and vessel dynamics.

## Heart Rate Estimation and Quality Control

Three heart rate estimates were compared for quality control: a quality-filtered one-minute ECG estimate, an earlier ten-minute ECG estimate, and a cardiac-frequency estimate derived from vessel diameter. The purpose was to corroborate the overall heart rate trajectory, evaluate the consequences of different ECG quality-control choices, and identify recordings in which heart rate was uncertain. Only the vessel-derived estimate came from a separate measurement modality.

### Quality-Filtered One-Minute ECG Heart Rate

The first estimate was calculated directly from the ECG. Heart rate was calculated from the interval between consecutive detected R peaks and expressed in beats per minute. Values outside the physiologically plausible range of 250 to 750 beats per minute were treated as invalid.

ECG quality was evaluated separately in each one-minute interval. A minute was rejected when at least half of its candidate R peaks had low confidence. For accepted minutes, the median valid heart rate was retained. Rejected minutes remained missing and were not filled by interpolation or extrapolation. Recordings with extensive missing coverage could then be excluded from the correlation analysis rather than having their poor-quality periods reconstructed.

This more conservative representation was used in the overlapping-window analysis because it retained one-minute temporal sampling while preventing low-confidence intervals from contributing artificial heart rate values.

### Ten-Minute ECG Heart Rate

The second estimate was the ECG heart rate summary produced by the earlier analysis workflow in ten-minute bins. It drew on the same R-peak detection and confidence assignments as the one-minute estimate and therefore was not an independent measurement. The two ECG estimates differed in bin duration and in the rules used to identify and handle low-quality periods or recordings.

The ten-minute ECG summary was retained as a corroborating calculation. Agreement between the one-minute and ten-minute ECG summaries tested whether the estimated heart rate trajectory was robust to those differences in binning and quality handling.

### Vessel-Derived Cardiac Frequency

The third comparison was the dominant cardiac-frequency peak in the vessel diameter power spectrum, converted from frequency to beats per minute. Unlike the two ECG summaries, this estimate was derived from the imaging signal. It therefore provided a cross-modal check that the cardiac rhythm identified in the ECG was also visible in vessel pulsation.

The vessel-derived estimate was used for heart rate validation only. It was not used as the heart rate measurement in the main heart rate-vessel correlation analysis because it is itself calculated from the vessel signal and would not constitute an independent predictor of vessel behavior.

Agreement among the three representations was examined within each recording. The comparison established whether they followed the same broad temporal pattern and helped identify recordings in which the estimates diverged enough to raise concern about signal quality. Agreement was interpreted as evidence about heart rate measurement reliability, not as evidence that heart rate controlled or was coupled to vessel diameter dynamics.

### ECG Heart Rate Used in the Five-Minute Correlation Analysis

The initial correlation analysis used the same beat-level ECG information summarized within five-minute windows. Candidate R peaks were retained when they passed the established confidence threshold and belonged to the set accepted during ECG quality review. The median heart rate was calculated for each window because it is less sensitive than the mean to isolated R-peak errors or unusually short or long beat intervals.

A five-minute window was considered poor quality when it contained too few valid beats or produced an implausible heart rate. When sufficient acceptable windows were available in a recording, poor-quality windows were replaced by linear interpolation or extrapolation from the acceptable windows. This analysis-specific five-minute summary was not treated as a fourth validation estimate; it was the form of ECG heart rate matched to the five-minute vessel features in the initial correlation analysis.

## Vessel Diameter Processing

The vessel diameter trace was aligned to the ECG time base using the recorded timing offset between the imaging and physiological acquisition systems. Diameter measurements were converted from pixels to micrometers before feature calculation.

Slow drift in baseline vessel diameter was removed by subtracting a ten-minute moving average from the original diameter trace. This step preserved fluctuations occurring on shorter time scales while reducing the influence of gradual changes in imaging position, vessel baseline, or recording conditions.

Three frequency-based vessel features were emphasized:

1. **Slow vasomotion amplitude.** The detrended diameter signal was filtered to retain slow fluctuations up to 0.3 Hz. The amplitude of these fluctuations within each analysis window was summarized by their interquartile range. This provided a robust estimate of the size of slow vessel oscillations.

2. **Slow vasomotion power.** The power spectrum of the detrended diameter signal was estimated within each analysis window. Spectral power was integrated from 0.01 to 0.3 Hz, providing a complementary measure of the strength of slow vasomotor activity.

3. **Cardiac-frequency vessel pulsation amplitude.** The detrended diameter signal was filtered within the cardiac range used for the dexmedetomidine recordings, 3.5 to 15 Hz. The interquartile range of the filtered signal was used to quantify the amplitude of vessel pulsation associated with the cardiac cycle.

The initial five-minute analysis also included mean vessel diameter. Mean diameter described the vessel's baseline caliber within a window, whereas the frequency-based features described the magnitude or power of ongoing diameter fluctuations.

## Time-Window Definitions

We used two related time-window strategies to test whether the result depended on temporal resolution.

The initial analysis summarized both heart rate and vessel features in five-minute windows advanced in one-minute steps. These windows began 14 minutes after the start of the recording. Median ECG heart rate was paired with slow vasomotion amplitude, slow vasomotion power, cardiac-frequency pulsation amplitude, and mean vessel diameter calculated over the same five-minute interval.

The later analysis increased the duration of the vessel window to ten minutes while retaining a one-minute step between consecutive windows. The longer window provided more data for estimating low-frequency vessel activity, particularly spectral power in the slow vasomotion band. The resulting windows overlapped substantially, producing a smoothly sampled description of how vessel dynamics changed over time.

In the implemented overlap analysis, each ten-minute vessel window was paired with the quality-filtered one-minute ECG heart rate value at the start of that window. Only windows with finite values for both measurements were included in a correlation. Because this pairing uses the heart rate at the window start rather than the average heart rate across the entire ten-minute vessel window, it should be described explicitly when the analysis is reported.

For both windowing approaches, correlations were evaluated over three recording periods:

- **Early:** times before 60 minutes.
- **Late:** times at or after 60 minutes.
- **Overall:** all available matched time points.

The early and late divisions were used to test whether heart rate-vessel relationships changed over the course of the recording rather than assuming that one relationship remained constant throughout the experiment.

## Within-Mouse Correlation Analysis

For each mouse, heart rate was correlated separately with each vessel feature across the matched time windows. Missing or rejected observations were removed pairwise, so a correlation used only time windows in which both heart rate and the vessel feature were available.

Spearman rank correlation was used as the main correlation measure. This approach tests whether heart rate and a vessel feature changed together in a consistent monotonic direction without requiring their relationship to be linear or the measurements to be normally distributed. It is also less sensitive than Pearson correlation to extreme values that may remain after physiological signal processing.

Pearson correlations were also calculated in the overlapping-window analysis as a complementary assessment of linear association. Comparisons between Pearson and Spearman results helped determine whether any apparent relationship was broadly consistent or depended strongly on the assumed form of the association.

The analyses produced one correlation coefficient for each mouse, vessel feature, and recording period. A positive coefficient indicated that the vessel feature tended to increase as heart rate increased. A negative coefficient indicated that the vessel feature tended to decrease as heart rate increased. The magnitude of the coefficient described the strength of the within-recording association.

## Recording-Level Quality Review

The quality-filtered overlap analysis required adequate ECG coverage. Recordings with extensive missing heart rate values could not provide a stable estimate of association and were excluded when more than half of the aligned heart rate observations were unavailable.

The diagnostic comparison of heart rate representations also identified recordings with substantial disagreement among the estimates. Six recordings were excluded from the final curated Spearman overlap analysis on this basis, leaving five recordings with the clearest agreement among heart rate measurements. Because this exclusion was informed by diagnostic review, the curated analysis should be presented together with the broader, less restrictive analyses as a sensitivity comparison rather than as the only view of the data.

## Group-Level Statistical Analysis

The mouse, rather than the individual time window, was used as the unit of inference. Time windows within a recording are repeated and strongly overlapping observations from the same animal; pooling them across animals would therefore overstate the effective sample size.

For each vessel feature and recording period, the distribution of mouse-level correlation coefficients was evaluated in two complementary ways. A Wilcoxon signed-rank test assessed whether the median correlation across mice differed from zero without assuming a normally distributed set of coefficients. A one-sample test was also performed after applying the Fisher transformation to the correlation coefficients. The transformed values were tested against zero, and the group mean and confidence interval were converted back to the correlation scale for presentation.

Group summary plots show one point per mouse, a group mean correlation, a 95% confidence interval, and a reference line at zero. This format displays both the average direction of the association and the degree of consistency or heterogeneity across animals.

The current analyses tested several vessel features across early, late, and overall periods. The resulting probability values were not adjusted for the number of comparisons. The correlation findings should therefore be described as exploratory unless a multiple-comparison strategy is selected before the final manuscript analysis.

## Planned Plots

### Main-Text Plots

1. **Analysis workflow schematic.** A compact diagram showing ECG R-peak detection, quality-filtered heart rate estimation, vessel detrending, separation into slow vasomotion and cardiac-frequency features, temporal window matching, within-mouse correlation, and group-level summary.

2. **Heart rate quality-control comparison.** Representative time-course overlays of the quality-filtered one-minute ECG heart rate, the related ten-minute ECG summary, and the vessel-derived cardiac-frequency estimate. A small group-level agreement panel could accompany the example recording to show that the calculations generally tracked the same physiological quantity while making clear that the two ECG summaries share the same underlying peak-quality pipeline.

3. **Representative within-mouse correlation plots.** Scatter plots from one high-quality recording showing heart rate against slow vasomotion amplitude, slow vasomotion power, and cardiac-frequency pulsation amplitude. Early and late observations would be distinguished by color, with a monotonic trend line added for visualization.

4. **Primary group-level correlation summary.** A dot-and-interval plot for the curated overlapping-window Spearman analysis. Each point would represent one mouse, organized by vessel feature and by early, late, and overall recording period. The plot would include the group estimate, 95% confidence interval, and zero-correlation reference line.

### Supplementary Plots

5. **Analysis sensitivity comparison.** Side-by-side group summaries comparing the five-minute and ten-minute window analyses, Pearson and Spearman correlations, and the broader versus curated recording sets. This would show which conclusions are stable across reasonable analysis choices.

6. **Recording-level correlation matrix.** A heat map showing the direction and magnitude of the correlation for every mouse, vessel feature, and recording period. This would make between-mouse heterogeneity visible without pooling the underlying time points.

7. **Data coverage and quality-control summary.** A plot showing the proportion of valid heart rate observations in each recording and indicating which recordings were retained or excluded from the curated analysis. This would make the quality-control process transparent.

## Points to Resolve Before Final Manuscript Use

The main remaining methodological decision is how heart rate should be paired with the ten-minute vessel windows. The current implementation uses the quality-filtered heart rate at the start of each vessel window. Averaging or taking the median heart rate over the same ten-minute interval would provide closer temporal matching and may be preferable for a final manuscript analysis.

The role of the curated five-recording cohort should also be stated clearly. A defensible presentation would treat the broader analysis as the primary description of the available recordings and the curated high-agreement cohort as a quality-focused sensitivity analysis, unless the exclusion criteria can be formalized independently of the observed correlation results.

Finally, the project team should decide whether the manuscript will designate one vessel feature and one recording period as the primary comparison or apply a correction across all feature-by-period tests. Until that decision is made, the group-level probability values should be interpreted as exploratory.
