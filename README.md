# Commbio_Hb_TAMRA_reporter_analysis
Matlab analysis associated with "Secreted hemophores drive early heme scavenging and interspecies competition for hemoglobin in pathogenic bacteria"

# MATLAB analysis code for [Manuscript Title]

This repository contains MATLAB scripts associated with:

[Jordan Ford, Jess Soule, Kylie Cheng, Brendan Mahoney, Clarissa Pham, John Meinert, Kat Ellis-Guardiola, and Robert T. Clubb]. "[Secreted hemophores drive early heme scavenging and interspecies competition for hemoglobin in pathogenic bacteria]." [Communications Biology].

## Contents

- rateanalysis.mlx — [accepts 3 replicate .csv files with paired OD and fluorescence data/baselines, performs analysis of timing parameters and rates with dipslays]
- analyze_one_rep.m — [module which performs most of the math, required]

## Usage

For each replicate, the OD baseline is fit using a cubic smoothing spline with smoothing parameter 0.2 and subtracted from the corresponding experimental OD trace. The fluorescence baseline is independently fit using a cubic smoothing spline with smoothing parameter 0.7 and subtracted from the experimental fluorescence trace. Baseline-corrected traces are subsequently smoothed using a three-point moving mean.

The temporal relationship between fluorescence and OD is evaluated using normalized rising-phase traces. Metrics include continuous-shift normalized cross-correlation, derivative cross-correlation, and the difference between the 50% rise times of the two signals.

Fluorescence normalized to baseline-corrected OD is analyzed over an interval in the time units of the input data. Rate estimates are obtained using a second-order Savitzky–Golay analysis and cubic smoothing-spline analysis. The Savitzky–Golay analysis uses a nominal five-point frame, subject to adjustment when required by the number of available points.

The scripts also retain processed baseline-corrected and normalized time series for visualization and calculation of replicate-average trajectories.

## Output files

Running rateanalysis produces replicate-level metrics, summary statistics, and numerical data corresponding to plotted analyses. Output CSV files include triplicate_metrics.csv, triplicate_summary.csv, FLC_OD_mean_plotdata.csv, and FLC_OD_derivative_plotdata.csv. When a Savitzky–Golay peak rate is identified, FLC_OD_peak_SG.csv is also generated.

## MATLAB requirements

The code requires MATLAB and uses functions from the Curve Fitting Toolbox, including csaps, fnval, and fnder, and the Signal Processing Toolbox, including sgolay, sgolayfilt, and xcorr.

The analysis was performed using MATLAB [R2023b].

## Usage

Place rateanalysis.mlx (or its plain-text .m equivalent), analyze_one_rep.m, and the three input CSV files in the same working directory.

Edit the files variable near the beginning of rateanalysis so that it contains the filenames of the three replicate datasets. Run rateanalysis in MATLAB. The script will analyze each replicate, calculate summary metrics, generate the associated plots, and write the numerical results to CSV files in the working directory.

## Citation

If using this code, please cite the associated manuscript:

[manuscript citation]

## License

Released under the MIT License.
