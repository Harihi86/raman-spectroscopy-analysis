# Methodology and evidence boundaries

Raw input, preprocessing output, fitted curves and interpretation are separate
artifacts. `processed.csv` preserves loaded raw intensities, estimated baseline
and the corrected difference. Input hashes, explicit missing-value decisions,
axis sorting and analysis parameters are recorded in `report.json`. No resampling,
normalization, smoothing or scientific classification happens implicitly.

AsLS minimizes a weighted squared data term plus an index-space second-difference
penalty, updating weights asymmetrically for positive peaks. The smoothing
parameter depends on sampling density and noise; it is not universally transferable
between instruments or spectra. A fixed number of iterations is recorded, not
reported as proof of convergence of every baseline choice.

The fit is a sum of positive Gaussians and an unconstrained linear residual
background. The intercept uses a centered axis for numerical conditioning. Inputs
need more points than fit parameters. The default loss is ordinary least squares;
no robust-loss covariance or noise-distribution guarantee is claimed.

For a Gaussian of amplitude A, center mu and width sigma:

- FWHM is `2 sqrt(2 ln 2) sigma`.
- Area is `sqrt(2 pi) A sigma`.
- Area variance uses the full A/sigma covariance block and a local delta-method
  gradient, rather than adding independent percentage errors.

When measurement standard deviations are supplied, weighted residuals and absolute
local covariance are used. Otherwise a residual-scaled independent equal-variance
noise model is assumed. The covariance is a linear approximation conditional on
the selected baseline, bounds, number of peaks and lineshape. Preprocessing creates
correlations not represented by this local matrix. Active constraints and nearly
unidentifiable peak combinations need further scientific analysis, even when the
optimizer returns a solution.

The program reports optimizer convergence, rank, condition number, active-bound
diagnostics, residual degrees of freedom and residual errors. PASS means these
configured numerical diagnostics passed. It is not evidence of correct band
assignment, calibrated uncertainty, phase transition, chemical identification or
experimental reproducibility across instruments.

Synthetic tests exercise known peaks across multiple counts and axis windows.
Local measured pressure spectra have now been exercised with full outcome and sensitivity reporting; see [measured-data evidence](real-validation.zh-CN.md). Physical count accuracy still requires independent reference labels.


Residual adequacy is reported separately from optimizer convergence. A lag-1
correlation magnitude above 0.35 (configurable with
`residual_correlation_limit`, or explicitly disabled with null) and a reduced
chi-square above 4 when measurement errors are supplied are heuristic reasons
for INCONCLUSIVE. These thresholds are diagnostic rules, not calibrated
statistical hypothesis tests. Near-zero numerical residuals are excluded from
the correlation rule to avoid interpreting floating-point noise.

During validation, an initial synthetic example with AsLS lambda=1e7 showed
correlated residuals and width bias. The released example explicitly uses
lambda=1e9 for its uniform synthetic grid. This parameter is not a universal
recommendation. Its saved truth, residual plots and report expose remaining
baseline/model bias instead of treating convergence as physical validation.

In manual positive-peak mode, flat or affine spectra within a relative intensity tolerance of 1e-12 do not establish a Gaussian peak: peak covariance is withheld and the result stays INCONCLUSIVE. Automatic mode explicitly includes a zero-peak background and can select it without asserting that a physical mode is absent. See [automatic model comparison](automatic.md).

An explicit AsLS `resample:true` option estimates only the baseline on a same-length uniform grid for irregular inputs; fitted observations stay on their original axes. The additional baseline uncertainty is not included in the local covariance. When supplied, fitted-sample standard errors are exported as `measurement_sigma` in `fit.csv`.
