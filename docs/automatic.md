# Automatic Gaussian proposals and component-count comparison

[中文操作说明](automatic.zh-CN.md) · [Configuration](configuration.md)

Use `fit.auto` instead of `fit.peaks`. The implementation compares a zero-component
linear background and 1–N positive Gaussian models on identical original samples.
Detection uses [SciPy find_peaks](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html)
and prominence-relative widths on a smoothed proposal signal. Widths are mapped
from fractional sample indices onto actual x coordinates. The optimizer always
uses unsmoothed observations. Local starts use detected proposals, unexplained
positive residuals and split broad components; this is a bounded search, not a
proof of a global optimum or a count of physical modes.

For N peaks, p = 3N + 2 fitted mean parameters (amplitude, center and sigma per
component, plus a residual-background intercept and slope). Let n be the number
of original fitted samples and RSS the unweighted residual sum of squares.

Without supplied measurement errors, estimate a common Gaussian variance, use
k = p + 1 and drop constants common to all candidates:

```
D = n log(max(RSS/n, smallest_positive_normal_float))
BIC = D + k log(n)
AICc = D + 2k + 2k(k+1)/(n-k-1)
```

With known independent standard errors, D = sum((residual_i / sigma_i)^2) and
k = p. The omitted log-error normalization is common because the same data and
errors enter every candidate. The same displayed AICc correction is used as a
ranking heuristic. It is not an exact finite-sample guarantee for nonlinear
Gaussian-component mixtures or known heteroscedastic errors. AICc is unavailable
when n <= k+1. Within each peak count, converged starts are compared using this
likelihood objective, including supplied weights.

References for the criteria: [Schwarz (1978)](https://doi.org/10.1214/aos/1176344136),
[Hurvich and Tsai (1989)](https://doi.org/10.1093/biomet/76.2.297).
Mixture boundaries, overlapping components, estimated baselines and correlated
residuals violate the simple regular-model interpretation. Score differences are
not calibrated peak-count probabilities or Bayes factors in this implementation.

All counts within `delta_threshold` of the minimum remain competitive. A best
candidate at `max_peaks`, a count without any converged start, a truncated or
weak component, unresolved separation or a selected-fit diagnostic prevents a
PASS. Numerically converged but diagnostically questionable fits remain visible
in the score table and cannot silently promote the result to trusted tracking.
The zero-peak candidate prevents forced discovery of peaks in featureless data.
A zero-component decision says nothing by itself about the existence of physical
Raman modes below resolution or under the background.

`search_complete` means every requested count had a converged attempt, with an
explicit affine-input exception, not that all parameter combinations were searched.
All local starts and diagnostics are saved in `model-selection.json`;
`model-selection.csv` provides a compact score table. Standard errors remain
conditional on the chosen model/baseline; count-selection uncertainty is not
included in those local errors.

Defaults: max count 5 (allowed 1–8), BIC, score margin 6, proposal SNR 4,
1.5-sample proposal smoothing, 3 starts (allowed 1–12). Physical sigma bounds
and separation may be declared; otherwise they derive from the actual sampling
and fit span. No material name, pressure table or experiment profile enters this
module. An analytic Gaussian Jacobian is tested against central finite differences
with and without measurement weighting.
