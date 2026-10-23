# Pressure-peak tracking method

This release tracks fitted peak estimates, not raw-signal maxima or inferred
chemical identities. Each spectrum is fitted independently. Only PASS fits with
finite pressure metadata update the tracker. All other fitted estimates remain
in the output with an exclusion event.

Records are processed by acquisition sequence within explicitly independent
series. Repeated pressures remain separate acquisitions. The optional linear
predictor uses the last accepted center and the most recent earlier accepted
point at a distinct pressure. Pressures equal within a 1e-12 absolute/relative
floating-point tolerance do not define a slope. `nearest` uses the last center.
The predictor is a local modeling assumption, not a general physical law.

For active track i and fitted observation j, the real-edge cost is the absolute
center-prediction residual plus `width_weight` times the absolute FWHM change.
Edges beyond `max_distance_cm1` are forbidden. A square augmented cost matrix adds
one gap option per track and one birth option per observation, each costing half
`max_distance_cm1`; the dummy-to-dummy block costs zero. This permits unmatched
observations and disappearing peaks instead of forcing every point into a track.

[SciPy linear_sum_assignment](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html)
solves each minimum-cost assignment. To assess ambiguity, every feasible track
identity and the birth option are forced in turn for each observation; the
remaining global assignment is re-solved. An identity remains competitive if its
optimal global cost is within `ambiguity_margin_cm1` of the original optimum.
Multiple competitive identities yield AMBIGUOUS, with no asserted track ID and
no predictor update. These margins are heuristic cost differences, not calibrated
probabilities or confidence intervals. A small exhaustive permutation oracle
checks the assignment/candidate logic independently in tests.

An unobserved acquisition creates a GAP with no measured center. A trajectory
expires after more than `max_gap` intervening acquisitions within its own series.
Interleaved acquisitions from a different series do not consume that budget.
After expiry, a returning peak starts a new ID. IDs are stable during the run but
do not represent universal peak names across independently initialized runs.

The synthetic crossing demonstration samples both sides of a frequency-order
reversal, with explicit fitting priors. Unit fixtures also exercise coincident
predictions, where identities remain ambiguous. A unique result under the chosen
predictor does not prove identity across a discontinuity, unresolved overlap,
mode mixing or pressure hysteresis. Series boundaries, model choice and physical
interpretation remain explicit scientific decisions.

Tracking figures connect adopted fit estimates in acquisition order and break
lines at explicit gaps. Unknown local errors remain missing; excluded estimates
are displayed separately. Lines guide the eye and do not add pressure samples.
No extrapolated centers are reported as measurements.
