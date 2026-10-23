# SHproject source audit and migration decisions

The owner identified SHproject as the source for this public project. Static
inspection found the following concrete restrictions; they are not inferred from
its name or README alone. Line numbers refer to the preserved original modules,
or to the private profile file where noted.

| Source location | Observed behavior | Public treatment |
|---|---|---|
| Original `profiles.py`, pressure lists and named EXP1/EXP2 profiles | Hard-coded pressure sequences, instrument/sample windows and comparison groups | Experimental values withheld; generic public templates supplied |
| Original `profiles.py`, heatmap conditional-fill rules | Selected pressure/frequency regions can be filled with zero, affecting visualization | No fill rule in supported API; missing/omitted points must not be represented as measured zero |
| `core_processing.py:159`, `trim_xy` | Truncates unequal arrays, drops invalid pairs, sorts and averages duplicate axis values | New import validates rows; dropping is opt-in and counted; duplicates rejected |
| `core_processing.py:489`, `detect_read_mode` | Even-column numeric tables can be classified as alternating X/Y pairs even if they are shared-X tables | New API requires explicit columns and never guesses layout |
| `core_processing.py:510`, `read_excel_spectra` | Excel-only import with header inference by numeric coercion | Explicit header or headerless CSV/TSV/TXT/XLSX import; units declared |
| `core_processing.py:572`, `als_baseline` | Index-space penalty; equality can assign zero weight on a flat signal | New AsLS requires uniform spacing, validates parameters and keeps positive equality weights |
| `core_processing.py:865` and `1045` | Single/double Gaussian routines, positive amplitude and local linear backgrounds | Manual 1..8 Gaussian/Lorentzian/Voigt models or automatic family/count comparison; bounds and ambiguity remain explicit |
| `core_processing.py:1276`, `fallback_peak_max` | Returns `success=True` for a maximum-based fallback while width/covariance remain NaN | New fit keeps FAIL/INCONCLUSIVE; no maximum fallback |
| `core_processing.py:2416`, `process_profile` | A pressure must exist for every spectrum, even for ordinary per-spectrum analysis | New API does not require pressure metadata |
| `core_processing.py:3120`, `_select_stable_segment` | Chooses favorable contiguous pressure segments for comparison | No automatic scientific trend/phase interpretation in public API |
| `photoactivity_processing.py:41`, `109` | Time-series Excel layout and configured excitation, exposure, expected frames, early/late windows | Source retained as experimental; no general time-series validation claimed |
| Legacy integration functions | Use `np.trapz`, removed from current NumPy | Full legacy workflows need a compatible NumPy environment or a separately tested migration |

The processing modules and original command-line scripts are preserved byte for
byte with a SHA-256 manifest. The original README's experiment-specific commands,
private profiles, raw spreadsheets, previous figures, archives and run logs are
not published. The new pipeline is an additive API; it does not alter the original
folder or assert that every historical analysis can be reproduced publicly.

The original dependency file is retained as source provenance, not a reproducible
lock. To run the unmodified full legacy workflow, use an isolated compatible
NumPy 1.x environment and validate it with authorized data. Only the legacy
nonflat AsLS differential check and listing interface are tested in this release.

Version 0.3 also exercises the new API on locally authorized measured pressure spectra. This does not run or certify the full historical workflow and does not use its fit-success flags or peak counts as reference truth. See [local validation scope](real-validation.zh-CN.md).
