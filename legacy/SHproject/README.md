# Historical SHproject implementation

The processing modules and command-line scripts are preserved unchanged from
the owner-provided SHproject source. `source-manifest.json` records SHA-256 values.
Experiment-specific pressure tables and photoactivity profiles are replaced by
public templates; measurement files, archives and historical outputs are excluded.

This is an **experimental legacy interface**, outside the supported wheel API.
Install the repository dependencies and `tqdm` to use it. From this directory,
`python run_pipeline.py --list-profiles` lists the generic template. Populate
explicit spectrum pressures locally, or use the supported `raman-analyze` entry
point, which does not require pressure data. The photoactivity template is empty
until locally configured. Do not infer validated general behavior from this archive.

Known issues and restrictions are documented in `../../docs/legacy-audit.md`.
