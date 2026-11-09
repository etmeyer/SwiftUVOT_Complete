# Swift UVOT Light-Curve Pipeline — Documentation

A step-by-step pipeline for multi-epoch light curves of point sources from Swift UVOT data.
It is the companion to the Swift XRT pipeline
([documentation](https://etmeyer.github.io/SwiftXRT_Complete/)).

> **Status (November 2026):** Steps 1–8 are done, tested on 3C 273 and on the faint, distant PKS 0537−286. Each step
> gets its own page, with the mechanism, inputs and outputs, gotchas and space for notes, as in
> the XRT documentation.

## Planned workflow

| Step | What it does |
| ---- | ------------ |
| [1. Setup](01-setup.md) | Install the pipeline and the UVOT CALDB; check the terminal with `swift_uvot_doctor.py` (one terminal, HEASoft, no CIAO) |
| [2. Download](02-download.md) | List the observations of a source with their UVOT exposure per filter; download sky images, exposure maps, raw images, event files and housekeeping; check them against the catalog (`swift_uvot_download.py`) |
| [3. Inventory](03-inventory.md) | Classify every exposure: filter, frame time, image or event mode, aspect correction, snapshot, whether the source is in the exposed field; one status each (`swift_uvot_inventory.py`) |
| [4. Positions and regions](04-positions.md) | Measure the source position, shape and concentration in each exposure; flag trailed and smeared images and sensitivity patches; Gaia field stars for each exposure's pointing and PSF; write source and background regions; contact sheets (`swift_uvot_positions.py`, `swift_uvot_viewer.py`) |
| [5. Photometry](05-photometry.md) | Run `uvotsource` on every exposure, keeping every number and flag; plot every exposure's rate against time, marked by what Step 6 will judge (`swift_uvot_photometry.py`) |
| [6. Master table](06-master-table.md) | Apply the quality rules, calibrated on 3C 273, and your overrides; one row per exposure with the rules that excluded it (`make_uvot_master_table.py`) |
| [7. Light curve](07-lightcurve.md) | Combine the included exposures into one point per observation and filter (weights from expected counts, background systematic): Poisson detection, KBN upper limits, sensitivity and saturated lower limits; magnitudes and flux densities; plot (`swift_uvot_lightcurve.py`) |
| [8. Flux densities](08-flux.md) | Flux densities corrected for Galactic dust (E(B−V) from IRSA, extinction folded through each filter) and, for a distant source, the intergalactic medium; a spectrum of your choice; spectral energy distributions (`swift_uvot_flux.py`) |

A per-exposure light curve (for short-timescale variability) is a planned extension.

## What the pipeline checks for

- **Saturation.** Bright sources can exceed the coincidence-loss limit (0.98 raw counts per
  frame) in full-frame exposures; those values are kept only as lower limits.
- **Pointing.** Exposures without aspect correction, trailed exposures taken while the
  spacecraft was still settling, and the spacecraft jitter of August 2023 – April 2024.
- **Detector sensitivity patches** flagged by the UVOT calibration.
- **Coverage.** Whether the source is inside the image at all.
