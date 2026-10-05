# Swift UVOT Light-Curve Pipeline — Documentation

A step-by-step pipeline for multi-epoch light curves of point sources from Swift UVOT data.
It is the companion to the Swift XRT pipeline
([documentation](https://etmeyer.github.io/SwiftXRT_Complete/)).

> **Status (October 2026):** Steps 1–3 are done; the other steps are being written. Each step
> gets its own page, with the mechanism, inputs and outputs, gotchas and space for notes, as in
> the XRT documentation.

## Planned workflow

| Step | What it does |
| ---- | ------------ |
| [1. Setup](01-setup.md) | Install the pipeline and the UVOT CALDB; check the terminal with `swift_uvot_doctor.py` (one terminal, HEASoft, no CIAO) |
| [2. Download](02-download.md) | List the observations of a source with their UVOT exposure per filter; download sky images, exposure maps, raw images, event files and housekeeping; check them against the catalog (`swift_uvot_download.py`) |
| [3. Inventory](03-inventory.md) | Classify every exposure: filter, frame time, image or event mode, aspect correction, snapshot, whether the source is in the exposed field; one status each (`swift_uvot_inventory.py`) |
| 4. Positions and regions | Measure the source position in each exposure; write source and background regions; diagnostics such as trailed images and sensitivity patches |
| 5. Photometry | Run `uvotsource` on every exposure, keeping every number and flag |
| 6. Master table | Apply the quality rules and your overrides; one row per exposure with the reason it is used or excluded |
| 7. Light curve | Combine the good exposures into one point per observation and filter; detections, upper limits and saturated lower limits; plots |

A per-exposure light curve (for short-timescale variability) is a planned extension.

## What the pipeline checks for

- **Saturation.** Bright sources can exceed the coincidence-loss limit (0.98 raw counts per
  frame) in full-frame exposures; those values are kept only as lower limits.
- **Pointing.** Exposures without aspect correction, trailed exposures taken while the
  spacecraft was still settling, and the spacecraft jitter of August 2023 – April 2024.
- **Detector sensitivity patches** flagged by the UVOT calibration.
- **Coverage.** Whether the source is inside the image at all.
