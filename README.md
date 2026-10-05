# Swift UVOT Light-Curve Pipeline

A Python pipeline, under development, for multi-epoch photometric light curves of point sources
(e.g. blazars and other AGN) from Swift Ultraviolet/Optical Telescope (UVOT) archival data, in
all six filters (V, B, U, UVW1, UVM2, UVW2).

It is the companion to the Swift XRT pipeline,
[SwiftXRT_Complete](https://github.com/etmeyer/SwiftXRT_Complete), and follows the same
approach: one script per step, visual inspection at each stage, explicit include/exclude tables,
and thresholds checked on real data.

**Author:** Eileen T. Meyer ([@etmeyer](https://github.com/etmeyer))  
**License:** MIT

> **Status (October 2026):** planning is complete; the scripts are not written yet.
> The planned workflow is in [docs/index.md](docs/index.md).

---

## Why a dedicated pipeline

Aperture photometry of a bright, variable point source with UVOT has several traps that the
standard tools do not guard against by default, among them:

- **Coincidence loss and saturation.** Above 0.98 raw counts per frame, `uvotsource` caps the
  rate and sets a 100 % error. Bright sources reach this in full-frame exposures (11 ms frames)
  but not in hardware-windowed ones (e.g. 3.6 ms frames), and both often appear in one
  observation.
- **One exposure per FITS extension.** A UVOT sky-image file holds every exposure in a filter
  as a separate extension; `uvotsource` measures only one extension per call.
- **Pointing problems.** Exposures without aspect correction, exposures that start while the
  spacecraft is still settling (trailed images), and the spacecraft jitter of
  August 2023 – April 2024.
- **Small-scale sensitivity patches** on the detector, flagged by recent HEASoft versions.
- **Non-detections** in some filters or epochs, which need upper limits rather than magnitudes.

The pipeline measures every exposure, records why each one is used or excluded, and combines
the good ones into one point per observation and filter.

## Planned workflow

| Step | What |
| ---- | ---- |
| 1 | Setup and environment check (HEASoft, UVOT CALDB) |
| 2 | Download UVOT data from the HEASARC archive |
| 3 | Inventory of every exposure (filter, frame time, aspect, field of view) |
| 4 | Source positions, regions and image diagnostics |
| 5 | Photometry of every exposure with `uvotsource` |
| 6 | Quality rules and overrides → master table |
| 7 | Light curve: combine per observation and filter, plot |

## Requirements

- [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) ≥ 6.30 (developed with 6.36)
- The Swift UVOT [CALDB](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/swift/) (index 20240201)
- Python 3 with astropy, numpy, scipy, matplotlib, requests
