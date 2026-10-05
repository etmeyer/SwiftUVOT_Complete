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

> **Status (October 2026):** Steps 1–2 (setup, download) are done; the other
> steps are being written. Documentation: <https://etmeyer.github.io/SwiftUVOT_Complete/>,
> built from [docs/](docs/index.md).

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

## Workflow

| Step | What | Script | Status |
| ---- | ---- | ------ | ------ |
| 1 | Setup and environment check (HEASoft, UVOT CALDB) — [docs](docs/01-setup.md) | `swift_uvot_doctor.py` | done |
| 2 | Download UVOT data from the HEASARC archive — [docs](docs/02-download.md) | `swift_uvot_download.py` | done |
| 3 | Inventory of every exposure (filter, frame time, aspect, field of view) | | planned |
| 4 | Source positions, regions and image diagnostics | | planned |
| 5 | Photometry of every exposure with `uvotsource` | | planned |
| 6 | Quality rules and overrides → master table | | planned |
| 7 | Light curve: combine per observation and filter, plot | | planned |

Everything runs in **one terminal** with HEASoft set up and CIAO not set up:

```bash
setup_swiftuvot; heainit
swift_uvot_doctor.py          # Step 1: check the terminal
swift_uvot_download.py --name "3C 273" --list-only                # Step 2: list ...
swift_uvot_download.py --name "3C 273" --outdir UVOT_input --nproc 4   # ... and download
```

See [docs/01-setup.md](docs/01-setup.md) for installation, the UVOT CALDB, the
shell setup and what the doctor checks.

### `swift_uvot_doctor.py`

Checks the terminal: pipeline on `PATH`, HEASoft (6.30 or newer) without
CIAO, the UVOT CALDB (every file `uvotsource` uses, its sensitivity-loss
version, a short enough `$CALDB`), Python packages and free disk. Exits
non-zero if any check fails.

```
swift_uvot_doctor.py [--quiet] [--no-color] [--data-dir DIR]
                     [--test-image FILE[+EXT] [--ra RA --dec DEC]]
```

`--test-image` also runs `uvotsource` on one exposure of a UVOT sky image.

### `swift_uvot_download.py`

Finds the Swift observations of a source (name, position, OBSID or a file of
OBSIDs; optional date window), lists their UVOT exposure per filter, and
downloads the chosen filters' sky images, exposure maps, raw images and event
files plus housekeeping and `auxil` into `UVOT_input/<OBSID>/`. Files are
downloaded under a temporary name and checked before they replace anything;
each observation is checked against the catalog's exposure per filter;
`download_report.txt` records the outcome. Re-running resumes and repairs.
Exits 1 if any file failed or an OBSID was not found.

```
swift_uvot_download.py (--name NAME | --ra RA --dec DEC | --obsid ID | --obsid-file FILE)
                       [--start-date YYYY-MM-DD] [--end-date YYYY-MM-DD] [--radius ARCMIN]
                       [--filters vv bb uu w1 m2 w2] [--skip-raw] [--products ...]
                       [--outdir UVOT_input] [--nproc N] [--list-only] [--max-obs N]
                       [--overwrite] [--verify]
```

### `swift_uvot_runner.py` and `swift_uvot_env.py`

Modules the scripts share: `swift_uvot_runner.py` runs every HEASoft call in
isolation (own working folder, private parameter files, no terminal prompts,
timeout, log) and judges success by the outputs, not only the exit status;
`swift_uvot_env.py` holds the environment checks.

## Requirements

- [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) ≥ 6.30 (developed with 6.36)
- The Swift UVOT [CALDB](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/swift/) (index 20240201)
- Python 3 with astropy, numpy, scipy, matplotlib, requests
