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

> **Status (October 2026):** Steps 1–5 (setup, download, inventory, positions
> and regions, photometry) are done; the other steps are being written. Documentation: <https://etmeyer.github.io/SwiftUVOT_Complete/>,
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
| 3 | Inventory of every exposure (filter, frame time, aspect, field of view) — [docs](docs/03-inventory.md) | `swift_uvot_inventory.py` | done |
| 4 | Source positions, shapes and regions; trailed and smeared images, sensitivity patches — [docs](docs/04-positions.md) | `swift_uvot_positions.py`, `swift_uvot_viewer.py` | done |
| 5 | Photometry of every exposure with `uvotsource` — [docs](docs/05-photometry.md) | `swift_uvot_photometry.py` | done |
| 6 | Quality rules and overrides → master table | | planned |
| 7 | Light curve: combine per observation and filter, plot | | planned |

Everything runs in **one terminal** with HEASoft set up and CIAO not set up:

```bash
setup_swiftuvot; heainit
swift_uvot_doctor.py          # Step 1: check the terminal
swift_uvot_download.py --name "3C 273" --list-only                # Step 2: list ...
swift_uvot_download.py --name "3C 273" --outdir UVOT_input --nproc 4   # ... and download
swift_uvot_inventory.py --ra 187.2779 --dec 2.0524                      # Step 3: every exposure
swift_uvot_positions.py --ra 187.2779 --dec 2.0524 --nproc 16           # Step 4: positions, regions
swift_uvot_viewer.py                                                    # ... and contact sheets
swift_uvot_photometry.py --nproc 16                                     # Step 5: photometry
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

### `swift_uvot_inventory.py`

Opens every extension (one per exposure) of every sky image in `UVOT_input`
and writes `UVOT_output/uvot_inventory.txt`: filter, data mode, frame time,
window, aspect correction, exposure and on-time, snapshot and whether the
exposure opens it, image/event duplicates, and how much of the 5″ source
circle and 27.5–35″ background annulus is exposed. Every exposure gets one
status (`candidate`, `not_covered`, `partial`, `settling`, `nonphot`,
`no_expmap`, `unreadable`). Also checks that every exposure in HEASARC's
UVOT exposure log has a sky image. Exits 1 if any file is unreadable.

```
swift_uvot_inventory.py --ra RA --dec DEC [--indir UVOT_input] [--outdir UVOT_output]
                        [--nproc 8] [--no-catalog]
```

### `swift_uvot_positions.py`

For every candidate exposure: the source's centroid and offset from your
position, its shape (axis ratio, sizes, position angle), the share of its
15″ counts inside the 5″ region (trailed and smeared images), its detector
position and the CALDB low-sensitivity flags at the LOW, MID and HIGH levels
(the same lookup as `uvotsource`). Writes fk5 source and background region
files: the 5″ circle (at the centroid when the source is clear and within
3″), and the 27.5–35″ annulus minus the sources `uvotdetect` finds, or an
offset circle when the annulus is not evenly exposed. Positions and
backgrounds set in `UVOT_output/uvot_region_overrides.txt` win. Output:
`UVOT_output/uvot_positions.txt`, `uvot_detections.txt`, `regions/`. Exits 1
if any exposure failed.

```
swift_uvot_positions.py --ra RA --dec DEC [--outdir UVOT_output] [--nproc 8]
```

### `swift_uvot_viewer.py`

A PDF of Step 4's results (`UVOT_output/uvot_positions.pdf`): distributions
of shape, concentration and offset, then contact sheets of the exposures
worth a look, worst first, with the regions drawn on; or every exposure of
the observations you name.

```
swift_uvot_viewer.py [--outdir UVOT_output] [--pdf FILE] [--per-category 48]
                     [--obsid OBSID ...]
```

### `swift_uvot_photometry.py`

Runs `uvotsource` on every exposure Step 4 passed, one extension per call
(with the exposure map, `apercorr=NONE`, `sigma=3`, `forcephot=no`,
`history=no`), and keeps everything: raw counts and areas, counts per frame
and the saturation flag, every correction factor, the corrected rate (also
on a low-sensitivity patch, where `uvotsource` gives no magnitude), limits,
magnitudes, flux densities and photometry flags. Checks the detector
position and the low-sensitivity flag against Step 4. Output:
`UVOT_output/uvot_photometry.txt`, `uvot_photometry.pdf` (every exposure's
rate against time, per filter) and `photometry/<OBSID>/` (each call's FITS
table and log). Exits 1 if any call failed.

```
swift_uvot_photometry.py [--outdir UVOT_output] [--nproc 8] [--plot-only]
```

### `swift_uvot_runner.py`, `swift_uvot_env.py`, `swift_uvot_tables.py`

Modules the scripts share: `swift_uvot_runner.py` runs every HEASoft call in
isolation (own working folder, private parameter files, no terminal prompts,
timeout, log) and judges success by the outputs, not only the exit status;
`swift_uvot_env.py` holds the environment checks; `swift_uvot_tables.py` reads
and writes the pipeline's aligned text tables.

## Requirements

- [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) ≥ 6.30 (developed with 6.36)
- The Swift UVOT [CALDB](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/swift/) (index 20240201)
- Python 3 with astropy, numpy, scipy, matplotlib, requests
