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

> **Status (November 2026):** Steps 1–8 (setup, download, inventory, positions
> and regions, photometry, quality rules, light curve, dereddened flux
> densities) are done; tested on 3C 273 and on the faint, distant
> PKS 0537−286. Documentation: <https://etmeyer.github.io/SwiftUVOT_Complete/>,
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
| 4 | Source positions, shapes and regions; trailed and smeared images, sensitivity patches; field stars (pointing, PSF) — [docs](docs/04-positions.md) | `swift_uvot_positions.py`, `swift_uvot_viewer.py` | done |
| 5 | Photometry of every exposure with `uvotsource` — [docs](docs/05-photometry.md) | `swift_uvot_photometry.py` | done |
| 6 | Quality rules and overrides → master table — [docs](docs/06-master-table.md) | `make_uvot_master_table.py` | done |
| 7 | Light curve: combine per observation and filter, plot — [docs](docs/07-lightcurve.md) | `swift_uvot_lightcurve.py` | done |
| 8 | Flux densities: dereddened, intergalactic absorption, your spectrum — [docs](docs/08-flux.md) | `swift_uvot_flux.py` | done |

Everything runs in **one terminal** with HEASoft set up and CIAO not set up:

```bash
setup_swiftuvot; heainit
swift_uvot_doctor.py          # Step 1: check the terminal
swift_uvot_download.py --name "3C 273" --list-only                # Step 2: list ...
swift_uvot_download.py --name "3C 273" --outdir UVOT_input --nproc 4   # ... and download
swift_uvot_inventory.py --ra 187.2779 --dec 2.0524                      # Step 3: every exposure
swift_uvot_positions.py --ra 187.2779 --dec 2.0524                     # Step 4: positions, regions
swift_uvot_viewer.py                                                    # ... and contact sheets
swift_uvot_photometry.py                                                # Step 5: photometry
make_uvot_master_table.py                                               # Step 6: quality rules
swift_uvot_lightcurve.py                                                # Step 7: light curve
swift_uvot_flux.py                                                      # Step 8: dereddened fluxes
```

See [docs/01-setup.md](docs/01-setup.md) for installation, the UVOT CALDB, the
shell setup and what the doctor checks.

### `swift_uvot_doctor.py`

Checks the terminal: pipeline on `PATH`, HEASoft (6.30 or newer) without
CIAO, the UVOT CALDB (every file `uvotsource` uses, its sensitivity-loss
version, a short enough `$CALDB`), Python packages, the Gaia archive and the
IRSA dust service, and free disk. Exits
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
                        [--nproc N] [--no-catalog]
```

### `swift_uvot_positions.py`

For every candidate exposure: the source's centroid and S/N; at S/N ≥ 10
its offset from your position and shape (axis ratio, sizes, position
angle); the share of its 15″ counts inside the 5″ region (trailed and
smeared images; only where that share is known to ±0.02); its detector
position and the CALDB low-sensitivity flags at the LOW, MID and HIGH levels
(the same lookup as `uvotsource`), and how much of the 5″ circle lies on a
LOW patch. Writes fk5 source and background region files: the 5″ circle (at
the centroid when the source is clear and within 3″), and the 27.5–35″
annulus minus the sources `uvotdetect` finds (each masked out to where its
PSF wing falls below 0.5 % of the sky; saturated ones also lose a read-out
strip) and minus any 45° sector that is too bright, or an offset circle
when too little of the annulus is left or it is not evenly exposed. Its own
background is the plain mean of what is left, as in `uvotsource`. Gaia DR3
field stars (fetched once from the ESA archive and cached) give each
exposure's pointing, which moves a faint source's regions when the image's
coordinates are off (or marks it not covered), and its PSF: the light a
broad or smeared image pushes out of the 5″ circle. Positions and
backgrounds set in `UVOT_output/uvot_region_overrides.txt` win. Output:
`UVOT_output/uvot_positions.txt`, `uvot_detections.txt`,
`uvot_field_stars.txt`, `uvot_gaia.txt`, `regions/`. Exits 1 if any
exposure failed.

```
swift_uvot_positions.py --ra RA --dec DEC [--outdir UVOT_output] [--nproc N]
                        [--no-field-stars] [--gaia-offline]
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
swift_uvot_photometry.py [--outdir UVOT_output] [--nproc N] [--plot-only]
```

### `make_uvot_master_table.py`

Joins Steps 3–5 into one row per exposure and decides by rule which ones the
light curve uses (`include` yes/no, with the rules that excluded each):
saturated, above 0.95 counts per frame (U: kept up to saturation but flagged above 0.90),
on a LOW low-sensitivity patch, source off the region centre, smeared
(concentration < 0.82), elongated (axis ratio > 1.5), no PSF check or
trusted pointing during the 2023–24 spacecraft jitter, duplicate exposure.
For a faint source the field stars' PSF loss is a flag (`psf_broad`;
`--psf-max-loss` to exclude). The limits were
calibrated on 3C 273. The table is regenerated, never edited; decisions go
in `UVOT_output/uvot_overrides.txt`. Lists the exposures whose `include`
changed since the last run. Output: `UVOT_output/uvot_master_table.txt` and
`uvot_master_table.pdf` (the rules against the data; the field stars' PSF
loss and pointing; every exposure's rate against time, included or
excluded).

```
make_uvot_master_table.py [--outdir UVOT_output] [--sss-level LOW|MID|HIGH]
                          [--coincidence-limit FILTER=VALUE ...]
                          [--psf-max-loss FRACTION]
```

### `swift_uvot_lightcurve.py`

Combines the included exposures into one point per observation and filter,
with weights from the counts each exposure is expected to have (unbiased for
faint sources; the same as inverse-variance weights for bright ones), a 1.5 %
floor on the combined rate (measured from 3C 273) and an additive background
systematic per filter (measured on blank sky). Detection from the Poisson
probability of the summed counts; upper limits by Kraft, Burrows & Nousek
(99.87 %); a 3σ sensitivity column. Statuses: detected, upper limit, lower
limit (all saturated), excluded, not covered. Count rates, Vega and AB
magnitudes and flux densities with `uvotsource`'s own conversion; flags
carried from the exposures; optionally corrected for the field stars' PSF
loss (`--psf-correct`). Output: `UVOT_output/uvot_lightcurve.txt` and
`uvot_lightcurve.pdf`.

```
swift_uvot_lightcurve.py [--outdir UVOT_output] [--floor 0.015] [--nsigma 3]
                         [--bkg-sys FILTER=VALUE,...] [--psf-correct]
```

### `swift_uvot_flux.py`

Turns the light curve's count rates into flux densities corrected for
Galactic dust and, with `--z`, for the intergalactic medium (Inoue et al.
2014 or Madau 1995). E(B−V) comes from IRSA (Schlafly & Finkbeiner 2011 by
default), cached once per analysis. The extinction (Fitzpatrick 1999 or
Cardelli et al. 1989) and the absorption are folded through each filter's
CALDB effective area. By default `uvotsource`'s factors are kept;
`--spectrum` folds a power law, blackbody or your own spectrum for the
factor too. Rows with strong intergalactic absorption or mostly red-leak
counts are flagged, not dropped. Needs no HEASoft. Output:
`UVOT_output/uvot_flux.txt` and `uvot_flux.pdf` (light curves in νF_ν and
spectral energy distributions).

```
swift_uvot_flux.py [--outdir UVOT_output] [--ebv E] [--ebv-map sf11|sfd]
                   [--law f99|ccm89] [--rv 3.1] [--z Z] [--igm inoue14|madau95|none]
                   [--spectrum powerlaw:ALPHA|blackbody:T|file:PATH]
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
