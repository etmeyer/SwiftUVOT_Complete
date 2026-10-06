# Step 5 — Photometry of every exposure

`swift_uvot_photometry.py` runs `uvotsource` on every exposure that
[Step 4](04-positions.md) passed, one exposure per call, with Step 4's source
and background regions and the exposure map. It keeps every number and flag
`uvotsource` reports. Nothing is excluded here: [Step 6](06-master-table.md) applies the
quality rules to these numbers, so a rule can be changed and re-applied
without running HEASoft again.

Run it in the HEASoft terminal (`setup_swiftuvot; heainit`), from the
analysis folder. For all of 3C 273 (2,891 exposures) it takes about 4½
minutes with `--nproc 16`.

## What runs

```bash
swift_uvot_photometry.py --nproc 16
```

| Flag | Purpose |
| ---- | ------- |
| `--outdir` | Step 3–4's output folder (default `UVOT_output`) |
| `--nproc` | Parallel `uvotsource` calls (default 8) |
| `--plot-only` | Only remake `uvot_photometry.pdf` from the tables (no HEASoft needed) |

## How it works

**One call per exposure.** Each call measures one extension of a sky image,
named explicitly (`image=sky.img.gz[uu387398273I]`, and the same for the
exposure map). Given no extension, `uvotsource` measures only the first
exposure in the file. Each call runs in its own folder through the
pipeline's runner. That matters because `uvotsource` appends to an existing
output table instead of replacing it. Each call's FITS table and log are
kept in `photometry/<OBSID>/`.

**The settings,** and why:

| Parameter | Value | Why |
| --------- | ----- | --- |
| `srcreg`, `bkgreg` | Step 4's regions | 5″ circle; 27.5–35″ annulus with exclusions, or an offset circle |
| `expfile` | the exposure map | `uvotsource` then refuses regions that are not fully or evenly exposed |
| `apercorr` | `NONE` | the 5″ aperture is the calibrated one |
| `sigma` | 3 | for the background-limited rate and magnitude |
| `syserr` | `no` | errors are statistical; the systematic (zero-point) error is kept in its own column |
| `forcephot`, `skipbad` | `no`, `no` | `uvotsource`'s own checks stay on, and a refusal becomes a failure with its reason |
| `centroid` | `no` | Step 4 placed the regions |
| `history` | `no` | `fthedit` can crash writing history keywords when CALDB paths are long |

**The corrections.** `uvotsource` turns the counts in the 5″ circle into a
rate, subtracts the background, and corrects in turn for:

- coincidence loss (`coi`);
- the large-scale sensitivity of the detector (`lss`);
- the small-scale sensitivity (`sss`);
- the loss of sensitivity since launch (`senscorr`).

The result is the corrected rate, from which come the magnitudes and flux
densities. The table keeps every factor, and the raw counts and areas too, so
that exposures can later be stacked from counts.

**Coincidence loss and saturation.** When more than about one photon per
frame arrives in the same spot, the detector counts them as one.
`uvotsource` corrects this from the raw counts per frame in the 5″ circle:
rate × (frame time − 0.174 ms dead time), column `counts_frame`. Above 0.98
counts per frame it caps the value, sets `saturated` = 1 and gives the rate a
100 % error. The rate is then only a lower limit. Exposures taken with a
hardware window have 3.6 ms frames instead of 11 ms, so a bright source stays
far below the limit in them.

**Low-sensitivity patches.** `uvotsource` reads the LOW small-scale
sensitivity map at the source's detector position, as Step 4 did for all
three levels. On a patch it reports no magnitude (`mag` 99; its own
`CORR_RATE` is −999). The rate is still computed, though: the table's `rate`
is `uvotsource`'s fully corrected `SENSCORR_RATE`, filled in either way.
Step 6 can then still compare such exposures with their neighbours. The
flag is in `photflags` (`BAD_SSS`).

**Statuses.** `ok` (measured) or `failed` (with `uvotsource`'s own error
lines and the log's path). Exposures Step 4 did not pass keep their Step 4
status and are not measured. Every Step 4 row comes back exactly once.

**Checks against Step 4.** `uvotsource` computes the detector position and
the low-sensitivity flag independently; both must agree with Step 4's.

For all of 3C 273:

```
[summary] 2891 exposures: ok 2891
          filter  measured  saturated  >0.95 counts/frame  LOW patch  nsigma<3
          V            487          0                   0         25         1
          B            196          0                  40          4         1
          U            504        176                 314         54         1
          UVW1         528          0                   3         99         1
          UVM2         565          0                   0        105         2
          UVW2         611          0                  18         98         1
[summary] photometry flags: BAD_SSS 385
[check] low-sensitivity flag: uvotsource agrees with Step 4 for 2891 of 2891 exposures
[check] detector position: uvotsource and Step 4 agree to 0.05 pixel (largest difference)
```

(Step 4 writes the detector position to 0.1 pixel.)

### The plot

`uvot_photometry.pdf` has a first page with the distribution of counts per
frame for each filter, against the saturation limit, and a bar chart of the
exposures by category. Then comes one page per filter: the rate of every
exposure against time. Each exposure is marked by the first of these that
applies: saturated, trailed or smeared (Step 4: concentration below 0.70 or
axis ratio above 1.3), faint at the source position (Step 4 S/N below 10),
on a LOW low-sensitivity patch, or none of these. Nothing has been excluded
yet, so this is where to see what Step 6 will remove.

## What to look for

- **Saturated U.** 176 of 3C 273's 358 full-frame U exposures are
  saturated, and none of the 146 with 3.6 ms frames. In the plot they form a
  flat band near the cap (370–420 counts/s), while unsaturated U exposures
  from the same years read up to 520. The median counts per frame
  is 0.96 in U, 0.88 in B, 0.80 in UVW1, 0.71 in UVW2, 0.64 in UVM2 and 0.55
  in V.
- **What each flag does to the rate:** see the table below.
- **Low significance.** Seven exposures have `nsigma` below 3, all of them
  the jitter-period exposures Step 4 found faint at the source position.
  Three have a negative rate. Negative rates are kept as they are: they
  matter for faint sources.
- **Failures.** None for 3C 273. A failure's reason quotes `uvotsource`,
  e.g. `error: low bkg.reg in FOV 0.000` for a background region off the
  image.

### What each flag does to the rate

A first look; Step 6 calibrates the rules. Each flagged exposure is compared
with the median rate of the clean exposures in the same observation and
filter. Clean means: not trailed or smeared, S/N ≥ 10, on no patch at any
level, not saturated, at most 0.95 counts per frame, and outside the jitter
period. The table gives the flagged exposure's rate as a fraction of that
median, for 3C 273:

| Exposures | Median rate / clean median |
| --------- | -------------------------- |
| On a LOW patch | 0.95 (UVM2 and UVW2 0.92–0.93; 29 % more than 10 % low) |
| On a MID patch but not LOW | 0.99 |
| On a HIGH patch but not MID | 1.00 |
| Concentration below 0.70 (Step 4) | 0.31 |
| Axis ratio above 1.3, concentration 0.70 or more | 0.92 |
| Saturated | 0.85 |
| U at 0.95–0.98 counts per frame, not saturated | 0.92 |

The U row compares against exposures below 0.85 counts per frame. These
numbers support excluding exposures that are on a LOW patch, trailed,
smeared or saturated, and a counts-per-frame limit for U below the
saturation cap. Step 6 sets the limits.

## Inputs and outputs

**Inputs:** Step 4's `uvot_positions.txt` and region files; Step 3's
`uvot_inventory.txt` (which file holds each exposure); the sky images and
exposure maps; the UVOT CALDB.

**Outputs** (in `UVOT_output`):

```
UVOT_output/
├── uvot_photometry.txt      one row per exposure (columns below)
├── uvot_photometry.pdf      counts per frame, categories, rate against time per filter
└── photometry/<OBSID>/
    ├── <EXTNAME>.fits       uvotsource's full output table (one row, 126 columns)
    └── <EXTNAME>.log        the call: command, exit status, complete output
```

For 3C 273, `photometry/` takes 126 MB.

| Column | What it is |
| ------ | ---------- |
| `exposure`, `frametime` | Exposure (s, dead-time corrected) and frame time (s) |
| `counts_frame`, `saturated` | Raw counts per frame in the 5″ circle; 1 if above 0.98 (the rate is then a lower limit) |
| `src_cnts`, `bkg_cnts`, `src_area`, `bkg_area` | Raw counts in the source and background regions, and their areas (arcsec²) |
| `coi`, `coi_bkg`, `lss`, `sss`, `senscorr` | Correction factors: coincidence loss (source, background), large- and small-scale sensitivity (`sss` −99.9 on a LOW patch), sensitivity loss |
| `rate`, `rate_err` | Corrected rate (counts/s) and its statistical error; filled also on a LOW patch |
| `rate_lim`, `mag_lim` | 3σ background-limited rate and Vega magnitude |
| `nsigma` | `uvotsource`'s significance |
| `mag`, `mag_err`, `mag_sys` | Vega magnitude, statistical error, systematic error (0.01 V, 0.02 B and U, 0.03 UV); 99 on a LOW patch or for a negative rate |
| `flux_aa`, `flux_aa_err` | Flux density (erg s⁻¹ cm⁻² Å⁻¹) and its statistical error |
| `photflags` | `uvotsource`'s photometry flags, e.g. `BAD_SSS` (`NO_QUALITY_MAP`, always set, is left out) |
| `phot_file` | The FITS table, with everything else (AB magnitudes, F_ν, background, limits) |
| `status`, `reason` | `ok` or `failed`, or Step 4's status; notes such as `saturated` or the error |

## Common variants

```bash
# Remake the plot after changing it, without HEASoft
swift_uvot_photometry.py --plot-only

# The saturated exposures
grep "lower limit" UVOT_output/uvot_photometry.txt | less -S

# Everything uvotsource said about one exposure
less UVOT_output/photometry/00035017124/uu387398273I.log
```

## Gotchas

1. **Re-run after Step 4.** The regions come from Step 4; after you change an
   override there, re-run this step.

2. **Re-running replaces** `uvot_photometry.txt` and each exposure's FITS
   table and log. An exposure Step 4 no longer passes keeps its old files in
   `photometry/`, but the table no longer points to them.

3. **`mag` 99 does not mean a non-detection.** It marks a LOW patch or a
   negative rate. Use `rate` and `rate_err`; upper limits are made in
   Step 7 from the combined rates.

4. **`nsigma` is `uvotsource`'s significance,** not a signed S/N. Step 7
   decides detections from `rate`/`rate_err` after combining exposures.

5. **A failed call keeps its working folder** in `/tmp` (named in the log)
   for inspection.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_
