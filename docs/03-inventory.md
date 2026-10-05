# Step 3 — Inventory of every exposure

`swift_uvot_inventory.py` opens every extension of every sky image from
[Step 2](02-download.md) and writes down what the later steps need to know
about it: what kind of exposure it is, where it sits in its observation,
whether it duplicates another one, and whether your source is inside its
exposed field. It gives every exposure a starting status — the start of the
pipeline's ledger, in which every exposure is accounted for until the light
curve. It changes nothing in `UVOT_input`.

It reads only FITS files, so it runs in either terminal (Python with
`astropy`; no HEASoft). All of 3C 273 — 2,930 exposures in 1,552 sky images
— takes about 40 seconds with `--nproc 8`.

## What runs

```bash
swift_uvot_inventory.py --ra 187.2779 --dec 2.0524
```

| Flag | Purpose |
| ---- | ------- |
| `--ra` / `--dec` | Your source position, J2000 decimal degrees (required) |
| `--indir` | Step 2's folder (default `UVOT_input`) |
| `--outdir` | Where the tables go (default `UVOT_output`) |
| `--nproc` | Files read at once (default 8) |
| `--no-catalog` | Skip the comparison with the Swift master catalog (no network needed) |

## How it works

**One extension per exposure.** Each sky image
`sw<OBSID>u<filter>_sk.img.gz` holds one extension per exposure in that
filter, named after the filter and the exposure's start, e.g.
`uu387398273I` (`I`: image mode; `E`: an image built from event-mode data).
The exposure map `*_ex.img.gz` has an extension of the same name. Exposures in
one file can differ in frame time, window, data mode and aspect correction:
in 3C 273's 2009–2016 observations, each filter has a short exposure with
11.0 ms frames followed by a longer one with 3.6 ms frames (a hardware
window), which is what keeps 3C 273 below saturation.

**What it records** for each extension (one row of `uvot_inventory.txt`):

| Column | What it is |
| ------ | ---------- |
| `obsid`, `filter`, `hdu`, `extname` | Which exposure (`hdu` is the extension number, as in `file.img.gz[2]`) |
| `mode`, `obs_mode` | `IMAGE`, `EVENT` or `IMAGEEVENT`; `POINTING` or `SETTLING` |
| `expid` | Exposure ID (shared by an image-mode and event-mode copy of one exposure) |
| `date_mid`, `mjd_mid` | Middle of the exposure (UTC, from `DATE-OBS`/`DATE-END`, as in the XRT pipeline) |
| `tstart`, `tstop` | Start and stop (Swift mission time, s) |
| `exposure`, `ontime` | Exposure with and without the dead-time correction (s) |
| `frametime`, `window`, `binning` | CCD frame time (s), science window (raw pixels), on-board binning |
| `aspcorr` | `DIRECT` if the pointing was corrected with field stars, `NONE` if not |
| `loss` | Time lost to telemetry/DPU stalls/filter blocking (s) |
| `snapshot`, `first` | Snapshot number (exposures less than 5 min apart) and whether this exposure opens it |
| `off_pnt` | Pointing offset from your source (arcmin) |
| `src_cover`, `bkg_cover` | Fraction of the 5″ source circle and of the 27.5–35″ background annulus that is fully exposed |
| `clearance` | Distance from the source to the nearest unexposed pixel (arcsec; `>60`) |
| `dup_of` | The other extension with the same `expid`, if any |
| `event` | Whether the event file is there (event-mode exposures) |
| `procver` | Archive processing version |
| `status`, `reason` | See below |

**Coverage.** A pixel counts as exposed if the exposure map is within 1 % of
its maximum near the source, the same tolerance `uvotsource` uses before it
refuses an uneven region. Pixels outside the image count as unexposed. In
pointings of other targets, or with small windows, your source can be partly
or wholly outside an exposure.

**Snapshots.** Swift observes in snapshots of up to ~30 minutes per orbit.
The first exposure of a snapshot can start while the spacecraft is still
settling; for event-mode exposures the archive's sky image is made with a
single pointing, so the drift shows as a trailed source (in 3C 273's
2009–2016 data, the UVW2 exposure that opens each snapshot). Step 4 looks at
the image shapes; here the exposure is only marked with `first`.

**Statuses.** Every extension gets exactly one:

| Status | Meaning |
| ------ | ------- |
| `candidate` | Goes on to positions and photometry (Steps 4–5) |
| `not_covered` | The source is outside the exposed field |
| `partial` | Part of the 5″ source circle is not fully exposed |
| `settling` | `OBS_MODE` is `SETTLING` |
| `nonphot` | WHITE, a grism, the magnifier or BLOCKED (only if you downloaded them) |
| `no_expmap` | No exposure-map extension with the same name |
| `unreadable` | The file or extension could not be read — re-download it (Step 2) |

Saturation, pointing errors, trailing and the detector's low-sensitivity
patches are not judged here; they need the photometry and images of
Steps 4–6. The inventory deliberately measures no brightness, so each
quantity comes from one place.

**The catalog check.** For each observation and filter, the summed on-time
(before dead-time correction, as the catalog counts it) is compared with the
Swift master catalog. A difference of more than 2 % usually means an exposure
the archive lists but has no sky image for. Step 2's check can't see that,
because the file exists. In 3C 273, six observation/filter pairs differ; for
00050900031 UVW1, HEASARC's exposure log lists a tenth exposure of 1,082 s
(extension name `UNDEF`) that is in no sky image.

The run prints one row per observation and filter (also saved as
`uvot_inventory_compact.txt`), then counts. From the 12-observation test set:

```
  obsid        date        filter  n_ext  n_cand  exp_s  cand_exp_s  ontime_s  cat_exp_s  frame_ms  aspcorr  statuses
  ...
  00035017124  2013-04-11  UVW2    2      2       145    145         152       152        3.6,11.0  N2       candidate:2
  00035017124  2013-04-11  V       2      2       145    145         152       152        3.6,11.0  D2       candidate:2
  00050900031  2026-01-15  UVW1    9      9       8706   8706        8846      9928       11.0      D9       candidate:9
  00091742011  2014-01-04  UVW2    2      1       2938   1691        2985      2985       11.0      D2       candidate:1,not_covered:1

[summary] 122 extensions: candidate 121, not_covered 1
          ...
[summary] 26 candidates open a snapshot; 0 extensions have an image/event duplicate
[note] 2 OBSID folder(s) have no sky images: 00035017005 00035017006
[check] on-time per OBSID and filter vs the catalog: 1 of 44 differ by more than 2% (and 10 s)
          00050900031 UVW1  inventory 8846 s, catalog 9928 s
          (usually an exposure the archive has no sky image for)
```

`aspcorr` counts the exposures: `D2` two `DIRECT`, `N2` two `NONE`. The two
OBSIDs without sky images are the BLOCKED-only and grism-only observations
Step 2 skipped as `no_data`.

## What to look for

- **`not_covered` and `partial`.** For 3C 273: 44 of 2,930 exposures — 42
  in pointings of the neighbouring AGN SDSS J122933+015810 (target IDs 91742
  and 32759, 3–12′ from 3C 273), one 3C 273 observation pointed 8′ away so
  that the source falls just outside the rotated field, and one 0.01 s
  exposure (`exposure only 0.01 s`).
- **Two frame times in one filter** (`frame_ms` `3.6,11.0`). Step 5 measures
  each exposure; for a bright source the 3.6 ms ones are the ones that stay
  below saturation.
- **`aspcorr` `N`.** 40 % of 3C 273's candidates (1,147 of 2,886) have no
  field-star correction, nearly all in the UV filters and U. Step 4 measures
  where the source actually is in each.
- **`first` = `yes` for event-mode exposures:** possible trailed images
  (189 for 3C 273).
- **`bkg_cover` below 1:** the standard background annulus runs off the
  exposed field (11 candidates for 3C 273); Step 4 chooses another
  background region there.
- **`[check]` lines:** exposures the archive lacks. Nothing to re-download.

## Inputs and outputs

**Inputs:** Step 2's `UVOT_input` (sky images and exposure maps), its
`download_report.txt` (OBSIDs whose download failed are listed), and your
source position.

**Outputs** (in `--outdir`):

```
UVOT_output/
├── uvot_inventory.txt           one row per exposure: the columns above
└── uvot_inventory_compact.txt   one row per observation and filter
```

Both are plain aligned text tables with `#` comment lines (command, source
position, pipeline version, counts).

## Common variants

```bash
# The usual run, from the analysis folder that holds UVOT_input
swift_uvot_inventory.py --ra 187.2779 --dec 2.0524

# Without network (no catalog comparison)
swift_uvot_inventory.py --ra 187.2779 --dec 2.0524 --no-catalog

# Look at one observation
grep 00035017124 UVOT_output/uvot_inventory.txt | less -S
```

## Gotchas

1. **Use the same `--ra`/`--dec` in every step.** Coverage is computed for
   that position.

2. **Re-run after changing `UVOT_input`** (more observations, a re-download):
   the inventory describes the files as they were when it ran.

3. **The exit status is 1 if anything was unreadable.** The `reason` column
   says what; re-run Step 2 for those observations (`--verify` replaces
   damaged files). Each `.gz` file is decompressed completely before it is
   read, because astropy shows a truncated file as one with fewer
   extensions — exposures would disappear without a word.

4. **The catalog comparison needs network.** Without it (`--no-catalog`) the
   `cat_exp_s` column stays empty.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_
