# Step 3 — Inventory of every exposure

`swift_uvot_inventory.py` opens every extension of every sky image from
[Step 2](02-download.md) and writes down what the later steps need to know
about it: what kind of exposure it is, where it sits in its observation,
whether it duplicates another one, and whether your source is inside its
exposed field. It gives every exposure a starting status — the start of the
pipeline's ledger, in which every exposure is accounted for until the light
curve. It changes nothing in `UVOT_input`.

It reads only FITS files, so it runs in either terminal (Python with
`astropy`; no HEASoft). All of 3C 273 — 2,935 exposures in 1,553 sky images
— takes under a minute with `--nproc 8`.

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
| `--no-catalog` | Skip the exposure-log check and the master catalog's exposures (no network needed) |

## How it works

**One extension per exposure.** Each sky image
`sw<OBSID>u<filter>_sk.img.gz` holds one extension per exposure in that
filter, named after the filter and the exposure's start, e.g.
`uu387398273I` (`I`: image mode; `E`: an image built from event-mode data).
The exposure map `*_ex.img.gz` has an extension of the same name. When the
archive split a long observation's event data into two event files
(`*_uf.evt.gz` and `*_uf_01.evt.gz`), there is a second sky image,
`*_sk_01.img.gz`, with its own exposure map; the inventory reads both
(3C 273 has one, in 00089771001: 5 more UVM2 exposures). Exposures in
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
| `loss` | Time lost to telemetry/DPU stalls/filter blocking (s); `-` if the keywords make no sense (one 3C 273 header has `STALLOSS` = 4×10²⁵²) |
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
After a slew the spacecraft reports itself settled while its pointing is
still converging by several arcseconds, so the first exposure of a snapshot
can catch that drift. Whether it shows depends on the data mode. In image
mode the instrument corrects the drift on board as it builds the image. In
event mode every photon is time-tagged, but the archive builds the sky image
with a single pointing for the whole exposure, so the drift is baked in as a
trail; the photons themselves are fine, and an image rebuilt from the event
file with the pointing at each photon's time is compact. Measured on all of
3C 273's candidates:

| Exposures | Number | Median axis ratio | Median offset | Trailed or offset* |
| --------- | ------ | ----------------- | ------------- | ------------------ |
| Event mode, opening a snapshot | 194 | 1.89 | 3.1″ | 80 % |
| Event mode, other | 531 | 1.07 | 0.6″ | 0.6 % |
| Image mode, opening a snapshot | 430 | 1.08 | 0.2″ | 5 % |
| Image mode, other | 1,736 | 1.07 | 0.3″ | 5 % |

\*axis ratio of the source above 1.3 or centroid more than 2″ from the
catalogue position (measured in [Step 4](04-positions.md)). In 2009–2016
each snapshot opened with a short UVW2 or UVW1 event-mode exposure. The
inventory marks openers with `first` and counts the event-mode ones; Step 4
measures every exposure's shape.

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

**The exposure-log check.** HEASARC's UVOT exposure log (`swiftuvlog`)
lists every exposure with its extension name. The inventory checks that each
logged exposure, in the filters you downloaded, has a sky image — Step 2's
check can't see a missing exposure, because the file exists. In 3C 273, 2,935
of 2,936 logged exposures have one; the exception, in 00050900031 UVW1, is
logged with the extension name `UNDEF` (1,082 s) and is in no sky or raw
image. The compact table also lists the master catalog's on-time per filter
(`cat_exp_s`). It is not a test: in 4 of 3C 273's 1,552 observation/filter
pairs it is 9–11 % above the sky images' on-time with nothing missing (once
because it counts the elapsed time of an `IMAGEEVENT` observation).

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
[summary] 26 candidates open a snapshot, 10 of them in event mode (often trailed; Step 4 measures them)
[summary] 0 extensions have an image/event duplicate
[note] 2 OBSID folder(s) have no sky images: 00035017005 00035017006
[check] exposure log (swiftuvlog): 122 of 123 logged exposures have a sky image
          missing: 00050900031 UVW1  UNDEF (1082 s), logged without an image
```

`aspcorr` counts the exposures: `D2` two `DIRECT`, `N2` two `NONE`. The two
OBSIDs without sky images are the BLOCKED-only and grism-only observations
Step 2 skipped as `no_data`.

## What to look for

- **`not_covered` and `partial`.** For 3C 273: 44 of 2,935 exposures — 42
  in pointings of SDSS J122933+015810, an X-ray source 8′ away (target IDs
  91742 and 32759; 3C 273 6–12′ off-axis), one 3C 273 observation pointed 8′
  away so that the source falls just outside the rotated field, and one
  0.01 s exposure (`exposure only 0.01 s`).
- **Two frame times in one filter** (`frame_ms` `3.6,11.0`). Step 5 measures
  each exposure; for a bright source the 3.6 ms ones are the ones that stay
  below saturation.
- **`aspcorr` `N`.** 40 % of 3C 273's candidates (1,151 of 2,891) have no
  field-star correction, nearly all in the UV filters and U. Step 4 measures
  where the source actually is in each.
- **`first` = `yes` with `mode` `EVENT`:** probably trailed (194 for
  3C 273, 80 % of them trailed or offset). Image-mode openers are fine.
- **`bkg_cover` below 1:** the standard background annulus runs off the
  exposed field (11 candidates for 3C 273); Step 4 chooses another
  background region there.
- **`missing:` lines:** exposures HEASARC logged that have no sky image.
  Nothing to re-download.

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

# Without network (no exposure-log check)
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

4. **The exposure-log check needs network.** Without it (`--no-catalog`) it
   is skipped and the `cat_exp_s` column stays empty.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_
