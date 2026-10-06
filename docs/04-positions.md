# Step 4 — Positions, shapes and regions

`swift_uvot_positions.py` looks at every candidate exposure from
[Step 3](03-inventory.md). It finds where your source actually is in the
exposure, what shape it has, how much of its light falls in the 5″ source
region, and whether it sits on one of the detector's low-sensitivity
patches. Then it writes the source and background regions that photometry
([Step 5](05-photometry.md)) uses. `swift_uvot_viewer.py` shows the exposures worth
looking at, with their regions drawn on.

Run both in the HEASoft terminal (`setup_swiftuvot; heainit`), from the
analysis folder that holds `UVOT_input` and `UVOT_output`. For all of
3C 273 (2,891 candidate exposures) the positions take about 4½ minutes with
`--nproc 16`, and the viewer about 20 seconds.

## What runs

```bash
swift_uvot_positions.py --ra 187.2779 --dec 2.0524 --nproc 16
swift_uvot_viewer.py                    # -> UVOT_output/uvot_positions.pdf
```

| Flag | Purpose |
| ---- | ------- |
| `--ra` / `--dec` | Your source position; must match the one Step 3 used |
| `--outdir` | Step 3's output folder (default `UVOT_output`) |
| `--nproc` | Parallel workers (default 8) |

`swift_uvot_viewer.py` takes `--outdir`, `--pdf`, `--per-category` (cutouts
per category, default 48) and `--obsid` (every exposure of the given
observations instead of the categories).

## How it works

**Position and shape.** In each exposure the script finds the brightest pixel
within 8″ of `--ra`/`--dec`, centroids the source around it (4″), and
measures its shape from the second moments within 6″: the axis ratio, the
sizes along the major and minor axes, and the major axis' position angle.
For 3C 273 the median axis ratio is 1.07 (95 % of ordinary exposures are
below 1.2), and 1.9 for the trailed event-mode snapshot openers of Step 3.

**Concentration.** The axis ratio sees only the central 6″. A trail longer
than that, or a smeared image with several lobes, can look round there and
still lose much of its light outside the 5″ source region. So the script also
divides the net counts in the 5″ region by the net counts in a 15″ circle
around it (column `conc`). For a point source this is about 0.85: the
median for 3C 273's ordinary exposures is 0.84–0.87 depending on the filter.
Strong coincidence loss flattens the core and lowers it a little: for U
exposures above 0.9 counts per frame the median is 0.83, and 5 % are below
0.77. Below 0.70 the reason column says
`trailed or smeared`. `conc` is left empty when even the 15″ circle holds no
clear source (S/N below 10 there).

**The source region** is a 5″ circle — the calibrated UVOT aperture — at the
centroid when the source is clearly detected (S/N ≥ 10) within 3″ of
`--ra`/`--dec`, otherwise at `--ra`/`--dec`. Centring on the source corrects
the pointing of exposures without aspect correction. For 3C 273, their
centroids are a median 0.77″ from the true position (90 % within 2.8″),
against 0.25″ (90 % within 0.45″) for aspect-corrected ones. A centroid
farther than 3″ away usually means a trailed image. The region then stays at
`--ra`/`--dec`; such an exposure is judged in Step 6, not rescued here.

**Low-sensitivity patches.** The script converts the source's sky position to
a detector position exactly as `uvotsource` does. It uses the sky image's
detector WCS, HEASoft's DET-to-RAW polynomial and the shift-and-add offset
recorded in newer processing (`UD_RAWX`/`UD_RAWY`). It then reads the
CALDB's small-scale-sensitivity maps at all three screening levels: LOW
(HEASoft's default), MID and HIGH. `-99.9` means the source sits on a patch
at that level. On 1,457 3C 273 exposures the detector positions agree with
`uvotsource`'s to better than 0.001 pixel, and the LOW flags agree in every
case. Which level to exclude at is decided in Step 6, from the data.

**The background region** is the 27.5–35″ annulus used in the UVOT
calibration (Poole et al. 2008). Every other source that `uvotdetect`
finds (2.5σ) in the deepest exposure of that observation and filter is
masked out with a circle sized by its magnitude, as in HEASoft's `uvotlc`:
20″ for magnitude ≤ 12, 15″ (≤ 14), 10″ (≤ 16), 7″ (≤ 18), otherwise 5″.
Your source itself is never masked (the brightest detection within 10″ of
`--ra`/`--dec`). `uvotsource` refuses a background that is not evenly
exposed, so the annulus must be fully exposed (within 1 %) and keep at least
half its area. Otherwise the script tries a 15″ circle 45″, 60″ or 90″ away,
wherever it is fully exposed and clear of sources. If none qualifies, the
status is `no_background` (you can set one by hand, below).

**Statuses.** Every candidate exposure gets one: `ok`, `no_background`, or
`failed` (with the reason). Nothing is excluded here; the flags and numbers
go to Step 6, which selects. Exposures that Step 3 did not pass keep their
Step 3 status.

For all of 3C 273:

```
[summary] 2891 candidate exposures: ok 2891
[summary] source region at: catalog 111, centroid 2780
[summary] background: annulus 2879, circle 12
          offsets, ASPCORR=DIRECT: median 0.25", 90% 0.45", 1 over 3"
          offsets, ASPCORR=NONE  : median 0.77", 90% 2.84", 110 over 3"
[summary] axis ratio > 1.3 (trailed or doubled?): 223, of which 146 event-mode snapshot openers
[summary] light outside the 5" region (5"/15" counts < 0.70; trailed or smeared): 135, of which 109 event-mode snapshot openers
[summary] S/N < 10 at the source: 7 (faint, or its light is elsewhere: see the viewer)
[summary] on a LOW  low-sensitivity patch: 385 (B 4, U 54, UVM2 105, UVW1 99, UVW2 98, V 25)
[summary] on a MID  low-sensitivity patch: 652 (B 22, U 154, UVM2 152, UVW1 156, UVW2 125, V 43)
[summary] on a HIGH low-sensitivity patch: 786 (B 28, U 200, UVM2 165, UVW1 184, UVW2 144, V 65)
```

1,059 of the 2,879 annuli have at least one source masked out (up to 10).
All 12 offset circles are in pointings of SDSS J122933+015810, 8′ from
3C 273, where 3C 273 is 7–12′ off-axis and the annulus is not evenly exposed.

### The viewer

`swift_uvot_viewer.py` writes one PDF. The summary page shows the axis ratio
and concentration of event-mode snapshot openers against everything else,
both against time with the 2023–24 jitter period marked, the offsets with and
without aspect correction, and the low-sensitivity flags per filter. Contact
sheets follow, worst first: the most light outside the 5″ region, the
faintest at the source position, the most elongated, the largest offsets,
non-standard backgrounds, exposures on a LOW patch, and a random sample of
typical exposures for comparison. Each 90″ cutout shows the source region
(green), `--ra`/`--dec` (red cross) and the background region (cyan; masked
sources dashed), with S/N, axis ratio (`ar`), concentration (`c`) and offset
in the title.

## What to look for

- **Trailed event-mode openers.** 146 of 3C 273's 223 exposures with an axis
  ratio above 1.3 open a snapshot in event mode (Step 3), and so do 109 of
  the 135 with less than 0.70 of their 15″ counts in the 5″ region.
- **The 2023–24 spacecraft jitter** (August 2023 – April 2024). Its 122
  3C 273 exposures have a median axis ratio of 1.12, against 1.07 at other
  times. 35 are above 1.3 and 37 have a concentration below 0.70. Many show
  the source doubled, or smeared over up to 30″, which the axis ratio can
  miss: 29 of the 135 low-concentration exposures have an
  axis ratio of 1.3 or less. Ten ordinary exposures earlier in 2023 are
  elongated too; in other years 0–5 % are, and 11 % in 2005 (8 of 75).
- **Low S/N for a bright source** means its light is not in the region.
  3C 273 has seven such exposures, all from the jitter period. In the second
  snapshot of 00031659123 (December 2023, image mode, no aspect correction)
  the pointing wandered during the exposures. In U, V and the UV filters the
  sources are curved streaks 30–150″ long. In B they are compact, but 3C 273
  is not within 8″ of where the image puts it. In one 1,635 s UVM2
  event-mode exposure of 00089771001 the image shows only faint curved
  tracks a few hundred arcseconds long. The pointing moved during the
  exposure, and the archive builds an event-mode sky image with one pointing
  (Step 3). For a faint source the same S/N could just mean a
  non-detection; telling the two apart needs field stars (planned for faint
  targets).
- **Offsets over 3″** (111): 99 event-mode openers, 9 jitter-period
  exposures and 3 just over the limit (3.2–3.35″).
- **Low-sensitivity patches.** For 3C 273 the LOW maps flag 13 % of the
  exposures, MID 23 % and HIGH 27 %, most in U and the UV filters. The
  summary page shows the fraction per filter.
- **Non-standard backgrounds** (`bkg_type` not `annulus`): check them in the
  viewer.

## Inputs and outputs

**Inputs:** Step 3's `UVOT_output/uvot_inventory.txt` and the files it lists;
the UVOT CALDB; optionally `UVOT_output/uvot_region_overrides.txt`.

**Outputs** (in `UVOT_output`):

```
UVOT_output/
├── uvot_positions.txt       one row per candidate exposure (columns below)
├── uvot_detections.txt      uvotdetect sources within 150" of your source
├── uvot_positions.pdf       the viewer's contact sheets
├── detections/              uvotdetect output and logs, per observation and filter
└── regions/<OBSID>/
    ├── <EXTNAME>_src.reg    fk5;circle(ra,dec,5.0")
    └── <EXTNAME>_bkg.reg    fk5;annulus(ra,dec,27.5",35.0") and fk5;-circle(...) exclusions,
                             or fk5;circle(ra,dec,15.0")
```

| Column | What it is |
| ------ | ---------- |
| `snr` | Source S/N in the 5″ circle |
| `offset`, `d_ra`, `d_dec` | Centroid minus `--ra`/`--dec` (arcsec; `d_ra` along the sky) |
| `axis_ratio`, `sig_major`, `sig_minor`, `major_pa` | Shape: ratio of the axes, their sizes (arcsec), the major axis' position angle (deg east of north) |
| `conc` | Net counts in the 5″ source region / in a 15″ circle (point source ~0.85) |
| `center` | Where the source region is: `centroid`, `catalog`, or `override:<line>` |
| `src_ra`, `src_dec` | The source region's centre |
| `detx`, `dety` | The source's detector position |
| `sss_low`, `sss_mid`, `sss_high` | CALDB small-scale-sensitivity value there (`-99.9`: on a patch) |
| `bkg_type`, `bkg_area`, `n_excl` | `annulus`, `circle`, `override:<line>` or `none`; area (arcsec²); masked sources |
| `src_reg`, `bkg_reg` | The region files |
| `status`, `reason` | `ok`, `no_background` or `failed`; why the region is not at the centroid, and flags such as `trailed or smeared` |

### Overrides

To set a position or a background by hand, create
`UVOT_output/uvot_region_overrides.txt` and re-run the step:

```
# obsid        filter  extname        what  value
00035017124    V       vv387397671I   src   187.27791 2.05238
00035017045    *       *              bkg   circle(187.2900,2.0300,15")
00091742011    UVW2    *              src   catalog
```

`*` matches anything. `src` takes an RA and Dec in degrees, or `catalog`;
`bkg` takes an fk5 region in `uvotsource`'s syntax. The `center` and
`bkg_type` columns name the override line used. The automatic choice is
recomputed every run, so the file is the one place your decisions live.

## Common variants

```bash
# One observation, every exposure, with regions
swift_uvot_viewer.py --obsid 00031659123 --pdf obs123.pdf

# Fewer cutouts per category
swift_uvot_viewer.py --per-category 24

# The exposures flagged as trailed or smeared
grep "trailed or smeared" UVOT_output/uvot_positions.txt | less -S
```

## Gotchas

1. **Same `--ra`/`--dec` as Step 3.** The script refuses a different
   position.

2. **HEASoft terminal.** `uvotdetect` and the CALDB lookups need it; the
   viewer does not.

3. **Run from the analysis folder.** Step 3 recorded the image paths relative
   to it.

4. **A neighbour within 15″ lowers `conc`.** Its counts land in the 15″
   circle. Look at the cutout before you blame trailing.

5. **The readout streak of a bright source is not masked** in the background
   annulus, which it crosses in every exposure. A to-do.

6. **Re-running replaces the region files** for every candidate; keep your
   own choices in the overrides file.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_
