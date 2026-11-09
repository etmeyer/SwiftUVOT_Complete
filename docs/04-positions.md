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
3C 273 (2,891 candidate exposures) the positions take about 13 minutes on
4 cores, and the viewer about 20 seconds.

## What runs

```bash
swift_uvot_positions.py --ra 187.2779 --dec 2.0524
swift_uvot_viewer.py                    # -> UVOT_output/uvot_positions.pdf
```

| Flag | Purpose |
| ---- | ------- |
| `--ra` / `--dec` | Your source position; must match the one Step 3 used |
| `--outdir` | Step 3's output folder (default `UVOT_output`) |
| `--nproc` | Parallel workers (default: the CPUs the run may use, at most 8; pin a shell with `taskset -c 17-20 bash` and runs follow) |
| `--no-field-stars` | Skip the Gaia field stars: no pointing correction, no PSF columns |
| `--gaia-offline` | Use only the cached Gaia list (`uvot_gaia.txt`); stop if there is none for this position |

`swift_uvot_viewer.py` takes `--outdir`, `--pdf`, `--per-category` (cutouts
per category, default 48) and `--obsid` (every exposure of the given
observations instead of the categories).

## How it works

**Position and shape.** In each exposure the script finds the brightest pixel
within 8″ of `--ra`/`--dec`, centroids the source around it (4″), and
measures its S/N in the 5″ circle. At S/N ≥ 10 it also measures the shape
from the second moments within 6″: the axis ratio, the sizes along the major
and minor axes, and the major axis' position angle, and the offset of the
centroid from `--ra`/`--dec`. Below S/N 10 these columns are `-`: the
moments of a noise-dominated image are noise (36 exposures of the faint
blazar PKS 0537−286 read axis ratios above 1.5 from noise alone). For
3C 273 the median axis ratio is 1.07 (95 % of ordinary exposures are below
1.2), and 1.9 for the trailed event-mode snapshot openers of Step 3.

**Concentration.** The axis ratio sees only the central 6″. A trail longer
than that, or a smeared image with several lobes, can look round there and
still lose much of its light outside the 5″ source region. So the script also
divides the net counts in the 5″ region by the net counts in a 15″ circle
around it (column `conc`), with its error (`conc_err`: the counts, plus 1 %
of the background as a systematic). For a point source this is 0.86–0.88
(median of 3C 273's clean exposures per filter; 5 % are below 0.84).
Strong coincidence loss flattens the core and lowers it a little. Below
0.70 the reason column says `trailed or smeared`. `conc` is reported only
when `conc_err` is at most 0.02, otherwise it is `-` and the reason says
why: a faint source's concentration is too noisy to judge (at S/N 10 the
error is 0.08–0.2, and a true point source would read `smeared` 10–18 % of
the time). For 3C 273, 2,790 of 2,891 exposures are measured; for
PKS 0537−286 none.

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
calibration (Poole et al. 2008), with three kinds of cut-out:

- *Neighbours.* Every other source that `uvotdetect` finds (2.5σ) in the
  deepest exposure of that observation and filter is masked with a circle
  that reaches out to where its PSF wing falls below 0.5 % of the sky. The
  wing is the CALDB's PSF profile (the encircled-energy file used for
  aperture corrections, with an exponential tail beyond its 30″ end); the
  source's brightness is `uvotdetect`'s coincidence-corrected rate; the sky
  is measured 35–100″ around your source. Radii run from 5″ to 40″ (3C 273:
  median 23″). The old magnitude bins of HEASoft's `uvotlc` (5–20″) were
  far too small for faint targets: the wings of a V = 13.2 star 38″ from
  PKS 0537−286 raised its whole annulus by 2–4 %, enough to make the target
  read 8 % (V) to 40 % (U) too faint.
- *Bright sources.* A source at 0.9 counts per frame or more is saturated or
  nearly so, and its measured rate understates its wings: it gets the 40″
  maximum and a strip 2″ either side of its read-out direction (detector
  Y), along which it leaves a streak across the image. That includes your
  own source when it is that bright (3C 273: 1,081 annuli lose a strip).
  Faint detections within 45″ of a source 100 times brighter are almost
  always fragments of its halo, which `uvotdetect` breaks into many
  "sources" (3C 273 has over a thousand such detections): they keep the
  `uvotlc` radius, or they would eat the annulus.
- *Bright sectors.* The annulus is cut into 8 sectors of 45°. A sector whose
  counts are too high for the others (Poisson, p < 0.001) is dropped whole,
  and the test repeats, leaving at least 4. Single pixels are never
  clipped: at a few counts per pixel that biases the background low.

The background is then the plain mean of the pixels left — the estimator
`uvotsource` uses. Step 4 computes it itself for its own S/N and
concentration (the median it used before misses 8–15 % of the background in
V, B and U at a few counts per pixel, and nearly all of it in the UV, which
made faint sources look smeared). The column `bkg5` gives it per 5″ circle
(counts/s); it agrees with Step 5's `uvotsource` value to 0.1 % (median;
5–95 % within 1.7 %) for 3C 273. All cut-outs go into the region file
(`-circle` and `-polygon`), so Step 5 measures exactly the same pixels.

Your source itself is never masked (the brightest detection within 10″ of
`--ra`/`--dec`). `uvotsource` refuses a background that is not evenly
exposed, so the annulus must be fully exposed (within 1 %) and keep at least
half its area. Otherwise the script tries a 15″ circle 45″, 60″ or 90″ away,
wherever it is fully exposed and clear of masks and strips. If none
qualifies, the status is `no_background` (you can set one by hand, below).

**Low-sensitivity patches at the edge.** Besides the centre, the script
samples the whole 5″ circle every 0.5″ and records what fraction lies on a
LOW patch (`sss_low_frac`); Step 6 flags `sss_edge` above 10 %.

**Statuses.** Every candidate exposure gets one: `ok`, `no_background`, or
`failed` (with the reason). Nothing is excluded here; the flags and numbers
go to Step 6, which selects. Exposures that Step 3 did not pass keep their
Step 3 status.

For all of 3C 273:

```
[summary] 2891 candidate exposures: ok 2891
[summary] source region at: catalog 111, centroid 2780
[summary] background: annulus 2838, circle 53
          offsets, ASPCORR=DIRECT: median 0.25", 90% 0.45", 1 over 3"
          offsets, ASPCORR=NONE  : median 0.77", 90% 2.59", 103 over 3"
[summary] axis ratio > 1.3 (trailed or doubled?): 222, of which 146 event-mode snapshot openers
[summary] light outside the 5" region (5"/15" counts < 0.70; trailed or smeared): 128, of which 106 event-mode snapshot openers
[summary] background: 1081 with read-out strips, 984 with bright sectors dropped (1391 sectors)
[summary] concentration measured (error <= 0.02): 2790 of 2891
[summary] more than 10% of the 5" circle on a LOW patch: 557
[summary] S/N < 10 at the source: 7 (faint, or its light is elsewhere: see the viewer)
[summary] on a LOW  low-sensitivity patch: 385 (B 4, U 54, UVM2 105, UVW1 99, UVW2 98, V 25)
[summary] on a MID  low-sensitivity patch: 652 (B 22, U 154, UVM2 152, UVW1 156, UVW2 125, V 43)
[summary] on a HIGH low-sensitivity patch: 786 (B 28, U 200, UVM2 165, UVW1 184, UVW2 144, V 65)
```

The 53 offset circles are mostly UVM2 and UVW2 (42), in pointings of other
targets where 3C 273 is far off-axis (00089771, 00091742, 00050900) and the
annulus is not evenly exposed or its neighbours' masks take more than half
of it. Dropped sectors are common for 3C 273 itself (984 exposures, mostly
U and the UV) because its own halo is lopsided across its annulus; there
they lower the background by 11 % (median) and its rate by 0.02 %. For the
faint PKS 0537−286, 5 of 363 exposures drop a sector, 23 use an offset
circle, and 125 lose a read-out strip from one of the G ≈ 10 stars 3–5′
away.

### Field stars: pointing and PSF

A faint source cannot say where it is or how sharp its image is: its own
centroid and shape are noise. The stars around it can. Step 4 fetches the
Gaia DR3 sources within 9′ of your source once per analysis (ESA archive,
asynchronous TAP; G < 21) and keeps them in `UVOT_output/uvot_gaia.txt`,
with the service, query, date and position in its header. Later runs read
that file. When the archive is busy (HTTP 429 or 503) Step 4 waits
30, 60, then 120 s, or as long as the archive asks, and retries. With no
network and no cache Step 4 stops and says so; `--no-field-stars` runs
without them.

**Pointing.** Without aspect correction an exposure's coordinates can be
wrong, sometimes by much more than the 5″ aperture. In each exposure the
60 brightest peaks within 4′ of your source are paired with every Gaia
source (G < 19.5) within 5′ of them. The offset shared by the most pairs
wins: random pairs spread over the whole range, while the true offset gets
one vote per star. The stars' centroids near that offset then refine it,
and at least 3 must agree (to 0.7″). Columns `pnt_dra`, `pnt_ddec` (image
minus Gaia, arcsec), `pnt_n` (stars) and `pnt_err`. An exposure with too
few stars of its own (short UVM2) borrows the median offset of the other
filters of the same snapshot, which share the spacecraft attitude
(`pnt_src = snapshot`).

When the offset is over 1″, Step 4 works in the corrected coordinates and
writes the regions where your source really is in that image (`center =
stars`; the reason gives the offset). If your source is bright, its own
centroid still decides (`center = centroid`). If the corrected position
falls outside the exposed field, the status is `not_covered`: Step 3 judged
the coverage with the wrong coordinates.

For aspect-corrected PKS 0537−286 exposures the stars agree with Gaia to
0.1–0.3″, so nothing moves. Its four unaspected exposures of 2005-11-24
(00035240001, second snapshot: V, UVW1, UVM2, UVW2) are off by 3.4′ (+54″,
+198″, all four alike), so PKS 0537−286 was not in the exposed field at
all. Until PR1b they were measured at an empty patch of sky 3.4′ away, and
read zero.

**PSF.** Each isolated Gaia star (G < 19.5, no source brighter than
G + 2.5 within 15″, RUWE ≤ 1.4, 15–300″ from your source) that reaches S/N
20 at ≤ 0.3 counts per frame gets its concentration measured, the same 5″/15″
ratio as your source's, with its error. Divided by the same star's median in
normal exposures (aspect-corrected, outside the jitter period), it says
how much a broad or smeared PSF pushes out of the 5″ circle: `psf_loss` =
1 − the weighted mean of the stars' relative concentrations (a star more
than 4σ from the others is dropped), with `psf_loss_err`. The stars' axis
ratio and width are listed too (`star_ar`, `star_sigma`). An exposure with
no star of its own borrows the snapshot's other filters, with their spread
as the error (`psf_src = snapshot`). Every star measurement is in
`uvot_field_stars.txt`.

On stars of known flux (24 Gaia stars measured with `uvotsource` in every
PKS 0537−286 exposure), `psf_loss` predicts the light the 5″ circle really
loses:

| `psf_loss` | Exposures | Star rate / its true rate | 1 − `psf_loss` |
| ---------- | --------- | ------------------------- | -------------- |
| below 0.02 | 245 | 1.003 | 0.999 |
| 0.02–0.05 | 9 | 0.995 | 0.976 |
| 0.05–0.10 | 8 | 0.931 | 0.929 |
| 0.10–0.20 | 6 | 0.916 | 0.837 |
| 0.20–0.30 | 5 | 0.788 | 0.776 |

Above 5 % loss it is right on average to a few per cent but scatters by
about 0.4 × the loss per exposure. Step 6 uses it as a flag, Step 7
optionally as a correction.

3C 273 lies at high Galactic latitude, where there are few stars: the
pointing is measured in 1,003 of 2,891 exposures (1,561 more borrow it from
their snapshot). Aspect-corrected exposures agree with Gaia to 0.24″
(median; 90 % within 0.39″), unaspected ones to 0.23″ (largest 2.7″). The
correction (over 1″) applies to 15 exposures. In all of them 3C 273's own
centroid then places the region. One is 00035017085's U exposure, which an
uncorrected 3.35″ offset had left at the catalogue position with half the
light missed: centred, it reads 396 instead of 251 counts/s (it is
saturated, so Step 6 still excludes it).

For PKS 0537−286 the PSF is measured in 358 of 363 exposures (22 borrowed
from the snapshot). 35 lose more than 3 %: the 2022-06-02 flare
(00030816013: V 7 %, UVW1/UVW2 10 %, UVM2 6 %; U and B normal) and the
2023 jitter-period observations (8–40 %; the second snapshot of
00030816018 is normal, and 00030816019 loses 3–7 % in V, B and U).

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
sources, read-out strips and dropped sectors dashed), with S/N, axis ratio (`ar`), concentration (`c`) and offset
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
  non-detection; the field stars tell the two apart (`pnt_*` and
  `psf_loss`).
- **Offsets over 3″** (104): mostly event-mode openers and jitter-period
  exposures.
- **Low-sensitivity patches.** For 3C 273 the LOW maps flag 13 % of the
  exposures, MID 23 % and HIGH 27 %, most in U and the UV filters. The
  summary page shows the fraction per filter.
- **Non-standard backgrounds** (`bkg_type` not `annulus`): check them in the
  viewer.
- **A faint source** has `-` in the shape columns and `conc`, with the reason.
  That is expected: Step 6 can't judge its PSF from the source itself, and
  says so.

## Inputs and outputs

**Inputs:** Step 3's `UVOT_output/uvot_inventory.txt` and the files it lists;
the UVOT CALDB; optionally `UVOT_output/uvot_region_overrides.txt`.

**Outputs** (in `UVOT_output`):

```
UVOT_output/
├── uvot_positions.txt       one row per candidate exposure (columns below)
├── uvot_detections.txt      uvotdetect sources within 150" of your source, with
│                            rate, counts per frame and mask radius
├── uvot_gaia.txt            Gaia DR3 sources within 9' (the cache; delete it to query again)
├── uvot_field_stars.txt     every PSF star in every exposure: offset, S/N, counts per
│                            frame, axis ratio, width, concentration (and relative)
├── uvot_positions.pdf       the viewer's contact sheets
├── detections/              uvotdetect output and logs, per observation and filter
└── regions/<OBSID>/
    ├── <EXTNAME>_src.reg    fk5;circle(ra,dec,5.0")
    └── <EXTNAME>_bkg.reg    fk5;annulus(ra,dec,27.5",35.0") and fk5;-circle(...) masks,
                             fk5;-polygon(...) read-out strips and dropped sectors,
                             or fk5;circle(ra,dec,15.0")
```

| Column | What it is |
| ------ | ---------- |
| `snr` | Source S/N in the 5″ circle |
| `offset`, `d_ra`, `d_dec` | Centroid minus `--ra`/`--dec` (arcsec; `d_ra` along the sky); only at S/N ≥ 10 |
| `axis_ratio`, `sig_major`, `sig_minor`, `major_pa` | Shape: ratio of the axes, their sizes (arcsec), the major axis' position angle (deg east of north); only at S/N ≥ 10 |
| `conc`, `conc_err` | Net counts in the 5″ source region / in a 15″ circle (point source 0.86–0.88), and its error; only where the error is ≤ 0.02 |
| `center` | Where the source region is: `centroid`, `catalog`, `stars` (catalogue position moved by the field stars' pointing offset), or `override:<line>` |
| `src_ra`, `src_dec` | The source region's centre, in the image's own coordinates (what `uvotsource` reads) |
| `detx`, `dety` | The source's detector position |
| `sss_low`, `sss_mid`, `sss_high` | CALDB small-scale-sensitivity value there (`-99.9`: on a patch) |
| `sss_low_frac` | Fraction of the 5″ circle on a LOW patch |
| `bkg_type`, `bkg_area`, `n_excl`, `n_strip`, `n_sect` | `annulus`, `circle`, `override:<line>` or `none`; area (arcsec²); masked sources, read-out strips and dropped sectors in it |
| `bkg5` | Step 4's background in the 5″ circle (counts/s); compare with Step 5's |
| `sky` | The sky (counts/s/arcsec²) the mask radii were set against |
| `pnt_dra`, `pnt_ddec`, `pnt_n`, `pnt_err`, `pnt_src` | Pointing offset from the field stars (image − Gaia, arcsec), how many stars, its uncertainty, `stars` or `snapshot` |
| `n_stars`, `star_conc`, `star_conc_err`, `star_ar`, `star_sigma` | PSF stars used; their concentration relative to normal exposures, its error; their median axis ratio and width (arcsec) |
| `psf_loss`, `psf_loss_err`, `psf_src` | Light the 5″ circle loses to a broad PSF (1 − `star_conc`), its error, `stars` or `snapshot` |
| `src_reg`, `bkg_reg` | The region files |
| `status`, `reason` | `ok`, `no_background`, `not_covered` (only after a pointing correction) or `failed`; why the region is where it is, and flags such as `trailed or smeared` |

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
`bkg` takes an fk5 region in `uvotsource`'s syntax. Both are in the image's
own coordinates (what you read off that image in ds9), also for an exposure
whose pointing Step 4 corrected; `src catalog` means the catalogue position
in corrected coordinates. The `center` and
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

5. **The read-out direction is detector Y.** In the PKS 0537−286 field the
   clearest streaks (G ≈ 10 stars in long U exposures) lie within 2° of it;
   a few weaker UV cases were ambiguous. If a streak in the viewer runs
   outside its dashed strip, tell us.

6. **Re-running replaces the region files** for every candidate; keep your
   own choices in the overrides file.

7. **The Gaia list is cached.** Delete `UVOT_output/uvot_gaia.txt` to query
   again (it is also queried again if `--ra`/`--dec` change). The field
   needs a few stars: in a sparse field (3C 273, high Galactic latitude)
   short exposures often have none, and borrow from their snapshot or stay
   unchecked.

8. **The PSF stars must be point sources.** Gaia DR3 contains some compact
   galaxies; each star is compared only with itself, so a galaxy's own
   concentration does not matter, but a variable or blended one adds
   scatter. The weighted mean drops a star more than 4σ from the others.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_
