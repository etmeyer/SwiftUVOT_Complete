# Development notes (out-of-scope findings for later sessions)

- **To do (idea, 2026-10-05): a "clean" option** that deletes the big,
  re-downloadable archive files of an analysis and keeps everything else.
  For 3C 273 the downloads are 4.9 GB (sky images 2.7, raw images 0.7,
  auxil 0.7, event files 0.5, exposure maps 0.2, housekeeping 0.1);
  everything the later steps write is small. Design notes:
  - Dry run by default: list what would go and how much space it frees;
    delete only with an explicit flag.
  - Write the exact command (and OBSID list) that brings the data back,
    e.g. `restore_data.sh`; the downloader resumes and fills in.
  - Keep a manifest of the archive files used (name, size, checksum,
    PROCVER) from the Step 3 inventory, so a later re-download can be
    compared: HEASARC sometimes reprocesses observations, and different
    files could change the results.
  - Also offer to remove regenerable intermediates (e.g. re-imaged event
    exposures, runner working folders); never touch user-edited files
    (overrides, notes) or the products.
  - After cleaning, Steps 6–7 (master table, light curve) still run, since
    they need only the photometry tables; Steps 3–5 need the data back.
  - Build it after Step 7, when the full set of products is known.

Logged during **Step 1** (setup, doctor, runner; branch `step1-setup`).

- **Things to fix in the XRT pipeline too** (found while planning to copy its
  downloader for Step 2; verified in `swift_xrt_download.py`, not yet fixed
  there):
  - `main()` always exits 0, even when files failed to download or the run was
    interrupted; the failure count is only printed.
  - `download_file()` writes straight to the destination, so an interrupted
    download leaves a truncated file. A re-run re-downloads it only if the
    server reports a size; when the HEAD request fails it keeps the file
    ("size unverifiable, skipping").
  - Only `requests.HTTPError` is caught per file; a dropped connection
    (`ConnectionError`, `ChunkedEncodingError`) ends the whole run.
  The UVOT copy (Step 2) will exit non-zero on any failure, download to a
  temporary name and rename, and treat every request error as a failed file.

- **HEASoft behaviour measured for the runner** (HEASoft 6.36, amorgos):
  - Without a controlling terminal and without `HEADASNOQUERY`/`HEADASPROMPT`,
    `uvotsource` prints `Task uvotsource 0.0 terminating with status 6` /
    `Unable to redirect prompts to the /dev/tty` and **exits 0** with no output
    file. `HEADASNOQUERY=1`, `HEADASPROMPT=/dev/null` or a pseudo-terminal each
    fix it; the runner sets both variables (no pty, so it is thread-safe) and
    also treats a "terminating with status N" line (N ≠ 0) as a failure. An
    empty `HEADASNOQUERY=` still counts as set.
  - `uvotsource history=yes` aborts in `fthedit` (`buffer overflow detected`,
    SIGABRT) when `$CALDB` is 78+ characters; 77 works. `uvotcoincidence`
    (inside `uvotdetect`, `uvotlc`) aborts the same way with a long `$CALDB`.
  - 384 parallel `uvotsource` calls (16 at a time) sharing one `PFILES` folder
    all succeeded with identical results; the folder stayed empty
    (`uvotsource` writes no `.par` files). Private `PFILES` are kept anyway
    for the other tools.
  - One `uvotsource` call: 0.5–0.9 s with the local CALDB, ~40 s with the
    remote CALDB (it downloads the 168 MB SSS map each time).
  - The runner's smoke tests (success, the silent `/dev/tty` failure, the
    long-`$CALDB` crash, timeout kills the process group, missing input,
    16-way == serial, and all of it under `setsid nohup`) are in
    `_dev_internal/step1_smoke/test_runner.py` (not published).

- **CALDB:** the UVOT tarball (`goodfiles_swift_uvota_20240201.tar.Z`) also
  contains `data/swift/mis` with index 20260901; amorgos has the 2023 `mis`
  (index 20230721). Only `data/swift/uvota` was extracted, so the XRT
  pipeline's environment is unchanged. Updating `mis` is a separate decision
  for the XRT side.

- **Shell setup:** `setup_swiftuvot` (docs/01-setup.md) has to be added to
  `/etc/bash.bashrc.local` by someone with root.

Logged during **Step 2** (download; branch `step2-download`).

- **Fixed in the UVOT copy of the downloader** (the three XRT issues above),
  plus two more that also affect `swift_xrt_download.py`:
  - The VO cone search returns only 11 columns unless asked for all with
    `VERB=3`, so `xrt_expo_pc` / `xrt_expo_wt` (and here `uvot_expo_*`) are
    never in the catalog rows; the XRT listing silently drops those columns.
  - The HEASARC year_month path was built from an MJD string as if it were a
    date (noted in the XRT doc_notes); the UVOT copy converts the MJD. The
    `heasarc_flat` candidate (`obs/<OBSID>/`) never exists and was dropped.
- With `astropy.io.votable`, HEASARC TAP results name the `obsid` column
  `DataLinkID` (its VOTable ID) unless `to_table(use_names_over_ids=True)`.
- UKSSDC's `archive/reproc/<OBSID>/uvot` files are byte-identical copies of
  HEASARC's (checked for 2013 and 2025 observations); the downloader uses
  HEASARC first and UKSSDC as fallback.
- Download tests (resume, truncated and corrupted files, unknown OBSID, 404,
  filter subset, Ctrl-C) are in `_dev_internal/step2_smoke/test_download.sh`;
  the 12-OBSID 3C 273 test set is `_dev_internal/step2_test_obsids.txt`.
  Note for such tests: non-interactive bash starts background jobs with
  SIGINT ignored; use `set -m` to send them Ctrl-C.
- Redirected to a file, Python block-buffers stdout, so a long download's log
  stopped at 35 of 368 observations while it was at ~70: it looked hung.
  The downloader now line-buffers stdout; later steps' scripts should too.
- A no-op re-run of all of 3C 273 took 4.5 min: one HEAD request per file
  (11,516), each on a new TLS connection. A `requests.Session` per thread
  (connection reuse) brought it to 1.7 min. The XRT downloader opens a new
  connection per request too. Decompressing every existing file to check it
  added ~5 min, so that is now `--verify`; new downloads are always checked.

Logged during **Step 3** (inventory; branch `step3-inventory`).

- **astropy and truncated .gz files:** `fits.open` on a truncated `.gz` sky
  image raised nothing and showed only the primary header, so an inventory
  that trusted `len(hdul)` would have dropped all its exposures without a
  word. The inventory now decompresses each `.gz` completely with `gzip`
  first (which raises on a truncation or a bad checksum) and parses the
  bytes; that was also twice as fast (37 s instead of 84 s for 3C 273). Later
  steps that open sky images with astropy should do the same (or run after
  the inventory has passed). CFITSIO-based HEASoft tools may behave
  differently; not checked.
- **The catalog's per-filter exposure is on-time** (before the dead-time
  correction): `ONTIME` sums match it, `EXPOSURE` sums are 1.6 % lower at
  11 ms frames and 4.8 % at 3.6 ms.
- **Archive gaps the file-level check can't see:** in 6 of 1,552 3C 273
  observation/filter pairs the catalog has more on-time than the sky images;
  for 00050900031 UVW1 HEASARC's exposure log (`swiftuvlog`) lists an
  exposure with extension `UNDEF` (1,082 s) that is in no sky image.
- **No image/event duplicates** (same `EXPID`) in any 3C 273 sky image, and no
  `OBS_MODE=SETTLING` extensions: the archive's sky images seem to hold one
  extension per exposure, pointing only. The checks stay (cheap) for other
  targets.
- **A 0.01 s exposure with a broken WCS** (00035017011, `bb156516976I`): the
  source maps to pixel (−20644, −15564). Classified `not_covered`, with the
  reason noting the tiny exposure.
- **Downloader, found while re-running Step 2 on all of 3C 273:** a no-op
  re-run that asked the server about every file (~100 requests/s with
  connection reuse) got HTTP 403 for 105 observations' folders. Fixed in
  this branch: files already on disk are kept without asking (this script
  only renames checked files into place; `--verify` still checks them),
  403/429/5xx are retried after 5, 15 and 45 s, and a folder HEASARC keeps
  refusing comes from the UKSSDC mirror. A no-op re-run now takes 56 s.
- Inventory tests (clean run with rows == extensions, truncated sky image and
  exposure map, missing exposure map, reasons read back) are in
  `_dev_internal/step3_smoke/test_inventory.sh`.
- **Which snapshot openers are trailed** (measured on all 2,886 3C 273
  candidates: second-moment axis ratio and centroid offset of the source):
  event-mode openers 82 % (axis ratio > 1.3 or offset > 2″; median ratio
  1.72, offset 2.9″, N=189); image-mode openers 5–7 %, the same as other
  exposures (2–6 %), consistent with UVOT's on-board drift correction
  ("shift-and-add") in image mode. So the inventory's summary counts the
  event-mode openers separately. Uncorrected aspect (`ASPCORR=NONE`) is
  almost entirely event mode (89 %) and UV images in the 5′ windows (84 %);
  full-frame image-mode exposures are 99 % corrected; most uncorrected
  exposures are still within ~1.3″ of the source (90 %).

Logged during **Step 4** (positions and regions; branch `step4-positions`).

- **Small-scale-sensitivity lookup ported from `uvotsource`** (UVOT::Source
  `updateDetectorPosition`/`applySmallScaleSensitivity`, UVOT::Calibration
  `estimateRAWfromDET`): DET from the sky image's `D` WCS (mm / 0.009075 +
  1100.5), the DET→RAW polynomial, clamp, subtract `UD_RAWX`/`UD_RAWY`
  ("best shift-and-add", set by `uvotimage` in newer processing), clamp. On
  1,457 3C 273 exposures DET agrees with `uvotsource` to < 0.001 pixel and
  the LOW flag agrees in all 1,457. `uvotsource` skips the `UD_RAW` shift if
  the environment variable `SSS_ADJUST_DISABLE` is set; the port always
  applies it. LOW is the default level (`SSS_TYPE || 'LOW'` in
  Calibration.pm).
- **`uvotlc`'s values, checked in HEASoft 6.36's UVOT/LCPar.pm:** background
  exclusion radius 20/15/10/7/5″ for magnitude ≤ 12/14/16/18/fainter, and
  `uvotdetect` threshold 2.5σ. (`uvotlc`'s own annulus is 12.5–25″; this
  pipeline uses Poole et al.'s 27.5–35″.) `uvotlc` also flags *source
  confusion* (a neighbour within 5″ plus a magnitude-dependent radius of
  the target); not done here yet — a to-do for crowded fields.
- **The target masked as its own neighbour:** in a trailed exposure
  (00035017175 UVW2) the target's `uvotdetect` position was more than 5″
  from the region centre, so its circle was cut out of the annulus. Now the
  brightest detection within 10″ of `--ra`/`--dec` is always the target.
- **Smeared images the axis ratio misses.** During the 2023–24 jitter,
  00089771001's 14 UVM2 event exposures (640–1,650 s, one pointing per sky
  image) show 3C 273 as a multi-lobed smear up to ~30″ across. In 13 of
  them the axis ratio within 6″ is 1.03–1.21 and S/N 130–220, but only
  33–60 % of the 15″ counts are in the 5″ circle, against ~86 % for a point
  source: photometry 1.4–2.6 times too faint, and nothing else flagged it.
  In the 14th only faint tracks remain (S/N 3.4). Hence
  the concentration column (5″/15″ net counts): ordinary exposures 0.84–0.87
  per filter; coincidence loss lowers it slightly (U above 0.9 counts/frame:
  median 0.83, 5 % below 0.77); flag below 0.70 (135 of 2,891: 109 event
  openers, 23 jitter-period image-mode exposures, 3 IMAGEEVENT exposures
  from 2005). 29 of the 135 have an axis ratio ≤ 1.3.
- **00031659123's second snapshot** (2023-12-05, image mode, ASPCORR=NONE):
  U, V, UVW1, UVW2, UVM2 sources are curved streaks 30–150″ long (on-board
  shift-and-add evidently could not follow the jitter); in B they are compact
  but 3C 273 is not within 8″ (S/N 1.2 where ~38 is expected). A
  peak-matching attempt to measure B's pointing error against the first
  snapshot did not converge (37 s exposure). All six have S/N < 10 at the
  source.
- **`uvotdetect` on a smeared exposure** detects lobes of the smeared target
  as separate sources, which are then masked in the annulus (00089771001
  UVM2, whose deepest aspect-corrected exposure is itself smeared). Harmless
  there, since those exposures are flagged; for other targets it would be
  better to detect on the deepest *compact* exposure (a second pass after
  the shapes are known). To do.
- **Found in Step 3 while checking these exposures** (fixed in this branch):
  - The inventory globbed `sw*_sk.img*`, which misses `sw…um2_sk_01.img.gz`:
    00089771001's event data were split into `_uf.evt.gz` and
    `_uf_01.evt.gz` (10.0 and 3.1 million events), each with its own sky
    image, exposure map and raw image. 5 UVM2 exposures (4,281 s) were
    silently missing from the ledger. The exposure-map name mapping
    (`_sk` → `_ex`) in the inventory, Step 4 and the doctor missed `_NN`
    too.
  - **Correction to the Step 3 note on archive gaps:** of the 6 catalog
    mismatches, one was the `_sk_01` file and one the `UNDEF` exposure; in
    the other 4, HEASARC's exposure log lists exactly the exposures in the
    sky images, and the master catalog's total is 9–11 % higher with nothing
    missing (00035017001 B: it counts the elapsed time of an `IMAGEEVENT`
    observation; the 3 others: unknown). So the inventory now checks the
    exposure log (`swiftuvlog`) by extension name instead: 2,935 of 2,936
    logged 3C 273 exposures have a sky image, the exception being the
    `UNDEF` one. Before the `_sk_01` fix this check would have reported the
    5 missing UVM2 exposures by name. The log lists event-mode exposures
    by their `…E` extension names too (from `_rw` images), plus rows for
    the event files without a name, which the check skips.
  - One header has `STALLOSS` = 3.99×10²⁵² (00089771001, `m2726326001E`);
    its `EXPOSURE` is unaffected (= `ONTIME` × the 11 ms dead-time factor).
    The inventory's `loss` column is now `-` with a note when the keywords
    are outside 0…(TSTOP−TSTART).
- **To do for faint targets:** a field-star check per exposure (positions
  and shapes of the brightest stars from the detection list): it tells a
  mis-pointed or smeared exposure from a genuine non-detection, which S/N at
  the target position cannot; it could also measure and correct the
  pointing of exposures without aspect correction.
- Step 4 tests (one row per candidate, trailed opener, jitter smear flagged
  by concentration, offset-circle background, SSS flags as `uvotsource`,
  region files accepted by `uvotsource`, overrides, position mismatch
  refused) are in `_dev_internal/step4_smoke/test_positions.sh`; Step 3's
  now include the `_sk_01` file and the exposure-log check.

Logged during **Step 5** (photometry; branch `step5-photometry`).

- **`uvotsource` 4.5 (HEASoft 6.36), read in UVOT/Source.pm and the task:**
  - Counts per frame for the coincidence correction =
    `RAW_STD_RATE` × (`FRAMTIME` − 6e-7 s × 290 rows), i.e. dead time
    0.174 ms per frame (1.6 % at 11 ms frames, 4.8 % at 3.6 ms). Above
    `FRAME_COUNT_LIMIT` = 0.98 the value is capped, `SATURATED` = 1 and
    the rate's error is set equal to the rate.
  - Count errors are binomial per frame (√(N(F−N)/F) for N counts in F
    frames, floor 1 count); the background's assumes an 80 arcsec²
    coincidence area.
  - It sets only these `PHOTFLAG` bits itself: `NO_EXPOSURE_MAP`,
    `NO_QUALITY_MAP` (always, without a quality map), `UNEVEN_EXPOSURE`,
    `BAD_SSS`, `BAD_LSS`. The other bits in UVOT::LCPar are `uvotlc`'s.
    `CORRFLAG` is always 0.
  - On a LOW patch it uses factor 1, then suppresses the magnitude:
    `CORR_RATE` and `SSS_RATE` −999, `MAG`/`AB_MAG` 99, but
    `SENSCORR_RATE` (the fully corrected rate) and `FLUX_*` stay filled.
    The pipeline's `rate` column is `SENSCORR_RATE` for that reason.
  - `ssstype` (LOW|MID|HIGH) is a hidden parameter, so the patch level can
    be changed in `uvotsource` itself; the pipeline records all three from
    Step 4 instead. Defaults worth overriding: `sigma` 5, `apercorr`
    CURVEOFGROWTH (no effect at 5″, but explicit is clearer).
  - `image=file.img.gz[EXTNAME]` works as well as `[N]`; the pipeline uses
    the extension name for the sky image and the exposure map, so the two
    files' extension order does not matter.
  - A refused region shows up as `error: ...` lines on stdout (e.g.
    `error: low bkg.reg in FOV 0.000`), then exit status 3 and
    "ERROR: No such process" on stderr; the reason column quotes the
    `error:` lines.
- **Run time:** 0.4 s per call alone; 2,891 calls 16 at a time took 4 m
  21 s (about 1.4 s each under load: every call decompresses its whole sky
  image). The output tables and logs take 126 MB for 3C 273.
- **First look at the flags' effect** (rate / median rate of clean exposures
  in the same observation and filter; numbers in docs/05-photometry.md):
  LOW patch 0.95 (29 % more than 10 % low), MID-only 0.99, HIGH-only 1.00
  — so exclude at LOW and check MID once more in Step 6; concentration
  < 0.70: 0.31; axis ratio > 1.3 with concentration ≥ 0.70: 0.92;
  saturated: 0.85; U at 0.95–0.98 counts/frame: 0.92. Near-saturation
  bins in B and the UV have too few same-observation references below
  0.85 counts/frame to judge; Step 6 should use the frame-time pairs
  (same snapshot, 11 ms and 3.6 ms) instead.
- **Tests no longer delete anything:** the smoke tests keep their work
  folders and print the command to remove them; deleting is left to the
  person running them. The Step 2 test truncates a file *inside the 3C 273
  test set* and lets the downloader repair it; ask before running it.
- Failed calls keep their working folders in `/tmp` (`swuvot_<tool>_*`) by
  design; they accumulate (tests provoke failures on purpose). To do:
  perhaps keep them under the analysis folder instead, where they are
  found and cleaned with the rest.
- Step 5 tests (one row per Step 4 row, LOW-patch and saturation handling,
  a refused background region, a carried `no_background`, serial ==
  parallel, `--plot-only` without HEASoft) are in
  `_dev_internal/step5_smoke/test_photometry.sh`.

Logged during **Step 6** (master table; branch `step6-master-table`).

- **Jitter dates, from the circulars:** GCN 34633 (2023-09-06) asks users to
  check UVOT data "starting ~ August 7" (gyro #1 noise, degraded attitude
  control). GCN 36033 (2024-04-05): two-gyro flight software, UVOT
  re-enabled on 4 April 2024; "the trailed point spread functions seen in
  some UVOT images during the last 6 months should be largely mitigated";
  and "the team is developing software updates to correct affected event
  mode data taken during this period" — worth checking before building
  event re-imaging (Step 3b). The rule uses 2023-08-07 to 2024-04-04.
- **Calibration** (3C 273; scripts in `_dev_internal/step6_calibration/`):
  same-snapshot 11 ms / 3.6 ms pairs (473) for coincidence loss; rate /
  median of clean same-observation exposures for shape, patches and
  centring. The numbers are in docs/06-master-table.md. Two surprises:
  UVW2 at 0.955–0.96 counts/frame (13 exposures, all 2005) reads 13–39 %
  high, and unsaturated full-frame U at 0.93–0.98 reads 8 % low against
  windowed U (which is itself corrected by a factor ~1.9, so the absolute
  U level is still open).
- **U in 2017–2022** is only full-frame, all at 0.95–0.98 counts/frame.
  First excluded by a 0.90 limit; at the user's request (2026-10-06) U is
  now kept up to saturation and flagged `high_coi` above 0.90 (about 8 %
  low). Readout-streak photometry (Page et al. 2013) is the clean fix; to
  do.
- **Significance needs both errors.** I first flagged same-snapshot pairs
  as 11–13σ discrepant using one exposure's error; with both errors
  combined they are 2.3–3.3σ (expected among ~470 pairs). The short
  full-frame member's binomial error may also be optimistic at high counts
  per frame. The `discrepant` flag uses the combined error and needs 10 %
  and 5σ; none for 3C 273.
- **Step 4 follow-up:** 3 non-trailed U exposures have centroids 3.2–3.35″
  from the source (pointing errors without aspect correction), so Step 4
  left their regions at the catalogue position and Step 6 excludes them
  (`off_centre`). Centring compact sources (concentration ≥ 0.78, axis
  ratio ≤ 1.5) up to ~5″ away would keep them. Minor; to do.
- Step 6 tests are in `_dev_internal/step6_smoke/test_master_table.sh`.
