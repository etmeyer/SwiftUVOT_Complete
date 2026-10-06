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
