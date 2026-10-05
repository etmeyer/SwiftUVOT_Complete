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
