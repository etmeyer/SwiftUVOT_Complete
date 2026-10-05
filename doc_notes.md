# Development notes (out-of-scope findings for later sessions)

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
