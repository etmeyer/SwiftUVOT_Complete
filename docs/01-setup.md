# Step 1 — Setup

Summary: install the pipeline scripts and the Swift UVOT calibration database
(CALDB), set up one terminal with HEASoft, and check it with
`swift_uvot_doctor.py`. Unlike the XRT pipeline, the UVOT pipeline needs only
**one terminal**: there is no spectral fitting, so CIAO is never needed, and
it must not be set up in that terminal.

## What runs

### Requirements

**Operating environment:**
- Python 3.9+
- [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) 6.30 or newer (tested
  with 6.36). 6.30 added the small-scale-sensitivity check to `uvotsource`,
  which the pipeline relies on. The pipeline runs `uvotsource`, `uvotinteg`,
  `uvotdetect`, `uvotimsum`, `quzcif`, `ftlist` and `fthedit`.

**CALDB:**
- The HEASoft CALDB with the **Swift UVOT** files (`data/swift/uvota`, index
  20240201). `uvotsource` reads the zero points, coincidence-loss correction,
  large- and small-scale sensitivity maps and the sensitivity loss with time
  from it every time it runs. See [Install the UVOT CALDB](#install-the-uvot-caldb).

**Python packages:**
- `astropy`, `numpy`, `scipy`, `matplotlib`
- `requests` (download only)
- `astroquery` (optional, fallback name resolver)

On amorgos the HEASoft terminal's Python (the `heasoft` conda env) has all of
them except the optional `astroquery`.

### Install the pipeline

```bash
# Example: install to /opt/swift-uvot-pipeline
git clone https://github.com/etmeyer/SwiftUVOT_Complete.git /opt/swift-uvot-pipeline
chmod +x /opt/swift-uvot-pipeline/*.py
```

### Shell setup

On a shared machine, define the setup command in `/etc/bash.bashrc.local`
next to `setup_swiftxrt` and `heainit`, so every user gets it:

```bash
# Puts the UVOT pipeline scripts on PATH. Touches nothing else.
setup_swiftuvot() {
    local pipedir="/opt/swift-uvot-pipeline"
    if [ ! -d "$pipedir" ]; then
        echo "ERROR: $pipedir not found."
        return 1
    fi
    case ":$PATH:" in
        *":$pipedir:"*) ;;
        *) export PATH="$pipedir:$PATH" ;;
    esac
    if command -v swift_uvot_doctor.py &>/dev/null; then
        echo "Swift UVOT Pipeline ready ($(ls "$pipedir"/*.py 2>/dev/null | wc -l) scripts in $pipedir)"
    else
        echo "WARNING: scripts found in $pipedir but not executable."
        echo "  Run: chmod +x $pipedir/*.py"
        return 1
    fi
}
```

`heainit` (HEASoft plus the HEASoft CALDB) is the same command the XRT
pipeline uses; on amorgos it is already defined in that file.

### Install the UVOT CALDB

The UVOT files go into the same CALDB as the XRT ones. On amorgos that is
`/opt/CALDB` (owned by `meyer`, so no root is needed), whose `caldb.config`
already lists `SWIFT UVOTA`:

```bash
cd /opt/CALDB
curl -O https://heasarc.gsfc.nasa.gov/FTP/caldb/data/swift/uvota/goodfiles_swift_uvota_20240201.tar.Z
gzip -dc goodfiles_swift_uvota_20240201.tar.Z | tar -xf - data/swift/uvota
rm goodfiles_swift_uvota_20240201.tar.Z
```

The download is 385 MB and the files take 1.3 GB (the large- and
small-scale sensitivity maps are 168 MB each). The tarball also carries a
newer `data/swift/mis`; naming `data/swift/uvota` in the `tar` command
extracts only the UVOT part and leaves the `mis` files the XRT pipeline uses
as they are. Check the result with the doctor.

### One terminal

```bash
setup_swiftuvot
heainit
swift_uvot_doctor.py
```

**Never run `ciao` in this terminal.** CIAO's `python3` resets `$HEADAS` and
`$CALDB` to CIAO's own trees inside every pipeline script, so the HEASoft
tools fail and the Swift calibration is missing (the XRT pipeline's
[Step 1](https://etmeyer.github.io/SwiftXRT_Complete/01-setup/#two-terminals)
explains it). Running `heainit` again does not undo it; open a new
terminal. The scripts refuse to start in a CIAO terminal.

## How it works

### Verify with the doctor

`swift_uvot_doctor.py` checks everything the pipeline assumes about the
terminal and then lists which steps the terminal can run. It checks that:

- the pipeline is on `PATH` (`setup_swiftuvot`),
- HEASoft is loaded, CIAO is not, and the HEASoft tools resolve,
- HEASoft is 6.30 or newer,
- `$CALDB` and `$CALDBCONFIG` are set and the configuration lists UVOT,
- the CALDB returns every file `uvotsource` uses. The doctor asks the CALDB
  itself with `quzcif`, through the same runner the pipeline uses for all
  HEASoft calls (below), so this also shows that HEASoft tools run from the
  pipeline,
- `$CALDB` is short enough (see [Gotchas](#gotchas)),
- the sensitivity-loss calibration is v007 or newer, and from when on it
  holds each filter's correction constant,
- the Python packages are importable,
- there is free disk where the data will go (`--data-dir`) and where the
  runner makes its working folders.

With `--test-image` it also measures one exposure of a UVOT sky image with
`uvotsource` (a 5″ circle at the image's target position, or `--ra`/`--dec`,
and a 27.5–35″ background annulus).

A healthy terminal on amorgos, with a 3C 273 image from Step 2:

```
$ swift_uvot_doctor.py --data-dir /media/drive2/meyer_swift_uvot \
      --test-image sw00035017124uvv_sk.img.gz+2
[OK] Pipeline on PATH: swift_uvot_doctor.py -> /opt/swift-uvot-pipeline/swift_uvot_doctor.py
[OK] HEASoft loaded ($HEADAS set, all UVOT tools on PATH)
       uvotsource   /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/uvotsource
       ...
[OK] HEASoft version 6.36
[OK] CALDB configured
       $CALDB=/opt/CALDB
       $CALDBCONFIG=/opt/CALDB/software/tools/caldb.config
[OK] UVOT CALDB caldb.indx20240201: every file uvotsource uses resolves
       COINCIDENCE  swucountcor20010101v103.fits  (coincidence-loss correction)
       COLORTABLE   swuphot20041120v108.fits  (zero points (Vega))
       ABCOLORTABLE swuphot20041120v108.fits  (zero points (AB))
       SKYFLAT      swulssens20041120v003.fits  (large-scale sensitivity (LSS) map)
       SKYFLAT_SSS  swusslsens20041120v001.fits  (small-scale sensitivity map, LOW)
       SKYFLAT_SSS  swussmsens20041120v001.fits  (small-scale sensitivity map, MID)
       SKYFLAT_SSS  swusshsens20041120v001.fits  (small-scale sensitivity map, HIGH)
       SENSCORR     swusenscorr20041120v007.fits  (sensitivity loss with time)
       REEF         swureef20041120v104.fits  (encircled energy (aperture corrections))
[OK] Sensitivity-loss calibration swusenscorr20041120v007.fits
       correction held constant from: V 2022, B 2026, U 2022, UVW1 2024, UVM2 2024, UVW2 2024, WHITE 2022
[OK] Python 3.12.13 (/opt/anaconda3/envs/heasoft/bin/python3)
[OK] astropy      7.2.0
[OK] numpy        2.4.2
[OK] scipy        1.15.2
[OK] matplotlib   3.10.8
[OK] requests     2.34.2
[WARN] astroquery not installed (optional; name resolution falls back to SIMBAD/NED/Sesame)
[OK] Disk free at /media/drive2/meyer_swift_uvot (data): 661.9 GB
[OK] Disk free at /tmp (working folders): 13.9 GB
[OK] uvotsource on sw00035017124uvv_sk.img.gz[2] (V, 0.5 s)
       5" circle at RA 187.27750 Dec 2.05219; exposure 131.1 s; frame time 0.0036 s
       rate 112.275 +/- 1.096 ct/s; Vega mag 12.764 +/- 0.011
       SATURATED=0  SSS_FACTOR=1.00  coincidence factor 1.193

This terminal can run:
  Step 2     download                           yes
  Steps 3-5  inventory, positions, photometry   yes
  Steps 6-7  master table, light curve          yes

16 checks: 15 ok, 1 warn, 0 fail
```

"Held constant from" is when the calibration stops changing a filter's
correction: the fit used data up to 2023-08-01, so observations after those
years are corrected as if the detector stopped losing sensitivity then.

The same doctor in a terminal where `ciao` was run:

```
[FAIL] CIAO is set up in this terminal
       CIAO's python replaced $HEADAS with /opt/ciao/ciao-4.16/spectral,
       so the HEASoft tools fail here. The UVOT pipeline needs no
       CIAO: open a new terminal and run setup_swiftuvot; heainit.
[FAIL] CALDB=/opt/ciao/ciao-4.16/CALDB is CIAO's Chandra calibration tree; it has no Swift UVOT files.
       The UVOT CALDB is installed at /opt/CALDB; point $CALDB there (heainit).
```

The exit status is 0 only if no check fails, so the doctor can start a
script or cron job.

### How the pipeline runs HEASoft tools

Every HEASoft call in the pipeline goes through one function,
`run_tool` in `swift_uvot_runner.py` (`run_many` runs many calls in
parallel). Each call gets:

| What | Why |
| ---- | --- |
| Its own working folder | `uvotsource` writes temporary files into the current folder, and a crashed call leaves them behind. Separate folders keep parallel calls apart. |
| Short names for its inputs | Each input file is linked into the working folder under a short name. HEASoft tools have fixed-size path buffers; the XRT pipeline's `xselect` failed beyond ~90 characters. |
| A private parameter-file folder (`$PFILES`) | Parallel calls sharing `$HOME/pfiles` can collide (the XRT pipeline's parallel extraction lost observations that way). In 384 parallel `uvotsource` calls on one shared folder nothing failed, because `uvotsource` writes no parameter files, but other tools do. |
| `HEADASNOQUERY=1`, `HEADASPROMPT=/dev/null` | Without them, a HEASoft tool started without a terminal (nohup, cron, a script started by another program) stops with `Task uvotsource 0.0 terminating with status 6` / `Unable to redirect prompts to the /dev/tty` — and `uvotsource` still exits with status 0. |
| stdin from `/dev/null`, a timeout, its own process group | A stuck call is killed together with everything it started. |
| A log, if asked for | The command, working folder, exit status, time and complete output. |

A call counts as successful only if the tool exits 0, prints no
`terminating with status N` line with N other than 0, and leaves every
expected output file. Then its outputs are moved to their destination and the
working folder is deleted; a failed call keeps its folder for inspection.

**Timing on amorgos:** one `uvotsource` call takes 0.5–0.9 s with the local
CALDB; 64 calls on 16 workers take about 4 s. A remote CALDB makes each call
about 40 s (see [Gotchas](#gotchas)).

## Inputs and outputs

**Inputs:** none from a previous step: the HEASoft installation (`$HEADAS`),
the CALDB with the UVOT files, and the pipeline scripts.

**Outputs:** a terminal with the pipeline on `PATH`, HEASoft and the
HEASoft CALDB set up, no CIAO, and a doctor run with no `[FAIL]`.

## Common variants

```bash
swift_uvot_doctor.py                       # full colored checklist
swift_uvot_doctor.py --quiet               # only what is broken
swift_uvot_doctor.py --no-color            # plain text for logs
swift_uvot_doctor.py --data-dir /media/drive2/meyer_swift_uvot   # free space there
swift_uvot_doctor.py --test-image sw00035017124uvv_sk.img.gz+2   # also run uvotsource
swift_uvot_doctor.py --test-image FILE+N --ra 187.2779 --dec 2.0524

# Working folders somewhere other than the system temporary directory
export SWIFT_UVOT_TMP=/media/drive2/meyer_swift_uvot/tmp

# Single-user install: put setup_swiftuvot in ~/.bashrc instead
```

## Gotchas

1. **Never run `ciao` in the UVOT terminal.** See [One terminal](#one-terminal).
   The UVOT pipeline does not need CIAO at all.

2. **Keep `$CALDB` short: at most 77 characters.** `uvotsource` writes the
   names of its CALDB files into history keywords, and the `fthedit` it calls
   aborts (`*** buffer overflow detected ***`) once `$CALDB` is 78 characters
   or longer (measured with HEASoft 6.36). `uvotcoincidence`, which
   `uvotdetect` runs, aborts the same way with a long `$CALDB`. `/opt/CALDB`
   is 10. The pipeline also runs `uvotsource` with `history=no` and writes the
   CALDB version into its own products instead.

3. **Don't use a remote CALDB** (`CALDB=https://heasarc.gsfc.nasa.gov/FTP/caldb`)
   for real work. It works, but every `uvotsource` call downloads the 168 MB
   small-scale-sensitivity map, about 40 s per call.

4. **Running HEASoft tools yourself in a script or cron job:** set
   `HEADASNOQUERY=1` (or `HEADASPROMPT=/dev/null`), or a tool started without
   a terminal stops with `Unable to redirect prompts to the /dev/tty`. With
   `uvotsource` the exit status is still 0, so check that the output file
   exists.

5. **`uvotsource` appends to its output file** unless `clobber=yes`. The
   runner gives each call a fresh output file; by hand, re-running into the
   same file adds duplicate rows.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_
