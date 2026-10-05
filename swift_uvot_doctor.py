#!/usr/bin/env python3
"""swift_uvot_doctor.py -- environment diagnostic for the Swift UVOT pipeline.

Prints a green/red checklist of everything the pipeline assumes about the
terminal (PATH, HEASoft, the UVOT CALDB, Python packages, disk) so a user
can verify readiness in one shot instead of watching the pipeline fail
midway. The UVOT pipeline needs one terminal: HEASoft set up, CIAO not set
up (see swift_uvot_env.py). The doctor ends by saying which steps this
terminal can run.

The CALDB checks ask the CALDB for every calibration file uvotsource uses,
through the pipeline's own runner (swift_uvot_runner.py), so they also
show that HEASoft tools run correctly from the pipeline. With --test-image
the doctor also runs uvotsource on one exposure of a UVOT sky image.

Exit status: 0 only if no check FAILs.

Usage:
    swift_uvot_doctor.py                    # full report (color on a terminal)
    swift_uvot_doctor.py --quiet            # only failing checks
    swift_uvot_doctor.py --no-color         # plain text, no ANSI
    swift_uvot_doctor.py --data-dir /media/drive2/meyer_swift_uvot
    swift_uvot_doctor.py --test-image sw00035017124uvv_sk.img.gz+2
"""

import argparse
import datetime
import importlib.metadata as _md
import importlib.util as _ilu
import os
import re
import shutil
import sys
import tempfile

import swift_uvot_env as uenv
from swift_uvot_runner import run_tool, WORK_ROOT_ENV

# Directory this script lives in; used to confirm PATH points at the pipeline.
PIPELINE_DIR = os.path.dirname(os.path.realpath(__file__))

# HEASoft tools the pipeline steps run (quzcif is used by the checks below;
# uvotsource runs uvotinteg, ftcreate and fthedit itself).
HEASOFT_TOOLS = ["uvotsource", "uvotinteg", "uvotdetect", "uvotimsum",
                 "quzcif", "ftlist", "fthedit"]

# Calibration files uvotsource reads at run time: (code name, filter,
# boundary expression, what it is). '-' means the file is not per filter.
CALDB_ITEMS = [
    ("COINCIDENCE", "-", "-", "coincidence-loss correction"),
    ("COLORTABLE", "-", "-", "zero points (Vega)"),
    ("ABCOLORTABLE", "-", "-", "zero points (AB)"),
    ("SKYFLAT", "V", "-", "large-scale sensitivity (LSS) map"),
    ("SKYFLAT_SSS", "V", "TYPE.eq.LOW", "small-scale sensitivity map, LOW"),
    ("SKYFLAT_SSS", "V", "TYPE.eq.MID", "small-scale sensitivity map, MID"),
    ("SKYFLAT_SSS", "V", "TYPE.eq.HIGH", "small-scale sensitivity map, HIGH"),
    ("SENSCORR", "V", "-", "sensitivity loss with time"),
    ("REEF", "V", "-", "encircled energy (aperture corrections)"),
]

# Oldest sensitivity-loss calibration the pipeline accepts: v007 (CALDB
# 20231208) refits all filters with data to 2023-08-01; v002 applied no
# correction at all.
MIN_SENSCORR_VERSION = 7

REQUIRED_PKGS = [
    ("astropy", "astropy"),
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("matplotlib", "matplotlib"),
    ("requests", "requests"),
]

OK, WARN, FAIL = "OK", "WARN", "FAIL"

_COLORS = {OK: "\033[32m", WARN: "\033[33m", FAIL: "\033[31m"}
_RESET = "\033[0m"

_results = []  # list of (status, message, [detail lines])
_use_color = False
_quiet = False
_resolved = {}  # (code name, boundary) -> CALDB file path, from check_uvot_caldb


def emit(status, message, details=None):
    """Record a check result and print it (unless suppressed by --quiet)."""
    _results.append((status, message, details or []))
    if _quiet and status != FAIL:
        return
    tag = "[%s]" % status
    if _use_color:
        tag = _COLORS[status] + tag + _RESET
    print("%s %s" % (tag, message))
    for line in (details or []):
        print("       %s" % line)


def _pkg_version(dist):
    try:
        return _md.version(dist)
    except Exception:
        return "?"


def _heasoft_ready():
    return bool(os.environ.get("HEADAS")) and not uenv.headas_is_ciao() \
        and all(shutil.which(t) for t in HEASOFT_TOOLS)


# --- individual checks ------------------------------------------------------

def check_path():
    resolved = shutil.which("swift_uvot_doctor.py")
    if resolved and os.path.realpath(resolved).startswith(PIPELINE_DIR):
        emit(OK, "Pipeline on PATH: swift_uvot_doctor.py -> %s" % resolved)
    elif resolved:
        emit(WARN, "swift_uvot_doctor.py resolves to %s, not %s "
             "(run setup_swiftuvot)" % (resolved, PIPELINE_DIR))
    else:
        emit(FAIL, "Pipeline not on PATH (no swift_uvot_doctor.py). "
             "Run setup_swiftuvot.")


def check_heasoft():
    headas = os.environ.get("HEADAS")
    if uenv.headas_is_ciao():
        emit(FAIL, "CIAO is set up in this terminal",
             ["CIAO's python replaced $HEADAS with %s," % headas,
              "so the HEASoft tools fail here. The UVOT pipeline needs no",
              "CIAO: open a new terminal and run setup_swiftuvot; heainit."])
        return
    if not headas:
        emit(FAIL, "HEASoft not loaded: $HEADAS unset. Run 'heainit'.")
        return
    missing, details = [], []
    for tool in HEASOFT_TOOLS:
        path = shutil.which(tool)
        if path:
            details.append("%-12s %s" % (tool, path))
        else:
            missing.append(tool)
            details.append("%-12s NOT FOUND" % tool)
    if missing:
        emit(FAIL, "HEASoft tools missing from PATH: %s (run 'heainit')"
             % ", ".join(missing), details)
    else:
        emit(OK, "HEASoft loaded ($HEADAS set, all UVOT tools on PATH)",
             details)


def check_heasoft_version():
    if uenv.headas_is_ciao() or not os.environ.get("HEADAS"):
        return  # check_heasoft already explained
    version = uenv.heasoft_version()
    if version is None:
        emit(WARN, "Could not parse HEASoft version (HEADAS=%s)"
             % (os.environ.get("HEADAS") or "unset"))
        return
    text = ".".join(str(v) for v in version)
    if version[:2] < uenv.MIN_HEASOFT:
        emit(FAIL, "HEASoft %s is older than %d.%d" % ((text,) + uenv.MIN_HEASOFT),
             ["uvotsource gained the small-scale-sensitivity check "
              "(SSS_FACTOR) in 6.30;", "the pipeline relies on it."])
    else:
        emit(OK, "HEASoft version %s" % text)


def check_caldb():
    caldb = os.environ.get("CALDB", "")
    cfg = os.environ.get("CALDBCONFIG", "")
    if not caldb:
        emit(FAIL, "$CALDB unset. Run 'heainit' (it sources caldbinit.sh).")
        return
    details = ["$CALDB=%s" % caldb]
    if not cfg:
        emit(FAIL, "$CALDBCONFIG unset (run 'heainit')", details)
        return
    details.append("$CALDBCONFIG=%s" % cfg)
    if not os.path.isfile(cfg):
        emit(FAIL, "$CALDBCONFIG=%s does not exist" % cfg, details)
        return
    with open(cfg) as f:
        lists_uvot = any(re.match(r"\s*SWIFT\s+UVOTA\s", line) for line in f)
    if not lists_uvot:
        emit(FAIL, "$CALDBCONFIG has no 'SWIFT UVOTA' line", details)
        return
    emit(OK, "CALDB configured", details)


def check_uvot_caldb():
    caldb = os.environ.get("CALDB", "")
    if not caldb:
        return  # check_caldb already reported it
    problem = uenv.caldb_problem(caldb)
    if problem:
        hint = []
        if uenv.caldb_has_uvot(uenv.DEFAULT_HEASOFT_CALDB) and \
                caldb != uenv.DEFAULT_HEASOFT_CALDB:
            hint = ["The UVOT CALDB is installed at %s; point $CALDB there "
                    "(heainit)." % uenv.DEFAULT_HEASOFT_CALDB]
        emit(FAIL, problem, hint)
        return
    index = uenv.uvot_caldb_index(caldb)
    if not _heasoft_ready():
        emit(WARN, "UVOT CALDB present (%s), but its files could not be "
             "looked up without HEASoft" % index)
        return
    details, failed = [], []
    for code, filt, expr, what in CALDB_ITEMS:
        res = run_tool("quzcif", ["SWIFT", "UVOTA", "-", filt, code,
                                  "now", "now", expr], timeout=60)
        found = res.stdout.split()[0] if res.ok and res.stdout.split() else ""
        if found and os.path.isfile(found):
            _resolved[(code, expr)] = found
            details.append("%-12s %s  (%s)" % (code, os.path.basename(found),
                                               what))
        else:
            failed.append(code)
            details.append("%-12s NOT FOUND (%s) %s" % (
                code, what, res.reason or res.tail.splitlines()[-1:]))
    if failed:
        emit(FAIL, "UVOT CALDB (%s) lacks files uvotsource needs: %s"
             % (index, ", ".join(sorted(set(failed)))), details)
    else:
        emit(OK, "UVOT CALDB %s: every file uvotsource uses resolves" % index,
             details)


def check_senscorr():
    path = _resolved.get(("SENSCORR", "-"))
    if not path:
        return  # check_uvot_caldb already reported
    m = re.search(r"v(\d+)\.fits", path)
    version = int(m.group(1)) if m else 0
    if _ilu.find_spec("astropy") is None:
        emit(WARN, "Sensitivity-loss file %s: astropy missing, coverage not "
             "checked" % os.path.basename(path))
        return
    from astropy.io import fits
    mjdref = 51910.0  # Swift MET zero point (2001-01-01)
    constant, extrapolated = [], []
    with fits.open(path) as hdul:
        for hdu in hdul[1:]:
            filt = hdu.header.get("FILTER", hdu.name)
            if filt == "MAGNIFIER":
                continue
            rows = sorted(zip(hdu.data["TIME"], hdu.data["SLOPE"]))
            if rows[-1][1] != 0:
                # the last row still has a slope: later data are extrapolated
                extrapolated.append("%s after %s" % (
                    filt, _met_date(rows[-1][0], mjdref)))
                continue
            # constant from the first row after the last one with a slope
            changing = [t for t, s in rows if s != 0]
            since = rows[0][0] if not changing else \
                min(t for t, s in rows if t > changing[-1])
            constant.append("%s %s" % (filt, _met_date(since, mjdref)[:4]))
    details = []
    if constant:
        details.append("correction held constant from: " + ", ".join(constant))
    if extrapolated:
        details.append("correction extrapolated with its last slope: "
                       + ", ".join(extrapolated))
    if version < MIN_SENSCORR_VERSION:
        emit(FAIL, "Sensitivity-loss calibration %s is older than v%03d: "
             "update the UVOT CALDB" % (os.path.basename(path),
                                        MIN_SENSCORR_VERSION), details)
    else:
        emit(OK, "Sensitivity-loss calibration %s" % os.path.basename(path),
             details)


def _met_date(met, mjdref):
    day = datetime.date(1858, 11, 17) + datetime.timedelta(
        days=mjdref + float(met) / 86400.0)
    return day.isoformat()


def check_python():
    v = sys.version_info
    msg = "Python %d.%d.%d (%s)" % (v.major, v.minor, v.micro, sys.executable)
    if (v.major, v.minor) < (3, 9):
        emit(FAIL, msg + " -- 3.9+ required")
    else:
        emit(OK, msg)


def check_required_pkgs():
    for imp, dist in REQUIRED_PKGS:
        if _ilu.find_spec(imp) is not None:
            emit(OK, "%-12s %s" % (imp, _pkg_version(dist)))
        elif imp == "requests":
            emit(WARN, "requests     not importable -- only Step 2 "
                 "(download) needs it")
        else:
            emit(FAIL, "%-12s NOT IMPORTABLE (pip install %s)" % (imp, dist))


def check_optional_pkgs():
    if _ilu.find_spec("astroquery") is not None:
        emit(OK, "astroquery   %s (optional)" % _pkg_version("astroquery"))
    else:
        emit(WARN, "astroquery not installed (optional; name resolution "
             "falls back to SIMBAD/NED/Sesame)")


def check_disk(data_dir):
    work_root = os.environ.get(WORK_ROOT_ENV) or tempfile.gettempdir()
    places = [("data", data_dir, 20), ("working folders", work_root, 2)]
    for label, path, min_gb in places:
        try:
            free_gb = shutil.disk_usage(path).free / (1024 ** 3)
        except OSError as exc:
            emit(WARN, "Could not check free space at %s: %s" % (path, exc))
            continue
        msg = "Disk free at %s (%s): %.1f GB" % (path, label, free_gb)
        if free_gb < min_gb:
            emit(WARN, msg + " (< %d GB)" % min_gb)
        else:
            emit(OK, msg)


def check_test_image(spec, ra=None, dec=None):
    """Run uvotsource on one exposure of a sky image, through the runner."""
    m = re.match(r"(.+?)(?:\+(\d+)|\[(\d+)\])?$", spec)
    path, ext = m.group(1), int(m.group(2) or m.group(3) or 1)
    if not os.path.isfile(path):
        emit(FAIL, "Test image %s not found" % path)
        return
    if not _heasoft_ready() or uenv.caldb_problem(os.environ.get("CALDB", "")):
        emit(FAIL, "Test image not measured: HEASoft or the UVOT CALDB is "
             "not ready (see above)")
        return
    from astropy.io import fits
    with fits.open(path) as hdul:
        hdr = hdul[ext].header
    ra = hdr.get("RA_OBJ") if ra is None else ra
    dec = hdr.get("DEC_OBJ") if dec is None else dec
    if ra is None or dec is None:
        emit(FAIL, "Test image %s[%d] has no RA_OBJ/DEC_OBJ; pass --ra and "
             "--dec" % (os.path.basename(path), ext))
        return
    work = tempfile.mkdtemp(prefix="swuvot_doctor_")
    with open(os.path.join(work, "src.reg"), "w") as f:
        f.write('fk5;circle(%.6f,%.6f,5")\n' % (ra, dec))
    with open(os.path.join(work, "bkg.reg"), "w") as f:
        f.write('fk5;annulus(%.6f,%.6f,27.5",35")\n' % (ra, dec))
    inputs = {"sky.img.gz": path, "src.reg": os.path.join(work, "src.reg"),
              "bkg.reg": os.path.join(work, "bkg.reg")}
    params = {"image": "sky.img.gz[%d]" % ext, "srcreg": "src.reg",
              "bkgreg": "bkg.reg", "sigma": 3, "apercorr": "NONE",
              "history": "no", "outfile": "phot.fits", "clobber": "yes",
              "chatter": 1}
    expmap = path.replace("_sk.img", "_ex.img")
    if expmap != path and os.path.isfile(expmap):
        inputs["ex.img.gz"] = expmap
        params["expfile"] = "ex.img.gz[%d]" % ext
    res = run_tool("uvotsource", params, inputs=inputs,
                   outputs={"phot.fits": os.path.join(work, "phot.fits")},
                   timeout=120)
    label = "uvotsource on %s[%d] (%s, %.1f s)" % (
        os.path.basename(path), ext, hdr.get("FILTER", "?"), res.elapsed)
    if not res.ok:
        emit(FAIL, label + ": " + res.reason, res.tail.splitlines())
        return
    row = fits.getdata(os.path.join(work, "phot.fits"), 1)[0]
    details = ["5\" circle at RA %.5f Dec %.5f; exposure %.1f s; frame time "
               "%.4f s" % (ra, dec, row["EXPOSURE"], row["FRAMTIME"]),
               "rate %.3f +/- %.3f ct/s; Vega mag %.3f +/- %.3f" % (
                   row["CORR_RATE"], row["CORR_RATE_ERR"], row["MAG"],
                   row["MAG_ERR"]),
               "SATURATED=%d  SSS_FACTOR=%.2f  coincidence factor %.3f" % (
                   row["SATURATED"], row["SSS_FACTOR"], row["COI_STD_FACTOR"])]
    if row["SATURATED"]:
        details.append("note: SATURATED -- over 0.98 raw counts per frame; "
                       "uvotsource caps the rate and sets a 100 % error")
    if row["SSS_FACTOR"] <= 0:
        details.append("note: the source sits on a low-sensitivity patch "
                       "(SSS_FACTOR %.1f); uvotsource sets MAG=99 and the "
                       "rate to -999" % row["SSS_FACTOR"])
    if "expfile" not in params:
        details.append("no exposure map found next to the image; field-of-"
                       "view check skipped")
    shutil.rmtree(work, ignore_errors=True)
    emit(OK, label, details)


def report_terminal_role():
    """Say which pipeline steps this terminal can run (not a check)."""
    def have(mod):
        return _ilu.find_spec(mod) is not None

    caldb_ok = uenv.caldb_problem(os.environ.get("CALDB", "")) is None
    rows = [
        ("Step 2", "download", have("requests"),
         "requests not importable here"),
        ("Steps 3-5", "inventory, positions, photometry",
         _heasoft_ready() and caldb_ok and have("astropy"),
         "needs HEASoft (no CIAO) and the UVOT CALDB: "
         "setup_swiftuvot; heainit"),
        ("Steps 6-7", "master table, light curve",
         all(have(m) for m in ("astropy", "numpy", "matplotlib")),
         "astropy/numpy/matplotlib missing"),
    ]
    print("\nThis terminal can run:")
    for step, what, ok, why in rows:
        print("  %-10s %-34s %s" % (step, what, "yes" if ok else "no -- " + why))


def main(argv=None):
    global _use_color, _quiet
    parser = argparse.ArgumentParser(
        description="Verify the Swift UVOT pipeline environment.")
    parser.add_argument("--quiet", action="store_true",
                        help="only print failing checks")
    parser.add_argument("--no-color", action="store_true",
                        help="disable ANSI color output")
    parser.add_argument("--data-dir", default=".",
                        help="where the UVOT data will live, for the "
                             "free-space check (default: current folder)")
    parser.add_argument("--test-image", metavar="FILE[+EXT]",
                        help="also run uvotsource on one exposure of this "
                             "UVOT sky image (*_sk.img.gz; default "
                             "extension 1)")
    parser.add_argument("--ra", type=float,
                        help="source RA for --test-image (default: RA_OBJ)")
    parser.add_argument("--dec", type=float,
                        help="source Dec for --test-image (default: DEC_OBJ)")
    args = parser.parse_args(argv)

    _quiet = args.quiet
    _use_color = sys.stdout.isatty() and not args.no_color \
        and os.environ.get("NO_COLOR") is None

    checks = [check_path, check_heasoft, check_heasoft_version, check_caldb,
              check_uvot_caldb, check_senscorr, check_python,
              check_required_pkgs, check_optional_pkgs,
              lambda: check_disk(args.data_dir)]
    if args.test_image:
        checks.append(lambda: check_test_image(args.test_image, args.ra,
                                               args.dec))
    for check in checks:
        try:
            check()
        except Exception as exc:  # a check must never crash the doctor
            emit(FAIL, "%s crashed: %s" % (getattr(check, "__name__", "check"),
                                           exc))

    report_terminal_role()

    n_ok = sum(1 for s, _, _ in _results if s == OK)
    n_warn = sum(1 for s, _, _ in _results if s == WARN)
    n_fail = sum(1 for s, _, _ in _results if s == FAIL)

    summary = "%d checks: %d ok, %d warn, %d fail" % (
        len(_results), n_ok, n_warn, n_fail)
    if _use_color:
        color = _COLORS[FAIL] if n_fail else (_COLORS[WARN] if n_warn
                                              else _COLORS[OK])
        summary = color + summary + _RESET
    print("\n" + summary)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
