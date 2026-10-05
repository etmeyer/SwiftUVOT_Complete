"""
swift_uvot_env.py

Shared environment checks for the Swift UVOT pipeline scripts.

The UVOT pipeline runs in one terminal, with HEASoft set up and CIAO not
set up:

    setup_swiftuvot; heainit

It never needs CIAO (there is no spectral fitting), and CIAO must not be
set up in that terminal: once it is, `python3` resolves to CIAO's wrapper
script, which sets HEADAS to CIAO's bundled spectral directory and CALDB
to the Chandra CALDB inside every Python process. The HEASoft tools a
pipeline script then starts cannot find their own files or the Swift
calibration. (The XRT pipeline found this the hard way; see its Step 1.)

uvotsource reads the UVOT calibration at run time -- zero points,
coincidence loss, large- and small-scale sensitivity, sensitivity loss --
so the CALDB must contain the UVOT files. Its path must also be short:
uvotsource writes the CALDB file names into history keywords, and the
fthedit it calls aborts ("*** buffer overflow detected ***") once $CALDB
is 78 characters or longer (measured with HEASoft 6.36; /opt/CALDB is
10). uvotcoincidence, which uvotdetect runs, aborts the same way.

The pipeline scripts import this module from their own directory.
"""

import os
import re
import shutil
import subprocess
import sys

# Where the UVOT calibration lives inside a HEASoft CALDB.
UVOT_CALDB_SUBDIR = os.path.join('data', 'swift', 'uvota')
UVOT_CALDB_INDEX = os.path.join(UVOT_CALDB_SUBDIR, 'caldb.indx')

# Conventional HEASoft CALDB location, suggested in error messages.
DEFAULT_HEASOFT_CALDB = '/opt/CALDB'

# Longest $CALDB that works when a tool writes the CALDB file names into
# history keywords (uvotsource history=yes, uvotdetect's uvotcoincidence).
# 78 characters crashed fthedit in HEASoft 6.36.
CALDB_PATH_MAX = 77

# Oldest HEASoft with the small-scale-sensitivity check in uvotsource
# (sssfile, SSS_FACTOR), which the pipeline relies on.
MIN_HEASOFT = (6, 30)

HEASOFT_SHELL_HINT = (
    "Run this step from a terminal with HEASoft set up and CIAO NOT set up:\n"
    "    setup_swiftuvot; heainit\n"
    "(elsewhere: source headas-init.sh and caldbinit.sh, but not ciao.sh).\n"
    "Running heainit again in a terminal where `ciao` was already run does\n"
    "not undo CIAO -- open a new terminal.")


def _is_under(path, root):
    """True if path is root or lies inside it (symlinks resolved)."""
    if not path or not root:
        return False
    path = os.path.realpath(path)
    root = os.path.realpath(root)
    return path == root or path.startswith(root + os.sep)


def ciao_install():
    """CIAO install root if CIAO is set up in this environment, else ''."""
    return os.environ.get('ASCDS_INSTALL', '')


def headas_is_ciao():
    """True if $HEADAS, as this process sees it, points into CIAO."""
    return _is_under(os.environ.get('HEADAS', ''), ciao_install())


def caldb_has_uvot(caldb):
    """
    True/False for a local CALDB directory; None when it cannot be
    checked (unset, or a remote CALDB URL).
    """
    if not caldb or '://' in caldb:
        return None
    return os.path.exists(os.path.join(caldb, UVOT_CALDB_INDEX))


def caldb_problem(caldb):
    """Describe why caldb can't serve the UVOT pipeline, or None if it can."""
    if not caldb:
        return "CALDB is not set."
    if '://' in caldb:
        return ("CALDB=%s is a remote CALDB. It works, but every uvotsource "
                "call then downloads the 168 MB small-scale-sensitivity map "
                "(about 40 s per call). Install the UVOT CALDB locally."
                % caldb)
    if _is_under(caldb, ciao_install()) or 'ciao' in caldb.lower():
        return ("CALDB=%s is CIAO's Chandra calibration tree; it has no "
                "Swift UVOT files." % caldb)
    if not caldb_has_uvot(caldb):
        return ("No Swift UVOT calibration under CALDB=%s (missing %s)."
                % (caldb, UVOT_CALDB_INDEX))
    if len(caldb) > CALDB_PATH_MAX:
        return ("CALDB=%s is %d characters long; HEASoft's fthedit aborts "
                "(\"buffer overflow detected\") when the UVOT tools write "
                "CALDB file names longer than this allows. Use a path of at "
                "most %d characters, e.g. a symbolic link such as %s."
                % (caldb, len(caldb), CALDB_PATH_MAX, DEFAULT_HEASOFT_CALDB))
    return None


def uvot_caldb_index(caldb):
    """
    The UVOT CALDB index file that caldb.indx points to, e.g.
    'caldb.indx20240201', or None.
    """
    if not caldb_has_uvot(caldb):
        return None
    return os.path.basename(os.path.realpath(
        os.path.join(caldb, UVOT_CALDB_INDEX)))


def heasoft_version(headas=None):
    """HEASoft version as a tuple of ints, from the $HEADAS path, or None."""
    headas = os.environ.get('HEADAS', '') if headas is None else headas
    candidates = [headas]
    tool = shutil.which('uvotsource')
    if tool:
        candidates.append(os.path.realpath(tool))
    for src in candidates:
        m = re.search(r"heasoft[-_]?(\d+)\.(\d+)(?:\.(\d+))?", src, re.I)
        if m:
            return tuple(int(p) for p in m.groups() if p)
    return None


def require_heasoft_shell(tools):
    """
    Exit with an actionable message unless HEASoft (and not CIAO) is set
    up for this process and the CALDB can serve UVOT. Returns
    (headas, caldb) on success.
    """
    headas = os.environ.get('HEADAS', '')
    caldb = os.environ.get('CALDB', '')
    problems = []

    if headas_is_ciao():
        problems.append(
            "CIAO is set up in this terminal, and CIAO's python has "
            "replaced $HEADAS with %s, so the HEASoft tools this step "
            "runs cannot find their files." % headas)
    elif not headas:
        problems.append("HEADAS is not set (HEASoft is not initialized).")

    problem = caldb_problem(caldb)
    if problem:
        problems.append(problem)

    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        problems.append("Not on PATH: %s" % ', '.join(missing))

    if problems:
        print("ERROR: this step needs a HEASoft environment without CIAO.",
              file=sys.stderr)
        for problem in problems:
            print("  - " + problem, file=sys.stderr)
        print(HEASOFT_SHELL_HINT, file=sys.stderr)
        sys.exit(1)

    return headas, caldb


def pipeline_commit():
    """Short git commit of the pipeline scripts (+'-dirty'), or 'unknown'."""
    here = os.path.dirname(os.path.realpath(__file__))
    try:
        commit = subprocess.run(
            ['git', '-C', here, 'rev-parse', '--short', 'HEAD'],
            capture_output=True, text=True, timeout=10).stdout.strip()
        dirty = subprocess.run(
            ['git', '-C', here, 'status', '--porcelain', '--untracked-files=no'],
            capture_output=True, text=True, timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return 'unknown'
    if not commit:
        return 'unknown'
    return commit + ('-dirty' if dirty else '')


def provenance():
    """
    What produced a product: HEASoft version, CALDB and its UVOT index,
    and the pipeline commit. Later steps write this into their outputs,
    because uvotsource runs with history=no (see CALDB_PATH_MAX) and so
    does not record its CALDB files itself.
    """
    caldb = os.environ.get('CALDB', '')
    version = heasoft_version()
    return {
        'HEASOFT': '.'.join(str(v) for v in version) if version else 'unknown',
        'CALDB': caldb or 'unset',
        'UVOTINDX': uvot_caldb_index(caldb) or 'unknown',
        'PIPECOMM': pipeline_commit(),
    }
