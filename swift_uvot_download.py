#!/usr/bin/env python3
"""
swift_uvot_download.py -- download Swift UVOT data from the HEASARC archive.

Finds Swift observations of a source (by name, position, observation ID or
a file of IDs), lists them with their UVOT exposure in each filter, and
downloads the UVOT files the pipeline needs into one folder per OBSID:

    <outdir>/<OBSID>/uvot/image/   sky images (*_sk), exposure maps (*_ex),
                                   raw images (*_rw)
    <outdir>/<OBSID>/uvot/event/   event-mode files (to rebuild trailed
                                   exposures later)
    <outdir>/<OBSID>/uvot/hk/      UVOT housekeeping
    <outdir>/<OBSID>/auxil/        spacecraft attitude and orbit files

Adapted from the XRT pipeline's swift_xrt_download.py: the name resolution
and the catalog query are the same code. The differences:

  * Only the chosen UVOT filters are downloaded (default: V, B, U, UVW1,
    UVM2, UVW2; WHITE, the grisms and the magnifier only on request).
  * Each file is downloaded under a temporary name and renamed when it is
    complete, so an interrupted run never leaves a truncated file under the
    real name. Every .gz file is decompressed once to check it. Network
    errors are retried, then counted as failed files instead of ending the
    run.
  * The HEASARC archive is reached through its year_month folders (from the
    catalog's start_time), with the identical UKSSDC mirror as fallback.
  * Afterwards, each OBSID's sky images are checked against the catalog's
    UVOT exposure in each filter, and download_report.txt in the output
    folder records the outcome for every OBSID.
  * The exit status is 1 if any file failed to download or an OBSID could
    not be found in the archive (re-running the same command retries them).

Requirements:
    requests; astropy (catalog lookup for --obsid / --obsid-file)
    astroquery (optional, fallback name resolver)

Usage:
    # List the observations of a source with their UVOT exposures
    swift_uvot_download.py --name "3C 273" --list-only

    # Download a date window
    swift_uvot_download.py --name "3C 273" --start-date 2009-01-01 \
        --end-date 2010-01-01 --outdir UVOT_input

    # By coordinates, or by observation ID(s)
    swift_uvot_download.py --ra 187.2779 --dec 2.0524 --outdir UVOT_input
    swift_uvot_download.py --obsid 00035017124 --outdir UVOT_input
    swift_uvot_download.py --obsid-file obsids.txt --outdir UVOT_input

    # Only some filters; four observations at a time
    swift_uvot_download.py --obsid-file obsids.txt --filters w2 m2 w1 --nproc 4
"""

from __future__ import annotations

import argparse
import gzip
import io
import math
import os
import re
import sys
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Name resolution (from swift_xrt_download.py)
# ---------------------------------------------------------------------------

def resolve_name_simbad(name: str) -> tuple[float, float]:
    """Resolve via the SIMBAD TAP service (robust, handles most names)."""
    url = "https://simbad.u-strasbg.fr/simbad/sim-id"
    params = {
        "Ident": name,
        "output.format": "votable",
        "output.params": "main_id,ra(d),dec(d)",
    }
    print(f"[info] Resolving '{name}' via SIMBAD …")
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    # Parse the VOTable response
    root = ET.fromstring(resp.text)
    ns_candidates = [
        "http://www.ivoa.net/xml/VOTable/v1.3",
        "http://www.ivoa.net/xml/VOTable/v1.2",
        "http://www.ivoa.net/xml/VOTable/v1.1",
        "",
    ]
    for ns_uri in ns_candidates:
        ns = {"v": ns_uri} if ns_uri else {}
        prefix = "v:" if ns_uri else ""
        fields = root.findall(f".//{prefix}FIELD", ns)
        if fields:
            break
    else:
        raise RuntimeError(f"SIMBAD returned no parseable VOTable for '{name}'")

    col_names = [f.attrib.get("name", f.attrib.get("ID", "")).lower()
                 for f in fields]
    tds = root.findall(f".//{prefix}TR/{prefix}TD", ns)
    if not tds or len(tds) < len(col_names):
        raise RuntimeError(f"SIMBAD returned no results for '{name}'")

    values = {c: td.text.strip() if td.text else ""
              for c, td in zip(col_names, tds)}

    # SIMBAD uses column names like "ra_d" / "dec_d" or "ra(d)" / "dec(d)"
    ra_str = values.get("ra_d") or values.get("ra(d)", "")
    dec_str = values.get("dec_d") or values.get("dec(d)", "")
    if not ra_str or not dec_str:
        raise RuntimeError(f"SIMBAD result for '{name}' missing coordinates "
                           f"(columns: {list(values.keys())})")

    ra, dec = float(ra_str), float(dec_str)
    print(f"[info] SIMBAD resolved to RA={ra:.5f}°, Dec={dec:.5f}°")
    return ra, dec


def resolve_name_ned(name: str) -> tuple[float, float]:
    """Resolve via the NASA/IPAC Extragalactic Database (NED)."""
    url = "https://ned.ipac.caltech.edu/srs/ObjectLookup"
    params = {"name": name, "of": "json"}
    print(f"[info] Resolving '{name}' via NED …")
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    # NED returns a "ResultCode" of 3 for success, with position in "Preferred"
    result_code = data.get("ResultCode", 0)
    if result_code != 3:
        raise RuntimeError(f"NED could not resolve '{name}' "
                           f"(ResultCode={result_code})")

    pos = data.get("Preferred", {}).get("Position", {})
    ra = pos.get("RA")
    dec = pos.get("Dec")
    if ra is None or dec is None:
        raise RuntimeError(f"NED result for '{name}' missing coordinates")

    ra, dec = float(ra), float(dec)
    print(f"[info] NED resolved to RA={ra:.5f}°, Dec={dec:.5f}°")
    return ra, dec


def resolve_name_sesame(name: str) -> tuple[float, float]:
    """Fallback resolver using CDS Sesame (queries SIMBAD, NED, VizieR)."""
    url = "https://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-ox/SNV"
    params = {"obj": name}
    print(f"[info] Resolving '{name}' via CDS Sesame …")
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    root = ET.fromstring(resp.text)
    ns = {"s": "http://vizier.u-strasbg.fr/xml/sesame_xml.xsd"}
    resolver = root.find(".//s:Resolver", ns)
    if resolver is None:
        raise RuntimeError(f"Sesame returned no result for '{name}'")

    jradeg = resolver.find("s:jradeg", ns)
    jdedeg = resolver.find("s:jdedeg", ns)
    if jradeg is None or jdedeg is None:
        raise RuntimeError(f"Sesame result for '{name}' missing coordinates")

    ra, dec = float(jradeg.text), float(jdedeg.text)
    print(f"[info] Sesame resolved to RA={ra:.5f}°, Dec={dec:.5f}°")
    return ra, dec


def resolve_name_astroquery(name: str) -> tuple[float, float]:
    """Last-resort resolver using astroquery (if installed)."""
    try:
        from astroquery.simbad import Simbad
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        print(f"[info] Resolving '{name}' via astroquery/SIMBAD …")
        result = Simbad.query_object(name)
        if result is None:
            raise RuntimeError(f"Could not resolve source name '{name}'")

        # astroquery ≥0.4.8 uses lowercase "ra"/"dec" (degrees);
        # older versions use uppercase "RA"/"DEC" (sexagesimal).
        colnames = [c.lower() for c in result.colnames]
        if "ra" in colnames and "dec" in colnames:
            ra_col = result.colnames[colnames.index("ra")]
            dec_col = result.colnames[colnames.index("dec")]
            ra_val, dec_val = result[ra_col][0], result[dec_col][0]
            # New astroquery returns degrees directly; old returns strings
            try:
                ra, dec = float(ra_val), float(dec_val)
            except (TypeError, ValueError):
                coord = SkyCoord(ra_val, dec_val, unit=(u.hourangle, u.deg))
                ra, dec = coord.ra.deg, coord.dec.deg
        else:
            raise RuntimeError(
                f"Unexpected columns in astroquery result: {result.colnames}")

        print(f"[info] astroquery resolved to RA={ra:.5f}°, Dec={dec:.5f}°")
        return ra, dec
    except ImportError:
        raise RuntimeError("astroquery is not installed")


# Ordered chain of resolvers – tried in sequence until one succeeds
_RESOLVERS = [
    ("SIMBAD", resolve_name_simbad),
    ("NED", resolve_name_ned),
    ("Sesame", resolve_name_sesame),
    ("astroquery", resolve_name_astroquery),
]


def resolve_name(name: str) -> tuple[float, float]:
    """Try each resolver in turn; return the first successful result."""
    errors: list[str] = []
    for label, func in _RESOLVERS:
        try:
            return func(name)
        except Exception as exc:
            print(f"[warn] {label} resolution failed: {exc}")
            errors.append(f"{label}: {exc}")
    raise RuntimeError(
        f"All resolvers failed for '{name}':\n  " + "\n  ".join(errors)
    )


# ---------------------------------------------------------------------------
# HEASARC query
# ---------------------------------------------------------------------------

HEASARC_CONESEARCH_URL = "https://heasarc.gsfc.nasa.gov/xamin/vo/cone"
HEASARC_TAP_URL = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
HEASARC_BROWSE_URL = "https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3query.pl"


def query_swift_master(ra: float, dec: float, radius_arcmin: float = 12.0) -> list[dict]:
    """
    Query the HEASARC swiftmastr catalog for Swift observations near (ra, dec).

    Tries multiple strategies in order:
      1. VO Cone Search (simple, standard, positional by design)
      2. W3Browse BatchDisplay text output
      3. HEASARC TAP with coordinate box

    Parameters
    ----------
    ra, dec : float
        Position in decimal degrees (J2000).
    radius_arcmin : float
        Search cone radius in arcminutes (default 12′, roughly the XRT FOV).

    Returns
    -------
    list of dict
        Each dict contains observation metadata fields.
    """
    strategies = [
        ("VO Cone Search",  lambda: _query_via_conesearch(ra, dec, radius_arcmin)),
        ("W3Browse (text)", lambda: _query_via_w3browse_text(ra, dec, radius_arcmin)),
        ("TAP (box)",       lambda: _query_via_tap_box(ra, dec, radius_arcmin)),
    ]
    for label, func in strategies:
        try:
            result = func()
            if result:
                return result
            print(f"[warn] {label} returned 0 results, trying next strategy …")
        except Exception as exc:
            print(f"[warn] {label} failed ({exc}), trying next strategy …")

    # All strategies exhausted
    return []


# --- Strategy 1: VO Cone Search (standard IVOA protocol) -------------------

def _query_via_conesearch(ra: float, dec: float,
                           radius_arcmin: float) -> list[dict]:
    """Query using the IVOA Simple Cone Search protocol.

    This is the most straightforward positional query — RA, Dec, and
    search radius are the only parameters.  Returns VOTable XML.
    """
    radius_deg = radius_arcmin / 60.0
    params = {
        "table": "swiftmastr",
        "RA": str(ra),
        "DEC": str(dec),
        "SR": str(radius_deg),
        # All columns: the default leaves out the per-filter exposures
        # (uvot_expo_vv, ...), which the listing and the check use.
        "VERB": "3",
    }
    print(f"[info] Querying HEASARC VO Cone Search (swiftmastr) within "
          f"{radius_arcmin}′ of RA={ra:.5f}, Dec={dec:.5f} …")
    resp = requests.get(HEASARC_CONESEARCH_URL, params=params, timeout=90)
    resp.raise_for_status()
    return _parse_votable(resp.text)


# --- Strategy 2: W3Browse text output (pipe / plus-delimited) ---------------

def _query_via_w3browse_text(ra: float, dec: float,
                              radius_arcmin: float) -> list[dict]:
    """Query using W3Browse with pipe-delimited text output."""
    # Format coordinates explicitly as decimal degrees for W3Browse.
    # The sign on Dec must be explicit; append 'd' to each value.
    dec_str = f"+{dec}" if dec >= 0 else f"{dec}"
    params = {
        "tablehead": "name=heasarc_swiftmastr&description=Swift Master Catalog",
        "Action": "Query",
        "Coordinates": f"{ra}d {dec_str}d",
        "Equinox": "2000",
        "Radius": str(radius_arcmin),
        "Radius_unit": "arcmin",
        "NR": "CheckCaches/GRB/SIMBAD/NED",
        "ResultMax": "10000",
        "displaymode": "BatchDisplay",
        "Fields": "All",
        "vession": "img",
        "gifsize": "0",
    }
    print(f"[info] Querying HEASARC W3Browse (text) within {radius_arcmin}′ of "
          f"RA={ra:.5f}, Dec={dec:.5f} …")
    resp = requests.get(HEASARC_BROWSE_URL, params=params, timeout=60)
    resp.raise_for_status()
    return _parse_batch_text(resp.text)


def _parse_batch_text(text: str) -> list[dict]:
    """Parse the pipe-delimited BatchDisplay output from W3Browse.

    Handles both separator styles::

        col1|col2|col3          (pipes only)
        ---|---|---

        |col1|col2|col3|        (pipes with leading/trailing)
        +----+----+----+        (plus-delimited separator)
    """
    lines = text.strip().splitlines()

    # Find the separator row.  It looks like one of:
    #   ---|---|---        (dashes and pipes)
    #   +------+------+   (dashes and plus signs)
    header_idx = None
    data_start = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        # A separator line is composed entirely of '-', '|', '+', and spaces
        if (stripped
                and all(c in "-|+ " for c in stripped)
                and "---" in stripped
                and i > 0):
            header_idx = i - 1
            data_start = i + 1
            break

    if header_idx is None or header_idx < 0:
        preview = "\n".join(lines[:10])
        raise RuntimeError(
            f"Could not find table header in W3Browse output:\n{preview}")

    col_names = [c.strip() for c in lines[header_idx].split("|") if c.strip()]

    observations = []
    for line in lines[data_start:]:
        line = line.strip()
        if not line or line.startswith("<") or line.startswith("Search"):
            continue
        # Stop at any trailing footer
        if line.startswith("BatchEnd") or line.startswith("***"):
            break
        values = [v.strip() for v in line.split("|")]
        # Strip empty tokens from leading/trailing pipes
        if values and values[0] == "":
            values = values[1:]
        if values and values[-1] == "":
            values = values[:-1]
        if len(values) < len(col_names):
            continue
        row = {c: values[j] for j, c in enumerate(col_names)}
        observations.append(row)

    print(f"[info] Found {len(observations)} observation(s).")
    return observations


# --- Strategy 3: TAP with coordinate box ------------------------------------

def _query_via_tap_box(ra: float, dec: float,
                        radius_arcmin: float) -> list[dict]:
    """Query HEASARC TAP using a simple coordinate-range WHERE clause."""
    radius_deg = radius_arcmin / 60.0
    cos_dec = math.cos(math.radians(dec))
    ra_range = radius_deg / max(cos_dec, 0.01)

    adql = (
        f"SELECT * FROM swiftmastr "
        f"WHERE ra BETWEEN {ra - ra_range} AND {ra + ra_range} "
        f"AND dec BETWEEN {dec - radius_deg} AND {dec + radius_deg}"
    )
    params = {
        "request": "doQuery",
        "version": "1.0",
        "lang": "ADQL",
        "format": "votable",
        "query": adql,
    }
    print(f"[info] Querying HEASARC TAP (box) (swiftmastr) within "
          f"{radius_arcmin}′ of RA={ra:.5f}, Dec={dec:.5f} …")
    resp = requests.get(HEASARC_TAP_URL, params=params, timeout=60)
    resp.raise_for_status()
    return _parse_votable(resp.text)


def _parse_votable(xml_text: str) -> list[dict]:
    """Minimal VOTable parser – extracts TABLEDATA rows."""
    # Guard against non-XML responses (HTML error pages, plain text, etc.)
    stripped = xml_text.strip()
    if not stripped.startswith("<?xml") and not stripped.startswith("<VOTABLE") \
       and not stripped.startswith("<vo:VOTABLE"):
        # Show a useful snippet for debugging
        preview = stripped[:300].replace("\n", " ")
        raise RuntimeError(
            f"Expected VOTable XML but got unexpected response: {preview!r}…"
        )

    root = ET.fromstring(xml_text)

    # Check for an INFO element with an error
    for ns_uri in ["http://www.ivoa.net/xml/VOTable/v1.3",
                    "http://www.ivoa.net/xml/VOTable/v1.2",
                    "http://www.ivoa.net/xml/VOTable/v1.1", ""]:
        ns = {"v": ns_uri} if ns_uri else {}
        prefix = "v:" if ns_uri else ""
        for info in root.findall(f".//{prefix}INFO", ns):
            if info.attrib.get("name") == "QUERY_STATUS" \
               and info.attrib.get("value") == "ERROR":
                msg = info.text.strip() if info.text else "unknown error"
                raise RuntimeError(f"HEASARC query error: {msg}")

    # VOTable namespace
    ns_candidates = [
        "http://www.ivoa.net/xml/VOTable/v1.3",
        "http://www.ivoa.net/xml/VOTable/v1.2",
        "http://www.ivoa.net/xml/VOTable/v1.1",
        "",
    ]
    for ns_uri in ns_candidates:
        ns = {"v": ns_uri} if ns_uri else {}
        prefix = "v:" if ns_uri else ""
        fields = root.findall(f".//{prefix}FIELD", ns)
        if fields:
            break
    else:
        print("[warn] No FIELD elements found in VOTable – "
              "the query may have returned no results.")
        return []

    col_names = [f.attrib.get("name", f"col{i}") for i, f in enumerate(fields)]
    rows_el = root.findall(f".//{prefix}TR", ns)
    observations = []
    for tr in rows_el:
        tds = tr.findall(f"{prefix}TD", ns)
        row = {}
        for cname, td in zip(col_names, tds):
            row[cname] = td.text.strip() if td.text else ""
        observations.append(row)

    print(f"[info] Found {len(observations)} observation(s).")
    return observations



# ---------------------------------------------------------------------------
# Catalog lookup by OBSID (for --obsid / --obsid-file)
# ---------------------------------------------------------------------------

def _cell(value) -> str:
    """One table cell as the string the cone-search parser would give."""
    if value is None or getattr(value, "mask", False) is True:
        return ""
    if isinstance(value, bytes):
        value = value.decode()
    return str(value).strip()


def query_swift_master_obsids(obsids: list[str]) -> dict[str, dict]:
    """
    swiftmastr rows for the given OBSIDs, from HEASARC TAP, keyed by OBSID.

    Returns what it could get (possibly {}) with a warning on failure; the
    download then falls back to the UKSSDC mirror and skips the
    per-filter check for OBSIDs without a row.
    """
    try:
        from astropy.io.votable import parse_single_table
    except ImportError:
        print("[warn] astropy is not importable here, so the catalog can't "
              "be asked about these OBSIDs: no listing, no per-filter check.")
        return {}
    rows: dict[str, dict] = {}
    for i in range(0, len(obsids), 100):
        chunk = obsids[i:i + 100]
        adql = ("SELECT * FROM swiftmastr WHERE obsid IN (%s)"
                % ",".join("'%s'" % o for o in chunk))
        try:
            resp = requests.get(HEASARC_TAP_URL, timeout=120, params={
                "REQUEST": "doQuery", "LANG": "ADQL", "QUERY": adql})
            resp.raise_for_status()
            # use_names_over_ids: the obsid column's VOTable ID is
            # 'DataLinkID'; its name is 'obsid'.
            table = parse_single_table(io.BytesIO(resp.content)).to_table(
                use_names_over_ids=True)
        except Exception as exc:  # network, server or parse error
            print(f"[warn] Catalog lookup of OBSIDs failed ({exc}); "
                  f"continuing without catalog rows for "
                  f"{len(obsids) - len(rows)} OBSID(s).")
            return rows
        for r in table:
            row = {c: _cell(r[c]) for c in table.colnames}
            rows[row.get("obsid", "")] = row
    missing = [o for o in obsids if o not in rows]
    if missing:
        print(f"[warn] {len(missing)} OBSID(s) not in the Swift master "
              f"catalog: {' '.join(missing[:10])}"
              f"{' ...' if len(missing) > 10 else ''}")
    return rows


# ---------------------------------------------------------------------------
# UVOT filters and file selection
# ---------------------------------------------------------------------------

# Filter codes as they appear in file names (sw<OBSID>u<code>...) and in the
# catalog's uvot_expo_<code> columns.
UVOT_FILTERS = {
    "vv": "V", "bb": "B", "uu": "U",
    "w1": "UVW1", "m2": "UVM2", "w2": "UVW2",
    "wh": "WHITE", "gu": "UGRISM", "gv": "VGRISM", "mg": "MAGNIFIER",
}
PHOTOMETRIC_FILTERS = ["vv", "bb", "uu", "w1", "m2", "w2"]

DEFAULT_PRODUCTS = ["uvot/image", "uvot/event", "uvot/hk", "auxil"]

# Products whose files carry a filter code in their name.
_PER_FILTER_PRODUCTS = ("uvot/image", "uvot/event")
_FILTER_IN_NAME = re.compile(r"^sw\d{11}u([a-z0-9]{2})")


def make_file_filter(filters: list[str], skip_raw: bool = False):
    """
    A predicate accept(product, filename): keep files of the chosen filters
    in uvot/image and uvot/event, optionally without raw images; keep every
    file of the other products (housekeeping, auxil).
    """
    wanted = set(filters)

    def accept(product: str, fname: str) -> bool:
        if product in _PER_FILTER_PRODUCTS:
            m = _FILTER_IN_NAME.match(fname)
            if not m or m.group(1) not in wanted:
                return False
            if skip_raw and "_rw.img" in fname:
                return False
        return True

    return accept


def uvot_exposure(row: dict, code: str | None = None) -> float:
    """Catalog UVOT exposure (s): total, or of one filter code."""
    key = "uvot_exposure" if code is None else f"uvot_expo_{code}"
    try:
        return float(row.get(key) or 0)
    except ValueError:
        return 0.0


# ---------------------------------------------------------------------------
# Archive locations
# ---------------------------------------------------------------------------

HEASARC_DATA_BASE = "https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs"
UKSSDC_DATA_BASE = "https://www.swift.ac.uk/archive/reproc"

_url_cache: dict[str, str] = {}
_cache_lock = threading.Lock()

# MJD 40587 == 1970-01-01 (the Unix epoch), used to turn an MJD into a date.
_MJD_UNIX_EPOCH = 40587


def _mjd_to_date(value: str) -> date | None:
    try:
        mjd = float(str(value).split()[0])
    except (ValueError, IndexError):
        return None
    if not 30000 <= mjd <= 90000:
        return None
    return date(1970, 1, 1) + timedelta(days=int(mjd - _MJD_UNIX_EPOCH))


def archive_candidates(obsid: str, meta: dict | None) -> list[str]:
    """
    Base URLs to try for an OBSID. The HEASARC archive files observations
    under the year and month they started (obs/2013_04/00035017124/); the
    catalog gives start_time as an MJD. The next month is tried too, for
    observations that start just before a month boundary. The UKSSDC
    mirror (identical files) needs no date.
    """
    urls = []
    d = _mjd_to_date((meta or {}).get("start_time", ""))
    if d:
        nxt = (d.replace(day=1) + timedelta(days=32)).replace(day=1)
        for day in (d, nxt):
            urls.append(f"{HEASARC_DATA_BASE}/{day:%Y_%m}/{obsid}/")
    urls.append(f"{UKSSDC_DATA_BASE}/{obsid}/")
    return urls


def _probe_url(url: str) -> bool:
    """Return True if *url* exists (HTTP 200 on HEAD or GET)."""
    try:
        resp = _session().head(url, timeout=15, allow_redirects=True)
        if resp.status_code == 200:
            return True
        if resp.status_code == 405:  # no HEAD on directories: try GET
            resp = _session().get(url, timeout=15, stream=True)
            resp.close()
            return resp.status_code == 200
        return False
    except requests.RequestException:
        return False


def resolve_data_url(obsid: str, meta: dict | None) -> str | None:
    """The archive folder of an OBSID (ending in '/'), or None."""
    with _cache_lock:
        if obsid in _url_cache:
            return _url_cache[obsid]
    for url in archive_candidates(obsid, meta):
        if _probe_url(url):
            with _cache_lock:
                _url_cache[obsid] = url
            return url
    return None


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------

class DownloadError(Exception):
    """A file could not be downloaded intact."""


_thread_state = threading.local()


def _session() -> requests.Session:
    """
    One HTTP session per thread, so archive requests reuse connections: a
    re-run checks every file's size with a HEAD request, and with a new
    TLS connection for each of 11,516 files (3C 273) that took 4.5 min.
    """
    session = getattr(_thread_state, "session", None)
    if session is None:
        session = requests.Session()
        _thread_state.session = session
    return session


def list_remote_dir(url: str) -> list[str]:
    """File names in a simple Apache/nginx directory listing."""
    resp = _session().get(url, timeout=60)
    resp.raise_for_status()
    names = []
    for href in re.findall(r'href="([^"]+)"', resp.text):
        if href.startswith(("/", "?", "#", "http:", "https:", "mailto:")) \
                or href.endswith("/") or href in (".", ".."):
            continue
        names.append(href)
    return names


def get_remote_size(url: str) -> int | None:
    """Return Content-Length from a HEAD request, or None if unavailable."""
    try:
        resp = _session().head(url, timeout=15, allow_redirects=True)
        resp.raise_for_status()
        cl = resp.headers.get("Content-Length")
        return int(cl) if cl else None
    except (requests.RequestException, ValueError):
        return None


def gz_intact(path: Path) -> bool:
    """True if a .gz file decompresses to the end; other files: True."""
    if not path.name.endswith(".gz"):
        return True
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except (OSError, EOFError):
        return False


def download_file(url: str, dest: Path, overwrite: bool = False,
                  verify: bool = False, retries: int = 3) -> str:
    """
    Download url to dest. Returns 'downloaded' or 'skipped'; raises
    DownloadError if the file could not be fetched intact.

    An existing file is kept if its size matches the server's (with
    verify, only if it also decompresses cleanly); if the server gives no
    size, only if it decompresses cleanly. New data go to <dest>.part and
    are renamed only when complete and intact.
    """
    if dest.exists() and not overwrite:
        remote = get_remote_size(url)
        local = dest.stat().st_size
        if remote == local and (not verify or gz_intact(dest)):
            return "skipped"
        if remote is None and gz_intact(dest):
            return "skipped"
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    problem = "unknown error"
    try:
        for attempt in range(1, retries + 1):
            try:
                with _session().get(url, stream=True, timeout=120) as resp:
                    resp.raise_for_status()
                    total = int(resp.headers.get("Content-Length") or 0)
                    with open(part, "wb") as fh:
                        for chunk in resp.iter_content(chunk_size=1 << 16):
                            fh.write(chunk)
                size = part.stat().st_size
                if total and size != total:
                    raise DownloadError(f"incomplete ({size} of {total} bytes)")
                if not gz_intact(part):
                    raise DownloadError("corrupt (does not decompress)")
                os.replace(part, dest)
                return "downloaded"
            except requests.HTTPError as exc:
                problem = f"HTTP {exc.response.status_code}"
                if exc.response.status_code == 404:
                    break  # no point retrying
            except (requests.RequestException, DownloadError, OSError) as exc:
                problem = str(exc) or exc.__class__.__name__
            if attempt < retries:
                time.sleep(2 * attempt)
    finally:
        # never leave a partial file, whatever stopped the download
        part.unlink(missing_ok=True)
    raise DownloadError(problem)


def download_obsid(obsid: str, outdir: Path, meta: dict | None,
                   products: list[str], accept, overwrite: bool,
                   stop: threading.Event | None = None,
                   verify: bool = False) -> dict:
    """
    Download one OBSID's files. Returns a dict with counts ('downloaded',
    'skipped', 'filtered'), 'failed' (list of (file, reason)), 'notes'
    (list of str), 'bytes', 'url' and 'found' (False if the OBSID is in
    no archive).
    """
    res = {"obsid": obsid, "downloaded": 0, "skipped": 0, "filtered": 0,
           "failed": [], "notes": [], "bytes": 0, "url": None, "found": False}
    base = resolve_data_url(obsid, meta)
    if base is None:
        res["failed"].append(("(observation folder)",
                              "not found at HEASARC or UKSSDC"))
        return res
    res["found"], res["url"] = True, base
    for prod in products:
        prod_url = base + prod + "/"
        try:
            names = list_remote_dir(prod_url)
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 404:
                # uvot/event exists only when an exposure was in event mode
                if prod != "uvot/event":
                    res["notes"].append(f"no {prod} folder in the archive")
            else:
                res["failed"].append((prod + "/", f"listing failed: {exc}"))
            continue
        except requests.RequestException as exc:
            res["failed"].append((prod + "/", f"listing failed: {exc}"))
            continue
        local_dir = outdir / obsid / prod
        for fname in names:
            if stop is not None and stop.is_set():
                res["failed"].append((f"{prod}/...", "interrupted"))
                return res
            if not accept(prod, fname):
                res["filtered"] += 1
                continue
            dest = local_dir / fname
            try:
                outcome = download_file(prod_url + fname, dest,
                                        overwrite=overwrite, verify=verify)
            except DownloadError as exc:
                res["failed"].append((f"{prod}/{fname}", str(exc)))
                continue
            res[outcome] += 1
            if outcome == "downloaded":
                res["bytes"] += dest.stat().st_size
    return res


# ---------------------------------------------------------------------------
# Check against the catalog
# ---------------------------------------------------------------------------

def check_against_catalog(obsid: str, outdir: Path, meta: dict | None,
                          filters: list[str]) -> tuple[list[str], list[str]]:
    """
    Compare an OBSID's sky images with the catalog's exposure per filter.
    Returns (expected filter codes, problems): a filter with catalog
    exposure but no sky image, or a sky image without its exposure map.
    """
    if not meta:
        return [], ["no catalog row: per-filter check skipped"]
    image_dir = outdir / obsid / "uvot" / "image"
    expected = [c for c in filters if uvot_exposure(meta, c) > 0]
    problems = []
    for code in filters:
        sky = image_dir / f"sw{obsid}u{code}_sk.img.gz"
        expmap = image_dir / f"sw{obsid}u{code}_ex.img.gz"
        if code in expected and not sky.exists():
            problems.append(f"{UVOT_FILTERS[code]}: catalog exposure "
                            f"{uvot_exposure(meta, code):.0f} s but no sky "
                            f"image")
        if sky.exists() and not expmap.exists():
            problems.append(f"{UVOT_FILTERS[code]}: sky image without "
                            f"exposure map")
    return expected, problems


# ---------------------------------------------------------------------------
# Listing
# ---------------------------------------------------------------------------

def _get_obsid(obs: dict) -> str:
    """Extract the observation ID from a catalog row."""
    for key in ("obsid", "obs_id", "OBSID", "OBS_ID"):
        val = (obs.get(key) or "").strip()
        if val:
            return val.zfill(11)
    return ""


def _offset_arcmin(obs: dict, ra: float | None, dec: float | None) -> str:
    try:
        ra0, dec0 = float(obs["ra"]), float(obs["dec"])
    except (KeyError, ValueError):
        return ""
    if ra is None or dec is None:
        return ""
    r1, d1, r2, d2 = map(math.radians, (ra0, dec0, ra, dec))
    cosang = (math.sin(d1) * math.sin(d2)
              + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2))
    return "%.1f" % (math.degrees(math.acos(max(-1.0, min(1.0, cosang)))) * 60)


def print_observations(observations: list[dict], ra=None, dec=None):
    """One row per observation: date, offset, target, UVOT exposure per filter."""
    head = ("  %-11s %-10s %6s  %-20s %7s  " % ("OBSID", "Date", "Off(')",
                                                 "Target", "UVOT(ks)")
            + " ".join("%5s" % UVOT_FILTERS[c] for c in PHOTOMETRIC_FILTERS)
            + "  other")
    print("\n" + head)
    print("  " + "-" * (len(head) - 2))
    for obs in observations:
        d = _mjd_to_date(obs.get("start_time", ""))
        cells = []
        for c in PHOTOMETRIC_FILTERS:
            e = uvot_exposure(obs, c)
            cells.append("%5.2f" % (e / 1e3) if e > 0 else "%5s" % "-")
        other = ", ".join("%s %.2f" % (UVOT_FILTERS[c], uvot_exposure(obs, c) / 1e3)
                          for c in ("wh", "gu", "gv", "mg")
                          if uvot_exposure(obs, c) > 0)
        print("  %-11s %-10s %6s  %-20s %7.2f  %s  %s" % (
            _get_obsid(obs), d.isoformat() if d else "?",
            _offset_arcmin(obs, ra, dec), (obs.get("name") or "")[:20],
            uvot_exposure(obs) / 1e3, " ".join(cells), other))
    print("  " + "-" * (len(head) - 2))
    print("  Exposures in ks per filter (catalog values). Off(') = pointing "
          "offset from the source position.\n")


# ---------------------------------------------------------------------------
# OBSID file and date window (from swift_xrt_download.py)
# ---------------------------------------------------------------------------

def parse_obsid_file(filepath: str) -> list[str]:
    """Read observation IDs from a text file (one per line).

    Blank lines and lines starting with '#' are ignored.  If a line contains
    whitespace-separated columns (e.g. "00035393001  Crab"), only the first
    token is taken as the obsid.
    """
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"Obsid file not found: {filepath}")

    obsids: list[str] = []
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            token = line.split()[0]
            if not token.isdigit():
                print(f"[warn] {filepath}:{lineno}: '{token}' doesn't look like "
                      f"an obsid, skipping.")
                continue
            obsids.append(token.zfill(11))

    if not obsids:
        raise ValueError(f"No valid observation IDs found in {filepath}")

    print(f"[info] Read {len(obsids)} observation ID(s) from {filepath}")
    return obsids


def filter_by_date(observations: list[dict],
                   start_date: date | None,
                   end_date: date | None) -> tuple[list[dict], int]:
    """Keep observations with start date in [start_date, end_date).

    Rows with no parseable date are dropped when a filter is active; the
    count of such rows is returned so the caller can warn about them.
    """
    if not start_date and not end_date:
        return observations, 0
    kept: list[dict] = []
    undated = 0
    for obs in observations:
        d = _mjd_to_date(obs.get("start_time", ""))
        if d is None:
            undated += 1
            continue
        if start_date and d < start_date:
            continue
        if end_date and d >= end_date:
            continue
        kept.append(obs)
    return kept, undated


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(path: Path, results: list[dict], args, filters):
    """download_report.txt: one line per OBSID, then the problems."""
    with open(path, "w") as fh:
        fh.write("# swift_uvot_download.py report, %s\n" % datetime.now(
            timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
        fh.write("# command: %s\n" % " ".join(sys.argv))
        fh.write("# filters: %s; products: %s\n" % (
            " ".join(UVOT_FILTERS[c] for c in filters), " ".join(args.products)))
        fh.write("# status: ok | failed (re-run to retry) | not_found | "
                 "check (files missing compared with the catalog) | "
                 "no_data (no exposure in the chosen filters)\n")
        fh.write("%-11s %-9s %5s %5s %5s  %-28s %s\n" % (
            "OBSID", "status", "new", "kept", "fail", "filters (catalog)",
            "source"))
        for r in results:
            fh.write("%-11s %-9s %5d %5d %5d  %-28s %s\n" % (
                r["obsid"], r["status"], r["downloaded"], r["skipped"],
                len(r["failed"]),
                " ".join(UVOT_FILTERS[c] for c in r["expected"]) or "-",
                "UKSSDC" if (r["url"] or "").startswith(UKSSDC_DATA_BASE)
                else ("HEASARC" if r["url"] else "-")))
        problems = [(r["obsid"], f"FAILED {name}: {why}")
                    for r in results for name, why in r["failed"]]
        problems += [(r["obsid"], f"CHECK {p}")
                     for r in results for p in r["problems"]]
        problems += [(r["obsid"], f"note: {n}")
                     for r in results for n in r["notes"]]
        if problems:
            fh.write("\n# Problems and notes\n")
            for obsid, text in problems:
                fh.write(f"{obsid}  {text}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _isodate(s: str) -> date:
    """argparse type converter: parse YYYY-MM-DD, error clearly otherwise."""
    try:
        return date.fromisoformat(s)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"invalid date '{s}': expected YYYY-MM-DD")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Download Swift UVOT data from the HEASARC archive.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--name", "-n", type=str,
                     help="Source name to resolve (e.g. '3C 273'), then a "
                          "cone search of the Swift master catalog. Every "
                          "pointing within --radius is returned, including "
                          "pointings of other targets.")
    grp.add_argument("--ra", type=float,
                     help="Right Ascension in decimal degrees (J2000). "
                          "Must also give --dec.")
    grp.add_argument("--obsid", type=str,
                     help="Download a single observation by its ID.")
    grp.add_argument("--obsid-file", type=str, metavar="FILE",
                     help="Text file with one observation ID per line.")

    p.add_argument("--dec", type=float, default=None,
                   help="Declination in decimal degrees (J2000).")
    p.add_argument("--radius", "-r", type=float, default=12.0,
                   help="Cone-search radius in arcminutes (default: 12).")
    p.add_argument("--start-date", type=_isodate, default=None,
                   metavar="YYYY-MM-DD",
                   help="Keep observations on or after this date (UTC, "
                        "inclusive; catalog queries only).")
    p.add_argument("--end-date", type=_isodate, default=None,
                   metavar="YYYY-MM-DD",
                   help="Keep observations strictly before this date "
                        "(UTC, exclusive; catalog queries only).")
    p.add_argument("--filters", nargs="+", default=PHOTOMETRIC_FILTERS,
                   choices=sorted(UVOT_FILTERS), metavar="CODE",
                   help="UVOT filters to download, by file-name code: "
                        "vv bb uu w1 m2 w2 (default, the photometric "
                        "filters), also wh (WHITE), gu gv (grisms), mg "
                        "(magnifier).")
    p.add_argument("--products", nargs="+", default=DEFAULT_PRODUCTS,
                   help="Archive folders to download (default: %s)."
                        % " ".join(DEFAULT_PRODUCTS))
    p.add_argument("--skip-raw", action="store_true",
                   help="Don't download raw images (*_rw.img.gz).")
    p.add_argument("--outdir", "-o", type=str, default="UVOT_input",
                   help="Output folder (default: UVOT_input).")
    p.add_argument("--nproc", type=int, default=1,
                   help="Observations to download at once (default: 1; "
                        "please keep it modest, e.g. 4).")
    p.add_argument("--list-only", "-l", action="store_true",
                   help="List the matching observations without "
                        "downloading.")
    p.add_argument("--max-obs", "-m", type=int, default=None,
                   help="Download at most this many observations.")
    p.add_argument("--overwrite", action="store_true",
                   help="Re-download files even if they are already there "
                        "and intact.")
    p.add_argument("--verify", action="store_true",
                   help="Also decompress files that are already there to "
                        "check them (slower: ~5 min for all of 3C 273); "
                        "by default a size match with the server is enough. "
                        "New downloads are always checked.")
    return p


def main():
    # One line at a time even when redirected to a file (nohup, logs):
    # block buffering made a running download look stuck.
    sys.stdout.reconfigure(line_buffering=True)
    parser = build_parser()
    args = parser.parse_args()
    filters = list(dict.fromkeys(args.filters))

    if (args.start_date or args.end_date) and (args.obsid or args.obsid_file):
        print("[warn] --start-date/--end-date are ignored for --obsid / "
              "--obsid-file.", file=sys.stderr)

    ra = dec = None
    if args.obsid or args.obsid_file:
        if args.obsid:
            if not args.obsid.strip().isdigit():
                parser.error(f"--obsid {args.obsid!r} is not an observation ID")
            obsid_list = [args.obsid.strip().zfill(11)]
        else:
            obsid_list = parse_obsid_file(args.obsid_file)
        meta_by_obsid = query_swift_master_obsids(obsid_list)
        observations = [meta_by_obsid[o] for o in obsid_list
                        if o in meta_by_obsid]
    else:
        if args.name:
            ra, dec = resolve_name(args.name)
        else:
            if args.dec is None:
                parser.error("--dec is required when using --ra")
            ra, dec = args.ra, args.dec
        observations = query_swift_master(ra, dec, radius_arcmin=args.radius)
        if not observations:
            print("[info] No observations found. Try increasing --radius.")
            sys.exit(0)
        if args.start_date or args.end_date:
            n_before = len(observations)
            observations, undated = filter_by_date(
                observations, args.start_date, args.end_date)
            print(f"[info] {len(observations)} of {n_before} observation(s) "
                  f"in the date window.")
            if undated:
                print(f"[warn] {undated} observation(s) had no parseable date "
                      f"and were dropped by the date filter.")
        observations.sort(key=lambda o: o.get("start_time", ""))
        no_uvot = [o for o in observations if uvot_exposure(o) <= 0]
        observations = [o for o in observations if uvot_exposure(o) > 0]
        if no_uvot:
            print(f"[info] {len(no_uvot)} observation(s) have no UVOT exposure "
                  f"and are left out.")
        obsid_list = [_get_obsid(o) for o in observations]
        meta_by_obsid = {_get_obsid(o): o for o in observations}

    if observations:
        print_observations(observations, ra, dec)
    if args.list_only:
        print(f"[info] {len(obsid_list)} observation(s). "
              "Run without --list-only to download.")
        sys.exit(0)
    if not obsid_list:
        print("[info] No observations to download.")
        sys.exit(0)
    if args.max_obs is not None:
        obsid_list = obsid_list[:args.max_obs]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    accept = make_file_filter(filters, skip_raw=args.skip_raw)
    print(f"[info] Downloading {len(obsid_list)} observation(s) to "
          f"{outdir.resolve()}")
    print(f"[info] Filters: {' '.join(UVOT_FILTERS[c] for c in filters)}; "
          f"folders: {' '.join(args.products)}"
          f"{'; no raw images' if args.skip_raw else ''}")

    results: list[dict] = []
    print_lock = threading.Lock()
    start = time.monotonic()

    stop = threading.Event()

    def work(obsid):
        meta = meta_by_obsid.get(obsid)
        t0 = time.monotonic()
        if meta and sum(uvot_exposure(meta, c) for c in filters) <= 0:
            res = {"obsid": obsid, "status": "no_data", "downloaded": 0,
                   "skipped": 0, "filtered": 0, "failed": [], "bytes": 0,
                   "url": None, "found": True, "expected": [], "problems": [],
                   "notes": ["no exposure in the chosen filters (catalog); "
                             "nothing downloaded"]}
        else:
            res = download_obsid(obsid, outdir, meta, args.products, accept,
                                 args.overwrite, stop, verify=args.verify)
            res["expected"], res["problems"] = ([], []) if not res["found"] \
                else check_against_catalog(obsid, outdir, meta, filters)
        if res.get("status") == "no_data":
            pass
        elif not res["found"]:
            res["status"] = "not_found"
        elif res["failed"]:
            res["status"] = "failed"
        elif res["problems"]:
            res["status"] = "check"
        else:
            res["status"] = "ok"
        with print_lock:
            results.append(res)
            n = len(results)
            elapsed = time.monotonic() - start
            eta = elapsed / n * (len(obsid_list) - n)
            print(f"[{n}/{len(obsid_list)}] {obsid}  {res['status']:<9} "
                  f"{res['downloaded']} new, {res['skipped']} kept, "
                  f"{len(res['failed'])} failed, "
                  f"{res['bytes'] / 1e6:.1f} MB, {time.monotonic() - t0:.0f} s"
                  f"  (~{eta / 60:.0f} min left)")
            for name, why in res["failed"]:
                print(f"      FAILED {name}: {why}")
            for p in res["problems"]:
                print(f"      CHECK  {p}")
        return res

    interrupted = False
    pool = ThreadPoolExecutor(max_workers=max(1, args.nproc))
    try:
        for fut in as_completed([pool.submit(work, o) for o in obsid_list]):
            fut.result()
    except KeyboardInterrupt:
        interrupted = True
        stop.set()
        print("\n[interrupted] Stopping after the files in progress. Re-run "
              "the same command to resume: complete files are kept, partial "
              "ones are deleted.")
    finally:
        pool.shutdown(wait=True, cancel_futures=True)

    results.sort(key=lambda r: r["obsid"])
    report = outdir / "download_report.txt"
    write_report(report, results, args, filters)

    count = {s: sum(r["status"] == s for r in results)
             for s in ("ok", "failed", "not_found", "check", "no_data")}
    n_files = sum(r["downloaded"] for r in results)
    n_failed = sum(len(r["failed"]) for r in results)
    print()
    print(f"[summary] {len(results)} of {len(obsid_list)} observation(s): "
          f"{count['ok']} ok, {count['check']} to check, {count['failed']} "
          f"with failed files, {count['not_found']} not found, "
          f"{count['no_data']} with no data in the chosen filters")
    print(f"[summary] Files: {n_files} downloaded "
          f"({sum(r['bytes'] for r in results) / 1e9:.2f} GB), "
          f"{sum(r['skipped'] for r in results)} already there, "
          f"{n_failed} failed")
    print(f"[summary] Report: {report}")
    if count["check"]:
        print("[summary] 'check' means the archive lacks files the catalog "
              "lists (see the report); re-running does not change that.")
    if interrupted:
        sys.exit(130)
    sys.exit(1 if (n_failed or count["not_found"]
                   or len(results) < len(obsid_list)) else 0)


if __name__ == "__main__":
    main()
