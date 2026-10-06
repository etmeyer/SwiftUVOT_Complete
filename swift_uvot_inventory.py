#!/usr/bin/env python3
"""
swift_uvot_inventory.py -- Step 3: an inventory of every UVOT exposure.

A UVOT sky image (sw<OBSID>u<filter>_sk.img.gz, plus _sk_01.img.gz and so
on when the archive split a long observation's event data) holds one FITS
extension per exposure, and the exposures in one file can differ in frame
time, window, data mode and aspect correction. This step opens every extension
of every sky image under --indir and records, for each:

  * what it is: filter, data mode (IMAGE, EVENT, IMAGEEVENT), OBS_MODE,
    exposure, frame time, window, binning, aspect correction, processing
    version, time (UTC midpoint);
  * where it sits in its observation: the snapshot (exposures closer than
    5 minutes) and whether it is the first exposure of that snapshot
    (those can be trailed: the spacecraft may still be settling);
  * whether it duplicates another extension (same EXPID: the image-mode
    and event-mode copy of one exposure);
  * whether the source -- your --ra/--dec -- is inside the exposed field:
    the fraction of the 5" source circle and of the 27.5-35" background
    annulus that is fully exposed (within 1 % of the local maximum of the
    exposure map, uvotsource's own tolerance), and the distance from the
    source to the nearest unexposed pixel.

Every extension gets exactly one status, the start of the pipeline's
ledger (each later step accounts for every extension it was given):

    candidate     goes on to positions and photometry (Steps 4-5)
    not_covered   the source is outside the exposed field
    partial       part of the 5" source circle is not fully exposed
    settling      OBS_MODE is SETTLING
    nonphot       WHITE, a grism, the magnifier or BLOCKED
    no_expmap     no exposure-map extension with the same name
    unreadable    the file or extension could not be read

The step reads only FITS files (no HEASoft). With network it also checks
that every exposure in HEASARC's UVOT exposure log has a sky image, and
lists the Swift master catalog's exposure per OBSID and filter.

Output (in --outdir):
    uvot_inventory.txt           one row per extension
    uvot_inventory_compact.txt   one row per OBSID and filter (also printed)

Exit status: 1 if any file or extension could not be read.

Usage:
    swift_uvot_inventory.py --ra 187.2779 --dec 2.0524
    swift_uvot_inventory.py --ra 187.2779 --dec 2.0524 --indir UVOT_input \
        --outdir UVOT_output --nproc 8
"""

import argparse
import glob
import gzip
import io
import os
import re
import sys
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning

import swift_uvot_env as uenv
from swift_uvot_tables import write_table

warnings.simplefilter('ignore', FITSFixedWarning)

# File-name filter codes (sw<OBSID>u<code>_sk.img.gz) and filter names.
FILTER_NAMES = {
    'vv': 'V', 'bb': 'B', 'uu': 'U', 'w1': 'UVW1', 'm2': 'UVM2', 'w2': 'UVW2',
    'wh': 'WHITE', 'gu': 'UGRISM', 'gv': 'VGRISM', 'mg': 'MAGNIFIER',
    'bl': 'BLOCKED',
}
PHOTOMETRIC = ('V', 'B', 'U', 'UVW1', 'UVM2', 'UVW2')

SRC_RADIUS = 5.0              # arcsec: the calibrated photometric aperture
BKG_INNER, BKG_OUTER = 27.5, 35.0   # arcsec: Poole et al. (2008) annulus
EXPOSED_FRACTION = 0.99       # pixel counts as exposed above 99 % of the max
CLEARANCE_MAX = 60.0          # arcsec searched for the nearest unexposed pixel
SNAPSHOT_GAP = 300.0          # s: a longer gap between exposures = new snapshot

STATUSES = ('candidate', 'not_covered', 'partial', 'settling', 'nonphot',
            'no_expmap', 'unreadable')

COLUMNS = ['obsid', 'filter', 'hdu', 'extname', 'mode', 'obs_mode', 'expid',
           'date_mid', 'mjd_mid', 'tstart', 'tstop', 'exposure', 'ontime',
           'frametime',
           'window', 'binning', 'aspcorr', 'loss', 'snapshot', 'first',
           'off_pnt', 'src_cover', 'bkg_cover', 'clearance', 'dup_of', 'event',
           'procver', 'file', 'status', 'reason']

# A long event-mode observation can have a second sky image, _sk_01.img.gz,
# made from a second event file (_uf_01.evt.gz).
_SKY_NAME = re.compile(r'^sw(\d{11})u([a-z0-9]{2})_sk(_\d+)?\.img(\.gz)?$')


def expmap_path(sky_path):
    """The exposure map of a sky image (..._ex[_NN].img.gz)."""
    head, name = os.path.split(sky_path)
    return os.path.join(head, re.sub(r'_sk(_\d+)?\.img', r'_ex\1.img', name))


def _loss(header):
    """Summed data-loss keywords (s), or '-' if they make no sense."""
    try:
        loss = sum(float(header.get(k, 0) or 0)
                   for k in ('TOSSLOSS', 'STALLOSS', 'BLOCLOSS'))
    except (TypeError, ValueError):
        return '-'
    # one 3C 273 header has STALLOSS = 4e252
    return ('%.1f' % loss if 0 <= loss <= header['TSTOP'] - header['TSTART']
            else '-')


def _coverage(expmap, header, src):
    """
    (src_cover, bkg_cover, clearance in arcsec or None if > CLEARANCE_MAX)
    of the source in one exposure-map extension, using the sky image's WCS.
    """
    x, y = WCS(header).world_to_pixel(src)
    x, y = float(x), float(y)
    scale = abs(header['CDELT1']) * 3600.0
    ny, nx = expmap.shape
    r = int(np.ceil(CLEARANCE_MAX / scale)) + 2
    xc, yc = int(round(x)), int(round(y))
    # Cutout around the source; pixels outside the image count as unexposed.
    cut = np.zeros((2 * r + 1, 2 * r + 1))
    x0, x1 = max(0, xc - r), min(nx, xc + r + 1)
    y0, y1 = max(0, yc - r), min(ny, yc + r + 1)
    if x0 < x1 and y0 < y1:
        cut[y0 - (yc - r):y1 - (yc - r), x0 - (xc - r):x1 - (xc - r)] = \
            np.nan_to_num(expmap[y0:y1, x0:x1])
    yy, xx = np.mgrid[yc - r:yc + r + 1, xc - r:xc + r + 1]
    dist = np.hypot(xx - x, yy - y) * scale
    peak = cut[dist <= BKG_OUTER].max() if (dist <= BKG_OUTER).any() else 0.0
    if peak <= 0:
        return 0.0, 0.0, 0.0
    exposed = cut >= EXPOSED_FRACTION * peak
    src_px = dist <= SRC_RADIUS
    bkg_px = (dist >= BKG_INNER) & (dist <= BKG_OUTER)
    unexposed = dist[(~exposed) & (dist <= CLEARANCE_MAX)]
    clearance = float(unexposed.min()) if unexposed.size else None
    return (float(exposed[src_px].mean()), float(exposed[bkg_px].mean()),
            clearance)


def _open_fits(path):
    """
    Open a FITS file, reading a .gz completely first: astropy silently
    shows a truncated .gz as having fewer extensions (a cut-off sky image
    looked like a file with no exposures at all), whereas gzip raises an
    error at the missing end or a bad checksum.
    """
    if path.endswith('.gz'):
        with gzip.open(path, 'rb') as fh:
            return fits.open(io.BytesIO(fh.read()))
    return fits.open(path)


def _date_mid(header):
    """(ISO UTC midpoint of DATE-OBS..DATE-END, MJD) or ('-', '-')."""
    try:
        t0 = Time(header['DATE-OBS'], format='isot', scale='utc')
        t1 = Time(header['DATE-END'], format='isot', scale='utc')
    except (KeyError, ValueError):
        return '-', '-'
    mid = t0 + (t1 - t0) / 2
    return mid.isot[:19], '%.5f' % mid.mjd


def inventory_file(sky_path, ra, dec):
    """Rows (dicts) for every extension of one sky image."""
    m = _SKY_NAME.match(os.path.basename(sky_path))
    obsid, code = m.group(1), m.group(2)
    filt = FILTER_NAMES.get(code, code.upper())
    base = {'obsid': obsid, 'filter': filt, 'file': sky_path}
    src = SkyCoord(ra, dec, unit='deg')
    ex_path = expmap_path(sky_path)
    rows = []
    try:
        sky = _open_fits(sky_path)
        procver = sky[0].header.get('PROCVER', '-')
        n_hdu = len(sky) - 1
    except Exception as exc:  # truncated or corrupt file
        return [dict(base, hdu='-', extname='-', status='unreadable',
                     reason='cannot read %s: %s' % (
                         os.path.basename(sky_path), exc))]
    if n_hdu < 1:
        return [dict(base, hdu='-', extname='-', status='unreadable',
                     reason='no extensions in %s' % os.path.basename(sky_path))]
    expmaps, ex_problem = {}, None
    if os.path.exists(ex_path):
        try:
            ex = _open_fits(ex_path)
            expmaps = {h.header.get('EXTNAME'): h for h in ex[1:]}
        except Exception as exc:
            ex_problem = 'cannot read %s: %s' % (os.path.basename(ex_path), exc)

    for idx in range(1, n_hdu + 1):
        row = dict(base, hdu=idx, procver=procver)
        try:
            hdu = sky[idx]
            h = hdu.header
            row.update({
                'extname': h.get('EXTNAME', '-'),
                'mode': h.get('DATAMODE', '-'),
                'obs_mode': h.get('OBS_MODE', '-'),
                'expid': h.get('EXPID', '-'),
                'tstart': '%.3f' % h['TSTART'],
                'tstop': '%.3f' % h['TSTOP'],
                'exposure': '%.2f' % h['EXPOSURE'],
                # before dead-time correction; the catalog sums this
                'ontime': '%.2f' % h.get('ONTIME', h.get('TELAPSE', 0.0)),
                'frametime': '%.6f' % h.get('FRAMTIME', 0.0110322),
                'window': '%sx%s' % (h.get('WINDOWDX', '?'),
                                     h.get('WINDOWDY', '?')),
                'binning': h.get('BINX', '-'),
                'aspcorr': str(h.get('ASPCORR', 'NONE')).split()[0],
                'loss': _loss(h),
            })
            row['date_mid'], row['mjd_mid'] = _date_mid(h)
            try:
                pnt = SkyCoord(h['RA_PNT'], h['DEC_PNT'], unit='deg')
                row['off_pnt'] = '%.2f' % pnt.separation(src).arcmin
            except KeyError:
                row['off_pnt'] = '-'

            reasons = []
            if row['loss'] == '-':
                reasons.append('data-loss keywords out of range')
            if filt not in PHOTOMETRIC:
                status = 'nonphot'
                reasons.append('%s is not a photometric filter' % filt)
            elif str(row['obs_mode']).upper() == 'SETTLING':
                status = 'settling'
                reasons.append('OBS_MODE=SETTLING')
            elif ex_problem:
                status = 'unreadable'
                reasons.append(ex_problem)
            elif row['extname'] not in expmaps:
                status = 'no_expmap'
                reasons.append('no extension %s in %s' % (
                    row['extname'], os.path.basename(ex_path)))
            elif hdu.data is None or getattr(hdu.data, 'ndim', 0) != 2:
                status = 'unreadable'
                reasons.append('no 2-D image in this extension')
            else:
                emap = expmaps[row['extname']].data
                if emap is None or emap.shape != hdu.data.shape:
                    status = 'no_expmap'
                    reasons.append('exposure map does not match the image')
                else:
                    sc, bc, clear = _coverage(emap, h, src)
                    row['src_cover'] = '%.2f' % sc
                    row['bkg_cover'] = '%.2f' % bc
                    row['clearance'] = ('>%d' % CLEARANCE_MAX if clear is None
                                        else '%.0f' % clear)
                    if sc == 0:
                        status = 'not_covered'
                        reasons.append('source outside the exposed field')
                    elif sc < 1:
                        status = 'partial'
                        reasons.append('%.0f%% of the 5" circle fully exposed'
                                       % (100 * sc))
                    else:
                        status = 'candidate'
                        if bc < 1:
                            reasons.append('background annulus %.0f%% exposed'
                                           % (100 * bc))
            if h['EXPOSURE'] < 1.0:
                # e.g. a 0.01 s exposure (3C 273 has one, with a broken WCS)
                reasons.append('exposure only %.2f s' % h['EXPOSURE'])
            row['status'] = status
            row['reason'] = '; '.join(reasons)
        except Exception as exc:
            row.update(status='unreadable', reason='cannot read: %s' % exc)
        rows.append(row)
    return rows


def annotate_observation(rows, event_codes):
    """Snapshot numbers, first-of-snapshot flags, duplicates, event files."""
    timed = sorted((r for r in rows if r.get('tstart') not in (None, '-')),
                   key=lambda r: float(r['tstart']))
    snapshot, last_stop = 0, None
    for r in timed:
        start, stop = float(r['tstart']), float(r['tstop'])
        first = last_stop is None or start - last_stop > SNAPSHOT_GAP
        if first:
            snapshot += 1
        r['snapshot'] = snapshot
        r['first'] = 'yes' if first else 'no'
        last_stop = stop if last_stop is None else max(last_stop, stop)
    by_expid = defaultdict(list)
    for r in rows:
        if r.get('expid') not in (None, '-'):
            by_expid[str(r['expid'])].append(r)
    for group in by_expid.values():
        if len(group) > 1:
            for r in group:
                r['dup_of'] = ','.join(o['extname'] for o in group if o is not r)
    for r in rows:
        if str(r.get('mode', '')).upper() in ('EVENT', 'IMAGEEVENT'):
            code = [c for c, n in FILTER_NAMES.items() if n == r['filter']]
            r['event'] = 'yes' if code and code[0] in event_codes else 'no'


def catalog_exposures(obsids):
    """{obsid: catalog row} from the Swift master catalog, or {} on failure."""
    try:
        from swift_uvot_download import query_swift_master_obsids
    except ImportError as exc:
        print('[warn] cannot query the catalog (%s)' % exc)
        return {}
    return query_swift_master_obsids(sorted(obsids))


def exposure_log(obsids):
    """HEASARC's UVOT exposure log for these OBSIDs, or None on failure."""
    try:
        from swift_uvot_download import query_uvot_exposure_log
    except ImportError as exc:
        print('[warn] cannot query the exposure log (%s)' % exc)
        return None
    return query_uvot_exposure_log(sorted(obsids))


def read_download_report(indir):
    """{obsid: status} from Step 2's download_report.txt, if present."""
    path = os.path.join(indir, 'download_report.txt')
    status = {}
    if not os.path.exists(path):
        return status
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2 and re.fullmatch(r'\d{11}', parts[0]) \
                    and parts[1] in ('ok', 'failed', 'not_found', 'check',
                                     'no_data'):
                status[parts[0]] = parts[1]
    return status


def compact_rows(rows, catalog):
    """One summary row per OBSID and filter."""
    groups = defaultdict(list)
    for r in rows:
        groups[(r['obsid'], r['filter'])].append(r)
    out = []
    for (obsid, filt), rs in sorted(groups.items()):
        expo = sum(float(r['exposure']) for r in rs
                   if r.get('exposure') not in (None, '-'))
        cand = [r for r in rs if r['status'] == 'candidate']
        frames = sorted({'%.1f' % (1e3 * float(r['frametime'])) for r in rs
                         if r.get('frametime') not in (None, '-')}, key=float)
        asp = Counter(r.get('aspcorr', '-') for r in rs)
        stat = Counter(r['status'] for r in rs)
        dates = [r['date_mid'] for r in rs
                 if r.get('date_mid') not in (None, '-')]
        row = {
            'obsid': obsid,
            'date': min(dates)[:10] if dates else '-',
            'filter': filt,
            'n_ext': len(rs),
            'n_cand': len(cand),
            'exp_s': '%.0f' % expo,
            'cand_exp_s': '%.0f' % sum(float(r['exposure']) for r in cand),
            'ontime_s': '%.0f' % sum(float(r['ontime']) for r in rs
                                     if r.get('ontime') not in (None, '-')),
            'frame_ms': ','.join(frames) or '-',
            'aspcorr': ','.join('%s%d' % (('D' if k == 'DIRECT' else
                                           'N' if k == 'NONE' else k[:1]), n)
                                for k, n in sorted(asp.items())),
            'statuses': ','.join('%s:%d' % (k, n) for k, n in sorted(stat.items())),
        }
        code = [c for c, n in FILTER_NAMES.items() if n == filt]
        cat = catalog.get(obsid)
        if cat and code:
            try:
                cexp = float(cat.get('uvot_expo_%s' % code[0]) or 0)
            except ValueError:
                cexp = 0.0
            row['cat_exp_s'] = '%.0f' % cexp
        out.append(row)
    return out


COMPACT_COLUMNS = ['obsid', 'date', 'filter', 'n_ext', 'n_cand', 'exp_s',
                   'cand_exp_s', 'ontime_s', 'cat_exp_s', 'frame_ms', 'aspcorr',
                   'statuses']


def main(argv=None):
    sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser(
        description='Step 3: inventory of every UVOT exposure.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--ra', type=float, required=True,
                        help='source RA, decimal degrees (J2000)')
    parser.add_argument('--dec', type=float, required=True,
                        help='source Dec, decimal degrees (J2000)')
    parser.add_argument('--indir', default='UVOT_input',
                        help='Step 2 output folder (default: UVOT_input)')
    parser.add_argument('--outdir', default='UVOT_output',
                        help='where to write the inventory (default: '
                             'UVOT_output)')
    parser.add_argument('--nproc', type=int, default=8,
                        help='files read at once (default: 8)')
    parser.add_argument('--no-catalog', action='store_true',
                        help="don't compare with HEASARC's master catalog "
                             "and exposure log (no network needed)")
    args = parser.parse_args(argv)

    sky_files = sorted(
        p for p in glob.glob(os.path.join(args.indir, '*', 'uvot', 'image',
                                          'sw*_sk*.img*'))
        if _SKY_NAME.match(os.path.basename(p)))
    if not sky_files:
        print('ERROR: no sky images (*/uvot/image/sw*_sk.img.gz) under %s'
              % args.indir, file=sys.stderr)
        return 1
    obsid_dirs = sorted(d for d in os.listdir(args.indir)
                        if re.fullmatch(r'\d{11}', d))
    print('[info] %d sky images in %d OBSID folders under %s'
          % (len(sky_files), len(obsid_dirs), args.indir))

    with ProcessPoolExecutor(max_workers=max(1, args.nproc)) as pool:
        per_file = list(pool.map(inventory_file, sky_files,
                                 [args.ra] * len(sky_files),
                                 [args.dec] * len(sky_files), chunksize=4))
    rows = [r for rs in per_file for r in rs]

    by_obsid = defaultdict(list)
    for r in rows:
        by_obsid[r['obsid']].append(r)
    for obsid, rs in by_obsid.items():
        event_dir = os.path.join(args.indir, obsid, 'uvot', 'event')
        codes = {m.group(1) for m in (re.match(r'sw\d{11}u([a-z0-9]{2})',
                                               os.path.basename(p))
                                      for p in glob.glob(event_dir + '/*'))
                 if m}
        annotate_observation(rs, codes)
    rows.sort(key=lambda r: (r['obsid'], float(r['tstart'])
                             if r.get('tstart') not in (None, '-') else 0))

    catalog = {} if args.no_catalog else catalog_exposures(by_obsid)
    log = None if args.no_catalog else exposure_log(by_obsid)
    compact = compact_rows(rows, catalog)
    downloads = read_download_report(args.indir)

    os.makedirs(args.outdir, exist_ok=True)
    counts = Counter(r['status'] for r in rows)
    prov = uenv.provenance()
    comments = [
        'swift_uvot_inventory.py, %s' % datetime.now(timezone.utc)
        .strftime('%Y-%m-%dT%H:%M:%SZ'),
        'command: %s' % ' '.join(sys.argv),
        'source: RA %.6f Dec %.6f; input %s' % (args.ra, args.dec,
                                                os.path.abspath(args.indir)),
        'pipeline %s' % prov['PIPECOMM'],
        '%d extensions in %d sky images: %s' % (
            len(rows), len(sky_files),
            ', '.join('%s %d' % (s, counts[s]) for s in STATUSES if counts[s])),
        'src_cover / bkg_cover: fraction of the 5" circle / 27.5-35" annulus '
        'fully exposed; clearance: arcsec to the nearest unexposed pixel',
        'first: first exposure of its snapshot (in event mode often trailed); '
        'dup_of: other extension with the same EXPID; event: event file '
        'present (event-mode exposures)',
    ]
    inv_path = os.path.join(args.outdir, 'uvot_inventory.txt')
    write_table(inv_path, rows, COLUMNS, comments)
    compact_path = os.path.join(args.outdir, 'uvot_inventory_compact.txt')
    write_table(compact_path, compact, COMPACT_COLUMNS,
                comments[:4] + ['one row per OBSID and filter; exposures in s '
                                '(exp_s dead-time corrected, ontime_s not); '
                                'cat_exp_s: Swift master catalog (on-time); '
                                'aspcorr: D=DIRECT N=NONE'])

    # --- report -------------------------------------------------------------
    with open(compact_path) as fh:
        lines = [ln.rstrip('\n') for ln in fh if not ln.startswith('#')]
    print()
    for ln in lines:
        print('  ' + ln)
    print()
    print('[summary] %d extensions: %s' % (
        len(rows), ', '.join('%s %d' % (s, counts[s]) for s in STATUSES
                             if counts[s])))
    by_filter = Counter((r['filter'], r['status']) for r in rows)
    for filt in PHOTOMETRIC:
        n = {s: by_filter[(filt, s)] for s in STATUSES}
        if sum(n.values()):
            print('          %-5s %s' % (filt, ', '.join(
                '%s %d' % (s, v) for s, v in n.items() if v)))
    # Only event-mode openers tend to be trailed: image mode is corrected
    # for drift on board, but the archive builds an event-mode sky image
    # with one pointing for the whole exposure (80 % of 3C 273's 194).
    openers = [r for r in rows if r['status'] == 'candidate'
               and r.get('first') == 'yes']
    event_openers = sum(1 for r in openers
                        if str(r.get('mode', '')).upper() == 'EVENT')
    dups = sum(1 for r in rows if r.get('dup_of') not in (None, '-'))
    print('[summary] %d candidates open a snapshot, %d of them in event mode '
          '(often trailed; Step 4 measures them)' % (len(openers),
                                                     event_openers))
    print('[summary] %d extensions have an image/event duplicate' % dups)
    missing = [o for o in obsid_dirs if o not in by_obsid]
    if missing:
        print('[note] %d OBSID folder(s) have no sky images: %s' % (
            len(missing), ' '.join(missing[:8]) + (' ...' if len(missing) > 8
                                                   else '')))
    bad_dl = sorted(o for o, s in downloads.items()
                    if s in ('failed', 'not_found', 'check'))
    if bad_dl:
        print('[warn] Step 2 reported problems for %d OBSID(s): %s -- see '
              '%s/download_report.txt' % (len(bad_dl), ' '.join(bad_dl[:8]),
                                          args.indir))
    if log is not None:
        # Every exposure HEASARC logged should have a sky image. (The master
        # catalog's per-filter totals are no test: they sometimes differ
        # with nothing missing, e.g. elapsed time in IMAGEEVENT mode.)
        have = {(r['obsid'], r['filter'], r['extname']) for r in rows}
        filters = {r['filter'] for r in rows}
        logged = [(o, f, ext, expo) for (o, f), entries in sorted(log.items())
                  if o in by_obsid and f in filters
                  for ext, expo in entries]
        missing = [e for e in logged if e[:3] not in have]
        print('[check] exposure log (swiftuvlog): %d of %d logged exposures '
              'have a sky image' % (len(logged) - len(missing), len(logged)))
        for o, f, ext, expo in missing[:20]:
            print('          missing: %s %-5s %s (%.0f s)%s' % (
                o, f, ext, expo, ', logged without an image'
                if ext == 'UNDEF' else ''))
    print('[summary] Written: %s, %s' % (inv_path, compact_path))
    if counts['unreadable']:
        print('ERROR: %d extension(s) could not be read; re-run Step 2 for '
              'those OBSIDs (see the reason column).' % counts['unreadable'],
              file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
