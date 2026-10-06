#!/usr/bin/env python3
"""
swift_uvot_positions.py -- Step 4: where the source is in each exposure,
and the source and background regions for photometry.

For every candidate exposure in Step 3's inventory:

  1. Position and shape. The brightest pixel within 8" of your --ra/--dec,
     a centroid around it, the offset of that centroid from --ra/--dec, and
     the second-moment shape of the source: axis ratio, sizes and position
     angle. Trailed exposures have axis ratios well above 1. Then the
     concentration: the net counts in the 5" source region divided by
     those in a 15" circle. A point source gives about 0.85; a long trail
     or a smeared image, which the axis ratio can miss, gives less (flagged
     below 0.70).
  2. The source region: a 5" circle at the centroid when the source is
     detected (S/N >= 10) within 3" of --ra/--dec -- this corrects the
     pointing of exposures without aspect correction -- otherwise at
     --ra/--dec.
  3. Small-scale sensitivity: whether the source's detector position falls
     on a low-sensitivity patch in each of the CALDB's LOW, MID and HIGH
     maps (the same lookup as uvotsource, which uses LOW by default).
  4. The background region: the 27.5-35" annulus (Poole et al. 2008), minus
     a circle around every source uvotdetect finds in the deepest exposure
     of that observation and filter (radius by magnitude, as HEASoft's
     uvotlc does). It must be fully and evenly exposed, as uvotsource
     requires, with at least half its area left; otherwise a 15" circle
     45-90" away is used, wherever it is fully exposed and clear.
  5. Overrides: positions and backgrounds you set in
     uvot_region_overrides.txt win over the automatic choice.

Region files (fk5, as uvotsource requires) go to
<outdir>/regions/<OBSID>/<EXTNAME>_src.reg and _bkg.reg. Every candidate
gets one status: ok, no_background (no clean background region; can be
set by hand), or failed (with the reason). Exposures that were not
candidates in Step 3 keep their Step 3 status.

Needs HEASoft (uvotdetect, quzcif) and the UVOT CALDB.

Output (in --outdir):
    uvot_positions.txt     one row per candidate exposure
    uvot_detections.txt    the sources uvotdetect found near your source
    regions/<OBSID>/       the region files

Exit status: 1 if any exposure failed.

Usage:
    swift_uvot_positions.py --ra 187.2779 --dec 2.0524
    swift_uvot_positions.py --ra 187.2779 --dec 2.0524 --nproc 16
"""

import argparse
import gzip
import io
import math
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
from astropy.wcs import WCS, FITSFixedWarning

import swift_uvot_env as uenv
from swift_uvot_runner import run_tool, run_many
from swift_uvot_tables import read_table, write_table
from swift_uvot_inventory import expmap_path

warnings.simplefilter('ignore', FITSFixedWarning)

SRC_RADIUS = 5.0            # arcsec, the calibrated photometric aperture
BKG_INNER, BKG_OUTER = 27.5, 35.0   # arcsec, Poole et al. (2008) annulus
SEARCH_RADIUS = 8.0         # arcsec around --ra/--dec for the brightest pixel
CENTROID_RADIUS = 4.0       # arcsec around the peak for the centroid
SHAPE_RADIUS = 6.0          # arcsec around the centroid for the shape
CONC_RADIUS = 15.0          # arcsec: the 5" counts are compared with this
CONC_MIN = 0.70             # 5"/15" counts below this: light outside 5"
MIN_SNR = 10.0              # the centroid is used only for a source this clear
MAX_SHIFT = 3.0             # arcsec: centroids farther away are not used
SAME_SOURCE = 5.0           # arcsec: a detection this close is your source
EXPOSED_TOLERANCE = 0.01    # uvotsource's expdeltaf: uneven exposure limit
MIN_BKG_FRACTION = 0.5      # of the annulus area left after exclusions
OFFSET_RADIUS = 15.0        # arcsec, fallback background circle
OFFSET_DISTANCES = (45.0, 60.0, 90.0)   # arcsec from the source
OFFSET_ANGLES = 24          # positions tried at each distance
DETECT_THRESHOLD = 2.5      # uvotdetect significance, as uvotlc
DETECTIONS_KEPT = 150.0     # arcsec: detections kept in uvot_detections.txt
# uvotlc's exclusion radius for a detected source, by magnitude
EXCLUSION_RADII = ((12.0, 20.0), (14.0, 15.0), (16.0, 10.0), (18.0, 7.0))
EXCLUSION_DEFAULT = 5.0

SSS_LEVELS = ('LOW', 'MID', 'HIGH')
OVERRIDES_FILE = 'uvot_region_overrides.txt'

COLUMNS = ['obsid', 'filter', 'hdu', 'extname', 'aspcorr', 'mode', 'first',
           'snr', 'offset', 'd_ra', 'd_dec', 'axis_ratio', 'sig_major',
           'sig_minor', 'major_pa', 'conc', 'center', 'src_ra', 'src_dec',
           'detx', 'dety', 'sss_low', 'sss_mid', 'sss_high', 'bkg_type',
           'bkg_area', 'n_excl', 'src_reg', 'bkg_reg', 'status', 'reason']

# DET -> RAW detector coordinates: the polynomial in HEASoft's
# UVOT::Calibration::estimateRAWfromDET (lib/perl/UVOT/Calibration.pm),
# which uvotsource uses for the small-scale-sensitivity lookup.
_KX = ((-78.649210, 0.015565185, -1.7515713e-07, 4.3219390e-10, -8.1065362e-13),
       (1.0244568, -6.4775283e-06, 2.5891338e-09, -3.8141409e-12, 9.0024773e-15),
       (2.5823715e-08, -1.7746245e-09, -6.1715464e-12, 7.1074961e-17, 2.1666564e-18),
       (9.0604263e-09, -2.9960087e-12, 1.7198618e-14, 4.0598908e-18, -7.6395990e-21),
       (-1.8713561e-12, -8.5996806e-16, 7.0224033e-19, -1.4290023e-21, -1.2595451e-25))
_KY = ((-75.469185, 0.98578166, -5.4731736e-06, 1.7460179e-08, -1.7065813e-12),
       (0.012609845, -4.0656442e-06, 9.0395420e-10, -3.6852265e-12, -6.8798446e-16),
       (4.1295693e-06, 2.8856846e-09, -5.1419799e-12, 1.8097329e-14, -1.5504342e-18),
       (-1.2398514e-09, -3.6408521e-12, -2.5093620e-15, 2.7977826e-18, 1.5104439e-21),
       (-2.8081675e-12, 1.1067681e-14, 7.2664194e-19, -5.8928808e-21, 1.6373034e-24))


def raw_from_det(detx, dety):
    """RAW detector position from DET, as uvotsource computes it."""
    x2, y2 = detx - 1023.5, dety - 1023.5
    xr = yr = 0.0
    for i in range(5):
        for j in range(5):
            xr += _KX[j][i] * x2 ** j * y2 ** i
            yr += _KY[j][i] * x2 ** j * y2 ** i
    return xr + 1023.5, yr + 1023.5


def detector_position(header, ra, dec):
    """
    (DETX, DETY, RAWX, RAWY) of a sky position in one exposure, as in
    UVOT::Source (updateDetectorPosition, applySmallScaleSensitivity): the
    sky image's 'D' WCS gives DET in mm; RAW is clamped to the detector and
    shifted by the on-board shift-and-add offset (UD_RAWX/Y, in newer
    processing only). Verified against uvotsource's DETX/DETY and
    SSS_FACTOR on 1,457 3C 273 exposures.
    """
    px, py = WCS(header).world_to_pixel_values(ra, dec)
    mmx, mmy = WCS(header, key='D').pixel_to_world_values(px, py)
    detx, dety = mmx / 0.009075 + 1100.5, mmy / 0.009075 + 1100.5
    rawx, rawy = raw_from_det(detx, dety)
    rawx = min(max(rawx, -0.5), 2047.49) - float(header.get('UD_RAWX', 0) or 0)
    rawy = min(max(rawy, -0.5), 2047.49) - float(header.get('UD_RAWY', 0) or 0)
    rawx = min(max(rawx, -0.5), 2047.49)
    rawy = min(max(rawy, -0.5), 2047.49)
    return float(detx), float(dety), rawx, rawy


def sss_map_paths():
    """{level: path} of the CALDB small-scale-sensitivity maps (quzcif)."""
    paths = {}
    for level in SSS_LEVELS:
        res = run_tool('quzcif', ['SWIFT', 'UVOTA', '-', 'V', 'SKYFLAT_SSS',
                                  'now', 'now', 'TYPE.eq.%s' % level],
                       timeout=60)
        words = res.stdout.split()
        if not res.ok or not words or not os.path.isfile(words[0]):
            raise RuntimeError('no %s small-scale-sensitivity map in the CALDB '
                               '(%s)' % (level, res.reason or res.tail))
        paths[level] = words[0]
    return paths


def exclusion_radius(mag):
    """uvotlc's background exclusion radius (arcsec) for a source magnitude."""
    for limit, radius in EXCLUSION_RADII:
        if mag <= limit:
            return radius
    return EXCLUSION_DEFAULT


def read_overrides(path):
    """
    Region overrides: lines of 'OBSID FILTER EXTNAME WHAT VALUE', with '*'
    as a wildcard in the first three fields. WHAT is 'src' (VALUE: 'RA DEC'
    in degrees, or 'catalog') or 'bkg' (VALUE: an fk5 region, e.g.
    circle(187.30,2.04,20")).
    """
    rules = []
    if not os.path.exists(path):
        return rules
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            parts = line.split(None, 4)
            if len(parts) != 5 or parts[3] not in ('src', 'bkg'):
                raise ValueError('%s:%d: expected OBSID FILTER EXTNAME src|bkg '
                                 'VALUE' % (path, lineno))
            rules.append(tuple(parts) + (lineno,))
    return rules


def matching_overrides(rules, row):
    """{what: (value, line number)} of the last rules matching a row."""
    found = {}
    for obsid, filt, ext, what, value, lineno in rules:
        if all(p == '*' or p.upper() == str(v).upper() for p, v in
               ((obsid, row['obsid']), (filt, row['filter']),
                (ext, row['extname']))):
            found[what] = (value, lineno)
    return found


def _open_fits(path):
    """Open a FITS file, reading a .gz completely first (see Step 3)."""
    if path.endswith('.gz'):
        with gzip.open(path, 'rb') as fh:
            return fits.open(io.BytesIO(fh.read()))
    return fits.open(path)


def measure_source(data, wcs, scale, x0, y0):
    """
    Peak, centroid, S/N and shape of the source near pixel (x0, y0).
    Returns a dict (pixel coordinates; sizes in arcsec), or None if nothing
    is there.
    """
    ny, nx = data.shape
    r = int(math.ceil((BKG_OUTER + 2) / scale))
    xc, yc = int(round(x0)), int(round(y0))
    xa, xb = max(0, xc - r), min(nx, xc + r + 1)
    ya, yb = max(0, yc - r), min(ny, yc + r + 1)
    if xa >= xb or ya >= yb:
        return None
    cut = np.nan_to_num(data[ya:yb, xa:xb].astype(float))
    yy, xx = np.mgrid[ya:yb, xa:xb]
    dist0 = np.hypot(xx - x0, yy - y0) * scale
    ring = (dist0 >= BKG_INNER) & (dist0 <= BKG_OUTER)
    bkg = float(np.median(cut[ring])) if ring.any() else 0.0
    search = dist0 <= SEARCH_RADIUS
    if not search.any():
        return None
    iy, ix = np.unravel_index(np.argmax(np.where(search, cut, -np.inf)),
                              cut.shape)
    cx, cy = float(xx[iy, ix]), float(yy[iy, ix])
    for _ in range(3):  # centroid around the current estimate
        m = np.hypot(xx - cx, yy - cy) * scale <= CENTROID_RADIUS
        w = np.clip(cut[m] - bkg, 0, None)
        if w.sum() <= 0:
            return None
        cx, cy = (xx[m] * w).sum() / w.sum(), (yy[m] * w).sum() / w.sum()
    m5 = np.hypot(xx - cx, yy - cy) * scale <= SRC_RADIUS
    total = cut[m5].sum()
    net = total - bkg * m5.sum()
    snr = net / math.sqrt(total) if total > 0 else 0.0
    ms = np.hypot(xx - cx, yy - cy) * scale <= SHAPE_RADIUS
    w = np.clip(cut[ms] - bkg, 0, None)
    dx, dy = xx[ms] - cx, yy[ms] - cy
    cov = np.array([[(w * dx * dx).sum(), (w * dx * dy).sum()],
                    [(w * dx * dy).sum(), (w * dy * dy).sum()]]) / w.sum()
    vals, vecs = np.linalg.eigh(cov)
    vmin, vmax = max(vals[0], 1e-9), max(vals[1], 1e-9)
    major = vecs[:, 1]
    # position angle of the major axis on the sky, east of north
    p1 = wcs.pixel_to_world(cx, cy)
    p2 = wcs.pixel_to_world(cx + major[0], cy + major[1])
    pa = p1.position_angle(p2).deg % 180.0
    return {'x': cx, 'y': cy, 'snr': snr, 'axis_ratio': math.sqrt(vmax / vmin),
            'sig_major': math.sqrt(vmax) * scale,
            'sig_minor': math.sqrt(vmin) * scale, 'major_pa': pa}


def concentration(data, scale, x, y):
    """
    (net counts in the 5" circle / net counts in the 15" circle, S/N in the
    15" circle) at pixel (x, y), with the median of the 27.5-35" annulus as
    background. A point source gives 0.83-0.87 (lower with strong
    coincidence loss); a trailed or smeared one less. None if the 15"
    circle holds no source light.
    """
    ny, nx = data.shape
    r = int(math.ceil((BKG_OUTER + 2) / scale))
    xc, yc = int(round(x)), int(round(y))
    xa, xb = max(0, xc - r), min(nx, xc + r + 1)
    ya, yb = max(0, yc - r), min(ny, yc + r + 1)
    if xa >= xb or ya >= yb:
        return None, 0.0
    cut = np.nan_to_num(data[ya:yb, xa:xb].astype(float))
    yy, xx = np.mgrid[ya:yb, xa:xb]
    dist = np.hypot(xx - x, yy - y) * scale
    ring = (dist >= BKG_INNER) & (dist <= BKG_OUTER)
    bkg = float(np.median(cut[ring])) if ring.any() else 0.0
    m5, m15 = dist <= SRC_RADIUS, dist <= CONC_RADIUS
    net5 = cut[m5].sum() - bkg * m5.sum()
    total15 = cut[m15].sum()
    net15 = total15 - bkg * m15.sum()
    if net15 <= 0 or total15 <= 0:
        return None, 0.0
    return net5 / net15, net15 / math.sqrt(total15)


def _uniform(expmap, mask, reference):
    """True if every masked pixel's exposure is within tolerance."""
    vals = expmap[mask]
    return vals.size > 0 and np.all(np.abs(vals - reference)
                                    <= EXPOSED_TOLERANCE * reference)


def choose_background(expmap, wcs, scale, cx, cy, exclusions):
    """
    (type, region lines, area in arcsec^2, n_excl, note) for the background:
    the annulus minus exclusions, else an offset circle, else none.
    exclusions: list of (x, y, radius_arcsec) in pixels of this image.
    """
    ny, nx = expmap.shape
    margin = 0.75 * scale  # uvotsource counts partly covered pixels
    reach = max(OFFSET_DISTANCES) + OFFSET_RADIUS + 2
    r = int(math.ceil(reach / scale))
    xc, yc = int(round(cx)), int(round(cy))
    xa, xb = max(0, xc - r), min(nx, xc + r + 1)
    ya, yb = max(0, yc - r), min(ny, yc + r + 1)
    emap = np.zeros((2 * r + 1, 2 * r + 1))
    emap[ya - (yc - r):yb - (yc - r), xa - (xc - r):xb - (xc - r)] = \
        np.nan_to_num(expmap[ya:yb, xa:xb])
    yy, xx = np.mgrid[yc - r:yc + r + 1, xc - r:xc + r + 1]
    dist = np.hypot(xx - cx, yy - cy) * scale
    ref = float(np.median(emap[dist <= SRC_RADIUS]))
    if ref <= 0:
        return 'none', [], 0.0, 0, 'no exposure at the source'

    def excluded(mask_xx, mask_yy, grow):
        out = np.zeros(mask_xx.shape, dtype=bool)
        for ex, ey, er in exclusions:
            out |= np.hypot(mask_xx - ex, mask_yy - ey) * scale <= er + grow
        return out

    annulus = (dist >= BKG_INNER - margin) & (dist <= BKG_OUTER + margin)
    keep = annulus & ~excluded(xx, yy, -margin)
    full_area = math.pi * (BKG_OUTER ** 2 - BKG_INNER ** 2)
    core = (dist >= BKG_INNER) & (dist <= BKG_OUTER) & ~excluded(xx, yy, 0.0)
    area = core.sum() * scale ** 2
    used = [(ex, ey, er) for ex, ey, er in exclusions
            if BKG_INNER - er <= math.hypot(ex - cx, ey - cy) * scale
            <= BKG_OUTER + er]
    note = ''
    if area >= MIN_BKG_FRACTION * full_area and _uniform(emap, keep, ref):
        return 'annulus', used, area, len(used), note
    note = ('annulus %.0f%% clear of sources' % (100 * area / full_area)
            if area < MIN_BKG_FRACTION * full_area
            else 'annulus not evenly exposed')
    for dmax in OFFSET_DISTANCES:
        for k in range(OFFSET_ANGLES):
            ang = 2 * math.pi * k / OFFSET_ANGLES
            ox = cx + dmax / scale * math.cos(ang)
            oy = cy + dmax / scale * math.sin(ang)
            if any(math.hypot(ox - ex, oy - ey) * scale <= er + OFFSET_RADIUS
                   for ex, ey, er in exclusions):
                continue
            circle = np.hypot(xx - ox, yy - oy) * scale <= OFFSET_RADIUS + margin
            if _uniform(emap, circle, ref):
                return ('circle', [(ox, oy)], math.pi * OFFSET_RADIUS ** 2, 0,
                        note + '; offset circle %.0f arcsec away' % dmax)
    return 'none', [], 0.0, 0, note + '; no clean offset circle either'


def process_file(args):
    """Measure every candidate exposure in one sky image."""
    sky_path, rows, src_ra, src_dec, detections, sss_paths, overrides, outdir \
        = args
    results = []
    try:
        sky = _open_fits(sky_path)
        ex = _open_fits(expmap_path(sky_path))
        ex_by_name = {h.header.get('EXTNAME'): h for h in ex[1:]}
        sss = {lvl: fits.open(p, memmap=True) for lvl, p in sss_paths.items()}
    except Exception as exc:
        return [dict(r, status='failed', reason='cannot read: %s' % exc)
                for r in rows]
    for row in rows:
        out = {k: row.get(k) for k in ('obsid', 'filter', 'hdu', 'extname',
                                       'aspcorr', 'mode', 'first')}
        reasons = []
        try:
            hdu = sky[int(row['hdu'])]
            h, data = hdu.header, hdu.data
            expmap = ex_by_name[row['extname']].data
            wcs = WCS(h)
            scale = abs(h['CDELT1']) * 3600.0
            x0, y0 = [float(v) for v in wcs.world_to_pixel_values(src_ra, src_dec)]
            meas = measure_source(data, wcs, scale, x0, y0)
            center_ra, center_dec, center = src_ra, src_dec, 'catalog'
            if meas is None:
                reasons.append('no source signal near the position')
            else:
                c = wcs.pixel_to_world(meas['x'], meas['y'])
                cat = SkyCoord(src_ra, src_dec, unit='deg')
                offset = c.separation(cat).arcsec
                out.update({
                    'snr': '%.1f' % meas['snr'],
                    'offset': '%.2f' % offset,
                    'd_ra': '%+.2f' % ((c.ra.deg - src_ra) * 3600.0
                                       * math.cos(math.radians(src_dec))),
                    'd_dec': '%+.2f' % ((c.dec.deg - src_dec) * 3600.0),
                    'axis_ratio': '%.2f' % meas['axis_ratio'],
                    'sig_major': '%.2f' % meas['sig_major'],
                    'sig_minor': '%.2f' % meas['sig_minor'],
                    'major_pa': '%.0f' % meas['major_pa'],
                })
                if meas['snr'] < MIN_SNR:
                    reasons.append('S/N %.1f: region at --ra/--dec'
                                   % meas['snr'])
                elif offset > MAX_SHIFT:
                    reasons.append('centroid %.1f arcsec away: region at '
                                   '--ra/--dec' % offset)
                else:
                    center_ra, center_dec, center = c.ra.deg, c.dec.deg, \
                        'centroid'
            ovr = matching_overrides(overrides, row)
            if 'src' in ovr:
                value, lineno = ovr['src']
                if value.lower() != 'catalog':
                    center_ra, center_dec = [float(v) for v in value.split()]
                else:
                    center_ra, center_dec = src_ra, src_dec
                center = 'override:%d' % lineno
            out.update({'center': center, 'src_ra': '%.6f' % center_ra,
                        'src_dec': '%.6f' % center_dec})
            cx, cy = [float(v) for v in wcs.world_to_pixel_values(center_ra,
                                                                  center_dec)]
            # How much of the source light near the region is inside it:
            # trailed and smeared images lose light the axis ratio can miss.
            conc, snr15 = concentration(data, scale, cx, cy)
            if conc is not None and snr15 >= MIN_SNR:
                out['conc'] = '%.2f' % conc
                if conc < CONC_MIN:
                    reasons.append('only %.0f%% of the counts within 15 arcsec '
                                   'are in the 5 arcsec circle: trailed or '
                                   'smeared' % (100 * conc))

            detx, dety, rawx, rawy = detector_position(h, center_ra, center_dec)
            out.update({'detx': '%.1f' % detx, 'dety': '%.1f' % dety})
            for lvl in SSS_LEVELS:
                val = sss[lvl]['SSSENS' + row['filter']].data[
                    int(rawy + 0.5), int(rawx + 0.5)]
                out['sss_%s' % lvl.lower()] = '%.1f' % val

            # background
            excl = []
            dets = detections or []
            centre = SkyCoord(center_ra, center_dec, unit='deg')
            catpos = SkyCoord(src_ra, src_dec, unit='deg')
            # Your source is the brightest detection within 10" of --ra/--dec:
            # in a trailed or mis-pointed exposure it can be several arcsec
            # from the region centre, and must not be masked as a neighbour.
            near = [(mag, i) for i, (ra_d, dec_d, mag) in enumerate(dets)
                    if SkyCoord(ra_d, dec_d, unit='deg').separation(catpos)
                    .arcsec <= SEARCH_RADIUS + 2]
            target = min(near)[1] if near else None
            for i, (ra_d, dec_d, mag) in enumerate(dets):
                if i == target or SkyCoord(ra_d, dec_d, unit='deg').separation(
                        centre).arcsec <= SAME_SOURCE:
                    continue
                ex_x, ex_y = [float(v) for v in wcs.world_to_pixel_values(ra_d,
                                                                          dec_d)]
                excl.append((ex_x, ex_y, exclusion_radius(mag)))
            if detections is None:
                reasons.append('no source list (uvotdetect failed): no '
                               'exclusions')
            btype, bgeom, area, n_excl, note = choose_background(
                expmap, wcs, scale, cx, cy, excl)
            if note:
                reasons.append(note)

            reg_dir = os.path.join(outdir, 'regions', row['obsid'])
            os.makedirs(reg_dir, exist_ok=True)
            src_reg = os.path.join(reg_dir, row['extname'] + '_src.reg')
            bkg_reg = os.path.join(reg_dir, row['extname'] + '_bkg.reg')
            with open(src_reg, 'w') as fh:
                fh.write('fk5;circle(%.6f,%.6f,%.1f")\n'
                         % (center_ra, center_dec, SRC_RADIUS))
            lines = []
            if 'bkg' in ovr:
                value, lineno = ovr['bkg']
                lines = ['fk5;' + value]
                btype, area, n_excl = 'override:%d' % lineno, float('nan'), 0
            elif btype == 'annulus':
                lines = ['fk5;annulus(%.6f,%.6f,%.1f",%.1f")'
                         % (center_ra, center_dec, BKG_INNER, BKG_OUTER)]
                for ex_x, ex_y, er in bgeom:
                    p = wcs.pixel_to_world(ex_x, ex_y)
                    lines.append('fk5;-circle(%.6f,%.6f,%.1f")'
                                 % (p.ra.deg, p.dec.deg, er))
            elif btype == 'circle':
                ox, oy = bgeom[0]
                p = wcs.pixel_to_world(ox, oy)
                lines = ['fk5;circle(%.6f,%.6f,%.1f")'
                         % (p.ra.deg, p.dec.deg, OFFSET_RADIUS)]
            if lines:
                with open(bkg_reg, 'w') as fh:
                    fh.write('\n'.join(lines) + '\n')
            elif os.path.exists(bkg_reg):
                os.remove(bkg_reg)
            out.update({'bkg_type': btype,
                        'bkg_area': '-' if area != area else '%.0f' % area,
                        'n_excl': n_excl,
                        'src_reg': os.path.relpath(src_reg, outdir),
                        'bkg_reg': os.path.relpath(bkg_reg, outdir)
                        if lines else '-'})
            out['status'] = 'ok' if lines else 'no_background'
        except Exception as exc:
            out['status'] = 'failed'
            reasons.append('%s: %s' % (exc.__class__.__name__, exc))
        out['reason'] = '; '.join(reasons)
        results.append(out)
    return results


def run_detections(groups, outdir, nproc):
    """
    uvotdetect on the deepest candidate exposure of each observation and
    filter (aspect-corrected ones preferred). Returns ({(obsid, filter):
    [(ra, dec, mag), ...] or None}, {(obsid, filter): reason}).
    """
    calls, keys = [], []
    det_dir = os.path.join(outdir, 'detections')
    os.makedirs(det_dir, exist_ok=True)
    for key, rows in sorted(groups.items()):
        best = max(rows, key=lambda r: (r['aspcorr'] == 'DIRECT',
                                        float(r['exposure'])))
        sky = best['file']
        ex = expmap_path(sky)
        dest = os.path.join(det_dir, '%s_%s_%s.fits' % (key[0], key[1],
                                                         best['extname']))
        calls.append(dict(
            tool='uvotdetect',
            params={'infile': 'sky.img.gz[%s]' % best['hdu'],
                    'outfile': 'det.fits',
                    'expfile': 'ex.img.gz[%s]' % best['hdu'],
                    'threshold': DETECT_THRESHOLD, 'calibrate': 'yes',
                    'clobber': 'yes', 'chatter': 1},
            inputs={'sky.img.gz': sky, 'ex.img.gz': ex},
            outputs={'det.fits': dest}, timeout=600,
            log=dest.replace('.fits', '.log')))
        keys.append((key, best['extname'], dest))
    results = run_many(calls, nproc=nproc)
    detections, problems = {}, {}
    for (key, extname, dest), res in zip(keys, results):
        if not res.ok:
            detections[key] = None
            problems[key] = res.reason
            continue
        tab = fits.getdata(dest, 'SOURCES')
        detections[key] = [(float(r['RA']), float(r['DEC']), float(r['MAG']))
                           for r in tab]
    return detections, problems


def main(argv=None):
    sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser(
        description='Step 4: source positions, shapes and regions.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--ra', type=float, required=True,
                        help='source RA, decimal degrees (J2000)')
    parser.add_argument('--dec', type=float, required=True,
                        help='source Dec, decimal degrees (J2000)')
    parser.add_argument('--outdir', default='UVOT_output',
                        help='Step 3 output folder (default: UVOT_output)')
    parser.add_argument('--nproc', type=int, default=8,
                        help='parallel workers (default: 8)')
    args = parser.parse_args(argv)

    uenv.require_heasoft_shell(['uvotdetect', 'quzcif'])
    inv_path = os.path.join(args.outdir, 'uvot_inventory.txt')
    if not os.path.exists(inv_path):
        print('ERROR: %s not found; run Step 3 (swift_uvot_inventory.py) '
              'first.' % inv_path, file=sys.stderr)
        return 1
    inventory = read_table(inv_path)
    with open(inv_path) as fh:
        header = ''.join(line for line in fh if line.startswith('#'))
    m = re.search(r'source: RA ([\d.+-]+) Dec ([\d.+-]+)', header)
    if m and SkyCoord(float(m.group(1)), float(m.group(2)), unit='deg') \
            .separation(SkyCoord(args.ra, args.dec, unit='deg')).arcsec > 1.0:
        print('ERROR: --ra/--dec differ from the inventory\'s source (RA %s '
              'Dec %s); use the same position in every step.'
              % (m.group(1), m.group(2)), file=sys.stderr)
        return 1
    candidates = [r for r in inventory if r['status'] == 'candidate']
    print('[info] %d candidate exposures of %d in the inventory'
          % (len(candidates), len(inventory)))

    try:
        sss_paths = sss_map_paths()
    except RuntimeError as exc:
        print('ERROR: %s' % exc, file=sys.stderr)
        return 1
    try:
        overrides = read_overrides(os.path.join(args.outdir, OVERRIDES_FILE))
    except ValueError as exc:
        print('ERROR: %s' % exc, file=sys.stderr)
        return 1
    if overrides:
        print('[info] %d region override(s) from %s'
              % (len(overrides), OVERRIDES_FILE))

    groups = defaultdict(list)
    for r in candidates:
        groups[(r['obsid'], r['filter'])].append(r)
    print('[info] uvotdetect on the deepest exposure of %d observation/filter '
          'pairs ...' % len(groups))
    detections, det_problems = run_detections(groups, args.outdir, args.nproc)

    by_file = defaultdict(list)
    for r in candidates:
        by_file[r['file']].append(r)
    jobs = [(path, rows, args.ra, args.dec,
             detections.get((rows[0]['obsid'], rows[0]['filter'])),
             sss_paths, overrides, args.outdir)
            for path, rows in sorted(by_file.items())]
    print('[info] measuring %d exposures in %d sky images ...'
          % (len(candidates), len(jobs)))
    with ProcessPoolExecutor(max_workers=max(1, args.nproc)) as pool:
        rows = [r for part in pool.map(process_file, jobs) for r in part]
    order = {(r['obsid'], r['filter'], r['extname']): i
             for i, r in enumerate(candidates)}
    # the ledger: every candidate must come back exactly once
    if sorted(order) != sorted((r['obsid'], r['filter'], r['extname'])
                               for r in rows):
        print('ERROR: internal: %d candidates in, %d rows out'
              % (len(candidates), len(rows)), file=sys.stderr)
        return 1
    rows.sort(key=lambda r: order[(r['obsid'], r['filter'], r['extname'])])

    # sources near the target, for the record
    src = SkyCoord(args.ra, args.dec, unit='deg')
    det_rows = []
    for (obsid, filt), dets in sorted(detections.items()):
        for ra_d, dec_d, mag in dets or []:
            sep = SkyCoord(ra_d, dec_d, unit='deg').separation(src).arcsec
            if sep <= DETECTIONS_KEPT:
                det_rows.append({'obsid': obsid, 'filter': filt,
                                 'ra': '%.6f' % ra_d, 'dec': '%.6f' % dec_d,
                                 'sep': '%.1f' % sep, 'mag': '%.2f' % mag,
                                 'excl_radius': exclusion_radius(mag)})

    counts = Counter(r['status'] for r in rows)
    prov = uenv.provenance()
    comments = [
        'swift_uvot_positions.py, %s' % datetime.now(timezone.utc)
        .strftime('%Y-%m-%dT%H:%M:%SZ'),
        'command: %s' % ' '.join(sys.argv),
        'source: RA %.6f Dec %.6f' % (args.ra, args.dec),
        'pipeline %s; HEASoft %s; CALDB %s (UVOT %s)' % (
            prov['PIPECOMM'], prov['HEASOFT'], prov['CALDB'], prov['UVOTINDX']),
        '%d candidate exposures: %s' % (len(rows), ', '.join(
            '%s %d' % (k, v) for k, v in sorted(counts.items()))),
        'offset, d_ra, d_dec: centroid - (--ra, --dec), arcsec; axis_ratio, '
        'sig_major/minor (arcsec), major_pa (deg E of N): source shape; '
        'conc: net counts in the 5" region / in 15" (point source ~0.85)',
        'center: where the 5" source region is (centroid, catalog or an '
        'override line); sss_*: CALDB small-scale-sensitivity factor at the '
        'source (-99.9 = low-sensitivity patch)',
        'bkg_type: annulus (27.5-35" minus n_excl source circles), circle '
        '(15" offset circle), override, or none; bkg_area in arcsec^2',
    ]
    pos_path = os.path.join(args.outdir, 'uvot_positions.txt')
    write_table(pos_path, rows, COLUMNS, comments)
    det_path = os.path.join(args.outdir, 'uvot_detections.txt')
    write_table(det_path, det_rows, ['obsid', 'filter', 'ra', 'dec', 'sep',
                                     'mag', 'excl_radius'],
                comments[:4] + ['uvotdetect sources within %.0f" of the source '
                                'in the deepest exposure of each observation '
                                'and filter' % DETECTIONS_KEPT])

    # --- report -------------------------------------------------------------
    print()
    print('[summary] %d candidate exposures: %s' % (len(rows), ', '.join(
        '%s %d' % (k, v) for k, v in sorted(counts.items()))))
    print('[summary] source region at: %s' % ', '.join(
        '%s %d' % (k, v) for k, v in sorted(Counter(
            r.get('center', '-').split(':')[0] for r in rows).items())))
    print('[summary] background: %s' % ', '.join(
        '%s %d' % (k, v) for k, v in sorted(Counter(
            str(r.get('bkg_type', '-')).split(':')[0] for r in rows).items())))
    for asp in ('DIRECT', 'NONE'):
        sel = [float(r['offset']) for r in rows
               if r.get('offset') not in (None, '-') and r['aspcorr'] == asp]
        if sel:
            print('          offsets, ASPCORR=%-6s: median %.2f", 90%% %.2f", '
                  '%d over %.0f"' % (asp, np.median(sel),
                                     np.percentile(sel, 90),
                                     sum(s > MAX_SHIFT for s in sel), MAX_SHIFT))
    trail = [r for r in rows if r.get('axis_ratio') not in (None, '-')
             and float(r['axis_ratio']) > 1.3]
    print('[summary] axis ratio > 1.3 (trailed or doubled?): %d, of which %d '
          'event-mode snapshot openers' % (len(trail), sum(
              1 for r in trail if r['mode'] == 'EVENT' and r['first'] == 'yes')))
    lost = [r for r in rows if r.get('conc') not in (None, '-')
            and float(r['conc']) < CONC_MIN]
    print('[summary] light outside the 5" region (5"/15" counts < %.2f; '
          'trailed or smeared): %d, of which %d event-mode snapshot openers'
          % (CONC_MIN, len(lost), sum(1 for r in lost if r['mode'] == 'EVENT'
                                      and r['first'] == 'yes')))
    faint = [r for r in rows if r['status'] != 'failed' and (
        r.get('snr') in (None, '-') or float(r['snr']) < MIN_SNR)]
    print('[summary] S/N < %.0f at the source: %d (faint, or its light is '
          'elsewhere: see the viewer)' % (MIN_SNR, len(faint)))
    for lvl in SSS_LEVELS:
        bad = Counter(r['filter'] for r in rows
                      if r.get('sss_%s' % lvl.lower()) not in (None, '-')
                      and float(r['sss_%s' % lvl.lower()]) <= 0)
        print('[summary] on a %-4s low-sensitivity patch: %d (%s)' % (
            lvl, sum(bad.values()), ', '.join('%s %d' % kv for kv in
                                              sorted(bad.items()))))
    if det_problems:
        print('[warn] uvotdetect failed for %d observation/filter pair(s): '
              'their backgrounds have no source exclusions' % len(det_problems))
        for key, why in list(det_problems.items())[:5]:
            print('          %s %s: %s' % (key[0], key[1], why))
    print('[summary] Written: %s, %s, %s/regions/' % (pos_path, det_path,
                                                     args.outdir))
    if counts['failed']:
        print('ERROR: %d exposure(s) failed; see the reason column.'
              % counts['failed'], file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
