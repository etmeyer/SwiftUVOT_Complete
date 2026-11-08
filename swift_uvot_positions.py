#!/usr/bin/env python3
"""
swift_uvot_positions.py -- Step 4: where the source is in each exposure,
and the source and background regions for photometry.

For every candidate exposure in Step 3's inventory:

  1. The background (item 4 below), then position and shape: the
     brightest pixel within 8" of your --ra/--dec, a centroid around it,
     its S/N, and -- only at S/N >= 10, where noise no longer dominates --
     the offset of that centroid from --ra/--dec and the second-moment
     shape: axis ratio, sizes and position angle. Trailed exposures have
     axis ratios well above 1. Then the concentration: the net counts in
     the 5" source region divided by those in a 15" circle, with its error
     (counts, plus 1 % of the background). A point source gives
     0.86-0.92; a long trail or a smeared image, which the axis ratio can
     miss, gives less (flagged below 0.70). It is reported only where its
     error is at most 0.02; a faint source gets '-' and the reason.
  2. The source region: a 5" circle at the centroid when the source is
     detected (S/N >= 10) within 3" of --ra/--dec -- this corrects the
     pointing of exposures without aspect correction -- otherwise at
     --ra/--dec.
  3. Small-scale sensitivity: whether the source's detector position falls
     on a low-sensitivity patch in each of the CALDB's LOW, MID and HIGH
     maps (the same lookup as uvotsource, which uses LOW by default), and
     what fraction of the 5" circle is on a LOW patch.
  4. The background region: the 27.5-35" annulus (Poole et al. 2008), minus
     a circle around every source uvotdetect finds in the deepest exposure
     of that observation and filter. Each circle reaches out to where the
     source's PSF wing (CALDB profile) falls below 0.5 % of the sky, 5-40";
     a source at 0.9 counts per frame or more, whose measured rate
     understates its wings, gets 40" and a 2" strip either side of its
     read-out direction. A faint detection within 45" of a source 100
     times brighter is most likely a fragment of that source's halo and
     keeps uvotlc's magnitude-binned radius (5-20"). Then the annulus is cut into 8 sectors, and a
     sector that is too bright against the others (Poisson p < 0.001) is
     dropped whole: never single pixels, which would bias a sky of a few
     counts per pixel low. The background is the plain mean of what is
     left -- uvotsource's own estimator. The annulus must be fully and
     evenly exposed, as uvotsource requires, with at least half its area
     left; otherwise a 15" circle 45-90" away is used, wherever it is fully
     exposed and clear of masks.
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
from scipy.stats import poisson

import swift_uvot_env as uenv
import swift_uvot_fieldstars as fs
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
CONC_ERR_MAX = 0.02         # concentration reported only if this precise
SKY_SYSTEMATIC = 0.01       # of the background, in the concentration error
MIN_SNR = 10.0              # centroid and shape only for a source this clear
MAX_SHIFT = 3.0             # arcsec: centroids farther away are not used
SAME_SOURCE = 5.0           # arcsec: a detection this close is your source
EXPOSED_TOLERANCE = 0.01    # uvotsource's expdeltaf: uneven exposure limit
MIN_BKG_FRACTION = 0.5      # of the annulus area left after exclusions
OFFSET_RADIUS = 15.0        # arcsec, fallback background circle
OFFSET_DISTANCES = (45.0, 60.0, 90.0)   # arcsec from the source
OFFSET_ANGLES = 24          # positions tried at each distance
DETECT_THRESHOLD = 2.5      # uvotdetect significance, as uvotlc
DETECTIONS_KEPT = 150.0     # arcsec: detections kept in uvot_detections.txt
# uvotlc's exclusion radius for a detected source, by magnitude: used only
# for the first, rough sky level that sizes the real masks
EXCLUSION_RADII = ((12.0, 20.0), (14.0, 15.0), (16.0, 10.0), (18.0, 7.0))
EXCLUSION_DEFAULT = 5.0
# Masks sized by flux: out to where a neighbour's PSF wing falls below this
# fraction of the sky surface brightness (CALDB REEF profile, exponential
# tail beyond its 30" end), between MASK_MIN and MASK_MAX.
WING_FRACTION = 0.005
MASK_MIN, MASK_MAX = 5.0, 40.0      # arcsec
TAIL_FIT = (20.0, 29.5)     # arcsec: the CALDB profile range the tail is fit on
SKY_RING = (35.0, 100.0)    # arcsec around the source: the sky for the masks
# A source at or above this many raw counts per frame is saturated or nearly
# so: its measured rate understates its wings. It gets MASK_MAX, plus a strip
# along the read-out direction (detector Y).
BRIGHT_CPF = 0.9
# A detection within HALO_REACH of a source HALO_RATIO times brighter (or of
# a bright one) is most likely a fragment of that source's halo or smoke
# ring, which uvotdetect breaks into many faint 'sources' around a
# saturated star: it keeps uvotlc's radius, or it would eat the annulus.
HALO_RATIO = 100.0
HALO_REACH = 45.0           # arcsec
STRIP_HALF_WIDTH = 2.0      # arcsec
# Background sectors: the annulus is cut into N_SECTORS; a sector whose
# counts are high at p < SECTOR_P (Poisson, against the other sectors) is
# dropped whole -- never single pixels, which would bias a sky of a few
# counts per pixel low.
N_SECTORS = 8
SECTOR_P = 1e-3
MIN_SECTORS_LEFT = 4
SSS_EDGE_FRACTION = 0.10    # of the 5" circle on a LOW patch: noted

SSS_LEVELS = ('LOW', 'MID', 'HIGH')
OVERRIDES_FILE = 'uvot_region_overrides.txt'

COLUMNS = ['obsid', 'filter', 'hdu', 'extname', 'aspcorr', 'mode', 'first',
           'snr', 'offset', 'd_ra', 'd_dec', 'axis_ratio', 'sig_major',
           'sig_minor', 'major_pa', 'conc', 'conc_err', 'center', 'src_ra',
           'src_dec', 'detx', 'dety', 'sss_low', 'sss_mid', 'sss_high',
           'sss_low_frac', 'bkg_type', 'bkg_area', 'n_excl', 'n_strip',
           'n_sect', 'bkg5', 'sky', 'pnt_dra', 'pnt_ddec', 'pnt_n', 'pnt_err',
           'pnt_src',
           'n_stars', 'star_conc', 'star_conc_err', 'star_ar', 'star_sigma',
           'psf_loss', 'psf_loss_err', 'psf_src', 'src_reg', 'bkg_reg',
           'status', 'reason']

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


def wing_profiles():
    """
    {filter: (radius grid, surface brightness)}: the PSF wing per arcsec^2,
    as a fraction of a source's 5" counts, from the CALDB encircled-energy
    file (REEF, normalised to 1 at 5"). Beyond the file's 30" the wing is an
    exponential fit to its last 10". Each value is the largest at that
    radius or beyond, so a radius whose value is below a level stays below
    it farther out.
    """
    res = run_tool('quzcif', ['SWIFT', 'UVOTA', '-', 'V', 'REEF', 'now', 'now',
                              '-'], timeout=60)
    words = res.stdout.split()
    if not res.ok or not words or not os.path.isfile(words[0]):
        raise RuntimeError('no PSF encircled-energy file (REEF) in the CALDB '
                           '(%s)' % (res.reason or res.tail))
    grid = np.arange(0.25, MASK_MAX + 0.01, 0.25)
    profiles = {}
    with fits.open(words[0]) as hdul:
        for hdu in hdul[1:]:
            name = str(hdu.header.get('EXTNAME', ''))
            if not name.startswith('REEF'):
                continue
            r = np.asarray(hdu.data['RADIUS'], float)
            ee = np.asarray(hdu.data['REEF'], float)
            mid = 0.5 * (r[1:] + r[:-1])
            sb = np.diff(ee) / (np.pi * np.diff(r ** 2))
            fit = (mid >= TAIL_FIT[0]) & (mid <= TAIL_FIT[1]) & (sb > 0)
            slope, icpt = np.polyfit(mid[fit], np.log(sb[fit]), 1)
            prof = np.where(grid <= mid[-1], np.interp(grid, mid, sb),
                            np.exp(icpt + slope * grid))
            prof = np.maximum.accumulate(np.clip(prof, 0, None)[::-1])[::-1]
            profiles[name[4:]] = (grid, prof)
    return profiles


def mask_radius(rate, sky, profile):
    """
    Mask radius (arcsec) for a neighbour of `rate` counts/s: where its wing
    falls below WING_FRACTION of the sky (counts/s/arcsec^2).
    """
    if not sky or sky <= 0:
        return MASK_MAX
    if not rate or rate != rate or rate <= 0:
        return MASK_MIN
    grid, prof = profile
    over = np.nonzero(rate * prof > WING_FRACTION * sky)[0]
    if not over.size:
        return MASK_MIN
    return float(min(MASK_MAX, max(MASK_MIN, grid[over[-1]] + 0.25)))


def _cutout(arrays, x, y, reach, scale):
    """
    Square cutouts of the arrays around pixel (x, y), `reach` arcsec on each
    side, with zeros off the image; and their pixel coordinates (xx, yy).
    """
    r = int(math.ceil(reach / scale))
    xc, yc = int(round(x)), int(round(y))
    ny, nx = arrays[0].shape
    xa, xb = max(0, xc - r), min(nx, xc + r + 1)
    ya, yb = max(0, yc - r), min(ny, yc + r + 1)
    out = []
    for a in arrays:
        c = np.zeros((2 * r + 1, 2 * r + 1))
        if xa < xb and ya < yb:
            c[ya - (yc - r):yb - (yc - r), xa - (xc - r):xb - (xc - r)] = \
                np.nan_to_num(a[ya:yb, xa:xb].astype(float))
        out.append(c)
    yy, xx = np.mgrid[yc - r:yc + r + 1, xc - r:xc + r + 1]
    return out, xx, yy


def readout_direction(header, x, y):
    """Unit pixel vector along detector Y (the read-out direction) at (x, y)."""
    wd = WCS(header, key='D')
    m0 = np.array(wd.pixel_to_world_values(x, y))
    jac = np.array([np.array(wd.pixel_to_world_values(x + 1, y)) - m0,
                    np.array(wd.pixel_to_world_values(x, y + 1)) - m0]).T
    v = np.linalg.solve(jac, [0.0, 1.0])
    return v / math.hypot(*v)


def excluded(xx, yy, scale, masks, strips, grow=0.0):
    """Pixels inside any mask circle (x, y, r) or strip (x, y, vx, vy, half)."""
    out = np.zeros(xx.shape, dtype=bool)
    for ex, ey, er in masks:
        out |= np.hypot(xx - ex, yy - ey) * scale <= er + grow
    for sx, sy, vx, vy, half in strips:
        dx, dy = xx - sx, yy - sy
        along = (dx * vx + dy * vy) * scale
        perp = (dy * vx - dx * vy) * scale
        out |= (np.abs(perp) <= STRIP_HALF_WIDTH + grow) & \
            (np.abs(along) <= half)
    return out


def sky_level(data, expmap, scale, x, y, masks):
    """
    Sky surface brightness (counts/s/arcsec^2) in the SKY_RING around (x, y),
    with the masks out and only fully exposed pixels: the level the mask
    radii are set against. None if the ring has no such pixel.
    """
    (cut, emap), xx, yy = _cutout((data, expmap), x, y, SKY_RING[1] + 2,
                                  scale)
    dist = np.hypot(xx - x, yy - y) * scale
    ring = (dist >= SKY_RING[0]) & (dist <= SKY_RING[1]) & (emap > 0)
    if not ring.any():
        return None
    ring &= emap >= 0.99 * np.median(emap[ring])
    ring &= ~excluded(xx, yy, scale, masks, [])
    if not ring.any():
        return None
    return float(cut[ring].sum() / (emap[ring].sum() * scale ** 2))


def _uniform(expmap, mask, reference):
    """True if every masked pixel's exposure is within tolerance."""
    vals = expmap[mask]
    return vals.size > 0 and np.all(np.abs(vals - reference)
                                    <= EXPOSED_TOLERANCE * reference)


def _poisson_high(n, mu):
    """P(N >= n) for a Poisson mean mu."""
    if mu <= 0:
        return 0.0 if n > 0 else 1.0
    return float(poisson.sf(round(n) - 1, mu))


def background(data, expmap, scale, cx, cy, masks, strips):
    """
    The background around pixel (cx, cy): the 27.5-35" annulus minus the
    masks and strips, minus whole sectors that are too bright; else an
    offset circle; else none. Returns a dict: type ('annulus', 'circle' or
    'none'), masks and strips that touch the region, rejected sectors, the
    offset circle centre, area (arcsec^2), mean counts per pixel and pixel
    count of the region (for 'none': of whatever annulus is left, for
    Step 4's own measurements), the exposure at the source, and a note.
    """
    reach = max(OFFSET_DISTANCES) + OFFSET_RADIUS + 2
    (cut, emap), xx, yy = _cutout((data, expmap), cx, cy, reach, scale)
    margin = 0.75 * scale  # uvotsource counts partly covered pixels
    dist = np.hypot(xx - cx, yy - cy) * scale
    ref = float(np.median(emap[dist <= SRC_RADIUS]))
    out = {'type': 'none', 'masks': [], 'strips': [], 'sectors': [],
           'offset': None, 'area': 0.0, 'mean': 0.0, 'npix': 0, 'ref': ref,
           'note': ''}
    if ref <= 0:
        out['note'] = 'no exposure at the source'
        return out
    ring = (dist >= BKG_INNER) & (dist <= BKG_OUTER)
    core = ring & ~excluded(xx, yy, scale, masks, strips) & (emap > 0)
    sector = ((np.degrees(np.arctan2(yy - cy, xx - cx)) % 360.0)
              // (360.0 / N_SECTORS)).astype(int)
    rejected = []
    while N_SECTORS - len(rejected) > MIN_SECTORS_LEFT:
        live = [k for k in range(N_SECTORS) if k not in rejected]
        worst = None
        for k in live:
            mine = core & (sector == k)
            rest = core & np.isin(sector, [j for j in live if j != k])
            if not mine.any() or not rest.any():
                continue
            p = _poisson_high(cut[mine].sum(),
                              cut[rest].sum() / rest.sum() * mine.sum())
            if worst is None or p < worst[0]:
                worst = (p, k)
        if worst is None or worst[0] >= SECTOR_P:
            break
        rejected.append(worst[1])
    keep = core & ~np.isin(sector, rejected)
    out.update(sectors=sorted(rejected), npix=int(keep.sum()),
               mean=float(cut[keep].mean()) if keep.any() else 0.0)
    full_area = math.pi * (BKG_OUTER ** 2 - BKG_INNER ** 2)
    area = keep.sum() * scale ** 2
    wide = ((dist >= BKG_INNER - margin) & (dist <= BKG_OUTER + margin)
            & ~excluded(xx, yy, scale, masks, strips, -margin)
            & ~np.isin(sector, rejected))
    if area >= MIN_BKG_FRACTION * full_area and _uniform(emap, wide, ref):
        out.update(type='annulus', area=area, masks=[
            m for m in masks if BKG_INNER - m[2] <= math.hypot(
                m[0] - cx, m[1] - cy) * scale <= BKG_OUTER + m[2]],
            strips=[s for s in strips
                    if (ring & excluded(xx, yy, scale, [], [s])).any()])
        return out
    out['note'] = ('annulus %.0f%% clear of sources' % (100 * area / full_area)
                   if area < MIN_BKG_FRACTION * full_area
                   else 'annulus not evenly exposed')
    blocked = excluded(xx, yy, scale, masks, strips)
    for dmax in OFFSET_DISTANCES:
        for k in range(OFFSET_ANGLES):
            ang = 2 * math.pi * k / OFFSET_ANGLES
            ox = cx + dmax / scale * math.cos(ang)
            oy = cy + dmax / scale * math.sin(ang)
            d_o = np.hypot(xx - ox, yy - oy) * scale
            circle = d_o <= OFFSET_RADIUS
            if (circle & blocked).any():
                continue
            if _uniform(emap, d_o <= OFFSET_RADIUS + margin, ref):
                out.update(type='circle', offset=(ox, oy), sectors=[],
                           area=math.pi * OFFSET_RADIUS ** 2,
                           mean=float(cut[circle].mean()),
                           npix=int(circle.sum()),
                           note=out['note'] + '; offset circle %.0f arcsec '
                           'away' % dmax)
                return out
    out['note'] += '; no clean offset circle either'
    return out


def measure_source(data, wcs, scale, x0, y0, bkg):
    """
    Peak, centroid, S/N and shape of the source near pixel (x0, y0), over a
    background of `bkg` counts per pixel. Returns a dict (pixel coordinates;
    sizes in arcsec), or None if nothing is there.
    """
    (cut,), xx, yy = _cutout((data,), x0, y0,
                             SEARCH_RADIUS + SHAPE_RADIUS + 2, scale)
    dist0 = np.hypot(xx - x0, yy - y0) * scale
    search = dist0 <= SEARCH_RADIUS
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


def concentration(data, scale, x, y, bkg, bkg_var):
    """
    (net counts in the 5" circle / net counts in the 15" circle, its error)
    at pixel (x, y), over a background of `bkg` counts per pixel with
    variance `bkg_var`. A point source gives 0.86-0.92 (lower with strong
    coincidence loss); a trailed or smeared one less. (None, None) if the
    15" circle holds no source light.
    """
    (cut,), xx, yy = _cutout((data,), x, y, CONC_RADIUS + 2, scale)
    dist = np.hypot(xx - x, yy - y) * scale
    m5, m15 = dist <= SRC_RADIUS, dist <= CONC_RADIUS
    t5, t15 = cut[m5].sum(), cut[m15].sum()
    a5, a15 = m5.sum(), m15.sum()
    n5 = t5 - bkg * a5
    ring = (t15 - t5) - bkg * (a15 - a5)   # net counts at 5-15"
    n15 = n5 + ring
    if n15 <= 0:
        return None, None
    v5 = max(t5, 0) + a5 ** 2 * bkg_var
    vr = max(t15 - t5, 0) + (a15 - a5) ** 2 * bkg_var
    cov = a5 * (a15 - a5) * bkg_var       # both share the background
    var = (ring ** 2 * v5 + n5 ** 2 * vr - 2 * n5 * ring * cov) / n15 ** 4
    return n5 / n15, math.sqrt(max(var, 0.0))


def sss_low_fraction(header, scale, cx, cy, sss_map):
    """
    Fraction of the 5" source circle at pixel (cx, cy) on a low-sensitivity
    patch of `sss_map` (sampled every 0.5"), with the same DET -> RAW
    mapping as detector_position.
    """
    g = np.arange(-SRC_RADIUS, SRC_RADIUS + 1e-9, 0.5)
    gx, gy = np.meshgrid(g, g)
    inside = np.hypot(gx, gy) <= SRC_RADIUS
    px, py = cx + gx[inside] / scale, cy + gy[inside] / scale
    mmx, mmy = WCS(header, key='D').pixel_to_world_values(px, py)
    rawx, rawy = raw_from_det(mmx / 0.009075 + 1100.5, mmy / 0.009075 + 1100.5)
    rawx = np.clip(np.clip(rawx, -0.5, 2047.49)
                   - float(header.get('UD_RAWX', 0) or 0), -0.5, 2047.49)
    rawy = np.clip(np.clip(rawy, -0.5, 2047.49)
                   - float(header.get('UD_RAWY', 0) or 0), -0.5, 2047.49)
    vals = sss_map[(rawy + 0.5).astype(int), (rawx + 0.5).astype(int)]
    return float(np.mean(vals <= 0))


def _polygon(wcs, points):
    """An fk5 polygon region (exclusion) from pixel points."""
    ra, dec = wcs.pixel_to_world_values(np.array([p[0] for p in points]),
                                        np.array([p[1] for p in points]))
    return 'fk5;-polygon(%s)' % ','.join('%.6f,%.6f' % (a, d)
                                         for a, d in zip(ra, dec))


def background_lines(wcs, scale, center_ra, center_dec, cx, cy, bkg):
    """The fk5 region lines of a background from background()."""
    if bkg['type'] == 'circle':
        p = wcs.pixel_to_world(*bkg['offset'])
        return ['fk5;circle(%.6f,%.6f,%.1f")'
                % (p.ra.deg, p.dec.deg, OFFSET_RADIUS)]
    if bkg['type'] != 'annulus':
        return []
    lines = ['fk5;annulus(%.6f,%.6f,%.1f",%.1f")'
             % (center_ra, center_dec, BKG_INNER, BKG_OUTER)]
    for ex, ey, er in bkg['masks']:
        p = wcs.pixel_to_world(ex, ey)
        lines.append('fk5;-circle(%.6f,%.6f,%.1f")'
                     % (p.ra.deg, p.dec.deg, er))
    for sx, sy, vx, vy, half in bkg['strips']:
        hl, hw = half / scale, STRIP_HALF_WIDTH / scale
        lines.append(_polygon(wcs, [
            (sx + a * hl * vx - b * hw * vy, sy + a * hl * vy + b * hw * vx)
            for a, b in ((1, 1), (1, -1), (-1, -1), (-1, 1))]))
    reach = (BKG_OUTER + 2) / scale
    width = 360.0 / N_SECTORS
    for k in bkg['sectors']:
        arc = np.radians(np.linspace(k * width, (k + 1) * width, 10))
        lines.append(_polygon(wcs, [(cx, cy)] + [
            (cx + reach * math.cos(a), cy + reach * math.sin(a))
            for a in arc]))
    return lines


def file_masks(data, expmap, header, wcs, scale, src_ra, src_dec, detections,
               profile):
    """
    The sky level and every detection's mask radius for one observation and
    filter, measured on its deepest exposure. Returns (sky or None,
    [(radius, bright), ...] in the order of detections).
    """
    x0, y0 = [float(v) for v in wcs.world_to_pixel_values(src_ra, src_dec)]
    prelim = []
    for d in detections:
        ex, ey = [float(v) for v in wcs.world_to_pixel_values(d['ra'],
                                                              d['dec'])]
        prelim.append((ex, ey, exclusion_radius(d['mag'])))
    sky = sky_level(data, expmap, scale, x0, y0, prelim)
    bright = [d['cpf'] >= BRIGHT_CPF or d['saturated'] for d in detections]
    pos = np.array([[p[0], p[1]] for p in prelim]) if prelim else None
    radii = []
    for i, d in enumerate(detections):
        rate = d['rate'] if d['rate'] == d['rate'] else 0.0
        halo = False
        if not bright[i]:
            dist = np.hypot(*(pos - pos[i]).T) * scale
            halo = any(dist[j] <= HALO_REACH and (
                bright[j] or detections[j]['rate'] >= HALO_RATIO * rate)
                for j in range(len(detections)) if j != i)
        if bright[i]:
            radii.append((MASK_MAX, 'yes'))
        elif halo or sky is None or profile is None:
            radii.append((exclusion_radius(d['mag']), 'halo' if halo else 'no'))
        else:
            radii.append((mask_radius(rate, sky, profile), 'no'))
    return sky, radii


def _key(row):
    return (row['obsid'], row['filter'], row['extname'])


def _true_wcs(wcs, pointing, dec):
    """The image WCS corrected by a field-star pointing offset larger than
    fs.RECENTRE_MIN, else the image WCS itself."""
    if pointing and math.hypot(pointing['d_ra'], pointing['d_dec']) \
            > fs.RECENTRE_MIN:
        return fs.corrected_wcs(wcs, pointing, dec)
    return wcs


def process_file(args):
    """Measure every candidate exposure in one sky image."""
    (sky_path, rows, src_ra, src_dec, detections, sss_paths, profiles,
     overrides, outdir, pointings) = args
    pointings = pointings or {}
    results = []
    try:
        sky = _open_fits(sky_path)
        ex = _open_fits(expmap_path(sky_path))
        ex_by_name = {h.header.get('EXTNAME'): h for h in ex[1:]}
        sss = {lvl: fits.open(p, memmap=True) for lvl, p in sss_paths.items()}
    except Exception as exc:
        return [dict(r, status='failed', reason='cannot read: %s' % exc)
                for r in rows], {}
    dets = detections or []
    catpos = SkyCoord(src_ra, src_dec, unit='deg')
    det_pos = SkyCoord([d['ra'] for d in dets], [d['dec'] for d in dets],
                       unit='deg') if dets else None
    # Your source is the brightest detection within 10" of --ra/--dec: in a
    # trailed or mis-pointed exposure it can be several arcsec from the
    # region centre, and must not be masked as a neighbour.
    near = [(d['mag'], i) for i, d in enumerate(dets)
            if det_pos[i].separation(catpos).arcsec <= SEARCH_RADIUS + 2]
    target = min(near)[1] if near else None

    # the masks: sized once per observation and filter, on its deepest
    # exposure (aspect-corrected preferred, as uvotdetect's)
    masks_info = {}
    radii = [(exclusion_radius(d['mag']), 'no') for d in dets]
    sky_sb = None
    try:
        deep = max(rows, key=lambda r: (r['aspcorr'] == 'DIRECT',
                                        float(r['exposure'])))
        hdu = sky[int(deep['hdu'])]
        sky_sb, radii = file_masks(
            hdu.data, ex_by_name[deep['extname']].data, hdu.header,
            _true_wcs(WCS(hdu.header), pointings.get(_key(deep)), src_dec),
            abs(hdu.header['CDELT1']) * 3600.0, src_ra,
            src_dec, dets, profiles.get(deep['filter']))
        masks_info = {'sky': sky_sb, 'radii': radii}
    except Exception as exc:
        masks_info = {'error': '%s: %s' % (exc.__class__.__name__, exc)}

    for row in rows:
        out = {k: row.get(k) for k in ('obsid', 'filter', 'hdu', 'extname',
                                       'aspcorr', 'mode', 'first')}
        reasons = []
        if 'error' in masks_info:
            reasons.append('mask sizes not measured (%s): uvotlc radii'
                           % masks_info['error'])
        try:
            hdu = sky[int(row['hdu'])]
            h, data = hdu.header, hdu.data
            expmap = ex_by_name[row['extname']].data
            # wcs gives true sky coordinates: the image's own, or corrected
            # by the offset its field stars show (unaspected exposures can
            # be minutes of arc off); wcs_img is the image's, which the
            # region files must use
            wcs_img = WCS(h)
            pointing = pointings.get(_key(row))
            wcs = _true_wcs(wcs_img, pointing, src_dec)
            shifted = wcs is not wcs_img
            if shifted:
                reasons.append('pointing corrected by %+.1f", %+.1f" (RA, '
                               'Dec; from %s)' % (
                                   pointing['d_ra'], pointing['d_dec'],
                                   '%d field stars' % pointing['n']
                                   if pointing['source'] == 'stars' else
                                   "the snapshot's other filters"))
            scale = abs(h['CDELT1']) * 3600.0
            x0, y0 = [float(v) for v in wcs.world_to_pixel_values(src_ra, src_dec)]
            det_xy = [[float(v) for v in wcs.world_to_pixel_values(d['ra'],
                                                                   d['dec'])]
                      for d in dets]

            def geometry(centre):
                """Mask circles and strips (pixels) around a centre."""
                masks, strips = [], []
                cx_, cy_ = [float(v) for v in wcs.world_to_pixel_values(
                    centre.ra.deg, centre.dec.deg)]
                for i, d in enumerate(dets):
                    sep = det_pos[i].separation(centre).arcsec
                    radius, bright = radii[i]
                    if bright == 'yes':
                        vx, vy = readout_direction(h, *det_xy[i])
                        strips.append((det_xy[i][0], det_xy[i][1], vx, vy,
                                       sep + max(OFFSET_DISTANCES)
                                       + OFFSET_RADIUS))
                    if i == target or sep <= SAME_SOURCE:
                        continue
                    masks.append((det_xy[i][0], det_xy[i][1], radius))
                return masks, strips

            masks, strips = geometry(catpos)
            bkg0 = background(data, expmap, scale, x0, y0, masks, strips)
            meas = measure_source(data, wcs, scale, x0, y0, bkg0['mean'])
            center_ra, center_dec = src_ra, src_dec
            center = 'stars' if shifted else 'catalog'
            if meas is None:
                reasons.append('no source signal near the position')
            else:
                out['snr'] = '%.1f' % meas['snr']
                c = wcs.pixel_to_world(meas['x'], meas['y'])
                offset = c.separation(catpos).arcsec
                if meas['snr'] < MIN_SNR:
                    reasons.append('S/N %.1f: no shape; region at --ra/--dec'
                                   % meas['snr'])
                else:
                    out.update({
                        'offset': '%.2f' % offset,
                        'd_ra': '%+.2f' % ((c.ra.deg - src_ra) * 3600.0
                                           * math.cos(math.radians(src_dec))),
                        'd_dec': '%+.2f' % ((c.dec.deg - src_dec) * 3600.0),
                        'axis_ratio': '%.2f' % meas['axis_ratio'],
                        'sig_major': '%.2f' % meas['sig_major'],
                        'sig_minor': '%.2f' % meas['sig_minor'],
                        'major_pa': '%.0f' % meas['major_pa'],
                    })
                    if offset > MAX_SHIFT:
                        reasons.append('centroid %.1f arcsec away: region at '
                                       '--ra/--dec%s' % (
                                           offset, ' (pointing-corrected)'
                                           if shifted else ''))
                    else:
                        center_ra, center_dec, center = c.ra.deg, c.dec.deg, \
                            'centroid'
            ovr = matching_overrides(overrides, row)
            if 'src' in ovr:
                value, lineno = ovr['src']
                if value.lower() != 'catalog':
                    # an override is in the image's own coordinates
                    ira, idec = [float(v) for v in value.split()]
                    center_ra, center_dec = [float(v) for v in
                                             wcs.pixel_to_world_values(
                                                 *wcs_img.world_to_pixel_values(
                                                     ira, idec))]
                else:
                    center_ra, center_dec = src_ra, src_dec
                center = 'override:%d' % lineno
            cx, cy = [float(v) for v in wcs.world_to_pixel_values(center_ra,
                                                                  center_dec)]
            # the region, in the image's own coordinates (what uvotsource
            # reads)
            img_ra, img_dec = [float(v) for v in
                               wcs_img.pixel_to_world_values(cx, cy)]
            out.update({'center': center, 'src_ra': '%.6f' % img_ra,
                        'src_dec': '%.6f' % img_dec})
            centre = SkyCoord(center_ra, center_dec, unit='deg')
            if (cx, cy) != (x0, y0):
                masks, strips = geometry(centre)
                bkg = background(data, expmap, scale, cx, cy, masks, strips)
            else:
                bkg = bkg0
            if bkg['note']:
                reasons.append(bkg['note'])
            if detections is None:
                reasons.append('no source list (uvotdetect failed): no '
                               'exclusions')

            # How much of the source light near the region is inside it:
            # trailed and smeared images lose light the axis ratio can miss.
            # Reported only where it is precise enough to judge.
            if bkg['npix']:
                b = bkg['mean']
                conc, conc_err = concentration(
                    data, scale, cx, cy, b,
                    b / bkg['npix'] + (SKY_SYSTEMATIC * b) ** 2)
            else:
                conc, conc_err = None, None
            if conc is not None and conc_err <= CONC_ERR_MAX:
                out['conc'] = '%.2f' % conc
                out['conc_err'] = '%.3f' % conc_err
                if conc < CONC_MIN:
                    reasons.append('only %.0f%% of the counts within 15 arcsec '
                                   'are in the 5 arcsec circle: trailed or '
                                   'smeared' % (100 * conc))
            elif conc is not None:
                reasons.append('concentration not measured (error %.2f > '
                               '%.2f)' % (conc_err, CONC_ERR_MAX))
            else:
                reasons.append('concentration not measured (no net counts '
                               'within 15 arcsec)')

            detx, dety, rawx, rawy = detector_position(h, img_ra, img_dec)
            out.update({'detx': '%.1f' % detx, 'dety': '%.1f' % dety})
            for lvl in SSS_LEVELS:
                val = sss[lvl]['SSSENS' + row['filter']].data[
                    int(rawy + 0.5), int(rawx + 0.5)]
                out['sss_%s' % lvl.lower()] = '%.1f' % val
            frac = sss_low_fraction(h, scale, cx, cy,
                                    sss['LOW']['SSSENS' + row['filter']].data)
            out['sss_low_frac'] = '%.2f' % frac

            reg_dir = os.path.join(outdir, 'regions', row['obsid'])
            os.makedirs(reg_dir, exist_ok=True)
            src_reg = os.path.join(reg_dir, row['extname'] + '_src.reg')
            bkg_reg = os.path.join(reg_dir, row['extname'] + '_bkg.reg')
            with open(src_reg, 'w') as fh:
                fh.write('fk5;circle(%.6f,%.6f,%.1f")\n'
                         % (img_ra, img_dec, SRC_RADIUS))
            btype, area = bkg['type'], bkg['area']
            n_excl, n_strip, n_sect = (len(bkg['masks']), len(bkg['strips']),
                                       len(bkg['sectors']))
            if 'bkg' in ovr:
                value, lineno = ovr['bkg']
                lines = ['fk5;' + value]
                btype, area = 'override:%d' % lineno, float('nan')
                n_excl = n_strip = n_sect = 0
            else:
                lines = background_lines(wcs_img, scale, img_ra, img_dec,
                                         cx, cy, bkg)
            if lines:
                with open(bkg_reg, 'w') as fh:
                    fh.write('\n'.join(lines) + '\n')
            elif os.path.exists(bkg_reg):
                os.remove(bkg_reg)
            # Step 4's background in the 5" circle, counts/s: compare with
            # Step 5's (uvotsource on the same regions)
            bkg5 = '-'
            if bkg['type'] in ('annulus', 'circle') and 'bkg' not in ovr:
                bkg5 = '%.5f' % (bkg['mean'] * math.pi * SRC_RADIUS ** 2
                                 / scale ** 2 / bkg['ref'])
            out.update({'bkg_type': btype,
                        'bkg_area': '-' if area != area else '%.0f' % area,
                        'n_excl': n_excl, 'n_strip': n_strip,
                        'n_sect': n_sect, 'bkg5': bkg5,
                        'sky': '-' if sky_sb is None else '%.3e' % sky_sb,
                        'src_reg': os.path.relpath(src_reg, outdir),
                        'bkg_reg': os.path.relpath(bkg_reg, outdir)
                        if lines else '-'})
            out['status'] = 'ok' if lines else 'no_background'
            if shifted:
                # Step 3 judged the coverage with the image's own (wrong)
                # coordinates: judge it again at the corrected position
                (emap,), ex_x, ex_y = _cutout((expmap,), cx, cy,
                                              SRC_RADIUS + 1, scale)
                inside = np.hypot(ex_x - cx, ex_y - cy) * scale <= SRC_RADIUS
                if not inside.any() or emap[inside].min() \
                        < 0.5 * float(np.max(expmap)):
                    out['status'] = 'not_covered'
                    reasons.append('with the corrected pointing the source '
                                   'is not (fully) in the exposed field')
        except Exception as exc:
            out['status'] = 'failed'
            reasons.append('%s: %s' % (exc.__class__.__name__, exc))
        out['reason'] = '; '.join(reasons)
        results.append(out)
    key = (rows[0]['obsid'], rows[0]['filter'])
    return results, {key: masks_info}


def run_detections(groups, outdir, nproc):
    """
    uvotdetect on the deepest candidate exposure of each observation and
    filter (aspect-corrected ones preferred). Returns ({(obsid, filter):
    [detection, ...] or None}, {(obsid, filter): reason}); a detection is a
    dict of ra, dec, mag, rate (coincidence-corrected counts/s), cpf (raw
    counts per frame) and saturated.
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
        keys.append((key, float(best.get('frametime') or 0.011032), dest))
    results = run_many(calls, nproc=nproc)
    detections, problems = {}, {}
    for (key, frametime, dest), res in zip(keys, results):
        if not res.ok:
            detections[key] = None
            problems[key] = res.reason
            continue
        tab = fits.getdata(dest, 'SOURCES')
        detections[key] = [
            {'ra': float(r['RA']), 'dec': float(r['DEC']),
             'mag': float(r['MAG']), 'rate': float(r['COI_RATE']),
             'cpf': float(r['RAW_TOT_RATE']) * frametime,
             'saturated': bool(r['TOT_SATURATED'])} for r in tab]
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
    parser.add_argument('--nproc', type=int,
                        default=uenv.default_nproc(),
                        help='parallel workers (default: the CPUs this run may use, '
                             'at most 8)')
    parser.add_argument('--no-field-stars', action='store_true',
                        help='skip the Gaia field stars (no pointing '
                             'correction, no PSF columns)')
    parser.add_argument('--gaia-offline', action='store_true',
                        help='use the cached Gaia list only; stop if there '
                             'is none for this position')
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
        profiles = wing_profiles()
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

    # field stars first: each exposure's pointing (used for the regions)
    # and PSF
    pointings, psf, star_rows, fs_note = {}, {}, [], ''
    if not args.no_field_stars:
        try:
            gaia, gaia_path = fs.gaia_stars(args.ra, args.dec, args.outdir,
                                            offline=args.gaia_offline)
        except Exception as exc:
            print('ERROR: Gaia field stars: %s\n(no network? run with '
                  '--no-field-stars to go on without the pointing and PSF '
                  'checks)' % exc, file=sys.stderr)
            return 1
        stars = fs.select(gaia, args.ra, args.dec)
        print('[info] Gaia DR3 (%s): %d sources; %d for registration, %d '
              'isolated PSF stars' % (gaia_path, len(gaia), len(stars[0]),
                                      len(stars[1])))
        print('[info] field stars in %d exposures ...' % len(candidates))
        with ProcessPoolExecutor(max_workers=max(1, args.nproc)) as pool:
            parts = list(pool.map(fs.measure_file, [
                (path, rows, stars, args.ra, args.dec)
                for path, rows in sorted(by_file.items())]))
        raw, meta, snap = {}, {}, {}
        for r in candidates:
            k = _key(r)
            snap[k] = (r['obsid'], r['snapshot'])
            meta[k] = {'snapshot': snap[k],
                       'normal': r['aspcorr'] == 'DIRECT' and not
                       fs.JITTER[0] <= r['date_mid'][:10] < fs.JITTER[1]}
        for part in parts:
            for key, pointing, srows, err in part:
                raw[key] = pointing
                star_rows.extend(srows)
                if err:
                    fs_note = err
        pointings = fs.pointing_fallback(raw, snap)
        psf = fs.normalise(star_rows, meta)
    jobs = [(path, rows, args.ra, args.dec,
             detections.get((rows[0]['obsid'], rows[0]['filter'])),
             sss_paths, profiles, overrides, args.outdir,
             {_key(r): pointings[_key(r)] for r in rows
              if _key(r) in pointings})
            for path, rows in sorted(by_file.items())]
    print('[info] measuring %d exposures in %d sky images ...'
          % (len(candidates), len(jobs)))
    with ProcessPoolExecutor(max_workers=max(1, args.nproc)) as pool:
        parts = list(pool.map(process_file, jobs))
    rows = [r for part, _ in parts for r in part]
    mask_sizes = {}
    for _, info in parts:
        mask_sizes.update(info)
    order = {(r['obsid'], r['filter'], r['extname']): i
             for i, r in enumerate(candidates)}
    # the ledger: every candidate must come back exactly once
    if sorted(order) != sorted((r['obsid'], r['filter'], r['extname'])
                               for r in rows):
        print('ERROR: internal: %d candidates in, %d rows out'
              % (len(candidates), len(rows)), file=sys.stderr)
        return 1
    rows.sort(key=lambda r: order[(r['obsid'], r['filter'], r['extname'])])
    for r in rows:
        k = _key(r)
        p = pointings.get(k)
        if p:
            r.update({'pnt_dra': '%+.2f' % p['d_ra'],
                      'pnt_ddec': '%+.2f' % p['d_dec'], 'pnt_n': p['n'],
                      'pnt_err': '%.2f' % p['err'], 'pnt_src': p['source']})
        q = psf.get(k)
        if q:
            r.update({'n_stars': q['n_stars'],
                      'psf_loss': '%.3f' % q['psf_loss'],
                      'psf_loss_err': '%.3f' % q['psf_loss_err'],
                      'psf_src': q['psf_source']})
            if 'star_conc' in q:
                r.update({'star_conc': '%.3f' % q['star_conc'],
                          'star_conc_err': '%.3f' % q['star_conc_err'],
                          'star_ar': '%.2f' % q['star_ar'],
                          'star_sigma': '%.2f' % q['star_sigma']})
        elif not args.no_field_stars:
            r['n_stars'] = 0

    # sources near the target, for the record
    src = SkyCoord(args.ra, args.dec, unit='deg')
    det_rows = []
    for (obsid, filt), dets in sorted(detections.items()):
        info = mask_sizes.get((obsid, filt), {})
        radii = info.get('radii') or [(exclusion_radius(d['mag']), 'no')
                                      for d in dets or []]
        for d, (radius, bright) in zip(dets or [], radii):
            sep = SkyCoord(d['ra'], d['dec'], unit='deg').separation(src).arcsec
            if sep <= DETECTIONS_KEPT:
                det_rows.append({'obsid': obsid, 'filter': filt,
                                 'ra': '%.6f' % d['ra'], 'dec': '%.6f' % d['dec'],
                                 'sep': '%.1f' % sep, 'mag': '%.2f' % d['mag'],
                                 'rate': '%.4g' % d['rate'],
                                 'cpf': '%.3f' % d['cpf'],
                                 'bright': bright,
                                 'mask_radius': '%.1f' % radius})

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
        'sig_major/minor (arcsec), major_pa (deg E of N): source shape, '
        'only at S/N >= %.0f' % MIN_SNR,
        'conc, conc_err: net counts in the 5" region / in 15" (point source '
        '0.86-0.92) and its error, only where the error is <= %.2f'
        % CONC_ERR_MAX,
        'center: where the 5" source region is (centroid, catalog or an '
        'override line); sss_*: CALDB small-scale-sensitivity factor at the '
        'source (-99.9 = low-sensitivity patch); sss_low_frac: fraction of '
        'the 5" circle on a LOW patch',
        'bkg_type: annulus (27.5-35" minus n_excl source circles, n_strip '
        'read-out strips and n_sect bright sectors), circle (15" offset '
        'circle), override, or none; bkg_area in arcsec^2; bkg5: that '
        'background in the 5" circle, counts/s; sky: counts/s/arcsec^2 '
        'the mask radii were set against',
        'pnt_dra, pnt_ddec, pnt_err: pointing offset (image - Gaia, arcsec) '
        'and its uncertainty, from pnt_n field stars or borrowed from the '
        'snapshot\'s other filters (pnt_src); applied to the regions '
        '(center = stars) when over %.0f"; status not_covered if the source '
        'is then outside the exposed field' % fs.RECENTRE_MIN,
        'n_stars, star_conc(_err), star_ar, star_sigma: isolated Gaia stars '
        '(S/N >= %.0f, <= %.1f counts/frame), their concentration relative '
        'to their median in normal exposures, axis ratio, width (arcsec); '
        'psf_loss(_err) = 1 - star_conc: light the 5" circle loses to a '
        'broad PSF (psf_src: stars, or snapshot)' % (fs.PSF_MIN_SNR,
                                                     fs.PSF_MAX_CPF),
    ]
    pos_path = os.path.join(args.outdir, 'uvot_positions.txt')
    write_table(pos_path, rows, COLUMNS, comments)
    det_path = os.path.join(args.outdir, 'uvot_detections.txt')
    write_table(det_path, det_rows, ['obsid', 'filter', 'ra', 'dec', 'sep',
                                     'mag', 'rate', 'cpf', 'bright',
                                     'mask_radius'],
                comments[:4] + [
                    'uvotdetect sources within %.0f" of the source in the '
                    'deepest exposure of each observation and filter'
                    % DETECTIONS_KEPT,
                    'rate: coincidence-corrected counts/s; cpf: raw counts per '
                    'frame; bright: yes if cpf >= %.1f or saturated (mask %.0f" '
                    'and a read-out strip), halo if within %.0f" of a source '
                    '%.0fx brighter (uvotlc radius); mask_radius: arcsec, '
                    'where the PSF wing falls below %.1f%% of the sky' % (
                        BRIGHT_CPF, MASK_MAX, HALO_REACH, HALO_RATIO,
                        100 * WING_FRACTION)])

    star_path = None
    if star_rows:
        star_path = os.path.join(args.outdir, 'uvot_field_stars.txt')
        write_table(star_path, star_rows, fs.STAR_COLUMNS, comments[:4] + [
            'Gaia field stars per exposure (PSF stars): offset from Gaia '
            '(image - Gaia, arcsec, proper motion applied), S/N, raw counts '
            'per frame in 5", axis ratio, width (arcsec), concentration '
            '(5"/15") with its error, and relative to the star\'s median in '
            'normal exposures'])

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
    print('[summary] background: %d with read-out strips, %d with bright '
          'sectors dropped (%d sectors)' % (
              sum(int(r.get('n_strip') or 0) > 0 for r in rows),
              sum(int(r.get('n_sect') or 0) > 0 for r in rows),
              sum(int(r.get('n_sect') or 0) for r in rows)))
    nconc = sum(r.get('conc') not in (None, '-') for r in rows)
    print('[summary] concentration measured (error <= %.2f): %d of %d'
          % (CONC_ERR_MAX, nconc, len(rows)))
    edge = [r for r in rows if r.get('sss_low_frac') not in (None, '-')
            and float(r['sss_low_frac']) > SSS_EDGE_FRACTION]
    print('[summary] more than %.0f%% of the 5" circle on a LOW patch: %d'
          % (100 * SSS_EDGE_FRACTION, len(edge)))
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
    if not args.no_field_stars:
        moved = Counter(r.get('pnt_src') for r in rows
                        if r.get('center') == 'stars')
        print('[summary] pointing measured from field stars: %d of %d (%d '
              'more from the snapshot); corrected (over %.0f"): %d, of which '
              '%d regions placed by it (%s) and %d not covered' % (
                  sum(r.get('pnt_src') == 'stars' for r in rows), len(rows),
                  sum(r.get('pnt_src') == 'snapshot' for r in rows),
                  fs.RECENTRE_MIN, sum('pointing corrected' in r['reason']
                                       for r in rows),
                  sum(moved.values()), ', '.join(
                      'from %s %d' % kv for kv in sorted(moved.items())) or
                  '-', sum(r['status'] == 'not_covered' for r in rows)))
        big = [r for r in rows if r.get('pnt_dra') not in (None, '-') and
               math.hypot(float(r['pnt_dra']), float(r['pnt_ddec'])) > 3]
        for r in big[:8]:
            print('          %s %s %s: %s", %s" (%s)' % (
                r['obsid'], r['filter'], r['extname'], r['pnt_dra'],
                r['pnt_ddec'], r['aspcorr']))
        if len(big) > 8:
            print('          ... %d more over 3"' % (len(big) - 8))
        lossy = [r for r in rows if r.get('psf_loss') not in (None, '-')
                 and float(r['psf_loss']) > 0.03]
        print('[summary] PSF from field stars: %d of %d (%d from the '
              'snapshot\'s other filters); psf_loss > 3 %%: %d' % (
                  sum(r.get('psf_loss') not in (None, '-') for r in rows),
                  len(rows), sum(r.get('psf_src') == 'snapshot' for r in rows),
                  len(lossy)))
        if fs_note:
            print('[warn] field stars: %s' % fs_note)
    if det_problems:
        print('[warn] uvotdetect failed for %d observation/filter pair(s): '
              'their backgrounds have no source exclusions' % len(det_problems))
        for key, why in list(det_problems.items())[:5]:
            print('          %s %s: %s' % (key[0], key[1], why))
    print('[summary] Written: %s, %s, %s%s/regions/' % (
        pos_path, det_path, star_path + ', ' if star_path else '',
        args.outdir))
    if counts['failed']:
        print('ERROR: %d exposure(s) failed; see the reason column.'
              % counts['failed'], file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
