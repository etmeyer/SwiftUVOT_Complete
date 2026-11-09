"""
swift_uvot_fieldstars.py -- Gaia field stars for Step 4: the pointing and
the PSF of every exposure.

A faint target cannot tell us where it is or how sharp the image is: its
own centroid and shape are noise. The stars around it can.

  * Pointing. Without aspect correction, or during the 2023-24 jitter, an
    exposure's sky coordinates can be off by several arcsec -- more than
    the 5" aperture (PKS 0537-286's unaspected exposures in 00035240001
    are 8-11" off). The Gaia stars are lined up with the image: the shift
    that puts the most starlight at their positions, refined by the
    centroids of the isolated stars. Step 4 moves a faint target's region
    by it.
  * PSF. Each isolated, unsaturated star's concentration (net counts in
    5" over 15") falls when the image is smeared. Divided by the same
    star's median in normal exposures, it measures the light the 5"
    aperture loses in that exposure (about 1:1 for stars): psf_loss =
    1 - relative concentration.

The Gaia DR3 list comes from the ESA archive (TAP, asynchronous) once per
analysis and is cached in the output folder (uvot_gaia.txt, with the
service, query, date and position in its header).

This module is imported by swift_uvot_positions.py; it is not a step of
its own.
"""

import math
import os
import time
from datetime import datetime, timezone

import numpy as np

from swift_uvot_tables import read_table, write_table

GAIA_TAP = 'https://gea.esac.esa.int/tap-server/tap'
GAIA_FILE = 'uvot_gaia.txt'
GAIA_RADIUS = 9.0           # arcmin around the target queried (registration
                            # pairs image peaks up to 4' away with stars up to 5' off)
GAIA_G_MAX = 21.0           # faintest Gaia source kept (masks, isolation)
GAIA_EPOCH = 2016.0         # Gaia DR3 reference epoch (Julian year)

# Swift's gyro noise: UVOT images degraded from about 2023-08-07 (GCN 34633)
# until UVOT was re-enabled in two-gyro mode on 2024-04-04 (GCN 36033)
JITTER = ('2023-08-07', '2024-04-04')

# Registration: image peaks are paired with every Gaia source within
# REG_MAX of them; the offset (image - Gaia) shared by the most pairs wins.
# Unaspected exposures can be minutes of arc off (PKS 0537-286's
# 00035240001 snapshot 2: 205"), so the search is wide; random pairs spread
# over the whole range and the true offset collects one vote per star.
REG_G = 19.5
REG_MAX = 300.0             # arcsec
REG_BIN = 2.0               # arcsec, vote histogram bin
REG_PEAKS = 60              # brightest image peaks used
REG_PEAK_RADIUS = 240.0     # arcsec around the target for image peaks
REG_PEAK_SIGMA = 6.0        # peak significance in a 4" box
REG_MIN_VOTES = 6           # pairs that must agree
REG_SEARCH = 3.0            # arcsec around the voted position for a centroid
REG_MIN_STARS = 3           # centroided stars that must agree for a shift
REG_SPREAD = 0.7            # arcsec: their offsets must agree this well (MAD)
REG_MIN_SNR = 5.0           # a star's centroid counts from this S/N
RECENTRE_MIN = 1.0          # arcsec: smaller shifts are not applied

# PSF stars: isolated, clear and not near coincidence loss
PSF_RADIUS = 300.0          # arcsec from the target
PSF_G_MAX = 19.5
PSF_ISOLATION = 15.0        # arcsec: no Gaia source brighter than G + 2.5
PSF_MIN_SNR = 20.0
PSF_MAX_CPF = 0.3           # raw counts per frame in 5"
PSF_CONC_ERR = 0.05         # a star's concentration used if this precise
PSF_MAX_RUWE = 1.4          # Gaia astrometric quality: binaries etc. out
PSF_SYSTEMATIC = 0.01       # star-to-star scatter of the relative value
PSF_REF_MIN = 3             # normal exposures a star needs for its median
MASK_RADIUS = 6.0           # arcsec around other Gaia sources in an annulus
TARGET_CLEAR = 15.0         # arcsec: stars this close to the target not used

STAR_COLUMNS = ['obsid', 'filter', 'extname', 'star', 'g', 'd_ra', 'd_dec',
                'snr', 'cpf', 'axis_ratio', 'sigma', 'conc', 'conc_err',
                'rel_conc']


# --- Gaia -------------------------------------------------------------------

def _tap_query(query, timeout=300):
    """CSV text of an ADQL query on the ESA Gaia archive (asynchronous TAP:
    the synchronous queue is often full)."""
    import requests
    r = requests.post(GAIA_TAP + '/async', timeout=60, allow_redirects=False,
                      data={'REQUEST': 'doQuery', 'LANG': 'ADQL',
                            'FORMAT': 'csv', 'PHASE': 'RUN', 'QUERY': query})
    job = r.headers.get('Location')
    if r.status_code not in (200, 303) or not job:
        raise RuntimeError('Gaia archive refused the query (HTTP %d)'
                           % r.status_code)
    start = time.time()
    while True:
        phase = requests.get(job + '/phase', timeout=30).text.strip()
        if phase == 'COMPLETED':
            break
        if phase in ('ERROR', 'ABORTED'):
            raise RuntimeError('Gaia query %s: %s' % (phase, requests.get(
                job + '/error', timeout=30).text[:300]))
        if time.time() - start > timeout:
            raise RuntimeError('Gaia query still %s after %d s'
                               % (phase, timeout))
        time.sleep(2)
    return requests.get(job + '/results/result', timeout=120).text


def gaia_stars(ra, dec, outdir, offline=False):
    """
    Gaia DR3 sources within GAIA_RADIUS of (ra, dec), G < GAIA_G_MAX, as a
    list of dicts (source_id, ra, dec, pmra, pmdec, g, ruwe). Read from
    <outdir>/uvot_gaia.txt if it was made for this position; otherwise
    queried and cached. offline: never query (RuntimeError if no cache).
    """
    path = os.path.join(outdir, GAIA_FILE)
    if os.path.exists(path):
        with open(path) as fh:
            header = [line for line in fh if line.startswith('#')]
        pos = [line for line in header if line.startswith('# position:')]
        if pos:
            p = pos[0].split()
            same = (abs(float(p[3]) - ra) * math.cos(math.radians(dec))
                    < 1 / 3600 and abs(float(p[5]) - dec) < 1 / 3600
                    and float(p[7]) >= GAIA_RADIUS)
            if same:
                return [dict(r, ra=float(r['ra']), dec=float(r['dec']),
                             pmra=_float(r['pmra'], 0.0),
                             pmdec=_float(r['pmdec'], 0.0),
                             g=_float(r['g'], 99.0),
                             ruwe=_float(r['ruwe'], 1.0))
                        for r in read_table(path)], path
    if offline:
        raise RuntimeError('no Gaia list for this position in %s, and '
                           'querying is off' % path)
    query = ("SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, ruwe "
             "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, "
             "dec), CIRCLE('ICRS', %.7f, %.7f, %.5f)) AND phot_g_mean_mag < "
             "%.1f" % (ra, dec, GAIA_RADIUS / 60.0, GAIA_G_MAX))
    text = _tap_query(query)
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith('source_id'):
        raise RuntimeError('unexpected answer from the Gaia archive: %r'
                           % text[:200])
    rows = []
    for ln in lines[1:]:
        sid, sra, sdec, pmra, pmdec, g, ruwe = ln.split(',')
        rows.append({'source_id': sid, 'ra': float(sra), 'dec': float(sdec),
                     'pmra': _float(pmra, 0.0), 'pmdec': _float(pmdec, 0.0),
                     'g': _float(g, 99.0), 'ruwe': _float(ruwe, 1.0)})
    out = [{'source_id': r['source_id'], 'ra': '%.8f' % r['ra'],
            'dec': '%.8f' % r['dec'], 'pmra': '%.3f' % r['pmra'],
            'pmdec': '%.3f' % r['pmdec'], 'g': '%.3f' % r['g'],
            'ruwe': '%.3f' % r['ruwe']} for r in rows]
    write_table(path, out, ['source_id', 'ra', 'dec', 'pmra', 'pmdec', 'g',
                            'ruwe'], [
        'Gaia DR3 sources near the target, for the field-star checks',
        'service: %s (asynchronous TAP), %s' % (
            GAIA_TAP, datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')),
        'query: %s' % query,
        'units: deg, mas/yr (pmra includes cos dec), G mag; epoch %.1f'
        % GAIA_EPOCH])
    with open(path) as fh:
        text = fh.read()
    with open(path, 'w') as fh:     # the position line, for reuse checks
        fh.write('# position: RA %.7f Dec %.7f radius %.2f arcmin\n'
                 % (ra, dec, GAIA_RADIUS) + text)
    return rows, path


def _float(value, default):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return default
    return default if v != v else v


def select(gaia, ra, dec):
    """
    (registration sources, PSF stars, mask sources) from the Gaia list:
    registration: G < REG_G anywhere; PSF stars: G < PSF_G_MAX within
    PSF_RADIUS of the target (not within TARGET_CLEAR), RUWE <= 1.4, no
    Gaia source brighter than G + 2.5 within PSF_ISOLATION; masks: all.
    """
    if not gaia:
        return [], [], []
    g = np.array([s['g'] for s in gaia])
    sra = np.array([s['ra'] for s in gaia])
    sdec = np.array([s['dec'] for s in gaia])
    cosd = math.cos(math.radians(dec))

    def sep(ra0, dec0):
        return 3600 * np.hypot((sra - ra0) * cosd, sdec - dec0)

    d_target = sep(ra, dec)
    reg = [s for s, gi in zip(gaia, g) if gi < REG_G]
    psf = []
    for i, s in enumerate(gaia):
        if not (g[i] < PSF_G_MAX and TARGET_CLEAR < d_target[i] < PSF_RADIUS
                and s['ruwe'] <= PSF_MAX_RUWE):
            continue
        d = sep(s['ra'], s['dec'])
        if np.any((d > 0.3) & (d < PSF_ISOLATION) & (g < g[i] + 2.5)):
            continue
        psf.append(s)
    return reg, psf, gaia


def at_epoch(stars, jyear):
    """(ra, dec) arrays at a Julian year, proper motion applied."""
    dt = jyear - GAIA_EPOCH
    dec = np.array([s['dec'] + s['pmdec'] * dt / 3.6e6 for s in stars])
    ra = np.array([s['ra'] + s['pmra'] * dt / 3.6e6
                   / math.cos(math.radians(s['dec'])) for s in stars])
    return ra, dec


# --- per exposure -----------------------------------------------------------

def _box_sums(data, half):
    """Sum of data in a (2 half + 1)^2 box around every pixel."""
    pad = np.pad(data, half + 1)
    c = pad.cumsum(0).cumsum(1)
    k = 2 * half + 1
    return (c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k])[
        :data.shape[0], :data.shape[1]]


def register(data, wcs, scale, target_xy, reg_ra, reg_dec):
    """
    Offset (image - Gaia) in pixels, (dx, dy), and the number of votes for
    it, from the brightest peaks within REG_PEAK_RADIUS of the target:
    each peak is paired with every Gaia source within REG_MAX, and the
    offset shared by the most pairs wins (2" bins, summed over 3x3).
    None if it does not win clearly (REG_MIN_VOTES, and twice the best
    offset more than 6" away).
    """
    from scipy.ndimage import maximum_filter
    if len(reg_ra) < REG_MIN_STARS:
        return None
    sums = _box_sums(data, max(1, int(round(2.0 / scale))))
    tx, ty = target_xy
    reach = int(REG_PEAK_RADIUS / scale)
    y0, y1 = max(0, int(ty) - reach), min(data.shape[0], int(ty) + reach)
    x0, x1 = max(0, int(tx) - reach), min(data.shape[1], int(tx) + reach)
    sub = sums[y0:y1, x0:x1]
    if sub.size == 0:
        return None
    pos = sub[sub > 0]
    med = float(np.median(pos)) if pos.size else 0.0
    sig = (sub - med) / math.sqrt(abs(med) + 1.0)
    peak = (sub == maximum_filter(sub, size=2 * int(4.0 / scale) + 1)) \
        & (sig > REG_PEAK_SIGMA)
    py, px = np.nonzero(peak)
    if len(px) < REG_MIN_STARS:
        return None
    order = np.argsort(sig[py, px])[::-1][:REG_PEAKS]
    px, py = px[order] + x0, py[order] + y0
    gx, gy = wcs.world_to_pixel_values(reg_ra, reg_dec)
    dx = (px[:, None] - gx[None, :]).ravel() * scale
    dy = (py[:, None] - gy[None, :]).ravel() * scale
    keep = (np.abs(dx) <= REG_MAX) & (np.abs(dy) <= REG_MAX)
    nb = int(2 * REG_MAX / REG_BIN)
    hist, xe, ye = np.histogram2d(dx[keep], dy[keep], bins=nb,
                                  range=[[-REG_MAX, REG_MAX]] * 2)
    smooth = sum(np.roll(np.roll(hist, i, 0), j, 1)
                 for i in (-1, 0, 1) for j in (-1, 0, 1))
    i, j = np.unravel_index(np.argmax(smooth), smooth.shape)
    votes = smooth[i, j]
    bx, by = (xe[i] + xe[i + 1]) / 2, (ye[j] + ye[j + 1]) / 2
    cx = (xe[:-1] + xe[1:]) / 2
    far = np.hypot(cx[:, None] - bx, cx[None, :] - by) > 6.0
    if votes < REG_MIN_VOTES or votes < 2 * smooth[far].max():
        return None
    near = keep & (np.hypot(dx - bx, dy - by) <= 1.5 * REG_BIN)
    return (float(np.median(dx[near])) / scale,
            float(np.median(dy[near])) / scale, int(votes))


def _annulus_background(data, xx, yy, x, y, scale, mask_xy):
    """Plain mean of the 27.5-35" annulus with other Gaia sources masked
    (as uvotsource's estimator), its pixel count; None if under half is
    left."""
    d = np.hypot(xx - x, yy - y) * scale
    ring = (d >= 27.5) & (d <= 35.0)
    keep = ring.copy()
    for mx, my in mask_xy:
        if abs(mx - x) * scale < 42 and abs(my - y) * scale < 42:
            keep &= np.hypot(xx - mx, yy - my) * scale > MASK_RADIUS
    if keep.sum() < 0.5 * ring.sum():
        return None, 0
    return float(data[yy[keep], xx[keep]].mean()), int(keep.sum())


def measure_star(data, scale, x0, y0, mask_xy, frametime, exposure):
    """
    A field star near pixel (x0, y0): centroid (within REG_SEARCH), S/N,
    raw counts per frame, axis ratio, width and concentration with its
    error. None if it is not on the image or not found.
    """
    from swift_uvot_positions import concentration, SKY_SYSTEMATIC
    reach = int(math.ceil(37.0 / scale))
    xc, yc = int(round(x0)), int(round(y0))
    ny, nx = data.shape
    if xc - reach < 0 or yc - reach < 0 or xc + reach >= nx or yc + reach >= ny:
        return None
    yy, xx = np.mgrid[yc - reach:yc + reach + 1, xc - reach:xc + reach + 1]
    bkg, npix = _annulus_background(data, xx, yy, x0, y0, scale, mask_xy)
    if bkg is None:
        return None
    cut = data[yc - reach:yc + reach + 1, xc - reach:xc + reach + 1]
    near = np.hypot(xx - x0, yy - y0) * scale <= REG_SEARCH
    iy, ix = np.unravel_index(np.argmax(np.where(near, cut, -np.inf)),
                              cut.shape)
    cx, cy = float(xx[iy, ix]), float(yy[iy, ix])
    for _ in range(5):
        m = np.hypot(xx - cx, yy - cy) * scale <= 4.0
        w = np.clip(cut[m] - bkg, 0, None)
        if w.sum() <= 0:
            return None
        cx, cy = (xx[m] * w).sum() / w.sum(), (yy[m] * w).sum() / w.sum()
    d = np.hypot(xx - cx, yy - cy) * scale
    m5 = d <= 5.0
    tot5 = cut[m5].sum()
    net5 = tot5 - bkg * m5.sum()
    snr = net5 / math.sqrt(tot5) if tot5 > 0 else 0.0
    ms = d <= 6.0
    w = np.clip(cut[ms] - bkg, 0, None)
    out = {'x': cx, 'y': cy, 'snr': snr,
           'cpf': tot5 / exposure * frametime if exposure > 0 else 0.0}
    if w.sum() > 0:
        dx, dy = xx[ms] - cx, yy[ms] - cy
        cov = np.array([[(w * dx * dx).sum(), (w * dx * dy).sum()],
                        [(w * dx * dy).sum(), (w * dy * dy).sum()]]) / w.sum()
        vals = np.linalg.eigvalsh(cov)
        vmin, vmax = max(vals[0], 1e-9), max(vals[1], 1e-9)
        out['axis_ratio'] = math.sqrt(vmax / vmin)
        out['sigma'] = math.sqrt((vmax + vmin) / 2) * scale
    conc, err = concentration(data, scale, cx, cy, bkg,
                              bkg / npix + (SKY_SYSTEMATIC * bkg) ** 2)
    out['conc'], out['conc_err'] = conc, err
    return out


def measure_exposure(data, wcs, scale, jyear, frametime, exposure, stars,
                     target_xy):
    """
    The field stars of one exposure. stars: (registration sources, PSF
    stars, mask sources) from select(). Returns (pointing, per-star rows):
    pointing = dict(d_ra, d_dec (arcsec, image minus Gaia: add to a
    catalogue position to get where it is in this image), n, spread) or
    None; per-star rows with the measurements of the PSF stars.
    """
    reg, psf, masks = stars
    data = np.nan_to_num(np.asarray(data, dtype=float))
    rra, rdec = at_epoch(reg, jyear)
    shift = register(data, wcs, scale, target_xy, rra, rdec)
    if shift is None:
        return None, []
    votes = shift[2]
    mra, mdec = at_epoch(masks, jyear)
    mx, my = wcs.world_to_pixel_values(mra, mdec)
    mask_xy = list(zip(mx, my))     # catalogue positions; shifted below
    pra, pdec = at_epoch(psf, jyear)
    offsets, rows = [], []
    for s, ra, dec in zip(psf, pra, pdec):
        x0, y0 = [float(v) for v in wcs.world_to_pixel_values(ra, dec)]
        own = [(a + shift[0], b + shift[1]) for a, b in mask_xy
               if math.hypot(a - x0, b - y0) * scale > 0.5]
        m = measure_star(data, scale, x0 + shift[0], y0 + shift[1], own,
                         frametime, exposure)
        if m is None or m['snr'] < REG_MIN_SNR:
            continue
        c = wcs.pixel_to_world(m['x'], m['y'])
        cosd = math.cos(math.radians(dec))
        d_ra = (c.ra.deg - ra) * 3600.0 * cosd
        d_dec = (c.dec.deg - dec) * 3600.0
        offsets.append((d_ra, d_dec))
        if m['snr'] >= PSF_MIN_SNR and m['cpf'] <= PSF_MAX_CPF \
                and m['conc'] is not None and m['conc_err'] <= PSF_CONC_ERR:
            rows.append({'star': s['source_id'], 'g': '%.2f' % s['g'],
                         'd_ra': '%+.2f' % d_ra, 'd_dec': '%+.2f' % d_dec,
                         'snr': '%.1f' % m['snr'], 'cpf': '%.3f' % m['cpf'],
                         'axis_ratio': '%.3f' % m.get('axis_ratio', np.nan),
                         'sigma': '%.3f' % m.get('sigma', np.nan),
                         'conc': '%.4f' % m['conc'],
                         'conc_err': '%.4f' % m['conc_err']})
    pointing = None
    if len(offsets) >= REG_MIN_STARS:
        o = np.array(offsets)
        med = np.median(o, axis=0)
        dist = np.hypot(*(o - med).T)
        spread = 1.4826 * np.median(dist)
        good = dist <= max(3 * spread, 0.5)
        if good.sum() >= REG_MIN_STARS:
            med = np.median(o[good], axis=0)
            pointing = {'d_ra': float(med[0]), 'd_dec': float(med[1]),
                        'n': int(good.sum()), 'spread': float(spread),
                        'votes': votes}
            if spread > REG_SPREAD:
                pointing['unreliable'] = True
    return pointing, rows


# --- after all exposures ----------------------------------------------------

def normalise(star_rows, rows_meta):
    """
    Relative concentrations and the PSF loss per exposure.

    star_rows: per-star rows (with obsid, filter, extname); rows_meta:
    {(obsid, filter, extname): dict(normal=bool, snapshot=key)}. A star's
    reference in a filter is its median concentration in normal exposures
    (aspect-corrected, outside the jitter period), from PSF_REF_MIN of
    them. Returns {(obsid, filter, extname): dict(n_stars, star_conc,
    star_conc_err, star_ar, star_sigma, psf_loss, psf_loss_err,
    psf_source)}, and fills rel_conc in star_rows.
    """
    ref = {}
    by_star = {}
    for r in star_rows:
        key = (r['obsid'], r['filter'], r['extname'])
        if rows_meta.get(key, {}).get('normal'):
            by_star.setdefault((r['star'], r['filter']), []).append(
                float(r['conc']))
    for k, v in by_star.items():
        if len(v) >= PSF_REF_MIN:
            ref[k] = float(np.median(v))
    per_exp = {}
    for r in star_rows:
        k = (r['star'], r['filter'])
        if k not in ref:
            r['rel_conc'] = '-'
            continue
        rel = float(r['conc']) / ref[k]
        err = math.hypot(float(r['conc_err']) / ref[k], PSF_SYSTEMATIC)
        r['rel_conc'] = '%.4f' % rel
        per_exp.setdefault((r['obsid'], r['filter'], r['extname']), []).append(
            (rel, err, float(r['axis_ratio']), float(r['sigma'])))
    out = {}
    for key, v in per_exp.items():
        rel = np.array([a[0] for a in v])
        err = np.array([a[1] for a in v])
        # one star that is itself odd (a neighbour, a cosmic ray) must not
        # decide: drop values more than 4 sigma from the median
        med = np.median(rel)
        keep = np.abs(rel - med) <= 4 * np.maximum(err, 1e-3) if len(v) > 2 \
            else np.ones(len(v), bool)
        w = 1 / err[keep] ** 2
        mean = float(np.sum(w * rel[keep]) / np.sum(w))
        sig = float(1 / math.sqrt(np.sum(w)))
        out[key] = {'n_stars': int(keep.sum()), 'star_conc': mean,
                    'star_conc_err': sig,
                    'star_ar': float(np.nanmedian([a[2] for a in v])),
                    'star_sigma': float(np.nanmedian([a[3] for a in v])),
                    'psf_loss': 1 - mean, 'psf_loss_err': sig,
                    'psf_source': 'stars'}
    # no star of its own (UVM2 often): the other filters of the same
    # snapshot share its pointing and jitter
    by_snap = {}
    for key, m in rows_meta.items():
        if key in out:
            by_snap.setdefault(m['snapshot'], []).append(out[key])
    for key, m in rows_meta.items():
        if key in out or m['snapshot'] not in by_snap:
            continue
        sib = by_snap[m['snapshot']]
        w = np.array([1 / s['psf_loss_err'] ** 2 for s in sib])
        losses = np.array([s['psf_loss'] for s in sib])
        loss = float(np.sum(w * losses) / w.sum())
        # the filters of one snapshot are taken one after another and the
        # jitter changes between them: their spread is the error
        err = max(float(1 / math.sqrt(w.sum())),
                  float(np.std(losses, ddof=1)) if len(sib) > 1 else 0.05)
        out[key] = {'n_stars': 0, 'psf_loss': loss, 'psf_loss_err': err,
                    'psf_source': 'snapshot'}
    return out


def measure_file(args):
    """
    Worker: the field stars of every candidate exposure in one sky image.
    args = (sky_path, rows, stars, src_ra, src_dec). Returns a list of
    (key, pointing, star rows, error) with key = (obsid, filter, extname).
    """
    from astropy.time import Time
    from astropy.wcs import WCS
    from swift_uvot_positions import _open_fits
    sky_path, rows, stars, src_ra, src_dec = args
    out = []
    try:
        sky = _open_fits(sky_path)
    except Exception as exc:
        return [((r['obsid'], r['filter'], r['extname']), None, [],
                 'cannot read: %s' % exc) for r in rows]
    for r in rows:
        key = (r['obsid'], r['filter'], r['extname'])
        try:
            hdu = sky[int(r['hdu'])]
            wcs = WCS(hdu.header)
            scale = abs(hdu.header['CDELT1']) * 3600.0
            txy = [float(v) for v in wcs.world_to_pixel_values(src_ra, src_dec)]
            jyear = Time(float(r['mjd_mid']), format='mjd').jyear
            pointing, srows = measure_exposure(
                hdu.data, wcs, scale, jyear, float(r['frametime']),
                float(r['exposure']), stars, txy)
            for s in srows:
                s.update(obsid=key[0], filter=key[1], extname=key[2])
            out.append((key, pointing, srows, ''))
        except Exception as exc:
            out.append((key, None, [], '%s: %s' % (exc.__class__.__name__,
                                                   exc)))
    return out


# arcsec: siblings must agree this well to lend their offset (unaspected
# filters of one snapshot differ by up to ~6"; what matters is whether the
# source is off by minutes of arc)
SNAPSHOT_SPREAD = 10.0


def pointing_fallback(pointings, snapshot_of):
    """
    Exposures without a pointing of their own (too few stars: short UVM2)
    borrow the median of the other filters of the same snapshot, which
    share the spacecraft attitude, if those agree within SNAPSHOT_SPREAD.
    pointings: {key: dict or None}; snapshot_of: {key: snapshot key}.
    Returns {key: pointing dict with 'source' = 'stars' or 'snapshot'}.
    """
    out = {}
    by_snap = {}
    for key, p in pointings.items():
        if p and not p.get('unreliable'):
            out[key] = dict(p, source='stars',
                            err=max(p['spread'] / math.sqrt(p['n']), 0.1))
            by_snap.setdefault(snapshot_of.get(key), []).append(p)
    for key, p in pointings.items():
        if key in out:
            continue
        sib = by_snap.get(snapshot_of.get(key))
        if not sib:
            continue
        o = np.array([(s['d_ra'], s['d_dec']) for s in sib])
        med = np.median(o, axis=0)
        spread = float(np.max(np.hypot(*(o - med).T)))
        if spread <= SNAPSHOT_SPREAD:
            out[key] = {'d_ra': float(med[0]), 'd_dec': float(med[1]),
                        'n': 0, 'spread': spread, 'source': 'snapshot',
                        'err': max(spread, 0.5)}
    return out


def corrected_wcs(wcs, pointing, dec):
    """A copy of the image WCS that gives true (Gaia-frame) coordinates:
    the image's coordinates minus the measured offset."""
    w = wcs.deepcopy()
    cra, cdec = w.wcs.crval
    w.wcs.crval = [cra - pointing['d_ra'] / 3600.0
                   / math.cos(math.radians(dec)),
                   cdec - pointing['d_dec'] / 3600.0]
    w.wcs.set()
    return w
