#!/usr/bin/env python3
"""
swift_uvot_lightcurve.py -- Step 7: the light curve.

Combines the included exposures of the Step 6 master table into one point
per observation (OBSID) and filter, and writes the light curve as a table
and a plot.

Combining. The rates of the included exposures are averaged with weights
from the counts each exposure is expected to have at the combined rate
(not the counts it happened to get, which bias a faint source's mean low;
Mighell 1999): uvotsource's error is rescaled to the expected counts, and
a systematic floor (--floor, default 1.5 % of the combined rate) is added
in quadrature. Included exposures of one 3C 273 observation scatter by
1.0-1.6 % more than their statistical errors (reduced chi^2 = 1 at V 1.0,
B 1.6, U 0.7, UVW1 1.1, UVM2 1.3, UVW2 1.4 %). An additive background
systematic (--bkg-sys, a fraction of the background in the aperture,
measured per filter on blank sky) is added once to the combined error. U
exposures flagged high_coi (probably about 8 % low) are used only when an
observation has no other U exposure. The time is the exposure-weighted
mean of the exposures' midpoints.

Each point gets one status:

    detected      the summed counts are at least 3 sigma (--nsigma) above
                  the background: Poisson probability, with the
                  background's uncertainty, as Gaussian sigma (det_sigma)
    upper_limit   below that; the limit is Kraft, Burrows & Nousek (1991)
                  at 99.87 % (3 sigma) on the summed counts (rate_limit)
    lower_limit   nothing included, but some exposures were excluded only
                  because they are saturated: the highest capped rate, a
                  lower limit on the true rate
    excluded      exposures were measured but none was included
    not_covered   no exposure had the source in its field

The measured rate and error are kept for every point, with the 3-sigma
sensitivity (rate_sens, mag_sens: what a source needs to be detected over
the background alone; not an upper limit).

Magnitudes and flux densities are converted from the combined rate with
uvotsource's own factors (Vega and AB zero points, erg/s/cm^2/A and mJy
per count/s), read from Step 5's outputs. mag_sys is the absolute
calibration (zero-point) error, kept separate from mag_err.

Output (in --outdir):
    uvot_lightcurve.txt    one row per observation and filter
    uvot_lightcurve.pdf    flux density against time, one panel per filter

Usage:
    swift_uvot_lightcurve.py
    swift_uvot_lightcurve.py --floor 0.02
    swift_uvot_lightcurve.py --bkg-sys UVM2=0.03
"""

import argparse
import math
import os
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits
from astropy.time import Time
from scipy.optimize import brentq
from scipy.special import gammainc, gammaincc
from scipy.stats import chi2 as chi2_dist, gamma as gamma_dist, norm

import swift_uvot_env as uenv
from swift_uvot_tables import read_table, write_table

FILTERS = ('V', 'B', 'U', 'UVW1', 'UVM2', 'UVW2')
DEFAULT_FLOOR = 0.015
DEFAULT_NSIGMA = 3.0
SCATTER_P = 0.001     # chi^2 probability below which a point is flagged
KEPT_FLAGS = ('high_coi', 'jitter_ok', 'discrepant', 'sss_mid', 'sss_high')
# Additive background systematic, as a fraction of the background in the
# 5" aperture, added once to each combined point. Measured on blank sky:
# 26 empty random positions in every included PKS 0537-286 exposure, run
# through this combination; the value makes the pulls of the combined
# points unit-width (the larger of the standard-deviation and MAD
# estimates). V, UVW1 and UVM2 need none.
BKG_SYS = {'V': 0.0, 'B': 0.018, 'U': 0.018, 'UVW1': 0.0, 'UVM2': 0.0,
           'UVW2': 0.026}
MAX_ITER = 50         # weight iterations (3C 273 needs at most 5)
N_PRIOR = 40          # quadrature nodes over the background's uncertainty

COLUMNS = ['obsid', 'filter', 'date_mid', 'mjd_mid', 'mjd_start', 'mjd_stop',
           'status', 'n_used', 'n_total', 'exposure', 'counts', 'bkg_counts',
           'det_sigma', 'rate', 'rate_err', 'snr', 'chi2_red', 'rate_limit',
           'rate_sens', 'mag', 'mag_err', 'ab_mag', 'ab_mag_err', 'mag_sys',
           'mag_sens', 'flux_aa', 'flux_aa_err', 'flux_mjy', 'flux_mjy_err',
           'flags', 'reason']


def _num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float('nan')


def conversion_factors(photometry, outdir):
    """
    {filter: dict(zp_vega, zp_ab, fcf_aa, fcf_mjy, mag_sys)} from up to 20
    unflagged Step 5 exposures per filter (uvotsource's own conversion).
    A filter with no positive rate (a faint source in the UV) takes them
    from uvotsource's limit columns instead (CORR_RATE_LIMIT with MAG_LIM,
    AB_MAG_LIM, FLUX_AA_LIM, FLUX_HZ_LIM): the same conversion.
    """
    factors = {}
    by_filter, spare = defaultdict(list), defaultdict(list)
    for r in photometry:
        if r['status'] != 'ok':
            continue
        if r['photflags'] == '-' and r['saturated'] == '0' \
                and _num(r['rate']) > 0:
            by_filter[r['filter']].append(r)
        else:
            spare[r['filter']].append(r)
    for filt in set(by_filter) | set(spare):
        rows = by_filter.get(filt)
        cols = ('CORR_RATE', 'MAG', 'AB_MAG', 'FLUX_AA', 'FLUX_HZ')
        if not rows:
            rows = spare[filt]
            cols = ('CORR_RATE_LIMIT', 'MAG_LIM', 'AB_MAG_LIM', 'FLUX_AA_LIM',
                    'FLUX_HZ_LIM')
        vals = defaultdict(list)
        for r in rows[:20]:
            with fits.open(os.path.join(outdir, r['phot_file'])) as hdul:
                d = hdul[1].data[0]
            rate = float(d[cols[0]])
            if not rate > 0:
                continue
            vals['zp_vega'].append(float(d[cols[1]]) + 2.5 * math.log10(rate))
            vals['zp_ab'].append(float(d[cols[2]]) + 2.5 * math.log10(rate))
            vals['fcf_aa'].append(float(d[cols[3]]) / rate)
            vals['fcf_mjy'].append(float(d[cols[4]]) / rate)
            sys_err = float(d['MAG_ERR_SYS'])
            if sys_err == sys_err:
                vals['mag_sys'].append(sys_err)
        if not vals['zp_vega']:
            continue
        factors[filt] = {k: float(np.median(v)) if v else float('nan')
                         for k, v in vals.items()}
        factors[filt].setdefault('mag_sys', float('nan'))
        spread = np.ptp(vals['zp_vega'])
        if spread > 0.002:
            print('[warn] %s: zero points from different exposures differ by '
                  '%.3f mag' % (filt, spread))
        if cols[0] == 'CORR_RATE_LIMIT':
            print('[info] %s: no positive rate; conversion factors from '
                  "uvotsource's limit columns" % filt)
    return factors


def exposure_inputs(r, p):
    """The numbers the combination needs for one included exposure: its
    rate and error (master table) and its counts, areas and corrections
    (Step 5 table)."""
    sss = _num(p['sss'])
    return dict(
        x=_num(r['rate']), e=_num(r['rate_err']), T=_num(r['exposure']),
        N=_num(p['src_cnts']), M=_num(p['bkg_cnts']),
        a_src=_num(p['src_area']), a_bkg=_num(p['bkg_area']),
        # corrected rate per raw count rate: coincidence loss, sensitivity
        # loss over time, small-scale sensitivity (-99.9 = not known),
        # large-scale sensitivity
        K=_num(p['coi']) * _num(p['senscorr']) * (sss if sss > 0 else 1.0)
        / _num(p['lss']))


def model_weights(ex, floor):
    """
    Weighted mean rate with weights from the counts each exposure is
    EXPECTED to have at the mean rate, not the counts it happened to get.

    Weights 1/error^2 with each exposure's error from its own counts give
    exposures that fluctuated low a smaller error and so a larger weight:
    the mean comes out low (about -1 count per exposure; Mighell 1999).
    Here the expected counts in the aperture are E_i = R T_i / K_i + B_i
    (B_i: the background expected from the observation's pooled background
    per arcsec^2 per s), and uvotsource's error is rescaled to them:
    sigma_i^2 = e_i^2 (E_i + V_i) / (max(N_i, 1) + V_i) + (floor R)^2,
    V_i the variance of the background subtracted. In uvotsource's
    coincidence-loss regime E_i ~ N_i, so its own error is kept. The floor
    is a fraction of the combined rate R, not of each exposure's rate
    (f |x_i| favours exposures that read low). Iterated from the
    exposure-weighted mean. Returns (R, error, chi2, converged).
    """
    x = np.array([d['x'] for d in ex])
    e = np.array([d['e'] for d in ex])
    T = np.array([d['T'] for d in ex])
    N = np.array([d['N'] for d in ex])
    M = np.array([d['M'] for d in ex])
    a_src = np.array([d['a_src'] for d in ex])
    a_bkg = np.array([d['a_bkg'] for d in ex])
    K = np.array([d['K'] for d in ex])
    B = M.sum() / np.sum(a_bkg * T) * a_src * T
    V = (a_src / a_bkg) ** 2 * np.maximum(M, 1.0)
    R = float(np.sum(T * x) / np.sum(T))
    converged = False
    for _ in range(MAX_ITER):
        Rp = max(R, 0.0)
        E = Rp * T / K + B
        w = 1.0 / (e ** 2 * (E + V) / (np.maximum(N, 1.0) + V)
                   + (floor * Rp) ** 2)
        R_new = float(np.sum(w * x) / np.sum(w))
        step, R = abs(R_new - R), R_new
        if step <= 1e-6 / math.sqrt(np.sum(w)):
            converged = True
            break
    if not converged:     # exposure weighting: unbiased, a little less precise
        R = float(np.sum(T * x) / np.sum(T))
        w = 1.0 / (e ** 2 + (floor * max(R, 0.0)) ** 2)
        err = float(math.sqrt(np.sum(T ** 2 / w)) / np.sum(T))
        return R, err, float(np.sum(w * (x - R) ** 2)), False
    return (R, float(1.0 / math.sqrt(np.sum(w))),
            float(np.sum(w * (x - R) ** 2)), True)


def background_prior(b_hat, var):
    """Nodes and weights over the background counts in the aperture: a
    gamma distribution with mean b_hat and variance var (the background is
    estimated, not known), by Gauss-Legendre in its quantiles."""
    if var <= 0 or b_hat <= 0:
        return np.array([max(b_hat, 0.0)]), np.array([1.0])
    u, w = np.polynomial.legendre.leggauss(N_PRIOR)
    shape, scale = b_hat ** 2 / var, var / b_hat
    return gamma_dist.ppf((u + 1) / 2, shape, scale=scale), w / 2


def detection_sigma(n, b_nodes, b_w):
    """Gaussian-equivalent significance of n summed counts over the
    background: the probability that background alone gives >= n counts
    (Poisson, continuous in n since the counts are sub-pixel sums,
    averaged over the background's uncertainty), as sigma."""
    if n <= 0:
        return -np.inf
    p = float(np.sum(b_w * gammainc(n, np.maximum(b_nodes, 1e-300))))
    return float(norm.isf(p)) if p > 0 else np.inf


def kbn_limit(n, b_nodes, b_w, cl):
    """Upper limit on the source counts at confidence cl (Kraft, Burrows &
    Nousek 1991: flat prior on the source counts >= 0), for n summed counts
    (continuous form) and a background known to its distribution: the
    posterior's CDF is sum_j w_j [Q(n+1, B_j) - Q(n+1, B_j + s)] /
    sum_j w_j Q(n+1, B_j), Q the regularised upper incomplete gamma."""
    n = max(n, 0.0)
    q0 = b_w * gammaincc(n + 1, b_nodes)
    if q0.sum() <= 0:       # far fewer counts than background: no limit
        return float('nan')

    def cdf(s):
        return np.sum(q0 - b_w * gammaincc(n + 1, b_nodes + s)) / q0.sum()

    hi = n + 50 * math.sqrt(n + float(np.max(b_nodes)) + 1) + 50
    while cdf(hi) < cl:
        hi *= 2
    return float(brentq(lambda s: cdf(s) - cl, 0.0, hi, xtol=1e-6))


def combine(rows, phot, floor, nsigma, eps):
    """One light-curve point (dict) from the master-table rows of one
    observation and filter. phot: Step 5 rows by (obsid, extname); eps:
    additive background systematic (fraction of the background in the
    aperture) for this filter."""
    used = [r for r in rows if r['include'] == 'yes']
    clean = [r for r in used if 'high_coi' not in r['flags']]
    if clean and len(clean) < len(used):
        used = clean     # biased high_coi U only when nothing else exists
    point = {'n_total': len(rows), 'n_used': len(used)}
    reasons = []
    if used:
        ex = [exposure_inputs(r, phot[(r['obsid'], r['extname'])])
              for r in used]
        rate, err_stat, c2, converged = model_weights(ex, floor)
        if not converged:
            reasons.append('weights did not converge: exposure-weighted mean')
        # summed counts: what detection and the upper limit are decided on
        n = sum(d['N'] for d in ex)
        b_hat = sum(d['M'] * d['a_src'] / d['a_bkg'] for d in ex)
        b_var = sum(d['M'] * (d['a_src'] / d['a_bkg']) ** 2 for d in ex)
        t_eff = sum(d['T'] / d['K'] for d in ex)
        # the background systematic does not average down within an
        # observation: it is added once, to the combined point
        b_sys = eps * b_hat
        err = math.hypot(err_stat, b_sys / t_eff)
        nodes, wts = background_prior(b_hat, b_var + b_sys ** 2)
        det = detection_sigma(n, nodes, wts)
        if len(used) > 1:
            point['chi2_red'] = '%.2f' % (c2 / (len(used) - 1))
            if chi2_dist.sf(c2, len(used) - 1) < SCATTER_P:
                reasons.append('exposures disagree (chi2 p < %g)' % SCATTER_P)
        point.update(rate=rate, rate_err=err, snr=rate / err, det_sigma=det,
                     counts=n, bkg_counts=b_hat,
                     rate_sens=nsigma * math.sqrt(b_hat + b_var + b_sys ** 2)
                     / t_eff,
                     status='detected' if det >= nsigma else 'upper_limit')
        if point['status'] == 'upper_limit':
            s = kbn_limit(n, nodes, wts, norm.cdf(nsigma))
            point['rate_limit'] = s / t_eff
        timed = used
        flags = sorted({f for r in used for f in r['flags'].split(',')
                        if f.split(':')[0] in KEPT_FLAGS})
        if reasons and reasons[-1].startswith('exposures disagree'):
            flags.append('scatter')
        point['flags'] = ','.join(flags) or '-'
    else:
        sat = [r for r in rows if r['rules'] == 'saturated']
        measured = [r for r in rows if r['rate'] not in (None, '-')]
        if sat:
            best = max(sat, key=lambda r: _num(r['rate']))
            point.update(status='lower_limit', rate=_num(best['rate']),
                         rate_err=float('nan'))
            timed = sat
            reasons.append('none included; %d exposure(s) excluded only for '
                           'saturation' % len(sat))
        elif measured:
            point['status'] = 'excluded'
            timed = []
            rules = sorted({x.split(':')[0] for r in rows
                            for x in r['rules'].split(',') if x != '-'})
            reasons.append('none included (%s)' % ', '.join(rules))
        else:
            point['status'] = 'not_covered'
            timed = []
            reasons.append(', '.join(sorted({r['rules'] for r in rows})))
    src = timed or rows
    exp = np.array([max(_num(r['exposure']), 1e-3) for r in src])
    mjd = np.array([_num(r['mjd_mid']) for r in src])
    ok = np.isfinite(mjd)
    if ok.any():
        mid = float(np.sum(exp[ok] * mjd[ok]) / np.sum(exp[ok]))
        point['mjd_mid'] = mid
        point['mjd_start'] = float(np.min(mjd[ok] - exp[ok] / 2 / 86400.0))
        point['mjd_stop'] = float(np.max(mjd[ok] + exp[ok] / 2 / 86400.0))
        point['date_mid'] = Time(mid, format='mjd').isot[:19]
    if timed:
        point['exposure'] = float(np.sum([_num(r['exposure']) for r in timed]))
    point['reason'] = '; '.join(reasons)
    return point


def to_units(point, fac):
    """Table strings for a point, with magnitudes and flux densities: of
    the rate if detected or a lower limit, of the upper limit otherwise."""
    out = {}
    status = point['status']
    rate, err = point.get('rate'), point.get('rate_err')
    for key in ('mjd_mid', 'mjd_start', 'mjd_stop'):
        if key in point:
            out[key] = '%.5f' % point[key]
    if 'exposure' in point:
        out['exposure'] = '%.1f' % point['exposure']
    if 'counts' in point:
        out['counts'] = '%.1f' % point['counts']
        out['bkg_counts'] = '%.1f' % point['bkg_counts']
        d = point['det_sigma']
        out['det_sigma'] = '%.1f' % d if np.isfinite(d) else \
            ('>37' if d > 0 else '-')
    if rate is None:
        return out
    out['rate'] = '%.6g' % rate
    if err == err:
        out['rate_err'] = '%.3g' % err
        out['snr'] = '%.1f' % point['snr']
    if point.get('rate_limit', float('nan')) == point.get('rate_limit'):
        out['rate_limit'] = '%.4g' % point['rate_limit']
    if 'rate_sens' in point:
        out['rate_sens'] = '%.4g' % point['rate_sens']
    if fac is None:
        return out
    value = point.get('rate_limit') if status == 'upper_limit' else rate
    out['mag_sys'] = '%.3f' % fac['mag_sys']
    if value is not None and value > 0:
        out['mag'] = '%.3f' % (fac['zp_vega'] - 2.5 * math.log10(value))
        out['ab_mag'] = '%.3f' % (fac['zp_ab'] - 2.5 * math.log10(value))
        out['flux_aa'] = '%.4e' % (value * fac['fcf_aa'])
        out['flux_mjy'] = '%.5g' % (value * fac['fcf_mjy'])
    if point.get('rate_sens', 0) > 0:
        out['mag_sens'] = '%.3f' % (fac['zp_vega']
                                    - 2.5 * math.log10(point['rate_sens']))
    if status == 'detected':
        merr = 2.5 / math.log(10) * err / rate
        out['mag_err'] = out['ab_mag_err'] = '%.3f' % merr
        out['flux_aa_err'] = '%.2e' % (err * fac['fcf_aa'])
        out['flux_mjy_err'] = '%.3g' % (err * fac['fcf_mjy'])
    return out


def plot_lightcurve(pdf_path, points, src_label):
    """Flux density against time, one panel per filter."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    ink, muted = '#1f4e79', '#7f8c99'
    filters = [f for f in FILTERS
               if any(p['filter'] == f and p.get('flux_mjy') not in (None, '-')
                      for p in points)]
    if not filters:
        return
    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(len(filters), 1, sharex=True,
                                 figsize=(11, 2.0 * len(filters) + 1.2))
        axes = np.atleast_1d(axes)
        for ax, f in zip(axes, filters):
            sel = [p for p in points if p['filter'] == f
                   and p.get('flux_mjy') not in (None, '-')]
            year = lambda p: Time(float(p['mjd_mid']), format='mjd').decimalyear
            det = [p for p in sel if p['status'] == 'detected']
            for flagged in (False, True):
                pts = [p for p in det if ('high_coi' in p['flags']) == flagged]
                if pts:
                    ax.errorbar([year(p) for p in pts],
                                [float(p['flux_mjy']) for p in pts],
                                yerr=[float(p['flux_mjy_err']) for p in pts],
                                fmt='o', ms=4, lw=0.8, color=ink,
                                mfc='white' if flagged else ink,
                                label='high_coi (U, ~8 % low)' if flagged
                                else 'detected')
            # limits as arrows from the limit, a quarter of its value long
            for status, label in (('lower_limit', 'lower limit (saturated)'),
                                  ('upper_limit', 'upper limit (KBN)')):
                pts = [p for p in sel if p['status'] == status]
                if pts:
                    y = np.array([float(p['flux_mjy']) for p in pts])
                    up = status == 'upper_limit'
                    ax.errorbar([year(p) for p in pts], y, yerr=0.25 * y,
                                uplims=up, lolims=not up, fmt='_', ms=6,
                                lw=0.8, color=muted, label=label)
            ax.set_ylabel('%s (mJy)' % f)
            ax.grid(axis='y', color='#e3e6ea', lw=0.6)
            for side in ('top', 'right'):
                ax.spines[side].set_visible(False)
            ax.legend(fontsize=7, loc='upper right', frameon=False)
        axes[-1].set_xlabel('year')
        fig.suptitle('UVOT light curve (%s): one point per observation and '
                     'filter' % src_label, fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        pdf.savefig(fig)
        plt.close(fig)


def main(argv=None):
    sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser(
        description='Step 7: the light curve.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--outdir', default='UVOT_output',
                        help='Steps 3-6 output folder (default: UVOT_output)')
    parser.add_argument('--floor', type=float, default=DEFAULT_FLOOR,
                        help='systematic error per exposure, as a fraction of '
                             'the combined rate (default: %.3f)'
                             % DEFAULT_FLOOR)
    parser.add_argument('--nsigma', type=float, default=DEFAULT_NSIGMA,
                        help='detection threshold in Gaussian sigma; also the '
                             'upper-limit confidence (3 = 99.87 %%) and the '
                             'sensitivity (default: 3)')
    parser.add_argument('--bkg-sys', default=None, metavar='VALUES',
                        help='additive background systematic as a fraction '
                             'of the background in the aperture: one number '
                             'for all filters, or e.g. V=0.01,UVM2=0.05 '
                             '(others keep their default); 0 to switch off '
                             '(default: %s)' % ','.join(
                                 '%s=%g' % kv for kv in BKG_SYS.items()))
    args = parser.parse_args(argv)
    if not 1.0 <= args.nsigma <= 8.0:
        print('ERROR: --nsigma must be between 1 and 8', file=sys.stderr)
        return 1
    bkg_sys = dict(BKG_SYS)
    if args.bkg_sys:
        try:
            if '=' in args.bkg_sys:
                for item in args.bkg_sys.split(','):
                    f, v = item.split('=')
                    if f.strip().upper() not in BKG_SYS:
                        raise ValueError(f)
                    bkg_sys[f.strip().upper()] = float(v)
            else:
                bkg_sys = {f: float(args.bkg_sys) for f in BKG_SYS}
        except ValueError:
            print('ERROR: --bkg-sys %s: give a number or FILTER=VALUE pairs '
                  '(filters %s)' % (args.bkg_sys, ', '.join(FILTERS)),
                  file=sys.stderr)
            return 1

    master_path = os.path.join(args.outdir, 'uvot_master_table.txt')
    phot_path = os.path.join(args.outdir, 'uvot_photometry.txt')
    for path, step in ((master_path, 'Step 6 (make_uvot_master_table.py)'),
                       (phot_path, 'Step 5 (swift_uvot_photometry.py)')):
        if not os.path.exists(path):
            print('ERROR: %s not found; run %s first.' % (path, step),
                  file=sys.stderr)
            return 1
    master = read_table(master_path)
    with open(master_path) as fh:
        header = ''.join(line for line in fh if line.startswith('#'))
    m = re.search(r'source: (RA [\d.+-]+ Dec [\d.+-]+)', header)
    src_label = m.group(1) if m else '?'
    photometry = read_table(phot_path)
    phot = {(r['obsid'], r['extname']): r for r in photometry}
    factors = conversion_factors(photometry, args.outdir)

    groups = defaultdict(list)
    for r in master:
        groups[(r['obsid'], r['filter'])].append(r)
    points = []
    for (obsid, filt), rows in sorted(groups.items()):
        p = combine(rows, phot, args.floor, args.nsigma,
                    bkg_sys.get(filt, 0.0))
        row = {'obsid': obsid, 'filter': filt, 'status': p['status'],
               'n_used': p['n_used'], 'n_total': p['n_total'],
               'date_mid': p.get('date_mid'), 'chi2_red': p.get('chi2_red'),
               'flags': p.get('flags', '-'), 'reason': p['reason']}
        row.update(to_units(p, factors.get(filt)))
        points.append(row)
    points.sort(key=lambda r: (FILTERS.index(r['filter'])
                               if r['filter'] in FILTERS else 99,
                               _num(r.get('mjd_mid'))))

    prov = uenv.provenance()
    comments = [
        'swift_uvot_lightcurve.py, %s' % datetime.now(timezone.utc)
        .strftime('%Y-%m-%dT%H:%M:%SZ'),
        'command: %s' % ' '.join(sys.argv),
        'source: %s' % src_label,
        'pipeline %s' % prov['PIPECOMM'],
        'one row per observation and filter; rate: weighted mean of the '
        'included exposures (counts/s), weights from the expected counts, '
        'floor %.1f %% of the rate; rate_err includes the background '
        'systematic (%s)' % (100 * args.floor, ','.join(
            '%s=%g' % (f, bkg_sys[f]) for f in FILTERS if f in bkg_sys)),
        'counts, bkg_counts: summed raw counts in the aperture and the '
        'background expected there; det_sigma: Poisson significance of '
        'counts over bkg_counts, as Gaussian sigma',
        'status: detected (det_sigma >= %g), upper_limit (mag/flux are '
        'rate_limit: Kraft-Burrows-Nousek at %.2f %%), lower_limit (only '
        'saturated exposures: highest capped rate), excluded, not_covered'
        % (args.nsigma, 100 * norm.cdf(args.nsigma)),
        'rate_sens, mag_sens: the rate a source needs to be %g sigma above '
        'the background alone (sensitivity, not an upper limit)'
        % args.nsigma,
        'mag/mag_err Vega, ab_mag AB; mag_sys: absolute calibration error '
        '(not in mag_err); flux_aa erg/s/cm^2/A, flux_mjy mJy; times MJD (UTC)',
    ]
    for f in FILTERS:
        if f in factors:
            fac = factors[f]
            comments.append('%s: Vega zero point %.3f, AB %.3f, %.4e '
                            'erg/s/cm^2/A and %.4g mJy per count/s'
                            % (f, fac['zp_vega'], fac['zp_ab'], fac['fcf_aa'],
                               fac['fcf_mjy']))
    lc_path = os.path.join(args.outdir, 'uvot_lightcurve.txt')
    write_table(lc_path, points, COLUMNS, comments)
    pdf_path = os.path.join(args.outdir, 'uvot_lightcurve.pdf')
    plot_lightcurve(pdf_path, points, src_label)

    print('[info] source %s; %d exposures in %d observation/filter pairs'
          % (src_label, len(master), len(points)))
    print()
    print('[summary] filter  detected  upper_limit  lower_limit  excluded  '
          'not_covered  high_coi  scatter')
    for f in FILTERS:
        sel = [p for p in points if p['filter'] == f]
        if sel:
            n = {s: sum(p['status'] == s for p in sel) for s in
                 ('detected', 'upper_limit', 'lower_limit', 'excluded',
                  'not_covered')}
            print('          %-6s  %8d  %11d  %11d  %8d  %11d  %8d  %7d' % (
                f, n['detected'], n['upper_limit'], n['lower_limit'],
                n['excluded'], n['not_covered'],
                sum('high_coi' in p['flags'] for p in sel),
                sum('scatter' in p['flags'] for p in sel)))
    print('[summary] Written: %s, %s' % (lc_path, pdf_path))
    return 0


if __name__ == '__main__':
    sys.exit(main())
