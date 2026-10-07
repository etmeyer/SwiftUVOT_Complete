#!/usr/bin/env python3
"""
swift_uvot_lightcurve.py -- Step 7: the light curve.

Combines the included exposures of the Step 6 master table into one point
per observation (OBSID) and filter, and writes the light curve as a table
and a plot.

Combining. The rates of the included exposures are averaged with
inverse-variance weights. Each exposure's error is its statistical error
plus a systematic floor (--floor, default 1.5 % of the rate), added in
quadrature: included exposures of one 3C 273 observation scatter by 1.0-1.6 %
more than their statistical errors (reduced chi^2 = 1 at V 1.0, B 1.6, U 0.7,
UVW1 1.1, UVM2 1.3, UVW2 1.4 %). U exposures flagged high_coi (probably
about 8 % low) are used only when an observation has no other U exposure.
The time is the exposure-weighted mean of the exposures' midpoints.

Each point gets one status:

    detected      combined rate / error >= 3 (--nsigma)
    upper_limit   below that; the limit is max(rate, 0) + 3 x error
    lower_limit   nothing included, but some exposures were excluded only
                  because they are saturated: the highest capped rate, a
                  lower limit on the true rate
    excluded      exposures were measured but none was included
    not_covered   no exposure had the source in its field

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
from scipy.stats import chi2 as chi2_dist

import swift_uvot_env as uenv
from swift_uvot_tables import read_table, write_table

FILTERS = ('V', 'B', 'U', 'UVW1', 'UVM2', 'UVW2')
DEFAULT_FLOOR = 0.015
DEFAULT_NSIGMA = 3.0
SCATTER_P = 0.001     # chi^2 probability below which a point is flagged
KEPT_FLAGS = ('high_coi', 'jitter_ok', 'discrepant', 'sss_mid', 'sss_high')

COLUMNS = ['obsid', 'filter', 'date_mid', 'mjd_mid', 'mjd_start', 'mjd_stop',
           'status', 'n_used', 'n_total', 'exposure', 'rate', 'rate_err',
           'snr', 'chi2_red', 'mag', 'mag_err', 'ab_mag', 'ab_mag_err',
           'mag_sys', 'flux_aa', 'flux_aa_err', 'flux_mjy', 'flux_mjy_err',
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
    """
    factors = {}
    by_filter = defaultdict(list)
    for r in photometry:
        if r['status'] == 'ok' and r['photflags'] == '-' \
                and r['saturated'] == '0' and _num(r['rate']) > 0:
            by_filter[r['filter']].append(r)
    for filt, rows in by_filter.items():
        vals = defaultdict(list)
        for r in rows[:20]:
            with fits.open(os.path.join(outdir, r['phot_file'])) as hdul:
                d = hdul[1].data[0]
            rate = float(d['CORR_RATE'])
            vals['zp_vega'].append(float(d['MAG']) + 2.5 * math.log10(rate))
            vals['zp_ab'].append(float(d['AB_MAG']) + 2.5 * math.log10(rate))
            vals['fcf_aa'].append(float(d['FLUX_AA']) / rate)
            vals['fcf_mjy'].append(float(d['FLUX_HZ']) / rate)
            vals['mag_sys'].append(float(d['MAG_ERR_SYS']))
        factors[filt] = {k: float(np.median(v)) for k, v in vals.items()}
        spread = np.ptp(vals['zp_vega'])
        if spread > 0.002:
            print('[warn] %s: zero points from different exposures differ by '
                  '%.3f mag' % (filt, spread))
    return factors


def combine(rows, floor, nsigma):
    """One light-curve point (dict) from the master-table rows of one
    observation and filter."""
    used = [r for r in rows if r['include'] == 'yes']
    clean = [r for r in used if 'high_coi' not in r['flags']]
    if clean and len(clean) < len(used):
        used = clean     # biased high_coi U only when nothing else exists
    point = {'n_total': len(rows), 'n_used': len(used)}
    reasons = []
    if used:
        x = np.array([_num(r['rate']) for r in used])
        e = np.hypot([_num(r['rate_err']) for r in used], floor * np.abs(x))
        w = 1.0 / e ** 2
        rate = float(np.sum(w * x) / np.sum(w))
        err = float(1.0 / math.sqrt(np.sum(w)))
        exp = np.array([_num(r['exposure']) for r in used])
        if len(used) > 1:
            c2 = float(np.sum(w * (x - rate) ** 2))
            point['chi2_red'] = '%.2f' % (c2 / (len(used) - 1))
            if chi2_dist.sf(c2, len(used) - 1) < SCATTER_P:
                reasons.append('exposures disagree (chi2 p < %g)' % SCATTER_P)
        snr = rate / err
        point.update(rate=rate, rate_err=err, snr=snr,
                     status='detected' if snr >= nsigma else 'upper_limit')
        timed = used
        flags = sorted({f for r in used for f in r['flags'].split(',')
                        if f.split(':')[0] in KEPT_FLAGS})
        if reasons:
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


def to_units(point, fac, nsigma):
    """Table strings for a point, with magnitudes and flux densities."""
    out = {}
    status = point['status']
    rate, err = point.get('rate'), point.get('rate_err')
    for key in ('mjd_mid', 'mjd_start', 'mjd_stop'):
        if key in point:
            out[key] = '%.5f' % point[key]
    if 'exposure' in point:
        out['exposure'] = '%.1f' % point['exposure']
    if rate is None or fac is None:
        return out
    if status == 'upper_limit':
        value = max(rate, 0.0) + nsigma * err
    else:
        value = rate
    out['rate'] = '%.6g' % rate
    if err == err:
        out['rate_err'] = '%.3g' % err
        out['snr'] = '%.1f' % point['snr']
    out['mag_sys'] = '%.3f' % fac['mag_sys']
    if value > 0:
        out['mag'] = '%.3f' % (fac['zp_vega'] - 2.5 * math.log10(value))
        out['ab_mag'] = '%.3f' % (fac['zp_ab'] - 2.5 * math.log10(value))
        out['flux_aa'] = '%.4e' % (value * fac['fcf_aa'])
        out['flux_mjy'] = '%.5g' % (value * fac['fcf_mjy'])
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
            for status, marker, label in (('lower_limit', '^',
                                           'lower limit (saturated)'),
                                          ('upper_limit', 'v', 'upper limit')):
                pts = [p for p in sel if p['status'] == status]
                if pts:
                    ax.scatter([year(p) for p in pts],
                               [float(p['flux_mjy']) for p in pts],
                               marker=marker, s=22, color=muted, label=label)
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
                             'its rate (default: %.3f)' % DEFAULT_FLOOR)
    parser.add_argument('--nsigma', type=float, default=DEFAULT_NSIGMA,
                        help='detection threshold and upper-limit sigma '
                             '(default: 3)')
    args = parser.parse_args(argv)

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
    factors = conversion_factors(read_table(phot_path), args.outdir)

    groups = defaultdict(list)
    for r in master:
        groups[(r['obsid'], r['filter'])].append(r)
    points = []
    for (obsid, filt), rows in sorted(groups.items()):
        p = combine(rows, args.floor, args.nsigma)
        row = {'obsid': obsid, 'filter': filt, 'status': p['status'],
               'n_used': p['n_used'], 'n_total': p['n_total'],
               'date_mid': p.get('date_mid'), 'chi2_red': p.get('chi2_red'),
               'flags': p.get('flags', '-'), 'reason': p['reason']}
        row.update(to_units(p, factors.get(filt), args.nsigma))
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
        'one row per observation and filter; rate: inverse-variance mean of '
        'the included exposures (counts/s), errors with a %.1f %% floor per '
        'exposure' % (100 * args.floor),
        'status: detected (snr >= %g), upper_limit (mag/flux are the %g-sigma '
        'limit), lower_limit (only saturated exposures: highest capped rate), '
        'excluded, not_covered' % (args.nsigma, args.nsigma),
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
