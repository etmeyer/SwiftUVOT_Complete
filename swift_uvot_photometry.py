#!/usr/bin/env python3
"""
swift_uvot_photometry.py -- Step 5: uvotsource photometry of every exposure.

For every exposure Step 4 passed (status ok), runs uvotsource on that one
extension, with Step 4's source and background regions and the exposure
map, and keeps everything it reports. Nothing is excluded here: Step 6
applies the quality rules to these numbers, so they can be changed and
re-applied without running HEASoft again.

The uvotsource settings, and why:

    image=<sky image>[EXTNAME]     one exposure per call (given no
    expfile=<exposure map>[EXTNAME]  extension, uvotsource measures only
                                   the first one in the file)
    srcreg, bkgreg                 Step 4's regions
    apercorr=NONE                  the 5" aperture is the calibrated one
    sigma=3                        for the background-limited magnitude
    syserr=no                      statistical errors; the systematic
                                   error is kept in its own column
    forcephot=no, skipbad=no       uvotsource's own checks stay on, and a
                                   refusal is a failure with its reason
    centroid=no                    Step 4 placed the regions
    history=no                     fthedit can crash writing history
                                   keywords with long CALDB paths

Each call runs in its own folder (uvotsource appends to an existing output
table); its FITS table and log are kept as photometry/<OBSID>/<EXTNAME>.fits
and .log.

Recorded per exposure (uvot_photometry.txt): the raw counts and areas, the
counts per frame in the 5" circle (uvotsource's coincidence-loss input;
above 0.98 it clamps the rate and sets SATURATED, so the rate is a lower
limit), every correction factor, the corrected rate and its error, the 3σ
background-limited rate and magnitude, the magnitude and flux density, and
uvotsource's photometry flags. The rate is uvotsource's fully corrected
rate (SENSCORR_RATE). On a LOW low-sensitivity patch uvotsource still
computes it but reports no magnitude (MAG 99, CORR_RATE -999), so Step 6
can still judge such exposures. Exposures Step 4 did not pass keep their
Step 4 status.

Also written: uvot_photometry.pdf, the rate of every exposure against time
per filter, marked by what Step 6 will look at (saturated, low-sensitivity
patch, trailed or smeared, faint at the source position).

Needs HEASoft (uvotsource) and the UVOT CALDB; --plot-only needs neither.

Exit status: 1 if any uvotsource call failed.

Usage:
    swift_uvot_photometry.py --nproc 16
    swift_uvot_photometry.py --plot-only
"""

import argparse
import os
import re
import sys
from collections import Counter
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits

import swift_uvot_env as uenv
from swift_uvot_inventory import expmap_path
from swift_uvot_runner import run_many
from swift_uvot_tables import read_table, write_table

# uvotsource's frame-count clamp and dead time per frame (UVOT::Source:
# FRAME_COUNT_LIMIT; VERTICAL_TRANSFER_TIME_s x VERTICAL_PIXELS_PER_IMAGE)
FRAME_COUNT_LIMIT = 0.98
DEAD_TIME_S = 6e-7 * 290

# PHOTFLAG bits in UVOT::LCPar order. uvotsource itself sets
# NO_EXPOSURE_MAP, NO_QUALITY_MAP, UNEVEN_EXPOSURE, BAD_SSS and BAD_LSS.
PHOTFLAG_NAMES = ('NO_EXPOSURE_MAP', 'NO_QUALITY_MAP', 'NO_ASPECT_CORRECTION',
                  'UNUSED', 'SOURCE_CONFUSION', 'SMALL_BKG_AREA',
                  'SHORT_EXPOSURE', 'UNEVEN_EXPOSURE', 'BIG_COI_FACTOR',
                  'BAD_SSS', 'BAD_LSS')
# always set: the archive has no quality maps
IGNORED_FLAGS = ('NO_QUALITY_MAP',)

UVOTSOURCE_PARAMS = {'sigma': 3, 'apercorr': 'NONE', 'syserr': 'no',
                     'forcephot': 'no', 'skipbad': 'no', 'centroid': 'no',
                     'history': 'no', 'chatter': 1}

# Step 4 limits used to mark exposures in the plot (Step 6 decides)
TRAIL_AXIS_RATIO = 1.3
FILTERS = ('V', 'B', 'U', 'UVW1', 'UVM2', 'UVW2')

COLUMNS = ['obsid', 'filter', 'extname', 'exposure', 'frametime',
           'counts_frame', 'saturated', 'src_cnts', 'bkg_cnts', 'src_area',
           'bkg_area', 'coi', 'coi_bkg', 'lss', 'sss', 'senscorr', 'rate',
           'rate_err', 'rate_lim', 'nsigma', 'mag', 'mag_err', 'mag_sys',
           'mag_lim', 'flux_aa', 'flux_aa_err', 'photflags', 'phot_file',
           'status', 'reason']


def counts_per_frame(raw_rate, frametime):
    """
    Raw counts per frame in the 5" circle, as uvotsource computes it for
    the coincidence-loss correction (rate x frame time x live fraction).
    """
    return raw_rate * (frametime - DEAD_TIME_S)


def photflag_names(value):
    """Names of the PHOTFLAG bits set, without the always-set ones."""
    return [name for i, name in enumerate(PHOTFLAG_NAMES)
            if value & (1 << i) and name not in IGNORED_FLAGS]


def _warnings(text):
    """Distinct warning lines from a tool's output."""
    seen = []
    for line in text.splitlines():
        line = line.strip()
        if 'warning' in line.lower() and line not in seen:
            seen.append(line)
    return seen


def _errors(text):
    """uvotsource's own error lines (it prints 'error: ...' on stdout)."""
    return [line.strip() for line in text.splitlines()
            if line.strip().lower().startswith('error:')
            and 'no such process' not in line.lower()]


def harvest(path, extname):
    """The one row of a uvotsource output table, checked."""
    with fits.open(path) as hdul:
        data = hdul[1].data
        n = 0 if data is None else len(data)
        if n != 1:
            raise ValueError('%d rows in the output table' % n)
        row = {name: data[name][0] for name in data.columns.names}
    if str(row['EXTNAME']).strip() != extname:
        raise ValueError('output is for %s' % str(row['EXTNAME']).strip())
    return row


def photometry_row(r, step4):
    """Table values (strings) and notes from one uvotsource output row."""
    frametime = float(r['FRAMTIME'])
    cpf = counts_per_frame(float(r['RAW_STD_RATE']), frametime)
    flags = photflag_names(int(r['PHOTFLAG']))
    out = {
        'exposure': '%.2f' % r['EXPOSURE'],
        'frametime': '%.6f' % frametime,
        'counts_frame': '%.3f' % cpf,
        'saturated': int(r['SATURATED']),
        'src_cnts': '%.2f' % r['RAW_TOT_CNTS'],
        'bkg_cnts': '%.2f' % r['RAW_BKG_CNTS'],
        'src_area': '%.2f' % r['SRC_AREA'],
        'bkg_area': '%.1f' % r['BKG_AREA'],
        'coi': '%.4f' % r['COI_STD_FACTOR'],
        'coi_bkg': '%.4f' % r['COI_BKG_FACTOR'],
        'lss': '%.4f' % r['LSS_FACTOR'],
        'sss': '%.4g' % r['SSS_FACTOR'],
        'senscorr': '%.4f' % r['SENSCORR_FACTOR'],
        'rate': '%.6g' % r['SENSCORR_RATE'],
        'rate_err': '%.3g' % r['SENSCORR_RATE_ERR'],
        'rate_lim': '%.3g' % r['CORR_RATE_LIMIT'],
        'nsigma': '%.1f' % r['NSIGMA'],
        'mag': '%.3f' % r['MAG'],
        'mag_err': '%.3f' % r['MAG_ERR_STAT'],
        'mag_sys': '%.3f' % r['MAG_ERR_SYS'],
        'mag_lim': '%.2f' % r['MAG_LIM'],
        'flux_aa': '%.4e' % r['FLUX_AA'],
        'flux_aa_err': '%.2e' % r['FLUX_AA_ERR_STAT'],
        'photflags': ','.join(flags) or '-',
    }
    notes = []
    if int(r['SATURATED']):
        notes.append('saturated (%.3f counts/frame): the rate is a lower '
                     'limit' % cpf)
    if 'BAD_SSS' in flags:
        notes.append('on a LOW low-sensitivity patch: uvotsource gives no '
                     'magnitude')
    if 'BAD_LSS' in flags:
        notes.append('no large-scale sensitivity value: uvotsource gives no '
                     'magnitude')
    # The same lookup as Step 4's, at the same position: they must agree.
    s4 = step4.get('sss_low')
    if s4 not in (None, '-') and (float(s4) <= 0) != (float(r['SSS_FACTOR'])
                                                       <= 0):
        notes.append('low-sensitivity flag differs from Step 4 (%s vs %s)'
                     % (s4, out['sss']))
    return out, notes


def run_photometry(rows, inv, outdir, nproc):
    """uvotsource on every row; returns the results in the same order."""
    calls = []
    for r in rows:
        sky = inv[(r['obsid'], r['filter'], r['extname'])]['file']
        ex = expmap_path(sky)
        suffix = '.img.gz' if sky.endswith('.gz') else '.img'
        dest = os.path.join(outdir, 'photometry', r['obsid'],
                            r['extname'] + '.fits')
        params = dict(UVOTSOURCE_PARAMS)
        params.update({'image': 'sky%s[%s]' % (suffix, r['extname']),
                       'expfile': 'ex%s[%s]' % (suffix, r['extname']),
                       'srcreg': 'src.reg', 'bkgreg': 'bkg.reg',
                       'outfile': 'phot.fits'})
        calls.append(dict(
            tool='uvotsource', params=params,
            inputs={'sky' + suffix: sky, 'ex' + suffix: ex,
                    'src.reg': os.path.join(outdir, r['src_reg']),
                    'bkg.reg': os.path.join(outdir, r['bkg_reg'])},
            outputs={'phot.fits': dest}, timeout=300,
            log=dest.replace('.fits', '.log')))
    step = max(1, len(calls) // 10)

    def progress(done, total):
        if done % step == 0 or done == total:
            print('[info] %d/%d exposures measured' % (done, total))

    return run_many(calls, nproc=nproc, progress=progress)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def _num(row, key):
    v = row.get(key)
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


def plot_photometry(pdf_path, rows, pos, inv, src_label):
    """Per-filter rate against time, marked by what Step 6 will judge."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from astropy.time import Time
    from swift_uvot_positions import CONC_MIN, MIN_SNR

    good = [r for r in rows if r['status'] == 'ok']
    for r in good:
        key = (r['obsid'], r['filter'], r['extname'])
        p, i = pos[key], inv[key]
        r['_year'] = Time(float(i['mjd_mid']), format='mjd').decimalyear
        # one category each, the one that matters most for the rate first
        if int(r['saturated']):
            r['_cat'] = 'saturated (rate is a lower limit)'
        elif _num(p, 'conc') < CONC_MIN or _num(p, 'axis_ratio') > \
                TRAIL_AXIS_RATIO:
            r['_cat'] = 'trailed or smeared (Step 4)'
        elif not _num(p, 'snr') >= MIN_SNR:
            r['_cat'] = 'faint at the position (Step 4 S/N < %.0f)' % MIN_SNR
        elif 'BAD_SSS' in r['photflags']:
            r['_cat'] = 'on a LOW low-sensitivity patch'
        else:
            r['_cat'] = 'none of these'
    cats = [('none of these', 'C0', 'o'),
            ('trailed or smeared (Step 4)', 'C4', 's'),
            ('faint at the position (Step 4 S/N < %.0f)' % MIN_SNR, 'k', 'v'),
            ('on a LOW low-sensitivity patch', 'C1', 'x'),
            ('saturated (rate is a lower limit)', 'C3', '^')]

    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
        ax = axes[0]
        bins = np.linspace(0, 1.0, 51)
        for f in FILTERS:
            v = [_num(r, 'counts_frame') for r in good if r['filter'] == f]
            if v:
                ax.hist(np.clip(v, 0, 0.999), bins, histtype='step', lw=1.3,
                        label='%s (%d)' % (f, len(v)))
        ax.axvline(FRAME_COUNT_LIMIT, color='k', ls=':', lw=0.8)
        ax.set_xlabel('raw counts per frame in the 5" circle '
                      '(dotted: saturation, %.2f)' % FRAME_COUNT_LIMIT)
        ax.set_ylabel('exposures')
        ax.legend(fontsize=8)
        ax = axes[1]
        x = np.arange(len(FILTERS))
        bottom = np.zeros(len(FILTERS))
        for name, color, _ in cats:
            n = np.array([sum(1 for r in good if r['filter'] == f
                              and r['_cat'] == name) for f in FILTERS])
            ax.bar(x, n, bottom=bottom, color=color, label=name)
            bottom += n
        ax.set_xticks(x)
        ax.set_xticklabels(FILTERS)
        ax.set_ylabel('exposures')
        ax.legend(fontsize=7)
        fig.suptitle('Step 5: photometry of %d exposures (%s)'
                     % (len(good), src_label))
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        pdf.savefig(fig)
        plt.close(fig)

        for f in FILTERS:
            sel = [r for r in good if r['filter'] == f]
            if not sel:
                continue
            fig, ax = plt.subplots(figsize=(11, 5.5))
            for name, color, marker in cats:
                pts = [r for r in sel if r['_cat'] == name]
                if not pts:
                    continue
                t = [r['_year'] for r in pts]
                y = [_num(r, 'rate') for r in pts]
                label = '%s (%d)' % (name, len(pts))
                if name.startswith('saturated'):
                    ax.scatter(t, y, s=14, color=color, marker=marker,
                               label=label, zorder=3)
                else:
                    ax.errorbar(t, y, yerr=[_num(r, 'rate_err') for r in pts],
                                fmt=marker, ms=3.5, color=color, lw=0.6,
                                mfc='none' if marker in 'sv' else color,
                                label=label, zorder=2)
            ys = np.array([_num(r, 'rate') for r in sel])
            ys = ys[np.isfinite(ys)]
            if ys.size:
                lo, hi = np.percentile(ys, [0.5, 99.5])
                if hi > 0:
                    ax.set_ylim(min(0.0, lo * 1.1), hi * 1.15)
            ax.set_xlabel('year')
            ax.set_ylabel('corrected count rate (counts/s)')
            ax.set_title('%s: rate of every exposure (Step 5; nothing '
                         'excluded yet)' % f)
            ax.legend(fontsize=8)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser(
        description='Step 5: uvotsource photometry of every exposure.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--outdir', default='UVOT_output',
                        help='Steps 3-4 output folder (default: UVOT_output)')
    parser.add_argument('--nproc', type=int, default=8,
                        help='parallel uvotsource calls (default: 8)')
    parser.add_argument('--plot-only', action='store_true',
                        help='only remake uvot_photometry.pdf from the tables')
    args = parser.parse_args(argv)

    pos_path = os.path.join(args.outdir, 'uvot_positions.txt')
    inv_path = os.path.join(args.outdir, 'uvot_inventory.txt')
    phot_path = os.path.join(args.outdir, 'uvot_photometry.txt')
    pdf_path = os.path.join(args.outdir, 'uvot_photometry.pdf')
    for path, step in ((inv_path, 'Step 3 (swift_uvot_inventory.py)'),
                       (pos_path, 'Step 4 (swift_uvot_positions.py)')):
        if not os.path.exists(path):
            print('ERROR: %s not found; run %s first.' % (path, step),
                  file=sys.stderr)
            return 1
    positions = read_table(pos_path)
    pos = {(r['obsid'], r['filter'], r['extname']): r for r in positions}
    inv = {(r['obsid'], r['filter'], r['extname']): r
           for r in read_table(inv_path)}
    with open(pos_path) as fh:
        header = ''.join(line for line in fh if line.startswith('#'))
    m = re.search(r'source: RA ([\d.+-]+) Dec ([\d.+-]+)', header)
    src_label = ('RA %s Dec %s' % (m.group(1), m.group(2))) if m else '?'

    if args.plot_only:
        if not os.path.exists(phot_path):
            print('ERROR: %s not found; run Step 5 first.' % phot_path,
                  file=sys.stderr)
            return 1
        plot_photometry(pdf_path, read_table(phot_path), pos, inv, src_label)
        print('[summary] Written: %s' % pdf_path)
        return 0

    uenv.require_heasoft_shell(['uvotsource'])
    todo = [r for r in positions if r['status'] == 'ok']
    print('[info] source %s; %d exposures with regions from Step 4 (of %d '
          'candidates)' % (src_label, len(todo), len(positions)))
    results = run_photometry(todo, inv, args.outdir, args.nproc)

    by_key = {}
    det_diffs = []
    for r, res in zip(todo, results):
        key = (r['obsid'], r['filter'], r['extname'])
        out = {'obsid': r['obsid'], 'filter': r['filter'],
               'extname': r['extname']}
        if not res.ok:
            why = _errors(res.stdout + '\n' + res.stderr)[:2] \
                or res.tail.splitlines()[-2:]
            out['status'] = 'failed'
            out['reason'] = 'uvotsource: %s; %s (log: %s)' % (
                res.reason, ' / '.join(why), res.log)
            by_key[key] = out
            continue
        try:
            row = harvest(res.outputs['phot.fits'], r['extname'])
            values, notes = photometry_row(row, r)
        except Exception as exc:  # unreadable or unexpected output table
            out['status'] = 'failed'
            out['reason'] = 'cannot read the uvotsource output: %s' % exc
            by_key[key] = out
            continue
        out.update(values)
        out['phot_file'] = os.path.relpath(res.outputs['phot.fits'],
                                           args.outdir)
        notes += ['uvotsource: ' + w for w in _warnings(res.stdout
                                                         + res.stderr)[:2]]
        out['status'] = 'ok'
        out['reason'] = '; '.join(notes)
        try:
            det_diffs.append(max(abs(float(row['DETX']) - float(r['detx'])),
                                 abs(float(row['DETY']) - float(r['dety']))))
        except (TypeError, ValueError):
            pass
        by_key[key] = out
    # exposures Step 4 did not pass keep their Step 4 status
    for r in positions:
        key = (r['obsid'], r['filter'], r['extname'])
        if key not in by_key:
            by_key[key] = {'obsid': r['obsid'], 'filter': r['filter'],
                           'extname': r['extname'], 'status': r['status'],
                           'reason': 'Step 4: %s' % r.get('reason', '')}
    # the ledger: every Step 4 row comes back exactly once, in its order
    rows = [by_key[(r['obsid'], r['filter'], r['extname'])] for r in positions]
    if len(by_key) != len(positions):
        print('ERROR: internal: %d exposures in, %d rows out'
              % (len(positions), len(by_key)), file=sys.stderr)
        return 1

    counts = Counter(r['status'] for r in rows)
    prov = uenv.provenance()
    comments = [
        'swift_uvot_photometry.py, %s' % datetime.now(timezone.utc)
        .strftime('%Y-%m-%dT%H:%M:%SZ'),
        'command: %s' % ' '.join(sys.argv),
        'source: %s' % src_label,
        'pipeline %s; HEASoft %s; CALDB %s (UVOT %s)' % (
            prov['PIPECOMM'], prov['HEASOFT'], prov['CALDB'], prov['UVOTINDX']),
        '%d exposures: %s' % (len(rows), ', '.join(
            '%s %d' % (k, v) for k, v in sorted(counts.items()))),
        'uvotsource %s' % ' '.join('%s=%s' % kv for kv in
                                   sorted(UVOTSOURCE_PARAMS.items())),
        'counts_frame: raw counts per frame in the 5" circle (uvotsource '
        'clamps above %.2f: saturated=1, rate a lower limit)'
        % FRAME_COUNT_LIMIT,
        'src_cnts, bkg_cnts: raw counts in src_area, bkg_area (arcsec^2); '
        'coi, coi_bkg, lss, sss, senscorr: correction factors',
        'rate, rate_err: corrected rate (counts/s; uvotsource SENSCORR_RATE, '
        'also on a low-sensitivity patch); rate_lim, mag_lim: 3-sigma '
        'background limits',
        'mag, mag_err, mag_sys: Vega magnitude, statistical and systematic '
        'error (99 on a LOW patch); flux_aa: erg/s/cm^2/A',
    ]
    write_table(phot_path, rows, COLUMNS, comments)

    # --- report -------------------------------------------------------------
    ok = [r for r in rows if r['status'] == 'ok']
    print()
    print('[summary] %d exposures: %s' % (len(rows), ', '.join(
        '%s %d' % (k, v) for k, v in sorted(counts.items()))))
    print('          filter  measured  saturated  >0.95 counts/frame  '
          'LOW patch  nsigma<3')
    for f in FILTERS:
        sel = [r for r in ok if r['filter'] == f]
        if sel:
            print('          %-6s  %8d  %9d  %18d  %9d  %8d' % (
                f, len(sel), sum(int(r['saturated']) for r in sel),
                sum(float(r['counts_frame']) > 0.95 for r in sel),
                sum('BAD_SSS' in r['photflags'] for r in sel),
                sum(float(r['nsigma']) < 3 for r in sel)))
    flags = Counter(f for r in ok for f in r['photflags'].split(',')
                    if f != '-')
    print('[summary] photometry flags: %s' % (', '.join(
        '%s %d' % kv for kv in sorted(flags.items())) or 'none'))
    differ = [r for r in ok if 'differs from Step 4' in r['reason']]
    n_cmp = sum(1 for r in ok if pos[(r['obsid'], r['filter'], r['extname'])]
                .get('sss_low') not in (None, '-'))
    print('[check] low-sensitivity flag: uvotsource agrees with Step 4 for '
          '%d of %d exposures' % (n_cmp - len(differ), n_cmp))
    if det_diffs:
        print('[check] detector position: uvotsource and Step 4 agree to '
              '%.2f pixel (largest difference)' % max(det_diffs))
    warned = [r for r in ok if 'uvotsource: ' in r['reason']]
    if warned:
        print('[note] uvotsource printed warnings for %d exposure(s), e.g. '
              '%s %s: %s' % (len(warned), warned[0]['obsid'],
                             warned[0]['extname'],
                             warned[0]['reason'].split('uvotsource: ')[1]))
    failed = [r for r in rows if r['status'] == 'failed']
    for r in failed[:10]:
        print('[fail] %s %s %s: %s' % (r['obsid'], r['filter'], r['extname'],
                                       r['reason']))

    plot_photometry(pdf_path, rows, pos, inv, src_label)
    print('[summary] Written: %s, %s, %s/photometry/' % (phot_path, pdf_path,
                                                        args.outdir))
    if failed:
        print('ERROR: %d exposure(s) failed; see the reason column.'
              % len(failed), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
