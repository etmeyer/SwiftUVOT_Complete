#!/usr/bin/env python3
"""
make_uvot_master_table.py -- Step 6: the quality rules and the master table.

Joins Steps 3-5 into one row per exposure -- every extension of every sky
image, so the ledger closes -- and decides by rule which exposures the light
curve uses: include = yes or no, with the rules that excluded it. The table
is regenerated, never edited. Your own decisions go in uvot_overrides.txt
in the output folder, which is applied every time:

    # obsid        filter  extname        action   note
    00031659123    *       *              exclude  streaked second snapshot
    00035017124    U       uu387398273I   include  checked by eye

(* matches anything; a later line wins over an earlier one.)

Exclusion rules (codes in the rules column). The limits were calibrated on
3C 273 (rates against the same observation's other exposures; see the
documentation):

    step3:<status>  not a candidate in Step 3 (not_covered, partial, ...)
    step4:<status>  no regions from Step 4 (no_background, failed)
    step5:failed    uvotsource failed
    saturated       above 0.98 counts per frame: the rate is a lower limit
    coincidence     above the filter's counts-per-frame limit (0.95; U
                    only at the saturation cap, but flagged high_coi above
                    0.90). UVW2 above 0.95 reads 13-39 % high; full-frame
                    U above 0.93 reads 8 % low against the windowed
                    exposure of the same snapshot.
    sss             on a low-sensitivity patch at --sss-level (default
                    LOW: 5 % low on average; MID-only patches 1 %)
    off_centre      the source was found more than 3" from the region
                    centre, so Step 4 left the region at --ra/--dec
    smeared         concentration (5"/15" counts) below 0.78
    elongated       axis ratio above 1.5
    jitter          taken during the 2023-24 spacecraft jitter with no
                    measurable PSF (otherwise the shape rules decide)
    duplicate       same EXPID as an included exposure
    override:N      excluded by line N of uvot_overrides.txt

Flags (noted, not excluding): sss_mid / sss_high (patch at a level the rule
does not exclude), jitter_ok (jitter period, PSF check passed), high_coi
(U above 0.90 counts per frame: about 8 % low), no_aspcorr, short
(exposure under 20 s), discrepant (differs from the rest of its
observation by more than 10 % and 5 sigma), override:N (included by line N,
despite the rules listed).

Output (in --outdir):
    uvot_master_table.txt   one row per exposure
    uvot_master_table.pdf   the rules against the data, then every
                            exposure's rate against time per filter

Rows whose include changed since the last run are listed.

Usage:
    make_uvot_master_table.py
    make_uvot_master_table.py --sss-level MID
    make_uvot_master_table.py --coincidence-limit U=0.95
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone

import numpy as np

import swift_uvot_env as uenv
from swift_uvot_tables import read_table, write_table

SATURATION = 0.98               # uvotsource caps counts per frame here
# Counts per frame above which an exposure is excluded; other filters:
# DEFAULT_COINCIDENCE_LIMIT. U is kept up to the saturation cap, because
# for bright sources full-frame U is often all there is, but flagged above
# COI_FLAG_LIMIT: it reads about 8 % low there (3C 273).
COINCIDENCE_LIMIT = {'U': SATURATION}
DEFAULT_COINCIDENCE_LIMIT = 0.95
COI_FLAG_LIMIT = {'U': 0.90}
CONC_MIN = 0.78                 # 5"/15" counts: 3 % low at 0.75-0.78
AXIS_RATIO_MAX = 1.5            # 5-12 % low above 1.5
SHORT_EXPOSURE = 20.0           # s
# an included exposure this far from the rest of its observation (and
# filter) is flagged 'discrepant'
DISCREPANT_FRACTION, DISCREPANT_SIGMA = 0.10, 5.0
# Swift's gyro noise: UVOT images degraded from about 2023-08-07 (GCN 34633)
# until UVOT was re-enabled in two-gyro mode on 2024-04-04 (GCN 36033)
JITTER = ('2023-08-07', '2024-04-04')
SSS_LEVELS = ('LOW', 'MID', 'HIGH')
OVERRIDES_FILE = 'uvot_overrides.txt'
FILTERS = ('V', 'B', 'U', 'UVW1', 'UVM2', 'UVW2')
# for the plot, an excluded exposure is shown under its first rule here
RULE_ORDER = ('saturated', 'coincidence', 'jitter', 'off_centre', 'smeared',
              'elongated', 'sss', 'duplicate', 'override', 'step5', 'step4')

COLUMNS = ['obsid', 'filter', 'extname', 'date_mid', 'mjd_mid', 'exposure',
           'frametime', 'mode', 'aspcorr', 'first', 'snr', 'offset',
           'axis_ratio', 'conc', 'center', 'sss_low', 'sss_mid', 'sss_high',
           'counts_frame', 'saturated', 'rate', 'rate_err', 'rate_lim', 'mag',
           'mag_err', 'mag_sys', 'flux_aa', 'flux_aa_err', 'include', 'rules',
           'flags']


def _num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float('nan')


def read_overrides(path):
    """[(line number, obsid, filter, extname, action, note)] from the file."""
    rules = []
    if not os.path.exists(path):
        return rules
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            text = line.split('#', 1)[0].strip()
            if not text:
                continue
            parts = text.split(None, 4)
            if len(parts) < 4 or parts[3] not in ('include', 'exclude'):
                raise ValueError('%s:%d: expected "OBSID FILTER EXTNAME '
                                 'include|exclude [note]"' % (path, lineno))
            rules.append((lineno, parts[0], parts[1], parts[2], parts[3],
                          parts[4] if len(parts) > 4 else ''))
    return rules


def matching_override(overrides, row):
    """The last override line matching this row, or None."""
    found = None
    for rule in overrides:
        _, obsid, filt, ext, _, _ = rule
        if all(pat == '*' or pat == val for pat, val in
               ((obsid, row['obsid']), (filt, row['filter']),
                (ext, row['extname']))):
            found = rule
    return found


def judge(inv, pos, phot, sss_level, limits):
    """
    (rules, flags) for one exposure; rules empty means include. limits:
    counts-per-frame limit per filter.
    """
    rules, flags = [], []
    if inv['status'] != 'candidate':
        return ['step3:%s' % inv['status']], flags
    if pos is None or pos['status'] != 'ok':
        return ['step4:%s' % (pos['status'] if pos else 'missing')], flags
    if phot is None or phot['status'] != 'ok':
        return ['step5:%s' % (phot['status'] if phot else 'missing')], flags

    cpf = _num(phot['counts_frame'])
    if phot['saturated'] == '1' or cpf > SATURATION:
        rules.append('saturated')
    elif cpf > limits.get(inv['filter'], DEFAULT_COINCIDENCE_LIMIT):
        rules.append('coincidence')
    elif cpf > COI_FLAG_LIMIT.get(inv['filter'], 2.0):
        flags.append('high_coi')
    worst = SSS_LEVELS.index(sss_level)
    on_patch = [lvl for lvl in SSS_LEVELS
                if _num(pos['sss_%s' % lvl.lower()]) <= 0]
    if any(SSS_LEVELS.index(lvl) <= worst for lvl in on_patch):
        rules.append('sss')
    elif on_patch:
        flags.append('sss_%s' % on_patch[0].lower())
    snr = _num(pos['snr'])
    if pos['center'] == 'catalog' and snr >= 10:
        rules.append('off_centre')
    conc, ar = _num(pos['conc']), _num(pos['axis_ratio'])
    if conc < CONC_MIN:
        rules.append('smeared')
    if ar > AXIS_RATIO_MAX:
        rules.append('elongated')
    if JITTER[0] <= inv['date_mid'][:10] < JITTER[1]:
        if conc != conc:    # no clear source within 15": no PSF to check
            rules.append('jitter')
        elif not rules:
            flags.append('jitter_ok')
    if inv['aspcorr'] == 'NONE':
        flags.append('no_aspcorr')
    if _num(inv['exposure']) < SHORT_EXPOSURE:
        flags.append('short')
    return rules, flags


def resolve_duplicates(rows):
    """Among included rows sharing an EXPID keep one: image mode, longest."""
    by_expid = defaultdict(list)
    for r in rows:
        if r['include'] == 'yes' and r.get('_expid') not in (None, '-', ''):
            by_expid[(r['obsid'], r['_expid'])].append(r)
    for group in by_expid.values():
        if len(group) < 2:
            continue
        group.sort(key=lambda r: (r['mode'] != 'IMAGE', -_num(r['exposure'])))
        for r in group[1:]:
            r['include'] = 'no'
            r['_rules'].append('duplicate')


def reference_ratios(rows):
    """
    For the plot and summary: each measured exposure's rate divided by the
    median rate of the other included exposures in its observation and
    filter (NaN without one).
    """
    groups = defaultdict(list)
    for r in rows:
        if r.get('rate') not in (None, '-'):
            groups[(r['obsid'], r['filter'])].append(r)
    for group in groups.values():
        for r in group:
            ref = [_num(o['rate']) for o in group
                   if o is not r and o['include'] == 'yes']
            r['_ratio'] = (_num(r['rate']) / np.median(ref)
                           if ref and np.median(ref) > 0 else float('nan'))


def rule_code(rule):
    """'override:12' -> 'override'; other rule codes as they are."""
    return 'override' if rule.startswith('override:') else rule


def flag_discrepant(rows):
    """
    Flag included exposures that differ from the error-weighted mean of
    the other included exposures of their observation and filter by more
    than DISCREPANT_FRACTION and DISCREPANT_SIGMA (both errors combined).
    """
    groups = defaultdict(list)
    for r in rows:
        if r['include'] == 'yes':
            groups[(r['obsid'], r['filter'])].append(r)
    for group in groups.values():
        for r in group:
            others = [(_num(o['rate']), _num(o['rate_err'])) for o in group
                      if o is not r and _num(o['rate_err']) > 0]
            if not others:
                continue
            w = np.array([1.0 / e ** 2 for _, e in others])
            mean = np.sum(w * [v for v, _ in others]) / np.sum(w)
            err = np.hypot(_num(r['rate_err']), 1.0 / np.sqrt(np.sum(w)))
            diff = _num(r['rate']) - mean
            if mean > 0 and abs(diff) > DISCREPANT_FRACTION * mean \
                    and abs(diff) > DISCREPANT_SIGMA * err:
                r['_flags'].append('discrepant')


def primary_rule(r):
    for code in RULE_ORDER:
        if any(rule.split(':')[0] == code for rule in r['_rules']):
            return code
    return r['_rules'][0].split(':')[0] if r['_rules'] else 'included'


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_master(pdf_path, rows, src_label, limits):
    """The rules against the data, then rate against time per filter."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from astropy.time import Time

    measured = [r for r in rows if r.get('rate') not in (None, '-')]
    for r in measured:
        r['_year'] = Time(_num(r['mjd_mid']), format='mjd').decimalyear
    colors = {'included': 'C0', 'saturated': 'C3', 'coincidence': 'C1',
              'jitter': 'k', 'off_centre': 'C5', 'smeared': 'C4',
              'elongated': 'C6', 'sss': 'C8', 'duplicate': 'C7',
              'override': 'C2', 'step5': 'C9', 'step4': 'C9'}

    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(1, 3, figsize=(11, 5))

        def panel(ax, key, own_rule, xlim, limit, xlabel):
            sel = [r for r in measured if np.isfinite(r['_ratio'])
                   and set(x.split(':')[0] for x in r['_rules'])
                   <= {own_rule}]
            for inc, color, label in ((True, 'C0', 'included'),
                                      (False, colors[own_rule],
                                       'excluded: ' + own_rule)):
                pts = [r for r in sel if (r['include'] == 'yes') == inc]
                ax.scatter([_num(r[key]) for r in pts],
                           [r['_ratio'] for r in pts], s=5, color=color,
                           label='%s (%d)' % (label, len(pts)))
            if limit is not None:
                ax.axvline(limit, color='k', ls=':', lw=0.8)
            ax.axhline(1.0, color='grey', lw=0.6)
            ax.set_xlim(*xlim)
            ax.set_ylim(0, 1.5)
            ax.set_xlabel(xlabel)
            ax.legend(fontsize=7, loc='lower right')

        panel(axes[0], 'conc', 'smeared', (0.3, 1.0), CONC_MIN,
              'concentration (5"/15" counts)')
        axes[0].set_ylabel('rate / median of the observation\'s other '
                           'included exposures')
        panel(axes[1], 'axis_ratio', 'elongated', (1.0, 3.0), AXIS_RATIO_MAX,
              'axis ratio')
        ax = axes[2]
        sel = [r for r in measured if np.isfinite(r['_ratio'])
               and set(x.split(':')[0] for x in r['_rules'])
               <= {'coincidence', 'saturated'}]
        for f in FILTERS:
            pts = [r for r in sel if r['filter'] == f]
            if pts:
                ax.scatter([_num(r['counts_frame']) for r in pts],
                           [r['_ratio'] for r in pts], s=5, label=f)
        for value in sorted(set(limits.values())):
            ax.axvline(value, color='k', ls=':', lw=0.8)
        ax.axvline(SATURATION, color='C3', ls=':', lw=0.8)
        ax.axhline(1.0, color='grey', lw=0.6)
        ax.set_xlim(0.2, 1.01)
        ax.set_ylim(0, 1.5)
        ax.set_xlabel('counts per frame (dotted: limits; red: saturation)')
        ax.legend(fontsize=7, loc='lower left', ncol=2)
        fig.suptitle('Step 6: the rules against the data (%s); each panel '
                     'shows exposures no other rule excludes' % src_label,
                     fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        pdf.savefig(fig)
        plt.close(fig)

        for f in FILTERS:
            sel = [r for r in measured if r['filter'] == f]
            if not sel:
                continue
            fig, ax = plt.subplots(figsize=(11, 5.5))
            groups = defaultdict(list)
            for r in sel:
                groups[primary_rule(r)].append(r)
            for name in ['included'] + [c for c in RULE_ORDER if c in groups]:
                pts = groups.get(name, [])
                if not pts:
                    continue
                t = [r['_year'] for r in pts]
                y = [_num(r['rate']) for r in pts]
                label = '%s (%d)' % ('included' if name == 'included'
                                     else 'excluded: ' + name, len(pts))
                if name == 'included':
                    ax.errorbar(t, y, yerr=[_num(r['rate_err']) for r in pts],
                                fmt='o', ms=3.5, lw=0.6, color=colors[name],
                                label=label, zorder=3)
                else:
                    ax.scatter(t, y, s=10, marker='^' if name == 'saturated'
                               else 'x', color=colors.get(name, 'C7'),
                               label=label, zorder=2)
            ys = np.array([_num(r['rate']) for r in sel])
            ys = ys[np.isfinite(ys)]
            if ys.size:
                lo, hi = np.percentile(ys, [0.5, 99.5])
                if hi > 0:
                    ax.set_ylim(min(0.0, lo * 1.1), hi * 1.15)
            ax.set_xlabel('year')
            ax.set_ylabel('corrected count rate (counts/s)')
            ax.set_title('%s: every exposure, included or excluded by its '
                         'first rule' % f)
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
        description='Step 6: quality rules and the master table.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--outdir', default='UVOT_output',
                        help='Steps 3-5 output folder (default: UVOT_output)')
    parser.add_argument('--sss-level', default='LOW', choices=SSS_LEVELS,
                        help='exclude exposures on a low-sensitivity patch at '
                             'this level or a stricter one (default: LOW)')
    parser.add_argument('--coincidence-limit', action='append', default=[],
                        metavar='FILTER=VALUE',
                        help='counts-per-frame limit for a filter (default: '
                             '0.95, U %.2f); can be repeated' % SATURATION)
    args = parser.parse_args(argv)
    limits = {f: COINCIDENCE_LIMIT.get(f, DEFAULT_COINCIDENCE_LIMIT)
              for f in FILTERS}
    for item in args.coincidence_limit:
        try:
            filt, value = item.split('=')
            if filt not in FILTERS or not 0 < float(value) <= SATURATION:
                raise ValueError
            limits[filt] = float(value)
        except ValueError:
            print('ERROR: --coincidence-limit %s: expected FILTER=VALUE with '
                  'FILTER one of %s and 0 < VALUE <= %.2f'
                  % (item, ' '.join(FILTERS), SATURATION), file=sys.stderr)
            return 1

    paths = {name: os.path.join(args.outdir, 'uvot_%s.txt' % name)
             for name in ('inventory', 'positions', 'photometry')}
    steps = {'inventory': 'Step 3 (swift_uvot_inventory.py)',
             'positions': 'Step 4 (swift_uvot_positions.py)',
             'photometry': 'Step 5 (swift_uvot_photometry.py)'}
    for name, path in paths.items():
        if not os.path.exists(path):
            print('ERROR: %s not found; run %s first.' % (path, steps[name]),
                  file=sys.stderr)
            return 1
    inventory = read_table(paths['inventory'])
    pos = {(r['obsid'], r['filter'], r['extname']): r
           for r in read_table(paths['positions'])}
    phot = {(r['obsid'], r['filter'], r['extname']): r
            for r in read_table(paths['photometry'])}
    with open(paths['positions']) as fh:
        header = ''.join(line for line in fh if line.startswith('#'))
    m = re.search(r'source: RA ([\d.+-]+) Dec ([\d.+-]+)', header)
    src_label = ('RA %s Dec %s' % (m.group(1), m.group(2))) if m else '?'
    try:
        overrides = read_overrides(os.path.join(args.outdir, OVERRIDES_FILE))
    except ValueError as exc:
        print('ERROR: %s' % exc, file=sys.stderr)
        return 1

    master_path = os.path.join(args.outdir, 'uvot_master_table.txt')
    previous = {}
    if os.path.exists(master_path):
        previous = {(r['obsid'], r['filter'], r['extname']): r['include']
                    for r in read_table(master_path)}

    rows = []
    for inv in inventory:
        key = (inv['obsid'], inv['filter'], inv['extname'])
        p, ph = pos.get(key), phot.get(key)
        row = {k: inv.get(k) for k in ('obsid', 'filter', 'extname',
                                       'date_mid', 'mjd_mid', 'exposure',
                                       'frametime', 'mode', 'aspcorr',
                                       'first')}
        row['_expid'] = inv.get('expid')
        if p:
            row.update({k: p.get(k) for k in ('snr', 'offset', 'axis_ratio',
                                              'conc', 'center', 'sss_low',
                                              'sss_mid', 'sss_high')})
        if ph and ph['status'] == 'ok':
            row.update({k: ph.get(k) for k in (
                'counts_frame', 'saturated', 'rate', 'rate_err', 'rate_lim',
                'mag', 'mag_err', 'mag_sys', 'flux_aa', 'flux_aa_err')})
            row['exposure'] = ph['exposure']
        rules, flags = judge(inv, p, ph, args.sss_level, limits)
        row['_rules'], row['_flags'] = rules, flags
        row['include'] = 'no' if rules else 'yes'
        rows.append(row)
    resolve_duplicates(rows)

    unusable = []
    for row in rows:
        ovr = matching_override(overrides, row)
        if ovr is None:
            continue
        lineno, action = ovr[0], ovr[4]
        if action == 'exclude':
            row['include'] = 'no'
            row['_rules'].append('override:%d' % lineno)
        elif row.get('rate') in (None, '-'):
            unusable.append((row, lineno))
        else:
            row['include'] = 'yes'
            row['_flags'].append('override:%d' % lineno)
    flag_discrepant(rows)
    for row in rows:
        row['rules'] = ','.join(row['_rules']) or '-'
        row['flags'] = ','.join(row['_flags']) or '-'
    reference_ratios(rows)

    # --- write -------------------------------------------------------------
    n_inc = sum(r['include'] == 'yes' for r in rows)
    prov = uenv.provenance()
    limit_text = ', '.join('%s %.2f' % (f, limits[f]) for f in FILTERS)
    comments = [
        'make_uvot_master_table.py, %s' % datetime.now(timezone.utc)
        .strftime('%Y-%m-%dT%H:%M:%SZ'),
        'command: %s' % ' '.join(sys.argv),
        'source: %s' % src_label,
        'pipeline %s' % prov['PIPECOMM'],
        '%d exposures: %d included, %d excluded' % (len(rows), n_inc,
                                                    len(rows) - n_inc),
        'rules: saturated > %.2f counts/frame; coincidence above %s; sss at '
        '%s; off_centre; smeared conc < %.2f; elongated axis ratio > %.1f; '
        'jitter %s to %s without a measurable PSF; duplicate EXPID'
        % (SATURATION, limit_text, args.sss_level, CONC_MIN, AXIS_RATIO_MAX,
           JITTER[0], JITTER[1]),
        'overrides: %s (%d line(s))' % (OVERRIDES_FILE, len(overrides)),
        'include: yes/no; rules: the rules that exclude the exposure (if an '
        'override includes it anyway, flags say so)',
    ]
    write_table(master_path, rows, COLUMNS, comments)

    # --- report -------------------------------------------------------------
    print('[info] source %s; %d exposures in the inventory' % (src_label,
                                                              len(rows)))
    if overrides:
        print('[info] %d override line(s) from %s' % (len(overrides),
                                                     OVERRIDES_FILE))
    for row, lineno in unusable:
        print('[warn] %s line %d includes %s %s %s, which has no photometry '
              '(%s): ignored' % (OVERRIDES_FILE, lineno, row['obsid'],
                                 row['filter'], row['extname'], row['rules']))
    print()
    print('[summary] %d exposures: %d included, %d excluded'
          % (len(rows), n_inc, len(rows) - n_inc))
    print('          filter  exposures  included  included exposure (s)')
    for f in FILTERS:
        sel = [r for r in rows if r['filter'] == f]
        if sel:
            inc = [r for r in sel if r['include'] == 'yes']
            print('          %-6s  %9d  %8d  %21.0f' % (
                f, len(sel), len(inc), sum(_num(r['exposure']) for r in inc)))
    counts = Counter(rule_code(rule) for r in rows for rule in r['_rules'])
    print('[summary] exposures per rule (one exposure can break several):')
    for rule, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        ratios = [r['_ratio'] for r in rows
                  if rule in [rule_code(x) for x in r['_rules']]
                  and np.isfinite(r.get('_ratio', np.nan))]
        effect = (', rate %.2f of the included median (%d with one)'
                  % (np.median(ratios), len(ratios))) if ratios else ''
        print('          %-24s %5d%s' % (rule, n, effect))
    flag_counts = Counter(flag.split(':')[0] for r in rows
                          for flag in r['_flags'] if r['include'] == 'yes')
    print('[summary] flags on included exposures: %s' % (', '.join(
        '%s %d' % kv for kv in sorted(flag_counts.items())) or 'none'))
    odd = [r for r in rows if 'discrepant' in r['_flags']]
    print('[check] included exposures that differ from the rest of their '
          'observation by more than %.0f %% and %.0f sigma: %d'
          % (100 * DISCREPANT_FRACTION, DISCREPANT_SIGMA, len(odd)))
    for r in odd[:10]:
        print('          %s %-5s %s rate %s +- %s' % (
            r['obsid'], r['filter'], r['extname'], r['rate'], r['rate_err']))
    if previous:
        changed = [r for r in rows if previous.get(
            (r['obsid'], r['filter'], r['extname'])) not in (None,
                                                            r['include'])]
        new = sum(1 for r in rows if (r['obsid'], r['filter'], r['extname'])
                  not in previous)
        print('[summary] since the last run: %d exposure(s) changed include, '
              '%d new' % (len(changed), new))
        for r in changed[:20]:
            print('          %s %-5s %s: %s -> %s (%s)' % (
                r['obsid'], r['filter'], r['extname'],
                previous[(r['obsid'], r['filter'], r['extname'])],
                r['include'], r['rules'] if r['include'] == 'no'
                else r['flags']))
    pdf_path = os.path.join(args.outdir, 'uvot_master_table.pdf')
    plot_master(pdf_path, rows, src_label, limits)
    print('[summary] Written: %s, %s' % (master_path, pdf_path))
    return 0


if __name__ == '__main__':
    sys.exit(main())
