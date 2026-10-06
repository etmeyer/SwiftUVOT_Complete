#!/usr/bin/env python3
"""
swift_uvot_viewer.py -- Step 4 (images): look at the exposures that need
looking at, with their regions.

Reads Step 4's uvot_positions.txt (and Step 3's inventory) and makes one
PDF:

  * a summary page: axis ratio and concentration (5"/15" counts) of the
    source for snapshot openers in event mode and for all other exposures,
    and both against time; centroid offsets with and without aspect
    correction; low-sensitivity flags per filter;
  * contact sheets, worst first: the sources with the most light outside
    the 5" region, the faintest (S/N < 10), the most elongated, the
    largest offsets, exposures whose background is not the standard
    annulus, ones on a LOW low-sensitivity patch, and a random sample of
    the rest for comparison. Each cutout (90" across) shows the 5" source
    region (green), your --ra/--dec (red cross), the background region
    (cyan) and its source exclusions (dashed);
  * with --obsid, every candidate exposure of those observations instead.

Usage:
    swift_uvot_viewer.py                     # UVOT_output/uvot_positions.pdf
    swift_uvot_viewer.py --obsid 00035017124 00035017045 --pdf two_obs.pdf
"""

import argparse
import gzip
import io
import os
import random
import re
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning

from swift_uvot_tables import read_table
from swift_uvot_positions import CONC_MIN, MIN_SNR

warnings.simplefilter('ignore', FITSFixedWarning)

CUTOUT = 45.0       # arcsec half-width of a cutout
JITTER = (2023.6, 2024.26)  # spacecraft jitter, Aug 2023 - Apr 2024 (GCN)
NCOL, NROW = 6, 4   # cutouts per page
FILTERS = ('V', 'B', 'U', 'UVW1', 'UVM2', 'UVW2')
_REGION = re.compile(r'(-?)(circle|annulus)\(([^)]*)\)')


def _open_fits(path):
    if path.endswith('.gz'):
        with gzip.open(path, 'rb') as fh:
            return fits.open(io.BytesIO(fh.read()))
    return fits.open(path)


def read_regions(path):
    """[(exclude?, shape, ra, dec, radii arcsec...)] from our region files."""
    shapes = []
    if not path or path == '-' or not os.path.exists(path):
        return shapes
    for line in open(path):
        m = _REGION.search(line)
        if m:
            vals = [v.strip().rstrip('"') for v in m.group(3).split(',')]
            shapes.append((m.group(1) == '-', m.group(2),
                           float(vals[0]), float(vals[1]),
                           [float(v) for v in vals[2:]]))
    return shapes


def draw_cutout(ax, data, wcs, scale, row, outdir, src_ra, src_dec):
    """One cutout centred on the source region, with the regions drawn."""
    cra, cdec = float(row['src_ra']), float(row['src_dec'])
    x0, y0 = [float(v) for v in wcs.world_to_pixel_values(cra, cdec)]
    half = CUTOUT / scale
    xa, xb = int(x0 - half), int(x0 + half) + 1
    ya, yb = int(y0 - half), int(y0 + half) + 1
    ny, nx = data.shape
    cut = np.zeros((yb - ya, xb - xa))
    sx0, sx1, sy0, sy1 = max(0, xa), min(nx, xb), max(0, ya), min(ny, yb)
    if sx0 < sx1 and sy0 < sy1:
        cut[sy0 - ya:sy1 - ya, sx0 - xa:sx1 - xa] = np.nan_to_num(
            data[sy0:sy1, sx0:sx1])
    lo, hi = np.percentile(cut, 5), max(np.percentile(cut, 99.8), 1e-6)
    shown = np.arcsinh((cut - lo) / max((hi - lo) / 50.0, 1e-9))
    ax.imshow(shown, origin='lower', cmap='gray_r',
              extent=(xa - 0.5, xb - 0.5, ya - 0.5, yb - 0.5))

    def pix(ra, dec):
        return [float(v) for v in wcs.world_to_pixel_values(ra, dec)]

    for excl, shape, ra, dec, radii in read_regions(
            os.path.join(outdir, row['bkg_reg']) if row['bkg_reg'] != '-'
            else None):
        px, py = pix(ra, dec)
        for rad in radii:
            ax.add_patch(Circle((px, py), rad / scale, fill=False, lw=0.8,
                                color='c', ls='--' if excl else '-'))
    sx, sy = pix(cra, cdec)
    ax.add_patch(Circle((sx, sy), 5.0 / scale, fill=False, lw=1.0, color='lime'))
    tx, ty = pix(src_ra, src_dec)
    ax.plot([tx], [ty], marker='+', color='r', ms=8, mew=1.0)
    ax.set_xlim(xa - 0.5, xb - 0.5)
    ax.set_ylim(ya - 0.5, yb - 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    flags = []
    if row['mode'] == 'EVENT':
        flags.append('ev')
    if row['first'] == 'yes':
        flags.append('1st')
    sss = [lvl[0] for lvl in ('low', 'mid', 'high')
           if row.get('sss_' + lvl) not in (None, '-')
           and float(row['sss_' + lvl]) <= 0]
    if sss:
        flags.append('SSS:' + ''.join(sss).upper())
    ax.set_title('%s %s %s\nS/N %s ar %s c %s off %s"\n%s %s' % (
        row['obsid'], row['filter'], row['extname'], row['snr'],
        row['axis_ratio'], row.get('conc', '-'), row['offset'],
        row['aspcorr'][:1], ' '.join(flags)), fontsize=5.5)


def contact_pages(pdf, title, rows, files, outdir, src_ra, src_dec):
    """Pages of cutouts for rows (already ordered)."""
    for start in range(0, len(rows), NCOL * NROW):
        chunk = rows[start:start + NCOL * NROW]
        fig, axes = plt.subplots(NROW, NCOL, figsize=(11, 8.5),
                                 gridspec_kw={'hspace': 0.6, 'wspace': 0.1})
        for ax in axes.flat:
            ax.axis('off')
        for ax, row in zip(axes.flat, chunk):
            ax.axis('on')
            try:
                hdul = files(row)
                hdu = hdul[int(row['hdu'])]
                wcs = WCS(hdu.header)
                draw_cutout(ax, hdu.data, wcs, abs(hdu.header['CDELT1']) * 3600,
                            row, outdir, src_ra, src_dec)
            except Exception as exc:  # show the problem instead of failing
                ax.text(0.5, 0.5, 'cannot show:\n%s' % exc, ha='center',
                        va='center', fontsize=6, transform=ax.transAxes)
        fig.suptitle('%s (%d-%d of %d)' % (title, start + 1,
                                            start + len(chunk), len(rows)),
                     fontsize=10)
        fig.text(0.5, 0.01, 'green: 5" source region; red +: --ra/--dec; '
                 'cyan: background (dashed: excluded sources); ar: axis ratio; '
                 'c: 5"/15" counts; off: centroid offset; D/N: ASPCORR; '
                 'ev: event mode; '
                 '1st: opens its snapshot; SSS: low-sensitivity patch at '
                 'L/M/H', ha='center', fontsize=6)
        fig.subplots_adjust(left=0.02, right=0.98, top=0.9, bottom=0.05)
        pdf.savefig(fig)
        plt.close(fig)


def summary_page(pdf, rows, inv):
    """Distributions of the Step 4 measurements."""
    def val(r, k):
        return float(r[k]) if r.get(k) not in (None, '-') else np.nan

    opener = [r['mode'] == 'EVENT' and r['first'] == 'yes' for r in rows]
    t = [Time(float(inv[(r['obsid'], r['filter'], r['extname'])]['mjd_mid']),
              format='mjd').decimalyear for r in rows]
    groups = (('all other', [not o for o in opener]),
              ('event mode, opens snapshot', opener))
    fig, axes = plt.subplots(2, 3, figsize=(11, 8.5))

    def histogram(ax, key, bins, lo, hi, label):
        for name, sel in groups:
            v = [val(r, key) for r, s in zip(rows, sel) if s]
            ax.hist(np.clip(v, lo, hi), bins, alpha=0.6,
                    label='%s (%d)' % (name, len(v)))
        ax.set_yscale('log')
        ax.set_xlabel(label)
        ax.legend(fontsize=7)

    def against_time(ax, key, lo, hi, label):
        for name, sel in groups:
            ax.scatter([x for x, s in zip(t, sel) if s],
                       [val(r, key) for r, s in zip(rows, sel) if s], s=3,
                       label=name)
        ax.axvspan(*JITTER, color='orange', alpha=0.2,
                   label='spacecraft jitter (Aug 2023 - Apr 2024)')
        ax.set_ylim(lo, hi)
        ax.set_xlabel('year')
        ax.set_ylabel(label)
        ax.legend(fontsize=6)

    histogram(axes[0, 0], 'axis_ratio', np.linspace(1.0, 3.0, 41), 1, 3,
              'axis ratio of the source')
    histogram(axes[0, 1], 'conc', np.linspace(0.0, 1.2, 49), 0, 1.2,
              'net counts in 5" / in 15" (dotted: flag limit)')
    axes[0, 1].axvline(CONC_MIN, color='k', ls=':', lw=0.8)
    ax = axes[0, 2]
    bins = np.linspace(0, 6, 49)
    for asp in ('DIRECT', 'NONE'):
        sel = [val(r, 'offset') for r in rows if r['aspcorr'] == asp]
        ax.hist(np.clip(sel, 0, 6), bins, alpha=0.6,
                label='ASPCORR=%s (%d)' % (asp, len(sel)))
    ax.axvline(3.0, color='k', ls=':', lw=0.8)
    ax.set_yscale('log')
    ax.set_xlabel('centroid offset from --ra/--dec (arcsec)')
    ax.legend(fontsize=7)
    against_time(axes[1, 0], 'axis_ratio', 0.95, 3.0, 'axis ratio')
    against_time(axes[1, 1], 'conc', 0.0, 1.2, '5"/15" counts')
    ax = axes[1, 2]
    width = 0.27
    x = np.arange(len(FILTERS))
    for i, lvl in enumerate(('low', 'mid', 'high')):
        n = [sum(1 for r in rows if r['filter'] == f
                 and r.get('sss_' + lvl) not in (None, '-')
                 and float(r['sss_' + lvl]) <= 0) for f in FILTERS]
        tot = [sum(1 for r in rows if r['filter'] == f) for f in FILTERS]
        ax.bar(x + (i - 1) * width, [100.0 * a / max(b, 1) for a, b in
                                     zip(n, tot)], width, label=lvl.upper())
    ax.set_xticks(x)
    ax.set_xticklabels(FILTERS)
    ax.set_ylabel('% of exposures on a low-sensitivity patch')
    ax.legend(fontsize=7)
    fig.suptitle('Step 4: positions, shapes and low-sensitivity patches '
                 '(%d candidate exposures)' % len(rows))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Step 4 images: contact sheets of the exposures to '
                    'check, with their regions.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--outdir', default='UVOT_output',
                        help='Step 4 output folder (default: UVOT_output)')
    parser.add_argument('--pdf', default=None,
                        help='output PDF (default: <outdir>/uvot_positions.pdf)')
    parser.add_argument('--per-category', type=int, default=48,
                        help='cutouts per category (default: 48)')
    parser.add_argument('--obsid', nargs='+',
                        help='show every candidate exposure of these '
                             'observations instead of the categories')
    args = parser.parse_args(argv)

    pos_path = os.path.join(args.outdir, 'uvot_positions.txt')
    if not os.path.exists(pos_path):
        print('ERROR: %s not found; run Step 4 (swift_uvot_positions.py) '
              'first.' % pos_path, file=sys.stderr)
        return 1
    rows = [r for r in read_table(pos_path)
            if r['status'] in ('ok', 'no_background')]
    inv = {(r['obsid'], r['filter'], r['extname']): r for r in
           read_table(os.path.join(args.outdir, 'uvot_inventory.txt'))}
    header = ''.join(l for l in open(pos_path) if l.startswith('#'))
    m = re.search(r'source: RA ([\d.+-]+) Dec ([\d.+-]+)', header)
    src_ra, src_dec = float(m.group(1)), float(m.group(2))
    pdf_path = args.pdf or os.path.join(args.outdir, 'uvot_positions.pdf')

    cache = {}

    def files(row):
        path = inv[(row['obsid'], row['filter'], row['extname'])]['file']
        if path not in cache:
            if len(cache) > 32:
                cache.clear()
            cache[path] = _open_fits(path)
        return cache[path]

    def num(r, k):
        return float(r[k]) if r.get(k) not in (None, '-') else -1.0

    with PdfPages(pdf_path) as pdf:
        if args.obsid:
            sel = [r for r in rows if r['obsid'] in set(args.obsid)]
            sel.sort(key=lambda r: float(inv[(r['obsid'], r['filter'],
                                               r['extname'])]['tstart']))
            contact_pages(pdf, 'Observations %s' % ' '.join(args.obsid),
                          sel, files, args.outdir, src_ra, src_dec)
            n_cut = len(sel)
        else:
            summary_page(pdf, rows, inv)
            k = args.per_category
            lost = [r for r in rows if r.get('conc') not in (None, '-')
                    and float(r['conc']) < CONC_MIN]
            faint = [r for r in rows if num(r, 'snr') < MIN_SNR]
            cats = [
                ('Light outside the 5" region (5"/15" counts < %.2f)'
                 % CONC_MIN, sorted(lost, key=lambda r: float(r['conc']))[:k]),
                ('Faint at the source position (S/N < %.0f)' % MIN_SNR,
                 sorted(faint, key=lambda r: num(r, 'snr'))[:k]),
                ('Most elongated sources', sorted(
                    rows, key=lambda r: -num(r, 'axis_ratio'))[:k]),
                ('Largest centroid offsets', sorted(
                    rows, key=lambda r: -num(r, 'offset'))[:k]),
                ('Background not the standard annulus', [
                    r for r in rows if r['bkg_type'] != 'annulus'][:k]),
                ('On a LOW low-sensitivity patch', [
                    r for r in rows if num(r, 'sss_low') <= 0
                    and r.get('sss_low') not in (None, '-')][:k]),
            ]
            rest = [r for r in rows if num(r, 'axis_ratio') <= 1.3
                    and num(r, 'conc') >= CONC_MIN
                    and num(r, 'offset') <= 1.0 and r['bkg_type'] == 'annulus'
                    and num(r, 'sss_low') > 0]
            random.Random(1).shuffle(rest)
            cats.append(('Typical exposures, for comparison (random sample)',
                         rest[:min(k, 24)]))
            n_cut = 0
            for title, sel in cats:
                if sel:
                    contact_pages(pdf, title, sel, files, args.outdir, src_ra,
                                  src_dec)
                    n_cut += len(sel)
    print('[summary] %d cutouts written to %s' % (n_cut, pdf_path))
    return 0


if __name__ == '__main__':
    sys.exit(main())
