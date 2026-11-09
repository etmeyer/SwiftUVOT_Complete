#!/usr/bin/env python3
"""
swift_uvot_flux.py -- Step 8: flux densities, dereddened.

Turns the count rates of the Step 7 light curve into flux densities
corrected for Galactic dust and, for a distant source, for the hydrogen of
the intergalactic medium (IGM). Steps 1-7 are not changed: the count rate
stays the measured quantity, and everything here depends on assumptions
(the spectrum, the dust law, the IGM model) written into the table header.

Why. A UVOT filter counts photons over a wide band, so the factor from
count rate to flux density depends on the shape of the spectrum across the
band; uvotsource's factors (Step 7) are for a smooth, GRB-like spectrum.
Dust dims the blue end of a broad filter more than the red end, so the
extinction of a filter also depends on the spectrum, most in the UV filters
that straddle the 2175 A bump. And above redshift ~0.3 the Lyman-alpha
forest and the Lyman limit remove light inside the UV filters; above z ~ 2
also inside U and B. For such a source a UV filter may count mostly the
red tail of its response (the red leak).

How. A spectrum (the source's own shape times the IGM transmission times
the Galactic extinction) is folded through the filter's effective area
from the CALDB (the newest SPECRESP file of caldb.indx) to give the count
rate it would produce; the measured rate fixes its normalisation. The flux
density is the intrinsic spectrum (without Galactic dust and IGM) at the
filter's reference wavelength (the one implied by uvotsource's two factors,
so the numbers compare with Step 7's directly). Folding reproduces
uvotsource's factors within 0.1-2.8 % for f_nu ~ nu^-alpha, alpha 0-1.

    default      keep uvotsource's factors; only correct for dust (and the
                 IGM if --z is given), with A and the transmission folded
                 for f_nu ~ nu^-1
    --spectrum   powerlaw:ALPHA (f_nu ~ nu^-ALPHA), blackbody:T (kelvin,
                 rest frame) or file:PATH (two columns, rest-frame
                 wavelength in A and f_lambda in any units): fold that
                 spectrum for the factor too

E(B-V) comes from the IRSA dust service (Schlafly & Finkbeiner 2011, i.e.
SFD x 0.86, by default; --ebv-map sfd for Schlegel et al. 1998), cached in
UVOT_output/uvot_ebv.txt; --ebv gives a value (0 to skip dereddening).
Its error (16 % of E(B-V), the SFD map's own; --ebv-err to change) is a
separate systematic, a_gal_err, not in the flux errors: it moves every
point of a filter together. Dust law: Fitzpatrick
(1999) by default, or Cardelli, Clayton & Mathis (1989); R_V 3.1.

IGM (--z): Inoue et al. (2014) by default, or Madau (1995). Both are mean
transmissions; real sightlines scatter around them, by a lot where the
transmission is low, so rows with igm_trans < 0.9 are flagged
igm_absorbed and their values are model-dependent. red_leak marks rows
where more than 30 % of the model counts come from beyond the long end of
the filter's main band (where the effective area falls below 10 % of its
peak): the measurement mostly tells about longer wavelengths.

Output (in --outdir):
    uvot_flux.txt    one row per light-curve row: reference wavelength,
                     A_gal, IGM transmission, red-leak share, the factor
                     relative to Step 7's flux, flux densities (erg/s/cm^2/A,
                     mJy, nu F_nu) and dereddened magnitudes; upper limits
                     converted the same way
    uvot_flux.pdf    dereddened nu F_nu against time per filter, and the
                     spectral energy distribution of every observation with
                     two or more filters

Usage:
    swift_uvot_flux.py
    swift_uvot_flux.py --ebv 0.05 --law ccm89
    swift_uvot_flux.py --z 3.104 --spectrum powerlaw:1.5
"""

import argparse
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits
from scipy.interpolate import CubicSpline

import swift_uvot_env as uenv
from swift_uvot_tables import read_table, write_table

FILTERS = ['V', 'B', 'U', 'UVW1', 'UVM2', 'UVW2']
H_CGS = 6.62607015e-27       # erg s
C_AA = 2.99792458e18         # speed of light, A/s
K_CGS = 1.380649e-16         # erg/K

IRSA_DUST = 'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust'
EBV_FILE = 'uvot_ebv.txt'
EBV_REL_ERR = 0.16           # SFD's quoted relative error
DEFAULT_RV = 3.1
DEFAULT_ALPHA = 1.0          # f_nu ~ nu^-1 for the default mode's folding
IGM_FLAG = 0.9               # transmission below this: igm_absorbed
RED_FLAG = 0.30              # red-leak share above this: red_leak
MAIN_BAND = 0.10             # main band: effective area >= 10 % of its peak

COLUMNS = ['obsid', 'filter', 'date_mid', 'mjd_mid', 'status', 'rate',
           'rate_err', 'rate_limit', 'lam_ref', 'a_gal', 'a_gal_err', 'igm_trans',
           'red_share', 'factor', 'flux_aa', 'flux_aa_err', 'flux_mjy',
           'flux_mjy_err', 'nufnu', 'nufnu_err', 'mag_dered', 'ab_mag_dered',
           'mag_err', 'flags', 'reason']


# --- dust laws --------------------------------------------------------------

def _fm_uv(x, rv):
    """Fitzpatrick (1999) UV curve, k = E(lambda-V)/E(B-V), x in 1/micron."""
    x0, gamma, c3, c4 = 4.596, 0.99, 3.23, 0.41
    c2 = -0.824 + 4.717 / rv
    c1 = 2.030 - 3.007 * c2
    d = x ** 2 / ((x ** 2 - x0 ** 2) ** 2 + (x * gamma) ** 2)
    y = np.clip(x - 5.9, 0, None)
    f = 0.5392 * y ** 2 + 0.05644 * y ** 3
    return c1 + c2 * x + c3 * d + c4 * f


def a_over_ebv_f99(lam, rv=DEFAULT_RV):
    """
    A(lambda)/E(B-V) of Fitzpatrick (1999): the FM90 form in the UV (x >=
    1/0.27 micron^-1), a natural cubic spline through his optical/IR anchor
    points (the UV ones from the FM90 form) below. lam in A.
    """
    x = 1e4 / np.asarray(lam, dtype=float)
    anchors_x = np.array([0.0, 1e4 / 26500, 1e4 / 12200, 1e4 / 6000,
                          1e4 / 5470, 1e4 / 4670, 1e4 / 4110, 1e4 / 2700,
                          1e4 / 2600])
    uv = _fm_uv(anchors_x[-2:], rv) + rv
    anchors_y = np.array([
        0.0, 0.26469 * rv / 3.1, 0.82925 * rv / 3.1,
        -0.422809 + 1.00270 * rv + 2.13572e-4 * rv ** 2,
        -5.13540e-2 + 1.00216 * rv - 7.35778e-5 * rv ** 2,
        0.700127 + 1.00184 * rv - 3.32598e-5 * rv ** 2,
        1.19456 + 1.01707 * rv - 5.46959e-3 * rv ** 2
        + 7.97809e-4 * rv ** 3 - 4.45636e-5 * rv ** 4,
        uv[0], uv[1]])
    spline = CubicSpline(anchors_x, anchors_y, bc_type='natural')
    out = spline(x)
    hi = x >= 1e4 / 2700
    out[hi] = _fm_uv(x[hi], rv) + rv
    return out


def a_over_ebv_ccm89(lam, rv=DEFAULT_RV):
    """A(lambda)/E(B-V) of Cardelli, Clayton & Mathis (1989); lam in A."""
    x = 1e4 / np.asarray(lam, dtype=float)
    a = np.zeros_like(x)
    b = np.zeros_like(x)
    ir = x < 1.1
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61
    opt = (x >= 1.1) & (x < 3.3)
    y = x[opt] - 1.82
    a[opt] = np.polyval([0.32999, -0.77530, 0.01979, 0.72085, -0.02427,
                         -0.50447, 0.17699, 1.0], y)
    b[opt] = np.polyval([-2.09002, 5.30260, -0.62251, -5.38434, 1.07233,
                         2.28305, 1.41338, 0.0], y)
    uv = x >= 3.3
    xu = x[uv]
    yy = np.clip(xu - 5.9, 0, None)
    fa = -0.04473 * yy ** 2 - 0.009779 * yy ** 3
    fb = 0.2130 * yy ** 2 + 0.1207 * yy ** 3
    a[uv] = 1.752 - 0.316 * xu - 0.104 / ((xu - 4.67) ** 2 + 0.341) + fa
    b[uv] = -3.090 + 1.825 * xu + 1.206 / ((xu - 4.62) ** 2 + 0.263) + fb
    return (a + b / rv) * rv


LAWS = {'f99': a_over_ebv_f99, 'ccm89': a_over_ebv_ccm89}


# --- IGM --------------------------------------------------------------------

# Inoue et al. (2014), Table 2: Lyman series wavelength (A), the Lyman-alpha
# forest (LAF) coefficients A_j1..3 and the damped-Lyman-alpha system (DLA)
# coefficients A_j1..2.
INOUE14 = np.array([
    (1215.670, 1.68976E-02, 2.35379E-03, 1.02611E-04, 1.61698E-04, 5.38995E-05),
    (1025.720, 4.69229E-03, 6.53625E-04, 2.84940E-05, 1.54539E-04, 5.15129E-05),
    (972.537, 2.23898E-03, 3.11884E-04, 1.35962E-05, 1.49767E-04, 4.99222E-05),
    (949.743, 1.31901E-03, 1.83735E-04, 8.00974E-06, 1.46031E-04, 4.86769E-05),
    (937.803, 8.70656E-04, 1.21280E-04, 5.28707E-06, 1.42893E-04, 4.76312E-05),
    (930.748, 6.17843E-04, 8.60640E-05, 3.75186E-06, 1.40159E-04, 4.67196E-05),
    (926.226, 4.60924E-04, 6.42055E-05, 2.79897E-06, 1.37714E-04, 4.59048E-05),
    (923.150, 3.56887E-04, 4.97135E-05, 2.16720E-06, 1.35495E-04, 4.51650E-05),
    (920.963, 2.84278E-04, 3.95992E-05, 1.72628E-06, 1.33452E-04, 4.44841E-05),
    (919.352, 2.31771E-04, 3.22851E-05, 1.40743E-06, 1.31561E-04, 4.38536E-05),
    (918.129, 1.92348E-04, 2.67936E-05, 1.16804E-06, 1.29785E-04, 4.32617E-05),
    (917.181, 1.62155E-04, 2.25878E-05, 9.84689E-07, 1.28117E-04, 4.27056E-05),
    (916.429, 1.38498E-04, 1.92925E-05, 8.41033E-07, 1.26540E-04, 4.21799E-05),
    (915.824, 1.19611E-04, 1.66615E-05, 7.26340E-07, 1.25041E-04, 4.16804E-05),
    (915.329, 1.04314E-04, 1.45306E-05, 6.33446E-07, 1.23614E-04, 4.12046E-05),
    (914.919, 9.17397E-05, 1.27791E-05, 5.57091E-07, 1.22248E-04, 4.07494E-05),
    (914.576, 8.12784E-05, 1.13219E-05, 4.93564E-07, 1.20938E-04, 4.03127E-05),
    (914.286, 7.25069E-05, 1.01000E-05, 4.40299E-07, 1.19681E-04, 3.98938E-05),
    (914.039, 6.50549E-05, 9.06198E-06, 3.95047E-07, 1.18469E-04, 3.94896E-05),
    (913.826, 5.86816E-05, 8.17421E-06, 3.56345E-07, 1.17298E-04, 3.90995E-05),
    (913.641, 5.31918E-05, 7.40949E-06, 3.23008E-07, 1.16167E-04, 3.87225E-05),
    (913.480, 4.84261E-05, 6.74563E-06, 2.94068E-07, 1.15071E-04, 3.83572E-05),
    (913.339, 4.42740E-05, 6.16726E-06, 2.68854E-07, 1.14011E-04, 3.80037E-05),
    (913.215, 4.06311E-05, 5.65981E-06, 2.46733E-07, 1.12983E-04, 3.76609E-05),
    (913.104, 3.73821E-05, 5.20723E-06, 2.27003E-07, 1.11972E-04, 3.73241E-05),
    (913.006, 3.45377E-05, 4.81102E-06, 2.09731E-07, 1.11002E-04, 3.70005E-05),
    (912.918, 3.19891E-05, 4.45601E-06, 1.94255E-07, 1.10051E-04, 3.66836E-05),
    (912.839, 2.97110E-05, 4.13867E-06, 1.80421E-07, 1.09125E-04, 3.63749E-05),
    (912.768, 2.76635E-05, 3.85346E-06, 1.67987E-07, 1.08220E-04, 3.60734E-05),
    (912.703, 2.58178E-05, 3.59636E-06, 1.56779E-07, 1.07337E-04, 3.57789E-05),
    (912.645, 2.41479E-05, 3.36374E-06, 1.46638E-07, 1.06473E-04, 3.54909E-05),
    (912.592, 2.26347E-05, 3.15296E-06, 1.37450E-07, 1.05629E-04, 3.52096E-05),
    (912.543, 2.12567E-05, 2.96100E-06, 1.29081E-07, 1.04802E-04, 3.49340E-05),
    (912.499, 1.99967E-05, 2.78549E-06, 1.21430E-07, 1.03991E-04, 3.46636E-05),
    (912.458, 1.88476E-05, 2.62543E-06, 1.14452E-07, 1.03198E-04, 3.43994E-05),
    (912.420, 1.77928E-05, 2.47850E-06, 1.08047E-07, 1.02420E-04, 3.41402E-05),
    (912.385, 1.68222E-05, 2.34330E-06, 1.02153E-07, 1.01657E-04, 3.38856E-05),
    (912.353, 1.59286E-05, 2.21882E-06, 9.67268E-08, 1.00908E-04, 3.36359E-05),
    (912.324, 1.50996E-05, 2.10334E-06, 9.16925E-08, 1.00168E-04, 3.33895E-05),
])
LYMAN_LIMIT = 911.8


def _inoue_lyman_series(lo, z):
    tau = np.zeros_like(lo)
    for lj, a1, a2, a3, d1, d2 in INOUE14:
        r = lo / lj
        m = lo < lj * (1 + z)
        laf = np.where(r < 2.2, a1 * r ** 1.2,
                       np.where(r < 5.7, a2 * r ** 3.7, a3 * r ** 5.5))
        dla = np.where(r < 3.0, d1 * r ** 2, d2 * r ** 3)
        tau += np.where(m, laf + dla, 0.0)
    return tau


def _inoue_lyman_continuum(lo, z):
    r = lo / LYMAN_LIMIT
    zz = 1.0 + z
    m = lo < LYMAN_LIMIT * zz
    rr = np.where(m, r, 1.0)    # placeholder where not absorbed
    if z < 1.2:
        laf = 0.3248 * (rr ** 1.2 - zz ** -0.9 * rr ** 2.1)
    elif z < 4.7:
        laf = np.where(rr < 2.2,
                       2.545e-2 * zz ** 1.6 * rr ** 2.1 + 0.3248 * rr ** 1.2
                       - 0.2496 * rr ** 2.1,
                       2.545e-2 * (zz ** 1.6 * rr ** 2.1 - rr ** 3.7))
    else:
        laf = np.where(
            rr < 2.2,
            5.221e-4 * zz ** 3.4 * rr ** 2.1 + 0.3248 * rr ** 1.2
            - 3.140e-2 * rr ** 2.1,
            np.where(rr < 5.7,
                     5.221e-4 * zz ** 3.4 * rr ** 2.1 + 0.2182 * rr ** 2.1
                     - 2.545e-2 * rr ** 3.7,
                     5.221e-4 * (zz ** 3.4 * rr ** 2.1 - rr ** 5.5)))
    if z < 2.0:
        dla = 0.2113 * zz ** 2 - 0.07661 * zz ** 2.3 * rr ** -0.3 \
            - 0.1347 * rr ** 2
    else:
        dla = np.where(
            rr < 3.0,
            0.6340 + 0.04696 * zz ** 3 - 0.01779 * zz ** 3.3 * rr ** -0.3
            - 0.1347 * rr ** 2 - 0.2905 * rr ** -0.3,
            0.04696 * zz ** 3 - 0.01779 * zz ** 3.3 * rr ** -0.3
            - 0.02916 * rr ** 3)
    return np.where(m, laf + dla, 0.0)


def igm_inoue14(lam_obs, z):
    """Mean IGM transmission of Inoue et al. (2014); observed wavelength in A."""
    lo = np.asarray(lam_obs, dtype=float)
    tau = _inoue_lyman_series(lo, z) + _inoue_lyman_continuum(lo, z)
    return np.exp(-tau)


def igm_madau95(lam_obs, z):
    """
    Mean IGM transmission of Madau (1995): Lyman alpha to delta (his eq. 15
    and Table 2) and the Lyman continuum (the approximation of his footnote
    3); observed wavelength in A.
    """
    lo = np.asarray(lam_obs, dtype=float)
    tau = np.zeros_like(lo)
    for lj, aj in ((1215.67, 3.6e-3), (1025.72, 1.7e-3), (972.537, 1.2e-3),
                   (949.743, 9.3e-4)):
        m = lo < lj * (1 + z)
        tau[m] += aj * (lo[m] / lj) ** 3.46
    zz = 1 + z
    m = lo < 912.0 * zz
    x = lo[m] / 912.0
    tau[m] += (0.25 * x ** 3 * (zz ** 0.46 - x ** 0.46)
               + 9.4 * x ** 1.5 * (zz ** 0.18 - x ** 0.18)
               - 0.7 * x ** 3 * (x ** -1.32 - zz ** -1.32)
               - 0.023 * (zz ** 1.68 - x ** 1.68))
    return np.exp(-tau)


IGM_MODELS = {'inoue14': igm_inoue14, 'madau95': igm_madau95}


# --- filters and spectra ----------------------------------------------------

def effective_areas(caldb=None):
    """
    {filter: (wavelength A, bin width A, area cm^2, file name)} from the
    newest good SPECRESP file of each filter in the UVOT caldb.indx (no
    HEASoft needed).
    """
    caldb = caldb or os.environ.get('CALDB') or uenv.DEFAULT_HEASOFT_CALDB
    index = os.path.join(caldb, uenv.UVOT_CALDB_INDEX)
    if not os.path.exists(index):
        raise FileNotFoundError('no UVOT CALDB index at %s (run heainit, or '
                                'set $CALDB)' % index)
    best = {}
    with fits.open(index) as hdul:
        for r in hdul[1].data:
            if r['CAL_CNAM'].strip() != 'SPECRESP' or r['CAL_QUAL'] != 0:
                continue
            filt = r['FILTER'].strip()
            if filt not in FILTERS:
                continue
            key = (r['CAL_VSD'], r['CAL_DATE'])
            if filt not in best or key > best[filt][0]:
                best[filt] = (key, os.path.join(caldb, r['CAL_DIR'].strip(),
                                                r['CAL_FILE'].strip()))
    out = {}
    for filt, (_, path) in best.items():
        with fits.open(path) as hdul:
            d = hdul[1].data
            lo, hi = d['WAVE_MIN'].astype(float), d['WAVE_MAX'].astype(float)
            lam = 0.5 * (lo + hi)
            order = np.argsort(lam)
            out[filt] = (lam[order], np.abs(hi - lo)[order],
                         d['SPECRESP'].astype(float)[order],
                         os.path.basename(path))
    return out


def parse_spectrum(text):
    """
    A function f_lambda(rest-frame wavelength in A), arbitrary units, and a
    label, from 'powerlaw:ALPHA', 'blackbody:T' or 'file:PATH'.
    """
    kind, _, value = text.partition(':')
    kind = kind.strip().lower()
    if kind == 'powerlaw':
        alpha = float(value)
        return (lambda lam: (np.asarray(lam, float) / 5000.0) ** (alpha - 2),
                'power law f_nu ~ nu^-%g' % alpha)
    if kind == 'blackbody':
        temp = float(value)
        if temp <= 0:
            raise ValueError('temperature must be positive')

        def bb(lam):
            lam = np.asarray(lam, float)
            arg = np.clip(H_CGS * C_AA / (lam * K_CGS * temp), None, 700)
            return (lam / 5000.0) ** -5 / np.expm1(arg)
        return bb, 'blackbody T = %g K' % temp
    if kind == 'file':
        data = np.loadtxt(value, usecols=(0, 1), comments='#', ndmin=2)
        data = data[np.argsort(data[:, 0])]
        if len(data) < 2 or np.any(data[:, 0] <= 0):
            raise ValueError('%s: need two or more rows of positive '
                             'wavelengths' % value)
        lam0, flam0 = data[:, 0], data[:, 1]

        def tab(lam):
            # zero outside the tabulated range
            return np.interp(lam, lam0, flam0, left=0.0, right=0.0)
        return tab, 'spectrum file %s (%.0f-%.0f A rest frame)' % (
            value, lam0[0], lam0[-1])
    raise ValueError("give powerlaw:ALPHA, blackbody:T or file:PATH")


def fold(ea, shape, z, law, ebv, rv, igm, lam_ref):
    """
    For one filter and one rest-frame spectrum shape: the counts per unit
    normalisation with and without the IGM and Galactic dust, the
    resulting A_gal (mag), IGM transmission, red-leak share, and the
    intrinsic flux density at lam_ref per unit normalisation.
    """
    lam, dl, area, _ = ea
    s = shape(lam / (1 + z))
    t = igm(lam, z) if (igm and z > 0) else np.ones_like(lam)
    ext = 10 ** (-0.4 * law(lam, rv) * ebv) if ebv > 0 else np.ones_like(lam)
    weight = lam / (H_CGS * C_AA) * area * dl     # photons per (f_lambda dlam)
    c_bare = np.sum(s * weight)
    c_igm = np.sum(s * t * weight)
    c_full = np.sum(s * t * ext * weight)
    a_gal = -2.5 * math.log10(c_full / c_igm) + 0.0 if c_igm > 0 \
        else float('nan')    # + 0.0: no '-0.000' when there is no dust
    igm_trans = c_igm / c_bare if c_bare > 0 else float('nan')
    peak = area.max()
    red_edge = lam[area >= MAIN_BAND * peak].max()
    counts = s * t * ext * weight
    red_share = counts[lam > red_edge].sum() / c_full if c_full > 0 \
        else float('nan')
    s_ref = float(shape(np.array([lam_ref / (1 + z)]))[0])
    return {'c_full': c_full, 'a_gal': a_gal, 'igm_trans': igm_trans,
            'red_share': red_share, 's_ref': s_ref}


# --- E(B-V) -----------------------------------------------------------------

def _irsa_value(text):
    m = re.search(r'([-+\d.eE]+)', text or '')
    return float(m.group(1)) if m else float('nan')


def query_irsa(ra, dec, timeout=60):
    """E(B-V) at (ra, dec) from the IRSA dust service (the map pixel at
    the position): {'sf11': value, 'sfd': value} in mag."""
    import requests
    resp = requests.get(IRSA_DUST, params={
        'locstr': '%.6f %+.6f equ j2000' % (ra, dec),
        'regSize': '2.0'}, timeout=timeout)    # degrees; IRSA's minimum
    resp.raise_for_status()
    root = ET.fromstring(resp.content)
    if root.attrib.get('status', 'ok') != 'ok':
        raise RuntimeError('IRSA: %s' % ''.join(root.itertext()).strip()[:200])
    out = {}
    for key, tag in (('sf11', 'SandF'), ('sfd', 'SFD')):
        val = float('nan')
        for el in root.iter():
            if el.tag == 'refPixelValue' + tag:
                val = _irsa_value(el.text)
        out[key] = val
    if not all(np.isfinite(v) for v in out.values()):
        raise RuntimeError('IRSA answer without E(B-V) values')
    return out


def ebv_for(ra, dec, outdir, which, timeout=60):
    """
    (E(B-V), source text) from the cache UVOT_output/uvot_ebv.txt, or IRSA
    (then cached).
    """
    path = os.path.join(outdir, EBV_FILE)
    vals = None
    if os.path.exists(path):
        rows = read_table(path)
        if rows and abs(float(rows[0]['ra']) - ra) < 1e-4 and \
                abs(float(rows[0]['dec']) - dec) < 1e-4:
            vals = {r['map']: float(r['ebv']) for r in rows}
            src = 'IRSA dust service (cached in %s)' % path
    if vals is None:
        vals = query_irsa(ra, dec, timeout)
        write_table(path, [{'map': k, 'ra': '%.6f' % ra, 'dec': '%.6f' % dec,
                            'ebv': '%.4f' % v} for k, v in vals.items()],
                    ['map', 'ra', 'dec', 'ebv'],
                    ['E(B-V) from %s, %s; sf11 = Schlafly & Finkbeiner 2011, '
                     'sfd = Schlegel et al. 1998; the map pixel at the source'
                     % (IRSA_DUST, datetime.now(timezone.utc)
                        .strftime('%Y-%m-%d'))])
        src = 'IRSA dust service (now cached in %s)' % path
    return vals[which], src


# --- light curve in, flux out -----------------------------------------------

def _num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float('nan')


def read_lightcurve(path):
    """Rows, (ra, dec), and Step 7's per-filter factors from its header."""
    rows = read_table(path)
    factors, radec = {}, None
    with open(path) as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            m = re.search(r'source: RA ([\d.+-]+) Dec ([\d.+-]+)', line)
            if m:
                radec = (float(m.group(1)), float(m.group(2)))
            m = re.match(r'# (\w+): Vega zero point ([\d.]+), AB ([\d.]+), '
                         r'([\d.eE+-]+) erg/s/cm\^2/A and ([\d.eE+-]+) mJy',
                         line)
            if m:
                factors[m.group(1)] = {
                    'zp_vega': float(m.group(2)), 'zp_ab': float(m.group(3)),
                    'fcf_aa': float(m.group(4)), 'fcf_mjy': float(m.group(5))}
    for fac in factors.values():
        # the wavelength at which uvotsource's two factors agree
        fac['lam_ref'] = math.sqrt(C_AA * fac['fcf_mjy'] * 1e-26
                                   / fac['fcf_aa'])
    return rows, radec, factors


def filter_terms(filters, areas, factors, shape, z, law, ebv, rv, igm,
                 use_shape):
    """Per filter: lam_ref, A_gal, IGM transmission, red-leak share and the
    factor from count rate to intrinsic f_lambda at lam_ref."""
    terms = {}
    for f in filters:
        if f not in areas or f not in factors:
            continue
        fac = factors[f]
        r = fold(areas[f], shape, z, law, ebv, rv, igm, fac['lam_ref'])
        if use_shape:
            fcf = r['s_ref'] / r['c_full']    # rate -> intrinsic f_lambda
        else:
            fcf = fac['fcf_aa'] * 10 ** (0.4 * r['a_gal']) / r['igm_trans']
        terms[f] = dict(r, lam_ref=fac['lam_ref'], fcf=fcf,
                        factor=fcf / fac['fcf_aa'])
    return terms


def convert(row, term, fac, a_err):
    """Table strings for one light-curve row."""
    out = {k: row.get(k) for k in ('obsid', 'filter', 'date_mid', 'mjd_mid',
                                   'status', 'rate', 'rate_err',
                                   'rate_limit')}
    flags = []
    if term is None:
        out.update(flags='-', reason='no conversion for this filter')
        return out
    lam = term['lam_ref']
    out.update(lam_ref='%.0f' % lam, a_gal='%.3f' % term['a_gal'],
               a_gal_err='%.3f' % a_err, igm_trans='%.3f' % term['igm_trans'],
               red_share='%.2f' % term['red_share'],
               factor='%.4f' % term['factor'])
    if term['igm_trans'] < IGM_FLAG:
        flags.append('igm_absorbed')
    if term['red_share'] > RED_FLAG:
        flags.append('red_leak')
    status = row.get('status')
    value = _num(row.get('rate_limit')) if status == 'upper_limit' \
        else _num(row.get('rate'))
    if status == 'lower_limit':
        value = _num(row.get('rate'))
    err = _num(row.get('rate_err')) if status == 'detected' else float('nan')
    reasons = []
    if term['igm_trans'] < IGM_FLAG:
        reasons.append('IGM keeps %.0f %% of the counts' %
                       (100 * term['igm_trans']))
    if term['red_share'] > RED_FLAG:
        reasons.append('%.0f %% of the counts from beyond %s\'s main band'
                       % (100 * term['red_share'], row.get('filter')))
    if status in ('detected', 'upper_limit', 'lower_limit') and value > 0:
        flam = value * term['fcf']
        fnu_mjy = flam * lam ** 2 / C_AA * 1e26
        out['flux_aa'] = '%.4e' % flam
        out['flux_mjy'] = '%.5g' % fnu_mjy
        out['nufnu'] = '%.4e' % (flam * lam)
        mag_shift = -2.5 * math.log10(term['factor'])
        out['mag_dered'] = '%.3f' % (fac['zp_vega'] - 2.5 * math.log10(value)
                                     + mag_shift)
        out['ab_mag_dered'] = '%.3f' % (fac['zp_ab'] - 2.5 * math.log10(value)
                                        + mag_shift)
        if np.isfinite(err):
            rel = err / value
            out['flux_aa_err'] = '%.2e' % (flam * rel)
            out['flux_mjy_err'] = '%.3g' % (fnu_mjy * rel)
            out['nufnu_err'] = '%.2e' % (flam * lam * rel)
            out['mag_err'] = '%.3f' % (2.5 / math.log(10) * rel)
        if status == 'upper_limit':
            reasons.insert(0, 'flux columns hold the upper limit')
        elif status == 'lower_limit':
            reasons.insert(0, 'flux columns hold the saturated lower limit')
    out['flags'] = ','.join(flags) if flags else '-'
    out['reason'] = '; '.join(reasons) if reasons else '-'
    return out


def plot(pdf_path, rows, label):
    """Dereddened nu F_nu against time per filter, then the SEDs."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    colors = {'V': '#2e7d32', 'B': '#1565c0', 'U': '#6a1b9a',
              'UVW1': '#c62828', 'UVM2': '#ef6c00', 'UVW2': '#4e342e'}
    with PdfPages(pdf_path) as pdf:
        filts = [f for f in FILTERS if any(
            r['filter'] == f and r.get('nufnu') not in (None, '-')
            for r in rows)]
        if filts:
            fig, axes = plt.subplots(len(filts), 1, sharex=True,
                                     figsize=(8.5, 1.8 + 1.6 * len(filts)),
                                     squeeze=False)
            for ax, f in zip(axes[:, 0], filts):
                for status, marker in (('detected', 'o'),
                                       ('upper_limit', 'v')):
                    pts = [r for r in rows if r['filter'] == f
                           and r['status'] == status
                           and r.get('nufnu') not in (None, '-')]
                    if not pts:
                        continue
                    t = [_num(r['mjd_mid']) for r in pts]
                    y = [_num(r['nufnu']) for r in pts]
                    e = [_num(r.get('nufnu_err')) for r in pts]
                    flagged = any(r.get('flags', '-') != '-' for r in pts)
                    face = 'none' if flagged else colors[f]
                    if status == 'detected':
                        ax.errorbar(t, y, yerr=e, fmt=marker, ms=3.5,
                                    color=colors[f], mfc=face, lw=0.8)
                    else:
                        ax.plot(t, y, marker, ms=4, color=colors[f], mfc=face,
                                alpha=0.6)
                ax.set_ylabel('%s\nnuFnu' % f, fontsize=8)
                ax.set_yscale('log')
            axes[-1, 0].set_xlabel('MJD')
            fig.suptitle('%s: dereddened nu F_nu (erg/s/cm^2); open: '
                         'flagged; triangles: upper limits' % label,
                         fontsize=9)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
        by_obs = {}
        for r in rows:
            if r['status'] == 'detected' and r.get('nufnu') not in (None, '-'):
                by_obs.setdefault(r['obsid'], []).append(r)
        seds = {k: v for k, v in by_obs.items() if len(v) >= 2}
        if seds:
            fig, ax = plt.subplots(figsize=(8.5, 6))
            cmap = plt.get_cmap('viridis')
            mjds = [np.mean([_num(r['mjd_mid']) for r in v])
                    for v in seds.values()]
            lo, hi = min(mjds), max(mjds)
            for (obsid, pts), mjd in zip(sorted(seds.items()), mjds):
                pts = sorted(pts, key=lambda r: _num(r['lam_ref']))
                nu = [C_AA / _num(r['lam_ref']) for r in pts]
                y = [_num(r['nufnu']) for r in pts]
                e = [_num(r.get('nufnu_err')) for r in pts]
                col = cmap((mjd - lo) / (hi - lo) if hi > lo else 0.5)
                ax.errorbar(nu, y, yerr=e, fmt='-o', ms=3, lw=0.6, color=col)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel('frequency (Hz, observed)')
            ax.set_ylabel('nu F_nu (erg/s/cm^2), dereddened')
            ax.set_title('%s: one line per observation with 2+ detected '
                         'filters, colour = time (MJD %.0f-%.0f)'
                         % (label, lo, hi), fontsize=9)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def main(argv=None):
    sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser(
        description='Step 8: flux densities, dereddened.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--outdir', default='UVOT_output',
                        help='Steps 3-7 output folder (default: UVOT_output)')
    parser.add_argument('--ebv', type=float, default=None,
                        help='Galactic E(B-V) in mag (default: from IRSA); '
                             '0 to skip dereddening')
    parser.add_argument('--ebv-err', type=float, default=None,
                        help='its error (default: 16 %% of E(B-V))')
    parser.add_argument('--ebv-map', choices=('sf11', 'sfd'), default='sf11',
                        help='IRSA map: Schlafly & Finkbeiner 2011 (sf11, '
                             'default) or Schlegel et al. 1998 (sfd)')
    parser.add_argument('--law', choices=sorted(LAWS), default='f99',
                        help='dust law: Fitzpatrick 1999 (f99, default) or '
                             'Cardelli et al. 1989 (ccm89)')
    parser.add_argument('--rv', type=float, default=DEFAULT_RV,
                        help='R_V (default: 3.1)')
    parser.add_argument('--z', type=float, default=0.0,
                        help='redshift of the source, for the IGM (default: '
                             '0, no IGM)')
    parser.add_argument('--igm', choices=sorted(IGM_MODELS) + ['none'],
                        default='inoue14',
                        help='IGM model when --z is given: Inoue et al. 2014 '
                             '(inoue14, default), Madau 1995 (madau95), none')
    parser.add_argument('--spectrum', default=None,
                        help='powerlaw:ALPHA (f_nu ~ nu^-ALPHA), '
                             'blackbody:T or file:PATH (rest-frame A, '
                             'f_lambda); default: keep uvotsource\'s '
                             'factors')
    parser.add_argument('--irsa-timeout', type=float, default=60,
                        help='seconds to wait for IRSA (default: 60)')
    args = parser.parse_args(argv)
    if args.z < 0 or args.z > 10:
        print('ERROR: --z must be between 0 and 10', file=sys.stderr)
        return 1
    if not 2.0 <= args.rv <= 6.0:
        print('ERROR: --rv must be between 2 and 6', file=sys.stderr)
        return 1

    lc_path = os.path.join(args.outdir, 'uvot_lightcurve.txt')
    if not os.path.exists(lc_path):
        print('ERROR: %s not found; run Step 7 (swift_uvot_lightcurve.py) '
              'first.' % lc_path, file=sys.stderr)
        return 1
    rows, radec, factors = read_lightcurve(lc_path)
    if not factors:
        print('ERROR: %s has no conversion factors in its header (re-run '
              'Step 7)' % lc_path, file=sys.stderr)
        return 1

    if args.ebv is not None:
        if args.ebv < 0:
            print('ERROR: --ebv must not be negative', file=sys.stderr)
            return 1
        ebv = args.ebv
        ebv_err = args.ebv_err if args.ebv_err is not None \
            else EBV_REL_ERR * ebv
        ebv_src = 'given (--ebv)'
    else:
        if radec is None:
            print('ERROR: no source position in %s; give --ebv'
                  % lc_path, file=sys.stderr)
            return 1
        try:
            ebv, ebv_src = ebv_for(radec[0], radec[1], args.outdir,
                                   args.ebv_map, args.irsa_timeout)
        except Exception as exc:
            print('ERROR: E(B-V) from IRSA failed (%s); give --ebv'
                  % exc, file=sys.stderr)
            return 1
        ebv_err = args.ebv_err if args.ebv_err is not None \
            else EBV_REL_ERR * ebv
        ebv_src += ', map %s' % args.ebv_map

    try:
        areas = effective_areas()
    except Exception as exc:
        print('ERROR: %s' % exc, file=sys.stderr)
        return 1
    if args.spectrum:
        try:
            shape, shape_label = parse_spectrum(args.spectrum)
        except (ValueError, OSError) as exc:
            print('ERROR: --spectrum %s: %s' % (args.spectrum, exc),
                  file=sys.stderr)
            return 1
    else:
        shape, _ = parse_spectrum('powerlaw:%g' % DEFAULT_ALPHA)
        shape_label = ("uvotsource's factors; corrections folded for "
                       "f_nu ~ nu^-%g" % DEFAULT_ALPHA)
    igm = None if args.igm == 'none' or args.z <= 0 else IGM_MODELS[args.igm]
    law = LAWS[args.law]
    filters = sorted({r['filter'] for r in rows},
                     key=lambda f: FILTERS.index(f) if f in FILTERS else 99)
    terms = filter_terms(filters, areas, factors, shape, args.z, law, ebv,
                         args.rv, igm, bool(args.spectrum))

    out = []
    for r in rows:
        f = r['filter']
        term = terms.get(f)
        a_err = term['a_gal'] / ebv * ebv_err if (term and ebv > 0) else 0.0
        out.append(convert(r, term, factors.get(f), a_err))

    prov = uenv.provenance()
    label = 'RA %.5f Dec %.5f' % radec if radec else '?'
    comments = [
        'swift_uvot_flux.py, %s' % datetime.now(timezone.utc)
        .strftime('%Y-%m-%dT%H:%M:%SZ'),
        'command: %s' % ' '.join(sys.argv),
        'source: %s' % label,
        'pipeline %s' % prov['PIPECOMM'],
        'input: %s (Step 7); rate, rate_err as there' % lc_path,
        'E(B-V) = %.4f +- %.4f mag, %s; dust law %s, R_V %.2f'
        % (ebv, ebv_err, ebv_src, args.law, args.rv),
        'spectrum: %s' % shape_label,
        'IGM: %s' % ('none' if igm is None else '%s at z = %g'
                     % (args.igm, args.z)),
        'lam_ref (A): wavelength of uvotsource\'s factors; flux_aa '
        '(erg/s/cm^2/A), flux_mjy (mJy), nufnu (erg/s/cm^2): intrinsic '
        'spectrum at lam_ref, without Galactic dust and IGM',
        'factor: flux_aa / Step 7 flux_aa; a_gal (mag): folded Galactic '
        'extinction, a_gal_err its systematic error (not in the flux '
        'errors); igm_trans: folded IGM transmission of the counts',
        'red_share: share of the model counts from beyond the filter\'s main '
        'band (effective area >= %g %% of peak); flags igm_absorbed '
        '(igm_trans < %g), red_leak (red_share > %g)'
        % (100 * MAIN_BAND, IGM_FLAG, RED_FLAG),
        'mag_dered, ab_mag_dered: Step 7 magnitudes minus the same correction; '
        'upper_limit rows: flux and magnitude columns hold the limit',
    ]
    for f in filters:
        if f in terms:
            t = terms[f]
            comments.append(
                '%s: %s, lam_ref %.0f A, A_gal %.3f mag (A/E(B-V) %s), '
                'igm_trans %.3f, red_share %.2f, factor %.4f'
                % (f, areas[f][3], t['lam_ref'], t['a_gal'],
                   '%.2f' % (t['a_gal'] / ebv) if ebv > 0 else '-',
                   t['igm_trans'], t['red_share'], t['factor']))
    flux_path = os.path.join(args.outdir, 'uvot_flux.txt')
    write_table(flux_path, out, COLUMNS, comments)
    pdf_path = os.path.join(args.outdir, 'uvot_flux.pdf')
    plot(pdf_path, out, label)

    print('[info] source %s; E(B-V) %.4f +- %.4f (%s)'
          % (label, ebv, ebv_err, ebv_src))
    print('[info] spectrum: %s; IGM: %s' % (
        shape_label, 'none' if igm is None else '%s, z = %g'
        % (args.igm, args.z)))
    print()
    print('[summary] filter  lam_ref  A_gal  A/E(B-V)  igm_trans  red_share  '
          'factor  rows  flagged')
    for f in filters:
        if f not in terms:
            print('          %-6s  no effective area or Step 7 factor' % f)
            continue
        t = terms[f]
        n = sum(r['filter'] == f for r in out)
        nf = sum(r['filter'] == f and r.get('flags', '-') != '-' for r in out)
        print('          %-6s  %7.0f  %5.3f  %8s  %9.3f  %9.2f  %6.3f  %4d  '
              '%7d' % (f, t['lam_ref'], t['a_gal'],
                       '%.2f' % (t['a_gal'] / ebv) if ebv > 0 else '-',
                       t['igm_trans'], t['red_share'], t['factor'], n, nf))
    print('[summary] Written: %s, %s' % (flux_path, pdf_path))
    return 0


if __name__ == '__main__':
    sys.exit(main())
