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

Host dust (--host-ebv, default 0): dust in the source's galaxy, applied at
the rest-frame wavelength inside the same folding: Gordon et al. (2003) SMC
bar (default; the usual choice for quasars) or LMC average, Fitzpatrick
(1999) R_V 3.1 (mw), or Calzetti et al. (2000) for star-forming hosts. The
curves end at 1000 A (Calzetti 1200 A) and are held constant below.

--fit-alpha: for each observation with 3 or more detected filters (not
flagged igm_absorbed or red_leak; --fit-filters to choose), fit
f_nu ~ nu^-alpha through the dust, the IGM and the folding to the measured
rates, each with rate_err and flux_sys combined, and convert that
observation's rows with the fitted alpha. uvot_alpha.txt holds the fits:
alpha, its error with and without flux_sys (the latter for changes between
observations), its shift when E(B-V) moves by its error, chi^2 and p
(poor_fit below 0.001, naming the worst filter), and the time the filters
span. Observations without a fit keep the default conversion.

flux_sys: the flux columns' calibration systematic (a fraction, not in the
errors): Step 7's mag_sys (uvotsource's zero-point error), the conversion's
own (--cal-sys; default 5 % in V, B, U and 4 % in the UV, from the star
check of known spectra) and, in UVW1 and UVW2, 0.25 and 0.6 x red_share
(the CALDB sensitivity-loss correction over-corrects the counts from the
filters' red tails), in quadrature.

Output (in --outdir):
    uvot_flux.txt    one row per light-curve row: reference wavelength,
                     A_gal, IGM transmission, red-leak share, the factor
                     relative to Step 7's flux, flux densities (erg/s/cm^2/A,
                     mJy, nu F_nu) and dereddened magnitudes; upper limits
                     converted the same way
    uvot_flux.pdf    dereddened nu F_nu against time per filter, the
                     spectral energy distribution of every observation with
                     two or more filters, and with --fit-alpha alpha against
                     time and brightness
    uvot_alpha.txt   with --fit-alpha: one row per fitted observation

Usage:
    swift_uvot_flux.py
    swift_uvot_flux.py --ebv 0.05 --law ccm89
    swift_uvot_flux.py --z 3.104 --spectrum powerlaw:1.5
    swift_uvot_flux.py --z 0.158 --fit-alpha --fit-filters V,B,UVW1,UVM2,UVW2
    swift_uvot_flux.py --z 0.5 --host-ebv 0.1 --host-law smc
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
from scipy.stats import chi2 as chi2_dist

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
ALPHA_GRID = np.round(np.arange(-3.0, 6.0001, 0.01), 2)   # --fit-alpha range
FIT_P_MIN = 1e-3             # fit p-value below this: poor_fit
FIT_MIN_FILTERS = 3
MAG_TO_REL = 2.5 / math.log(10)
# Calibration systematic of the flux conversion (fraction), from the Step 8b
# star check: hot white dwarfs read 3.5-6.5 % low in V and B, Gaia XP field
# stars up to 6 % high in B; U 5 %; UV 2-5 % (doc_notes.md, Step 8 colour
# term). Added in quadrature to Step 7's mag_sys into flux_sys.
CAL_SYS = {'V': 0.05, 'B': 0.05, 'U': 0.05, 'UVW1': 0.04, 'UVM2': 0.04,
           'UVW2': 0.04}
# Red-tail term, times red_share: SENSCORR over-corrects the counts from the
# UV filters' red tails (solar-type stars read up to ~30 % high in UVW2 after
# 2015, ~1.00 in 2009-10; weaker in UVW1, none in UVM2).
RED_TAIL_SYS = {'UVW1': 0.25, 'UVW2': 0.6}

COLUMNS = ['obsid', 'filter', 'date_mid', 'mjd_mid', 'status', 'rate',
           'rate_err', 'rate_limit', 'lam_ref', 'a_gal', 'a_gal_err', 'a_host',
           'igm_trans', 'red_share', 'alpha', 'alpha_src', 'factor', 'flux_aa', 'flux_aa_err', 'flux_mjy',
           'flux_mjy_err', 'nufnu', 'nufnu_err', 'flux_sys', 'mag_dered',
           'ab_mag_dered', 'mag_err', 'flags', 'reason']
ALPHA_COLUMNS = ['obsid', 'mjd_mid', 'n_filt', 'filters', 'alpha', 'alpha_err',
                 'alpha_err_stat', 'alpha_sys_ebv', 'chi2', 'dof', 'p_value',
                 'dt_span_h', 'nufnu_5000', 'flags', 'reason']


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


# --- host-galaxy dust (rest frame) -------------------------------------------

def _gordon03(lam, rv, uv, optnir_lam, optnir_axav):
    """
    A(lambda)/E(B-V) of a Gordon et al. (2003) average curve, built the way
    Fitzpatrick (1999) builds his: the FM90 form (uv = c1, c2, c3, c4, x0,
    gamma) for x >= 1/0.27 micron^-1, a natural cubic spline through the
    optical/IR anchors (A/A(V) at optnir_lam, micron) and two UV points
    below. Held constant shortward of 1000 A, where the curves end. lam in A.
    """
    c1, c2, c3, c4, x0, gamma = uv
    x = np.minimum(1e4 / np.asarray(lam, dtype=float), 10.0)

    def fm(xx):
        d = xx ** 2 / ((xx ** 2 - x0 ** 2) ** 2 + (xx * gamma) ** 2)
        y = np.clip(xx - 5.9, 0, None)
        k = c1 + c2 * xx + c3 * d + c4 * (0.5392 * y ** 2 + 0.05644 * y ** 3)
        return k / rv + 1.0                            # A/A(V)
    anchors_x = np.concatenate([[0.0], 1.0 / np.asarray(optnir_lam),
                                1e4 / np.array([2700.0, 2600.0])])
    anchors_y = np.concatenate([[0.0], optnir_axav,
                                fm(1e4 / np.array([2700.0, 2600.0]))])
    order = np.argsort(anchors_x)
    spline = CubicSpline(anchors_x[order], anchors_y[order], bc_type='natural')
    out = spline(x)
    hi = x >= 1e4 / 2700
    out[hi] = fm(x[hi])
    return out * rv


def a_over_ebv_smc(lam):
    """Gordon et al. (2003) SMC bar, R_V 2.74 (smoothed anchors of Gordon
    et al. 2016, as in the dust_extinction package)."""
    return _gordon03(lam, 2.74, (-4.959, 2.264, 0.389, 0.461, 4.6, 1.0),
                     [2.198, 1.65, 1.25, 0.81, 0.65, 0.55, 0.44, 0.37],
                     [0.11, 0.169, 0.25, 0.567, 0.801, 1.00, 1.374, 1.672])


def a_over_ebv_lmc(lam):
    """Gordon et al. (2003) LMC average, R_V 3.41."""
    return _gordon03(lam, 3.41, (-0.890, 0.998, 2.719, 0.400, 4.579, 0.934),
                     [2.198, 1.65, 1.25, 0.55, 0.44, 0.37],
                     [0.10, 0.186, 0.257, 1.000, 1.293, 1.518])


def a_over_ebv_calzetti(lam):
    """Calzetti et al. (2000) starburst attenuation k(lambda), R_V 4.05
    (their eq. 4); held constant outside 0.12-2.2 micron."""
    m = np.clip(np.asarray(lam, dtype=float) / 1e4, 0.12, 2.2)
    blue = 2.659 * (-2.156 + 1.509 / m - 0.198 / m ** 2 + 0.011 / m ** 3)
    red = 2.659 * (-1.857 + 1.040 / m)
    return np.where(m < 0.63, blue, red) + 4.05


HOST_LAWS = {'smc': a_over_ebv_smc, 'lmc': a_over_ebv_lmc,
             'mw': a_over_ebv_f99, 'calzetti': a_over_ebv_calzetti}


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


def fold(ea, shape, z, law, ebv, rv, igm, lam_ref, host=None, host_ebv=0.0):
    """
    For one filter and one rest-frame spectrum shape: the counts per unit
    normalisation with and without the IGM, host dust and Galactic dust,
    the resulting A_gal and A_host (mag), IGM transmission, red-leak share,
    and the intrinsic flux density at lam_ref per unit normalisation.
    """
    lam, dl, area, _ = ea
    s = shape(lam / (1 + z))
    t = igm(lam, z) if (igm and z > 0) else np.ones_like(lam)
    h = 10 ** (-0.4 * host(lam / (1 + z)) * host_ebv) \
        if (host is not None and host_ebv > 0) else np.ones_like(lam)
    ext = 10 ** (-0.4 * law(lam, rv) * ebv) if ebv > 0 else np.ones_like(lam)
    weight = lam / (H_CGS * C_AA) * area * dl     # photons per (f_lambda dlam)
    c_bare = np.sum(s * weight)
    c_igm = np.sum(s * t * weight)
    c_host = np.sum(s * t * h * weight)
    c_full = np.sum(s * t * h * ext * weight)
    # + 0.0: no '-0.000' when there is no dust
    a_gal = -2.5 * math.log10(c_full / c_host) + 0.0 if c_host > 0 \
        else float('nan')
    a_host = -2.5 * math.log10(c_host / c_igm) + 0.0 if c_igm > 0 \
        else float('nan')
    igm_trans = c_igm / c_bare if c_bare > 0 else float('nan')
    peak = area.max()
    red_edge = lam[area >= MAIN_BAND * peak].max()
    counts = s * t * h * ext * weight
    red_share = counts[lam > red_edge].sum() / c_full if c_full > 0 \
        else float('nan')
    s_ref = float(shape(np.array([lam_ref / (1 + z)]))[0])
    return {'c_full': c_full, 'a_gal': a_gal, 'a_host': a_host,
            'igm_trans': igm_trans, 'red_share': red_share, 's_ref': s_ref}


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
                 use_shape, host=None, host_ebv=0.0):
    """Per filter: lam_ref, A_gal, A_host, IGM transmission, red-leak share
    and the factor from count rate to intrinsic f_lambda at lam_ref."""
    terms = {}
    for f in filters:
        if f not in areas or f not in factors:
            continue
        fac = factors[f]
        r = fold(areas[f], shape, z, law, ebv, rv, igm, fac['lam_ref'],
                 host, host_ebv)
        if use_shape:
            fcf = r['s_ref'] / r['c_full']    # rate -> intrinsic f_lambda
        else:
            fcf = fac['fcf_aa'] * 10 ** (0.4 * (r['a_gal'] + r['a_host'])) \
                / r['igm_trans']
        terms[f] = dict(r, lam_ref=fac['lam_ref'], fcf=fcf,
                        factor=fcf / fac['fcf_aa'])
    return terms


def flux_sys(filt, mag_sys, red_share, cal_sys):
    """Relative calibration systematic of a filter's flux: Step 7's
    mag_sys (zero point), the conversion's cal_sys and the red-tail term
    (RED_TAIL_SYS x red_share) in quadrature. 0 cal_sys drops the latter two."""
    zp = mag_sys / MAG_TO_REL if np.isfinite(mag_sys) else 0.0
    cal = cal_sys.get(filt, 0.0)
    red = RED_TAIL_SYS.get(filt, 0.0) * red_share \
        if (cal > 0 and np.isfinite(red_share)) else 0.0
    return float(np.sqrt(zp ** 2 + cal ** 2 + red ** 2))


# --- --fit-alpha --------------------------------------------------------------

def alpha_grid_counts(filters, areas, z, law, ebv, rv, igm, host, host_ebv):
    """
    {filter: counts per unit normalisation for f_nu ~ nu^-alpha, one value
    per ALPHA_GRID entry}, the power law normalised to f_lambda = 1 at
    5000 A in the rest frame, through IGM, host and Galactic dust.
    """
    out = {}
    for f in filters:
        if f not in areas:
            continue
        lam, dl, area, _ = areas[f]
        t = igm(lam, z) if (igm and z > 0) else np.ones_like(lam)
        h = 10 ** (-0.4 * host(lam / (1 + z)) * host_ebv) \
            if (host is not None and host_ebv > 0) else np.ones_like(lam)
        ext = 10 ** (-0.4 * law(lam, rv) * ebv) if ebv > 0 \
            else np.ones_like(lam)
        base = t * h * ext * lam / (H_CGS * C_AA) * area * dl
        logl = np.log(lam / (1 + z) / 5000.0)
        out[f] = np.exp(np.outer(ALPHA_GRID - 2, logl)) @ base
    return out


def _chi2_curve(y, sig, model):
    """Chi^2 against ALPHA_GRID with the normalisation fitted analytically
    (model: filters x grid); returns chi2 and the normalisations."""
    w = 1.0 / sig[:, None] ** 2
    norm = np.sum(model * y[:, None] * w, 0) / np.sum(model ** 2 * w, 0)
    return np.sum((y[:, None] - norm * model) ** 2 * w, 0), norm


def _minimum(chi):
    """Best alpha (parabola through the grid minimum), the Delta chi^2 = 1
    half-width (crossings interpolated between grid points), and whether
    the minimum sits at the edge of the grid."""
    i = int(np.argmin(chi))
    step = ALPHA_GRID[1] - ALPHA_GRID[0]
    edge = i in (0, len(chi) - 1)
    best, cmin = ALPHA_GRID[i], chi[i]
    if not edge:
        c0, c1, c2 = chi[i - 1], chi[i], chi[i + 1]
        den = c0 - 2 * c1 + c2
        if den > 0:
            shift = 0.5 * (c0 - c2) / den
            best = best + shift * step
            cmin = c1 - 0.25 * (c0 - c2) * shift
    level = cmin + 1.0

    def crossing(direction):
        j = i
        while 0 <= j + direction < len(chi) and chi[j + direction] < level:
            j += direction
        k = j + direction
        if not 0 <= k < len(chi):
            return ALPHA_GRID[j]
        frac = (level - chi[j]) / (chi[k] - chi[j]) if chi[k] != chi[j] \
            else 0.0
        return ALPHA_GRID[j] + direction * frac * step
    return best, 0.5 * (crossing(1) - crossing(-1)), edge


def fit_alphas(rows, grids, grids_ebv_hi, eligible, z, min_filters,
               sys_of=None):
    """
    One power-law fit per observation: the detected rates of the eligible
    filters, each with rate_err and its flux_sys (sys_of(row), a fraction;
    default mag_sys alone) combined. Returns
    {obsid: table row (strings) plus '_alpha' (float or None if not usable)}.
    """
    if sys_of is None:
        def sys_of(r):
            return flux_sys(r['filter'], _num(r.get('mag_sys')), 0.0, {})
    by_obs = {}
    for r in rows:
        if r.get('status') != 'detected' or (r['obsid'], r['filter']) \
                not in eligible or r['filter'] not in grids:
            continue
        rate, err = _num(r.get('rate')), _num(r.get('rate_err'))
        if not (rate > 0 and err > 0):
            continue
        by_obs.setdefault(r['obsid'], []).append(r)
    results = {}
    for obsid, pts in by_obs.items():
        if len(pts) < min_filters:
            continue
        pts.sort(key=lambda r: FILTERS.index(r['filter']))
        y = np.array([_num(r['rate']) for r in pts])
        stat = np.array([_num(r['rate_err']) for r in pts])
        fsys = np.array([sys_of(r) for r in pts])
        sig = np.hypot(stat, y * fsys)
        model = np.array([grids[r['filter']] for r in pts])
        chi, norm = _chi2_curve(y, sig, model)
        alpha, err, edge = _minimum(chi)
        chi_stat, _ = _chi2_curve(y, stat, model)
        _, err_stat, _ = _minimum(chi_stat)
        chi_hi, _ = _chi2_curve(y, sig, np.array(
            [grids_ebv_hi[r['filter']] for r in pts]))
        alpha_hi, _, _ = _minimum(chi_hi)
        i = int(np.argmin(chi))
        dof = len(pts) - 2
        chi_min = float(chi[i])
        p = float(chi2_dist.sf(chi_min, dof)) if dof > 0 else float('nan')
        resid = (y - norm[i] * model[:, i]) / sig
        worst = int(np.argmax(np.abs(resid)))
        flags, reasons = [], []
        if edge:
            flags.append('alpha_at_edge')
            reasons.append('best alpha at the end of the range %g..%g; '
                           'fluxes use the default conversion'
                           % (ALPHA_GRID[0], ALPHA_GRID[-1]))
        if dof == 0:
            flags.append('two_filters')
            reasons.append('two filters: alpha is exact, the fit cannot be '
                           'checked')
        elif p < FIT_P_MIN:
            flags.append('poor_fit')
        if dof > 0:
            reasons.append('worst filter %s: %+.1f sigma (measured/model '
                           '%.3f)' % (pts[worst]['filter'], resid[worst],
                                      y[worst] / (norm[i] * model[worst, i])))
        # intrinsic nu F_nu of the fitted power law at 5000 A observed
        lam5 = 5000.0
        nufnu5 = norm[i] * (lam5 / (1 + z) / 5000.0) ** (ALPHA_GRID[i] - 2) \
            * lam5
        starts = [_num(r.get('mjd_start')) for r in pts]
        stops = [_num(r.get('mjd_stop')) for r in pts]
        span = (np.nanmax(stops) - np.nanmin(starts)) * 24.0 \
            if np.isfinite(starts).any() and np.isfinite(stops).any() \
            else float('nan')
        results[obsid] = {
            'obsid': obsid,
            'mjd_mid': '%.5f' % np.mean([_num(r['mjd_mid']) for r in pts]),
            'n_filt': str(len(pts)),
            'filters': ','.join(r['filter'] for r in pts),
            'alpha': '%.3f' % alpha, 'alpha_err': '%.3f' % err,
            'alpha_err_stat': '%.3f' % err_stat,
            'alpha_sys_ebv': '%.3f' % abs(alpha_hi - alpha),
            'chi2': '%.2f' % chi_min, 'dof': str(dof),
            'p_value': '%.3g' % p if np.isfinite(p) else '-',
            'dt_span_h': '%.2f' % span if np.isfinite(span) else '-',
            'nufnu_5000': '%.4e' % nufnu5,
            'flags': ','.join(flags) if flags else '-',
            'reason': '; '.join(reasons) if reasons else '-',
            '_alpha': None if edge else round(float(alpha), 2),
            '_flags': flags,
        }
    return results


def convert(row, term, fac, a_err, alpha='-', alpha_src='default',
            fit_flags=(), cal_sys=None):
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
               a_gal_err='%.3f' % a_err, a_host='%.3f' % term['a_host'],
               igm_trans='%.3f' % term['igm_trans'],
               red_share='%.2f' % term['red_share'], alpha=alpha,
               alpha_src=alpha_src, factor='%.4f' % term['factor'],
               flux_sys='%.3f' % flux_sys(row.get('filter'),
                                          _num(row.get('mag_sys')),
                                          term['red_share'],
                                          cal_sys or {}))
    if term['igm_trans'] < IGM_FLAG:
        flags.append('igm_absorbed')
    if term['red_share'] > RED_FLAG:
        flags.append('red_leak')
    flags.extend(f for f in fit_flags if f != 'alpha_at_edge')
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
    if 'poor_fit' in fit_flags:
        reasons.append('the power law fits this observation poorly (p < '
                       '%g; uvot_alpha.txt)' % FIT_P_MIN)
    if 'two_filters' in fit_flags:
        reasons.append('alpha from two filters only')
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


def plot(pdf_path, rows, label, alpha_rows=()):
    """Dereddened nu F_nu against time per filter, the SEDs, and the
    fitted alpha against time and brightness."""
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
        if alpha_rows:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 9))
            for sel, mfc in ((lambda a: not a['_flags'], 'k'),
                             (lambda a: bool(a['_flags']), 'none')):
                pts = [a for a in alpha_rows if sel(a)]
                if not pts:
                    continue
                t = [_num(a['mjd_mid']) for a in pts]
                al = [_num(a['alpha']) for a in pts]
                e = [_num(a['alpha_err_stat']) for a in pts]
                fl = [_num(a['nufnu_5000']) for a in pts]
                ax1.errorbar(t, al, yerr=e, fmt='o', ms=3, color='k', mfc=mfc,
                             lw=0.6)
                ax2.errorbar(fl, al, yerr=e, fmt='o', ms=3, color='k',
                             mfc=mfc, lw=0.6)
            ax1.set_xlabel('MJD')
            ax1.set_ylabel('alpha (f_nu ~ nu^-alpha)')
            ax2.set_xscale('log')
            ax2.set_xlabel('fitted nu F_nu at 5000 A observed (erg/s/cm^2), '
                           'dereddened')
            ax2.set_ylabel('alpha')
            ax1.set_title('%s: fitted alpha per observation (error without '
                          'flux_sys); open: flagged' % label, fontsize=9)
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
    parser.add_argument('--host-ebv', type=float, default=0.0,
                        help='E(B-V) of dust in the host galaxy, applied in '
                             'the rest frame (default: 0)')
    parser.add_argument('--host-law', choices=sorted(HOST_LAWS),
                        default='smc',
                        help='host dust law: Gordon et al. 2003 SMC bar (smc, '
                             'default) or LMC average (lmc), Fitzpatrick 1999 '
                             'R_V 3.1 (mw), Calzetti et al. 2000 (calzetti)')
    parser.add_argument('--fit-alpha', action='store_true',
                        help='fit f_nu ~ nu^-alpha per observation to the '
                             'detected filters and convert that '
                             'observation with it')
    parser.add_argument('--fit-filters', default=None,
                        help='comma-separated filters the fit may use '
                             '(default: all; rows flagged igm_absorbed or '
                             'red_leak are always left out)')
    parser.add_argument('--fit-min-filters', type=int,
                        default=FIT_MIN_FILTERS,
                        help='fewest detected filters for a fit (default: '
                             '3; 2 gives an exact, unchecked alpha, flagged '
                             'two_filters)')
    parser.add_argument('--cal-sys', default=None, metavar='VALUES',
                        help='calibration systematic of the conversion, a '
                             'fraction: one number for all filters, or e.g. '
                             'V=0.03,UVW2=0.08; 0 leaves only mag_sys '
                             '(default: %s, plus %s x red_share)'
                             % (','.join('%s=%g' % kv for kv in
                                         CAL_SYS.items()),
                                ','.join('%s %g' % kv for kv in
                                         RED_TAIL_SYS.items())))
    parser.add_argument('--irsa-timeout', type=float, default=60,
                        help='seconds to wait for IRSA (default: 60)')
    args = parser.parse_args(argv)
    if args.host_ebv < 0:
        print('ERROR: --host-ebv must not be negative', file=sys.stderr)
        return 1
    if args.fit_min_filters < 2:
        print('ERROR: --fit-min-filters must be 2 or more', file=sys.stderr)
        return 1
    fit_filters = FILTERS
    if args.fit_filters:
        fit_filters = [f.strip().upper() for f in args.fit_filters.split(',')
                       if f.strip()]
        bad = [f for f in fit_filters if f not in FILTERS]
        if bad or len(fit_filters) < args.fit_min_filters:
            print('ERROR: --fit-filters: give %d or more of %s'
                  % (args.fit_min_filters, ', '.join(FILTERS)),
                  file=sys.stderr)
            return 1
    cal_sys = dict(CAL_SYS)
    if args.cal_sys:
        try:
            if '=' in args.cal_sys:
                for item in args.cal_sys.split(','):
                    f, v = item.split('=')
                    if f.strip().upper() not in FILTERS:
                        raise ValueError(f)
                    cal_sys[f.strip().upper()] = float(v)
            else:
                cal_sys = {f: float(args.cal_sys) for f in CAL_SYS}
        except ValueError:
            print('ERROR: --cal-sys %s: give a number or FILTER=VALUE pairs '
                  '(filters %s)' % (args.cal_sys, ', '.join(FILTERS)),
                  file=sys.stderr)
            return 1
        if any(v < 0 for v in cal_sys.values()):
            print('ERROR: --cal-sys must not be negative', file=sys.stderr)
            return 1
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
    host = HOST_LAWS[args.host_law] if args.host_ebv > 0 else None
    filters = sorted({r['filter'] for r in rows},
                     key=lambda f: FILTERS.index(f) if f in FILTERS else 99)
    terms = filter_terms(filters, areas, factors, shape, args.z, law, ebv,
                         args.rv, igm, bool(args.spectrum), host,
                         args.host_ebv)

    fits_by_obs = {}
    if args.fit_alpha:
        # eligible: filters asked for, and not flagged with the default terms
        eligible = {(r['obsid'], r['filter']) for r in rows
                    if r['filter'] in fit_filters and r['filter'] in terms
                    and terms[r['filter']]['igm_trans'] >= IGM_FLAG
                    and terms[r['filter']]['red_share'] <= RED_FLAG}
        grids = alpha_grid_counts(filters, areas, args.z, law, ebv, args.rv,
                                  igm, host, args.host_ebv)
        grids_hi = alpha_grid_counts(filters, areas, args.z, law,
                                     ebv + ebv_err, args.rv, igm, host,
                                     args.host_ebv)
        def sys_of(r):
            return flux_sys(r['filter'], _num(r.get('mag_sys')),
                            terms[r['filter']]['red_share'], cal_sys)
        fits_by_obs = fit_alphas(rows, grids, grids_hi, eligible, args.z,
                                 args.fit_min_filters, sys_of)
    alpha_terms = {}
    given_alpha, given_src = '-', 'uvotsource'
    if args.spectrum:
        given_src = 'spectrum'
        m = re.match(r'\s*powerlaw:\s*([-+\d.eE]+)', args.spectrum)
        if m:
            given_alpha = '%g' % float(m.group(1))

    out = []
    for r in rows:
        f = r['filter']
        term = terms.get(f)
        alpha, src, fflags = given_alpha, given_src, ()
        fit = fits_by_obs.get(r['obsid'])
        if fit is not None and term is not None:
            fflags = fit['_flags']
            if fit['_alpha'] is not None:
                a = fit['_alpha']
                if a not in alpha_terms:
                    sh, _ = parse_spectrum('powerlaw:%g' % a)
                    alpha_terms[a] = filter_terms(
                        filters, areas, factors, sh, args.z, law, ebv,
                        args.rv, igm, True, host, args.host_ebv)
                term = alpha_terms[a].get(f, term)
                alpha, src = fit['alpha'], 'fit'
        a_err = term['a_gal'] / ebv * ebv_err if (term and ebv > 0) else 0.0
        out.append(convert(r, term, factors.get(f), a_err, alpha, src,
                           fflags, cal_sys))

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
        'host dust: %s' % ('none' if host is None else
                           'E(B-V) = %.4f mag, law %s, rest frame (z = %g)'
                           % (args.host_ebv, args.host_law, args.z)),
        'alpha (f_nu ~ nu^-alpha): %s' % (
            'fitted per observation (--fit-alpha, %d+ detected filters of '
            '%s; uvot_alpha.txt); alpha_src fit, else the default above'
            % (args.fit_min_filters, ','.join(fit_filters))
            if args.fit_alpha else 'not fitted'),
        'lam_ref (A): wavelength of uvotsource\'s factors; flux_aa '
        '(erg/s/cm^2/A), flux_mjy (mJy), nufnu (erg/s/cm^2): intrinsic '
        'spectrum at lam_ref, without Galactic dust, host dust and IGM',
        'factor: flux_aa / Step 7 flux_aa; a_gal (mag): folded Galactic '
        'extinction, a_gal_err its systematic error (not in the flux '
        'errors); a_host (mag): folded host-galaxy extinction; igm_trans: '
        'folded IGM transmission of the counts',
        'red_share: share of the model counts from beyond the filter\'s main '
        'band (effective area >= %g %% of peak); flags igm_absorbed '
        '(igm_trans < %g), red_leak (red_share > %g)'
        % (100 * MAIN_BAND, IGM_FLAG, RED_FLAG),
        'flux_sys: relative calibration systematic of the flux columns (not '
        'in the errors): Step 7 mag_sys, the conversion (--cal-sys: %s) and '
        '%s x red_share in quadrature'
        % (','.join('%s=%g' % kv for kv in cal_sys.items()),
           ','.join('%s %g' % kv for kv in RED_TAIL_SYS.items())
           if any(cal_sys.values()) else 'no red-tail term'),
        'mag_dered, ab_mag_dered: Step 7 magnitudes minus the same correction; '
        'upper_limit rows: flux and magnitude columns hold the limit',
    ]
    if args.fit_alpha:
        comments.append('per-filter values below: for the default spectrum; '
                        'rows with alpha_src fit have their own')
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
    alpha_rows = sorted(fits_by_obs.values(), key=lambda a: _num(a['mjd_mid']))
    alpha_path = None
    if args.fit_alpha:
        alpha_path = os.path.join(args.outdir, 'uvot_alpha.txt')
        write_table(alpha_path, alpha_rows, ALPHA_COLUMNS, comments[:10] + [
            'one row per observation with %d+ detected filters among %s, '
            'none flagged igm_absorbed or red_leak'
            % (args.fit_min_filters, ','.join(fit_filters)),
            'model: f_nu ~ nu^-alpha (rest frame) through IGM, host and '
            'Galactic dust, folded through each filter; normalisation fitted',
            'errors per filter: rate_err and flux_sys combined; alpha_err: '
            'Delta chi^2 = 1; alpha_err_stat: the same with rate_err only '
            '(for changes of alpha between observations)',
            'alpha_sys_ebv: change of alpha when Galactic E(B-V) moves up by '
            'its error; dt_span_h: first start to last stop of the filters '
            'used (hours)',
            'nufnu_5000: the fitted intrinsic nu F_nu at 5000 A observed '
            '(erg/s/cm^2)',
            'flags: poor_fit (p < %g), two_filters, alpha_at_edge (alpha '
            'outside %g..%g; not used for the fluxes)'
            % (FIT_P_MIN, ALPHA_GRID[0], ALPHA_GRID[-1])])
    pdf_path = os.path.join(args.outdir, 'uvot_flux.pdf')
    plot(pdf_path, out, label, alpha_rows)

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
    if args.fit_alpha:
        good = [a for a in alpha_rows if a['_alpha'] is not None]
        av = np.array([_num(a['alpha']) for a in good])
        npoor = sum('poor_fit' in a['_flags'] for a in alpha_rows)
        if len(av):
            print('[summary] alpha fits: %d observations (%d used for the '
                  'fluxes); median alpha %.2f (16-84 %%: %.2f..%.2f); '
                  'poor_fit %d' % (len(alpha_rows), len(good), np.median(av),
                                   *np.percentile(av, [16, 84]), npoor))
        else:
            print('[summary] alpha fits: none (no observation with %d+ '
                  'eligible detected filters)' % args.fit_min_filters)
    print('[summary] Written: %s, %s%s' % (
        flux_path, pdf_path, ', ' + alpha_path if alpha_path else ''))
    return 0


if __name__ == '__main__':
    sys.exit(main())
