# Step 8: Flux densities, dereddened

`swift_uvot_flux.py` turns the count rates of the [Step 7](07-lightcurve.md)
light curve into flux densities corrected for Galactic dust and, for a
distant source, for the hydrogen between the galaxies (the intergalactic
medium, IGM).

It changes nothing in Steps 1–7. The count rate stays the measured quantity.
Everything here depends on assumptions (the spectrum, the dust law, the IGM
model), and all of them are written into the table header. Step 8 reads only
tables and the CALDB, needs no HEASoft, and takes a few seconds.

## What runs

```bash
swift_uvot_flux.py                                   # dereddened, uvotsource's factors
swift_uvot_flux.py --z 3.104 --spectrum powerlaw:1.5 # a distant source: IGM too
swift_uvot_flux.py --z 0.158 --fit-alpha             # fit the spectral slope per observation
```

| Flag | Purpose |
| ---- | ------- |
| `--outdir` | Steps 3–7's output folder (default `UVOT_output`) |
| `--ebv` | Galactic E(B−V) in mag. Default: from the IRSA dust service. `0` skips dereddening |
| `--ebv-err` | Its error (default 16 % of E(B−V)) |
| `--ebv-map` | IRSA map: `sf11` (Schlafly & Finkbeiner 2011, default) or `sfd` (Schlegel et al. 1998) |
| `--law` | Dust law: `f99` (Fitzpatrick 1999, default) or `ccm89` (Cardelli, Clayton & Mathis 1989) |
| `--rv` | R_V, the dust's total-to-selective extinction ratio (default 3.1) |
| `--z` | Redshift of the source, for the IGM (default 0: no IGM) |
| `--igm` | IGM model when `--z` is given: `inoue14` (Inoue et al. 2014, default), `madau95` (Madau 1995) or `none` |
| `--spectrum` | `powerlaw:ALPHA` (f_ν ∝ ν^−α), `blackbody:T` (kelvin, rest frame) or `file:PATH` (rest-frame wavelength in Å, f_λ in any units). Default: keep `uvotsource`'s factors |
| `--host-ebv` | E(B−V) of dust in the host galaxy, applied in the rest frame (default 0) |
| `--host-law` | Host dust law: `smc` (Gordon et al. 2003 SMC bar, default), `lmc` (Gordon et al. 2003 LMC average), `mw` (Fitzpatrick 1999, R_V 3.1) or `calzetti` (Calzetti et al. 2000) |
| `--fit-alpha` | Fit a power law f_ν ∝ ν^−α to each observation's filters and convert that observation with it |
| `--fit-filters` | Filters the fit may use, e.g. `V,B,UVW1,UVM2,UVW2` (default: all) |
| `--fit-min-filters` | Fewest detected filters for a fit (default 3; 2 is allowed and flagged) |
| `--cal-sys` | Calibration systematic of the conversion, a fraction: one number, or e.g. `V=0.03,UVW2=0.08`; `0` leaves only `mag_sys` (default: 0.05 in V, B, U; 0.04 in UVW1, UVM2, UVW2; plus the red-tail term) |
| `--irsa-timeout` | Seconds to wait for IRSA (default 60) |

## How it works

**The problem.** A UVOT filter counts photons over a wide band of
wavelengths. Step 7 turns its count rate into a flux density with one fixed
factor per filter, `uvotsource`'s, which assumes a smooth, GRB-like spectrum.
That goes wrong in three ways:

1. **Another spectral shape needs another factor.** The difference is a few
   per cent for smooth spectra and more for steep ones.
2. **Dust dims the source, and the dimming depends on the spectrum.** Dust
   dims the blue end of a broad filter more than the red end, so a filter's
   extinction depends on where its photons come from. The effect is largest
   in the UV filters, which straddle the dust's strong absorption at 2175 Å.
3. **Light is removed inside the band.** For a distant source, hydrogen along
   the line of sight absorbs light bluer than its Lyman-α line (1216 Å in the
   rest frame, the "Lyman-α forest"). It absorbs everything bluer than
   912 Å (the Lyman limit). Above redshift z ≈ 0.3 this happens inside the
   UV filters; above z ≈ 2 also inside U and B. A UV filter may then count
   mostly the far red tail of its response, the **red leak**.

**The fix: folding.** Step 8 builds a model spectrum: the source's own shape,
times the IGM transmission, times the Galactic extinction. It then computes
the count rate that model would give through the filter's effective area
(the newest `SPECRESP` file of each filter in the UVOT CALDB). This is called
"folding". The measured rate fixes the model's normalisation. The flux
density is the *intrinsic* spectrum (without dust and IGM) at the filter's
reference wavelength `lam_ref`. That is the wavelength at which
`uvotsource`'s two factors (erg/s/cm²/Å and mJy) agree, so the numbers
compare directly with Step 7's.

- **Default** (no `--spectrum`): `uvotsource`'s factors are kept and only
  the corrections are applied. The Galactic extinction and, with `--z`, the
  IGM transmission are folded for f_ν ∝ ν^−1.
- **`--spectrum`**: the factor itself comes from folding that spectrum.

**E(B−V)** comes from IRSA's dust service, for the map pixel at the source.
It is cached in `UVOT_output/uvot_ebv.txt`, so IRSA is asked once per
analysis. Its error is 16 %, the SFD map's own, kept as a separate column
(`a_gal_err`) and not added to the flux errors. Like `mag_sys`, it moves
every point of a filter together.

**Calibration systematic** (`flux_sys`, a fraction of the flux, not added
to the flux errors). The check against stars of known spectra (Validation,
below) shows that the UVOT calibration itself is good to about 5 % for one
source, depending on its colour. `flux_sys` combines in quadrature:

- Step 7's `mag_sys`: `uvotsource`'s zero-point error (1–3 %);
- the conversion's own error, `--cal-sys`: 5 % in V, B and U, 4 % in
  UVW1, UVM2 and UVW2 (default);
- a **red-tail term** in UVW1 and UVW2: 0.25 and 0.6 × `red_share`. The
  CALDB correction for the slow loss of sensitivity (SENSCORR) is one
  factor per filter, measured on hot stars. The red tail of the UV filters
  seems to have lost little sensitivity, so the correction over-corrects
  the counts from there, by more every year (Validation). For a power law
  with α = 1.5 the term is about 4 % in UVW2; for a Sun-like star 40 %.

`--cal-sys 0` leaves `mag_sys` alone (the behaviour before this option).

**Flags** (flag-only, nothing is dropped):

- `igm_absorbed`: the IGM keeps less than 90 % of the counts. Both IGM
  models are *mean* transmissions, and single sightlines scatter widely
  around them where the transmission is low, so these values depend on the
  model.
- `red_leak`: more than 30 % of the model's counts come from beyond the long
  end of the filter's main band (where the effective area falls below 10 %
  of its peak). The measurement then mostly tells about longer wavelengths.

**Host-galaxy dust** (`--host-ebv`). Dust in the source's own galaxy dims
it at the wavelength the light had when it left, λ/(1+z). The host
extinction is folded in the same way as the Galactic one and written as
`a_host`. Laws:

- `smc`, the SMC bar curve of Gordon et al. (2003), is the default. Quasar
  reddening looks SMC-like (Hopkins et al. 2004), and it has no 2175 Å bump.
- `lmc` is the LMC average of Gordon et al. (2003).
- `mw` is the Milky Way curve (Fitzpatrick 1999, R_V 3.1).
- `calzetti` is the attenuation curve of star-forming galaxies (Calzetti
  et al. 2000).

The curves end at 1000 Å (Calzetti: 1200 Å) and are held constant below,
where the IGM removes most of the light anyway. The SMC and LMC curves
match the `dust_extinction` package to 10⁻¹⁴.

**Fitting the spectral slope** (`--fit-alpha`). The slope α of
f_ν ∝ ν^−α is a quantity of its own (for a blazar, "bluer when
brighter"), and it sets the folding. With `--fit-alpha`, each observation
(OBSID) with three or more detected filters gets its own power law:

1. The model is the power law, times the IGM transmission, the host and
   the Galactic dust, folded through each filter. Its normalisation and α
   are fitted to the measured count rates.
2. Each filter's error is `rate_err` and its calibration systematic
   `flux_sys` combined.
3. That observation's rows (detections and limits) are then converted
   with the fitted α. Observations without a fit keep the default
   conversion.

Rows flagged `igm_absorbed` or `red_leak` are never fitted; `--fit-filters`
leaves out more. The fits are listed in `uvot_alpha.txt`:

- **Two errors.** `alpha_err` contains `flux_sys`: use it for the absolute
  α. But `flux_sys` is (nearly) the same at every epoch, so it does not blur changes
  of α from one observation to the next: for those use `alpha_err_stat`.
- **Goodness of fit.** `poor_fit` marks fits with p < 0.001; `reason`
  names the filter that deviates most. A power law cannot follow emission
  lines or a break.
- **Simultaneity.** `dt_span_h` is the time from the first filter's start
  to the last filter's end. A fast-varying source can change in between.

How much the fluxes change: for an unreddened power law, the folded factor
differs from `uvotsource`'s by ≤ 3 % for α between −0.5 and 2, and by 8 %
(UVW2) at α = 3. For most sources the gain is α itself, and dust and IGM
corrections that match the spectrum.

## Validation

| Check | Result |
| ----- | ------ |
| Folding the CALDB effective areas vs `uvotsource`'s factors, f_ν ∝ ν^−α, α = 0–1 | within 0.1–1.1 % in V, B, U, UVM2, UVW2; 2.8 % in UVW1 |
| Dust laws vs the `extinction` package (Barbary), 1600–8000 Å, R_V 2.5–4 | identical (≤ 10⁻¹⁵) |
| Inoue 2014 Lyman-α forest optical depth vs observed (Becker et al. 2013) | z = 3: 0.400 vs 0.402; z = 4: 0.911 vs 0.887 |
| IRSA E(B−V) for PKS 0537−286 | S&F11 0.0208, SFD 0.0242 |
| Host-dust laws (SMC, LMC: Gordon et al. 2003) vs the `dust_extinction` package, 1200–12000 Å | identical (≤ 10⁻⁸); Calzetti k(1500 Å) = 10.33, k(5500 Å) = 4.05 |
| `--fit-alpha` on synthetic count rates (α = −0.3, 0.8, 1.7; with Galactic and SMC host dust at z = 0.5) | input α recovered; pull spread 0.8–1.2; a 20 % excess in U gives `poor_fit` naming U |
| Folding vs measured count rates of stars with known spectra | see **Stars with known spectra** below: within 3 % for some filters and colours, up to 6 % (V, B) and 23 % (UVW2, red stars) for others |

**Folded Galactic extinction**, A_filter / E(B−V) for f_ν ∝ ν^−1
(Fitzpatrick 1999, R_V = 3.1), compared with the curve's value at the single
wavelength `lam_ref`:

| | V | B | U | UVW1 | UVM2 | UVW2 |
| - | - | - | - | ---- | ---- | ---- |
| folded | 3.10 | 4.09 | 4.94 | 6.87 | 8.44 | 8.07 |
| at `lam_ref` only | 3.11 | 4.12 | 4.89 | 6.63 | 9.30 | 8.82 |

The difference reaches 0.8 × E(B−V) in UVM2 and UVW2: 0.02 mag for
PKS 0537−286 (E(B−V) = 0.021), but 0.25 mag in a field with E(B−V) = 0.3.

**PKS 0537−286** (z = 3.104, f_ν ∝ ν^−1.5), IGM transmission of the counts
and the resulting factor relative to Step 7:

| IGM model | V | B | U | UVW1 | UVM2 | UVW2 |
| --------- | - | - | - | ---- | ---- | ---- |
| Inoue 2014 | 1.00 (×1.06) | 0.73 (×1.48) | 0.36 (×3.1) | 0.13 (×8.6) | 0.09 (×13) | 0.10 (×11) |
| Madau 1995 | 1.00 (×1.06) | 0.68 (×1.59) | 0.24 (×4.6) | 0.05 (×22) | 0.02 (×60) | 0.03 (×38) |

V is unaffected. In the UV the two models differ by 2–3 times, so these
values are bounds at best, and the rows are flagged. With Madau 1995, 60 %
of the UVW1 counts come from its red tail (flag `red_leak`). In the 2005
observations (00035240001/2), νF_ν is:

- **Without the IGM correction:** U is 7 times below V.
- **With Inoue 2014:** B sits 4–17 % above V, and U is still 2.5–2.7
  times below it.

This is the expected picture when the U band lies mostly beyond the Lyman
limit at this redshift: a mean model cannot predict one sightline's Lyman
continuum.

**3C 273** (z = 0.158, E(B−V) = 0.018): no IGM effect in any filter
(Lyman-α is at 1408 Å). The dereddened factors are 1.05 (V) to 1.15 (UVM2).
The spectral energy distributions (second page of the PDF) show the rise
into the UV (the "big blue bump"). U sits a median 19 % above the line from
B to UVW1 (92 observations); Mg II (2798 Å, observed at 3240 Å) and the Balmer continuum fall in U, a
case where a line-free power law is not the right spectrum.
With `--fit-alpha`, the median α is 0.22 (16–84 %: 0.14–0.30). A power law
fits all filters except U: with U the median χ²/dof is 7.7 and almost every
`poor_fit` names U as the worst filter; without U (`--fit-filters
V,B,UVW1,UVM2,UVW2`) it is 0.2.

### Stars with known spectra

The checks above compare Step 8 with `uvotsource`, which uses the same
calibration. This one tests the whole chain (Steps 4–7 photometry,
coincidence loss, aperture correction, effective areas, folding): for stars
whose spectrum is known, fold the spectrum and compare the predicted count
rate with the rate the pipeline measured. Ratio = median measured /
predicted, over all detected exposures without `senscorr_frozen` or
saturation flags.

- **Field stars:** 9 isolated Gaia stars in the 3C 273 field (G 13.4–15.2,
  BP−RP 0.75–1.61, i.e. G and K stars), with their Gaia XP spectra
  (336–1020 nm: V and B only; U is only 60 % covered).
- **HST CALSPEC standards:** three hot white dwarfs (WD 1657+343, GD 50,
  GD 153) and two Sun-like G stars (P041C, P177D). The G stars' STIS spectra
  start at 2222 Å; below that the CALSPEC model is used, so UVM2 (24 % of
  the counts from the model) is not a test for them.

| Filter | Field stars (9) | White dwarfs | G-star standards |
| ------ | --------------- | ------------ | ---------------- |
| V | 1.005 (0.990–1.025) | 0.945–0.965 | – |
| B | 1.043 (1.016–1.064) | 0.934–0.978 | – |
| U | – | 0.990, 1.054 | 0.953, 0.980 |
| UVW1 | – | 0.960–0.969 | 1.000, 1.049 |
| UVM2 | – | 0.978, 0.984 | (1.04, model) |
| UVW2 | – | 0.973 | 1.16, 1.23 |

Typical statistical error of a median: 0.2–1 %. The pattern is the same in
V, B, UVW1 and UVW2: **blue stars measure a few % low, red stars a few %
high**. That points to the shape of the CALDB effective-area curves (a
colour term), not to the folding:

- Folding the CALSPEC spectrum of Vega reproduces the UVOT zero points to
  0.3–3.2 % in all six filters.
- It is not coincidence loss: WD 1657+343 is faint in V (3.7 counts/s) and
  still 3.5 % low. Among the field stars, the brightest (100 counts/s in B)
  has the lowest ratio in V and B, a smaller effect of 1–3 %.
- It is not the Gaia XP calibration: XP photometry standardised to Landolt
  stars differs from the raw XP spectra by ≤ 1.2 % in B and V, mostly in the
  direction that would widen the gap.
- The white dwarf's ratios are steady over time (WD 1657+343, 2005–2026).
  **The G stars' UV ratios are not:** in UVW2, P041C reads 0.94–0.98 in
  2009, 1.08–1.10 in 2012 and 1.19–1.26 from 2017 on; P177D 1.01 in 2010
  and 1.23–1.33 in 2019–20 (UVW1: 0.95–1.00 to 1.03–1.08; UVM2: no
  trend). At the time of the UV calibration (Breeveld et al. 2011, which
  used these stars) they agreed. About 75 % of a G star's UVW2 counts come
  from the filter's red tail (> 2600 Å; a white dwarf's 3 %). A model in
  which that tail keeps its sensitivity while SENSCORR corrects all counts
  by the white dwarfs' loss fits the trend (χ² 39 for 52 degrees of
  freedom, against 151 for a constant per star). This was checked against
  coincidence loss, background, detector position, exposure times and the
  stars' own variability; none explains it. Not followed up (the CALDB
  note lists these two stars among its monitoring stars, yet shows no such
  drift); the red-tail term of `flux_sys` covers it.

**What this means:** absolute flux densities are good to about ±5 % in V
and B, depending on the source's colour, and to 3–5 % in U, UVW1 and UVM2;
`flux_sys` carries this. In UVW2 a red source may read 20–30 % high after
2015 (flag `red_leak` when most counts come from the red tail). For a
blazar-like power law the UVW2 excess is 2–4 % after 2015 and grows slowly
with time, so it can mimic a slow brightening of that size. For
`--fit-alpha` the offsets can shift the absolute α; changes of α between
observations are affected only through the slow UVW2 drift.

## Outputs

`UVOT_output/uvot_flux.txt`, one row per light-curve row:

| Column | Meaning |
| ------ | ------- |
| `obsid`, `filter`, `date_mid`, `mjd_mid`, `status`, `rate`, `rate_err`, `rate_limit` | as in Step 7 |
| `lam_ref` | reference wavelength (Å) |
| `a_gal`, `a_gal_err` | folded Galactic extinction (mag) and its systematic error |
| `a_host` | folded host-galaxy extinction (mag; 0 without `--host-ebv`) |
| `igm_trans` | folded IGM transmission of the counts (1 without `--z`) |
| `red_share` | share of the model's counts from beyond the filter's main band |
| `alpha`, `alpha_src` | the slope used: `fit` (`--fit-alpha`), `spectrum` (`--spectrum`; α if a power law) or `uvotsource` (default factors) |
| `factor` | `flux_aa` / Step 7's `flux_aa` |
| `flux_aa`, `flux_aa_err` | intrinsic f_λ at `lam_ref` (erg/s/cm²/Å) |
| `flux_mjy`, `flux_mjy_err` | the same as f_ν (mJy) |
| `nufnu`, `nufnu_err` | ν F_ν (erg/s/cm²), for spectral energy distributions |
| `flux_sys` | calibration systematic of the flux columns, a fraction (not in the errors): `mag_sys`, `--cal-sys` and the red-tail term in quadrature |
| `mag_dered`, `ab_mag_dered`, `mag_err` | Step 7's Vega and AB magnitudes minus the same correction |
| `flags`, `reason` | `igm_absorbed`, `red_leak`, and from the fit `poor_fit`, `two_filters`, with the numbers |

Upper-limit rows hold the limit in the flux and magnitude columns, without
errors. The header lists, per filter, the effective-area file, A_gal and
A/E(B−V), the IGM transmission, the red-leak share and the factor.

`UVOT_output/uvot_flux.pdf`:

- dereddened νF_ν against time, one panel per filter (open symbols:
  flagged; triangles: upper limits);
- the spectral energy distribution of every observation with two or more
  detected filters, coloured by time;
- with `--fit-alpha`: α against time, and α against the fitted νF_ν at
  5000 Å (open symbols: flagged; error bars without `flux_sys`).

`UVOT_output/uvot_alpha.txt` (with `--fit-alpha`), one row per fitted
observation:

| Column | Meaning |
| ------ | ------- |
| `obsid`, `mjd_mid` | observation, mean mid-time of the filters used |
| `n_filt`, `filters` | the filters fitted |
| `alpha`, `alpha_err`, `alpha_err_stat` | the slope; its error (Δχ² = 1) with and without `flux_sys` |
| `alpha_sys_ebv` | the change of α when the Galactic E(B−V) moves up by its error |
| `chi2`, `dof`, `p_value` | goodness of fit |
| `dt_span_h` | hours from the first filter's start to the last filter's end |
| `nufnu_5000` | the fitted intrinsic νF_ν at 5000 Å observed (erg/s/cm²) |
| `flags`, `reason` | `poor_fit` (p < 0.001), `two_filters`, `alpha_at_edge` (α outside −3…6, not used for the fluxes); the worst filter's deviation |

`UVOT_output/uvot_ebv.txt`: the cached E(B−V) values of both maps.

## Examples

```bash
# Default: E(B-V) from IRSA, uvotsource's factors, dereddened
swift_uvot_flux.py

# A given E(B-V) (offline), Cardelli et al. law, R_V 3.1
swift_uvot_flux.py --ebv 0.05 --law ccm89

# A z = 3.1 quasar, with the older IGM model for comparison
swift_uvot_flux.py --z 3.104 --spectrum powerlaw:1.5
swift_uvot_flux.py --z 3.104 --spectrum powerlaw:1.5 --igm madau95

# Fit alpha per observation; leave out U, where 3C 273's emission lines fall
swift_uvot_flux.py --z 0.158 --fit-alpha --fit-filters V,B,UVW1,UVM2,UVW2

# A reddened quasar: SMC-like dust in the host, E(B-V) = 0.1
swift_uvot_flux.py --z 0.5 --host-ebv 0.1

# Your own spectrum (rest frame, e.g. a quasar template with emission lines)
swift_uvot_flux.py --z 0.158 --spectrum file:template.txt
```

## Gotchas

1. **Count rates are the measurement.** Every flux here depends on the
   assumed spectrum. Quote the assumptions (they are in the header) with
   the numbers.

2. **The IGM models are averages.** Single sightlines vary, most where the
   transmission is low (U and the UV filters above z ≈ 2.5). Treat
   `igm_absorbed` rows as bounds or model-dependent values.

3. **Emission lines and breaks.** A power law has none. For a source with a
   strong line or a break in a band, give its spectrum with `--spectrum
   file:`. Rest-frame wavelengths: the file is shifted by `--z`.

4. **Dust in the host galaxy** is corrected only when you give
   `--host-ebv`. Its value has to come from elsewhere (e.g. the Balmer
   decrement, or an X-ray column). Step 8 does not fit it: with three to
   six broad filters it is degenerate with α.

5. **Colour-dependent calibration.** Measured rates of stars with known
   spectra differ from the folded prediction by up to ~6 % in V and B
   (blue stars low, red stars high) and by 20–30 % in UVW2 for red stars
   after 2015 (see Validation). `flux_sys` carries this; add it in
   quadrature when comparing with other instruments. Trust differences
   between epochs more than absolute fluxes or the absolute α, except for
   slow trends of a few % in UVW2.

6. **`poor_fit` is rarer with `flux_sys`.** A 5 % systematic per filter
   absorbs moderate deviations: on 3C 273, U (emission lines, ~15 % above
   the power law) is still the worst filter in most fits, but only 3 of 239
   fits are `poor_fit` (with `--cal-sys 0`, over 100). Read `reason`, or
   leave the filter out with `--fit-filters`.

7. **Re-run after Step 7.** Step 8 reads `uvot_lightcurve.txt`.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_
