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

**Flags** (flag-only, nothing is dropped):

- `igm_absorbed`: the IGM keeps less than 90 % of the counts. Both IGM
  models are *mean* transmissions, and single sightlines scatter widely
  around them where the transmission is low, so these values depend on the
  model.
- `red_leak`: more than 30 % of the model's counts come from beyond the long
  end of the filter's main band (where the effective area falls below 10 %
  of its peak). The measurement then mostly tells about longer wavelengths.

## Validation

| Check | Result |
| ----- | ------ |
| Folding the CALDB effective areas vs `uvotsource`'s factors, f_ν ∝ ν^−α, α = 0–1 | within 0.1–1.1 % in V, B, U, UVM2, UVW2; 2.8 % in UVW1 |
| Dust laws vs the `extinction` package (Barbary), 1600–8000 Å, R_V 2.5–4 | identical (≤ 10⁻¹⁵) |
| Inoue 2014 Lyman-α forest optical depth vs observed (Becker et al. 2013) | z = 3: 0.400 vs 0.402; z = 4: 0.911 vs 0.887 |
| IRSA E(B−V) for PKS 0537−286 | S&F11 0.0208, SFD 0.0242 |

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

## Outputs

`UVOT_output/uvot_flux.txt`, one row per light-curve row:

| Column | Meaning |
| ------ | ------- |
| `obsid`, `filter`, `date_mid`, `mjd_mid`, `status`, `rate`, `rate_err`, `rate_limit` | as in Step 7 |
| `lam_ref` | reference wavelength (Å) |
| `a_gal`, `a_gal_err` | folded Galactic extinction (mag) and its systematic error |
| `igm_trans` | folded IGM transmission of the counts (1 without `--z`) |
| `red_share` | share of the model's counts from beyond the filter's main band |
| `factor` | `flux_aa` / Step 7's `flux_aa` |
| `flux_aa`, `flux_aa_err` | intrinsic f_λ at `lam_ref` (erg/s/cm²/Å) |
| `flux_mjy`, `flux_mjy_err` | the same as f_ν (mJy) |
| `nufnu`, `nufnu_err` | ν F_ν (erg/s/cm²), for spectral energy distributions |
| `mag_dered`, `ab_mag_dered`, `mag_err` | Step 7's Vega and AB magnitudes minus the same correction |
| `flags`, `reason` | `igm_absorbed`, `red_leak`, with the numbers |

Upper-limit rows hold the limit in the flux and magnitude columns, without
errors. The header lists, per filter, the effective-area file, A_gal and
A/E(B−V), the IGM transmission, the red-leak share and the factor.

`UVOT_output/uvot_flux.pdf`:

- dereddened νF_ν against time, one panel per filter (open symbols:
  flagged; triangles: upper limits);
- the spectral energy distribution of every observation with two or more
  detected filters, coloured by time.

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

4. **Dust in the host galaxy** is not corrected. Only Galactic dust and the
   IGM are.

5. **Re-run after Step 7.** Step 8 reads `uvot_lightcurve.txt`.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_
