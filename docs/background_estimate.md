# Analytic estimate of the GEM background rates (SoLID DDVCS configuration)

Script: `scripts/estimate_background_rates.py`. Output table: `config/gem_background_rates_analytic.txt`,
usable directly with `-Pgem:background_file=config/gem_background_rates_analytic.txt`.

This replaces the hand-set placeholder in `BackgroundModel` with a first-principles **direct-source
floor**: every particle that can travel from the target region straight to a GEM plane, computed from
the cross sections in the form used by the Geant4 Physics Reference Manual (PRM) and the PDG, folded
with the actual beam current, target, solenoid field at the target and GEM geometry. It deliberately
does **not** contain the albedo of the downstream beam pipe, collimators and endcap; that term needs
a Geant4 beam-on-target run and is discussed at the end.

## Beam and target

| Quantity | Value |
|---|---|
| Beam energy | 11 GeV |
| Target | LH2, 15 cm, ρ = 0.0708 g/cm³, X₀ = 890 cm, Al cell wall 0.18 mm |
| Areal density | 6.35×10²³ electrons (= protons) / cm² |
| Current for L = 1.2×10³⁷ | 3.0 µA (matches the J/ψ configuration) |
| Current for L = 10³⁸ / 10³⁹ with 15 cm | 25 µA / 253 µA |
| Alternative for 10³⁹ | 40 cm LH2 at 95 µA (PVDIS-like) |
| Target thickness | 0.017 X₀ (15 cm), 0.045 X₀ (40 cm) |
| B on axis at the target (z = −315 cm) | 0.28 T (model, no yoke); 1.5 T at the coil centre |

## Sources and formulas

1. **Møller scattering** e⁻e⁻ → e⁻e⁻ on target electrons (PRM "Møller–Bhabha ionisation"):
   dσ/dT with the full Møller expression; recoil kinetic energy vs lab angle T(θ) = 2mₑp₀²cos²θ /
   ((E₀+mₑ)² − p₀²cos²θ), i.e. T ≈ 2mₑcot²θ. Integrated above T = 1 MeV: **3.1×10¹² electrons/s at 3 µA**,
   uniform in solid angle at small angles. Whether one reaches a GEM plane is decided by integrating
   its helix from the target in the on-axis field: electrons emitted above **25°** are reflected by the
   magnetic mirror (sin²θ · B_max/B_target > 1); those below 25° have gyroradii of a few cm that shrink
   as ∝1/√B toward the coil, so their transverse excursion stays inside the GEM inner radii.
   **Result: zero direct Møller hits on all six planes.** This is the SoLID design point and the
   estimate reproduces it quantitatively.
2. **Bremsstrahlung of the Møller electrons** while leaving the target (PRM "Seltzer–Berger /
   Bethe–Heitler", complete screening, thin-target yield N_γ ≈ (4/3) t [ln(T/k_min) − 1 + k_min/T]):
   path in LH2 up to the cell wall (1.9 cm / sinθ) plus the Al wall; k_min = 0.1 MeV. The photons keep
   the electron direction and are not bent, so they hit the planes at the Møller angle. Weighted with
   the photon interaction probability in one GEM plane.
3. **Compton scattering of the primary-beam bremsstrahlung** inside the target (PRM "Klein–Nishina"):
   the beam radiates **4.7×10¹² photons/s (k > 0.1 MeV, 560 W)** at 3 µA, all within 0.05 mrad of the
   beam; a fraction Compton-scatters in the remaining target (average path L_t/2) to the GEM angles with
   the Klein–Nishina angular distribution. Scattered photon energy k′ = k/(1 + (k/mₑ)(1 − cosθ)).
4. **Hadronic photo/electro-production**: Weizsäcker–Williams equivalent-photon flux of the electron
   (Q²_max ≈ 1 GeV²) plus real bremsstrahlung photons interacting in the target, times σ_γp(k) (Δ region
   350 µb, second/third resonances 200 µb, Donnachie–Landshoff fit above 2 GeV). σ_eff ≈ 84 µb/proton →
   **1.0×10⁹ events/s at 3 µA**. Secondaries: 1.5 charged (assumed unconfined for p > 0.1 GeV) and one
   π⁰ photon per event, with dN/dΩ ∝ 1/(θ² + (20°)²). This is the most model-dependent term.

Photon → GEM hit: interaction probability in a 0.26 g/cm² plane (Klein–Nishina Compton with Z/A = 0.5,
photoelectric on the 17% Cu mass fraction, pair production above ~20 MeV), times 0.5 for the interaction
to leave a signal in the drift/transfer gaps. Rates per unit area on a plane at distance d use
dA = d² dΩ / cos³θ.

## Results at L = 1.2×10³⁷ (3.0 µA, 15 cm LH2)

Hit rate densities in kHz/mm² (columns: at the inner radius, mid-plane, outer radius); the last
column is the integrated rate over the plane compared with the placeholder table used for the first
luminosity scan.

| Plane | R_in → mid → R_out | Møller-brems γ | Compton γ | hadrons | **total** | placeholder | integrated (analytic / placeholder) |
|---|---|---|---|---|---|---|---|
| 1 | 39 → 91 → 135 cm | 0.17 / 0.10 / 0.07 | 0.21 / 0.10 / 0.05 | 0.33 / 0.09 / 0.04 | **0.71 / 0.30 / 0.16** | 1.05 / 0.19 / 0.09 | 1.61 / 1.25 GHz |
| 2 | 21 → 83 → 135 cm | 0.23 / 0.09 / 0.06 | 0.19 / 0.10 / 0.05 | 0.37 / 0.11 / 0.04 | **0.80 / 0.29 / 0.15** | 3.63 / 0.23 / 0.09 | 1.63 / 1.89 GHz |
| 3 | 25 → 85 → 135 cm | 0.16 / 0.07 / 0.05 | 0.14 / 0.08 / 0.05 | 0.26 / 0.10 / 0.04 | **0.56 / 0.25 / 0.14** | 2.05 / 0.18 / 0.07 | 1.34 / 1.37 GHz |
| 4 | 32 → 88 → 135 cm | 0.10 / 0.05 / 0.04 | 0.09 / 0.06 / 0.04 | 0.16 / 0.08 / 0.04 | **0.35 / 0.19 / 0.12** | 0.94 / 0.12 / 0.05 | 1.00 / 0.87 GHz |
| 5 | 42 → 73 → 99 cm | 0.06 / 0.04 / 0.03 | 0.05 / 0.04 / 0.04 | 0.10 / 0.08 / 0.06 | **0.21 / 0.15 / 0.13** | 0.45 / 0.15 / 0.08 | 0.39 / 0.43 GHz |
| 6 | 53 → 91 → 122 cm | 0.04 / 0.02 / 0.02 | 0.03 / 0.03 / 0.02 | 0.06 / 0.05 / 0.04 | **0.13 / 0.10 / 0.08** | 0.23 / 0.08 / 0.04 | 0.37 / 0.34 GHz |

Fitted parametrisation ρ(R) = ρ₀ (40 cm / R)^α at L = 1.2×10³⁷ (this is what the config file contains):

| Plane | ρ₀ [kHz/mm²] | α | placeholder ρ₀, α |
|---|---|---|---|
| 1 | 0.82 | 1.29 | 1.0, 2.0 |
| 2 | 0.61 | 1.10 | 1.0, 2.0 |
| 3 | 0.49 | 0.99 | 0.8, 2.0 |
| 4 | 0.35 | 0.85 | 0.6, 2.0 |
| 5 | 0.22 | 0.60 | 0.5, 2.0 |
| 6 | 0.15 | 0.59 | 0.4, 2.0 |

Two conclusions:

* **The integrated rates agree with the placeholder to within ±30%** (0.4 to 1.6 GHz per plane at
  3 µA). The placeholder normalisation was not the problem.
* **The radial profile is much flatter than 1/R².** The Møller-brems and Compton photon terms are
  nearly uniform in solid angle, so the hit density falls only as ∝ cos³θ/d². The placeholder therefore
  overestimated the inner-edge densities of planes 2–3 by a factor 3–4.5 and underestimated the outer
  edges by a factor ~2. Since strip occupancy is driven by the densest part of the strip, this matters
  for the strip-readout conclusion (see below).

The three direct sources are of comparable size. The hadronic term dominates at the inner radii but
carries the largest model uncertainty (assumed multiplicity and angular distribution, factor ~3).

## Scaling with luminosity and target length

At fixed target the direct sources scale linearly with current, hence with L. If a higher luminosity is
reached with a longer target, the Compton term and the real-photon part of the hadronic term scale as
(target length)² × current, because both the photon yield and the interaction path grow with the length:

| L = 10³⁹ reached with | Current | Plane-2 total at R_in | Plane-2 integrated | Strip occupancy at R_in (100 ns) |
|---|---|---|---|---|
| 15 cm LH2 | 253 µA | 66 kHz/mm² | 136 GHz | 2.7 |
| 40 cm LH2 (PVDIS-like) | 95 µA | 105 kHz/mm² | 219 GHz | 4.2 |
| placeholder (any) | — | 302 kHz/mm² | 158 GHz | 12 |

So the "dense/long target amplifies the photon background beyond linear luminosity scaling" effect is
real, about a factor 1.6 for 40 cm versus 15 cm at the same luminosity. Note also that 253 µA on a
15 cm LH2 cell is not a realistic CEBAF/target combination; 10³⁹ in practice means a long target.

## Occupancy consequences (100 ns window, inner radius of plane 2)

| L | Readout | Occupancy, placeholder | Occupancy, analytic (15 cm) | Hit survival e^(−2·occ), analytic |
|---|---|---|---|---|
| 1.2×10³⁷ | strip (0.4 mm × 100 cm) | 0.145 | 0.032 | 0.94 |
| 10³⁸ | strip | 1.2 | 0.27 | 0.58 |
| 10³⁹ | strip | 12 | 2.7 | 0.005 |
| 10³⁹ | 1 mm² pixel | 0.030 | 0.0066 | 0.99 |

Per plane and per hit, the analytic profile makes the strip readout at 10³⁸ look "degraded" (58%
survival at the worst point) rather than "dead" (9%). End to end the improvement is much smaller,
because the survival compounds over the five planes a forward track crosses and over the three
leptons, and because the u/v ghost density, which scales as the square of the hit density times the
module area, still floods the pattern recognition. The full chain rerun with this table
(`results_analytic/`) gives for the strip readout:

| L | events | with ≥4 GEM hits on e′, μ⁺, μ⁻ | fully reconstructed | placeholder table |
|---|---|---|---|---|
| 1.2×10³⁷ | 10 000 | 978 | 991 (98% of accepted) | 934 (92%) |
| 10³⁸ | 2 000 | 45 | 19 | 12 |
| 10³⁹ | 2 000 | 0 | 0 | 0 |

So the conclusion for long strips does not change: usable at the J/ψ luminosity, marginal to unusable
at 10³⁸, unusable at 10³⁹. Pixels remain comfortable throughout (see `results_analytic/RESULTS.md`).

## What is not included, and how to get it

The 3×10¹² Møller electrons/s (tens of MeV each, a few W) and the 560 W of bremsstrahlung photons are
swept down the beam pipe. Wherever they stop — beam pipe walls, the target-chamber exit window,
collimators, the endcap nose — they shower, and low-energy photons and neutrons from those showers can
reach the GEMs from behind or from the side. This **albedo** term depends entirely on the beamline
geometry and shielding, cannot be estimated with cross sections alone, and in past SoLID studies was
a significant part of the total. It should be obtained from a Geant4 beam-on-target simulation
(Møller + bremsstrahlung generator, or `solid_gemc`'s `bggen`, through the full geometry with the
beam pipe and shielding) and added to this table; the script accepts `--flat-albedo-kHz-mm2` as a
crude way to add a uniform term meanwhile.

### Low-energy photons: the estimate's weakest point

The default photon cutoff is k_min = 100 keV, and the photon sources are treated as "generated at the
target, then a single interaction in the GEM". Both choices matter, because bremsstrahlung and Compton
spectra are ∝ 1/k, so every decade of photon energy carries a comparable number of photons, and the GEM
response to soft photons is *larger*, not smaller: a 10–30 keV photon is absorbed photoelectrically in
the 5 µm Cu electrodes with high probability, and the photoelectron (range 1–3 µm in Cu) often escapes
into the gas and deposits its full energy, several times a MIP signal. Lowering the cutoff in the
script shows the sensitivity (plane 2, L = 1.2×10³⁷):

| k_min | Compton term at R_in [kHz/mm²] | plane-2 total at R_in | plane-2 integrated |
|---|---|---|---|
| 100 keV (default) | 0.19 | 0.80 | 1.6 GHz |
| 30 keV | 0.54 | 1.20 | 3.0 GHz |
| 10 keV | 1.30 | 2.00 | 5.8 GHz |

The 10 keV row overstates the effect because the script has no attenuation between the target and the
GEM gas: the 0.18 mm Al cell wall alone is about one mean free path at 10 keV, the 50 µm Al GEM entrance
window another 0.3, and a few metres of air about 0.4, so photons below ~15–20 keV mostly never reach
the gas, while 30–100 keV photons largely do. What the script also lacks is the *degradation cascade*:
MeV photons that Compton-scatter repeatedly in the target, the cell, the beam pipe and the detector
frames end up as tens-of-keV photons arriving from all directions, plus K-shell fluorescence from Cu
(8 keV), Fe (6.4 keV, yoke) and Pb (75 keV, shielding). In SoLID's own Geant4 background studies this
soft-photon component, not the direct MeV photons, dominates the GEM hit rate. **Treat the soft-photon
contribution as a factor 2–3 upward uncertainty on the photon terms**, and as the main reason the
direct-source floor here is a floor. It is also the part of the problem where an analytic estimate is
least trustworthy and a Geant4 run with full low-energy EM physics (Livermore or Penelope, production
cuts at a few keV) is required.

Other approximations to keep in mind: the GEM integration window is taken as 100 ns (an APV25-based
readout with a 3 mm drift gap is realistically 100–200 ns, so occupancies could be up to 2× higher);
the field is the analytic coil model without the iron yoke (B at the target could differ by ~30%,
which moves the mirror angle but not the confinement conclusion); photons are assumed to keep the
direction of the radiating electron; hadronic multiplicities and angular distributions are assumed.
