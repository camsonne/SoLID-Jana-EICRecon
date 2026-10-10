# SoLID-Jana-EICRecon

JANA2-based simulation and reconstruction chain for **SoLID** in the **DDVCS** configuration
(`e p → e' p' γ* → e' p' μ⁺μ⁻`, 11 GeV beam, 15 cm LH2 target at z = −315 cm), built to
quantify the **missing-mass resolution of the undetected recoil proton** with full tracking at
luminosities of **10³⁸ and 10³⁹ cm⁻² s⁻¹** (and the 1.2×10³⁷ J/ψ reference point).

The results of the luminosity scan are in [`results/RESULTS.md`](results/RESULTS.md).

## What is in the repository

| Path | Content |
|---|---|
| `geometry/solid_ddvcs.xml`, `geometry/compact/` | DD4hep compact description of the DDVCS configuration: LH2 target, DDVCS GEM tracker (6 planes, planes 1–4 extended to R = 135 cm), LGC/HGC, LAEC/FAEC, SPD, MRPC, DDVCS muon detector (iron + scintillator). Detector plugins (`SoLID_*`) come from [JeffersonLab/solid_dd4hep](https://github.com/JeffersonLab/solid_dd4hep); positions follow `solid_gemc` (`gem_moved/solid_SIDIS_DDVCS_gem`, `prototype/muon_pro`). |
| `scripts/run_ddsim_ddvcs.py` | DDG4 steering for the full Geant4 simulation (HepMC3 in, EDM4hep/podio out). |
| `src/generator/` | DDVCS kinematics generator (`solid_ddvcs_gen`, HepMC3 ASCII output) — phase-space with importance sampling in the SoLID DDVCS region; weight 1, no amplitude. `HepMC3Reader` reads physics events from [EpIC](https://github.com/pawelsznajder/epic) (DDVCS + BH cross section), `solid_hepmc_to_solid` converts them for ddsim, and `DDVCSKinematicsCalc.h` computes Q², x_B, Q′², t, the Trento φ and the lepton-pair angles θ_ℓ, φ_ℓ with EpIC's conventions, for truth and reconstruction alike. |
| `config/epic/` | EpIC scenario for SoLID DDVCS (11 GeV e⁻ on LH2, μ⁺μ⁻ pair). |
| `src/services/` | `SoLIDGeometry_service` (DDVCS geometry parameters) and `SoLIDField_service` (CLEO-II solenoid: Biot–Savart model of the coil packs from `solid_gemc/field/CLEOv9.am`, or the `solenoid_CLEOv9.dat` map when available). |
| `src/fastsim/` | `SoLIDFastSim_source`: generator + RK4 transport in the solenoid field, multiple scattering (target, air, GEM planes), GEM hit smearing/efficiency, calorimeter clusters, muon-detector acceptance and a **luminosity-scaled background model** (`BackgroundModel`: rate density per plane vs radius, GEM time window, strip vs pixel readout, cluster-overlap losses, u/v ghost hits). |
| `src/acts/` | ACTS tracking (`tracking:engine=acts`, default when built with `-DUSE_ACTS=ON`): `SoLIDActs_service` (ACTS tracking geometry of the GEM planes and the SoLID field as an ACTS field provider), `ActsTrackFinder_factory` (ACTS triplet seeding, seed pre-fit, Combinatorial Kalman Filter, greedy ambiguity resolution, beam-line constraint). |
| `src/reco/` | Native tracker (`tracking:engine=native`): `TrackPropagator` (RK4), `TrackFitter` (global χ² fit of (q/p, θ, φ, z_v) with the SoLID beam-line vertex constraint and a full multiple-scattering covariance), `TrackFinder_factory` (calorimeter-seeded progressive track finding). `DDVCSReco_factory` (e', μ⁺, μ⁻ selection and missing mass) serves both engines. |
| `src/analysis/` | `DDVCSAnalysis_processor`: per-event CSV and JSON summary (acceptance, efficiency, resolutions). |
| `src/io/SoLIDEDM4hep_source` | Reader for ddsim output (EDM4hep), digitises `GEMTrackerHits`, clusters the calorimeter hits, flags muon-detector hits. Built with `-DUSE_EDM4HEP=ON`. |
| `scripts/run_ddvcs_lumi_scan.sh`, `scripts/plot_ddvcs_results.py` | Luminosity scan driver and plots/tables. |
| `config/gem_background_rates*.txt` | Background-rate tables for `gem:background_file`: the placeholder and the analytic estimate. |
| `scripts/estimate_background_rates.py`, `docs/background_estimate.md` | First-principles estimate of the direct GEM background from cross sections. |

## Building

Requires CMake ≥ 3.16, a C++17 compiler and [JANA2](https://github.com/JeffersonLab/JANA2) (v2.4+;
the `Input<>/Output<>/Parameter<>` component API). ROOT is **not** required.

```bash
# JANA2 (if not already available)
git clone https://github.com/JeffersonLab/JANA2 && cmake -S JANA2 -B JANA2/build -DCMAKE_INSTALL_PREFIX=$PWD/jana2-install \
    -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF && cmake --build JANA2/build -j4 && cmake --install JANA2/build

# this package
cmake -S . -B build -DCMAKE_PREFIX_PATH=$PWD/jana2-install -DCMAKE_INSTALL_PREFIX=$PWD/install
cmake --build build -j4
export PATH=$PWD/jana2-install/bin:$PATH JANA_PLUGIN_PATH=$PWD/build/src
ctest --test-dir build
```

Add `-DUSE_EDM4HEP=ON` (needs podio and EDM4hep ≥ 0.99) to build the ddsim reader.

**ACTS tracking** (`-DUSE_ACTS=ON`): needs [ACTS](https://github.com/acts-project/acts) ≥ 46 (Core
only; the code is developed against 46.8.1, the version in eic-shell) and switches the build to
C++20. ACTS is then the default tracking engine; `-Ptracking:engine=native` selects the built-in
tracker.

```bash
# inside eic-shell ACTS is already installed; otherwise (needs Eigen3 >= 3.4 and Boost >= 1.78)
git clone --branch v46.8.1 --depth 1 https://github.com/acts-project/acts
cmake -S acts -B acts/build -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=$PWD/acts-install \
      -DACTS_BUILD_EXAMPLES=OFF -DACTS_BUILD_UNITTESTS=OFF -DACTS_BUILD_INTEGRATIONTESTS=OFF
cmake --build acts/build -j4 && cmake --install acts/build

cmake -S . -B build -DUSE_ACTS=ON -DCMAKE_PREFIX_PATH="$PWD/jana2-install;$PWD/acts-install"
```

## Running the fast-simulation + reconstruction chain

```bash
# one configuration
jana -Pplugins=solid_ddvcs -Pnthreads=4 -Pfastsim:nevents=10000 \
     -Psolid:luminosity=1e38 -Pgem:readout=strip \
     -Panalysis:output=L1e38_strip_events.csv -Panalysis:summary=L1e38_strip_summary.json

# full scan (1.2e37, 1e38, 1e39) x (strip, pixel) + plots
scripts/run_ddvcs_lumi_scan.sh results 10000 4
```

Useful parameters (all have `-Pname=value` form; `jana -Pplugins=solid_ddvcs -Pjana:parameter_verbosity=2` lists them):

| Parameter | Default | Meaning |
|---|---|---|
| `solid:luminosity` | 1.2e37 | Luminosity used by the background model [cm⁻² s⁻¹] |
| `gem:readout` | strip | `strip` (u/v strips, 0.4 mm pitch, length = module radial size ≤ 100 cm) or `pixel` (1 mm² pads) |
| `gem:time_window_ns` | 100 | GEM hit integration window |
| `gem:strip_length_cm`, `gem:pad_area_mm2` | 0 (auto), 1.0 | Readout granularity |
| `gem:background_file` | — | Per-plane rate parametrisation (see `config/gem_background_rates.txt`) |
| `gem:resolution_um`, `gem:efficiency`, `gem:x0_fraction` | 70, 0.97, 0.007 | GEM plane properties |
| `field:mode`, `field:map`, `field:central_tesla` | analytic, —, 1.5 | Solenoid model (`map` reads `solenoid_CLEOv9.dat`, GEMC ASCII `r z Br Bz` format) |
| `tracking:min_hits`, `tracking:chi2ndf_max` | 4, 10 | Track selection |
| `solid:muon_pmin` | 2.0 | Minimum momentum to traverse FAEC + 108 cm Fe to the last muon scintillator [GeV] |
| `ddvcs:Q2_min/max`, `ddvcs:Qp2_min/max`, `ddvcs:xB_min/max`, `ddvcs:tabs_max` | 1–6, 2–9, 0.10–0.55, 1.5 | Generator ranges |
| `fastsim:hepmc_out` | — | Also write the generated events as HepMC3 |
| `fastsim:hepmc_in` | — | Read the events from a HepMC3 ASCII file (EpIC) instead of the internal generator |
| `fastsim:hepmc_vertex` | target | `target`: vertex uniform along the LH2 cell (EpIC writes none); `file`: use the file's vertex |

## Physics events from EpIC

[EpIC](https://github.com/pawelsznajder/epic) (PARTONS-based; `DDVCSGeneratorService`, model
`DDVCSProcessDMSW22`) generates unweighted DDVCS + Bethe–Heitler events with beam helicity, and is the
event generator for the physics analysis. The internal generator remains for pure acceptance and
resolution studies.

```bash
# 1. generate (EpIC installation or the partons/epic Docker image); one file per beam helicity
epic --scenario=config/epic/epic_scenario_SoLID_DDVCS.xml --seed=1        # -> epic_SoLID_DDVCS.hepmc

# 2a. fast simulation + reconstruction on the EpIC events
jana -Pplugins=solid_ddvcs -Pfastsim:hepmc_in=epic_SoLID_DDVCS.hepmc -Psolid:luminosity=1e38 -Pgem:readout=pixel \
     -Panalysis:output=epic_L1e38_pixel_events.csv -Panalysis:summary=epic_L1e38_pixel_summary.json

# 2b. or convert for the Geant4 path (beam along +z, vertex in the target) and run ddsim on it
build/src/solid_hepmc_to_solid -i epic_SoLID_DDVCS.hepmc -o ddvcs_solid.hepmc
```

What the reader does with EpIC files:

* **Frame.** EpIC puts the lepton beam along −z (collider convention). Such events are rotated by π
  around y, so the beam points along +z as in SoLID; kinematics and angles are unchanged.
* **Fixed target.** Use `hadron_energy = Mp` (as in the provided scenario). The EpIC example's
  `1.00001*Mp` corresponds to an 11.05 GeV fixed-target beam and leaves the proton moving at 4 MeV/c.
  The source warns when the target is not at rest or the beam energy differs from `solid:beam_energy`.
* **Event record.** e, p, e′, p′, γ*, ℓ⁻, ℓ⁺ are identified from the vertex structure (status 4 beams,
  status 1 final state, ℓ pair produced by the outgoing γ*). Additional final-state particles, such as
  photons from EpIC radiative corrections, are kept. The pair must be μ⁺μ⁻ (`decay_type = mu-`).
* **Vertex.** EpIC writes none, so the vertex is placed uniformly along the 15 cm target.
* **Cross section.** EpIC events are unweighted. The integrated cross section is written after the
  last event and logged at the end of the run; the event weight goes to `GenKinematics::weight`, the
  `weight` CSV column and `weight_sum` in the JSON summary.
* **Angles.** `phi` (Trento, target rest frame), `costhl` and `phil` (ℓ⁻ relative to the recoil
  proton in the ℓℓ rest frame) in the analysis CSV are exactly EpIC's `phi`, `cos(thetaL)` and
  `phiL`. `tests/test_hepmc_reader.cc` checks this to 1e-12 on events built with EpIC's own
  construction (`tests/make_epic_like_hepmc.py`). EpIC itself is not run in the tests.

## Full Geant4 simulation path

```bash
export SOLID_DD4HEP=/path/to/solid_dd4hep DETECTOR_PATH=$PWD/geometry     # inside eic-shell
build/src/solid_ddvcs_gen -n 10000 -o ddvcs.hepmc
python3 scripts/run_ddsim_ddvcs.py --compact geometry/solid_ddvcs.xml -i ddvcs.hepmc -n 10000 -o ddvcs_sim
jana -Pplugins=solid_ddvcs -Psolid:source=edm4hep -Psolid:luminosity=1e38 ddvcs_sim.edm4hep.root
```

The same tracking, DDVCS reconstruction and analysis run on the Geant4 hits; the background
model is applied on top of the simulated hits in the same way as in the fast simulation.
The compact files have been written against the `solid_dd4hep` plugin interfaces but could not
be executed in the environment used to develop this package (no DD4hep/Geant4 available), so
expect to iterate on them the first time they are loaded.

## Physics and detector model

* **Geometry** (all z relative to the solenoid centre): target 15 cm LH2 at −315 cm (1.89 cm radius Al cell);
  GEM planes at z = −175, −150, −119, −68, 5, 92 cm with R_in = 39, 21, 25, 32, 42, 53 cm and
  R_out = 135, 135, 135, 135, 99, 122 cm; LAEC front face z = −65 cm (R 83–137 cm); FAEC front
  face z = 425 cm (R 98–230 cm); forward muon detector: 3 × 36 cm Fe with scintillator annuli
  (R 80–285 cm) at z = 661, 725, 789 cm; large-angle muon scintillator cylinder R = 290 cm, z = 209–620 cm.
  Forward tracks (θ ≈ 7.5–16°) cross planes 2–6, large-angle tracks (θ ≈ 16–28°) planes 1–4.
* **Field**: Biot–Savart field of the CLEO-II coil packs (r = 152.5/154.1 cm, z = ±173.8 cm) scaled
  to 1.5 T at the centre. The iron yoke is not modelled, so the fringe field at the target differs
  from the real map; use `field:mode=map field:map=solenoid_CLEOv9.dat` when the map is available.
* **Tracking**: hits are 2D points with 70 µm resolution; the fit propagates with RK4 in the field
  and uses a measurement covariance that includes the multiple scattering in the target (LH2 +
  Al wall), the air between planes and the GEM planes (0.7 % X₀ each), so χ²/ndf ≈ 1 and the
  parameter errors are realistic. Tracks are seeded from a calorimeter cluster + a GEM hit + the
  vertex constraint (as in the SoLID tracking strategy), extended plane by plane inside a window
  given by the predicted-position uncertainty, refitted and cleaned of outliers.
* **Background**: hit rate density ρ_i(R) = ρ₀,ᵢ (40 cm / R)² · L/1.2×10³⁷ with ρ₀ ≈ 0.4–1 kHz/mm²
  (order of magnitude of the SoLID pCDR estimates; factor 2–3 uncertainty — replace with numbers
  from a dedicated `solid_gemc` background run through `gem:background_file`). Hits within the
  GEM time window are generated in a ±3 cm road around each true crossing. For strip readout the
  channel occupancy o = ρ·T·pitch·length removes true hits with probability 1 − e^(−2o) and
  produces u/v ghost hits of density (ρT)²·A_module; pixel readout uses o = ρ·T·A_pad and no ghosts.
* **PID**: electrons must reach a calorimeter, muons the muon detector behind ≥ 108 cm of iron
  (p > 2 GeV); the e/μ/π separation itself is assumed perfect (truth association) — the
  Cherenkov, calorimeter and muon-detector responses are not simulated. Background in the
  calorimeters is not simulated.
* **Not included**: bremsstrahlung/energy loss of the leptons in the target and detectors, the
  physical DDVCS + Bethe–Heitler cross section in the internal generator (all events have weight 1;
  use EpIC events for it), and pile-up in the
  calorimeters.

## Results

Full tables and figures: [`results/RESULTS.md`](results/RESULTS.md). Six configurations were run
(luminosities 1.2×10³⁷, 10³⁸, 10³⁹ cm⁻² s⁻¹ × strip/pixel GEM readout, 10 000 events each, 2000 for
the saturated strip configurations):

* **Pixel readout (1 mm² pads)** tracks the DDVCS missing mass with a Gaussian-core resolution of
  **56–60 MeV²** essentially flat from 1.2×10³⁷ to 10³⁹ cm⁻² s⁻¹ (σ(M_X) ≈ 30 MeV). The fraction of
  geometrically-accepted events with e′, μ⁺ and μ⁻ all reconstructed falls from 0.82 to 0.73 as GEM
  occupancy grows with luminosity.
* **Strip readout (0.4 mm pitch, up to 100 cm length)** is background-limited in this model: the
  channel occupancy at the inner radius of GEM plane 2 is already ≈ 0.15 at 1.2×10³⁷ (comparable
  resolution to pixel, 56 MeV²) but reaches ≈ 1.2 at 10³⁸ and ≈ 12 at 10³⁹, so cluster overlap
  destroys the real hits faster than tracks can be found: efficiency drops to 0.04 at 10³⁸ and 0 at
  10³⁹ cm⁻² s⁻¹. At the highest SoLID DDVCS luminosities this points to finer segmentation
  (shorter strips or pixels) rather than the SIDIS-style long strips.

These numbers use the placeholder background-rate table compiled into `BackgroundModel`. A
first-principles estimate of the direct background (Møller, bremsstrahlung, Compton, hadronic
photoproduction, with the solenoid confinement of the Møller electrons) is in
[`docs/background_estimate.md`](docs/background_estimate.md), produced by
`scripts/estimate_background_rates.py`. It agrees with the placeholder totals to ±30% but has a much
flatter radial profile; `config/gem_background_rates_analytic.txt` carries that profile and the scan
with it is in [`results_analytic/RESULTS.md`](results_analytic/RESULTS.md), with a side-by-side table in
[`results_analytic/COMPARISON.md`](results_analytic/COMPARISON.md). The two tables give the same
physics conclusions: pixel efficiency and resolution agree within statistics at every luminosity
(0.72 to 0.82, 56 to 60 MeV²), and long strips go from usable at 1.2×10³⁷ to 6% efficiency at 10³⁸
and zero at 10³⁹. Neither table contains the beamline/collimator albedo, which requires a Geant4
beam-on-target simulation.
