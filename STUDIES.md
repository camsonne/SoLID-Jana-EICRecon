# SoLID DDVCS: studies done so far

This page summarises every study done with this repository: what was asked, how it was done,
what came out, what is still open, and how to reproduce it. Each section links to the detailed
write-up. Start here, then pick an open item at the end to contribute.

**Physics goal.** Measure DDVCS, e p → e′ p′ γ* → e′ p′ μ⁺μ⁻, with SoLID: 11 GeV electrons on a
15 cm LH2 target, at luminosities up to 10³⁸–10³⁹ cm⁻² s⁻¹. The recoil proton is not detected; it
is reconstructed through the missing mass MM² = (k + P − k′ − μ⁺ − μ⁻)².

## Status at a glance

| # | Study | Status | Main result | Details |
|---|---|---|---|---|
| 1 | Missing-mass resolution vs luminosity (fast simulation) | done | pixel GEM readout: σ(MM²) ≈ 56–60 MeV² up to 10³⁹; long strips unusable from 10³⁸ | [results/RESULTS.md](results/RESULTS.md) |
| 2 | First-principles GEM background estimate | done (direct sources only) | agrees with the placeholder table within ±30%, flatter in radius; same conclusions | [docs/background_estimate.md](docs/background_estimate.md), [results_analytic/COMPARISON.md](results_analytic/COMPARISON.md) |
| 3 | EpIC physics generator and DDVCS angles | done, tested with an EpIC-like fixture only | EpIC HepMC3 input; reconstructed φ, θ_ℓ, φ_ℓ match EpIC's conventions to 1e-12 | [README](README.md#physics-events-from-epic) |
| 4 | ACTS tracking | done | same resolution as the native tracker, equal or better efficiency, 3–7× faster | [docs/acts_tracking.md](docs/acts_tracking.md) |
| 5 | GNN tracking (Interaction Network, ONNX) | done | same resolution; fewest fake tracks at 10³⁹ (0.2% against 3.4% for ACTS) | [docs/gnn_tracking.md](docs/gnn_tracking.md) |
| 6 | Full Geant4 simulation | first run done | acceptance 12.3% (fast sim 14.2%), σ(MM²) 83–88 MeV², reconstruction efficiency of accepted events 68–73% | [docs/geant4_simulation.md](docs/geant4_simulation.md) |

All work is on branch `ccr-ca3416b1-vmzs77`.

## Software chain

```
 generator                     detector simulation                    reconstruction (JANA2 plugin solid_ddvcs)
 ---------                     -------------------                    -----------------------------------------
 internal phase space ──┐      fast simulation (RK4 transport,  ──┐   tracking: acts | gnn | native
 (solid_ddvcs_gen)      ├──►   multiple scattering, GEM          ├──► DDVCS reconstruction (e', mu+, mu-, MM², Q²,
 EpIC (DDVCS+BH, HepMC3)┘      background model)                 │    x_B, Q'², t, phi, theta_l, phi_l)
                               DD4hep/Geant4 (ddsim) → EDM4hep ──┘    analysis: per-event CSV + JSON summary
```

| Component | Where | Selected with |
|---|---|---|
| Geometry (DD4hep compact) | `geometry/` | `$DETECTOR_PATH/solid_ddvcs.xml` |
| Geometry numbers for fast sim and reconstruction | `src/services/SoLIDGeometry_service` | — |
| Magnetic field | `src/services/SoLIDField_service` | `field:mode=analytic` (coil only) or `map` (CLEO-II `solenoid_CLEOv9.dat`) |
| Event source | `src/fastsim/`, `src/io/` | `solid:source=fastsim` (`fastsim:hepmc_in=` for EpIC) or `edm4hep` |
| Tracking | `src/reco/`, `src/acts/`, `src/gnn/` | `tracking:engine=acts\|gnn\|native` |
| Reconstruction and analysis | `src/reco/DDVCSReco_factory`, `src/analysis/` | `analysis:output`, `analysis:summary` |

Build options: `-DUSE_ACTS=ON` (ACTS ≥ 46), `-DUSE_ONNX=ON` (ONNX Runtime ≥ 1.17),
`-DUSE_EDM4HEP=ON` (podio ≥ 1.0). The first and last switch the build to C++20. See the
[README](README.md#building).

---

## 1. Missing-mass resolution vs luminosity (fast simulation)

**Question.** Can SoLID reconstruct the DDVCS missing mass at 10³⁸–10³⁹ cm⁻² s⁻¹, and does the
GEM readout (long u/v strips or 1 mm² pads) matter?

**Method.**
* Internal phase-space generator; fast simulation with RK4 transport in the solenoid field, and
  multiple scattering in the target, air and GEM planes.
* GEM hit efficiency, plus a luminosity-scaled background model: hit rate density per plane vs
  radius in a 100 ns window. Strips lose hits by cluster overlap and gain u/v ghost hits.
* Native tracker: calorimeter-seeded, global χ² fit with the beam-line constraint.
* 10,000 events per point (2000 for saturated strips) at 1.2×10³⁷, 10³⁸ and 10³⁹, for each readout.

**Results** (`results/RESULTS.md`):

| L [cm⁻²s⁻¹] | readout | ε(reco \| geom. acc.) | σ(MM²) core [MeV²] | GEM plane 2 occupancy at R_in |
|---|---|---|---|---|
| 1.2×10³⁷ | pixel | 0.82 | 59 | 0.000 |
| 10³⁸ | pixel | 0.79 | 56 | 0.003 |
| 10³⁹ | pixel | 0.73 | 60 | 0.030 |
| 1.2×10³⁷ | strip | 0.62 | 56 | 0.145 |
| 10³⁸ | strip | 0.04 | — | 1.21 |
| 10³⁹ | strip | 0 | — | 12.1 |

**Conclusion.**
* Pixel readout keeps σ(MM²) ≈ 56–60 MeV² (σ(M_X) ≈ 30 MeV) flat up to 10³⁹.
* Long strips are background-limited from 10³⁸, so the highest DDVCS luminosities need finer
  segmentation.

**Caveats.** The background rates are a placeholder table (factor 2–3 uncertainty). PID is truth
based, and there is no energy loss or radiation. The fast simulation turned out to be optimistic
compared with Geant4 (study 6).

**Reproduce.** `scripts/run_ddvcs_lumi_scan.sh results 10000 4` (plots with
`scripts/plot_ddvcs_results.py`).

## 2. First-principles GEM background estimate

**Question.** Is the placeholder background table credible?

**Method.** Analytic rates (`scripts/estimate_background_rates.py`) from cross sections folded
with beam, target, field and GEM geometry:
* Møller scattering, including confinement of the Møller electrons in the solenoid;
* bremsstrahlung of Møller electrons;
* Compton scattering of the beam bremsstrahlung;
* hadronic photo- and electro-production.

The scan was then rerun with the analytic table (`config/gem_background_rates_analytic.txt`).

**Results.**
* Møller electrons never reach the GEMs (magnetic mirror and small gyroradii).
* Totals agree with the placeholder within ±30%, with a much flatter radial profile.
* The scan gives the same conclusions:
  * Pixel: ε 0.72–0.82 and σ(MM²) 56–60 MeV² with either table.
  * Strips: 6% efficiency at 10³⁸, zero at 10³⁹ (`results_analytic/COMPARISON.md`).

**Caveats.**
* The albedo from the beamline, collimators and endcap is not included; it needs a Geant4
  beam-on-target simulation.
* The soft-photon component dominates and is sensitive to the low-energy cutoff (see the doc).

## 3. EpIC physics generator and DDVCS angular variables

**Question.** Can EpIC (PARTONS-based DDVCS + Bethe–Heitler, with beam helicity) feed the chain,
and are the angles needed for asymmetries reconstructed?

**Method.**
* `HepMC3Reader` reads EpIC HepMC3 files:
  * EpIC's beam along −z is rotated to SoLID's +z;
  * particles are identified from the vertex structure;
  * the vertex is placed in the target;
  * the integrated cross section is read from the end of the file.
* `DDVCSKinematicsCalc.h` computes Q², x_B, Q′², t, φ (Trento), θ_ℓ and φ_ℓ, identically for
  truth and reconstruction.
* `config/epic/epic_scenario_SoLID_DDVCS.xml` is a ready-made EpIC scenario: 11 GeV, target at
  rest, μ⁺μ⁻ decay, SoLID kinematic range.
* `solid_hepmc_to_solid` converts EpIC files for ddsim.

**Results.**
* On events built step by step with EpIC's own construction code (`tests/make_epic_like_hepmc.py`),
  the reconstructed truth reproduces EpIC's φ, θ_L and φ_L to 1e-12.
* Angular resolutions in the fast simulation: σ(φ) ≈ 0.024 rad, σ(cos θ_ℓ) ≈ 0.013,
  σ(φ_ℓ) ≈ 0.023 rad.

**Caveats.**
* **EpIC itself has not been run here** (it needs PARTONS and ROOT); the test uses an EpIC-like
  fixture.
* EpIC's example uses `hadron_energy = 1.00001*Mp`, which corresponds to an 11.05 GeV beam. The
  SoLID scenario uses `Mp`.

## 4. ACTS tracking

**Question.** Can SoLID use the EIC-standard ACTS tracking?

**Method.** ACTS 46.8.1 (the eic-shell version), see `src/acts/`:
* Geometry: GEM planes as planar rings with material in an air volume.
* Seeding: Seeding2 triplet seeding, with the collision region at the target.
* Seed pre-fit, then a Combinatorial Kalman Filter from the beam line, smoothing and greedy
  ambiguity resolution.
* SoLID vertex treatment: beam-line constraint and target multiple scattering.
* Calorimeter match: the track must point to a calorimeter cluster.

**Results** (fast simulation, identical events, `docs/acts_tracking.md`):
* Same resolution as native: σ(MM²) 56–58 MeV², δp/p 0.77%.
* Efficiency equal or better: 97.5% vs 94.0% at 10³⁹ pixel, 98.8% vs 90.8% at 1.2×10³⁷ strip.
* 3–7× faster.
* Three SoLID-specific fixes were needed:
  * The seed z0 window must allow for helix bending.
  * Seeds need a pre-fit, otherwise the CKF picks wrong hits.
  * The surfaces need a margin beyond the active edges.

**Caveats.**
* At 10³⁹, ACTS keeps 3.4% fake tracks even with the calorimeter match (native: 0.5%).
* In Geant4 the calorimeter match must be off (study 6).

## 5. GNN tracking

**Question.** Is a graph neural network tracker (Exa.TrkX style) better?

**Method.** See `src/gnn/` and `scripts/gnn/`:
* Hit graph: edges between neighbouring planes, or skipping one, with geometric cuts that keep
  99.9% of lepton segments.
* Interaction Network edge classifier (160k parameters), trained in PyTorch on 16k fast-simulation
  events and exported to ONNX (`config/gnn/solid_gem_in.onnx`).
* Walkthrough track building, then the native fitter.

**Results** (`docs/gnn_tracking.md`):
* Validation edge efficiency and purity: 0.9995 / 0.9998 at 1.2×10³⁷, 0.967 / 0.817 at 10³⁹.
* Same resolution as the other engines.
* Efficiency equal at low occupancy; 95.9% at 10³⁹ pixel (native 94.0%, ACTS 97.5%; the
  differences are within about 1% statistical error).
* **Fewest fake tracks at 10³⁹: 0.20%**, against 0.51% (native) and 3.4% (ACTS).
* Fastest below 10³⁹.
* Walkthrough building was essential: mutual-best chaining gave only 75% at 10³⁹.

**Verdict.** Not better in resolution or clearly in efficiency; better in purity at high
occupancy. It is trained on the fast simulation only.

## 6. Full Geant4 simulation

**Question.** Does the full simulation work, and how optimistic is the fast simulation?

**Method.**
* DD4hep 1.38 with Geant4 11.4 and EDM4hep output, plus `solid_dd4hep`.
* CLEO-II field map. Stack recipe in `scripts/geant4/setup_stack.sh`.
* 2000 generated events, run through both Geant4 and the fast simulation, with each tracking
  engine.

**First run.** The geometry had never been loaded. Fixed:
* compact syntax errors;
* **no magnet in the geometry:** the field map contains the yoke return field, but there was no
  iron, so large-angle muons bent back in air (now `geometry/compact/magnet_ddvcs.xml`, converted
  from `solid_gemc`);
* the large-angle muon scintillator not sensitive;
* SPD/MRPC overlaps;
* reconstruction mismatches: field model, calorimeter cluster threshold (100 MeV, while a muon MIP
  deposits 36 MeV) and cluster depth.

**Results** (same events; `docs/geant4_simulation.md`):

| | fast simulation | Geant4 |
|---|---|---|
| geometric acceptance | 14.2% | 12.3% (muons range out in the iron) |
| σ(MM²) core | 56–60 MeV² | 83–88 MeV² |
| σ(p)/p | 0.80% | ≈ 1.05% |
| ε(reco \| accepted): ACTS / GNN / native | 98% / 98% / 94% | 73% / 68% / 59% |

The Geant4 efficiency losses, per lepton:
* about 9% of electrons are not tracked (bremsstrahlung);
* 8–12% of muons fail identification, which extrapolates the track instead of using the
  muon-detector hits;
* the calorimeter match costs a further 10%.

The MM² peak sits 0.03 GeV² high, because energy loss in the target, Cherenkovs and windows is
not corrected.

**Known geometry issue.** The SIDIS-layout LGC/HGC tanks of `solid_dd4hep` overlap the long-endcap
nose cone by 8.0 and 2.4 cm.

---

## Consolidated conclusions so far

1. The **missing-mass method works**: Geant4 gives σ(MM²) ≈ 85 MeV², σ(M_X) ≈ 45 MeV. The fast
   simulation's 56–60 MeV² is optimistic by about 1.5×.
2. **GEM readout**: long strips fail from 10³⁸; pixels or short strips are needed for the
   high-luminosity running.
3. **Acceptance** is limited by the muon detector: muons below about 3 GeV stop in the iron, and
   large-angle muons must cross the yoke.
4. **Tracking**: ACTS is the recommended default (standard, efficient, fast). Use it without the
   calorimeter match on realistic simulation. The GNN is the best option if fake tracks become the
   limitation.
5. The current reconstruction loses about 30% of accepted Geant4 events, mostly in muon
   identification and electron tracking. These are the most valuable things to fix (below).

## Open items: where to contribute

Ordered by expected impact. Each item names the files involved and how to judge success.
Compare against the Geant4 numbers in study 6, on the same 2000-event sample (see "Reproducing").

| # | Item | Files | Done when |
|---|---|---|---|
| A | **Muon identification from muon-detector hits**: match the track to the simulated muon-detector hits (and calorimeter MIP) instead of the extrapolation and the 2 GeV cut | `src/io/SoLIDEDM4hep_source.cc`, `src/reco/TrackUtils.cc`, `src/reco/DDVCSReco_factory.cc` | muon "tracked + identified" in the analysis JSON ≥ 97% on Geant4 |
| B | **Electron fit with ACTS's Gaussian-sum filter** (bremsstrahlung) | `src/acts/ActsTrackFinder_factory.cc` | electron tracked fraction ≥ 97%, smaller σ(p)/p for e′ |
| C | **Vertex energy-loss correction**: target, Cherenkov and window material in the ACTS geometry, or a dE/dx correction at the vertex | `src/acts/`, `src/reco/` | ⟨MM²⟩ within 0.005 GeV² of Mp² on Geant4 |
| D | **DDVCS layout of LGC/HGC** in `solid_dd4hep` (or decide whether DDVCS runs without them) | upstream `JeffersonLab/solid_dd4hep`, `geometry/solid_ddvcs.xml` | `checkOverlaps` clean |
| E | **Realistic background**: Geant4 beam-on-target run (with albedo) → `gem:background_file` table, and mixing of background hits into Geant4 events | `scripts/`, `config/` | rate table per plane vs radius from simulation |
| F | **Run EpIC** with `config/epic/epic_scenario_SoLID_DDVCS.xml`, both helicities; check the reader on real output; first beam-spin-asymmetry projection | `src/generator/HepMC3Reader.cc`, analysis | rates and BSA statistical errors vs (Q², Q′², t, φ) |
| G | **Retrain the GNN on Geant4 hits**, including showers and secondaries | `scripts/gnn/`, `-Pgnn:dump` on EDM4hep input | GNN efficiency on Geant4 ≥ ACTS |
| H | **ACTS geometry from DD4hep** (ACTS DD4hep plugin, as EICrecon does) instead of building it in code | `src/acts/SoLIDActs_service.cc` | same performance, geometry read from `solid_ddvcs.xml` |
| I | **Cherenkov and calorimeter PID response** (currently truth PID; optical photons available with `--optical`) | `scripts/run_ddsim_ddvcs.py`, `src/io/`, `src/reco/` | e/π and μ/π separation, and π→μ misidentification background |
| J | **Physics analysis**: binning, unfolding, asymmetry extraction, systematics | new `src/analysis/` code, `docs/` | projected observables |

## How to contribute

1. **Branch** from the current work, `git checkout -b <you>/<topic> origin/ccr-ca3416b1-vmzs77`,
   and open a pull request against it, or against `main` once it is merged.
2. **Build** the configuration you need (see the README). The fast-simulation chain needs only
   JANA2. ACTS, ONNX Runtime and EDM4hep are optional, and are available inside eic-shell.
3. **Test**: `ctest --test-dir build` runs smoke tests of every engine that is built, the
   generator, the EpIC reader and the GNN graph export. Add a test for new code.
4. **Keep the interfaces**: all tracking engines produce the same `Track` objects
   (`src/datamodel/SoLIDDataModel.h`), and truth matching and PID flags go through
   `src/reco/TrackUtils`. This keeps engines interchangeable and comparisons fair.
5. **Report results reproducibly**:
   * State the commit, the command line, the number of events and the seeds.
   * Commit the JSON summary (small), not the per-event CSVs (git-ignored).
   * Add a `docs/<topic>.md`, and a row in "Status at a glance" above.
   * Compare on identical events, using `fastsim:hepmc_in` or the same EDM4hep files.
6. **Conventions**: cm, GeV, ns and T inside JANA (ACTS uses mm). Parameters are
   `-P<group>:<name>`; `jana -Pplugins=solid_ddvcs -Pjana:parameter_verbosity=2` lists them.

## Reproducing the Geant4 comparison

```bash
source scripts/geant4/env.sh                       # stack from scripts/geant4/setup_stack.sh or eic-shell
build/src/solid_ddvcs_gen -n 2000 -s 2026 -o ddvcs2000.hepmc
for j in 0 1 2 3; do                               # 4 x 500 events, ~6 min on 4 cores
  python3 scripts/run_ddsim_ddvcs.py --compact $DETECTOR_PATH/solid_ddvcs.xml -i ddvcs2000.hepmc \
          -n 500 -s $((j*500)) --seed $((1000+j)) -o g4_part$j &
done; wait
FM="-Pfield:mode=map -Pfield:map=$SOLID_FIELD_MAP"
# Geant4, ACTS (built with -DUSE_ACTS=ON -DUSE_EDM4HEP=ON)
jana -Pplugins=solid_ddvcs -Psolid:source=edm4hep $FM -Ptracking:engine=acts -Pacts:energy_loss=1 \
     -Pacts:ec_match_cm=0 -Panalysis:summary=g4_acts.json g4_part*.edm4hep.root
# the same events in the fast simulation
jana -Pplugins=solid_ddvcs $FM -Ptracking:engine=acts -Pfastsim:hepmc_in=ddvcs2000.hepmc \
     -Panalysis:summary=fs_acts.json
```

Versions used: JANA2 master (2026-10), ACTS 46.8.1, ONNX Runtime 1.31, DD4hep 1.38, Geant4 11.4.2,
ROOT 6.36.14, podio 1.8, EDM4hep 1.2, PyTorch 2.14 (GNN training).
