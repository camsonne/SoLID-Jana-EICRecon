# Full Geant4 simulation of the SoLID DDVCS configuration

First end-to-end run of the full simulation path: DDVCS events → DD4hep/Geant4 (`ddsim`-style
DDG4 steering) → EDM4hep → JANA2 reconstruction with each tracking engine. The run was compared with
the fast simulation on **the same generated events**.

## Software

| Package | Version | Source |
|---|---|---|
| Geant4 | 11.4.2 (with data) | conda-forge |
| ROOT | 6.36.14 | conda-forge |
| HepMC3 | conda-forge | conda-forge |
| podio | 1.8 | source |
| EDM4hep | 1.2 | source |
| DD4hep | 1.38 (Geant4, HepMC3 and EDM4hep output) | source |
| solid_dd4hep | main | JeffersonLab |
| Field map | CLEO-II `solenoid_CLEOv9.dat` | JeffersonLab/solid_gemc on GitHub |

`scripts/geant4/setup_stack.sh` reproduces the stack outside eic-shell, and
`scripts/geant4/env.sh` sets up the environment. Inside eic-shell only `solid_dd4hep` needs
building. The script works around four build problems:
* conda-forge's Geant4 exports `-std=c++17`, which overrides DD4hep's C++20 and breaks podio 1.x
  headers.
* DD4hep refuses that Geant4; the check is downgraded to a warning.
* `solid_dd4hep` asks for ACTS, which its sources don't use, and misses `<fmt/format.h>`.
* ROOT's CMake config requires the exact Python version it was built with.

```bash
source scripts/geant4/env.sh                       # after setup_stack.sh, in the run directory
build/src/solid_ddvcs_gen -n 2000 -s 2026 -o ddvcs.hepmc
python3 scripts/run_ddsim_ddvcs.py --compact $DETECTOR_PATH/solid_ddvcs.xml -i ddvcs.hepmc -n 500 -s 0 -o g4_part0
jana -Pplugins=solid_ddvcs -Psolid:source=edm4hep -Pfield:mode=map -Pfield:map=$SOLID_FIELD_MAP \
     -Ptracking:engine=acts -Pacts:energy_loss=1 -Pacts:ec_match_cm=0 g4_part*.edm4hep.root
```

Geant4 takes about 1 s per event, with 40 s to load the geometry and the field. A 500-event
EDM4hep file is about 100 MB.

## What the first run found and fixed

The compact files had never been loaded before. In the order found:

1. **Compact syntax**: `<limit>` needs a separate `unit` attribute. `DD4hep_TubeSegment` needs a
   `<rotation>` and `DD4hep_BoxSegment` reads the material from a `<material>` child.
   `DD4hep_CylinderTracker` does not exist, so the muon cylinders are now sensitive tube segments.
   The disk-tracker readouts needed the `barrel` and `slice` volume-ID fields.
2. **Field map**: `solid_dd4hep` downloads the map from jlabsvn. Here it comes from `solid_gemc` on
   GitHub instead, stored in the loader's hash-named cache layout.
3. **No magnet in the geometry**: every large-angle muon left the world without reaching the
   large-angle muon detector. The CLEO map includes the yoke return field (Bz ≈ −1.3 T at
   r = 2 m, Br ≈ 1.5 T at the endcap), but the geometry had no iron, so muons bent back in air
   instead of crossing about 1 m of steel. `geometry/compact/magnet_ddvcs.xml` now carries the
   coil, cryostat and CLEO iron. It is converted from `solid_gemc` `geometry/magnet_moved` by
   `scripts/gemc_magnet_to_dd4hep.py`, with the octagons as polycones at the flat-face radius.
4. **Muon scintillator not sensitive**: `DD4hep_TubeSegment` needs `sensitive="yes"` on the
   detector element; the `<sensitive>` child only sets the type. Before this fix, 0% of muons at
   17–28° reached the detector; after it, 58% do.
5. **Overlaps** (`checkOverlaps`, 0.1 mm tolerance): the SPD and MRPC planes overlapped the
   `solid_dd4hep` calorimeters by 2–3 cm. They now sit in the actual gaps: LASPD upstream of GEM
   plane 4 (the GEM4–LAEC gap is only 4 mm), MRPC and FASPD between the HGC (ends at 421 cm) and
   the FAEC (425 cm). **Still open, upstream in solid_dd4hep**: its SIDIS-layout LGC and HGC tanks
   overlap the long-endcap nose cone by 8.0 and 2.4 cm, and the LGC mirrors extrude its gas volume.
6. **DDG4 steering**: EDM4hep output through `Geant4Output2EDM4hep`. The Cherenkov readout used
   npsim's `PhotoMultiplierSDAction`, which is not in DD4hep; optical photons are now opt-in
   (`--optical`).
7. **Reconstruction reading Geant4**:
   * The track fit must use the same field as the simulation: with the coil-only analytic model,
     momenta came out 6% high; with `field:mode=map` the remaining bias is −0.6%.
   * The calorimeter cluster threshold was 100 MeV, while a muon deposits a sampled MIP signal of
     36 MeV in the FAEC. The threshold is now `ec:cluster_emin_gev`, default 20 MeV.
   * Cluster positions are energy-weighted shower centroids at their own depth (not the front
     face), and tracks are projected to that depth. Without this, inclined tracks missed their
     cluster by several cm, and the native tracker, which seeds from calorimeter clusters, found
     nothing.

## Results: 2000 events, internal phase-space generator, 1.2×10³⁷, strip readout, background model on

The same generated events go through the fast simulation and through Geant4, both with the CLEO
field map. "Accepted" means the e′ reaches a calorimeter and both muons reach the muon detector
(last forward layer or the large-angle cylinder), with ≥ 4 GEM hits each. The columns
"tracked + identified" give, for accepted events, the fraction of each lepton with a matched track
that also passes the identification extrapolation.

| simulation | tracking | geometric acceptance | ε(reco \| accepted) | e / μ⁺ / μ⁻ tracked + identified | ⟨MM²⟩ [GeV²] | σ(MM²) core [MeV²] | σ₆₈(MM²) [MeV²] | σ(p)/p μ | σ(p)/p e | σ(θ) μ [mrad] | σ(z_v) [cm] | fake tracks |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fast sim | native | 0.142 | 0.94 | — | 0.874 | 60 | 63 | 0.79% | 0.81% | 0.47 | 0.25 | 0.02% |
| Geant4 | native | 0.123 | 0.59 | 0.84 / 0.86 / 0.83 | 0.909 | 88 | 106 | 1.07% | 1.22% | 0.46 | 0.33 | 0.11% |
| fast sim | ACTS (calo match 5 cm) | 0.142 | 0.98 | — | 0.880 | 56 | 61 | 0.80% | 0.82% | 0.47 | 0.26 | 0.05% |
| Geant4 | ACTS (calo match 5 cm) | 0.123 | 0.48 | 0.80 / 0.81 / 0.77 | 0.915 | 88 | 105 | 1.07% | 1.20% | 0.46 | 0.32 | 0.06% |
| Geant4 | ACTS (no calo match) | 0.123 | 0.73 | 0.91 / 0.92 / 0.89 | 0.912 | 83 | 104 | 1.07% | 1.21% | 0.47 | 0.33 | 0.17% |
| fast sim | GNN | 0.142 | 0.98 | — | 0.874 | 58 | 61 | 0.80% | 0.82% | 0.47 | 0.25 | 0.02% |
| Geant4 | GNN | 0.123 | 0.68 | 0.87 / 0.91 / 0.88 | 0.914 | 86 | 99 | 1.04% | 1.20% | 0.47 | 0.33 | 0.08% |
| Geant4 | GNN (calo match 5 cm) | 0.123 | 0.48 | 0.79 / 0.80 / 0.77 | 0.914 | 88 | 98 | 1.05% | 1.21% | 0.46 | 0.32 | 0.03% |

**What Geant4 changes:**

* **Acceptance 14.2% → 12.3%.** Muons below about 3 GeV range out in the 108 cm of iron in front of
  the last forward muon layer, which the fast simulation's p > 2 GeV threshold does not capture.
  Large-angle muons must also cross the yoke.
* **Resolution σ(MM²) 56–60 → 83–88 MeV²** (σ(p)/p 0.80% → 1.05%, σ(z_v) 0.25 → 0.33 cm). The
  sources are energy-loss straggling, electron bremsstrahlung, and the full material (GEM modules,
  Cherenkov gas, windows) that the fast simulation lumps into 0.7% X₀ per plane. The MM² peak
  sits at 0.91 GeV², 0.03 GeV² above Mp², because the energy lost in the target, Cherenkovs and
  windows is not corrected.
* **Reconstruction efficiency for accepted events 98% → 68–73%** (ACTS or GNN without the
  calorimeter match; native 59%). Per lepton, the losses are:
  * **Electron tracking:** 91% (ACTS), 87% (GNN). Bremsstrahlung kinks need a Gaussian-sum filter.
  * **Muon identification:** muons are tracked at 99–100%, but 8–12% fail the identification
    extrapolation, which ignores the iron and uses a 2 GeV range cut.
  * **Calorimeter match:** requiring it costs a further 10% per lepton, because only about 60% of
    muons leave a MIP cluster above 20 MeV. Use `-Pacts:ec_match_cm=0` on Geant4 input; the match
    was tuned on the fast simulation.
* **Fake tracks stay below 0.2%** at this luminosity.

**Next steps, in order of impact:**
1. Muon identification from the simulated muon-detector hits matched to the track, instead of
   extrapolation plus a momentum cut.
2. Electron fitting with the ACTS Gaussian-sum filter.
3. Target, Cherenkov and window material in the ACTS geometry, or an energy-loss correction at
   the vertex, to centre the MM² peak.
4. A DDVCS layout of the LGC/HGC in `solid_dd4hep`.
5. Retrain the GNN on Geant4 hits.
6. Feed the background model from a Geant4 beam-on-target run instead of the parametrised table.
