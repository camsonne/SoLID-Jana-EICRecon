# ACTS tracking for SoLID DDVCS

`tracking:engine=acts` (default when built with `-DUSE_ACTS=ON`) runs ACTS 46 (the eic-shell
version, 46.8.1) on the SoLID GEM tracker. The code is in `src/acts/`.

## Chain

1. **Geometry** (`SoLIDActs_service`): ACTS Gen1 tracking geometry built from
   `SoLIDGeometry_service`, the same numbers the fast simulation uses. It is one air-filled
   cylindrical volume from upstream of the target to behind GEM plane 6. Each GEM plane is a
   sensitive planar ring: `PlaneSurface` with `EllipseBounds`, widened by `acts:surface_margin_mm`
   so that tracks predicted just outside an edge are still followed. Each ring carries the plane
   material. The field is `SoLIDField_service` (analytic coil model or field map) wrapped as an
   `Acts::MagneticFieldProvider`. When the DD4hep compact files can be loaded, this service is
   the place to switch to the ACTS DD4hep plugin, as EICrecon does.
2. **Measurements**: each GEM hit becomes a 2D (local x, y) measurement and a seeding space point.
3. **Seeding**: ACTS `Seeding2` triplet seeding (k-d tree, doublet and triplet finders, broad seed
   filter), the API EICrecon uses with ACTS 46. The collision region is the LH2 target widened by
   `acts:seed_z0_margin_mm` (300 mm), because the straight r–z extrapolation of a bending track
   overshoots upstream. The impact-parameter cut is loose (150 mm): the seed helix uses the
   field inside the solenoid, while the track starts in the weak fringe field at the target.
4. **Seed pre-fit**: seed parameters are transported to the beam line and refined by a CKF pass
   that is offered only the three seed measurements. Without this, the loose seed covariance
   lets the CKF accept a wrong hit on a plane where the track has none.
5. **CKF** from the beam line through the GEM planes, then smoothing and extrapolation to the
   perigee.
6. **Greedy ambiguity resolution** on shared measurements.
7. **SoLID vertex treatment**: a Kalman update with the beam-line constraint d0 = 0
   (`acts:beamline_sigma_mm`, 0.1 mm), then the multiple scattering in the LH2 target and Al
   cell added to the angular covariance. The scattering uses the same model as the fast
   simulation. It happens at the vertex, so it is not part of the ACTS geometry.
8. **Calorimeter match** (`acts:ec_match_cm`, 5 cm): the track must point to a FAEC or LAEC
   cluster. The native tracker gets the same effect from seeding on calorimeter clusters.
   Without this cut, combinations of background hits form fake tracks at high occupancy.
   Tighter windows trade efficiency for purity: at 10³⁹ with pixel readout, 3 cm gives 1.5%
   fakes and 4% efficiency loss, and 2 cm gives 0.8% fakes and 30% efficiency loss.

The output is the same `Track` object as the native tracker, so `DDVCSReco_factory` and the
analysis work unchanged with either engine. Energy loss in the fit is off by default
(`acts:energy_loss`) because the fast simulation has none; turn it on for Geant4 input.

## Native vs ACTS on the same events (fast simulation, placeholder background table)

| L [cm⁻²s⁻¹] | readout | engine | events | ε(reco \| geom. acc) | ε(reco \| ≥4 GEM hits) | σ(MM²) core [MeV²] | σ(p)/p μ | σ(θ) μ [mrad] | σ(z_v) [cm] | tracks w/ bkg hit | fake tracks | median χ²/ndf | time / event [ms, 4 threads] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.2e+37 | pixel | native | 4000 | 0.825 | 0.998 | 58 | 0.76% | 0.45 | 0.24 | 0.000 | 0.0000 | 0.90 | 8.8 |
| 1.2e+37 | pixel | acts | 4000 | 0.822 | 0.994 | 58 | 0.77% | 0.45 | 0.24 | 0.000 | 0.0000 | 0.78 | 2.6 |
| 1.0e+38 | pixel | native | 4000 | 0.784 | 0.981 | 56 | 0.76% | 0.44 | 0.25 | 0.001 | 0.0000 | 0.89 | 10.8 |
| 1.0e+38 | pixel | acts | 4000 | 0.787 | 0.985 | 56 | 0.77% | 0.45 | 0.25 | 0.001 | 0.0000 | 0.78 | 2.6 |
| 1.0e+39 | pixel | native | 4000 | 0.766 | 0.940 | 59 | 0.77% | 0.44 | 0.24 | 0.023 | 0.0051 | 0.89 | 80.8 |
| 1.0e+39 | pixel | acts | 4000 | 0.793 | 0.975 | 58 | 0.78% | 0.44 | 0.24 | 0.047 | 0.0342 | 0.80 | 13.0 |
| 1.2e+37 | strip | native | 4000 | 0.614 | 0.908 | 52 | 0.76% | 0.45 | 0.25 | 0.006 | 0.0001 | 0.89 | 22.5 |
| 1.2e+37 | strip | acts | 4000 | 0.668 | 0.988 | 53 | 0.77% | 0.45 | 0.25 | 0.006 | 0.0001 | 0.77 | 3.1 |
| 1.0e+38 | strip | native | 1000 | 0.043 | 0.312 | — | 0.76% | 0.70 | 0.28 | 0.549 | 0.4414 | 2.46 | 525.0 |
| 1.0e+38 | strip | acts | 1000 | 0.029 | 0.250 | — | 0.77% | 0.55 | 0.25 | 0.895 | 0.8828 | 2.24 | 208.9 |
| 1.0e+39 | strip | native | 200 | 0.000 | 0.000 | — | — | — | — | 1.000 | 1.0000 | 5.51 | 557.1 |
| 1.0e+39 | strip | acts | 200 | 0.000 | 0.000 | — | — | — | — | 1.000 | 1.0000 | 2.74 | 1132.6 |

Event samples are identical between the two engines at each point (same seed). Times are wall
clock per event with 4 threads and include the fast simulation.

* **Resolution is the same with both engines**: σ(MM²) 56–59 MeV², σ(p)/p 0.77%, σ(θ) 0.45 mrad,
  σ(z_v) 0.24 cm.
* **Efficiency**: ACTS equals the native tracker at 1.2×10³⁷ and 10³⁸ with pixel readout. It is
  better at 10³⁹ pixel (97.5% vs 94.0% of events with three leptons of at least 4 GEM hits) and
  with strip readout at 1.2×10³⁷ (98.8% vs 90.8%).
* **Speed**: ACTS is 3–7× faster where the readout is usable.
* **Fakes**: at 10³⁹ pixel, ACTS keeps more tracks containing background hits (3.4% fakes vs
  0.5%) despite the calorimeter match. They don't affect the DDVCS selection with truth PID;
  with realistic PID, see the window trade-off above.
* **Saturated strip readout** (10³⁸ and above) is unusable with either engine, as before.
