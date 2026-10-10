// ACTS track finding and fitting for the SoLID GEM tracker (tracking:engine=acts).
//
// Chain, following EICrecon (ePIC) and the ACTS examples:
//   1. GEM hits -> 2D measurements (local x, y) on the ACTS plane surfaces of SoLIDActs_service,
//      and space points for seeding.
//   2. Triplet seeding (ACTS Seeding2: k-d tree, doublet/triplet finders, broad seed filter), with
//      the collision region set to the LH2 target widened by acts:seed_z0_margin_mm.
//   3. Track parameters from each seed (helix in the local field), transported back to the beam
//      line (perigee), and refined by a pre-fit in which the CKF is offered only the three seed
//      measurements. Without it the loose seed covariance lets the CKF pick a wrong hit on a plane
//      where the track has none, and pull the track off.
//   4. Combinatorial Kalman Filter from the beam line through the GEM planes, smoothing, and
//      extrapolation of the smoothed track to the perigee.
//   5. Greedy ambiguity resolution on shared measurements.
//   6. SoLID beam-line constraint (Kalman update d0 = 0 with acts:beamline_sigma_mm) and the
//      multiple scattering in the LH2 target and Al cell added to the angular covariance; then
//      the same Track object as the native tracker (q/p, theta, phi, z_v and covariance), truth
//      matching and the calorimeter / muon-detector extrapolation flags used for PID.
//   7. Calorimeter match (acts:ec_match_cm): the track must point to a FAEC or LAEC cluster.
#pragma once

#include <JANA/JFactory.h>
#include <memory>
#include "acts/SoLIDActs_service.h"
#include "datamodel/SoLIDDataModel.h"
#include "reco/TrackPropagator.h"
#include "services/SoLIDField_service.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

struct ActsTrackingImpl;   // ACTS objects (propagators, CKF), defined in the .cc

class ActsTrackFinder_factory : public JFactory {
public:
  ActsTrackFinder_factory();
  ~ActsTrackFinder_factory() override;
  void Init() override;
  void Process(const JEvent& event) override;
  void Finish() override;

private:
  bool matchesCluster(const Track& trk);

public:

  struct Stats {
    long events = 0, hits = 0, seeds = 0, ckf_tracks = 0, tracks = 0, failed_seeds = 0;
  };

private:
  Input<GEMHit> m_hits_in{this};
  Input<ECCluster> m_ec_in{this};
  Output<Track> m_tracks_out{this};

  Service<SoLIDGeometry_service> m_geo{this};
  Service<SoLIDField_service> m_field{this};
  Service<SoLIDActs_service> m_acts{this};

  // seeding (ACTS units: mm, GeV)
  Parameter<double> m_seed_minpt{this, "acts:seed_min_pt", 0.2, "Minimum seed pT [GeV]"};
  Parameter<double> m_seed_impact{this, "acts:seed_impact_max_mm", 150.0,
      "Maximum seed impact parameter [mm] (loose: the seed helix uses the field inside the solenoid, "
      "the track starts in the weak fringe field at the target)"};
  Parameter<double> m_seed_z0_margin{this, "acts:seed_z0_margin_mm", 300.0,
      "Margin around the target for the seed z0 (straight r-z extrapolation of a helix overshoots upstream) [mm]"};
  Parameter<double> m_seed_dphi{this, "acts:seed_delta_phi_max", 0.6, "Maximum phi difference between seed space points [rad]"};
  Parameter<double> m_seed_drmin{this, "acts:seed_delta_r_min_mm", 10.0, "Minimum radial distance between seed space points [mm]"};
  Parameter<double> m_seed_drmax{this, "acts:seed_delta_r_max_mm", 1200.0, "Maximum radial distance between seed space points [mm]"};
  Parameter<double> m_seed_cotmax{this, "acts:seed_cot_theta_max", 9.5, "Maximum cot(theta) of a seed (theta > 6 deg)"};
  Parameter<int> m_seed_max_per_middle{this, "acts:seed_max_per_middle_sp", 5, "Seeds kept per middle space point"};
  // CKF
  Parameter<double> m_chi2_cut{this, "acts:ckf_chi2_cut", 15.0, "CKF measurement chi2 cut"};
  Parameter<int> m_max_branches{this, "acts:ckf_max_measurements_per_surface", 1,
      "Measurements followed per surface (>1 lets the CKF branch)"};
  Parameter<int> m_min_meas{this, "acts:min_measurements", 4, "Minimum number of GEM measurements per track"};
  Parameter<int> m_max_shared{this, "acts:max_shared_hits", 1, "Maximum measurements shared with a better track"};
  Parameter<double> m_chi2ndf_max{this, "acts:chi2ndf_max", 10.0, "Maximum chi2/ndf of accepted tracks"};
  Parameter<bool> m_eloss{this, "acts:energy_loss", false,
      "Energy-loss correction in the fit (off: the fast simulation has no energy loss; on for Geant4 input)"};
  // SoLID vertex treatment
  Parameter<double> m_beamline_sigma{this, "acts:beamline_sigma_mm", 0.1,
      "Beam-line constraint on the transverse vertex position [mm] (<= 0 disables it)"};
  Parameter<bool> m_target_ms{this, "acts:target_scattering", true,
      "Add the multiple scattering in the LH2 target and Al cell to the vertex angular covariance"};
  Parameter<double> m_ec_match{this, "acts:ec_match_cm", 5.0,
      "Keep only tracks pointing to a FAEC/LAEC cluster within this distance [cm] (<= 0 disables). "
      "The native tracker gets the same effect from its calorimeter seeding; without it, combinations "
      "of background hits form fake tracks at high occupancy"};
  Parameter<std::string> m_log_level{this, "acts:log_level", "WARNING", "ACTS log level: VERBOSE DEBUG INFO WARNING ERROR"};
  Parameter<double> m_muon_pmin{this, "solid:muon_pmin", 2.0,
      "Minimum momentum for a muon to reach the last FA muon scintillator [GeV]"};

  std::unique_ptr<ActsTrackingImpl> m_impl;
  std::unique_ptr<TrackPropagator> m_prop;   // for the PID extrapolation (same as the native tracker)
  Stats m_stats;
};

}  // namespace solid
