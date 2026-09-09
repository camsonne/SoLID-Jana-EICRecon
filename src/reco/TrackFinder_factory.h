// Track finding + fitting for the SoLID GEM tracker.
//
// Seeding follows the SoLID strategy: the beam-line vertex constraint plus a calorimeter
// cluster (FAEC or LAEC, ~1 cm position resolution) define a one-parameter family of
// trajectories; a hit on one GEM plane fixes the seed. Hits are then collected on all planes
// inside a road, the best candidates are fitted with the global chi2 fit (RK4 propagation,
// multiple-scattering covariance), re-collected with a tight window and refitted with
// outlier rejection. Without calorimeter clusters a GEM-pair seeding is used instead.
#pragma once

#include <JANA/JFactory.h>
#include <memory>
#include "datamodel/SoLIDDataModel.h"
#include "reco/TrackFitter.h"
#include "services/SoLIDField_service.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

class TrackFinder_factory : public JFactory {
public:
  TrackFinder_factory();
  void Init() override;
  void Process(const JEvent& event) override;

private:
  struct Candidate {
    std::vector<FitHit> hits;   // GEM hits (plane>0) and optionally the calorimeter point (plane 0)
    std::vector<int> indices;   // GEMHit indices for the GEM hits
    TrackParams seed;
    double seed_score = 0;      // sum of seed-prediction residuals
    FitResult fit;
  };

  /// Progressive extension: starting from the fitted seed (one or two GEM hits, optional calorimeter
  /// point, vertex prior) add the nearest hit plane by plane inside a window derived from the
  /// predicted position uncertainty, refitting after each addition.
  bool extend(const FitResult& seedfit, const std::vector<FitHit>& seed_hits, const std::vector<std::vector<int>>& by_plane,
              const std::vector<const GEMHit*>& hits, const std::vector<char>& used, const ECCluster* ec, Candidate& cand);
  /// Predicted position and (1D) uncertainty at plane `plane` for a fitted track.
  bool predictAt(const TrackParams& par, const std::array<double, 16>& cov, const std::vector<FitHit>& gem_hits, int plane,
                 Vec3& pos, double& sigma);
  bool refine(Candidate& cand, const std::vector<std::vector<int>>& by_plane, const std::vector<const GEMHit*>& hits,
              const std::vector<char>& used, const ECCluster* ec);
  Track* makeTrack(const Candidate& cand, const std::vector<const GEMHit*>& hits);
  void setAcceptanceFlags(Track& trk);

  Input<GEMHit> m_hits_in{this};
  Input<ECCluster> m_ec_in{this};
  Output<Track> m_tracks_out{this};

  Service<SoLIDGeometry_service> m_geo{this};
  Service<SoLIDField_service> m_field{this};

  Parameter<int> m_min_hits{this, "tracking:min_hits", 4, "Minimum number of GEM hits per track"};
  Parameter<double> m_seed_window{this, "tracking:seed_window_cm", 5.0, "Maximum road half-width around the seed prediction [cm]"};
  Parameter<double> m_window_margin{this, "tracking:window_margin_cm", 0.15, "Additive safety margin of the search window [cm]"};
  Parameter<double> m_tight_window{this, "tracking:tight_window_cm", 0.3, "Minimum half-width of the refit window [cm]"};
  Parameter<double> m_nsigma{this, "tracking:window_nsigma", 5.0, "Refit window in units of the predicted position error"};
  Parameter<double> m_chi2ndf_max{this, "tracking:chi2ndf_max", 10.0, "Maximum chi2/ndf of accepted tracks"};
  Parameter<double> m_b_eff{this, "tracking:b_eff_tesla", 1.0, "Effective uniform field for the analytic seed [T]"};
  Parameter<double> m_zv_sigma{this, "tracking:zv_prior_cm", 6.0, "Vertex-z prior width (target length / sqrt(12) inflated) [cm]"};
  Parameter<double> m_pmin{this, "tracking:pmin", 0.3, "Minimum seed momentum [GeV]"};
  Parameter<int> m_max_hits_plane{this, "tracking:max_hits_per_plane", 600, "Planes with more hits than this (saturated readout) are not used for seeding"};
  Parameter<int> m_max_seeds{this, "tracking:max_seeds_per_cluster", 40, "Seed-plane hits tried per calorimeter cluster (closest to the cluster-target straight line first)"};
  Parameter<int> m_max_fits{this, "tracking:max_fits_per_seed_cluster", 6, "Number of best candidates fitted per calorimeter cluster"};
  Parameter<bool> m_use_ec{this, "tracking:use_ec_seed", true, "Seed with calorimeter clusters when available"};
  Parameter<double> m_step{this, "tracking:rk_step_cm", 2.0, "RK4 step [cm]"};
  Parameter<double> m_muon_pmin{this, "solid:muon_pmin", 2.0, "Minimum momentum for a muon to reach the last FA muon scintillator [GeV]"};

  std::unique_ptr<TrackFitter> m_fitter;
};

}  // namespace solid
