// GNN track finding for the SoLID GEM tracker (tracking:engine=gnn).
//
//   1. Hit graph (GNNGraph): GEM hits as nodes, candidate segments between neighbouring planes
//      (or skipping one plane) as edges, with geometric cuts.
//   2. Edge classification by a graph neural network (Interaction Network trained with
//      scripts/gnn/train_gnn.py on fast-simulation events) run through ONNX Runtime.
//   3. Track candidates from the edges above gnn:edge_threshold (gnn:track_builder):
//      walkthrough (default): all paths from hits without an incoming edge, ranked by length and
//      summed edge score, fitted and accepted greedily without shared hits;
//      mutual_best: chains of edges that are the best outgoing edge of their inner hit and the best
//      incoming edge of their outer hit (faster, but a background hit next to a true hit that wins
//      one comparison cuts the track).
//   4. Fit of each chain with at least gnn:min_hits hits by the native global chi2 fitter (RK4 in the
//      field, multiple-scattering covariance, beam-line constraint), with removal of the worst hit
//      while chi2/ndf is too large.
//   5. The same Track object as the other engines: truth matching and calorimeter / muon-detector
//      flags for PID; optional calorimeter match (gnn:ec_match_cm).
#pragma once

#include <JANA/JFactory.h>
#include <memory>
#include <string>
#include "datamodel/SoLIDDataModel.h"
#include "gnn/GNNEdgeClassifier.h"
#include "gnn/GNNGraphDump_processor.h"
#include "reco/TrackFitter.h"
#include "services/SoLIDField_service.h"
#include "services/SoLIDGeometry_service.h"

#ifndef SOLID_GNN_DEFAULT_MODEL
#define SOLID_GNN_DEFAULT_MODEL "solid_gem_in.onnx"
#endif

namespace solid {

class GNNTrackFinder_factory : public JFactory {
public:
  GNNTrackFinder_factory();
  void Init() override;
  void Process(const JEvent& event) override;
  void Finish() override;

private:
  bool fitCandidate(const std::vector<int>& hit_indices, const std::vector<const GEMHit*>& hits, FitResult& fit,
                    std::vector<int>& kept);
  bool matchesCluster(const Track& trk);

  Input<GEMHit> m_hits_in{this};
  Input<ECCluster> m_ec_in{this};
  Output<Track> m_tracks_out{this};
  Service<SoLIDGeometry_service> m_geo{this};
  Service<SoLIDField_service> m_field{this};

  GNNGraphParameters m_graph{this};
  Parameter<std::string> m_model{this, "gnn:model", SOLID_GNN_DEFAULT_MODEL, "ONNX model of the edge classifier"};
  Parameter<double> m_threshold{this, "gnn:edge_threshold", 0.02,
      "Minimum edge score for the track building (low: the walkthrough and the fit reject wrong paths)"};
  Parameter<std::string> m_builder{this, "gnn:track_builder", "walkthrough",
      "Track building from the scored graph: walkthrough (all paths, fitted and accepted greedily) | mutual_best"};
  Parameter<int> m_max_paths{this, "gnn:max_paths_per_start", 64, "Walkthrough: maximum paths followed from one hit"};
  Parameter<int> m_min_hits{this, "gnn:min_hits", 4, "Minimum number of GEM hits per track"};
  Parameter<double> m_chi2ndf_max{this, "gnn:chi2ndf_max", 10.0, "Maximum chi2/ndf of accepted tracks"};
  Parameter<double> m_b_eff{this, "gnn:b_eff_tesla", 1.0, "Effective uniform field for the analytic fit seed [T]"};
  Parameter<double> m_zv_sigma{this, "gnn:zv_prior_cm", 6.0, "Vertex-z prior width of the fit [cm]"};
  Parameter<double> m_ec_match{this, "gnn:ec_match_cm", 0.0,
      "Keep only tracks pointing to a FAEC/LAEC cluster within this distance [cm] (<= 0 disables)"};
  Parameter<double> m_muon_pmin{this, "solid:muon_pmin", 2.0,
      "Minimum momentum for a muon to reach the last FA muon scintillator [GeV]"};

  std::unique_ptr<GNNEdgeClassifier> m_classifier;
  std::unique_ptr<TrackFitter> m_fitter;
  long m_events = 0, m_nodes = 0, m_edges = 0, m_kept_edges = 0, m_candidates = 0, m_tracks = 0;
};

}  // namespace solid
