// Training-data export for the GNN edge classifier (-Pgnn:dump=<prefix>): writes the hit graphs
// built by buildGNNGraph, with their features and truth labels, to <prefix>_nodes.csv and
// <prefix>_edges.csv (scripts/gnn/train_gnn.py reads them), and reports the graph-construction
// efficiency: the fraction of true track segments (consecutive hits of one particle) that are
// edges of the graph.
#pragma once

#include <JANA/JEventProcessor.h>
#include <fstream>
#include <string>
#include "datamodel/SoLIDDataModel.h"
#include "gnn/GNNGraph.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

/// gnn:* graph-construction parameters shared by the dump processor and the GNN track finder.
struct GNNGraphParameters {
  explicit GNNGraphParameters(jana::components::JComponent* c)
      : max_plane_gap(c, "gnn:max_plane_gap", 2, "Edges from plane i to planes i+1 .. i+gap"),
        max_dphi(c, "gnn:max_dphi", 0.25, "Maximum azimuthal difference between connected hits [rad]"),
        z0_window(c, "gnn:z0_window_cm", 50.0, "Maximum |z0 - z_target| of the r-z line through the two hits [cm]"),
        max_curv(c, "gnn:max_dphi_per_dr", 1.2, "Maximum |dphi/dr| (minimum pT) [rad/m]"),
        max_edges(c, "gnn:max_edges", 200000, "Maximum number of graph edges per event") {}
  GNNGraphConfig config() {
    GNNGraphConfig cfg;
    cfg.max_plane_gap = max_plane_gap();
    cfg.max_dphi = max_dphi();
    cfg.z0_window_cm = z0_window();
    cfg.max_dphi_per_dr = max_curv();
    cfg.max_edges = max_edges();
    return cfg;
  }
  jana::components::JComponent::Parameter<int> max_plane_gap;
  jana::components::JComponent::Parameter<double> max_dphi, z0_window, max_curv;
  jana::components::JComponent::Parameter<int> max_edges;
};

class GNNGraphDump_processor : public JEventProcessor {
public:
  GNNGraphDump_processor();
  void Init() override;
  void ProcessSequential(const JEvent& event) override;
  void Finish() override;

private:
  Input<GEMHit> m_hits_in{this};
  Service<SoLIDGeometry_service> m_geo{this};
  Parameter<std::string> m_prefix{this, "gnn:dump", "gnn_graphs", "Output prefix of the GNN training data"};
  GNNGraphParameters m_graph{this};

  std::ofstream m_nodes, m_edges;
  long m_events = 0, m_edges_total = 0, m_true_edges = 0, m_true_segments = 0, m_true_segments_in_graph = 0;
  long m_truncated = 0;
};

}  // namespace solid
