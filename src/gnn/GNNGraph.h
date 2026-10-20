// Hit graph for GNN track finding in the SoLID GEM tracker (Exa.TrkX-style edge classification).
//
// Nodes are GEM hits. Directed edges go from a hit on plane i to a hit on plane j = i+1 or i+2
// (the second allows one missed plane), outward in radius, with geometric cuts that keep the
// true segments of tracks from the target:
//   * delta phi between the two hits (bending in the solenoid),
//   * z0: the straight r-z line through the two hits, extrapolated to r = 0, near the target,
//   * |delta phi| / delta r (a curvature proxy, i.e. a minimum pT).
// The default cuts keep 99.9% of the e' / mu+- segments of DDVCS events and ~80% of those of the
// soft, strongly bending recoil proton (fast simulation, 1.2e37).
// Node and edge features are computed here, for both the training-data export
// (GNNGraphDump_processor) and the inference (GNNTrackFinder_factory), so that training sees
// exactly what the reconstruction computes. Feature layout (all in cm or rad, scaled to O(1)):
//   node: x/100, y/100, r/100, z/100
//   edge: dr/100, dphi, dz/100, (z0 - z_target)/100, dphi/dr*100 [rad/m], plane gap - 1
#pragma once

#include <cstdint>
#include <vector>
#include "datamodel/SoLIDDataModel.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

struct GNNGraphConfig {
  int max_plane_gap = 2;          // edges to the next plane and the one after
  double max_dphi = 0.25;         // rad
  double z0_window_cm = 50.0;     // |z0 - z_target| (straight r-z extrapolation of a bending track overshoots)
  double max_dphi_per_dr = 1.2;   // rad/m
  int max_edges = 200000;         // safety bound for saturated planes
};

struct GNNGraph {
  static constexpr int kNodeFeatures = 4;
  static constexpr int kEdgeFeatures = 6;
  std::vector<int> hit;                 // node -> index in the GEMHit collection
  std::vector<float> node_features;     // N x kNodeFeatures, row-major
  std::vector<std::int64_t> src, dst;   // edges (node indices), src on the inner plane
  std::vector<float> edge_features;     // E x kEdgeFeatures, row-major
  bool truncated = false;               // max_edges reached
  std::size_t numNodes() const { return hit.size(); }
  std::size_t numEdges() const { return src.size(); }
};

GNNGraph buildGNNGraph(const std::vector<const GEMHit*>& hits, const SoLIDGeometry_service& geo,
                       const GNNGraphConfig& cfg);

/// Truth label of each edge: both hits from the same particle and no hit of that particle on a
/// plane in between (so an edge skipping a plane is true only where the particle missed it).
std::vector<int> gnnEdgeTruth(const GNNGraph& g, const std::vector<const GEMHit*>& hits);

}  // namespace solid
