#include "GNNGraph.h"
#include <cmath>
#include <map>
#include <set>

namespace solid {

GNNGraph buildGNNGraph(const std::vector<const GEMHit*>& hits, const SoLIDGeometry_service& geo,
                       const GNNGraphConfig& cfg) {
  GNNGraph g;
  const int nplanes = static_cast<int>(geo.gemPlanes().size());
  const double ztgt = geo.target().z_center;
  std::vector<std::vector<int>> by_plane(nplanes + 1);   // node indices per plane
  for (size_t i = 0; i < hits.size(); ++i) {
    const GEMHit& h = *hits[i];
    if (h.plane < 1 || h.plane > nplanes) continue;
    const int node = static_cast<int>(g.hit.size());
    g.hit.push_back(static_cast<int>(i));
    g.node_features.insert(g.node_features.end(), {static_cast<float>(h.x / 100.0), static_cast<float>(h.y / 100.0),
                                                   static_cast<float>(h.r() / 100.0), static_cast<float>(h.z / 100.0)});
    by_plane[h.plane].push_back(node);
  }
  for (int pa = 1; pa <= nplanes; ++pa) {
    for (int pb = pa + 1; pb <= std::min(nplanes, pa + cfg.max_plane_gap); ++pb) {
      for (int na : by_plane[pa]) {
        const GEMHit& a = *hits[g.hit[na]];
        const double ra = a.r(), phia = a.phi();
        for (int nb : by_plane[pb]) {
          const GEMHit& b = *hits[g.hit[nb]];
          const double rb = b.r();
          const double dr = rb - ra;
          if (dr <= 0) continue;
          const double dphi = std::remainder(b.phi() - phia, 2 * M_PI);
          if (std::fabs(dphi) > cfg.max_dphi) continue;
          const double dz = b.z - a.z;
          const double z0 = a.z - ra * dz / dr;
          if (std::fabs(z0 - ztgt) > cfg.z0_window_cm) continue;
          const double curv = dphi / dr * 100.0;   // rad/m
          if (std::fabs(curv) > cfg.max_dphi_per_dr) continue;
          if (static_cast<int>(g.src.size()) >= cfg.max_edges) { g.truncated = true; return g; }
          g.src.push_back(na);
          g.dst.push_back(nb);
          g.edge_features.insert(g.edge_features.end(),
                                 {static_cast<float>(dr / 100.0), static_cast<float>(dphi), static_cast<float>(dz / 100.0),
                                  static_cast<float>((z0 - ztgt) / 100.0), static_cast<float>(curv),
                                  static_cast<float>(pb - pa - 1)});
        }
      }
    }
  }
  return g;
}

std::vector<int> gnnEdgeTruth(const GNNGraph& g, const std::vector<const GEMHit*>& hits) {
  std::map<int, std::set<int>> planes_of;   // particle -> planes with a hit
  for (const auto* h : hits)
    if (h->mc_index >= 0) planes_of[h->mc_index].insert(h->plane);
  std::vector<int> label(g.numEdges(), 0);
  for (size_t e = 0; e < g.numEdges(); ++e) {
    const GEMHit& a = *hits[g.hit[g.src[e]]];
    const GEMHit& b = *hits[g.hit[g.dst[e]]];
    if (a.mc_index < 0 || a.mc_index != b.mc_index) continue;
    const auto& pl = planes_of[a.mc_index];
    bool gap_free = true;
    for (int p = a.plane + 1; p < b.plane; ++p)
      if (pl.count(p)) { gap_free = false; break; }
    label[e] = gap_free ? 1 : 0;
  }
  return label;
}

}  // namespace solid
