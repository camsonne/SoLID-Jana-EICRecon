#include "GNNGraphDump_processor.h"
#include <JANA/JEvent.h>
#include <JANA/JLogger.h>
#include <iomanip>
#include <map>
#include <set>

namespace solid {

GNNGraphDump_processor::GNNGraphDump_processor() {
  SetTypeName(NAME_OF_THIS);
  SetCallbackStyle(CallbackStyle::ExpertMode);
  m_hits_in.SetDatabundleName("GEMHits");
}

void GNNGraphDump_processor::Init() {
  m_nodes.open(m_prefix() + "_nodes.csv");
  m_edges.open(m_prefix() + "_edges.csv");
  m_nodes << "event,node,plane";
  for (int i = 0; i < GNNGraph::kNodeFeatures; ++i) m_nodes << ",f" << i;
  m_nodes << ",mc\n";
  m_edges << "event,src,dst";
  for (int i = 0; i < GNNGraph::kEdgeFeatures; ++i) m_edges << ",f" << i;
  m_edges << ",label\n";
  m_nodes << std::setprecision(6);
  m_edges << std::setprecision(6);
}

void GNNGraphDump_processor::ProcessSequential(const JEvent& event) {
  const auto& hits = m_hits_in();
  const GNNGraph g = buildGNNGraph(hits, m_geo(), m_graph.config());
  const std::vector<int> label = gnnEdgeTruth(g, hits);
  const auto ev = event.GetEventNumber();
  ++m_events;
  if (g.truncated) ++m_truncated;
  for (size_t n = 0; n < g.numNodes(); ++n) {
    const GEMHit& h = *hits[g.hit[n]];
    m_nodes << ev << "," << n << "," << h.plane;
    for (int i = 0; i < GNNGraph::kNodeFeatures; ++i) m_nodes << "," << g.node_features[n * GNNGraph::kNodeFeatures + i];
    m_nodes << "," << h.mc_index << "\n";
  }
  for (size_t e = 0; e < g.numEdges(); ++e) {
    m_edges << ev << "," << g.src[e] << "," << g.dst[e];
    for (int i = 0; i < GNNGraph::kEdgeFeatures; ++i) m_edges << "," << g.edge_features[e * GNNGraph::kEdgeFeatures + i];
    m_edges << "," << label[e] << "\n";
    m_true_edges += label[e];
  }
  m_edges_total += static_cast<long>(g.numEdges());

  // true segments: consecutive hits of one particle (no hit of it on a plane in between)
  std::map<int, std::set<int>> planes_of;
  for (const auto* h : hits)
    if (h->mc_index >= 0) planes_of[h->mc_index].insert(h->plane);
  for (const auto& [mc, planes] : planes_of) {
    if (planes.size() < 2) continue;
    m_true_segments += static_cast<long>(planes.size()) - 1;
  }
  for (int l : label) m_true_segments_in_graph += l;
}

void GNNGraphDump_processor::Finish() {
  m_nodes.close();
  m_edges.close();
  LOG_INFO(GetLogger()) << "GNN graphs: " << m_events << " events, " << double(m_edges_total) / std::max(1L, m_events)
                        << " edges/event, true-edge fraction " << double(m_true_edges) / std::max(1L, m_edges_total)
                        << ", graph-construction efficiency " << double(m_true_segments_in_graph) / std::max(1L, m_true_segments)
                        << " (" << m_true_segments_in_graph << "/" << m_true_segments << " true segments), "
                        << m_truncated << " truncated graphs";
}

}  // namespace solid
