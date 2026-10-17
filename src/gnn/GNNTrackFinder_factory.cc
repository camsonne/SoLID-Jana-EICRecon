#include "GNNTrackFinder_factory.h"
#include <JANA/JEvent.h>
#include <JANA/JLogger.h>
#include <algorithm>
#include <cmath>
#include "reco/TrackUtils.h"

namespace solid {

GNNTrackFinder_factory::GNNTrackFinder_factory() {
  SetTypeName(NAME_OF_THIS);
  m_hits_in.SetDatabundleName("GEMHits");
  m_ec_in.SetDatabundleName("ECClusters");
  m_ec_in.SetOptional(true);
  m_tracks_out.SetShortName("Tracks");
}

void GNNTrackFinder_factory::Init() {
  m_classifier = std::make_unique<GNNEdgeClassifier>(m_model());
  m_fitter = std::make_unique<TrackFitter>(&m_geo(), &m_field(), 2.0);
  LOG_INFO(GetLogger()) << "GNN tracking: model " << m_model() << ", edge threshold " << m_threshold() << ", min "
                        << m_min_hits() << " hits";
}

namespace {
FitHit toFitHit(const GEMHit& h, int index) {
  FitHit f;
  f.plane = h.plane; f.x = h.x; f.y = h.y; f.z = h.z; f.sx = h.sigma_x; f.sy = h.sigma_y; f.index = index;
  return f;
}
}  // namespace

// Fit a chain of hits (ordered by plane); drop the hit with the largest pull while chi2/ndf is
// too large and enough hits remain.
bool GNNTrackFinder_factory::fitCandidate(const std::vector<int>& hit_indices, const std::vector<const GEMHit*>& hits,
                                          FitResult& fit, std::vector<int>& kept) {
  const double ztgt = m_geo().target().z_center;
  kept = hit_indices;
  while (static_cast<int>(kept.size()) >= m_min_hits()) {
    std::vector<FitHit> fh;
    for (int idx : kept) fh.push_back(toFitHit(*hits[idx], idx));
    TrackParams seed;
    if (!m_fitter->seedFromTwoHits(fh.front(), fh.back(), m_b_eff(), seed)) return false;
    // the effective-field seed can misassign the charge of stiff tracks: keep the better fit
    FitResult best;
    for (int sign : {+1, -1}) {
      TrackParams s = seed;
      if (sign < 0) s.qoverp = -s.qoverp;
      FitResult f = m_fitter->fit(s, fh, ztgt, m_zv_sigma(), 10);
      if (f.ok && f.ndf > 0 && (!best.ok || f.chi2 < best.chi2)) best = f;
    }
    if (!best.ok) return false;
    fit = best;
    if (fit.chi2 / fit.ndf <= m_chi2ndf_max()) return true;
    if (static_cast<int>(kept.size()) == m_min_hits()) return false;
    // remove the hit with the largest pull
    int worst = -1;
    double wp = -1;
    for (size_t i = 0; i < kept.size(); ++i) {
      const double p = std::hypot(fit.pulls[2 * i], fit.pulls[2 * i + 1]);
      if (p > wp) { wp = p; worst = static_cast<int>(i); }
    }
    kept.erase(kept.begin() + worst);
  }
  return false;
}

bool GNNTrackFinder_factory::matchesCluster(const Track& trk) {
  const auto& geo = m_geo();
  for (int calo = 0; calo < 2; ++calo) {
    const auto& ec = calo == 0 ? geo.faec() : geo.laec();
    TrackState s;
    s.pos = {0, 0, trk.params.zv};
    s.mom = trk.params.momentum();
    s.charge = trk.params.charge();
    if (!m_fitter->propagator().propagateToZ(s, ec.z)) continue;
    for (const auto* c : m_ec_in()) {
      if (c->calo != calo) continue;
      if (std::hypot(c->x - s.pos.x, c->y - s.pos.y) < m_ec_match()) return true;
    }
  }
  return false;
}

void GNNTrackFinder_factory::Process(const JEvent& event) {
  const auto& hits = m_hits_in();
  std::vector<Track*> out;
  ++m_events;
  const GNNGraph g = buildGNNGraph(hits, m_geo(), m_graph.config());
  m_nodes += static_cast<long>(g.numNodes());
  m_edges += static_cast<long>(g.numEdges());
  if (g.numEdges() == 0) { m_tracks_out() = out; return; }
  const std::vector<float> score = m_classifier->score(g);

  const size_t N = g.numNodes();
  std::vector<std::vector<int>> candidates;   // hit-index lists ordered by plane
  if (m_builder() == "mutual_best") {
    // edges that are the best outgoing edge of their inner hit and the best incoming edge of their
    // outer hit form disjoint chains
    std::vector<int> best_out(N, -1), best_in(N, -1);
    for (size_t e = 0; e < g.numEdges(); ++e) {
      if (score[e] < m_threshold()) continue;
      const auto a = g.src[e], b = g.dst[e];
      if (best_out[a] < 0 || score[e] > score[best_out[a]]) best_out[a] = static_cast<int>(e);
      if (best_in[b] < 0 || score[e] > score[best_in[b]]) best_in[b] = static_cast<int>(e);
    }
    std::vector<int> next(N, -1);
    std::vector<char> has_prev(N, 0);
    for (size_t a = 0; a < N; ++a) {
      const int e = best_out[a];
      if (e < 0 || best_in[g.dst[e]] != e) continue;
      next[a] = static_cast<int>(g.dst[e]);
      has_prev[g.dst[e]] = 1;
      ++m_kept_edges;
    }
    for (size_t a = 0; a < N; ++a) {
      if (has_prev[a] || next[a] < 0) continue;
      std::vector<int> chain;
      for (int n = static_cast<int>(a); n >= 0; n = next[n]) chain.push_back(g.hit[n]);
      if (static_cast<int>(chain.size()) >= m_min_hits()) candidates.push_back(chain);
    }
  } else {
    // walkthrough: every path along edges above threshold, from hits without an incoming
    // above-threshold edge; ranked by length and mean edge score, then fitted and accepted
    // greedily without shared hits. Unlike mutual-best chains, a background hit that wins one
    // comparison next to a true hit cannot cut the track.
    std::vector<std::vector<std::pair<float, int>>> adj(N);   // (score, dst node)
    std::vector<char> has_in(N, 0);
    for (size_t e = 0; e < g.numEdges(); ++e) {
      if (score[e] < m_threshold()) continue;
      adj[g.src[e]].emplace_back(score[e], static_cast<int>(g.dst[e]));
      has_in[g.dst[e]] = 1;
      ++m_kept_edges;
    }
    for (auto& v : adj) std::sort(v.begin(), v.end(), [](const auto& x, const auto& y) { return x.first > y.first; });
    struct Path { std::vector<int> nodes; double sum = 0; };
    std::vector<Path> paths;
    const int max_paths = m_max_paths();
    for (size_t a = 0; a < N; ++a) {
      if (has_in[a] || adj[a].empty()) continue;
      int found = 0;
      std::vector<Path> stack{{{static_cast<int>(a)}, 0.0}};
      while (!stack.empty() && found < max_paths) {
        Path p = std::move(stack.back());
        stack.pop_back();
        const auto& out_edges = adj[p.nodes.back()];
        if (out_edges.empty()) {
          if (static_cast<int>(p.nodes.size()) >= m_min_hits()) { paths.push_back(std::move(p)); ++found; }
          continue;
        }
        for (auto it = out_edges.rbegin(); it != out_edges.rend(); ++it) {   // best edge explored first
          Path q = p;
          q.nodes.push_back(it->second);
          q.sum += it->first;
          stack.push_back(std::move(q));
        }
      }
    }
    std::sort(paths.begin(), paths.end(), [](const Path& x, const Path& y) {
      if (x.nodes.size() != y.nodes.size()) return x.nodes.size() > y.nodes.size();
      return x.sum > y.sum;
    });
    for (const auto& p : paths) {
      std::vector<int> c;
      for (int n : p.nodes) c.push_back(g.hit[n]);
      candidates.push_back(std::move(c));
    }
  }

  std::vector<char> used(hits.size(), 0);
  for (const auto& cand : candidates) {
    std::vector<int> free;
    for (int idx : cand) if (!used[idx]) free.push_back(idx);
    if (static_cast<int>(free.size()) < m_min_hits()) continue;
    ++m_candidates;
    FitResult fit;
    std::vector<int> kept;
    if (!fitCandidate(free, hits, fit, kept)) continue;
    auto* t = new Track();
    t->params = fit.params;
    t->cov = fit.cov;
    t->chi2 = fit.chi2;
    t->ndf = fit.ndf;
    t->hit_indices = kept;
    for (int idx : kept) t->planes.push_back(hits[idx]->plane);
    t->nhits = static_cast<int>(kept.size());
    matchTrackToTruth(*t, hits);
    setAcceptanceFlags(*t, m_geo(), m_fitter->propagator(), m_muon_pmin());
    if (m_ec_match() > 0 && !matchesCluster(*t)) { delete t; continue; }
    for (int idx : kept) used[idx] = 1;
    out.push_back(t);
  }
  m_tracks += static_cast<long>(out.size());
  m_tracks_out() = out;
}

void GNNTrackFinder_factory::Finish() {
  const double ev = std::max(1L, m_events);
  LOG_INFO(GetLogger()) << "GNN tracking: " << m_events << " events, " << m_nodes / ev << " hits/event, " << m_edges / ev
                        << " graph edges/event, " << m_kept_edges / ev << " kept edges/event, " << m_candidates / ev
                        << " candidates/event, " << m_tracks / ev << " tracks/event";
}

}  // namespace solid
