#include "TrackFinder_factory.h"
#include <JANA/JEvent.h>
#include <algorithm>
#include <cmath>
#include <map>

namespace solid {

TrackFinder_factory::TrackFinder_factory() {
  SetTypeName(NAME_OF_THIS);
  m_hits_in.SetDatabundleName("GEMHits");
  m_ec_in.SetDatabundleName("ECClusters");
  m_ec_in.SetOptional(true);
  m_tracks_out.SetShortName("Tracks");
}

void TrackFinder_factory::Init() {
  m_fitter = std::make_unique<TrackFitter>(&m_geo(), &m_field(), m_step());
}

namespace {
FitHit toFitHit(const GEMHit& h, int index) {
  FitHit f;
  f.plane = h.plane; f.x = h.x; f.y = h.y; f.z = h.z; f.sx = h.sigma_x; f.sy = h.sigma_y; f.index = index;
  return f;
}
FitHit toFitHit(const ECCluster& c) {
  FitHit f;
  f.plane = 0; f.x = c.x; f.y = c.y; f.z = c.z; f.sx = f.sy = c.sigma_xy; f.index = -1;
  return f;
}
}  // namespace

bool TrackFinder_factory::predictAt(const TrackParams& par, const std::array<double, 16>& cov,
                                    const std::vector<FitHit>& gem_hits, int plane, Vec3& pos, double& sigma) {
  const auto& g = m_geo->gemPlane(plane);
  // measurement list sorted in z: the GEM hits already on the track (material) plus the probe
  std::vector<FitHit> list = gem_hits;
  FitHit probe; probe.plane = plane; probe.z = g.z; probe.index = -2;
  list.push_back(probe);
  std::sort(list.begin(), list.end(), [](const FitHit& a, const FitHit& b) { return a.z < b.z; });
  std::vector<TrackState> states;
  if (!m_fitter->predict(par, list, states)) return false;
  size_t ip = 0;
  for (; ip < list.size(); ++ip) if (list[ip].index == -2) break;
  pos = states[ip].pos;
  std::vector<double> V;
  m_fitter->measurementCovariance(par, list, states, V);
  const int n = 2 * static_cast<int>(list.size());
  double sig2 = 0.5 * (V[(2 * ip) * n + 2 * ip] + V[(2 * ip + 1) * n + 2 * ip + 1]);
  // parameter-uncertainty contribution
  const double dpar[4] = {std::max(1e-5, 1e-3 * std::fabs(par.qoverp)), 1e-4, 1e-4, 0.1};
  double J[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  for (int k = 0; k < 4; ++k) {
    TrackParams pp = par;
    (k == 0 ? pp.qoverp : k == 1 ? pp.theta : k == 2 ? pp.phi : pp.zv) += dpar[k];
    TrackState s; s.pos = {0, 0, pp.zv}; s.mom = pp.momentum(); s.charge = pp.charge();
    if (!m_fitter->propagator().propagateToZ(s, g.z)) return false;
    J[k] = (s.pos.x - pos.x) / dpar[k];
    J[4 + k] = (s.pos.y - pos.y) / dpar[k];
  }
  double var_par = 0;
  for (int k = 0; k < 4; ++k)
    for (int l = 0; l < 4; ++l) var_par += 0.5 * (J[k] * J[l] + J[4 + k] * J[4 + l]) * cov[k * 4 + l];
  sigma = std::sqrt(std::max(0.0, sig2 + var_par));
  return true;
}

// Progressive extension of a fitted seed.
bool TrackFinder_factory::extend(const FitResult& seedfit, const std::vector<FitHit>& seed_hits,
                                 const std::vector<std::vector<int>>& by_plane, const std::vector<const GEMHit*>& hits,
                                 const std::vector<char>& used, const ECCluster* ec, Candidate& cand) {
  const double zv0 = m_geo->target().z_center;
  cand.hits.clear(); cand.indices.clear(); cand.seed_score = 0;
  std::vector<FitHit> gem;   // GEM hits on the candidate, in z order
  double z_ref = 0;
  for (const auto& h : seed_hits) {
    if (h.plane > 0) { gem.push_back(h); cand.indices.push_back(h.index); z_ref = h.z; }
  }
  cand.fit = seedfit;
  cand.seed = seedfit.params;
  // planes to visit, nearest to the seed plane first
  std::vector<int> todo;
  for (const auto& g : m_geo->gemPlanes()) {
    bool have = false;
    for (const auto& h : gem) if (h.plane == g.id) have = true;
    if (!have) todo.push_back(g.id);
  }
  std::sort(todo.begin(), todo.end(), [&](int a, int b) {
    return std::fabs(m_geo->gemPlane(a).z - z_ref) < std::fabs(m_geo->gemPlane(b).z - z_ref);
  });
  for (int pl : todo) {
    Vec3 pos; double sigma;
    if (!predictAt(cand.fit.params, cand.fit.cov, gem, pl, pos, sigma)) continue;
    const auto& g = m_geo->gemPlane(pl);
    const double r = pos.perp();
    const double win = std::min(m_seed_window(), std::max(m_tight_window(), m_nsigma() * sigma + m_window_margin()));
    if (r < g.rin - win || r > g.rout + win) continue;
    int best = -1; double best_d2 = win * win;
    for (int idx : by_plane[pl]) {
      if (used[idx]) continue;
      const double dx = hits[idx]->x - pos.x, dy = hits[idx]->y - pos.y;
      const double d2 = dx * dx + dy * dy;
      if (d2 < best_d2) { best_d2 = d2; best = idx; }
    }
    if (best < 0) continue;
    std::vector<FitHit> trial = gem;
    trial.push_back(toFitHit(*hits[best], best));
    std::sort(trial.begin(), trial.end(), [](const FitHit& a, const FitHit& b) { return a.z < b.z; });
    std::vector<FitHit> meas = trial;
    if (ec) meas.push_back(toFitHit(*ec));
    FitResult fr = m_fitter->fit(cand.fit.params, meas, zv0, m_zv_sigma(), 4);
    if (!fr.ok) continue;
    gem = trial;
    cand.indices.push_back(best);
    cand.fit = fr;
    cand.seed_score += std::sqrt(best_d2) / std::max(sigma, 1e-3);
  }
  cand.hits = gem;
  if (ec) cand.hits.push_back(toFitHit(*ec));
  return static_cast<int>(cand.indices.size()) >= m_min_hits();
}

// Fit, re-collect hits in a tight window around the fitted trajectory, refit, reject outliers.
bool TrackFinder_factory::refine(Candidate& cand, const std::vector<std::vector<int>>& by_plane,
                                 const std::vector<const GEMHit*>& hits, const std::vector<char>& used,
                                 const ECCluster* ec) {
  const double zv0 = m_geo->target().z_center;
  if (!cand.fit.ok) return false;

  for (int pass = 0; pass < 2; ++pass) {
    // re-collect: predicted positions and errors on every plane
    std::vector<FitHit> probes;
    for (const auto& g : m_geo->gemPlanes()) { FitHit f; f.plane = g.id; f.z = g.z; probes.push_back(f); }
    std::vector<TrackState> states;
    std::vector<FitHit> reachable;
    {
      TrackState s; s.pos = {0, 0, cand.fit.params.zv}; s.mom = cand.fit.params.momentum(); s.charge = cand.fit.params.charge();
      for (auto& f : probes) {
        if (!m_fitter->propagator().propagateToZ(s, f.z)) break;
        const auto& g = m_geo->gemPlane(f.plane);
        const double r = s.pos.perp();
        if (r < g.rin || r > g.rout) continue;
        states.push_back(s); reachable.push_back(f);
      }
    }
    // predicted position uncertainty from the parameter covariance (numerical derivative) + MS
    std::vector<double> V;
    std::vector<FitHit> dummy = reachable;
    m_fitter->measurementCovariance(cand.fit.params, dummy, states, V);
    const int n = 2 * static_cast<int>(reachable.size());
    std::vector<FitHit> newhits; std::vector<int> newidx;
    for (size_t i = 0; i < reachable.size(); ++i) {
      double sig2 = 0.5 * (V[(2 * i) * n + 2 * i] + V[(2 * i + 1) * n + 2 * i + 1]);
      // parameter-uncertainty contribution: dx/dpar via finite differences
      double var_par = 0;
      {
        const double dpar[4] = {std::max(1e-5, 1e-3 * std::fabs(cand.fit.params.qoverp)), 1e-4, 1e-4, 0.1};
        std::vector<double> J(8, 0.0);
        for (int k = 0; k < 4; ++k) {
          TrackParams pp = cand.fit.params;
          (k == 0 ? pp.qoverp : k == 1 ? pp.theta : k == 2 ? pp.phi : pp.zv) += dpar[k];
          TrackState s; s.pos = {0, 0, pp.zv}; s.mom = pp.momentum(); s.charge = pp.charge();
          if (!m_fitter->propagator().propagateToZ(s, reachable[i].z)) continue;
          J[k] = (s.pos.x - states[i].pos.x) / dpar[k];
          J[4 + k] = (s.pos.y - states[i].pos.y) / dpar[k];
        }
        for (int k = 0; k < 4; ++k)
          for (int l = 0; l < 4; ++l) var_par += 0.5 * (J[k] * J[l] + J[4 + k] * J[4 + l]) * cand.fit.cov[k * 4 + l];
      }
      const double win = std::max(m_tight_window(), m_nsigma() * std::sqrt(std::max(0.0, sig2 + var_par)));
      int best = -1; double best_d2 = win * win;
      for (int idx : by_plane[reachable[i].plane]) {
        if (used[idx]) continue;
        const double dx = hits[idx]->x - states[i].pos.x, dy = hits[idx]->y - states[i].pos.y;
        const double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; best = idx; }
      }
      if (best >= 0) { newhits.push_back(toFitHit(*hits[best], best)); newidx.push_back(best); }
    }
    if (static_cast<int>(newidx.size()) < m_min_hits()) return false;
    if (ec) newhits.push_back(toFitHit(*ec));
    const bool same = (newidx == cand.indices);
    cand.hits = newhits; cand.indices = newidx;
    FitResult fr = m_fitter->fit(cand.fit.params, cand.hits, zv0, m_zv_sigma());
    if (!fr.ok) return false;
    cand.fit = fr;
    if (same) break;
  }
  // outlier rejection: drop the GEM hit with the largest pull while chi2/ndf is too large
  while (cand.fit.ndf > 0 && cand.fit.chi2 / cand.fit.ndf > m_chi2ndf_max() &&
         static_cast<int>(cand.indices.size()) > m_min_hits()) {
    int worst = -1; double worst_pull = 0;
    for (size_t i = 0; i < cand.indices.size(); ++i) {
      double pull = std::hypot(cand.fit.pulls[2 * i], cand.fit.pulls[2 * i + 1]);
      if (pull > worst_pull) { worst_pull = pull; worst = static_cast<int>(i); }
    }
    if (worst < 0) break;
    cand.hits.erase(cand.hits.begin() + worst);
    cand.indices.erase(cand.indices.begin() + worst);
    FitResult fr = m_fitter->fit(cand.fit.params, cand.hits, zv0, m_zv_sigma());
    if (!fr.ok) return false;
    cand.fit = fr;
  }
  return cand.fit.ndf > 0 && cand.fit.chi2 / cand.fit.ndf <= m_chi2ndf_max();
}

Track* TrackFinder_factory::makeTrack(const Candidate& cand, const std::vector<const GEMHit*>& hits) {
  auto* t = new Track();
  t->params = cand.fit.params;
  t->cov = cand.fit.cov;
  t->chi2 = cand.fit.chi2;
  t->ndf = cand.fit.ndf;
  t->nhits = static_cast<int>(cand.indices.size());
  t->hit_indices = cand.indices;
  std::map<int, int> votes;
  for (int idx : cand.indices) {
    t->planes.push_back(hits[idx]->plane);
    if (hits[idx]->mc_index >= 0) votes[hits[idx]->mc_index]++;
  }
  int best = -1, nbest = 0;
  for (auto& kv : votes) if (kv.second > nbest) { nbest = kv.second; best = kv.first; }
  if (nbest * 2 > t->nhits) { t->mc_index = best; t->n_bkg_hits = t->nhits - nbest; t->purity = double(nbest) / t->nhits; }
  else { t->mc_index = -1; t->n_bkg_hits = t->nhits - nbest; t->purity = double(nbest) / t->nhits; }
  setAcceptanceFlags(*t);
  return t;
}

void TrackFinder_factory::setAcceptanceFlags(Track& trk) {
  const auto& geo = m_geo();
  const auto& prop = m_fitter->propagator();
  TrackState s; s.pos = {0, 0, trk.params.zv}; s.mom = trk.params.momentum(); s.charge = trk.params.charge();
  // LAEC
  TrackState t = s;
  if (prop.propagateToZ(t, geo.laec().z)) {
    const double r = t.pos.perp();
    trk.reaches_laec = (r >= geo.laec().rin && r <= geo.laec().rout);
    s = t;
    if (prop.propagateToZ(t, geo.faec().z)) {
      const double r2 = t.pos.perp();
      trk.reaches_faec = (r2 >= geo.faec().rin && r2 <= geo.faec().rout);
      if (trk.reaches_faec && trk.p() > m_muon_pmin()) {
        bool ok = true;
        for (const auto& sc : geo.muon().fa_scint) {
          if (!prop.propagateToZ(t, sc.z)) { ok = false; break; }
          const double r3 = t.pos.perp();
          if (r3 < sc.rin || r3 > sc.rout) { ok = false; break; }
        }
        trk.reaches_fa_muon = ok;
        trk.p_at_muon = trk.p() - 0.0114 * geo.muon().fa_iron_total - 0.4;  // ~11.4 MeV/cm Fe + FAEC
      }
    }
  }
  if (trk.reaches_laec && !trk.reaches_fa_muon) {
    TrackState u = s;
    const auto& mu = geo.muon();
    if (prop.propagateToR(u, mu.la_radius, mu.la_zmax) && u.pos.z >= mu.la_zmin && u.pos.z <= mu.la_zmax &&
        trk.p() > 0.5 * m_muon_pmin())
      trk.reaches_la_muon = true;
  }
}

void TrackFinder_factory::Process(const JEvent& event) {
  const auto& hits = m_hits_in();
  std::vector<std::vector<int>> by_plane(7);
  for (size_t i = 0; i < hits.size(); ++i) {
    if (hits[i]->plane >= 1 && hits[i]->plane <= 6) by_plane[hits[i]->plane].push_back(static_cast<int>(i));
  }
  std::vector<char> used(hits.size(), 0);
  std::vector<Track*> out;
  const double ztgt = m_geo->target().z_center;

  auto processSeeds = [&](std::vector<Candidate>& cands, const ECCluster* ec) {
    // rank: more hits first, then smaller seed residuals
    std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b) {
      if (a.indices.size() != b.indices.size()) return a.indices.size() > b.indices.size();
      return a.seed_score < b.seed_score;
    });
    Candidate best; bool have = false;
    int nfit = 0;
    for (auto& c : cands) {
      if (nfit >= m_max_fits()) break;
      ++nfit;
      if (!refine(c, by_plane, hits, used, ec)) continue;
      auto better = [](const Candidate& a, const Candidate& b) {
        if (a.indices.size() != b.indices.size()) return a.indices.size() > b.indices.size();
        return a.fit.chi2 / std::max(1, a.fit.ndf) < b.fit.chi2 / std::max(1, b.fit.ndf);
      };
      if (!have || better(c, best)) { best = c; have = true; }
    }
    if (have) {
      for (int idx : best.indices) used[idx] = 1;
      out.push_back(makeTrack(best, hits));
    }
    return have;
  };

  const bool ec_seeding = m_use_ec() && !m_ec_in().empty();
  if (ec_seeding) {
    // order clusters by energy (electron first)
    std::vector<const ECCluster*> clusters(m_ec_in().begin(), m_ec_in().end());
    std::sort(clusters.begin(), clusters.end(), [](const ECCluster* a, const ECCluster* b) { return a->energy > b->energy; });
    for (const ECCluster* ec : clusters) {
      const FitHit ecf = toFitHit(*ec);
      const double phi_ec = std::atan2(ec->y, ec->x);
      const double r_ec = std::hypot(ec->x, ec->y);
      const std::vector<int> seed_planes = (ec->calo == 1) ? std::vector<int>{3, 2, 4, 1} : std::vector<int>{4, 3, 5, 2, 6};
      bool found = false;
      for (int pl : seed_planes) {
        if (static_cast<int>(by_plane[pl].size()) > m_max_hits_plane()) continue;  // saturated plane
        std::vector<Candidate> cands;
        // seed-plane hits passing the prefilters, ranked by the distance to the straight line
        // from the target centre to the calorimeter cluster (a cheap road); only the best
        // tracking:max_seeds_per_cluster are tried, which bounds the cost on saturated planes
        const double zpl = m_geo->gemPlane(pl).z;
        const double fline = (zpl - ztgt) / (ec->z - ztgt);
        const double xl = ec->x * fline, yl = ec->y * fline;
        std::vector<std::pair<double, int>> ranked;
        for (int idx : by_plane[pl]) {
          if (used[idx]) continue;
          const GEMHit& h = *hits[idx];
          // prefilters: azimuthal compatibility and (R,z) line pointing back to the target
          double dphi = std::remainder(h.phi() - phi_ec, 2 * M_PI);
          if (std::fabs(dphi) > 1.2) continue;
          const double r_h = h.r();
          if (std::fabs(r_ec - r_h) < 1e-3) continue;
          const double zv = h.z - r_h * (ec->z - h.z) / (r_ec - r_h);
          if (std::fabs(zv - ztgt) > 80.0) continue;
          ranked.emplace_back(std::hypot(h.x - xl, h.y - yl), idx);
        }
        std::sort(ranked.begin(), ranked.end());
        if (static_cast<int>(ranked.size()) > m_max_seeds()) ranked.resize(m_max_seeds());
        for (const auto& rk : ranked) {
          const int idx = rk.second;
          const GEMHit& h = *hits[idx];
          for (int sign : {+1, -1}) {
            TrackParams seed;
            if (!m_fitter->seedFromTwoHits(toFitHit(h, idx), ecf, m_b_eff(), seed)) continue;
            if (sign < 0) {
              // the effective-field seed can misassign the charge for low curvature: try both
              seed.qoverp = -seed.qoverp;
              if (std::fabs(seed.qoverp) < 0.2) continue;  // only for p > 5 GeV where charge is ambiguous
            }
            if (seed.p() < m_pmin() || seed.p() > 30.0) continue;
            // refine the analytic seed with a 2-point fit (GEM hit + calorimeter point + vertex prior)
            // so that the extension does not depend on the effective-field approximation
            std::vector<FitHit> two = {toFitHit(h, idx), ecf};
            FitResult sf = m_fitter->fit(seed, two, ztgt, m_zv_sigma(), 6);
            if (!sf.ok) continue;
            if (sf.params.p() < m_pmin() || sf.params.p() > 30.0) continue;
            Candidate c;
            if (extend(sf, two, by_plane, hits, used, ec, c)) cands.push_back(std::move(c));
          }
        }
        if (processSeeds(cands, ec)) { found = true; break; }
      }
      (void)found;
    }
  } else {
    // GEM-pair seeding fallback: pairs of planes, vertex-line prefilter
    const std::vector<std::pair<int, int>> pairs = {{2, 3}, {3, 4}, {1, 2}, {2, 4}, {4, 5}, {1, 3}};
    for (const auto& pr : pairs) {
      if (static_cast<int>(by_plane[pr.first].size()) > m_max_hits_plane() ||
          static_cast<int>(by_plane[pr.second].size()) > m_max_hits_plane()) continue;
      std::vector<Candidate> cands;
      for (int ia : by_plane[pr.first]) {
        if (used[ia]) continue;
        for (int ib : by_plane[pr.second]) {
          if (used[ib]) continue;
          const GEMHit& a = *hits[ia];
          const GEMHit& b = *hits[ib];
          if (std::fabs(std::remainder(a.phi() - b.phi(), 2 * M_PI)) > 0.5) continue;
          const double ra = a.r(), rb = b.r();
          if (rb <= ra) continue;
          const double zv = a.z - ra * (b.z - a.z) / (rb - ra);
          if (std::fabs(zv - ztgt) > 40.0) continue;
          TrackParams seed;
          if (!m_fitter->seedFromTwoHits(toFitHit(a, ia), toFitHit(b, ib), m_b_eff(), seed)) continue;
          if (seed.p() < m_pmin() || seed.p() > 30.0) continue;
          std::vector<FitHit> two = {toFitHit(a, ia), toFitHit(b, ib)};
          FitResult sf = m_fitter->fit(seed, two, ztgt, m_zv_sigma(), 6);
          if (!sf.ok) continue;
          Candidate c;
          if (extend(sf, two, by_plane, hits, used, nullptr, c)) cands.push_back(std::move(c));
        }
      }
      // fit the best candidates repeatedly until none is left
      while (!cands.empty()) {
        std::vector<Candidate> alive;
        for (auto& c : cands) {
          bool ok = true;
          for (int idx : c.indices) if (used[idx]) { ok = false; break; }
          if (ok) alive.push_back(c);
        }
        if (alive.empty()) break;
        if (!processSeeds(alive, nullptr)) break;
        cands.swap(alive);
      }
    }
  }
  m_tracks_out() = out;
}

}  // namespace solid
