#include "TrackUtils.h"
#include <map>

namespace solid {

void matchTrackToTruth(Track& trk, const std::vector<const GEMHit*>& hits) {
  std::map<int, int> votes;
  for (int idx : trk.hit_indices)
    if (hits[idx]->mc_index >= 0) votes[hits[idx]->mc_index]++;
  int best = -1, nbest = 0;
  for (auto& kv : votes) if (kv.second > nbest) { nbest = kv.second; best = kv.first; }
  trk.mc_index = (nbest * 2 > trk.nhits) ? best : -1;
  trk.n_bkg_hits = trk.nhits - nbest;
  trk.purity = trk.nhits > 0 ? double(nbest) / trk.nhits : 0;
}

void setAcceptanceFlags(Track& trk, const SoLIDGeometry_service& geo, const TrackPropagator& prop, double muon_pmin) {
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
      if (trk.reaches_faec && trk.p() > muon_pmin) {
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
        trk.p() > 0.5 * muon_pmin)
      trk.reaches_la_muon = true;
  }
}

}  // namespace solid
