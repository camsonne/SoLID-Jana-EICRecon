#include "DDVCSReco_factory.h"
#include <JANA/JEvent.h>
#include <cmath>
#include "generator/DDVCSKinematicsCalc.h"

namespace solid {

DDVCSReco_factory::DDVCSReco_factory() {
  SetTypeName(NAME_OF_THIS);
  m_tracks_in.SetDatabundleName("Tracks");
  m_mc_in.SetDatabundleName("MCParticles");
  m_out.SetShortName("DDVCSEvents");
}

void DDVCSReco_factory::Process(const JEvent& event) {
  auto* ev = new DDVCSEvent();
  const auto& tracks = m_tracks_in();
  const auto& mcs = m_mc_in();
  const double Ebeam = m_geo->beamEnergy();
  // nominal beam and target used for the reconstruction
  const Vec4 k_nom(Ebeam, {0, 0, std::sqrt(Ebeam * Ebeam - kElectronMass * kElectronMass)});
  const Vec4 P_nom(kProtonMass, {0, 0, 0});

  // truth: beams and the first e', mu-, mu+ in the MC record
  int i_e = -1, i_mum = -1, i_mup = -1;
  Vec4 k_true = k_nom, P_true = P_nom;
  for (const auto* mc : mcs) {
    if (mc->status == 4) {
      if (std::abs(mc->pdg) == 11) k_true = mc->p4; else P_true = mc->p4;
      continue;
    }
    if (mc->status != 1) continue;
    if (mc->pdg == 11 && i_e < 0) i_e = mc->index;
    else if (mc->pdg == 13 && i_mum < 0) i_mum = mc->index;
    else if (mc->pdg == -13 && i_mup < 0) i_mup = mc->index;
  }
  if (i_e >= 0 && i_mum >= 0 && i_mup >= 0) {
    const auto& e = *mcs[i_e]; const auto& mm = *mcs[i_mum]; const auto& mp = *mcs[i_mup];
    const DDVCSKinematics kt = computeDDVCSKinematics(k_true, P_true, e.p4, mm.p4, mp.p4);
    ev->mm2_true = kt.mm2; ev->Q2_true = kt.Q2; ev->xB_true = kt.xB; ev->Qp2_true = kt.Qp2; ev->t_true = kt.t;
    ev->phi_true = kt.phi_trento; ev->cos_theta_l_true = kt.cos_theta_l; ev->phi_l_true = kt.phi_l;
    ev->geom_accepted = (e.hits_faec || e.hits_laec) && (mm.hits_fa_muon || mm.hits_la_muon) &&
                        (mp.hits_fa_muon || mp.hits_la_muon);
    ev->accepted_true = ev->geom_accepted && e.n_gem_hits >= 4 && mm.n_gem_hits >= 4 && mp.n_gem_hits >= 4;
  }

  // candidate selection
  auto better = [](const Track* a, const Track* b) {
    if (a->nhits != b->nhits) return a->nhits > b->nhits;
    return a->chi2 / std::max(1, a->ndf) < b->chi2 / std::max(1, b->ndf);
  };
  const Track* te = nullptr; const Track* tmup = nullptr; const Track* tmum = nullptr;
  int ie = -1, imup = -1, imum = -1;
  for (size_t i = 0; i < tracks.size(); ++i) {
    const Track* t = tracks[i];
    const bool is_muon_like = t->reaches_fa_muon || t->reaches_la_muon;
    const bool is_calo = t->reaches_faec || t->reaches_laec;
    if (m_truth_pid()) {
      if (t->mc_index < 0) continue;
      const int pdg = mcs[t->mc_index]->pdg;
      if (pdg == 11 && is_calo && (!te || better(t, te))) { te = t; ie = static_cast<int>(i); }
      else if (pdg == 13 && is_muon_like && (!tmum || better(t, tmum))) { tmum = t; imum = static_cast<int>(i); }
      else if (pdg == -13 && is_muon_like && (!tmup || better(t, tmup))) { tmup = t; imup = static_cast<int>(i); }
    } else {
      // charge + acceptance only (ambiguous for the two negative tracks: take the muon-like one as mu-)
      if (t->charge() > 0 && is_muon_like && (!tmup || better(t, tmup))) { tmup = t; imup = static_cast<int>(i); }
      else if (t->charge() < 0 && is_muon_like && (!tmum || better(t, tmum))) { tmum = t; imum = static_cast<int>(i); }
      else if (t->charge() < 0 && is_calo && (!te || better(t, te))) { te = t; ie = static_cast<int>(i); }
    }
  }
  if (te && tmup && tmum) {
    ev->found = true;
    ev->electron_track = ie; ev->muplus_track = imup; ev->muminus_track = imum;
    ev->electron = Vec4::fromPM(te->params.momentum(), kElectronMass);
    ev->muplus = Vec4::fromPM(tmup->params.momentum(), kMuonMass);
    ev->muminus = Vec4::fromPM(tmum->params.momentum(), kMuonMass);
    const DDVCSKinematics k = computeDDVCSKinematics(k_nom, P_nom, ev->electron, ev->muminus, ev->muplus);
    ev->mm2 = k.mm2; ev->Q2 = k.Q2; ev->xB = k.xB; ev->Qprime2 = k.Qp2; ev->t = k.t; ev->W = k.W;
    ev->phi = k.phi_trento; ev->cos_theta_l = k.cos_theta_l; ev->phi_l = k.phi_l;
    ev->missing = k_nom + P_nom - ev->electron - ev->muplus - ev->muminus;
    ev->all_matched = (te->mc_index == i_e) && (tmup->mc_index == i_mup) && (tmum->mc_index == i_mum);
  }
  m_out().push_back(ev);
}

}  // namespace solid
