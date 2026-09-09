#include "DDVCSAnalysis_processor.h"
#include <JANA/JEvent.h>
#include <JANA/JLogger.h>
#include <algorithm>
#include <cmath>
#include <iomanip>

namespace solid {

double DDVCSAnalysis_processor::Stats::mean() const {
  if (v.empty()) return 0;
  double s = 0; for (double x : v) s += x; return s / v.size();
}
double DDVCSAnalysis_processor::Stats::rms() const {
  if (v.size() < 2) return 0;
  double m = mean(), s = 0; for (double x : v) s += (x - m) * (x - m); return std::sqrt(s / (v.size() - 1));
}
double DDVCSAnalysis_processor::Stats::quantile(double q) const {
  if (v.empty()) return 0;
  std::vector<double> s = v; std::sort(s.begin(), s.end());
  double pos = q * (s.size() - 1); size_t i = static_cast<size_t>(pos); double f = pos - i;
  if (i + 1 >= s.size()) return s.back();
  return s[i] * (1 - f) + s[i + 1] * f;
}
double DDVCSAnalysis_processor::Stats::gaussCoreSigma(double& mean_out) const {
  if (v.size() < 10) { mean_out = mean(); return rms(); }
  double m = quantile(0.5), s = sigma68();
  for (int it = 0; it < 10; ++it) {
    double sum = 0, sum2 = 0; long n = 0;
    for (double x : v) {
      if (std::fabs(x - m) < 2.5 * s) { sum += x; sum2 += x * x; ++n; }
    }
    if (n < 5) break;
    double nm = sum / n, ns = std::sqrt(std::max(0.0, sum2 / n - nm * nm));
    // correct for the 2.5 sigma truncation of a Gaussian (variance fraction 0.9838 -> sigma factor)
    ns /= 0.9738;
    if (std::fabs(ns - s) < 1e-4 * s) { m = nm; s = ns; break; }
    m = nm; s = ns;
  }
  mean_out = m;
  return s;
}

DDVCSAnalysis_processor::DDVCSAnalysis_processor() {
  SetTypeName(NAME_OF_THIS);
  SetCallbackStyle(CallbackStyle::ExpertMode);
  m_ev_in.SetDatabundleName("DDVCSEvents");
  m_gen_in.SetDatabundleName("GenKinematics");
  m_gen_in.SetOptional(true);
  m_tracks_in.SetDatabundleName("Tracks");
  m_mc_in.SetDatabundleName("MCParticles");
}

void DDVCSAnalysis_processor::Init() {
  m_out.open(m_csv());
  m_out << "event,Q2_true,xB_true,Qp2_true,t_true,geom_accepted,accepted,found,all_matched,mm2,mm2_true,Q2,Qp2,t,"
           "pe_true,pe_reco,the_true,the_reco,pmup_true,pmup_reco,pmum_true,pmum_reco,"
           "ntracks,nhits_e,nhits_mup,nhits_mum,nbkg_e,nbkg_mup,nbkg_mum,chi2ndf_e,chi2ndf_mup,chi2ndf_mum\n";
  m_out << std::setprecision(7);
}

void DDVCSAnalysis_processor::ProcessSequential(const JEvent& event) {
  const auto& evs = m_ev_in();
  if (evs.empty()) return;
  const DDVCSEvent& ev = *evs[0];
  const auto& tracks = m_tracks_in();
  const auto& mcs = m_mc_in();
  if (!m_gen_in().empty()) m_lumi = m_gen_in()[0]->luminosity;

  ++m_n_events;
  if (ev.geom_accepted) ++m_n_geom;
  if (ev.geom_accepted && ev.found) ++m_n_found_geom;
  if (ev.accepted_true) ++m_n_accepted;
  if (ev.found) ++m_n_found;
  if (ev.found && ev.accepted_true) ++m_n_found_accepted;
  if (ev.found && ev.all_matched) ++m_n_all_matched;

  for (const auto* t : tracks) {
    ++m_n_tracks;
    if (t->n_bkg_hits > 0) ++m_n_tracks_with_bkg;
    if (t->mc_index < 0) ++m_n_fake_tracks;
    m_chi2ndf.add(t->chi2 / std::max(1, t->ndf));
    if (t->mc_index >= 0) {
      const auto& mc = *mcs[t->mc_index];
      const double pt = mc.p4.p.mag(), pr = t->p();
      const double dpp = (pr - pt) / pt;
      m_dzv.add(t->params.zv - mc.vertex.z);
      if (std::abs(mc.pdg) == 11) { m_dpp_e.add(dpp); m_dth_e.add(t->params.theta - mc.p4.p.theta()); }
      if (std::abs(mc.pdg) == 13) {
        m_dpp_mu.add(dpp); m_dth_mu.add(t->params.theta - mc.p4.p.theta());
        m_dph_mu.add(std::remainder(t->params.phi - mc.p4.p.phi(), 2 * M_PI));
      }
    }
  }

  const Track* te = ev.electron_track >= 0 ? tracks[ev.electron_track] : nullptr;
  const Track* tp = ev.muplus_track >= 0 ? tracks[ev.muplus_track] : nullptr;
  const Track* tm = ev.muminus_track >= 0 ? tracks[ev.muminus_track] : nullptr;
  const MCParticle* me = nullptr; const MCParticle* mp = nullptr; const MCParticle* mm = nullptr;
  for (const auto* mc : mcs) {
    if (mc->status != 1) continue;
    if (mc->pdg == 11) me = mc; else if (mc->pdg == -13) mp = mc; else if (mc->pdg == 13) mm = mc;
  }
  if (ev.found) {
    m_mm2.add(ev.mm2);
    if (ev.mm2 > 0) m_mm.add(std::sqrt(ev.mm2));
    if (ev.all_matched && te->n_bkg_hits == 0 && tp->n_bkg_hits == 0 && tm->n_bkg_hits == 0) m_mm2_clean.add(ev.mm2);
    m_dQ2.add(ev.Q2 - ev.Q2_true);
    m_dQp2.add(ev.Qprime2 - ev.Qp2_true);
    m_dt.add(ev.t - ev.t_true);
  }

  auto tv = [](const Track* t, auto f, double def = 0) { return t ? f(t) : def; };
  m_out << event.GetEventNumber() << "," << ev.Q2_true << "," << ev.xB_true << "," << ev.Qp2_true << "," << ev.t_true << ","
        << ev.geom_accepted << "," << ev.accepted_true << "," << ev.found << "," << ev.all_matched << "," << ev.mm2 << "," << ev.mm2_true << ","
        << ev.Q2 << "," << ev.Qprime2 << "," << ev.t << ","
        << (me ? me->p4.p.mag() : 0) << "," << tv(te, [](const Track* t) { return t->p(); }) << ","
        << (me ? me->p4.p.theta() : 0) << "," << tv(te, [](const Track* t) { return t->params.theta; }) << ","
        << (mp ? mp->p4.p.mag() : 0) << "," << tv(tp, [](const Track* t) { return t->p(); }) << ","
        << (mm ? mm->p4.p.mag() : 0) << "," << tv(tm, [](const Track* t) { return t->p(); }) << ","
        << tracks.size() << ","
        << tv(te, [](const Track* t) { return double(t->nhits); }) << "," << tv(tp, [](const Track* t) { return double(t->nhits); }) << ","
        << tv(tm, [](const Track* t) { return double(t->nhits); }) << ","
        << tv(te, [](const Track* t) { return double(t->n_bkg_hits); }) << "," << tv(tp, [](const Track* t) { return double(t->n_bkg_hits); }) << ","
        << tv(tm, [](const Track* t) { return double(t->n_bkg_hits); }) << ","
        << tv(te, [](const Track* t) { return t->chi2 / std::max(1, t->ndf); }) << ","
        << tv(tp, [](const Track* t) { return t->chi2 / std::max(1, t->ndf); }) << ","
        << tv(tm, [](const Track* t) { return t->chi2 / std::max(1, t->ndf); }) << "\n";
}

void DDVCSAnalysis_processor::Finish() {
  m_out.close();
  std::ofstream js(m_json());
  js << std::setprecision(6);
  auto stat = [&](const char* name, Stats& s, bool last = false) {
    double gm = 0, gs = s.gaussCoreSigma(gm);
    js << "  \"" << name << "\": {\"n\": " << s.n() << ", \"mean\": " << s.mean() << ", \"rms\": " << s.rms()
       << ", \"median\": " << s.quantile(0.5) << ", \"sigma68\": " << s.sigma68() << ", \"sigma95\": " << s.sigma95()
       << ", \"gauss_mean\": " << gm << ", \"gauss_sigma\": " << gs << "}" << (last ? "\n" : ",\n");
  };
  js << "{\n";
  js << "  \"label\": \"" << m_label() << "\",\n";
  js << "  \"luminosity\": " << m_lumi << ",\n";
  js << "  \"n_events\": " << m_n_events << ",\n";
  js << "  \"n_geom_accepted\": " << m_n_geom << ",\n";
  js << "  \"n_found_geom_accepted\": " << m_n_found_geom << ",\n";
  js << "  \"geometric_acceptance\": " << (m_n_events ? double(m_n_geom) / m_n_events : 0) << ",\n";
  js << "  \"efficiency_vs_geometric_acceptance\": " << (m_n_geom ? double(m_n_found_geom) / m_n_geom : 0) << ",\n";
  js << "  \"n_accepted_true\": " << m_n_accepted << ",\n";
  js << "  \"n_found\": " << m_n_found << ",\n";
  js << "  \"n_found_and_accepted\": " << m_n_found_accepted << ",\n";
  js << "  \"n_found_all_matched\": " << m_n_all_matched << ",\n";
  js << "  \"acceptance\": " << (m_n_events ? double(m_n_accepted) / m_n_events : 0) << ",\n";
  js << "  \"reco_efficiency_in_acceptance\": " << (m_n_accepted ? double(m_n_found_accepted) / m_n_accepted : 0) << ",\n";
  js << "  \"n_tracks\": " << m_n_tracks << ",\n";
  js << "  \"frac_tracks_with_bkg_hit\": " << (m_n_tracks ? double(m_n_tracks_with_bkg) / m_n_tracks : 0) << ",\n";
  js << "  \"frac_fake_tracks\": " << (m_n_tracks ? double(m_n_fake_tracks) / m_n_tracks : 0) << ",\n";
  stat("mm2", m_mm2);
  stat("mm", m_mm);
  stat("mm2_clean_tracks", m_mm2_clean);
  stat("dp_over_p_electron", m_dpp_e);
  stat("dp_over_p_muon", m_dpp_mu);
  stat("dtheta_electron", m_dth_e);
  stat("dtheta_muon", m_dth_mu);
  stat("dphi_muon", m_dph_mu);
  stat("dzv", m_dzv);
  stat("dQ2", m_dQ2);
  stat("dQprime2", m_dQp2);
  stat("dt", m_dt);
  stat("chi2ndf", m_chi2ndf, true);
  js << "}\n";
  js.close();
  double gm = 0, gs = m_mm2.gaussCoreSigma(gm);
  LOG_INFO(GetLogger()) << "DDVCS summary (L=" << m_lumi << "): events=" << m_n_events << " accepted=" << m_n_accepted
                        << " found=" << m_n_found << " eff(in acc)=" << (m_n_accepted ? double(m_n_found_accepted) / m_n_accepted : 0)
                        << "  MM2: sigma68=" << m_mm2.sigma68() << " GeV^2, Gaussian core sigma=" << gs << " GeV^2 (mean " << gm
                        << ");  MM: sigma68=" << m_mm.sigma68() << " GeV";
}

}  // namespace solid
