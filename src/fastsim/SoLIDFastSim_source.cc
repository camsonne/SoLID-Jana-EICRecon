#include "SoLIDFastSim_source.h"
#include <JANA/JEvent.h>
#include <JANA/JLogger.h>
#include <cmath>

namespace solid {

SoLIDFastSim_source::SoLIDFastSim_source() : JEventSource() {
  SetTypeName(NAME_OF_THIS);
  SetCallbackStyle(CallbackStyle::ExpertMode);
}

void SoLIDFastSim_source::Init() {
  DDVCSConfig cfg;
  cfg.beam_energy = m_geo->beamEnergy();
  cfg.target_z = m_geo->target().z_center;
  cfg.target_half_length = m_geo->target().half_length;
  cfg.Q2_min = m_Q2_min(); cfg.Q2_max = m_Q2_max();
  cfg.xB_min = m_xB_min(); cfg.xB_max = m_xB_max();
  cfg.Qp2_min = m_Qp2_min(); cfg.Qp2_max = m_Qp2_max();
  cfg.t_slope = m_t_slope(); cfg.tabs_max = m_tabs_max();
  m_gen = std::make_unique<DDVCSGenerator>(cfg, static_cast<unsigned>(m_seed()));
  m_rng.seed(static_cast<unsigned>(m_seed()) + 7919u);
  if (!m_hepmc_in().empty()) {
    m_reader = std::make_unique<HepMC3Reader>(m_hepmc_in());
    if (!m_reader->good()) throw JException("Cannot open HepMC3 input %s", m_hepmc_in().c_str());
    if (m_hepmc_vertex() != "target" && m_hepmc_vertex() != "file")
      throw JException("fastsim:hepmc_vertex must be 'target' or 'file'");
    LOG_INFO(GetLogger()) << "Reading generated events from " << m_hepmc_in() << " (vertex: " << m_hepmc_vertex() << ")";
  }

  BackgroundConfig bcfg;
  bcfg.luminosity = m_luminosity();
  bcfg.readout = m_readout();
  bcfg.time_window_ns = m_time_window();
  bcfg.strip_pitch_mm = m_strip_pitch();
  bcfg.strip_length_cm = m_strip_length();
  bcfg.pad_area_mm2 = m_pad_area();
  bcfg.file = m_bkg_file();
  m_bkg = std::make_unique<BackgroundModel>(bcfg, &m_geo());
  m_prop = std::make_unique<TrackPropagator>(&m_field(), 2.0);
  if (!m_hepmc_out().empty()) m_hepmc = std::make_unique<HepMC3Writer>(m_hepmc_out());

  LOG_INFO(GetLogger()) << m_bkg->summary();
}

// Transport one charged particle through the detector, producing true GEM hits and recording
// the plane crossings (for background generation) and the acceptance flags.
void SoLIDFastSim_source::transport(MCParticle& mc, std::vector<GEMHit>& hits, std::vector<ECCluster>& clusters,
                                    std::vector<std::pair<int, Vec3>>& crossings) {
  auto makeCluster = [&](int calo, const Vec3& pos) {
    ECCluster c;
    c.calo = calo;
    std::normal_distribution<double> Gc(0.0, m_ec_sigma());
    c.x = pos.x + Gc(m_rng);
    c.y = pos.y + Gc(m_rng);
    c.z = pos.z;
    c.sigma_xy = m_ec_sigma();
    const bool em = std::abs(mc.pdg) == 11;
    const double E = em ? mc.p4.E : 0.3;  // muons/hadrons: MIP-like deposit
    std::normal_distribution<double> Ge(0.0, m_ec_eres() * std::sqrt(E));
    c.energy = std::max(0.0, E + Ge(m_rng));
    c.mc_index = mc.index;
    clusters.push_back(c);
  };
  const auto& geo = m_geo();
  const auto& tg = geo.target();
  const bool ms = m_ms();
  const double p = mc.p4.p.mag();
  const double beta = p / mc.p4.E;
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::normal_distribution<double> G(0.0, 1.0);

  TrackState s;
  s.pos = mc.vertex;
  s.mom = mc.p4.p;
  s.charge = static_cast<int>(mc.charge);
  if (s.dir().z <= 0.05) return;  // backward / very large angle: outside the spectrometer

  // --- target: LH2 path + Al wall or window, one lumped scatter at the exit point
  {
    const double sinth = std::max(1e-6, s.dir().perp()), costh = std::max(1e-6, s.dir().z);
    const double L_side = tg.radius / sinth, L_end = (tg.z_center + tg.half_length - s.pos.z) / costh;
    double L, xx0;
    if (L_side < L_end) { L = L_side; xx0 = L / tg.x0_lh2 + tg.wall_thickness / (tg.x0_al * sinth); }
    else { L = L_end; xx0 = L / tg.x0_lh2 + tg.window_thickness / (tg.x0_al * costh); }
    const double z_exit = s.pos.z + L * costh;
    if (!m_prop->propagateToZ(s, z_exit)) return;
    if (ms) TrackPropagator::scatter(s, highlandTheta0(p, beta, xx0), m_rng);
  }

  const double air_x0 = ms ? geo.airX0() : 0.0;
  const bool is_muon = std::abs(mc.pdg) == 13;
  const bool is_electron = std::abs(mc.pdg) == 11;
  bool stopped = false;

  // --- GEM planes
  for (const auto& g : geo.gemPlanes()) {
    if (g.z < s.pos.z) continue;
    // large-angle calorimeter sits between plane 4 and plane 5: absorbs electrons/hadrons
    const auto& laec = geo.laec();
    if (!stopped && laec.z > s.pos.z && laec.z < g.z) {
      TrackState t = s;
      if (!m_prop->propagateToZ(t, laec.z, ms ? &m_rng : nullptr, air_x0)) return;
      const double r = t.pos.perp();
      if (r >= laec.rin && r <= laec.rout) {
        mc.hits_laec = true;
        makeCluster(1, t.pos);
        if (!is_muon) stopped = true;
      }
      s = t;
    }
    if (stopped) break;
    if (!m_prop->propagateToZ(s, g.z, ms ? &m_rng : nullptr, air_x0)) return;
    const double r = s.pos.perp();
    if (r < g.rin || r > g.rout) continue;
    crossings.emplace_back(g.id, s.pos);
    // hit efficiency: intrinsic x occupancy survival
    const double eff = g.efficiency * m_bkg->survivalProbability(g.id, r);
    if (U(m_rng) < eff) {
      GEMHit h;
      h.plane = g.id;
      h.x = s.pos.x + g.resolution * G(m_rng);
      h.y = s.pos.y + g.resolution * G(m_rng);
      h.z = g.z;
      h.sigma_x = h.sigma_y = g.resolution;
      h.mc_index = mc.index;
      hits.push_back(h);
      mc.n_gem_hits++;
    }
    if (ms) TrackPropagator::scatter(s, highlandTheta0(p, beta, g.x0_fraction / std::max(0.05, s.dir().z)), m_rng);
  }
  if (stopped) return;

  // --- forward calorimeter
  const auto& faec = geo.faec();
  if (faec.z > s.pos.z) {
    TrackState t = s;
    if (!m_prop->propagateToZ(t, faec.z, ms ? &m_rng : nullptr, air_x0)) return;
    const double r = t.pos.perp();
    if (r >= faec.rin && r <= faec.rout) {
      mc.hits_faec = true;
      makeCluster(0, t.pos);
    }
    s = t;
  }
  if (is_electron || !is_muon) return;

  // --- muons: forward-angle muon detector (behind the FAEC + iron), large-angle cylinder
  if (mc.hits_faec && p > m_muon_pmin()) {
    bool ok = true;
    TrackState t = s;
    for (const auto& sc : geo.muon().fa_scint) {
      if (!m_prop->propagateToZ(t, sc.z)) { ok = false; break; }
      const double r = t.pos.perp();
      if (r < sc.rin || r > sc.rout) { ok = false; break; }
    }
    mc.hits_fa_muon = ok;
  }
  if (!mc.hits_fa_muon && mc.hits_laec) {
    TrackState t = s;
    const auto& mu = geo.muon();
    if (m_prop->propagateToR(t, mu.la_radius, mu.la_zmax) && t.pos.z >= mu.la_zmin && t.pos.z <= mu.la_zmax &&
        p > 0.5 * m_muon_pmin()) {
      mc.hits_la_muon = true;
    }
  }
}

// Next generated event, from the HepMC3 file or the internal generator. Returns false at the
// end of the input file.
bool SoLIDFastSim_source::nextEvent(std::vector<MCParticle>& parts) {
  if (!m_reader) {
    DDVCSKinematics kin;
    if (!m_gen->generate(parts, kin)) throw JException("DDVCS generator failed to produce an event");
    return true;
  }
  try {
    if (!m_reader->next(parts)) {
      const auto& run = m_reader->runInfo();
      LOG_INFO(GetLogger()) << "End of HepMC3 input after " << m_emitted << " events"
                            << (run.tool.empty() ? "" : " (" + run.tool + ")")
                            << (run.has_xsec ? ", integrated cross section " + std::to_string(run.xsec) + " +- " +
                                                   std::to_string(run.xsec_err) + " (generator units)"
                                             : "");
      return false;
    }
  } catch (const std::exception& e) {
    throw JException("%s", e.what());
  }
  if (std::abs(parts[5].pdg) != 13 || std::abs(parts[6].pdg) != 13)
    throw JException("HepMC3 event %d: lepton pair is %d %d, the SoLID DDVCS chain expects mu- mu+ "
                     "(EpIC: set decay_type to mu-)", m_reader->eventNumber(), parts[5].pdg, parts[6].pdg);
  if (m_emitted == 0) {
    const double Eb = parts[0].p4.E, Eg = m_geo->beamEnergy();
    if (std::fabs(Eb - Eg) > 1e-3 * Eg)
      LOG_WARN(GetLogger()) << "HepMC3 beam energy " << Eb << " GeV differs from solid geometry beam energy " << Eg
                            << " GeV used by the reconstruction";
    if (parts[1].p4.p.mag() > 1e-3)
      LOG_WARN(GetLogger()) << "HepMC3 target is not at rest (|p| = " << parts[1].p4.p.mag() * 1e3
                            << " MeV); the reconstruction assumes a target at rest (EpIC: set hadron_energy to Mp)";
    if (m_reader->flipped()) LOG_INFO(GetLogger()) << "HepMC3 beam along -z: events rotated by pi around y (SoLID frame)";
  }
  if (m_hepmc_vertex() == "target" || !m_reader->vertexInFile()) {
    const auto& tg = m_geo->target();
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    const Vec3 vtx(0, 0, tg.z_center + tg.half_length * U(m_rng));
    for (auto& mc : parts) mc.vertex = vtx;
  }
  return true;
}

JEventSource::Result SoLIDFastSim_source::Emit(JEvent& event) {
  if (m_nevents() > 0 && m_emitted >= m_nevents()) return Result::FailureFinished;

  std::vector<MCParticle> parts;
  if (!nextEvent(parts)) return Result::FailureFinished;
  const double weight = m_reader ? m_reader->weight() : 1.0;
  const DDVCSKinematics kin = computeDDVCSKinematics(parts[0].p4, parts[1].p4, parts[2].p4, parts[5].p4, parts[6].p4);
  if (m_hepmc) m_hepmc->write(m_emitted + 1, parts, weight);

  std::vector<GEMHit> hits;
  std::vector<ECCluster> clusters;
  std::vector<std::pair<int, Vec3>> crossings;
  for (auto& mc : parts) {
    if (mc.status != 1 || mc.charge == 0) continue;
    transport(mc, hits, clusters, crossings);
  }
  // background in a road around every plane crossing of a charged final-state particle
  const double w = m_bkg_road();
  for (const auto& c : crossings) m_bkg->generateInRoad(c.first, c.second.x, c.second.y, w, m_rng, hits);

  std::vector<MCParticle*> mc_out;
  for (const auto& mc : parts) mc_out.push_back(new MCParticle(mc));
  std::vector<GEMHit*> hit_out;
  for (const auto& h : hits) hit_out.push_back(new GEMHit(h));
  std::vector<ECCluster*> ec_out;
  for (const auto& c : clusters) ec_out.push_back(new ECCluster(c));
  auto* gk = new GenKinematics();
  gk->Q2 = kin.Q2; gk->xB = kin.xB; gk->nu = kin.nu; gk->W = kin.W; gk->Qp2 = kin.Qp2; gk->t = kin.t;
  gk->tmin = kin.tmin; gk->phi_trento = kin.phi_trento; gk->y = kin.y; gk->luminosity = m_luminosity();
  gk->cos_theta_l = kin.cos_theta_l; gk->phi_l = kin.phi_l;
  gk->weight = weight; gk->source = m_reader ? 1 : 0;

  event.SetEventNumber(++m_emitted);
  event.SetRunNumber(1);
  event.Insert(mc_out, "MCParticles");
  event.Insert(hit_out, "GEMHits");
  event.Insert(ec_out, "ECClusters");
  event.Insert(gk, "GenKinematics");
  return Result::Success;
}

}  // namespace solid
