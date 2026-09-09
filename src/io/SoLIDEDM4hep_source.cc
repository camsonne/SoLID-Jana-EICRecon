#include "SoLIDEDM4hep_source.h"
#include <JANA/JEvent.h>
#include <JANA/JLogger.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <cmath>
#include <map>

namespace solid {

SoLIDEDM4hep_source::SoLIDEDM4hep_source(std::string resource_name, JApplication* app) : JEventSource(resource_name, app) {
  SetTypeName(NAME_OF_THIS);
  SetCallbackStyle(CallbackStyle::ExpertMode);
}

SoLIDEDM4hep_source::~SoLIDEDM4hep_source() = default;

void SoLIDEDM4hep_source::Init() {
  BackgroundConfig bcfg;
  bcfg.luminosity = m_luminosity();
  bcfg.readout = m_readout();
  bcfg.time_window_ns = m_time_window();
  bcfg.file = m_bkg_file();
  m_bkg = std::make_unique<BackgroundModel>(bcfg, &m_geo());
  m_rng.seed(static_cast<unsigned>(m_seed()));
  LOG_INFO(GetLogger()) << m_bkg->summary();
}

void SoLIDEDM4hep_source::Open() {
  m_reader = std::make_unique<podio::ROOTReader>();
  m_reader->openFile(GetResourceName());
  m_entries = m_reader->getEntries("events");
  m_current = 0;
  LOG_INFO(GetLogger()) << "Opened " << GetResourceName() << " with " << m_entries << " events";
}

void SoLIDEDM4hep_source::Close() { m_reader.reset(); }

JEventSource::Result SoLIDEDM4hep_source::Emit(JEvent& event) {
  if (m_current >= m_entries) return Result::FailureFinished;
  auto frame = podio::Frame(m_reader->readNextEntry("events"));
  ++m_current;
  const auto& geo = m_geo();

  // --- MC particles
  const auto& mcs = frame.get<edm4hep::MCParticleCollection>("MCParticles");
  std::vector<MCParticle*> mc_out;
  std::map<unsigned, int> id_to_index;
  int idx = 0;
  for (const auto& mc : mcs) {
    auto* p = new MCParticle();
    p->index = idx;
    p->pdg = mc.getPDG();
    p->status = mc.getGeneratorStatus();
    p->charge = mc.getCharge();
    p->mass = mc.getMass();
    const auto mom = mc.getMomentum();
    p->p4 = Vec4(mc.getEnergy(), Vec3(mom.x, mom.y, mom.z));
    const auto v = mc.getVertex();
    p->vertex = Vec3(v.x * 0.1, v.y * 0.1, v.z * 0.1);  // mm -> cm
    p->time = mc.getTime();
    if (mc.parents_size() > 0) {
      auto it = id_to_index.find(mc.getParents(0).getObjectID().index);
      p->parent = (it != id_to_index.end()) ? it->second : -1;
    }
    id_to_index[mc.getObjectID().index] = idx;
    mc_out.push_back(p);
    ++idx;
  }
  auto mcIndexOf = [&](const edm4hep::MCParticle& mc) -> int {
    if (!mc.isAvailable()) return -1;
    auto it = id_to_index.find(mc.getObjectID().index);
    return it == id_to_index.end() ? -1 : it->second;
  };

  // --- GEM hits: digitise sim hits (one hit per particle per plane), then add background
  std::vector<GEMHit> hits;
  std::vector<std::pair<int, Vec3>> crossings;
  std::normal_distribution<double> G(0.0, 1.0);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  {
    const auto& sim = frame.get<edm4hep::SimTrackerHitCollection>("GEMTrackerHits");
    std::map<std::pair<int, int>, std::pair<Vec3, double>> merged;  // (mc, plane) -> (edep-weighted pos, edep)
    for (const auto& h : sim) {
      const auto pos = h.getPosition();
      const Vec3 x(pos.x * 0.1, pos.y * 0.1, pos.z * 0.1);
      int plane = 0;
      for (const auto& g : geo.gemPlanes()) if (std::fabs(x.z - g.z) < 2.0) plane = g.id;
      if (plane == 0) continue;
      const int mc = mcIndexOf(h.getParticle());
      auto& m = merged[{mc, plane}];
      const double e = h.getEDep();
      m.first = m.first + x * e;
      m.second += e;
    }
    for (auto& kv : merged) {
      if (kv.second.second * 1e6 < m_edep_min()) continue;  // GeV -> keV
      const auto& g = geo.gemPlane(kv.first.second);
      Vec3 x = kv.second.first * (1.0 / kv.second.second);
      crossings.emplace_back(g.id, x);
      if (U(m_rng) > g.efficiency * m_bkg->survivalProbability(g.id, x.perp())) continue;
      GEMHit hit;
      hit.plane = g.id;
      hit.x = x.x + g.resolution * G(m_rng);
      hit.y = x.y + g.resolution * G(m_rng);
      hit.z = g.z;
      hit.sigma_x = hit.sigma_y = g.resolution;
      hit.edep = kv.second.second;
      hit.mc_index = kv.first.first;
      hits.push_back(hit);
      if (kv.first.first >= 0) mc_out[kv.first.first]->n_gem_hits++;
    }
    for (const auto& c : crossings) m_bkg->generateInRoad(c.first, c.second.x, c.second.y, m_bkg_road(), m_rng, hits);
  }

  // --- calorimeter clusters: greedy proximity clustering of the shower hits, energy-weighted position
  std::vector<ECCluster> clusters;
  auto clusterCalo = [&](const char* name, int calo, double zface) {
    if (!frame.getAvailableCollections().empty()) {
      const auto* collp = frame.get(name);
      if (!collp) return;
    }
    const auto& coll = frame.get<edm4hep::SimCalorimeterHitCollection>(name);
    struct Cl { double x = 0, y = 0, e = 0; std::map<int, double> mc; };
    std::vector<Cl> cls;
    for (const auto& h : coll) {
      const auto pos = h.getPosition();
      const double x = pos.x * 0.1, y = pos.y * 0.1, e = h.getEnergy();
      if (e <= 0) continue;
      Cl* best = nullptr;
      for (auto& c : cls) {
        if (std::hypot(c.x / c.e - x, c.y / c.e - y) < m_ec_cluster_radius()) { best = &c; break; }
      }
      if (!best) { cls.emplace_back(); best = &cls.back(); }
      best->x += x * e; best->y += y * e; best->e += e;
      for (unsigned i = 0; i < h.contributions_size(); ++i) {
        const auto contrib = h.getContributions(i);
        int mc = mcIndexOf(contrib.getParticle());
        // walk up to the primary
        while (mc >= 0 && mc_out[mc]->parent >= 0 && mc_out[mc_out[mc]->parent]->status != 4) mc = mc_out[mc]->parent;
        best->mc[mc] += contrib.getEnergy();
      }
    }
    for (const auto& c : cls) {
      if (c.e < 0.1) continue;  // 100 MeV threshold
      ECCluster ec;
      ec.calo = calo;
      ec.x = c.x / c.e + m_ec_sigma() * G(m_rng);
      ec.y = c.y / c.e + m_ec_sigma() * G(m_rng);
      ec.z = zface;
      ec.sigma_xy = m_ec_sigma();
      ec.energy = c.e;
      int best = -1; double ebest = 0;
      for (const auto& kv : c.mc) if (kv.second > ebest) { ebest = kv.second; best = kv.first; }
      ec.mc_index = best;
      clusters.push_back(ec);
      if (best >= 0) { if (calo == 0) mc_out[best]->hits_faec = true; else mc_out[best]->hits_laec = true; }
    }
  };
  clusterCalo("FAEC_ShHits", 0, geo.faec().z);
  clusterCalo("LAEC_ShHits", 1, geo.laec().z);

  // --- muon detector flags
  {
    const auto& mu = frame.get<edm4hep::SimTrackerHitCollection>("MuonHits");
    std::map<int, std::pair<int, bool>> layers;  // mc -> (deepest FA layer, LA/barrel hit)
    for (const auto& h : mu) {
      const int mc = mcIndexOf(h.getParticle());
      if (mc < 0) continue;
      const auto pos = h.getPosition();
      const double z = pos.z * 0.1, r = std::hypot(pos.x, pos.y) * 0.1;
      auto& l = layers[mc];
      if (r > 286.0) l.second = true;
      else {
        int layer = 0;
        for (size_t i = 0; i < geo.muon().fa_scint.size(); ++i)
          if (std::fabs(z - geo.muon().fa_scint[i].z) < 5.0) layer = static_cast<int>(i) + 1;
        l.first = std::max(l.first, layer);
      }
    }
    for (const auto& kv : layers) {
      if (kv.second.first == static_cast<int>(geo.muon().fa_scint.size())) mc_out[kv.first]->hits_fa_muon = true;
      if (kv.second.second) mc_out[kv.first]->hits_la_muon = true;
    }
  }

  std::vector<GEMHit*> hit_out;
  for (const auto& h : hits) hit_out.push_back(new GEMHit(h));
  std::vector<ECCluster*> ec_out;
  for (const auto& c : clusters) ec_out.push_back(new ECCluster(c));
  auto* gk = new GenKinematics();
  gk->luminosity = m_luminosity();

  event.SetEventNumber(m_current);
  event.SetRunNumber(1);
  event.Insert(mc_out, "MCParticles");
  event.Insert(hit_out, "GEMHits");
  event.Insert(ec_out, "ECClusters");
  event.Insert(gk, "GenKinematics");
  return Result::Success;
}

}  // namespace solid
