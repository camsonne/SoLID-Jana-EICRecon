#include "HepMC3Reader.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace solid {

namespace {
bool isChargedLepton(int pdg) { int a = std::abs(pdg); return a == 11 || a == 13 || a == 15; }
double leptonCharge(int pdg) { return pdg > 0 ? -1.0 : +1.0; }   // l- has positive PDG code
double chargeOf(int pdg) {
  if (isChargedLepton(pdg)) return leptonCharge(pdg);
  switch (std::abs(pdg)) {
    case 2212: case 211: case 321: return pdg > 0 ? 1.0 : -1.0;
    default: return 0.0;
  }
}
std::vector<int> parseList(const std::string& s) {   // "[1,2]" or "[1]"
  std::vector<int> out;
  std::string body = s;
  body.erase(std::remove(body.begin(), body.end(), '['), body.end());
  body.erase(std::remove(body.begin(), body.end(), ']'), body.end());
  std::stringstream ss(body);
  std::string item;
  while (std::getline(ss, item, ',')) if (!item.empty()) out.push_back(std::stoi(item));
  return out;
}
}  // namespace

HepMC3Reader::HepMC3Reader(const std::string& filename) : m_in(filename) {}

void HepMC3Reader::parseRunLine(const std::string& line) {
  std::istringstream is(line);
  std::string tag;
  is >> tag;
  if (tag == "T") {
    std::string rest;
    std::getline(is, rest);
    const auto b = rest.find_first_not_of(' ');
    rest = b == std::string::npos ? "" : rest.substr(b);
    if (m_run.tool.find(rest) != std::string::npos) return;   // run information repeated at the end of the file
    if (!m_run.tool.empty()) m_run.tool += "; ";
    m_run.tool += rest;
  } else if (tag == "A") {
    std::string name, value;
    is >> name;
    std::getline(is, value);
    const auto b = value.find_first_not_of(' ');
    value = b == std::string::npos ? "" : value.substr(b);
    m_run.attributes[name] = value;
    try {
      if (name == "integrated_cross_section_value") { m_run.xsec = std::stod(value); m_run.has_xsec = true; }
      if (name == "integrated_cross_section_uncertainty") m_run.xsec_err = std::stod(value);
    } catch (...) {}
  }
}

bool HepMC3Reader::readEvent() {
  m_parts.clear();
  m_vertices.clear();
  m_weight = 1.0;
  m_event_xsec = 0;
  m_vertex_in_file = false;
  std::string line;
  bool in_event = false;
  if (m_have_pending) { line = m_pending; m_have_pending = false; in_event = true; }
  auto parseE = [&](const std::string& l) {
    std::istringstream is(l);
    std::string tag;
    is >> tag >> m_event_number;
  };
  if (in_event) parseE(line);
  while (std::getline(m_in, line)) {
    if (line.empty()) continue;
    const char c = line[0];
    if (c == 'E') {
      if (in_event) { m_pending = line; m_have_pending = true; return true; }
      in_event = true;
      parseE(line);
      continue;
    }
    if (!in_event) {               // header / run information
      parseRunLine(line);
      continue;
    }
    std::istringstream is(line);
    std::string tag;
    is >> tag;
    if (tag == "U") {
      std::string mu, lu;
      is >> mu >> lu;
      m_mom_unit = (mu == "MEV") ? 1e-3 : 1.0;
      m_len_unit = (lu == "CM") ? 1.0 : 0.1;
    } else if (tag == "W") {
      // event weight; a "W <names>" line is run information (EpIC rewrites it after the last event)
      std::string w;
      if (is >> w) { try { m_weight = std::stod(w); } catch (...) {} }
    } else if (tag == "P") {
      RawParticle p;
      is >> p.id >> p.parent >> p.pdg >> p.px >> p.py >> p.pz >> p.e >> p.m >> p.status;
      m_parts.push_back(p);
    } else if (tag == "V") {
      int id = 0, status = 0;
      std::string list;
      is >> id >> status >> list;
      RawVertex v;
      v.in = parseList(list);
      std::string at;
      if (is >> at && at == "@") {
        double t = 0;
        is >> v.x >> v.y >> v.z >> t;
        v.has_pos = true;
      }
      m_vertices[id] = v;
    } else if (tag == "A" || tag == "T" || tag == "N") {
      // event attributes (e.g. GenCrossSection) or run-info lines written mid-file
      if (tag != "A") parseRunLine(line);
      else {
        // "A <id> <name> <value>" for event attributes; "A <name> <value>" for run attributes
        std::string second;
        is >> second;
        if (!second.empty() && !std::isdigit(static_cast<unsigned char>(second[0])) && second[0] != '-') {
          parseRunLine(line);
        } else {
          std::string name;
          double value = 0;
          if (is >> name >> value && name == "GenCrossSection") m_event_xsec = value;
        }
      }
    } else if (line.rfind("HepMC::", 0) == 0) {
      // end-of-file marker
    }
  }
  return in_event;
}

bool HepMC3Reader::next(std::vector<MCParticle>& particles) {
  if (!readEvent()) return false;

  // incoming particle ids of the production vertex of every particle
  auto incoming = [&](const RawParticle& p) -> std::vector<int> {
    if (p.parent > 0) return {p.parent};
    auto it = m_vertices.find(p.parent);
    return it == m_vertices.end() ? std::vector<int>{} : it->second.in;
  };
  auto byId = [&](int id) -> const RawParticle* {
    for (const auto& p : m_parts) if (p.id == id) return &p;
    return nullptr;
  };
  auto has = [](const std::vector<int>& v, int id) { return std::find(v.begin(), v.end(), id) != v.end(); };

  const RawParticle* beam = nullptr; const RawParticle* target = nullptr;
  for (const auto& p : m_parts) {
    if (p.status != 4) continue;
    if (isChargedLepton(p.pdg)) beam = &p; else target = &p;
  }
  if (!beam || !target) throw std::runtime_error("HepMC3Reader: event " + std::to_string(m_event_number) +
                                                 " has no lepton beam and target (status 4)");
  const RawParticle* scat = nullptr; const RawParticle* recoil = nullptr;
  const RawParticle* lm = nullptr; const RawParticle* lp = nullptr;
  std::vector<const RawParticle*> extra;
  for (const auto& p : m_parts) {
    if (p.status != 1) continue;
    if (p.pdg == target->pdg && !recoil) { recoil = &p; continue; }
    if (isChargedLepton(p.pdg)) {
      const auto in = incoming(p);
      const bool from_beam = has(in, beam->id);
      if (p.pdg == beam->pdg && (from_beam || in.empty()) && !scat) { scat = &p; continue; }
      if (!from_beam) {
        if (p.pdg > 0 && !lm) { lm = &p; continue; }
        if (p.pdg < 0 && !lp) { lp = &p; continue; }
      }
    }
    extra.push_back(&p);
  }
  if (!scat || !recoil || !lm || !lp)
    throw std::runtime_error("HepMC3Reader: event " + std::to_string(m_event_number) +
                             " does not have the DDVCS topology (l' p' l- l+)");
  // gamma* : parent of the lepton pair when present, otherwise the pair sum
  const RawParticle* gstar = nullptr;
  {
    const auto in = incoming(*lm);
    if (in.size() == 1) gstar = byId(in[0]);
  }

  // frame: beam along +z
  m_flipped = beam->pz < 0;
  const double sx = m_flipped ? -1.0 : 1.0, sz = m_flipped ? -1.0 : 1.0;  // (+ 0.0 below avoids writing -0)
  auto p4 = [&](const RawParticle& r) {
    return Vec4(r.e * m_mom_unit, {sx * r.px * m_mom_unit + 0.0, r.py * m_mom_unit, sz * r.pz * m_mom_unit + 0.0});
  };

  // vertex: production vertex of the scattered lepton, if the file has positions
  Vec3 vtx;
  if (scat->parent < 0) {
    auto it = m_vertices.find(scat->parent);
    if (it != m_vertices.end() && it->second.has_pos) {
      m_vertex_in_file = true;
      vtx = {sx * it->second.x * m_len_unit, it->second.y * m_len_unit, sz * it->second.z * m_len_unit};
    }
  }

  particles.clear();
  auto add = [&](const RawParticle* r, int pdg, int status, int parent, const Vec4& v) {
    MCParticle mc;
    mc.index = static_cast<int>(particles.size());
    mc.pdg = pdg;
    mc.status = status;
    mc.parent = parent;
    mc.charge = chargeOf(pdg);
    mc.p4 = v;
    mc.mass = r ? r->m * m_mom_unit : std::sqrt(std::max(0.0, v.m2()));
    mc.vertex = vtx;
    particles.push_back(mc);
  };
  add(beam, beam->pdg, 4, -1, p4(*beam));
  add(target, target->pdg, 4, -1, p4(*target));
  add(scat, scat->pdg, 1, 0, p4(*scat));
  add(recoil, recoil->pdg, 1, 1, p4(*recoil));
  const Vec4 pair = p4(*lm) + p4(*lp);
  add(nullptr, 22, 11, 0, gstar ? p4(*gstar) : pair);
  add(lm, lm->pdg, 1, 4, p4(*lm));
  add(lp, lp->pdg, 1, 4, p4(*lp));
  for (const auto* r : extra) add(r, r->pdg, 1, -1, p4(*r));
  return true;
}

}  // namespace solid
