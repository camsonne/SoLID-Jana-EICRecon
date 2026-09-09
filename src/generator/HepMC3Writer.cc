#include "HepMC3Writer.h"
#include <iomanip>

namespace solid {

HepMC3Writer::HepMC3Writer(const std::string& filename) : m_out(filename) {
  m_out << "HepMC::Version 3.02.05\nHepMC::Asciiv3-START\n";
  m_out << std::setprecision(10);
}

HepMC3Writer::~HepMC3Writer() {
  if (m_out.is_open()) {
    m_out << "HepMC::Asciiv3-END\n";
    m_out.close();
  }
}

void HepMC3Writer::write(int event_number, const std::vector<MCParticle>& p, double weight) {
  // vertices: -1 primary (in: 1,2 ; out: 3,4,5), -2 gamma* decay (in: 5 ; out: 6,7)
  m_out << "E " << event_number << " 2 " << p.size() << "\n";
  m_out << "U GEV MM\n";
  m_out << "W " << weight << "\n";
  auto P = [&](int id, int vtx, const MCParticle& mc) {
    m_out << "P " << id << " " << vtx << " " << mc.pdg << " " << mc.p4.p.x << " " << mc.p4.p.y << " " << mc.p4.p.z
          << " " << mc.p4.E << " " << mc.mass << " " << mc.status << "\n";
  };
  P(1, 0, p[0]);
  P(2, 0, p[1]);
  const Vec3& v = p[2].vertex;  // cm -> mm
  m_out << "V -1 0 [1,2] @ " << v.x * 10 << " " << v.y * 10 << " " << v.z * 10 << " 0\n";
  P(3, -1, p[2]);
  P(4, -1, p[3]);
  P(5, -1, p[4]);
  m_out << "V -2 0 [5]\n";
  P(6, -2, p[5]);
  P(7, -2, p[6]);
}

}  // namespace solid
