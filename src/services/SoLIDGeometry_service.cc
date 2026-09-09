#include "SoLIDGeometry_service.h"
#include <sstream>

namespace solid {

SoLIDGeometry_service::SoLIDGeometry_service() {
  SetTypeName(NAME_OF_THIS);
}

void SoLIDGeometry_service::Init() {
  m_beam_energy_val = m_beam_energy();
  m_target.z_center = m_target_z();
  m_target.half_length = 0.5 * m_target_length();

  const double res = m_gem_resolution() * 1e-4;  // um -> cm
  const double x0 = m_gem_x0();
  const double eff = m_gem_efficiency();

  // z positions and inner radii from solid_gemc geometry/gem_moved parameters (long endcap).
  // Outer radii: DDVCS variant (solid_SIDIS_DDVCS_gem) extends planes 1-4 to 135 cm.
  const bool ddvcs = (m_variant() != "sidis");
  m_gem = {
      {1, -175.0, 39.0, ddvcs ? 135.0 : 87.0, x0, res, eff, 30},
      {2, -150.0, 21.0, ddvcs ? 135.0 : 98.0, x0, res, eff, 24},
      {3, -119.0, 25.0, ddvcs ? 135.0 : 112.0, x0, res, eff, 21},
      {4, -68.0, 32.0, 135.0, x0, res, eff, 21},
      {5, 5.0, 42.0, 99.0, x0, res, eff, 29},
      {6, 92.0, 53.0, 122.0, x0, res, eff, 25},
  };

  // Forward-angle muon detector (prototype/muon_pro/solid_DDVCS_muon_forwardangle_geometry.pl):
  // iron slabs of 36 cm centred at z = 638, 702, 766 cm; scintillator annuli R = 80-285 cm
  // located 23 cm downstream of each slab centre.
  m_muon.fa_scint = {{661.0, 80.0, 285.0}, {725.0, 80.0, 285.0}, {789.0, 80.0, 285.0}};
  m_muon.fa_iron_total = 3 * 36.0;

  LOG_INFO(GetLogger()) << describe();
}

std::string SoLIDGeometry_service::describe() const {
  std::ostringstream os;
  os << "SoLID DDVCS geometry: E_beam=" << m_beam_energy_val << " GeV, target z=" << m_target.z_center
     << " cm (L=" << 2 * m_target.half_length << " cm LH2)\n";
  for (const auto& g : m_gem) {
    os << "  GEM" << g.id << ": z=" << g.z << " cm  R=[" << g.rin << "," << g.rout << "] cm  X0=" << g.x0_fraction
       << "  sigma=" << g.resolution * 1e4 << " um  eff=" << g.efficiency << "\n";
  }
  os << "  LAEC front z=" << m_laec.z << " R=[" << m_laec.rin << "," << m_laec.rout << "]\n";
  os << "  FAEC front z=" << m_faec.z << " R=[" << m_faec.rin << "," << m_faec.rout << "]\n";
  os << "  FA muon scint z=" << m_muon.fa_scint.back().z << " R=[" << m_muon.fa_scint.back().rin << ","
     << m_muon.fa_scint.back().rout << "] behind " << m_muon.fa_iron_total << " cm Fe";
  return os.str();
}

}  // namespace solid
