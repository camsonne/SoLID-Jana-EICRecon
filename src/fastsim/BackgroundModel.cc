#include "BackgroundModel.h"
#include <cmath>
#include <fstream>
#include <sstream>

namespace solid {

BackgroundModel::BackgroundModel(const BackgroundConfig& cfg, const SoLIDGeometry_service* geo) : m_cfg(cfg), m_geo(geo) {
  if (!m_cfg.file.empty()) {
    std::ifstream in(m_cfg.file);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream iss(line);
      int plane; double rho0, alpha;
      if (iss >> plane >> rho0 >> alpha && plane >= 1 && plane <= 6) {
        m_cfg.rho0_kHz_mm2[plane - 1] = rho0;
        m_cfg.alpha[plane - 1] = alpha;
      }
    }
  }
}

double BackgroundModel::rateDensity(int plane, double R) const {
  const auto& g = m_geo->gemPlane(plane);
  const double Rc = std::max(R, g.rin);
  const double scale = m_cfg.luminosity / m_cfg.ref_luminosity;
  return m_cfg.rho0_kHz_mm2[plane - 1] * 1e3 * std::pow(m_cfg.R0_cm / Rc, m_cfg.alpha[plane - 1]) * scale;
}

double BackgroundModel::hitDensity(int plane, double R) const {
  return rateDensity(plane, R) * m_cfg.time_window_ns * 1e-9 * 100.0;  // Hz/mm^2 * s * (100 mm^2/cm^2)
}

double BackgroundModel::moduleArea(int plane) const {
  const auto& g = m_geo->gemPlane(plane);
  return M_PI * (g.rout * g.rout - g.rin * g.rin) / g.nmodules;
}

double BackgroundModel::stripLength(int plane) const {
  if (m_cfg.strip_length_cm > 0) return m_cfg.strip_length_cm;
  const auto& g = m_geo->gemPlane(plane);
  return std::min(100.0, g.rout - g.rin);
}

double BackgroundModel::occupancy(int plane, double R) const {
  const double n_per_mm2 = hitDensity(plane, R) / 100.0;
  if (m_cfg.readout == "pixel") return n_per_mm2 * m_cfg.pad_area_mm2;
  return n_per_mm2 * m_cfg.strip_pitch_mm * stripLength(plane) * 10.0;
}

double BackgroundModel::survivalProbability(int plane, double R) const {
  return std::exp(-2.0 * occupancy(plane, R));
}

double BackgroundModel::ghostDensity(int plane, double R) const {
  if (m_cfg.readout == "pixel") return 0;
  const double rho = hitDensity(plane, R);
  return rho * rho * moduleArea(plane);
}

void BackgroundModel::generateInRoad(int plane, double x, double y, double w, std::mt19937& rng,
                                     std::vector<GEMHit>& out) const {
  const double R = std::hypot(x, y);
  const auto& g = m_geo->gemPlane(plane);
  const double area = 4 * w * w;
  const double n_real = hitDensity(plane, R) * area;
  const double n_ghost = std::min(ghostDensity(plane, R) * area, 500.0);
  std::uniform_real_distribution<double> U(-w, w);
  auto make = [&](int n_mean, bool ghost) {
    std::poisson_distribution<int> P(n_mean);
    int n = P(rng);
    for (int i = 0; i < n; ++i) {
      GEMHit h;
      h.plane = plane;
      h.x = x + U(rng);
      h.y = y + U(rng);
      h.z = g.z;
      h.sigma_x = h.sigma_y = g.resolution;
      h.mc_index = -1;
      h.ghost = ghost;
      const double r = h.r();
      if (r < g.rin || r > g.rout) continue;
      out.push_back(h);
    }
  };
  make(n_real, false);
  make(n_ghost, true);
}

std::string BackgroundModel::summary() const {
  std::ostringstream os;
  os << "GEM background: L=" << m_cfg.luminosity << " cm^-2 s^-1, window=" << m_cfg.time_window_ns << " ns, readout="
     << m_cfg.readout << "\n";
  for (int pl = 1; pl <= 6; ++pl) {
    const auto& g = m_geo->gemPlane(pl);
    os << "  GEM" << pl << " @Rin: rate=" << rateDensity(pl, g.rin) / 1e3 << " kHz/mm^2, hits=" << hitDensity(pl, g.rin)
       << "/cm^2, occupancy=" << occupancy(pl, g.rin) << ", survival=" << survivalProbability(pl, g.rin)
       << ", ghosts=" << ghostDensity(pl, g.rin) << "/cm^2\n";
  }
  return os.str();
}

}  // namespace solid
