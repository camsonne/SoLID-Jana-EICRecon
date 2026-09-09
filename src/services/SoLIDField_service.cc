#include "SoLIDField_service.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <JANA/JException.h>

namespace solid {

SoLIDField_service::SoLIDField_service() { SetTypeName(NAME_OF_THIS); }

void SoLIDField_service::Init() {
  m_uniform = (m_mode() == "uniform");
  m_uniform_val = m_uniform_bz();
  m_scale_val = m_scale();
  if (m_mode() == "uniform") {
    LOG_INFO(GetLogger()) << "Uniform solenoid field Bz=" << m_uniform_bz() * m_scale() << " T";
    return;
  }
  if (m_mode() == "map") {
    if (!loadMap(m_map_file())) {
      throw JException("SoLIDField_service: cannot read field map '%s'", m_map_file().c_str());
    }
    return;
  }
  buildAnalyticGrid();
  double Br, Bz;
  BrBz(0, 0, Br, Bz);
  double Bt, Bzt;
  BrBz(0, -315.0, Bt, Bzt);
  LOG_INFO(GetLogger()) << "Analytic CLEO-II solenoid: Bz(0,0)=" << Bz << " T, Bz(0,target)=" << Bzt << " T";
}

// Complete elliptic integrals K(m), E(m) with parameter m = k^2, via Abramowitz & Stegun 17.3.34/36.
void SoLIDField_service::ellipticKE(double m, double& K, double& E) {
  double m1 = 1.0 - m;
  if (m1 < 1e-15) m1 = 1e-15;
  const double lnm1 = std::log(1.0 / m1);
  K = (1.38629436112 + m1 * (0.09666344259 + m1 * (0.03590092383 + m1 * (0.03742563713 + m1 * 0.01451196212)))) +
      (0.5 + m1 * (0.12498593597 + m1 * (0.06880248576 + m1 * (0.03328355346 + m1 * 0.00441787012)))) * lnm1;
  E = (1.0 + m1 * (0.44325141463 + m1 * (0.06260601220 + m1 * (0.04757383546 + m1 * 0.01736506451)))) +
      (m1 * (0.24998368310 + m1 * (0.09200180037 + m1 * (0.04069697526 + m1 * 0.00526449639)))) * lnm1;
}

// Field of a circular current loop (radius a [cm], at z0 [cm], current I [A]) at (r,z) [cm], in Tesla.
void SoLIDField_service::loopField(double a, double z0, double I, double r, double z, double& Br, double& Bz) {
  const double mu0_over_2pi = 2.0e-7;  // T m / A
  const double zz = (z - z0) * 0.01, aa = a * 0.01, rr = r * 0.01;  // to metres
  if (rr < 1e-9) {
    Br = 0;
    Bz = mu0_over_2pi * M_PI * I * aa * aa / std::pow(aa * aa + zz * zz, 1.5);
    return;
  }
  const double alpha2 = aa * aa + rr * rr + zz * zz - 2 * aa * rr;
  const double beta2 = aa * aa + rr * rr + zz * zz + 2 * aa * rr;
  const double beta = std::sqrt(beta2);
  const double k2 = 1.0 - alpha2 / beta2;
  double K, E;
  ellipticKE(k2, K, E);
  const double C = mu0_over_2pi * I;
  Bz = C / (alpha2 * beta) * ((aa * aa - rr * rr - zz * zz) * E + alpha2 * K);
  Br = C * zz / (alpha2 * beta * rr) * ((aa * aa + rr * rr + zz * zz) * E - alpha2 * K);
}

void SoLIDField_service::buildAnalyticGrid() {
  // Coil packs from CLEOv9.am (POISSON input): two layers, three z-segments each.
  struct Sheet { double a, z1, z2, jz; };  // radius, z-range [cm], current per length [A/cm]
  const double j_end = 542156.0 / 88.3, j_cen = 1009194.0 / 171.0;
  std::vector<Sheet> sheets = {
      {152.5, -173.8, -85.5, j_end}, {152.5, -85.5, 85.5, j_cen}, {152.5, 85.5, 173.8, j_end},
      {154.1, -173.8, -85.5, j_end}, {154.1, -85.5, 85.5, j_cen}, {154.1, 85.5, 173.8, j_end},
  };
  // discretise sheets into loops of 1 cm pitch
  struct Loop { double a, z0, I; };
  std::vector<Loop> loops;
  for (const auto& s : sheets) {
    int n = static_cast<int>(std::lround(s.z2 - s.z1));
    double dz = (s.z2 - s.z1) / n;
    for (int i = 0; i < n; ++i) loops.push_back({s.a, s.z1 + (i + 0.5) * dz, s.jz * dz});
  }
  // scale so that Bz(0,0) = central
  double Br0 = 0, Bz0 = 0;
  for (const auto& l : loops) {
    double br, bz;
    loopField(l.a, l.z0, l.I, 0, 0, br, bz);
    Bz0 += bz;
  }
  const double scale = m_central() / Bz0;
  (void)Br0;

  m_nr = static_cast<int>(std::lround((m_rmax - m_rmin) / m_dr)) + 1;
  m_nz = static_cast<int>(std::lround((m_zmax - m_zmin) / m_dz)) + 1;
  m_Br.assign(static_cast<size_t>(m_nr) * m_nz, 0.f);
  m_Bz.assign(static_cast<size_t>(m_nr) * m_nz, 0.f);
  for (int ir = 0; ir < m_nr; ++ir) {
    const double r = m_rmin + ir * m_dr;
    for (int iz = 0; iz < m_nz; ++iz) {
      const double z = m_zmin + iz * m_dz;
      double sBr = 0, sBz = 0;
      for (const auto& l : loops) {
        // skip the singular point on the sheet itself
        if (std::fabs(r - l.a) < 0.6 && std::fabs(z - l.z0) < 0.6) continue;
        double br, bz;
        loopField(l.a, l.z0, l.I, r, z, br, bz);
        sBr += br;
        sBz += bz;
      }
      m_Br[ir * m_nz + iz] = static_cast<float>(sBr * scale);
      m_Bz[ir * m_nz + iz] = static_cast<float>(sBz * scale);
    }
  }
  m_have_grid = true;
}

bool SoLIDField_service::loadMap(const std::string& file) {
  std::ifstream in(file);
  if (!in) return false;
  struct Row { double r, z, br, bz; };
  std::vector<Row> rows;
  std::string line;
  double rmin = 1e9, rmax = -1e9, zmin = 1e9, zmax = -1e9;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '<' || line[0] == '#') continue;
    std::istringstream iss(line);
    Row row;
    if (!(iss >> row.r >> row.z >> row.br >> row.bz)) continue;
    rows.push_back(row);
    rmin = std::min(rmin, row.r); rmax = std::max(rmax, row.r);
    zmin = std::min(zmin, row.z); zmax = std::max(zmax, row.z);
  }
  if (rows.size() < 4) return false;
  // infer grid steps from the two smallest distinct coordinates
  double dr = 1e9, dz = 1e9;
  for (const auto& row : rows) {
    if (row.r > rmin + 1e-6) dr = std::min(dr, row.r - rmin);
    if (row.z > zmin + 1e-6) dz = std::min(dz, row.z - zmin);
  }
  m_rmin = rmin; m_rmax = rmax; m_dr = dr;
  m_zmin = zmin; m_zmax = zmax; m_dz = dz;
  m_nr = static_cast<int>(std::lround((rmax - rmin) / dr)) + 1;
  m_nz = static_cast<int>(std::lround((zmax - zmin) / dz)) + 1;
  m_Br.assign(static_cast<size_t>(m_nr) * m_nz, 0.f);
  m_Bz.assign(static_cast<size_t>(m_nr) * m_nz, 0.f);
  const double gauss = 1e-4;
  for (const auto& row : rows) {
    int ir = static_cast<int>(std::lround((row.r - rmin) / dr));
    int iz = static_cast<int>(std::lround((row.z - zmin) / dz));
    if (ir < 0 || ir >= m_nr || iz < 0 || iz >= m_nz) continue;
    m_Br[ir * m_nz + iz] = static_cast<float>(row.br * gauss);
    m_Bz[ir * m_nz + iz] = static_cast<float>(row.bz * gauss);
  }
  m_have_grid = true;
  LOG_INFO(GetLogger()) << "Loaded field map " << file << ": r=[" << rmin << "," << rmax << "] dr=" << dr
                        << "  z=[" << zmin << "," << zmax << "] dz=" << dz << "  (" << rows.size() << " rows)";
  return true;
}

void SoLIDField_service::BrBz(double r, double z, double& Br, double& Bz) const {
  if (m_uniform) {
    Br = 0;
    Bz = m_uniform_val * m_scale_val;
    return;
  }
  if (!m_have_grid || r < m_rmin || r >= m_rmax || z < m_zmin || z >= m_zmax) {
    Br = Bz = 0;
    return;
  }
  const double fr = (r - m_rmin) / m_dr, fz = (z - m_zmin) / m_dz;
  int ir = static_cast<int>(fr), iz = static_cast<int>(fz);
  if (ir >= m_nr - 1) ir = m_nr - 2;
  if (iz >= m_nz - 1) iz = m_nz - 2;
  const double tr = fr - ir, tz = fz - iz;
  auto at = [&](const std::vector<float>& v, int i, int j) { return static_cast<double>(v[i * m_nz + j]); };
  auto bil = [&](const std::vector<float>& v) {
    return (1 - tr) * (1 - tz) * at(v, ir, iz) + tr * (1 - tz) * at(v, ir + 1, iz) + (1 - tr) * tz * at(v, ir, iz + 1) +
           tr * tz * at(v, ir + 1, iz + 1);
  };
  Br = bil(m_Br) * m_scale_val;
  Bz = bil(m_Bz) * m_scale_val;
}

Vec3 SoLIDField_service::B(const Vec3& pos) const {
  const double r = pos.perp();
  double Br, Bz;
  BrBz(r, pos.z, Br, Bz);
  if (r < 1e-9) return {0, 0, Bz};
  return {Br * pos.x / r, Br * pos.y / r, Bz};
}

}  // namespace solid
