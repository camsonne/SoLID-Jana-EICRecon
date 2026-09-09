#include "TrackPropagator.h"
#include <cmath>

namespace solid {

namespace {
constexpr double kLight = 0.00299792458;  // GeV / (T cm) : dp/ds = q * kLight * (u x B)
}

double highlandTheta0(double p, double beta, double x_over_x0) {
  if (x_over_x0 <= 0 || p <= 0) return 0;
  double t = 0.0136 / (beta * p) * std::sqrt(x_over_x0) * (1.0 + 0.038 * std::log(x_over_x0));
  return t > 0 ? t : 0;
}

void TrackPropagator::transverseBasis(const Vec3& u, Vec3& e1, Vec3& e2) {
  Vec3 ref = std::fabs(u.z) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
  e1 = u.cross(ref).unit();
  e2 = u.cross(e1).unit();
}

void TrackPropagator::scatter(TrackState& s, double theta0, std::mt19937& rng) {
  if (theta0 <= 0) return;
  std::normal_distribution<double> g(0.0, theta0);
  Vec3 u = s.dir(), e1, e2;
  transverseBasis(u, e1, e2);
  double a1 = g(rng), a2 = g(rng);
  Vec3 nu = (u + e1 * a1 + e2 * a2).unit();
  s.mom = nu * s.mom.mag();
}

void TrackPropagator::rk4Step(TrackState& s, double h) const {
  // state y = (pos, mom); dy/ds = (u, q k (u x B))
  const double q = s.charge;
  auto deriv = [&](const Vec3& pos, const Vec3& mom, Vec3& dpos, Vec3& dmom) {
    Vec3 u = mom.unit();
    Vec3 B = m_field->B(pos);
    dpos = u;
    dmom = u.cross(B) * (q * kLight);
  };
  Vec3 k1p, k1m, k2p, k2m, k3p, k3m, k4p, k4m;
  deriv(s.pos, s.mom, k1p, k1m);
  deriv(s.pos + k1p * (h / 2), s.mom + k1m * (h / 2), k2p, k2m);
  deriv(s.pos + k2p * (h / 2), s.mom + k2m * (h / 2), k3p, k3m);
  deriv(s.pos + k3p * h, s.mom + k3m * h, k4p, k4m);
  s.pos = s.pos + (k1p + k2p * 2 + k3p * 2 + k4p) * (h / 6);
  s.mom = s.mom + (k1m + k2m * 2 + k3m * 2 + k4m) * (h / 6);
  s.path += h;
}

bool TrackPropagator::propagateToZ(TrackState& s, double z_target, std::mt19937* rng, double air_x0,
                                   double rmax) const {
  if (z_target < s.pos.z) return false;
  const double p = s.mom.mag();
  int guard = 0;
  while (s.pos.z < z_target - 1e-7) {
    if (++guard > 100000) return false;
    Vec3 u = s.dir();
    if (u.z <= 1e-6) return false;
    if (s.pos.perp() > rmax) return false;
    double remaining = (z_target - s.pos.z) / u.z;  // path length estimate
    double h = std::min(m_step, remaining);
    if (h < 1e-9) break;
    rk4Step(s, h);
    if (rng && air_x0 > 0) scatter(s, highlandTheta0(p, 1.0, h / air_x0), *rng);
  }
  // final exact landing on the plane (linear correction of the tiny overshoot / undershoot)
  Vec3 u = s.dir();
  if (u.z > 1e-9) {
    double dz = z_target - s.pos.z;
    s.pos = s.pos + u * (dz / u.z);
    s.path += dz / u.z;
  }
  return true;
}

bool TrackPropagator::propagateToR(TrackState& s, double r_target, double z_stop, double rmax) const {
  int guard = 0;
  while (s.pos.perp() < r_target) {
    if (++guard > 100000) return false;
    if (s.pos.z > z_stop || s.pos.perp() > rmax) return false;
    if (s.dir().z <= 1e-6 && s.pos.perp() < 1e-3) return false;
    rk4Step(s, m_step);
  }
  return true;
}

}  // namespace solid
