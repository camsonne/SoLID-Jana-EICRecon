#include "DDVCSGenerator.h"
#include <cmath>

namespace solid {

namespace {
Vec3 rotateZ(const Vec3& v, double a) {
  return {v.x * std::cos(a) - v.y * std::sin(a), v.x * std::sin(a) + v.y * std::cos(a), v.z};
}
}  // namespace

double DDVCSGenerator::missingMass2(double Ebeam, const Vec4& ep, const Vec4& mup, const Vec4& mum) {
  Vec4 k(Ebeam, {0, 0, std::sqrt(Ebeam * Ebeam - kElectronMass * kElectronMass)});
  Vec4 P(kProtonMass, {0, 0, 0});
  Vec4 miss = k + P - ep - mup - mum;
  return miss.m2();
}

bool DDVCSGenerator::generate(std::vector<MCParticle>& particles, DDVCSKinematics& kin) {
  for (int attempt = 0; attempt < 10000; ++attempt) {
    if (tryGenerate(particles, kin)) return true;
  }
  return false;
}

bool DDVCSGenerator::tryGenerate(std::vector<MCParticle>& particles, DDVCSKinematics& kin) {
  std::uniform_real_distribution<double> U(0.0, 1.0);
  const double M = kProtonMass, E = m_cfg.beam_energy;

  // --- leptonic vertex
  const double Q2 = m_cfg.Q2_min * std::pow(m_cfg.Q2_max / m_cfg.Q2_min, U(m_rng));
  const double xB = m_cfg.xB_min + (m_cfg.xB_max - m_cfg.xB_min) * U(m_rng);
  const double nu = Q2 / (2 * M * xB);
  const double Ep = E - nu;
  if (Ep < m_cfg.Ep_min) return false;
  const double cos_te = 1.0 - Q2 / (2 * E * Ep);
  if (cos_te <= -1 || cos_te >= 1) return false;
  const double W2 = M * M + 2 * M * nu - Q2;
  if (W2 <= m_cfg.W_min * m_cfg.W_min) return false;
  const double W = std::sqrt(W2);

  // --- dilepton mass
  const double Qp2 = m_cfg.Qp2_min * std::pow(m_cfg.Qp2_max / m_cfg.Qp2_min, U(m_rng));
  const double Mll = std::sqrt(Qp2);
  if (W <= M + Mll) return false;

  // --- lab four-vectors of the leptons
  const double sin_te = std::sqrt(1 - cos_te * cos_te);
  const double phi_e = 2 * M_PI * U(m_rng);
  Vec4 k(E, {0, 0, E});
  Vec4 kp(Ep, rotateZ({Ep * sin_te, 0, Ep * cos_te}, phi_e));
  Vec4 P(M, {0, 0, 0});
  Vec4 q = k - kp;

  // --- gamma* p centre of mass
  Vec4 tot = q + P;
  Vec3 beta = tot.p * (1.0 / tot.E);
  Vec4 q_cm = q.boosted(beta * -1.0);
  Vec4 k_cm = k.boosted(beta * -1.0);
  const double Eq = q_cm.E, pq = q_cm.p.mag();
  const double Eqp = (W2 + Qp2 - M * M) / (2 * W);
  const double pqp2 = Eqp * Eqp - Qp2;
  if (pqp2 <= 0) return false;
  const double pqp = std::sqrt(pqp2);

  // t range
  const double tmin = -Q2 + Qp2 - 2 * (Eq * Eqp - pq * pqp);   // cos theta* = +1 (closest to zero)
  const double tmax = -Q2 + Qp2 - 2 * (Eq * Eqp + pq * pqp);   // cos theta* = -1
  // sample |t'| = tmin - t exponentially
  const double tp_max = std::min(tmin - tmax, m_cfg.tabs_max - std::fabs(tmin));
  if (tp_max <= 0) return false;
  const double b = m_cfg.t_slope;
  const double u = U(m_rng);
  const double tp = -std::log(1.0 - u * (1.0 - std::exp(-b * tp_max))) / b;
  const double t = tmin - tp;
  double cos_ts = (Eq * Eqp - (-Q2 + Qp2 - t) / 2.0) / (pq * pqp);
  if (cos_ts > 1) cos_ts = 1;
  if (cos_ts < -1) cos_ts = -1;
  const double sin_ts = std::sqrt(1 - cos_ts * cos_ts);
  const double phi = 2 * M_PI * U(m_rng);

  // basis in the CM: ez along q, ex in the leptonic plane
  Vec3 ez = q_cm.p.unit();
  Vec3 ex = (k_cm.p - ez * k_cm.p.dot(ez)).unit();
  Vec3 ey = ez.cross(ex);
  Vec3 qp_cm_p = (ex * (sin_ts * std::cos(phi)) + ey * (sin_ts * std::sin(phi)) + ez * cos_ts) * pqp;
  Vec4 qp_cm(Eqp, qp_cm_p);
  Vec4 Pp_cm(W - Eqp, qp_cm_p * -1.0);
  Vec4 qp = qp_cm.boosted(beta);
  Vec4 Pp = Pp_cm.boosted(beta);

  // --- gamma* -> mu+ mu- in its rest frame, 1 + cos^2 theta_l
  double cos_tl;
  do {
    cos_tl = 2 * U(m_rng) - 1;
  } while (U(m_rng) * 2.0 > 1.0 + cos_tl * cos_tl);
  const double sin_tl = std::sqrt(1 - cos_tl * cos_tl);
  const double phi_l = 2 * M_PI * U(m_rng);
  const double pl = std::sqrt(Qp2 / 4.0 - kMuonMass * kMuonMass);
  // decay axis: gamma* direction in the lab (helicity frame)
  Vec3 lz = qp.p.unit();
  Vec3 lx = std::fabs(lz.z) < 0.99 ? Vec3(0, 0, 1).cross(lz).unit() : Vec3(1, 0, 0).cross(lz).unit();
  Vec3 ly = lz.cross(lx);
  Vec3 lm_rest = (lx * (sin_tl * std::cos(phi_l)) + ly * (sin_tl * std::sin(phi_l)) + lz * cos_tl) * pl;
  Vec4 mum_rest(Mll / 2, lm_rest);
  Vec4 mup_rest(Mll / 2, lm_rest * -1.0);
  Vec3 beta_g = qp.p * (1.0 / qp.E);
  Vec4 mum = mum_rest.boosted(beta_g);
  Vec4 mup = mup_rest.boosted(beta_g);

  // --- vertex
  std::normal_distribution<double> G(0.0, 1.0);
  Vec3 vtx(m_cfg.beam_sigma_xy > 0 ? m_cfg.beam_sigma_xy * G(m_rng) : 0.0,
           m_cfg.beam_sigma_xy > 0 ? m_cfg.beam_sigma_xy * G(m_rng) : 0.0,
           m_cfg.target_z + m_cfg.target_half_length * (2 * U(m_rng) - 1));

  particles.clear();
  auto add = [&](int pdg, int status, int parent, double charge, double mass, const Vec4& p4) {
    MCParticle mc;
    mc.index = static_cast<int>(particles.size());
    mc.pdg = pdg; mc.status = status; mc.parent = parent; mc.charge = charge; mc.mass = mass;
    mc.p4 = p4; mc.vertex = vtx;
    particles.push_back(mc);
  };
  add(11, 4, -1, -1, kElectronMass, k);
  add(2212, 4, -1, +1, kProtonMass, P);
  add(11, 1, 0, -1, kElectronMass, kp);
  add(2212, 1, 1, +1, kProtonMass, Pp);
  add(22, 11, 0, 0, Mll, qp);        // virtual photon (intermediate)
  add(13, 1, 4, -1, kMuonMass, mum);
  add(-13, 1, 4, +1, kMuonMass, mup);

  kin.Q2 = Q2; kin.xB = xB; kin.nu = nu; kin.W = W; kin.Qp2 = Qp2; kin.t = t; kin.tmin = tmin;
  kin.phi_trento = phi; kin.y = nu / E; kin.cos_theta_l = cos_tl; kin.phi_l = phi_l;
  return true;
}

}  // namespace solid
