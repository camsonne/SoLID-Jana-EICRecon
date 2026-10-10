// DDVCS event kinematics from four-vectors:  e(k) + p(P) -> e'(k') + p'(P') + gamma*(q'),  gamma* -> l- l+
//
// Used for the generator truth (internal generator or EpIC input) and for the reconstructed
// events, so that truth and reconstruction follow exactly the same conventions:
//   Q2 = -q^2 with q = k - k',  xB = Q2 / (2 P.q),  y = P.q / P.k,  W^2 = (q + P)^2
//   Q'2 = (l- + l+)^2,  t = (q - q')^2,  t_min at cos(theta*) = 1 in the gamma* p centre of mass
//   phi      Trento angle, in the target rest frame, between the lepton plane (k, k') and the
//            production plane (q, q'):
//              cos(phi) ~ (q x k).(q x q'),  sin(phi) ~ (k x q').q_hat,  in [0, 2 pi)
//   theta_l  angle between l- and the recoil-proton direction in the dilepton rest frame
//   phi_l    azimuth of l- around that axis, measured from n x z' with n = P x P' (both in the
//            dilepton rest frame) and z' the recoil-proton direction; in [0, 2 pi)
// phi, theta_l and phi_l are the angles EpIC (DDVCSKinematicDefault) uses to build the event, so
// EpIC truth read back with these functions reproduces its generated phi, thetaL and phiL
// (verified by tests/test_hepmc_reader.cc).
#pragma once

#include <cmath>
#include "datamodel/SoLIDDataModel.h"

namespace solid {

struct DDVCSKinematics {
  double Q2 = 0, xB = 0, nu = 0, W = 0, Qp2 = 0, t = 0, tmin = 0, phi_trento = 0, y = 0;
  double cos_theta_l = 0, phi_l = 0;
  double mm2 = 0;  // (k + P - k' - l- - l+)^2
};

inline double wrapTwoPi(double a) {
  a = std::fmod(a, 2 * M_PI);
  return a < 0 ? a + 2 * M_PI : a;
}

/// All kinematic variables from the beam (k), target (P), scattered lepton (k') and the lepton
/// pair (lm = negative lepton, lp = positive lepton). The recoil proton is taken as the missing
/// four-momentum, so the same function serves truth and reconstruction.
inline DDVCSKinematics computeDDVCSKinematics(const Vec4& k, const Vec4& P, const Vec4& kp, const Vec4& lm,
                                              const Vec4& lp) {
  DDVCSKinematics r;
  const Vec4 q = k - kp;
  const Vec4 qp = lm + lp;
  const Vec4 Pp = k + P - kp - lm - lp;
  r.mm2 = Pp.m2();
  r.Q2 = -q.m2();
  const double Pq = P.dot(q);
  r.xB = Pq != 0 ? r.Q2 / (2 * Pq) : 0;
  const double M = std::sqrt(std::max(0.0, P.m2()));
  r.nu = M > 0 ? Pq / M : 0;
  r.y = P.dot(k) != 0 ? Pq / P.dot(k) : 0;
  const Vec4 tot = q + P;
  const double W2 = tot.m2();
  r.W = W2 > 0 ? std::sqrt(W2) : 0;
  r.Qp2 = qp.m2();
  r.t = (q - qp).m2();

  // t_min in the gamma* p centre of mass (recoil mass = target mass)
  if (W2 > 0) {
    const Vec3 beta = tot.p * (-1.0 / tot.E);
    const Vec4 q_cm = q.boosted(beta);
    const double Eq = q_cm.E, pq = q_cm.p.mag();
    const double Eqp = (W2 + r.Qp2 - M * M) / (2 * r.W);
    const double pqp2 = Eqp * Eqp - r.Qp2;
    if (pqp2 > 0) r.tmin = -r.Q2 + r.Qp2 - 2 * (Eq * Eqp - pq * std::sqrt(pqp2));
  }

  // Trento phi, in the target rest frame
  {
    const Vec3 bP = P.p * (-1.0 / P.E);
    const Vec3 kr = k.boosted(bP).p, qr = q.boosted(bP).p, qpr = qp.boosted(bP).p;
    const Vec3 qh = qr.unit();
    const Vec3 a = qh.cross(kr), b = qh.cross(qpr);
    const double na = a.mag(), nb = b.mag();
    if (na > 0 && nb > 0) {
      const double c = a.dot(b) / (na * nb);
      const double s = kr.cross(qpr).dot(qh) / (na * nb);
      r.phi_trento = wrapTwoPi(std::atan2(s, c));
    }
  }

  // lepton angles in the dilepton rest frame
  if (r.Qp2 > 0) {
    const Vec3 beta = qp.p * (-1.0 / qp.E);
    const Vec3 l = lm.boosted(beta).p;
    const Vec3 zp = Pp.boosted(beta).p.unit();
    const Vec3 n = P.boosted(beta).p.cross(zp).unit();
    const Vec3 xp = n.cross(zp);
    const double lmag = l.mag();
    if (lmag > 0) r.cos_theta_l = l.dot(zp) / lmag;
    r.phi_l = wrapTwoPi(std::atan2(l.dot(n), l.dot(xp)));
  }
  return r;
}

}  // namespace solid
