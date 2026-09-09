// DDVCS kinematics generator:  e(k) + p(P) -> e'(k') + p'(P') + gamma*(q') ,  gamma* -> mu+ mu-
//
// This is a phase-space style generator with importance sampling in the SoLID DDVCS region
// (11 GeV beam, LH2 target). It is intended for detector-response studies (acceptance,
// tracking resolution, missing-mass resolution) and does not implement the DDVCS/BH
// amplitude; every event carries weight 1.  Sampling densities:
//    Q^2   log-uniform in [Q2_min, Q2_max]
//    x_B   uniform     in [xB_min, xB_max]
//    Q'^2  log-uniform in [Qp2_min, Qp2_max]      (dilepton mass squared)
//    t     t_min - |t'| with |t'| ~ exp(-b |t'|), truncated at |t| < tabs_max and physical t_max
//    phi   (Trento) uniform, dilepton decay 1 + cos^2(theta_l) in the gamma* rest frame.
#pragma once

#include <random>
#include <vector>
#include "datamodel/SoLIDDataModel.h"

namespace solid {

struct DDVCSConfig {
  double beam_energy = 11.0;   // GeV
  double Q2_min = 1.0, Q2_max = 6.0;
  double xB_min = 0.10, xB_max = 0.55;
  double Qp2_min = 2.0, Qp2_max = 9.0;   // M_mumu between ~1.4 and 3.0 GeV
  double t_slope = 1.5;                  // GeV^-2
  double tabs_max = 1.5;                 // GeV^2
  double W_min = 2.0;                    // GeV
  double target_z = -315.0, target_half_length = 7.5;  // cm
  double beam_sigma_xy = 0.0;            // cm (rastered beam position is assumed known)
  double Ep_min = 0.5;                   // GeV, scattered electron energy floor
};

struct DDVCSKinematics {
  double Q2 = 0, xB = 0, nu = 0, W = 0, Qp2 = 0, t = 0, tmin = 0, phi_trento = 0, y = 0;
  double cos_theta_l = 0, phi_l = 0;
};

class DDVCSGenerator {
public:
  explicit DDVCSGenerator(const DDVCSConfig& cfg, unsigned seed = 12345) : m_cfg(cfg), m_rng(seed) {}

  /// Generate one event. Fills particles (beam e-, target p, e', p', gamma*, mu-, mu+) and kinematics.
  /// Returns false only if no physical configuration was found after many attempts.
  bool generate(std::vector<MCParticle>& particles, DDVCSKinematics& kin);

  const DDVCSConfig& config() const { return m_cfg; }
  std::mt19937& rng() { return m_rng; }

  /// Missing mass squared of (k + P - k' - l+ - l-) for the given four-vectors.
  static double missingMass2(double Ebeam, const Vec4& ep, const Vec4& mup, const Vec4& mum);

private:
  bool tryGenerate(std::vector<MCParticle>& particles, DDVCSKinematics& kin);
  DDVCSConfig m_cfg;
  std::mt19937 m_rng;
};

}  // namespace solid
