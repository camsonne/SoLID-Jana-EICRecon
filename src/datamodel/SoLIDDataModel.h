// SoLID lightweight data model for the JANA2 reconstruction chain.
// These are plain JObjects so that the chain can run without PODIO/EDM4hep.
// When ddsim (DD4hep/Geant4) output is available, SoLIDEDM4hep_source converts
// edm4hep::MCParticle / edm4hep::SimTrackerHit into these same objects.
#pragma once

#include <JANA/JObject.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace solid {

// Units used throughout: cm, GeV, ns, Tesla.

struct Vec3 {
  double x = 0, y = 0, z = 0;
  Vec3() = default;
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  double mag() const { return std::sqrt(x * x + y * y + z * z); }
  double perp() const { return std::sqrt(x * x + y * y); }
  double theta() const { return std::atan2(perp(), z); }
  double phi() const { return std::atan2(y, x); }
  Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
  Vec3 cross(const Vec3& o) const { return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x}; }
  Vec3 unit() const { double m = mag(); return m > 0 ? (*this) * (1.0 / m) : *this; }
};

struct Vec4 {
  double E = 0;
  Vec3 p;
  Vec4() = default;
  Vec4(double E_, const Vec3& p_) : E(E_), p(p_) {}
  static Vec4 fromPM(const Vec3& p_, double m) { return Vec4(std::sqrt(p_.mag() * p_.mag() + m * m), p_); }
  double m2() const { return E * E - p.dot(p); }
  double m() const { double v = m2(); return v >= 0 ? std::sqrt(v) : -std::sqrt(-v); }
  Vec4 operator+(const Vec4& o) const { return {E + o.E, p + o.p}; }
  Vec4 operator-(const Vec4& o) const { return {E - o.E, p - o.p}; }
  double dot(const Vec4& o) const { return E * o.E - p.dot(o.p); }
  // Boost this four-vector by velocity beta (of the new frame relative to the current one, i.e.
  // to go to the rest frame of a particle with momentum P use beta = -P.p/P.E).
  Vec4 boosted(const Vec3& beta) const {
    double b2 = beta.dot(beta);
    if (b2 <= 0) return *this;
    double gamma = 1.0 / std::sqrt(1.0 - b2);
    double bp = beta.dot(p);
    double gamma2 = (gamma - 1.0) / b2;
    Vec3 pp = p + beta * (gamma2 * bp + gamma * E);
    return {gamma * (E + bp), pp};
  }
};

/// Generated (truth) particle.
struct MCParticle : public JObject {
  JOBJECT_PUBLIC(MCParticle)
  int index = -1;          // position in the event record
  int pdg = 0;
  int status = 1;          // 1 = final state, 4 = beam, 11 = intermediate (decayed)
  int parent = -1;         // index of parent (-1 for beam/none)
  double charge = 0;       // in units of e
  double mass = 0;         // GeV
  Vec4 p4;                 // GeV
  Vec3 vertex;             // cm
  double time = 0;         // ns
  // Fast-simulation truth bookkeeping (filled by SoLIDFastSim_source)
  int n_gem_hits = 0;      // GEM planes with a recorded hit from this particle
  bool hits_faec = false, hits_laec = false, hits_fa_muon = false, hits_la_muon = false;

  void Summarize(JObjectSummary& s) const override {
    s.add(index, NAME_OF(index), "%d");
    s.add(pdg, NAME_OF(pdg), "%d");
    s.add(status, NAME_OF(status), "%d");
    s.add(p4.p.mag(), "p", "%.4f", "GeV");
    s.add(p4.p.theta() * 180 / M_PI, "theta", "%.3f", "deg");
    s.add(p4.p.phi() * 180 / M_PI, "phi", "%.3f", "deg");
    s.add(vertex.z, "vz", "%.2f", "cm");
  }
};

/// A reconstructed-level GEM hit (2D point on a GEM plane). Both the fast simulation
/// and the EDM4hep digitizer produce these. mc_index < 0 marks a background hit,
/// mc_index >= 0 points to the MCParticle that produced it.
struct GEMHit : public JObject {
  JOBJECT_PUBLIC(GEMHit)
  int plane = 0;            // 1..6
  double x = 0, y = 0, z = 0;  // cm (z = nominal plane z)
  double sigma_x = 0.007;   // cm  (70 um)
  double sigma_y = 0.007;   // cm
  double time = 0;          // ns
  double edep = 0;          // GeV (0 for fast-sim)
  int mc_index = -1;        // truth link, -1 = background/ghost
  bool ghost = false;       // ghost combination from strip readout

  double r() const { return std::sqrt(x * x + y * y); }
  double phi() const { return std::atan2(y, x); }

  void Summarize(JObjectSummary& s) const override {
    s.add(plane, NAME_OF(plane), "%d");
    s.add(x, NAME_OF(x), "%.4f", "cm");
    s.add(y, NAME_OF(y), "%.4f", "cm");
    s.add(z, NAME_OF(z), "%.2f", "cm");
    s.add(mc_index, NAME_OF(mc_index), "%d");
    s.add(ghost, NAME_OF(ghost), "%d");
  }
};

/// Calorimeter cluster (FAEC or LAEC) used to seed the tracking and for e/mu identification.
struct ECCluster : public JObject {
  JOBJECT_PUBLIC(ECCluster)
  int calo = 0;             // 0 = FAEC, 1 = LAEC
  double x = 0, y = 0, z = 0;   // cm, cluster position at the calorimeter front face
  double sigma_xy = 1.0;    // cm
  double energy = 0;        // GeV
  int mc_index = -1;
  void Summarize(JObjectSummary& s) const override {
    s.add(calo, NAME_OF(calo), "%d");
    s.add(x, NAME_OF(x), "%.2f", "cm");
    s.add(y, NAME_OF(y), "%.2f", "cm");
    s.add(energy, NAME_OF(energy), "%.3f", "GeV");
    s.add(mc_index, NAME_OF(mc_index), "%d");
  }
};

/// Track parameters at the vertex: charge/momentum, direction, and vertex z on the beam line.
struct TrackParams {
  double qoverp = 0;   // e/GeV
  double theta = 0;    // rad
  double phi = 0;      // rad
  double zv = 0;       // cm (x=y=0 on the beam line)
  double p() const { return 1.0 / std::fabs(qoverp); }
  int charge() const { return qoverp > 0 ? +1 : -1; }
  Vec3 momentum() const {
    double pp = p();
    return {pp * std::sin(theta) * std::cos(phi), pp * std::sin(theta) * std::sin(phi), pp * std::cos(theta)};
  }
};

/// Fitted charged track.
struct Track : public JObject {
  JOBJECT_PUBLIC(Track)
  TrackParams params;
  std::array<double, 16> cov{};   // 4x4 covariance of (q/p, theta, phi, zv), row-major
  double chi2 = 0;
  int ndf = 0;
  int nhits = 0;
  std::vector<int> hit_indices;  // indices into the GEMHit collection used
  std::vector<int> planes;       // planes with a hit
  int mc_index = -1;             // majority truth match (-1 = fake)
  int n_bkg_hits = 0;            // hits on this track not belonging to mc_index
  double purity = 0;
  // Extrapolations used for PID/acceptance
  bool reaches_faec = false, reaches_laec = false, reaches_fa_muon = false, reaches_la_muon = false;
  double p_at_muon = 0;          // momentum after the FAEC/iron absorber (GeV)

  double p() const { return params.p(); }
  int charge() const { return params.charge(); }
  double sigma_p_over_p() const { return std::sqrt(cov[0]) * p(); }  // sigma(q/p) * p
  double sigma_theta() const { return std::sqrt(cov[5]); }
  double sigma_phi() const { return std::sqrt(cov[10]); }
  double sigma_zv() const { return std::sqrt(cov[15]); }

  void Summarize(JObjectSummary& s) const override {
    s.add(charge(), "q", "%d");
    s.add(p(), "p", "%.4f", "GeV");
    s.add(params.theta * 180 / M_PI, "theta", "%.3f", "deg");
    s.add(params.phi * 180 / M_PI, "phi", "%.3f", "deg");
    s.add(params.zv, "zv", "%.2f", "cm");
    s.add(chi2 / std::max(1, ndf), "chi2/ndf", "%.2f");
    s.add(nhits, NAME_OF(nhits), "%d");
    s.add(mc_index, NAME_OF(mc_index), "%d");
  }
};

/// DDVCS candidate: e' mu+ mu- with the recoil proton reconstructed as the missing particle.
struct DDVCSEvent : public JObject {
  JOBJECT_PUBLIC(DDVCSEvent)
  bool found = false;          // all three leptons reconstructed
  int electron_track = -1, muplus_track = -1, muminus_track = -1;  // indices into Track collection
  Vec4 electron, muplus, muminus;  // reconstructed four-vectors
  Vec4 missing;                // k + P - e' - mu+ - mu-
  double mm2 = 0;              // missing mass squared (GeV^2)
  double Q2 = 0, xB = 0, Qprime2 = 0, t = 0, W = 0;  // reconstructed kinematics
  double mm2_true = 0;         // from truth four-vectors (should equal Mp^2)
  double Q2_true = 0, xB_true = 0, Qp2_true = 0, t_true = 0;
  bool accepted_true = false;  // all three true leptons within acceptance with >= 4 GEM hits each
  bool geom_accepted = false;  // all three true leptons reach their PID detectors (no GEM-hit requirement)
  bool all_matched = false;    // the three reco tracks are matched to the right truth particles

  void Summarize(JObjectSummary& s) const override {
    s.add(found, NAME_OF(found), "%d");
    s.add(mm2, NAME_OF(mm2), "%.4f", "GeV^2");
    s.add(Q2, NAME_OF(Q2), "%.3f", "GeV^2");
    s.add(Qprime2, NAME_OF(Qprime2), "%.3f", "GeV^2");
    s.add(t, NAME_OF(t), "%.3f", "GeV^2");
  }
};

/// Generator-level kinematics of the DDVCS event.
struct GenKinematics : public JObject {
  JOBJECT_PUBLIC(GenKinematics)
  double Q2 = 0, xB = 0, nu = 0, W = 0, Qp2 = 0, t = 0, tmin = 0, phi_trento = 0, y = 0;
  double weight = 1.0;
  double luminosity = 0;   // cm^-2 s^-1 used for the background model
  void Summarize(JObjectSummary& s) const override {
    s.add(Q2, NAME_OF(Q2), "%.3f", "GeV^2");
    s.add(xB, NAME_OF(xB), "%.3f");
    s.add(Qp2, NAME_OF(Qp2), "%.3f", "GeV^2");
    s.add(t, NAME_OF(t), "%.3f", "GeV^2");
  }
};

constexpr double kProtonMass = 0.938272088;   // GeV
constexpr double kElectronMass = 0.000510999;  // GeV
constexpr double kMuonMass = 0.1056583755;     // GeV

}  // namespace solid
