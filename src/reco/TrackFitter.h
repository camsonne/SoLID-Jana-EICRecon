// Global chi-square track fit with RK4 propagation in the SoLID field.
//
// Track model: (q/p, theta, phi, zv) at the beam line x=y=0 (SoLID vertex constraint).
// Measurements: (x,y) on each GEM plane. The covariance contains the intrinsic GEM resolution
// plus the fully correlated multiple-scattering contributions from the target, the air gaps and
// the GEM planes themselves (linearised kink model), so the chi-square is well defined and the
// parameter covariance is realistic.
#pragma once

#include <vector>
#include "TrackPropagator.h"
#include "datamodel/SoLIDDataModel.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

struct FitHit {
  int plane = 0;   // GEM plane 1..6, or 0 for a calorimeter point (no material, only position)
  double x = 0, y = 0, z = 0;
  double sx = 0.007, sy = 0.007;
  int index = -1;  // index in the input GEMHit collection
};

struct FitResult {
  bool ok = false;
  TrackParams params;
  std::array<double, 16> cov{};
  double chi2 = 0;
  int ndf = 0;
  std::vector<double> pulls;  // per measurement (x,y interleaved), uncorrelated approximation
  int iterations = 0;
};

class TrackFitter {
public:
  TrackFitter(const SoLIDGeometry_service* geo, const SoLIDField_service* field, double step_cm = 2.0);

  /// Predicted crossing points at the z of each measurement (GEM plane or calorimeter point).
  bool predict(const TrackParams& par, const std::vector<FitHit>& hits, std::vector<TrackState>& states) const;

  /// Seed from the beam line + two hits assuming an effective uniform field b_eff [T].
  bool seedFromTwoHits(const FitHit& a, const FitHit& b, double b_eff, TrackParams& out) const;

  /// Fit. zv_prior_sigma <= 0 disables the vertex-z prior.
  FitResult fit(const TrackParams& seed, const std::vector<FitHit>& hits, double zv_prior = 0,
                double zv_prior_sigma = 0, int max_iter = 10) const;

  /// Full 2N x 2N measurement covariance (intrinsic + multiple scattering) for a given track.
  void measurementCovariance(const TrackParams& par, const std::vector<FitHit>& hits,
                             const std::vector<TrackState>& states, std::vector<double>& V) const;

  const TrackPropagator& propagator() const { return m_prop; }

private:
  struct Scatterer { double z, path, theta0; Vec3 u; };
  void buildScatterers(const TrackParams& par, const std::vector<FitHit>& hits, const std::vector<TrackState>& states,
                       std::vector<Scatterer>& out) const;

  const SoLIDGeometry_service* m_geo;
  const SoLIDField_service* m_field;
  TrackPropagator m_prop;
};

/// In-place inversion of a dense n x n matrix (row-major). Returns false if singular.
bool invertMatrix(std::vector<double>& A, int n);

}  // namespace solid
