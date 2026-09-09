// Runge-Kutta (RK4) charged-particle propagation in the SoLID solenoid field.
#pragma once

#include <random>
#include "datamodel/SoLIDDataModel.h"
#include "services/SoLIDField_service.h"

namespace solid {

struct TrackState {
  Vec3 pos;   // cm
  Vec3 mom;   // GeV
  int charge = -1;
  double path = 0;  // accumulated path length [cm]
  Vec3 dir() const { return mom.unit(); }
};

/// Highland formula: RMS projected scattering angle [rad] for momentum p [GeV], beta, thickness x/X0.
double highlandTheta0(double p, double beta, double x_over_x0);

class TrackPropagator {
public:
  explicit TrackPropagator(const SoLIDField_service* field, double step_cm = 2.0) : m_field(field), m_step(step_cm) {}

  /// Propagate to the plane z = z_target (forward in z only). Returns false if the track turns
  /// around (pz <= 0) or leaves the fiducial radius before reaching the plane.
  /// If rng != nullptr, multiple scattering in air (radiation length air_x0 [cm]) is applied per step.
  bool propagateToZ(TrackState& s, double z_target, std::mt19937* rng = nullptr, double air_x0 = 0,
                    double rmax = 400.0) const;

  /// Propagate until the cylindrical radius reaches r_target (used for the barrel/large-angle muon
  /// detector). Stops at z_stop.
  bool propagateToR(TrackState& s, double r_target, double z_stop, double rmax = 400.0) const;

  /// Apply a random multiple-scattering kick with RMS projected angle theta0 to the momentum direction.
  static void scatter(TrackState& s, double theta0, std::mt19937& rng);

  /// Orthonormal basis (e1, e2) perpendicular to direction u.
  static void transverseBasis(const Vec3& u, Vec3& e1, Vec3& e2);

  double step() const { return m_step; }

private:
  void rk4Step(TrackState& s, double h) const;
  const SoLIDField_service* m_field;
  double m_step;
};

}  // namespace solid
