// Luminosity-scaled GEM background model.
//
// The background hit rate density on GEM plane i at radius R is parametrised as
//     rho_i(R) = rho0_i * (R0 / R)^alpha  * (L / L_ref)          [Hz/mm^2]
// with rho0_i the rate at R0 = 40 cm for the reference luminosity L_ref = 1.2e37 cm^-2 s^-1
// (SoLID J/psi configuration: 15 cm LH2, 3 uA). The default rho0_i values are of the order of
// the SoLID pCDR background estimates (~1 kHz/mm^2 at the inner edge of the first planes) and
// carry a factor ~2-3 uncertainty; override them with gem:background_file (one line per plane:
// "plane rho0_kHz_per_mm2 alpha") when results from a dedicated solid_gemc background run are
// available. Hits within a GEM time window T are integrated: n/mm^2 = rho * T.
//
// Two readout models:
//   strip : 2D u/v strips of pitch p and length L. Strip occupancy o = rho*T*p*L; the true hit
//           survives with probability exp(-2 o) (no overlapping cluster in either coordinate) and
//           uncorrelated u/v combinations produce ghost hits of density (rho*T)^2 * A_module.
//   pixel : pads of area A_pad; occupancy o = rho*T*A_pad, no ghosts.
#pragma once

#include <random>
#include <string>
#include <vector>
#include "datamodel/SoLIDDataModel.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

struct BackgroundConfig {
  double luminosity = 1.2e37;       // cm^-2 s^-1
  double ref_luminosity = 1.2e37;   // cm^-2 s^-1
  double time_window_ns = 100.0;
  std::string readout = "strip";    // strip | pixel
  double strip_pitch_mm = 0.4;
  double strip_length_cm = 0;       // 0 -> (rout - rin) capped at 100 cm
  double pad_area_mm2 = 1.0;
  double R0_cm = 40.0;
  std::vector<double> rho0_kHz_mm2 = {1.0, 1.0, 0.8, 0.6, 0.5, 0.4};  // per plane at R0, L_ref
  std::vector<double> alpha = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
  std::string file;                 // optional override file
};

class BackgroundModel {
public:
  BackgroundModel(const BackgroundConfig& cfg, const SoLIDGeometry_service* geo);

  /// Hit rate density [Hz/mm^2] at the configured luminosity.
  double rateDensity(int plane, double R) const;
  /// Integrated hit density [1/cm^2] within the time window.
  double hitDensity(int plane, double R) const;
  /// Channel occupancy (dimensionless).
  double occupancy(int plane, double R) const;
  /// Survival probability of a true hit against cluster overlap.
  double survivalProbability(int plane, double R) const;
  /// Ghost hit density [1/cm^2] (strip readout only).
  double ghostDensity(int plane, double R) const;

  /// Generate background (and ghost) hits in a square road of half-width w [cm] around (x,y).
  void generateInRoad(int plane, double x, double y, double w, std::mt19937& rng, std::vector<GEMHit>& out) const;

  const BackgroundConfig& config() const { return m_cfg; }
  std::string summary() const;

private:
  double moduleArea(int plane) const;   // cm^2
  double stripLength(int plane) const;  // cm
  BackgroundConfig m_cfg;
  const SoLIDGeometry_service* m_geo;
};

}  // namespace solid
