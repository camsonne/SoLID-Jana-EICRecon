// Magnetic field of the SoLID (CLEO-II) solenoid.
//
// Two back-ends:
//  * "map": a (r,z) field map in the solid_gemc ASCII format used by solid_dd4hep's
//    SoLID_FieldMapB (columns r[cm] z[cm] Br[G] Bz[G], any header lines starting with '<'
//    are skipped), e.g. solenoid_CLEOv9.dat.  Set field:map=/path/to/solenoid_CLEOv9.dat.
//  * "analytic" (default): Biot-Savart field of the CLEO-II coil packs taken from
//    solid_gemc/field/CLEOv9.am (two layers at r=152.5/154.1 cm, z in [-173.8,173.8] cm,
//    5902 A/cm centre and 6140 A/cm ends) computed with complete elliptic integrals and
//    tabulated on an (r,z) grid.  The iron yoke is not modelled; the current is scaled so
//    that Bz(0,0) matches field:central_tesla (default 1.5 T).
#pragma once

#include <JANA/JService.h>
#include <JANA/Components/JComponent.h>
#include <string>
#include <vector>
#include "datamodel/SoLIDDataModel.h"

namespace solid {

class SoLIDField_service : public JService {
public:
  SoLIDField_service();
  void Init() override;

  /// Field in Tesla at position (cm).
  Vec3 B(const Vec3& pos) const;
  /// Cylindrical components (Br, Bz) in Tesla at (r, z) in cm.
  void BrBz(double r, double z, double& Br, double& Bz) const;

  bool isUniform() const { return m_uniform; }

private:
  void buildAnalyticGrid();
  bool loadMap(const std::string& file);
  static void ellipticKE(double k2, double& K, double& E);
  static void loopField(double a, double z0, double I, double r, double z, double& Br, double& Bz);

  Parameter<std::string> m_mode{this, "field:mode", "analytic", "analytic | map | uniform"};
  Parameter<std::string> m_map_file{this, "field:map", "", "GEMC-format r z Br Bz map (cm, Gauss)"};
  Parameter<double> m_central{this, "field:central_tesla", 1.5, "Bz at the magnet centre for the analytic model [T]"};
  Parameter<double> m_uniform_bz{this, "field:uniform_bz", 1.5, "Bz for field:mode=uniform [T]"};
  Parameter<double> m_scale{this, "field:scale", 1.0, "Overall multiplicative scale applied to the field"};

  bool m_uniform = false;
  double m_uniform_val = 1.5, m_scale_val = 1.0;
  // grid
  double m_rmin = 0, m_rmax = 300, m_dr = 1.0;
  double m_zmin = -450, m_zmax = 850, m_dz = 2.0;
  int m_nr = 0, m_nz = 0;
  std::vector<float> m_Br, m_Bz;  // [ir*m_nz + iz]
  bool m_have_grid = false;
};

}  // namespace solid
