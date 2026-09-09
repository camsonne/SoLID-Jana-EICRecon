// SoLID DDVCS (J/psi-type) configuration geometry used by the fast simulation and the
// reconstruction. The numbers mirror geometry/solid_ddvcs.xml (DD4hep compact) and the
// solid_gemc "moved" (long endcap) geometry:
//   target        LH2, 15 cm long, centred at z = -315 cm
//   GEM planes    6 planes, DDVCS variant (Rout of planes 1-4 extended to 135 cm)
//   LAEC / FAEC   large- and forward-angle electromagnetic calorimeters
//   muon detector iron/scintillator layers behind the FAEC and around the coil
// Coordinates: z along the beam, origin at the solenoid centre, units cm.
#pragma once

#include <JANA/JService.h>
#include <JANA/Components/JComponent.h>
#include <string>
#include <vector>

namespace solid {

struct GEMPlane {
  int id;
  double z;      // cm
  double rin;    // cm
  double rout;   // cm
  double x0_fraction;   // material of the plane in units of X0
  double resolution;    // cm, single-coordinate spatial resolution
  double efficiency;    // intrinsic hit efficiency
  int nmodules;         // azimuthal modules (for the strip-readout ghost model)
};

struct Annulus {
  double z, rin, rout;
};

struct TargetGeometry {
  double z_center = -315.0;   // cm
  double half_length = 7.5;   // cm
  double radius = 1.8872;     // cm (inner radius of the cell)
  double wall_thickness = 0.0178;   // cm Al (1.905-1.8872)
  double window_thickness = 0.0051; // cm Al entrance/exit windows
  double x0_lh2 = 890.4;       // cm
  double x0_al = 8.897;        // cm
};

struct MuonDetector {
  // forward-angle: three iron absorbers (36 cm each) interleaved with scintillator annuli
  std::vector<Annulus> fa_scint;   // z of scint planes with radial acceptance
  double fa_iron_total = 108.0;    // cm Fe traversed to the last scintillator
  // large-angle: scintillator cylinder outside the endcap yoke
  double la_radius = 290.0, la_zmin = 209.0, la_zmax = 620.0;
  // barrel: scintillator outside the coil yoke
  double barrel_radius = 285.0, barrel_zmin = -266.0, barrel_zmax = 182.0;
};

class SoLIDGeometry_service : public JService {
public:
  SoLIDGeometry_service();
  void Init() override;

  const TargetGeometry& target() const { return m_target; }
  const std::vector<GEMPlane>& gemPlanes() const { return m_gem; }
  const GEMPlane& gemPlane(int id) const { return m_gem.at(id - 1); }
  const Annulus& laec() const { return m_laec; }
  const Annulus& faec() const { return m_faec; }
  const MuonDetector& muon() const { return m_muon; }
  double beamEnergy() const { return m_beam_energy_val; }

  // material between planes is air (X0 = 30390 cm) unless configured otherwise
  double airX0() const { return 30390.0; }
  double z_min() const { return m_target.z_center - m_target.half_length - 1.0; }
  double z_max() const { return m_muon.fa_scint.back().z + 1.0; }

  std::string describe() const;

private:
  Parameter<double> m_beam_energy{this, "solid:beam_energy", 11.0, "Beam energy [GeV]"};
  Parameter<double> m_target_z{this, "solid:target_z", -315.0, "Target centre z [cm]"};
  Parameter<double> m_target_length{this, "solid:target_length", 15.0, "Target length [cm]"};
  Parameter<double> m_gem_resolution{this, "gem:resolution_um", 70.0, "GEM single-coordinate resolution [um]"};
  Parameter<double> m_gem_x0{this, "gem:x0_fraction", 0.007, "Material per GEM plane [X0]"};
  Parameter<double> m_gem_efficiency{this, "gem:efficiency", 0.97, "Intrinsic GEM hit efficiency"};
  Parameter<std::string> m_variant{this, "solid:gem_variant", "ddvcs", "GEM radii variant: ddvcs (Rout 135 cm planes 1-4) or sidis"};

  double m_beam_energy_val = 11.0;
  TargetGeometry m_target;
  std::vector<GEMPlane> m_gem;
  Annulus m_laec{-65.0, 83.0, 137.0};   // preshower front face (LAEC z_shower = -35 cm)
  Annulus m_faec{425.0, 98.0, 230.0};   // preshower front face (FAEC z_shower = 433/467 cm)
  MuonDetector m_muon;
};

}  // namespace solid
