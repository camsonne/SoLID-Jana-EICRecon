// Fast simulation event source: DDVCS generator + charged-particle transport through the SoLID
// DDVCS geometry in the solenoid field (RK4), multiple scattering in the target, air and GEM
// planes, GEM hit smearing/efficiency, and luminosity-scaled background/ghost hits.
//
// Output collections:  MCParticle "MCParticles",  GEMHit "GEMHits",  ECCluster "ECClusters",  GenKinematics "GenKinematics"
#pragma once

#include <JANA/JEventSource.h>
#include <memory>
#include <random>
#include "fastsim/BackgroundModel.h"
#include "generator/DDVCSGenerator.h"
#include "generator/HepMC3Writer.h"
#include "reco/TrackPropagator.h"
#include "services/SoLIDField_service.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

class SoLIDFastSim_source : public JEventSource {
public:
  SoLIDFastSim_source();
  ~SoLIDFastSim_source() override = default;

  void Init() override;
  Result Emit(JEvent& event) override;
  static std::string GetDescription() { return "SoLID DDVCS fast simulation (generator + tracking-level detector response)"; }

private:
  void transport(MCParticle& mc, std::vector<GEMHit>& hits, std::vector<ECCluster>& clusters,
                 std::vector<std::pair<int, Vec3>>& crossings);

  Service<SoLIDGeometry_service> m_geo{this};
  Service<SoLIDField_service> m_field{this};

  Parameter<int> m_nevents{this, "fastsim:nevents", 10000, "Number of events to generate (0 = unlimited)"};
  Parameter<int> m_seed{this, "fastsim:seed", 20250908, "Random seed"};
  Parameter<double> m_luminosity{this, "solid:luminosity", 1.2e37, "Luminosity for the background model [cm^-2 s^-1]"};
  Parameter<std::string> m_readout{this, "gem:readout", "strip", "GEM readout model: strip | pixel"};
  Parameter<double> m_time_window{this, "gem:time_window_ns", 100.0, "GEM hit integration window [ns]"};
  Parameter<double> m_strip_pitch{this, "gem:strip_pitch_mm", 0.4, "Strip pitch [mm]"};
  Parameter<double> m_strip_length{this, "gem:strip_length_cm", 0.0, "Strip length [cm] (0 = module radial size, max 100)"};
  Parameter<double> m_pad_area{this, "gem:pad_area_mm2", 1.0, "Pad area for pixel readout [mm^2]"};
  Parameter<std::string> m_bkg_file{this, "gem:background_file", "", "Override file for the background rate parametrisation"};
  Parameter<double> m_bkg_road{this, "fastsim:bkg_road_cm", 3.0, "Half-width of the road in which background hits are generated [cm]"};
  Parameter<bool> m_ms{this, "fastsim:multiple_scattering", true, "Apply multiple scattering"};
  Parameter<double> m_ec_sigma{this, "ec:position_sigma_cm", 1.0, "Calorimeter cluster position resolution [cm]"};
  Parameter<double> m_ec_eres{this, "ec:energy_resolution", 0.05, "Calorimeter sampling term sigma_E/E = a/sqrt(E)"};
  Parameter<double> m_muon_pmin{this, "solid:muon_pmin", 2.0, "Minimum momentum for a muon to reach the last FA muon scintillator [GeV]"};
  Parameter<std::string> m_hepmc_out{this, "fastsim:hepmc_out", "", "If set, also write the generated events as HepMC3 ASCII"};
  // generator
  Parameter<double> m_Q2_min{this, "ddvcs:Q2_min", 1.0, "GeV^2"};
  Parameter<double> m_Q2_max{this, "ddvcs:Q2_max", 6.0, "GeV^2"};
  Parameter<double> m_xB_min{this, "ddvcs:xB_min", 0.10, ""};
  Parameter<double> m_xB_max{this, "ddvcs:xB_max", 0.55, ""};
  Parameter<double> m_Qp2_min{this, "ddvcs:Qp2_min", 2.0, "GeV^2 (dilepton mass squared)"};
  Parameter<double> m_Qp2_max{this, "ddvcs:Qp2_max", 9.0, "GeV^2"};
  Parameter<double> m_t_slope{this, "ddvcs:t_slope", 1.5, "GeV^-2"};
  Parameter<double> m_tabs_max{this, "ddvcs:tabs_max", 1.5, "GeV^2"};

  std::unique_ptr<DDVCSGenerator> m_gen;
  std::unique_ptr<BackgroundModel> m_bkg;
  std::unique_ptr<TrackPropagator> m_prop;
  std::unique_ptr<HepMC3Writer> m_hepmc;
  std::mt19937 m_rng;
  int m_emitted = 0;
};

}  // namespace solid
