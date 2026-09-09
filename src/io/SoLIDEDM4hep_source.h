// Event source reading DDG4 / ddsim output (EDM4hep in podio ROOT files) and converting it to
// the SoLID lightweight data model:
//   edm4hep::MCParticle          "MCParticles"      -> solid::MCParticle  "MCParticles"
//   edm4hep::SimTrackerHit       "GEMTrackerHits"   -> solid::GEMHit      "GEMHits"   (digitised)
//   edm4hep::SimCalorimeterHit   FAEC_*/LAEC_*Hits  -> solid::ECCluster   "ECClusters"
//   edm4hep::SimTrackerHit       "MuonHits"         -> acceptance flags on the MCParticles
// Requires podio and EDM4hep (>= 0.99 for SimTrackerHit::getParticle). Built with -DUSE_EDM4HEP=ON.
#pragma once

#include <JANA/JEventSource.h>
#include <memory>
#include <random>
#include <string>
#include "fastsim/BackgroundModel.h"
#include "services/SoLIDGeometry_service.h"

namespace podio { class ROOTReader; }

namespace solid {

class SoLIDEDM4hep_source : public JEventSource {
public:
  SoLIDEDM4hep_source(std::string resource_name, JApplication* app);
  ~SoLIDEDM4hep_source() override;

  void Init() override;
  void Open() override;
  void Close() override;
  Result Emit(JEvent& event) override;
  static std::string GetDescription() { return "EDM4hep (ddsim) reader for the SoLID DDVCS geometry"; }

private:
  Service<SoLIDGeometry_service> m_geo{this};
  Parameter<double> m_luminosity{this, "solid:luminosity", 1.2e37, "Luminosity for the background model [cm^-2 s^-1]"};
  Parameter<std::string> m_readout{this, "gem:readout", "strip", "GEM readout model: strip | pixel"};
  Parameter<double> m_time_window{this, "gem:time_window_ns", 100.0, "GEM hit integration window [ns]"};
  Parameter<std::string> m_bkg_file{this, "gem:background_file", "", "Override file for the background rate parametrisation"};
  Parameter<double> m_bkg_road{this, "fastsim:bkg_road_cm", 3.0, "Half-width of the road in which background hits are generated [cm]"};
  Parameter<double> m_ec_sigma{this, "ec:position_sigma_cm", 1.0, "Calorimeter cluster position resolution [cm]"};
  Parameter<double> m_ec_cluster_radius{this, "ec:cluster_radius_cm", 15.0, "Clustering radius for calorimeter hits [cm]"};
  Parameter<double> m_edep_min{this, "gem:edep_min_keV", 0.5, "Minimum energy deposit for a GEM hit [keV]"};
  Parameter<int> m_seed{this, "fastsim:seed", 20250908, "Random seed for digitisation and background"};

  std::unique_ptr<podio::ROOTReader> m_reader;
  std::unique_ptr<BackgroundModel> m_bkg;
  std::mt19937 m_rng;
  unsigned m_entries = 0, m_current = 0;
};

}  // namespace solid
