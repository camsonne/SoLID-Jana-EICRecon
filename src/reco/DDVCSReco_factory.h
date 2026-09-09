// DDVCS event reconstruction: identify e', mu+, mu- among the fitted tracks and compute the
// missing mass of the undetected recoil proton,  MM^2 = (k + P - k' - l+ - l-)^2.
//
// Particle identification uses the detector-acceptance flags of the tracks combined with the
// truth association (the fast simulation does not model the Cherenkov / calorimeter / muon
// detector responses beyond acceptance): an electron candidate must reach a calorimeter, a muon
// candidate must reach the muon detector. Reconstructed four-vectors come only from the tracks.
#pragma once

#include <JANA/JFactory.h>
#include "datamodel/SoLIDDataModel.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

class DDVCSReco_factory : public JFactory {
public:
  DDVCSReco_factory();
  void Process(const JEvent& event) override;

private:
  Input<Track> m_tracks_in{this};
  Input<MCParticle> m_mc_in{this};
  Output<DDVCSEvent> m_out{this};
  Service<SoLIDGeometry_service> m_geo{this};
  Parameter<bool> m_truth_pid{this, "ddvcs:truth_pid", true,
                              "Use truth association for e/mu identification (acceptance flags still required)"};
};

}  // namespace solid
