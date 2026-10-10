// Helpers shared by the native and the ACTS track finders: truth matching of a track from its
// GEM hits, and the extrapolation to the calorimeters / muon detector used for e/mu identification.
#pragma once

#include <vector>
#include "datamodel/SoLIDDataModel.h"
#include "reco/TrackPropagator.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

/// Fill mc_index (majority of the hits, -1 if no particle owns more than half), n_bkg_hits and
/// purity from trk.hit_indices.
void matchTrackToTruth(Track& trk, const std::vector<const GEMHit*>& hits);

/// Propagate the track from its vertex to the LAEC, FAEC and muon detectors and set the
/// reaches_* flags and p_at_muon.
void setAcceptanceFlags(Track& trk, const SoLIDGeometry_service& geo, const TrackPropagator& prop, double muon_pmin);

}  // namespace solid
