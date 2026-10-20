// ACTS geometry and magnetic field for the SoLID DDVCS GEM tracker.
//
// Tracking geometry (ACTS Gen1, units mm): one cylindrical volume from upstream of the LH2 target
// to downstream of the last GEM plane, filled with air (homogeneous volume material), containing
// one sensitive PlaneSurface per GEM plane at the plane z: a ring (EllipseBounds rin..rout, widened
// by acts:surface_margin_mm) with local x/y = global x/y, carrying the plane material
// (gem:x0_fraction, path-length corrected by ACTS). The numbers come from
// SoLIDGeometry_service, the same source as the fast simulation and the native tracker, so the
// two trackers and the simulation describe the same detector. The LH2 target and its Al cell are
// not part of the ACTS geometry: the multiple scattering there happens at the vertex and is added
// to the perigee covariance after the fit (ActsTrackFinder_factory), together with the beam-line
// constraint.
//
// Field: SoLIDField_service (analytic CLEO-II coil model or the solenoid_CLEOv9.dat map) wrapped
// as an Acts::MagneticFieldProvider.
//
// When the DD4hep geometry (geometry/solid_ddvcs.xml) can be loaded, this service is the place
// to switch to Acts::DD4hepDetector (ACTS PluginDD4hep), as EICrecon does for ePIC.
#pragma once

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>
#include <JANA/JService.h>
#include <map>
#include <memory>
#include <vector>
#include "services/SoLIDField_service.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

/// SoLIDField_service seen by ACTS (positions in mm, field in ACTS units).
class SoLIDActsField final : public Acts::MagneticFieldProvider {
public:
  struct Cache {
    explicit Cache(const Acts::MagneticFieldContext& /*mctx*/) {}
  };
  explicit SoLIDActsField(const SoLIDField_service* field) : m_field(field) {}
  Acts::MagneticFieldProvider::Cache makeCache(const Acts::MagneticFieldContext& mctx) const override {
    return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }
  Acts::Result<Acts::Vector3> getField(const Acts::Vector3& position,
                                       Acts::MagneticFieldProvider::Cache& /*cache*/) const override;

private:
  const SoLIDField_service* m_field;
};

class SoLIDActs_service : public JService {
public:
  SoLIDActs_service();
  void Init() override;

  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const { return m_tgeo; }
  std::shared_ptr<const Acts::MagneticFieldProvider> field() const { return m_bfield; }
  /// Sensitive surface of GEM plane `plane` (1..6).
  const Acts::Surface& gemSurface(int plane) const { return *m_surfaces.at(plane - 1); }
  /// GEM plane (1..6) of a sensitive surface, 0 if unknown.
  int planeOf(const Acts::GeometryIdentifier& id) const {
    auto it = m_plane_of.find(id);
    return it == m_plane_of.end() ? 0 : it->second;
  }

  const Acts::GeometryContext& geoContext() const { return m_gctx; }
  const Acts::MagneticFieldContext& fieldContext() const { return m_mctx; }
  const Acts::CalibrationContext& calibContext() const { return m_cctx; }

private:
  Parameter<double> m_margin{this, "acts:surface_margin_mm", 30.0,
      "ACTS GEM surfaces extend this far beyond the active radii, so that tracks predicted just outside an "
      "edge are still followed there (hits only exist inside the active area) [mm]"};

  std::shared_ptr<const Acts::TrackingGeometry> m_tgeo;
  std::shared_ptr<const Acts::MagneticFieldProvider> m_bfield;
  std::vector<std::shared_ptr<const Acts::Surface>> m_surfaces;
  std::map<Acts::GeometryIdentifier, int> m_plane_of;
  Acts::GeometryContext m_gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  Acts::MagneticFieldContext m_mctx;
  Acts::CalibrationContext m_cctx;
};

}  // namespace solid
