#include "SoLIDActs_service.h"
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/CylinderVolumeBounds.hpp>
#include <Acts/Geometry/PlaneLayer.hpp>
#include <Acts/Geometry/LayerArrayCreator.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Material/HomogeneousVolumeMaterial.hpp>
#include <Acts/Material/Material.hpp>
#include <Acts/Material/MaterialSlab.hpp>
#include <Acts/Surfaces/EllipseBounds.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/SurfaceArray.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <JANA/JApplication.h>
#include <sstream>

namespace solid {

using namespace Acts::UnitLiterals;

Acts::Result<Acts::Vector3> SoLIDActsField::getField(const Acts::Vector3& position,
                                                     Acts::MagneticFieldProvider::Cache& /*cache*/) const {
  const Vec3 b = m_field->B({position.x() / 10.0, position.y() / 10.0, position.z() / 10.0});   // mm -> cm, T
  return Acts::Result<Acts::Vector3>::success(Acts::Vector3(b.x, b.y, b.z) * Acts::UnitConstants::T);
}

SoLIDActs_service::SoLIDActs_service() { SetTypeName(NAME_OF_THIS); }

void SoLIDActs_service::Init() {
  auto app = GetApplication();
  auto geo = app->GetService<SoLIDGeometry_service>();
  auto field = app->GetService<SoLIDField_service>();
  m_bfield = std::make_shared<SoLIDActsField>(field.get());

  const auto& tg = geo->target();
  const auto& planes = geo->gemPlanes();
  const double zmin = (tg.z_center - tg.half_length) * 10.0 - 100.0;   // mm
  const double zmax = planes.back().z * 10.0 + 100.0;
  double rmax = 0;
  for (const auto& g : planes) rmax = std::max(rmax, g.rout * 10.0);
  rmax += m_margin() + 50.0;

  // GEM plane material: Kapton-like slab with thickness x0_fraction * X0 (only X0 matters here;
  // energy loss is off by default, consistent with the fast simulation).
  const Acts::Material gemMat = Acts::Material::fromMassDensity(28.6_cm, 70.9_cm, 12.7, 6.4, 1.42_g / 1_cm3);
  const Acts::Material air = Acts::Material::fromMassDensity(geo->airX0() * 1_cm, 747.6_m, 14.4, 7.3, 1.205e-3_g / 1_cm3);

  Acts::LayerVector layers;
  m_surfaces.clear();
  for (const auto& g : planes) {
    const Acts::Transform3 trafo(Acts::Translation3(0., 0., g.z * 10.0));
    // planar ring (local x, y = global x, y): Cartesian measurements, no phi wrap-around in the fit
    const double rin = std::max(0.0, g.rin * 10.0 - m_margin()), rout = g.rout * 10.0 + m_margin();
    auto bounds = std::make_shared<const Acts::EllipseBounds>(rin, rin, rout, rout);
    auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, bounds);
    surface->assignIsSensitive(true);   // free surface (no detector element): mark it as a measurement surface
    surface->assignSurfaceMaterial(std::make_shared<Acts::HomogeneousSurfaceMaterial>(
        Acts::MaterialSlab(gemMat, g.x0_fraction * gemMat.X0())));
    auto layer = Acts::PlaneLayer::create(trafo, bounds, std::make_unique<Acts::SurfaceArray>(surface), 5_mm);
    surface->associateLayer(*layer);
    layers.push_back(layer);
    m_surfaces.push_back(surface);
  }

  auto logger = Acts::getDefaultLogger("SoLIDActsGeometry", Acts::Logging::WARNING);
  Acts::LayerArrayCreator layArrCreator(Acts::LayerArrayCreator::Config{}, logger->clone());
  std::unique_ptr<const Acts::LayerArray> layArr(
      layArrCreator.layerArray(m_gctx, layers, zmin, zmax, Acts::BinningType::arbitrary, Acts::AxisDirection::AxisZ));

  const Acts::Transform3 trafoVol(Acts::Translation3(0., 0., 0.5 * (zmin + zmax)));
  auto boundsVol = std::make_shared<Acts::CylinderVolumeBounds>(0.0, rmax, 0.5 * (zmax - zmin));
  auto volume = std::make_shared<Acts::TrackingVolume>(trafoVol, boundsVol,
                                                       std::make_shared<Acts::HomogeneousVolumeMaterial>(air),
                                                       std::move(layArr), nullptr,
                                                       Acts::MutableTrackingVolumeVector{}, "SoLIDTracker");
  m_tgeo = std::make_shared<Acts::TrackingGeometry>(volume);

  m_plane_of.clear();
  std::ostringstream os;
  os << "ACTS tracking geometry: volume z=[" << zmin << "," << zmax << "] mm, R<" << rmax << " mm, air";
  for (size_t i = 0; i < m_surfaces.size(); ++i) {
    m_plane_of[m_surfaces[i]->geometryId()] = static_cast<int>(i) + 1;
    os << "\n  GEM" << i + 1 << " " << m_surfaces[i]->geometryId() << " z=" << planes[i].z * 10.0 << " mm";
  }
  LOG_INFO(GetLogger()) << os.str();
}

}  // namespace solid
