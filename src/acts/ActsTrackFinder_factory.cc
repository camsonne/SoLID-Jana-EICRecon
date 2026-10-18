#include "ActsTrackFinder_factory.h"

#include <Acts/AmbiguityResolution/GreedyAmbiguityResolution.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/BoundTrackParameters.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/EventData/SeedContainer2.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/SpacePointContainer2.hpp>
#include <Acts/EventData/TrackContainer.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Propagator/SympyStepper.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Seeding2/BroadTripletSeedFilter.hpp>
#include <Acts/Seeding2/CylindricalSpacePointKDTree.hpp>
#include <Acts/Seeding2/DoubletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeeder.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>
#include <Acts/TrackFinding/TrackStateCreator.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/TrackHelpers.hpp>
#include <JANA/JEvent.h>
#include <JANA/JLogger.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <optional>
#include "reco/TrackUtils.h"

namespace solid {

using namespace Acts::UnitLiterals;

namespace {

using TrackContainer = Acts::TrackContainer<Acts::VectorTrackContainer, Acts::VectorMultiTrajectory, std::shared_ptr>;
using Propagator = Acts::Propagator<Acts::SympyStepper, Acts::Navigator>;
using CKF = Acts::CombinatorialKalmanFilter<Propagator, TrackContainer>;

/// A GEM hit as ACTS measurement: 2D (local x, y) on the plane surface.
struct Measurement {
  Acts::GeometryIdentifier gid;
  Acts::Vector2 loc;
  Acts::SquareMatrix2 cov;
  int hit = -1;   // index in the GEMHit collection
};

/// Source link: surface identifier + index into the measurement vector.
struct SoLIDSourceLink {
  Acts::GeometryIdentifier gid;
  std::uint32_t index = 0;
  friend bool operator==(const SoLIDSourceLink& a, const SoLIDSourceLink& b) {
    return a.gid == b.gid && a.index == b.index;
  }
};

struct GidLess {
  bool operator()(const SoLIDSourceLink& a, const Acts::GeometryIdentifier& b) const { return a.gid < b; }
  bool operator()(const Acts::GeometryIdentifier& a, const SoLIDSourceLink& b) const { return a < b.gid; }
};

struct SourceLinkAccessor {
  using BaseIterator = std::vector<SoLIDSourceLink>::const_iterator;
  using Iterator = Acts::SourceLinkAdapterIterator<BaseIterator>;
  const std::vector<SoLIDSourceLink>* links = nullptr;   // sorted by geometry identifier
  std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const {
    auto [b, e] = std::equal_range(links->begin(), links->end(), surface.geometryId(), GidLess{});
    return {Iterator{b}, Iterator{e}};
  }
};

struct Calibrator {
  const std::vector<Measurement>* meas = nullptr;
  void calibrate(const Acts::GeometryContext& /*gctx*/, const Acts::CalibrationContext& /*cctx*/,
                 const Acts::SourceLink& sourceLink, Acts::VectorMultiTrajectory::TrackStateProxy trackState) const {
    trackState.setUncalibratedSourceLink(Acts::SourceLink{sourceLink});
    const Measurement& m = meas->at(sourceLink.get<SoLIDSourceLink>().index);
    trackState.allocateCalibrated(m.loc, m.cov);
    trackState.setProjectorSubspaceIndices(std::array{Acts::eBoundLoc0, Acts::eBoundLoc1});
  }
};

struct SurfaceAccessor {
  const Acts::TrackingGeometry* geometry = nullptr;
  const Acts::Surface* operator()(const Acts::SourceLink& sl) const {
    return geometry->findSurface(sl.get<SoLIDSourceLink>().gid);
  }
};

/// Highland projected scattering angle [rad].
double highland(double p, double beta, double xx0) {
  if (xx0 <= 0) return 0;
  return 0.0136 / (beta * p) * std::sqrt(xx0) * (1 + 0.038 * std::log(xx0 / (beta * beta)));
}

}  // namespace

struct ActsTrackingImpl {
  std::unique_ptr<const Acts::Logger> logger;
  std::shared_ptr<const Acts::Surface> perigee;
  std::unique_ptr<Propagator> propagator;
  std::unique_ptr<CKF> ckf;
  Acts::GainMatrixUpdater updater;
  std::unique_ptr<Acts::MeasurementSelector> selector, preSelector;
  std::unique_ptr<Acts::GreedyAmbiguityResolution> ambiguity;
  Acts::TripletSeeder seeder;
  std::unique_ptr<Acts::DoubletSeedFinder> bottomFinder, topFinder;
  std::unique_ptr<Acts::TripletSeedFinder> tripletFinder;
  Acts::BroadTripletSeedFilter::Config filterConfig;
  Acts::Experimental::CylindricalSpacePointKDTree::Options lhOptions, hlOptions;
  Acts::Range1D<float> rMiddle{0.f, 0.f};
  Acts::Range1D<float> zMiddle{0.f, 0.f};

  explicit ActsTrackingImpl(std::unique_ptr<const Acts::Logger> l)
      : logger(std::move(l)), seeder(logger->cloneWithSuffix("Seeder")) {}
};

ActsTrackFinder_factory::ActsTrackFinder_factory() {
  SetTypeName(NAME_OF_THIS);
  m_hits_in.SetDatabundleName("GEMHits");
  m_ec_in.SetDatabundleName("ECClusters");
  m_ec_in.SetOptional(true);
  m_tracks_out.SetShortName("Tracks");
}

ActsTrackFinder_factory::~ActsTrackFinder_factory() = default;

void ActsTrackFinder_factory::Init() {
  const auto& acts = m_acts();
  const auto& geo = m_geo();
  m_prop = std::make_unique<TrackPropagator>(&m_field(), 2.0);
  const std::map<std::string, Acts::Logging::Level> levels = {
      {"VERBOSE", Acts::Logging::VERBOSE}, {"DEBUG", Acts::Logging::DEBUG}, {"INFO", Acts::Logging::INFO},
      {"WARNING", Acts::Logging::WARNING}, {"ERROR", Acts::Logging::ERROR}};
  const auto lvl = levels.count(m_log_level()) ? levels.at(m_log_level()) : Acts::Logging::WARNING;
  m_impl = std::make_unique<ActsTrackingImpl>(Acts::getDefaultLogger("SoLIDActs", lvl));
  auto& I = *m_impl;
  I.perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(0., 0., 0.));

  Acts::Navigator::Config ncfg{acts.trackingGeometry()};
  ncfg.resolvePassive = false;
  ncfg.resolveMaterial = true;
  ncfg.resolveSensitive = true;
  Acts::Navigator navigator(ncfg, I.logger->cloneWithSuffix("Navigator"));
  I.propagator = std::make_unique<Propagator>(Acts::SympyStepper(acts.field()), navigator,
                                              I.logger->cloneWithSuffix("Propagator"));
  I.ckf = std::make_unique<CKF>(*I.propagator, I.logger->cloneWithSuffix("CKF"));

  Acts::MeasurementSelectorCuts cuts;
  cuts.chi2CutOff = {m_chi2_cut()};
  cuts.numMeasurementsCutOff = {static_cast<std::size_t>(std::max(1, m_max_branches()))};
  I.selector = std::make_unique<Acts::MeasurementSelector>(cuts);
  Acts::MeasurementSelectorCuts preCuts;   // seed pre-fit: only the seed's own measurements are offered
  preCuts.chi2CutOff = {1e4};
  preCuts.numMeasurementsCutOff = {1};
  I.preSelector = std::make_unique<Acts::MeasurementSelector>(preCuts);

  Acts::GreedyAmbiguityResolution::Config gcfg;
  gcfg.maximumSharedHits = static_cast<std::uint32_t>(m_max_shared());
  gcfg.nMeasurementsMin = static_cast<std::size_t>(m_min_meas());
  I.ambiguity = std::make_unique<Acts::GreedyAmbiguityResolution>(gcfg, I.logger->cloneWithSuffix("Ambiguity"));

  // --- seeding configuration. The collision region is the LH2 target; the field used for the
  // seed helix is the solenoid field in the GEM region.
  const auto& tg = geo.target();
  const float zcol_min = static_cast<float>((tg.z_center - tg.half_length) * 10.0 - m_seed_z0_margin());
  const float zcol_max = static_cast<float>((tg.z_center + tg.half_length) * 10.0 + m_seed_z0_margin());
  const auto& planes = geo.gemPlanes();
  double rmax = 0;
  for (const auto& g : planes) rmax = std::max(rmax, g.rout * 10.0);
  const Vec3 bmid = m_field().B(Vec3(0.0, 0.0, 0.5 * (planes.front().z + planes.back().z)));
  const float bz = static_cast<float>(bmid.z * Acts::UnitConstants::T);

  Acts::DoubletSeedFinder::Config dcfg;
  dcfg.spacePointsSortedByRadius = false;
  dcfg.candidateDirection = Acts::Direction::Backward();
  dcfg.deltaRMin = static_cast<float>(m_seed_drmin());
  dcfg.deltaRMax = static_cast<float>(m_seed_drmax());
  dcfg.impactMax = static_cast<float>(m_seed_impact());
  dcfg.interactionPointCut = false;
  dcfg.collisionRegionMin = zcol_min;
  dcfg.collisionRegionMax = zcol_max;
  dcfg.cotThetaMax = static_cast<float>(m_seed_cotmax());
  dcfg.minPt = static_cast<float>(m_seed_minpt() * Acts::UnitConstants::GeV);
  I.bottomFinder = Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(dcfg, bz));
  dcfg.candidateDirection = Acts::Direction::Forward();
  I.topFinder = Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(dcfg, bz));

  Acts::TripletSeedFinder::Config tcfg;
  tcfg.useStripInfo = false;
  tcfg.sortedByCotTheta = true;
  tcfg.minPt = dcfg.minPt;
  tcfg.impactMax = dcfg.impactMax;
  tcfg.radLengthPerSeed = 0.03f;   // three GEM planes and the air between them
  I.tripletFinder = Acts::TripletSeedFinder::create(Acts::TripletSeedFinder::DerivedConfig(tcfg, bz));

  I.filterConfig.deltaRMin = dcfg.deltaRMin;
  I.filterConfig.maxSeedsPerSpM = static_cast<unsigned int>(m_seed_max_per_middle());
  I.filterConfig.impactWeightFactor = 0.1f;   // the impact parameter is only approximate here

  I.lhOptions.rMax = static_cast<float>(rmax + 10.0);
  I.lhOptions.zMin = zcol_min;
  I.lhOptions.zMax = static_cast<float>(planes.back().z * 10.0 + 50.0);
  I.lhOptions.deltaRMin = dcfg.deltaRMin;
  I.lhOptions.deltaRMax = dcfg.deltaRMax;
  I.lhOptions.collisionRegionMin = zcol_min;
  I.lhOptions.collisionRegionMax = zcol_max;
  I.lhOptions.cotThetaMax = dcfg.cotThetaMax;
  I.lhOptions.deltaPhiMax = static_cast<float>(m_seed_dphi());
  I.hlOptions = I.lhOptions;
  I.rMiddle = Acts::Range1D<float>(0.f, static_cast<float>(rmax + 10.0));
  I.zMiddle = Acts::Range1D<float>(static_cast<float>(planes.front().z * 10.0 - 50.0),
                                   static_cast<float>(planes.back().z * 10.0 + 50.0));

  LOG_INFO(GetLogger()) << "ACTS tracking: seeding Bz=" << bmid.z << " T, collision region z=[" << zcol_min << ","
                        << zcol_max << "] mm, CKF chi2<" << m_chi2_cut() << ", min " << m_min_meas()
                        << " measurements, beam-line constraint " << m_beamline_sigma() << " mm";
}

void ActsTrackFinder_factory::Process(const JEvent& event) {
  const auto& acts = m_acts();
  const auto& geo = m_geo();
  auto& I = *m_impl;
  const auto& gctx = acts.geoContext();
  const auto& mctx = acts.fieldContext();
  const auto& cctx = acts.calibContext();
  const auto hits = m_hits_in();
  std::vector<Track*> out;
  ++m_stats.events;

  // --- measurements, source links and space points
  std::vector<Measurement> meas;
  std::vector<SoLIDSourceLink> links;
  Acts::SpacePointContainer2 sps(Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
                                 Acts::SpacePointColumns::Phi | Acts::SpacePointColumns::VarianceZ |
                                 Acts::SpacePointColumns::VarianceR | Acts::SpacePointColumns::CopyFromIndex);
  Acts::Experimental::CylindricalSpacePointKDTreeBuilder kdBuilder;
  meas.reserve(hits.size());
  for (size_t i = 0; i < hits.size(); ++i) {
    const GEMHit& h = *hits[i];
    if (h.plane < 1 || h.plane > static_cast<int>(geo.gemPlanes().size())) continue;
    const Acts::Surface& surf = acts.gemSurface(h.plane);
    const Acts::Vector3 glob(h.x * 10.0, h.y * 10.0, surf.center(gctx).z());
    auto loc = surf.globalToLocal(gctx, glob, Acts::Vector3::UnitZ());
    if (!loc.ok()) continue;
    Measurement m;
    m.gid = surf.geometryId();
    m.loc = *loc;
    m.cov = Acts::SquareMatrix2::Zero();
    m.cov(0, 0) = std::pow(h.sigma_x * 10.0, 2);
    m.cov(1, 1) = std::pow(h.sigma_y * 10.0, 2);
    m.hit = static_cast<int>(i);
    links.push_back({m.gid, static_cast<std::uint32_t>(meas.size())});

    const Acts::SpacePointIndex2 spIndex = sps.size();
    auto sp = sps.createSpacePoint();
    sp.xy() = std::array<float, 2>{static_cast<float>(glob.x()), static_cast<float>(glob.y())};
    sp.zr() = std::array<float, 2>{static_cast<float>(glob.z()), static_cast<float>(std::hypot(glob.x(), glob.y()))};
    sp.phi() = static_cast<float>(std::atan2(glob.y(), glob.x()));
    sp.varianceZ() = 0.01f;   // (0.1 mm)^2: the plane z is known
    sp.varianceR() = static_cast<float>(m.cov(0, 0));
    sp.copyFromIndex() = static_cast<Acts::SpacePointIndex2>(meas.size());   // -> measurement index
    kdBuilder.insert(spIndex, sp.phi(), sp.zr()[1], sp.zr()[0]);
    meas.push_back(m);
  }
  m_stats.hits += static_cast<long>(meas.size());
  std::sort(links.begin(), links.end(), [](const SoLIDSourceLink& a, const SoLIDSourceLink& b) {
    return a.gid < b.gid || (a.gid == b.gid && a.index < b.index);
  });

  // --- triplet seeding
  Acts::SeedContainer2 seeds;
  if (meas.size() >= 3) {
    auto kdTree = kdBuilder.build();
    Acts::BroadTripletSeedFilter::State filterState;
    Acts::BroadTripletSeedFilter::Cache filterCache;
    auto filterLogger = I.logger->cloneWithSuffix("SeedFilter");
    Acts::BroadTripletSeedFilter seedFilter(I.filterConfig, filterState, filterCache, *filterLogger);
    Acts::TripletSeeder::Cache cache;
    Acts::Experimental::CylindricalSpacePointKDTree::Candidates candidates;
    seeds.assignSpacePointContainer(sps);
    for (const auto& middle : kdTree) {
      const auto spM = sps.at(middle.second).asConst();
      const float rM = spM.zr()[1], zM = spM.zr()[0];
      if (rM < I.rMiddle.min() || rM > I.rMiddle.max() || zM < I.zMiddle.min() || zM > I.zMiddle.max()) continue;
      candidates.clear();
      kdTree.validTuples(I.lhOptions, I.hlOptions, spM, 0, candidates);
      auto bottom = sps.subset(candidates.bottom_lh_v).asConst();
      auto top = sps.subset(candidates.top_lh_v).asConst();
      I.seeder.createSeedsFromGroup(cache, *I.bottomFinder, *I.topFinder, *I.tripletFinder, seedFilter, sps, bottom,
                                    spM, top, seeds);
      bottom = sps.subset(candidates.bottom_hl_v).asConst();
      top = sps.subset(candidates.top_hl_v).asConst();
      I.seeder.createSeedsFromGroup(cache, *I.bottomFinder, *I.topFinder, *I.tripletFinder, seedFilter, sps, bottom,
                                    spM, top, seeds);
    }
  }
  m_stats.seeds += static_cast<long>(seeds.size());

  // --- CKF set-up
  Calibrator calibrator{&meas};
  SourceLinkAccessor slAccessor{&links};
  using TrackStateCreator = Acts::TrackStateCreator<SourceLinkAccessor::Iterator, TrackContainer>;
  TrackStateCreator stateCreator;
  stateCreator.sourceLinkAccessor.template connect<&SourceLinkAccessor::range>(&slAccessor);
  stateCreator.calibrator.template connect<&Calibrator::calibrate>(&calibrator);
  stateCreator.measurementSelector.template connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
      I.selector.get());
  Acts::CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
  extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&I.updater);
  extensions.createTrackStates.template connect<&TrackStateCreator::createTrackStates>(&stateCreator);

  // seed pre-fit: the same CKF, offered only the three seed measurements
  std::vector<SoLIDSourceLink> seedLinks;
  SourceLinkAccessor seedAccessor{&seedLinks};
  TrackStateCreator preCreator;
  preCreator.sourceLinkAccessor.template connect<&SourceLinkAccessor::range>(&seedAccessor);
  preCreator.calibrator.template connect<&Calibrator::calibrate>(&calibrator);
  preCreator.measurementSelector.template connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
      I.preSelector.get());
  Acts::CombinatorialKalmanFilterExtensions<TrackContainer> preExtensions;
  preExtensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&I.updater);
  preExtensions.createTrackStates.template connect<&TrackStateCreator::createTrackStates>(&preCreator);

  Acts::PropagatorPlainOptions popts(gctx, mctx);
  popts.maxSteps = 10000;
  popts.direction = Acts::Direction::Forward();
  Acts::CombinatorialKalmanFilterOptions<TrackContainer> ckfOptions(gctx, mctx, cctx, extensions, popts);
  ckfOptions.energyLoss = m_eloss();
  ckfOptions.multipleScattering = true;
  Acts::CombinatorialKalmanFilterOptions<TrackContainer> preOptions(gctx, mctx, cctx, preExtensions, popts);
  preOptions.energyLoss = m_eloss();

  using ExtrapolatorOptions =
      Propagator::template Options<Acts::ActorList<Acts::MaterialInteractor, Acts::EndOfWorldReached>>;
  ExtrapolatorOptions extOptions(gctx, mctx);

  TrackContainer tracks(std::make_shared<Acts::VectorTrackContainer>(), std::make_shared<Acts::VectorMultiTrajectory>());
  TrackContainer tracksTemp(std::make_shared<Acts::VectorTrackContainer>(),
                            std::make_shared<Acts::VectorMultiTrajectory>());

  Acts::EstimateTrackParamCovarianceConfig covCfg;
  covCfg.initialSigmas = {2_mm, 2_mm, 1_degree, 1_degree, 0.1 / 1_GeV, 1_ns};
  covCfg.initialSigmaPtRel = 0.3;   // the seed helix uses a uniform field
  covCfg.initialVarInflation = {4., 4., 4., 4., 1., 1.};
  const auto hypothesis = Acts::ParticleHypothesis::muon();

  for (auto seed : seeds) {
    const auto idx = seed.spacePointIndices();
    if (idx.size() != 3) continue;
    std::array<Acts::Vector3, 3> pos;
    std::array<const Measurement*, 3> m3{};
    for (int k = 0; k < 3; ++k) {
      const auto sp = sps.at(idx[k]);
      pos[k] = Acts::Vector3(sp.xy()[0], sp.xy()[1], sp.zr()[0]);
      m3[k] = &meas[sp.copyFromIndex()];
    }
    const Acts::Surface* s0 = acts.trackingGeometry()->findSurface(m3[0]->gid);
    if (!s0) continue;
    // field at the middle space point for the seed helix
    Acts::Vector3 bfield = Acts::Vector3::Zero();
    {
      auto bcache = acts.field()->makeCache(mctx);
      auto b = acts.field()->getField(pos[1], bcache);
      if (b.ok()) bfield = *b;
    }
    auto bound = Acts::estimateTrackParamsFromSeed(gctx, *s0, pos[0], 0.0, pos[1], pos[2], bfield);
    if (!bound.ok()) { ++m_stats.failed_seeds; continue; }
    Acts::BoundVector params = *bound;
    params[Acts::eBoundTime] = 0;
    const Acts::BoundMatrix cov = Acts::estimateTrackParamCovariance(covCfg, params, false);
    Acts::BoundTrackParameters seedPars(s0->getSharedPtr(), params, cov, hypothesis);

    // transport the seed to the beam line, start the CKF there (as in EICrecon)
    Propagator::Options<> backOpts(gctx, mctx);
    backOpts.direction = Acts::Direction::Backward();
    backOpts.maxSteps = 10000;
    auto back = I.propagator->propagate(seedPars, *I.perigee, backOpts);
    if (!back.ok() || !back->endParameters.has_value()) { ++m_stats.failed_seeds; continue; }
    const Acts::BoundTrackParameters& seedStart = *back->endParameters;

    // pre-fit of the three seed measurements: perigee parameters with a realistic covariance, so
    // that the full CKF does not pick up a wrong hit on a plane where this track has none
    seedLinks.clear();
    for (const auto* m : m3) seedLinks.push_back({m->gid, static_cast<std::uint32_t>(m - meas.data())});
    std::sort(seedLinks.begin(), seedLinks.end(),
              [](const SoLIDSourceLink& a, const SoLIDSourceLink& b) { return a.gid < b.gid; });
    tracksTemp.clear();
    std::optional<Acts::BoundTrackParameters> refined;
    if (auto pre = I.ckf->findTracks(seedStart, preOptions, tracksTemp); pre.ok()) {
      for (auto& t : *pre) {
        if (t.nMeasurements() < 3) continue;
        auto c = tracksTemp.makeTrack();
        c.copyFrom(t);
        if (!Acts::smoothTrack(gctx, c, *I.logger).ok()) continue;
        if (!Acts::extrapolateTrackToReferenceSurface(c, *I.perigee, *I.propagator, extOptions,
                                                      Acts::TrackExtrapolationStrategy::first, *I.logger)
                 .ok())
          continue;
        refined = c.createParametersAtReference();
        break;
      }
    }
    if (!refined) { ++m_stats.failed_seeds; continue; }
    const Acts::BoundTrackParameters& start = *refined;

    tracksTemp.clear();
    auto res = I.ckf->findTracks(start, ckfOptions, tracksTemp);
    if (!res.ok()) { ++m_stats.failed_seeds; continue; }
    for (auto& t : *res) {
      if (t.nMeasurements() < static_cast<unsigned>(m_min_meas())) continue;
      auto cand = tracksTemp.makeTrack();
      cand.copyFrom(t);
      if (!Acts::smoothTrack(gctx, cand, *I.logger).ok()) continue;
      if (!Acts::extrapolateTrackToReferenceSurface(cand, *I.perigee, *I.propagator, extOptions,
                                                    Acts::TrackExtrapolationStrategy::first, *I.logger)
               .ok())
        continue;
      Acts::calculateTrackQuantities(cand);
      auto dest = tracks.makeTrack();
      dest.copyFrom(cand);
      ++m_stats.ckf_tracks;
    }
  }

  // --- ambiguity resolution on shared measurements
  Acts::GreedyAmbiguityResolution::State state;
  I.ambiguity->computeInitialState(
      tracks, state,
      [](const Acts::SourceLink& sl) {
        const auto& s = sl.get<SoLIDSourceLink>();
        return std::hash<std::uint64_t>{}(s.gid.value()) ^ (std::hash<std::uint32_t>{}(s.index) << 1);
      },
      [](const Acts::SourceLink& a, const Acts::SourceLink& b) {
        return a.get<SoLIDSourceLink>() == b.get<SoLIDSourceLink>();
      });
  I.ambiguity->resolve(state);

  // --- SoLID vertex treatment and conversion to Track
  const auto& tg = geo.target();
  for (auto iTrack : state.selectedTracks) {
    auto track = tracks.getTrack(state.trackTips.at(iTrack));
    if (!track.hasReferenceSurface()) continue;
    const int ndf = static_cast<int>(track.nDoF()) - 5;
    if (ndf <= 0 || track.chi2() / ndf > m_chi2ndf_max()) continue;
    Acts::BoundVector x = track.parameters();
    Acts::BoundMatrix C = track.covariance();

    // beam-line constraint: Kalman update with the pseudo-measurement d0 = 0
    if (m_beamline_sigma() > 0) {
      const double s2 = m_beamline_sigma() * m_beamline_sigma();
      const double denom = C(Acts::eBoundLoc0, Acts::eBoundLoc0) + s2;
      const Acts::BoundVector K = C.col(Acts::eBoundLoc0) / denom;
      x -= K * x[Acts::eBoundLoc0];
      C -= K * C.row(Acts::eBoundLoc0);
    }

    const double theta = x[Acts::eBoundTheta], phi = x[Acts::eBoundPhi], qop = x[Acts::eBoundQOverP];
    const double zv = x[Acts::eBoundLoc1] / 10.0;   // cm
    const double p = 1.0 / std::fabs(qop);

    // multiple scattering in the target between the vertex and the cell wall (side) or end window,
    // the same model as the fast simulation
    if (m_target_ms()) {
      const double sinth = std::max(1e-6, std::sin(theta)), costh = std::max(1e-6, std::cos(theta));
      const double L_side = tg.radius / sinth, L_end = (tg.z_center + tg.half_length - zv) / costh;
      const double xx0 = (L_side < L_end) ? L_side / tg.x0_lh2 + tg.wall_thickness / (tg.x0_al * sinth)
                                          : std::max(0.0, L_end) / tg.x0_lh2 + tg.window_thickness / (tg.x0_al * costh);
      const double m = hypothesis.mass();
      const double th0 = highland(p, p / std::sqrt(p * p + m * m), xx0);
      C(Acts::eBoundTheta, Acts::eBoundTheta) += th0 * th0;
      C(Acts::eBoundPhi, Acts::eBoundPhi) += th0 * th0 / (sinth * sinth);
    }

    auto* t = new Track();
    t->params.qoverp = qop;
    t->params.theta = theta;
    t->params.phi = phi;
    t->params.zv = zv;
    const int map[4] = {Acts::eBoundQOverP, Acts::eBoundTheta, Acts::eBoundPhi, Acts::eBoundLoc1};
    const double scale[4] = {1.0, 1.0, 1.0, 0.1};   // z0 mm -> cm
    for (int a = 0; a < 4; ++a)
      for (int b = 0; b < 4; ++b) t->cov[4 * a + b] = C(map[a], map[b]) * scale[a] * scale[b];
    t->chi2 = track.chi2();
    // ACTS nDoF is the sum of the measurement dimensions; subtract the five fitted parameters and
    // add the beam-line constraint, as the native tracker counts it
    t->ndf = static_cast<int>(track.nDoF()) - 5 + (m_beamline_sigma() > 0 ? 1 : 0);
    for (auto ts : track.trackStatesReversed()) {
      if (!ts.typeFlags().isMeasurement()) continue;
      const SoLIDSourceLink sl = ts.getUncalibratedSourceLink().get<SoLIDSourceLink>();
      const int hit = meas[sl.index].hit;
      t->hit_indices.push_back(hit);
      t->planes.push_back(hits[hit]->plane);
    }
    std::reverse(t->hit_indices.begin(), t->hit_indices.end());
    std::reverse(t->planes.begin(), t->planes.end());
    t->nhits = static_cast<int>(t->hit_indices.size());
    matchTrackToTruth(*t, hits);
    setAcceptanceFlags(*t, geo, *m_prop, m_muon_pmin());
    if (m_ec_match() > 0 && !matchesCluster(*t)) { delete t; continue; }
    out.push_back(t);
  }
  m_stats.tracks += static_cast<long>(out.size());
  m_tracks_out() = out;
}

bool ActsTrackFinder_factory::matchesCluster(const Track& trk) {
  // project the track to the depth of each cluster (cluster centroids are inside the shower)
  for (const auto* c : m_ec_in()) {
    TrackState s;
    s.pos = {0, 0, trk.params.zv};
    s.mom = trk.params.momentum();
    s.charge = trk.params.charge();
    if (!m_prop->propagateToZ(s, c->z)) continue;
    if (std::hypot(c->x - s.pos.x, c->y - s.pos.y) < m_ec_match()) return true;
  }
  return false;
}

void ActsTrackFinder_factory::Finish() {
  LOG_INFO(GetLogger()) << "ACTS tracking: " << m_stats.events << " events, " << m_stats.hits << " measurements, "
                        << m_stats.seeds << " seeds (" << m_stats.failed_seeds << " failed), " << m_stats.ckf_tracks
                        << " CKF tracks, " << m_stats.tracks << " tracks after ambiguity resolution and selection";
}

}  // namespace solid
