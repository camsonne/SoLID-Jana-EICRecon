// JANA2 plugin entry point: solid_ddvcs
//
//   jana -Pplugins=solid_ddvcs -Psolid:luminosity=1e38 -Pgem:readout=strip -Pfastsim:nevents=20000
//   jana -Pplugins=solid_ddvcs -Ptracking:engine=native ...        (built-in tracker instead of ACTS)
//
// Registers the geometry and field services, the fast-simulation event source (unless an input
// file is given on the command line), the tracking and DDVCS factories and the analysis processor.
#include <JANA/JApplication.h>
#include <JANA/JFactoryGenerator.h>
#include "analysis/DDVCSAnalysis_processor.h"
#include "gnn/GNNGraphDump_processor.h"
#ifdef SOLID_HAVE_ONNX
#include "gnn/GNNTrackFinder_factory.h"
#endif
#include "fastsim/SoLIDFastSim_source.h"
#include "reco/DDVCSReco_factory.h"
#include "reco/TrackFinder_factory.h"
#include "services/SoLIDField_service.h"
#include "services/SoLIDGeometry_service.h"
#ifdef SOLID_HAVE_EDM4HEP
#include "io/SoLIDEDM4hep_source.h"
#endif
#ifdef SOLID_HAVE_ACTS
#include "acts/ActsTrackFinder_factory.h"
#include "acts/SoLIDActs_service.h"
#endif

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);
  app->ProvideService(std::make_shared<solid::SoLIDGeometry_service>());
  app->ProvideService(std::make_shared<solid::SoLIDField_service>());

  std::string source = "fastsim";
  app->SetDefaultParameter("solid:source", source, "Event source: fastsim | edm4hep (input files on the command line)");
  if (source == "fastsim") {
    app->Add(new solid::SoLIDFastSim_source());
  }
#ifdef SOLID_HAVE_EDM4HEP
  else if (source == "edm4hep") {
    app->Add(new JEventSourceGeneratorT<solid::SoLIDEDM4hep_source>());
  }
#endif
  // GNN training-data export: graphs only, no tracking or analysis
  std::string gnn_dump;
  app->SetDefaultParameter("gnn:dump", gnn_dump,
                           "If set, only export the GNN hit graphs to <prefix>_nodes.csv / <prefix>_edges.csv");
  if (!gnn_dump.empty()) {
    app->Add(new solid::GNNGraphDump_processor());
    return;
  }

  std::string engine =
#ifdef SOLID_HAVE_ACTS
      "acts";
#else
      "native";
#endif
  app->SetDefaultParameter("tracking:engine", engine,
                           "Track finding/fitting: acts (ACTS seeding + CKF, needs -DUSE_ACTS=ON) | "
                           "gnn (graph neural network edge classifier, needs -DUSE_ONNX=ON) | native");
  if (engine == "acts") {
#ifdef SOLID_HAVE_ACTS
    app->ProvideService(std::make_shared<solid::SoLIDActs_service>());
    app->Add(new JFactoryGeneratorT<solid::ActsTrackFinder_factory>());
#else
    throw JException("tracking:engine=acts requested but this build has no ACTS (configure with -DUSE_ACTS=ON)");
#endif
  } else if (engine == "gnn") {
#ifdef SOLID_HAVE_ONNX
    app->Add(new JFactoryGeneratorT<solid::GNNTrackFinder_factory>());
#else
    throw JException("tracking:engine=gnn requested but this build has no ONNX Runtime (configure with -DUSE_ONNX=ON)");
#endif
  } else if (engine == "native") {
    app->Add(new JFactoryGeneratorT<solid::TrackFinder_factory>());
  } else {
    throw JException("tracking:engine must be acts, gnn or native, got '%s'", engine.c_str());
  }
  app->Add(new JFactoryGeneratorT<solid::DDVCSReco_factory>());
  app->Add(new solid::DDVCSAnalysis_processor());
}
}
