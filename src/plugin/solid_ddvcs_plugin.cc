// JANA2 plugin entry point: solid_ddvcs
//
//   jana -Pplugins=solid_ddvcs -Psolid:luminosity=1e38 -Pgem:readout=strip -Pfastsim:nevents=20000
//
// Registers the geometry and field services, the fast-simulation event source (unless an input
// file is given on the command line), the tracking and DDVCS factories and the analysis processor.
#include <JANA/JApplication.h>
#include <JANA/JFactoryGenerator.h>
#include "analysis/DDVCSAnalysis_processor.h"
#include "fastsim/SoLIDFastSim_source.h"
#include "reco/DDVCSReco_factory.h"
#include "reco/TrackFinder_factory.h"
#include "services/SoLIDField_service.h"
#include "services/SoLIDGeometry_service.h"
#ifdef SOLID_HAVE_EDM4HEP
#include "io/SoLIDEDM4hep_source.h"
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
  app->Add(new JFactoryGeneratorT<solid::TrackFinder_factory>());
  app->Add(new JFactoryGeneratorT<solid::DDVCSReco_factory>());
  app->Add(new solid::DDVCSAnalysis_processor());
}
}
