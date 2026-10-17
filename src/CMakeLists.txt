# Core reconstruction library: datamodel, services, generator, fast simulation, tracking, DDVCS reco.
set(SOLID_SOURCES
  services/SoLIDGeometry_service.cc
  services/SoLIDField_service.cc
  generator/DDVCSGenerator.cc
  generator/HepMC3Writer.cc
  generator/HepMC3Reader.cc
  fastsim/SoLIDFastSim_source.cc
  fastsim/BackgroundModel.cc
  reco/TrackPropagator.cc
  reco/TrackFitter.cc
  reco/TrackFinder_factory.cc
  reco/TrackUtils.cc
  gnn/GNNGraph.cc
  gnn/GNNGraphDump_processor.cc
  reco/DDVCSReco_factory.cc
  analysis/DDVCSAnalysis_processor.cc
)

add_library(solid_recon SHARED ${SOLID_SOURCES})
target_include_directories(solid_recon PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/SoLID>)
target_link_libraries(solid_recon PUBLIC JANA::jana2_shared_lib)
target_compile_options(solid_recon PRIVATE -Wall -Wextra -Wno-unused-parameter)

if(USE_EDM4HEP)
  target_sources(solid_recon PRIVATE io/SoLIDEDM4hep_source.cc)
  target_link_libraries(solid_recon PUBLIC podio::podio podio::podioRootIO EDM4HEP::edm4hep)
  target_compile_definitions(solid_recon PUBLIC SOLID_HAVE_EDM4HEP=1)
endif()

if(USE_ACTS)
  target_sources(solid_recon PRIVATE acts/SoLIDActs_service.cc acts/ActsTrackFinder_factory.cc)
  target_link_libraries(solid_recon PUBLIC Acts::Core)
  target_compile_definitions(solid_recon PUBLIC SOLID_HAVE_ACTS=1)
endif()

if(USE_ONNX)
  target_sources(solid_recon PRIVATE gnn/GNNEdgeClassifier.cc gnn/GNNTrackFinder_factory.cc)
  target_link_libraries(solid_recon PUBLIC onnxruntime::onnxruntime)
  target_compile_definitions(solid_recon PUBLIC SOLID_HAVE_ONNX=1
    SOLID_GNN_DEFAULT_MODEL="${PROJECT_SOURCE_DIR}/config/gnn/solid_gem_in.onnx")
endif()

# JANA plugin: solid_ddvcs.so
add_library(solid_ddvcs SHARED plugin/solid_ddvcs_plugin.cc)
set_target_properties(solid_ddvcs PROPERTIES PREFIX "" SUFFIX ".so")
target_link_libraries(solid_ddvcs PUBLIC solid_recon)

# Stand-alone generator executable (writes HepMC3 ASCII for ddsim / DDG4)
add_executable(solid_ddvcs_gen generator/ddvcs_gen_main.cc)
target_link_libraries(solid_ddvcs_gen PUBLIC solid_recon)

# Converter for external generators (EpIC HepMC3 -> SoLID frame and target vertex, for ddsim / DDG4)
add_executable(solid_hepmc_to_solid generator/hepmc_to_solid_main.cc)
target_link_libraries(solid_hepmc_to_solid PUBLIC solid_recon)

install(TARGETS solid_recon solid_ddvcs_gen solid_hepmc_to_solid
        LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
install(TARGETS solid_ddvcs LIBRARY DESTINATION lib/SoLID/plugins)
install(DIRECTORY ./ DESTINATION include/SoLID FILES_MATCHING PATTERN "*.h")
