# Core reconstruction library: datamodel, services, generator, fast simulation, tracking, DDVCS reco.
set(SOLID_SOURCES
  services/SoLIDGeometry_service.cc
  services/SoLIDField_service.cc
  generator/DDVCSGenerator.cc
  generator/HepMC3Writer.cc
  fastsim/SoLIDFastSim_source.cc
  fastsim/BackgroundModel.cc
  reco/TrackPropagator.cc
  reco/TrackFitter.cc
  reco/TrackFinder_factory.cc
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

# JANA plugin: solid_ddvcs.so
add_library(solid_ddvcs SHARED plugin/solid_ddvcs_plugin.cc)
set_target_properties(solid_ddvcs PROPERTIES PREFIX "" SUFFIX ".so")
target_link_libraries(solid_ddvcs PUBLIC solid_recon)

# Stand-alone generator executable (writes HepMC3 ASCII for ddsim / DDG4)
add_executable(solid_ddvcs_gen generator/ddvcs_gen_main.cc)
target_link_libraries(solid_ddvcs_gen PUBLIC solid_recon)

install(TARGETS solid_recon solid_ddvcs_gen
        LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
install(TARGETS solid_ddvcs LIBRARY DESTINATION lib/SoLID/plugins)
install(DIRECTORY ./ DESTINATION include/SoLID FILES_MATCHING PATTERN "*.h")
