# Smoke test: run 300 fast-sim events through the full chain and check the summary is produced.
add_test(NAME solid_ddvcs_smoke
  COMMAND ${CMAKE_COMMAND} -E env JANA_PLUGIN_PATH=${CMAKE_BINARY_DIR}/src
          jana -Pplugins=solid_ddvcs -Pfastsim:nevents=300 -Pnthreads=2
               -Panalysis:output=${CMAKE_BINARY_DIR}/smoke_events.csv
               -Panalysis:summary=${CMAKE_BINARY_DIR}/smoke_summary.json)
add_test(NAME solid_ddvcs_generator
  COMMAND $<TARGET_FILE:solid_ddvcs_gen> -n 200 -o ${CMAKE_BINARY_DIR}/smoke.hepmc)

# EpIC HepMC3 input: frame conversion and kinematic conventions (phi, thetaL, phiL) against an
# EpIC-like file, then the full chain reading that file.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  add_executable(test_hepmc_reader test_hepmc_reader.cc)
  target_link_libraries(test_hepmc_reader PRIVATE solid_recon)
  add_test(NAME solid_epic_fixture
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/make_epic_like_hepmc.py
            ${CMAKE_BINARY_DIR}/epic_like.hepmc ${CMAKE_BINARY_DIR}/epic_like_truth.csv 300 7)
  set_tests_properties(solid_epic_fixture PROPERTIES FIXTURES_SETUP epic_file)
  add_test(NAME solid_hepmc_reader
    COMMAND test_hepmc_reader ${CMAKE_BINARY_DIR}/epic_like.hepmc ${CMAKE_BINARY_DIR}/epic_like_truth.csv)
  add_test(NAME solid_epic_chain
    COMMAND ${CMAKE_COMMAND} -E env JANA_PLUGIN_PATH=${CMAKE_BINARY_DIR}/src
            jana -Pplugins=solid_ddvcs -Pnthreads=2 -Pfastsim:hepmc_in=${CMAKE_BINARY_DIR}/epic_like.hepmc
                 -Panalysis:output=${CMAKE_BINARY_DIR}/epic_events.csv
                 -Panalysis:summary=${CMAKE_BINARY_DIR}/epic_summary.json)
  add_test(NAME solid_hepmc_to_solid
    COMMAND $<TARGET_FILE:solid_hepmc_to_solid> -i ${CMAKE_BINARY_DIR}/epic_like.hepmc
            -o ${CMAKE_BINARY_DIR}/epic_like_solid.hepmc)
  set_tests_properties(solid_hepmc_reader solid_epic_chain solid_hepmc_to_solid PROPERTIES FIXTURES_REQUIRED epic_file)
endif()

# Built-in tracker (the default engine when ACTS is not built)
add_test(NAME solid_ddvcs_smoke_native
  COMMAND ${CMAKE_COMMAND} -E env JANA_PLUGIN_PATH=${CMAKE_BINARY_DIR}/src
          jana -Pplugins=solid_ddvcs -Ptracking:engine=native -Pfastsim:nevents=300 -Pnthreads=2
               -Panalysis:output=${CMAKE_BINARY_DIR}/smoke_native_events.csv
               -Panalysis:summary=${CMAKE_BINARY_DIR}/smoke_native_summary.json)
