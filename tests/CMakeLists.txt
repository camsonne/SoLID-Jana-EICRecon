# Smoke test: run 300 fast-sim events through the full chain and check the summary is produced.
add_test(NAME solid_ddvcs_smoke
  COMMAND ${CMAKE_COMMAND} -E env JANA_PLUGIN_PATH=${CMAKE_BINARY_DIR}/src
          jana -Pplugins=solid_ddvcs -Pfastsim:nevents=300 -Pnthreads=2
               -Panalysis:output=${CMAKE_BINARY_DIR}/smoke_events.csv
               -Panalysis:summary=${CMAKE_BINARY_DIR}/smoke_summary.json)
add_test(NAME solid_ddvcs_generator
  COMMAND $<TARGET_FILE:solid_ddvcs_gen> -n 200 -o ${CMAKE_BINARY_DIR}/smoke.hepmc)
