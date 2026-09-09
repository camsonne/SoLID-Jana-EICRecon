cmake_minimum_required(VERSION 3.16)
project(SoLIDJanaEICRecon VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# JANA2 (https://github.com/JeffersonLab/JANA2). Point JANA_ROOT or CMAKE_PREFIX_PATH at the install.
find_package(JANA REQUIRED)

option(USE_EDM4HEP "Build the EDM4hep/podio event source for reading ddsim output" OFF)
if(USE_EDM4HEP)
  find_package(podio REQUIRED)
  find_package(EDM4HEP REQUIRED)
endif()

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_INSTALL_PREFIX}/lib/SoLID/plugins;${JANA_DIR}/../../")
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

add_subdirectory(src)

# Install run-time configuration and scripts
install(DIRECTORY config/   DESTINATION share/SoLID/config)
install(DIRECTORY geometry/ DESTINATION share/SoLID/geometry)
install(PROGRAMS scripts/run_ddvcs_lumi_scan.sh scripts/plot_ddvcs_results.py scripts/run_ddsim_ddvcs.py
        DESTINATION bin)

enable_testing()
add_subdirectory(tests)
