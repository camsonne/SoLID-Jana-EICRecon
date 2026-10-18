#!/bin/bash
# Build the DD4hep / Geant4 / EDM4hep stack for the SoLID DDVCS full simulation outside eic-shell
# (inside eic-shell everything except solid_dd4hep is already installed).
#
# conda-forge (micromamba) provides ROOT 6.36, Geant4 11.4 with data, HepMC3, compilers and tools;
# podio, EDM4hep and DD4hep (with Geant4, HepMC3 and EDM4hep output) are built from source, then the
# solid_dd4hep detector plugins. Tested on Linux x86_64 (4 cores: about 1 h, 8 GB of disk).
#
#   PREFIX=/opt/g4env SRC=/opt/src scripts/geant4/setup_stack.sh
#   source scripts/geant4/env.sh            # afterwards, in every shell
set -euo pipefail
PREFIX=${PREFIX:-/opt/g4env}
SRC=${SRC:-/opt/src}
MAMBA=${MAMBA:-micromamba}
mkdir -p "$SRC"
cd "$SRC"

# 1. conda-forge environment (python pinned to the version ROOT was built with: ROOT's CMake
#    config requires it exactly)
$MAMBA create -y -p "$PREFIX" -c conda-forge python=3.12.14 "root=6.36" "geant4=11.4.2" hepmc3 xerces-c \
  libboost-devel cmake make ninja cxx-compiler c-compiler fmt jinja2 pyyaml tbb-devel expat nlohmann_json eigen
# conda-forge's Geant4 exports -std=c++17 in Geant4_CXX_FLAGS, which overrides DD4hep's C++20
sed -i 's/ -std=c++17")/ -std=c++20")/' "$PREFIX/lib/cmake/Geant4/Geant4Config.cmake"

# 2. sources
[ -d podio ]    || git clone --depth 1 --branch v01-08 https://github.com/AIDASoft/podio
[ -d EDM4hep ]  || git clone --depth 1 --branch v01-02 https://github.com/key4hep/EDM4hep
[ -d DD4hep ]   || git clone --depth 1 --branch v01-38 https://github.com/AIDASoft/DD4hep
[ -d solid_dd4hep ] || git clone --depth 1 https://github.com/JeffersonLab/solid_dd4hep
if [ ! -f solid_gemc/field/solenoid_CLEOv9.dat ]; then   # field map (jlabsvn is not always reachable)
  mkdir -p solid_gemc && cd solid_gemc && git init -q && git remote add origin https://github.com/JeffersonLab/solid_gemc
  git config core.sparseCheckout true && echo "field/solenoid_CLEOv9.dat" > .git/info/sparse-checkout
  git fetch -q --depth 1 --filter=blob:none origin HEAD && git checkout -q FETCH_HEAD && cd ..
fi
# DD4hep refuses a Geant4 that reports C++17; its libraries work from C++20
sed -i 's/    message(FATAL_ERROR "Geant4 was compiled with C++${Geant4_CXX_STANDARD}, but DD4hep requires C++${CMAKE_CXX_STANDARD}")/    message(WARNING "Geant4 C++${Geant4_CXX_STANDARD}, DD4hep C++${CMAKE_CXX_STANDARD}")/' DD4hep/CMakeLists.txt
# solid_dd4hep: find_package(Acts) is not used by its sources; fmt::format needs <fmt/format.h>
sed -i 's/^find_package(Acts REQUIRED COMPONENTS Core PluginTGeo PluginDD4hep )/#&/' solid_dd4hep/CMakeLists.txt
for f in $(grep -ln "fmt::" solid_dd4hep/src/* | xargs grep -L "fmt/format.h"); do
  sed -i '0,/#include/s//#include <fmt\/format.h>\n#include/' "$f"
done

# 3. build (inside the environment)
cat > build.sh <<EOF
set -e
P=$PREFIX
cmake -S podio -B podio/build -G Ninja -DCMAKE_INSTALL_PREFIX=\$P -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=20 \\
  -DBUILD_TESTING=OFF -DENABLE_SIO=OFF -DENABLE_RNTUPLE=OFF
cmake --build podio/build -j4 && cmake --install podio/build
cmake -S EDM4hep -B EDM4hep/build -G Ninja -DCMAKE_INSTALL_PREFIX=\$P -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=20 \\
  -DBUILD_TESTING=OFF -DUSE_EXTERNAL_CATCH2=OFF
cmake --build EDM4hep/build -j4 && cmake --install EDM4hep/build
cmake -S DD4hep -B DD4hep/build -G Ninja -DCMAKE_INSTALL_PREFIX=\$P -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=20 \\
  -DDD4HEP_USE_GEANT4=ON -DDD4HEP_USE_XERCESC=ON -DDD4HEP_USE_EDM4HEP=ON -DDD4HEP_USE_HEPMC3=ON -DDD4HEP_USE_LCIO=OFF \\
  -DDD4HEP_IGNORE_GEANT4_TLS=ON -DBUILD_TESTING=OFF -DBUILD_DOCS=OFF -DCMAKE_PREFIX_PATH=\$P \\
  -DEXPAT_LIBRARY=\$P/lib/libexpat.so -DEXPAT_INCLUDE_DIR=\$P/include \\
  -DDD4HEP_BUILD_PACKAGES="DDRec DDDetectors DDCond DDAlign DDCAD DDG4 UtilityApps DDEve"
cmake --build DD4hep/build -j4 && cmake --install DD4hep/build
source \$P/bin/thisdd4hep.sh
cmake -S solid_dd4hep -B solid_dd4hep/build -G Ninja -DCMAKE_INSTALL_PREFIX=\$P -DCMAKE_BUILD_TYPE=Release \\
  -DCMAKE_CXX_STANDARD=20 -DCMAKE_PREFIX_PATH=\$P
cmake --build solid_dd4hep/build -j4 && cmake --install solid_dd4hep/build
EOF
$MAMBA run -p "$PREFIX" bash build.sh
echo "done: source scripts/geant4/env.sh (PREFIX=$PREFIX SRC=$SRC)"
