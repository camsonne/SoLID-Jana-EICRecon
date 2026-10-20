# Environment for the SoLID DDVCS full simulation (source after scripts/geant4/setup_stack.sh), e.g.
#   eval "$(micromamba shell hook -s bash)"; micromamba activate /opt/g4env
#   PREFIX=/opt/g4env SRC=/opt/src source scripts/geant4/env.sh
PREFIX=${PREFIX:-/opt/g4env}
SRC=${SRC:-/opt/src}
source "$PREFIX/bin/thisdd4hep.sh"
export SOLID_DD4HEP=$SRC/solid_dd4hep
export DETECTOR_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../geometry" && pwd)
export SOLID_FIELD_MAP=$SRC/solid_gemc/field/solenoid_CLEOv9.dat
# the solid_dd4hep field loader expects fieldmaps/<hash> in the working directory, with
# fieldmaps/solenoid_CLEOv9.dat pointing to it
mkdir -p fieldmaps
[ -e fieldmaps/c37b8743f2d1b205 ] || cp "$SOLID_FIELD_MAP" fieldmaps/c37b8743f2d1b205
[ -e fieldmaps/solenoid_CLEOv9.dat ] || ln -s c37b8743f2d1b205 fieldmaps/solenoid_CLEOv9.dat
