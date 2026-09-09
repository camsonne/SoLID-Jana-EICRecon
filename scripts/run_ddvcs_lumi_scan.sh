#!/usr/bin/env bash
# Luminosity scan for the DDVCS missing-mass study.
#   run_ddvcs_lumi_scan.sh [OUTDIR] [NEVENTS] [NTHREADS]
# Requires `jana` in PATH and JANA_PLUGIN_PATH pointing at the directory containing solid_ddvcs.so.
set -euo pipefail
OUTDIR=${1:-results}
NEVENTS=${2:-20000}
NTHREADS=${3:-4}
LUMIS=${LUMIS:-"1.2e37 1e38 1e39"}
READOUTS=${READOUTS:-"strip pixel"}
EXTRA=${EXTRA:-"-Pjana:timeout=600 -Pjana:warmup_timeout=600"}
# Strip readout above ~5e37 is saturated (hit survival < 10 %, hundreds of ghost hits per plane):
# the tracking efficiency is ~0 and the events are slow, so fewer events are simulated there.
NEVENTS_SATURATED=${NEVENTS_SATURATED:-2000}
mkdir -p "$OUTDIR"
for ro in $READOUTS; do
  for L in $LUMIS; do
    label="L${L}_${ro}"
    echo "=== $label"
    n=$NEVENTS
    if [ "$ro" = "strip" ] && awk "BEGIN{exit !($L > 5e37)}"; then n=$NEVENTS_SATURATED; fi
    jana -Pplugins=solid_ddvcs -Pnthreads=$NTHREADS -Pfastsim:nevents=$n \
         -Psolid:luminosity=$L -Pgem:readout=$ro -Panalysis:label=$label \
         -Panalysis:output="$OUTDIR/${label}_events.csv" -Panalysis:summary="$OUTDIR/${label}_summary.json" \
         $EXTRA > "$OUTDIR/${label}.log" 2>&1
    grep "DDVCS summary" "$OUTDIR/${label}.log" | sed 's/.*DDVCS summary/DDVCS summary/'
  done
done
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python3 "$SCRIPT_DIR/plot_ddvcs_results.py" "$OUTDIR"
