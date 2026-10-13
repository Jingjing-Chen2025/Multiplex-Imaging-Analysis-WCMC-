#!/bin/bash
# Segment all slides in 3 parallel R workers (2 slides each), then combine.
# Per-slide checkpoints (cells_raw_<sid>.rds) are reused on the combine pass.
set -u
cd "$(dirname "$0")"
OUT=${OUT:-"/Volumes/T9/Case Studies/multiplex_analysis_20261013_PDGFRA_v2"}
LOG="$OUT/logs"; mkdir -p "$LOG"
run_worker() {
  Rscript -e "OUT<-'$OUT'; RUN_SLIDES<-c($1); source('02_segment_measure.R')" > "$LOG/02_worker_$2.log" 2>&1
}
run_worker '"PM0T32","PM423"'      A &
run_worker '"PM677_A8","PM12_B3"'  B &
run_worker '"PM159T5","PM1024_B"'  C &
wait
Rscript -e "OUT<-'$OUT'; REUSE_CHECKPOINTS<-TRUE; source('02_segment_measure.R')" > "$LOG/02_combine.log" 2>&1
echo "02 done, exit $?"
