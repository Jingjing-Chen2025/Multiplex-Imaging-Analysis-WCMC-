#!/bin/bash
# Segment all slides in 3 parallel R workers (2 slides each), then combine.
# Per-slide checkpoints (cells_raw_<sid>.rds) make this safe to re-run.
set -u
cd "$(dirname "$0")"
OUT=${OUT:-/Volumes/T9/results_ASPC_NEPC}
LOG="$OUT/logs"; mkdir -p "$LOG"
run_worker() {
  Rscript -e "OUT<-'$OUT'; RUN_SLIDES<-c($1); COMBINE<-FALSE; source('02_segment_measure.R')" > "$LOG/02_worker_$2.log" 2>&1
}
run_worker '"PM677_A8","PM423"'  1 &
run_worker '"PM0T32","PM12_B3"'  2 &
run_worker '"PM159T5","PM1024_B"' 3 &
wait
Rscript -e "OUT<-'$OUT'; source('02_segment_measure.R')" > "$LOG/02_combine.log" 2>&1
echo "02 done, exit $?"
