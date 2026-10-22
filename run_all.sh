#!/bin/bash
# 2026-10-22: full re-run with the v11.3 tissue-mask / nuclear-threshold fix in 02_segment_measure.R.
# Populations as in the 2026-10-13 PDGFRA_v2 run (03/04 unchanged). All six slides segmented from scratch.
# Previous header:
#   PI16 = PDGFRA+ CEBPD+ ; EPI_CA = AR+ FOLH1+ panCK+ ; NE_CA = DLL3+ panCK+
# All six slides segmented from scratch in this run.
set -u
export PATH=/usr/local/bin:/opt/homebrew/bin:$PATH
cd "$(dirname "$0")"
OUT="/Volumes/T9/Case Studies/multiplex_analysis_20261022_maskfix"
LOG="$OUT/logs"
d() { echo "$(date) $*" >> "$LOG/driver.log"; }
d start
Rscript -e "source('00_config.R'); source('01_channels.R')" > "$LOG/00_01.log" 2>&1 || { d "00/01 failed"; exit 1; }
d "00/01 ok"
w() { Rscript -e "RUN_SLIDES<-c($1); source('02_segment_measure.R')" > "$LOG/02_worker_$2.log" 2>&1; d "worker $2 exit $?"; }
w '"PM0T32","PM423"' A &
w '"PM677_A8","PM12_B3"' B &
w '"PM159T5","PM1024_B"' C &
wait
for s in PM0T32 PM677_A8 PM12_B3 PM1024_B PM423 PM159T5; do
  [ -f "$OUT/cells_raw_$s.rds" ] || { d "missing $s, stopping"; exit 1; }
done
Rscript -e "REUSE_CHECKPOINTS<-TRUE; source('02_segment_measure.R')" > "$LOG/02_combine.log" 2>&1 || { d "02 combine failed"; exit 1; }
d "02 combine ok"
Rscript -e "source('03_classify.R')" > "$LOG/03.log" 2>&1 || { d "03 failed"; exit 1; }
d "03 ok"
Rscript -e "source('04_abundance_proximity.R')" > "$LOG/04.log" 2>&1 || { d "04 failed"; exit 1; }
d "ALL DONE"
