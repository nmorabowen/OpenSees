#!/bin/bash
# bash -lc "~/wp144/deck/k5_submit.sh <leg> <footing_ab_ns.py args...>"   (Slurm 24.11: needs a login shell for the PATH)
# Compute nodes lost libpython3.10.so.1.0 and libgfortran.so.5 in the 2026-11 reinstall: copies live in ~/wp144/lib.
set -u
LEG=$1; shift
D=$HOME/wp144/deck
mkdir -p $D/runs/$LEG/logs
export LD_LIBRARY_PATH=$HOME/wp144/lib
JID=$(sbatch --parsable --time=${TIME:-12:00:00} --partition=computes,new --exclude=node5 --mem=6G --job-name=ns_$LEG --export=ALL,LEG=$LEG,BIN=${BIN:-pdl},ARGS="$*" $D/leg.sbatch)
echo "$(date -Is) $LEG job $JID bin ${BIN:-pdl} args: $*" | tee -a $D/JOBS.txt
