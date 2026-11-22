#!/bin/bash
# Non-softening control (2026-11-13, O15): UW DruckerPrager, psi 0, perfectly plastic, on the Siddiquee pit-box geometry.
set -eu
cd $HOME/wp144/deck
export LD_LIBRARY_PATH=$HOME/wp144/lib
BIN=$HOME/wp144/bin/k5_0d38e57b4
S="--mat dp --sand toyoura --dr 0.8342 --dp-phi 38 --dp-g 30000 --target 0.20 --ckpt-every 50 --bproto 0.5 --gamma 15.60 --dom-w 14 --dom-d 8"
sub() {  # leg mesh time wall
  mkdir -p runs/$1/logs
  local a="$S --mesh $2 --wall $4"
  local jid=$(sbatch --parsable --time=$3 --partition=computes,new --exclude=node5 --mem=6G --job-name=ns_$1 \
    --export=ALL,LEG=$1,BIN=$BIN,ARGS="$a" leg.sbatch)
  echo "$(date -Is) $1 job $jid bin $BIN args: $a" | tee -a JOBS.txt
}
sub K5dp_sid_b4 b4 1:00:00 3000
sub K5dp_sid_b8 b8 1:00:00 3000
sub K5dp_sid_b16 b16 4:00:00 13000
