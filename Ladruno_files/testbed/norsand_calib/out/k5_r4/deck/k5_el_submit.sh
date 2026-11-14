#!/bin/bash
# O15 element-technology sensitivity (2026-11-13): associated Drucker-Prager phi 38, rough strip, Siddiquee pit box,
# target s/B 10 %, meshes b4 and b8; deck copy footing_ab_ns_k5.py, switch ELEM=<tech>.  Legs K5el_<tech>_<mesh>.
set -eu
cd $HOME/wp144/deck
export LD_LIBRARY_PATH=$HOME/wp144/lib
BIN=$HOME/wp144/bin/k5_0d38e57b4
S="--mat dp --sand toyoura --dr 0.8342 --dp-phi 38 --dp-g 30000 --target 0.10 --ckpt-every 50 --bproto 0.5 --gamma 15.60 --dom-w 14 --dom-d 8"
sub() {  # leg mesh time wall envs
  mkdir -p runs/$1/logs
  local a="$S --mesh $2 --wall $4"
  local jid=$(sbatch --parsable --time=$3 --partition=computes,new --exclude=node5 --mem=8G --job-name=el_$1 \
    --export=ALL,$5,LEG=$1,BIN=$BIN,ARGS="$a" leg_k5.sbatch)
  echo "$(date -Is) $1 job $jid env $5 args: $a" | tee -a JOBS.txt
}
for e in qbbar qstd qssp qeas cst lst bez bezb; do
  sub K5el_${e}_b4 b4 1:00:00 3300 DP_ASSOC=1,ELEM=$e
  sub K5el_${e}_b8 b8 2:00:00 6800 DP_ASSOC=1,ELEM=$e
done
