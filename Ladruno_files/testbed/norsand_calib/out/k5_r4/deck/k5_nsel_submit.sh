#!/bin/bash
# O15/L23 element control with NorSand K3 (2026-11-13): Siddiquee pit box (14 B by 8 B), rough, target s/B 0.20, as K5sp_sid_*_k3
# but with ELEM=qeas (LadrunoQuad EAS) or ELEM=lst (LadrunoLST, 2 per cell, node grid of the quads at the next mesh).
# LST caveat: the per-row prescribed K0 field is not in nodal equilibrium on the T6 (midside body load); one gravity step re-balances it.
set -eu
cd $HOME/wp144/deck
export LD_LIBRARY_PATH=$HOME/wp144/lib
BIN=$HOME/wp144/bin/k5_0d38e57b4
S="--mat norsand --sand toyoura --energy HAR --target 0.20 --ckpt-every 50 --bproto 0.5 --dr 0.8342 --gamma 15.60 --dom-w 14 --dom-d 8 --ns-params $HOME/wp144/deck/ns_params_toyoura_k3.json"
sub() {  # leg elem mesh time wall
  mkdir -p runs/$1/logs
  local a="$S --mesh $3 --wall $5"
  local jid=$(sbatch --parsable --time=$4 --partition=computes,new --exclude=node5 --mem=8G --job-name=ns_$1 \
    --export=ALL,ELEM=$2,LEG=$1,BIN=$BIN,ARGS="$a" leg_k5.sbatch)
  echo "$(date -Is) $1 job $jid env ELEM=$2 args: $a" | tee -a JOBS.txt
}
for e in qeas lst; do
  sub K5nse_${e}_b4 $e b4 3:00:00 10000
  sub K5nse_${e}_b8 $e b8 12:00:00 42000
  sub K5nse_${e}_b16 $e b16 2-00:00:00 170000
done
