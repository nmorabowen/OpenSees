#!/bin/bash
# Siddiquee et al. (1999) footing with the test pit as the box (2026-11-12): pit 7 m long, about 4 m deep, B 0.5 m
# -> total width 14 B, depth 8 B (2001 paper p. 61, as reported by the SANISAND TIM session). Run with bash -lc on Esmeralda.
set -eu
cd $HOME/wp144/deck
export LD_LIBRARY_PATH=$HOME/wp144/lib
BIN=$HOME/wp144/bin/k5_0d38e57b4
S="--mat norsand --sand toyoura --energy HAR --target 0.20 --ckpt-every 50 --bproto 0.5 --dr 0.8342 --gamma 15.60 --dom-w 14 --dom-d 8"
sub() {  # leg params mesh time wall
  mkdir -p runs/$1/logs
  local a="$S --ns-params $HOME/wp144/deck/$2 --mesh $3 --wall $5"
  local jid=$(sbatch --parsable --time=$4 --partition=computes,new --exclude=node5 --mem=6G --job-name=ns_$1 \
    --export=ALL,LEG=$1,BIN=$BIN,ARGS="$a" leg.sbatch)
  echo "$(date -Is) $1 job $jid bin $BIN args: $a" | tee -a JOBS.txt
}
sub K5sp_sid_b4_k3 ns_params_toyoura_k3.json b4 1:00:00 3000
sub K5sp_sid_b8_k3 ns_params_toyoura_k3.json b8 1:00:00 3000
sub K5sp_sid_b16_k3 ns_params_toyoura_k3.json b16 4:00:00 13000
for r in 1.0 1.5 2.0; do sub K5sp_sid_b8_r$r ns_params_toyoura_k3_pi0r$r.json b8 1:00:00 3000; done
sub K5sp_sid_b8_r1.5_G3 ns_params_toyoura_k3_pi0r1.5_G0x3.json b8 2:00:00 6500
