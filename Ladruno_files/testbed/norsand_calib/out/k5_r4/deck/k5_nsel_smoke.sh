#!/bin/bash
set -eu
cd $HOME/wp144/deck
export LD_LIBRARY_PATH=$HOME/wp144/lib
BIN=$HOME/wp144/bin/k5_0d38e57b4
S="--mat norsand --sand toyoura --energy HAR --target 0.01 --ckpt-every 50 --bproto 0.5 --dr 0.8342 --gamma 15.60 --dom-w 14 --dom-d 8 --ns-params $HOME/wp144/deck/ns_params_toyoura_k3.json --mesh b4 --wall 1500"
for e in qeas lst; do L=K5smoke_ns_$e; rm -rf runs/$L; mkdir -p runs/$L/logs
  sbatch --parsable --time=0:30:00 --partition=computes,new --exclude=node5 --mem=6G --job-name=sm_$L --export=ALL,ELEM=$e,LEG=$L,BIN=$BIN,ARGS="$S" leg_k5.sbatch; done
