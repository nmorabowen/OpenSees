#!/bin/bash
# O15 separating controls (2026-11-13): Drucker-Prager phi 38 on the Siddiquee pit box, deck copy footing_ab_ns_k5.py (leg_k5.sbatch)
#   K5dpa_* : ASSOCIATED (psi = phi, plane-strain matched), rough footing      -> DP_ASSOC=1
#   K5dps_* : non-associated (psi 0), SMOOTH footing (vertical DOF coupled)    -> FOOT_SMOOTH=1
set -eu
cd $HOME/wp144/deck
export LD_LIBRARY_PATH=$HOME/wp144/lib
BIN=$HOME/wp144/bin/k5_0d38e57b4
S="--mat dp --sand toyoura --dr 0.8342 --dp-phi 38 --dp-g 30000 --target 0.20 --ckpt-every 50 --bproto 0.5 --gamma 15.60 --dom-w 14 --dom-d 8"
sub() {  # leg mesh time wall envs
  mkdir -p runs/$1/logs
  local a="$S --mesh $2 --wall $4"
  local jid=$(sbatch --parsable --time=$3 --partition=computes,new --exclude=node5 --mem=6G --job-name=ns_$1 \
    --export=ALL,$5,LEG=$1,BIN=$BIN,ARGS="$a" leg_k5.sbatch)
  echo "$(date -Is) $1 job $jid env $5 args: $a" | tee -a JOBS.txt
}
for m in b4 b8; do
  sub K5dpa_sid_$m $m 1:00:00 3000 DP_ASSOC=1
  sub K5dps_sid_$m $m 1:00:00 3000 FOOT_SMOOTH=1
done
sub K5dpa_sid_b16 b16 4:00:00 13000 DP_ASSOC=1
sub K5dps_sid_b16 b16 4:00:00 13000 FOOT_SMOOTH=1
