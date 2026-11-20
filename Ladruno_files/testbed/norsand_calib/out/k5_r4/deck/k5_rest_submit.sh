#!/bin/bash
# K5, the rest (2026-11-11). Run with bash -lc on Esmeralda. Argument: ns | san | smoke
set -eu
cd $HOME/wp144/deck
export LD_LIBRARY_PATH=$HOME/wp144/lib
BIN=$HOME/wp144/bin/k5_0d38e57b4
NS="--mat norsand --sand toyoura --energy HAR --ns-params $HOME/wp144/deck/ns_params_toyoura_k3.json --target 0.20 --ckpt-every 50"
ns() {  # leg sbatchfile time args...
  local leg=$1 sb=$2 tm=$3; shift 3
  mkdir -p runs/$leg/logs
  local jid=$(sbatch --parsable --time=$tm --partition=computes,new --exclude=node5 --mem=6G --job-name=ns_$leg \
    --export=ALL,LEG=$leg,BIN=$BIN,ARGS="$NS $*" $sb)
  echo "$(date -Is) $leg job $jid bin $BIN args: $NS $*" | tee -a JOBS.txt
}
san() {  # leg kimura meshargs
  local leg=$1 kim=$2; shift 2
  mkdir -p runs/$leg/logs
  local jid=$(sbatch --parsable --time=26:00:00 --partition=computes,new --exclude=node5 --mem=6G --job-name=ns_$leg \
    --export=ALL,LEG=$leg,BIN=$BIN,KIM=$kim,MESHARGS="$*" leg_k5san.sbatch)
  echo "$(date -Is) $leg job $jid bin $BIN sanisand baseline, kimura $kim, $*" | tee -a JOBS.txt
}
case $1 in
smoke)
  $HOME/wp144/venv/bin/python k5_patch.py
  $HOME/wp144/venv/bin/python -m py_compile footing_ab_ns_k5.py
  sed s/footing_ab_ns.py/footing_ab_ns_k5.py/ launch.py > launch_k5.py
  sed s/launch.py/launch_k5.py/ leg.sbatch > leg_k5.sbatch
  ns K5r4_kimura_V75.1_b8_shear15 leg_k5.sbatch 1:00:00 --kimura V75.1 --mesh b8 --mesh-perturb shear:15 --wall 3000
  ;;
ns)
  for k in V85.6 V75.1 V64.6; do
    ns K5r4_kimura_${k}_b4_F leg.sbatch 1:00:00 --kimura $k --mesh b4 --wall 3000
    ns K5r4_kimura_${k}_b4_shear15 leg_k5.sbatch 1:00:00 --kimura $k --mesh b4 --mesh-perturb shear:15 --wall 3000
    [ $k = V75.1 ] || ns K5r4_kimura_${k}_b8_shear15 leg_k5.sbatch 1:00:00 --kimura $k --mesh b8 --mesh-perturb shear:15 --wall 3000
    ns K5r4_kimura_${k}_b16_shear15 leg_k5.sbatch 4:00:00 --kimura $k --mesh b16 --mesh-perturb shear:15 --wall 12000
  done
  ;;
sid)
  # Siddiquee et al. (1999): B 0.5 m, e0 0.66 (Dr 0.8342 with e_max 0.977, e_min 0.597), gamma = 2.64 x 9.81 / 1.66 = 15.60
  S="--bproto 0.5 --dr 0.8342 --gamma 15.60"
  ns K5r4_sid_B50_b8_F leg.sbatch 1:00:00 $S --mesh b8 --wall 3000
  ns K5r4_sid_B50_b8_Fhalf leg.sbatch 1:00:00 $S --mesh b8 --wall 3000 --pmin 0.2525
  ns K5r4_sid_B50_b8_shear15 leg_k5.sbatch 1:00:00 $S --mesh b8 --mesh-perturb shear:15 --wall 3000
  ns K5r4_sid_B50_b4_F leg.sbatch 1:00:00 $S --mesh b4 --wall 3000
  ns K5r4_sid_B50_b16_F leg.sbatch 4:00:00 $S --mesh b16 --wall 12000
  ;;
san)
  for k in V85.6 V75.1 V64.6; do
    san K5san_kimura_${k}_b8 $k --mesh b8
    san K5san_kimura_${k}_b8_shear15 $k --mesh b8 --mesh-perturb shear:15
  done
  ;;
esac
