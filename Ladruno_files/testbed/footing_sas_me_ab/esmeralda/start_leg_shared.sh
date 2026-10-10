#!/bin/bash
# Like start_leg.sh, but a NON-exclusive 4-CPU srun step pinned to a node, so
# several single-threaded legs (MKL/OMP threads = 1) can share one node.
#   bash start_leg_shared.sh <node> <leg> <footing_ab args...>
set -u
NODE=$1; shift
LEG=$1; shift
export PATH=/opt/slurm/bin:$PATH
export TMPDIR=$HOME/ladruno_build_test/tmp
D=$HOME/ladruno_wp138/deck
mkdir -p $D/runs/$LEG/logs
cd $D
setsid nohup srun --job-name=wp138_$LEG --nodes=1 --ntasks=1 --cpus-per-task=4 --nodelist=$NODE \
  $HOME/ladruno_build_test/conan_venv/bin/python -u $D/launch.py --leg $LEG -- "$@" \
  > $D/runs/$LEG/logs/srun.out 2>&1 < /dev/null &
sleep 4
JID=$(squeue -h -u $USER -n wp138_$LEG -o "%i %N %T" | head -1)
echo "$(date -Is) $LEG job/node: $JID (shared node, 4 cpus) args: $*" | tee -a $D/JOBS.txt
