#!/bin/bash
# Start ONE WP-138 leg as an exclusive single-node srun step, detached from the
# ssh session (setsid nohup), and record the SLURM job id.
#   bash start_leg.sh <leg> <footing_ab args...>
set -u
LEG=$1; shift
export PATH=/opt/slurm/bin:$PATH
export TMPDIR=$HOME/ladruno_build_test/tmp
D=$HOME/ladruno_wp138/deck
mkdir -p $D/runs/$LEG/logs
cd $D
setsid nohup srun --job-name=wp138_$LEG --nodes=1 --ntasks=1 --cpus-per-task=32 --exclusive \
  $HOME/ladruno_build_test/conan_venv/bin/python -u $D/launch.py --leg $LEG -- "$@" \
  > $D/runs/$LEG/logs/srun.out 2>&1 < /dev/null &
sleep 4
JID=$(squeue -h -u $USER -n wp138_$LEG -o "%i %N" | head -1)
echo "$(date -Is) $LEG job/node: $JID args: $*" | tee -a $D/JOBS.txt
