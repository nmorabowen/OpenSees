#!/bin/bash
#SBATCH --exclude=node5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --partition=computes
#SBATCH --output=/mnt/deadmanschest/nmorabowen/wp144/deck/runs/%x-%j.out
# One strip leg on one shared compute node (house style of ~/ladruno_wp152/deck_*/start_leg_shared.sh: non-exclusive, 4 CPUs,
# 6 GB, single-threaded MKL/OMP). Do not call directly: use submit_leg.sh <leg> <footing_ab_ns.py args...>
set -u
D=$HOME/wp144/deck
mkdir -p $HOME/wp144/tmp
cd $D
exec $HOME/wp144/venv/bin/python -u $D/launch_k5.py --leg "$LEG" --bin "${BIN:-pdl}" -- $ARGS
