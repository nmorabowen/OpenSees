#!/bin/bash
#SBATCH --exclude=node5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --partition=computes
#SBATCH --output=/mnt/deadmanschest/nmorabowen/wp144/deck/runs/%x-%j.out
# K5 comparator: LadrunoSANISAND, DM04 Toyoura set, the TIMs baseline flags (SAS-ME scheme 129, TolR 1e-4, R1 set),
# on the NorSand deck in Kimura mode. Env: LEG, BIN, KIM (test id), MESHARGS (e.g. "--mesh b8 --mesh-perturb shear:15").
set -u
D=$HOME/wp144/deck
mkdir -p $HOME/wp144/tmp
cd $D
exec $HOME/wp144/venv/bin/python -u $D/launch.py --leg "$LEG" --bin "$BIN" -- \
  --mat sanisand --scheme 129 --tolr 1e-4 --sand toyoura --kimura "$KIM" --target 0.20 --wall 80000 --ckpt-every 25 \
  $MESHARGS --extra "-sasHFloor 1.0 -sasReseatHyst 1.0 -sasSoftCap 0.5 -sasTensionCutoff 0.5 1.0"
