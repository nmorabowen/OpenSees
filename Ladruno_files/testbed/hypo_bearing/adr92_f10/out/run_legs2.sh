#!/bin/sh
cd "$(dirname "$0")/.."
for L in I J K; do
  echo "########## leg $L start $(date +%H:%M:%S)"
  python3.12 -u f10_selfweight_wall.py leg --leg $L --h0 0.5 --wall 400 --out out > out/leg_$L.log 2>&1
  echo "########## leg $L done rc=$? $(date +%H:%M:%S)"
  grep -E "MODE =|s/B reached|refusals:|guards:" out/leg_$L.log
done
for L in A C H; do
  echo "########## census $L start $(date +%H:%M:%S)"
  python3.12 -u f10_selfweight_wall.py census --leg $L --h0 0.5 --census-steps 3 --census-ds 1e-3 --out out > out/census_$L.log 2>&1
  echo "########## census $L done rc=$?"
done
echo "########## census B start $(date +%H:%M:%S)"
python3.12 -u f10_selfweight_wall.py census --leg B --h0 0.5 --census-steps 3 --census-ds 1e-3 --out out > out/census_B.log 2>&1
echo "########## census B done rc=$?"
echo "########## ALL DONE 2"
