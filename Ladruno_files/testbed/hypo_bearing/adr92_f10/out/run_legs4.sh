#!/bin/sh
cd "$(dirname "$0")/.."
for L in N N1; do
  echo "########## leg $L start $(date +%H:%M:%S)"
  python3.12 -u f10_selfweight_wall.py leg --leg $L --h0 0.5 --wall 400 --out out > out/leg_$L.log 2>&1
  echo "########## leg $L done rc=$? $(date +%H:%M:%S)"
  grep -E "MODE =|s/B reached|refusals:|guards:" out/leg_$L.log
done
echo "########## ALL DONE 4"
