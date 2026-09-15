#!/bin/sh
cd "$(dirname "$0")/.."
for DS in 8e-5 4e-5 2e-5 1e-5 5e-6 2e-6 5e-7; do
  echo "##### probe B s=0.0085 ds=$DS $(date +%H:%M:%S)"
  python3.12 -u f10_selfweight_wall.py probe --leg B --h0 0.5 --probe-s 0.0085 \
      --census-ds $DS --wall 600 --out out > out/probe_B_ds$DS.log 2>&1
  grep -E "walked|PROBE|worst" out/probe_B_ds$DS.log | head -4
done
echo "##### PROBE DONE"
