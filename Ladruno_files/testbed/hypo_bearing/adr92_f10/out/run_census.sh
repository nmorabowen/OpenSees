#!/bin/sh
cd "$(dirname "$0")/.."
for L in B A C H; do
  for DS in 2e-5 4e-5 8e-5 2e-4; do
    echo "##### census $L ds=$DS $(date +%H:%M:%S)"
    python3.12 -u f10_selfweight_wall.py census --leg $L --h0 0.5 --census-steps 3 \
        --census-ds $DS --out out > out/census_${L}_ds${DS}.log 2>&1
    grep -E "ORIGIN: eta|^  -- step 3|error > 0.05|error > 0.1 |error:|>0.05 pop" out/census_${L}_ds${DS}.log | tail -8
  done
done
echo "##### CENSUS SWEEP DONE"
