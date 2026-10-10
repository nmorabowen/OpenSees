#!/bin/bash
# Pull logs + steps.csv + summaries (NOT the ckpt/replay field dumps) back to SCR.
S="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$S/runs"
ssh -o BatchMode=yes esmeralda 'cd ~/ladruno_wp138/deck && tar czf - JOBS.txt $(ls -d runs/*/logs runs/*/steps.csv runs/*/summary.json runs/*/census_last_converged.csv census 2>/dev/null)' | tar xzf - -C "$S"
ssh -o BatchMode=yes esmeralda 'export PATH=/opt/slurm/bin:$PATH; squeue -u $USER -o "%.8i %.16j %.2t %.10M %N"' > "$S/squeue_last.txt"
date -Is >> "$S/pulls.txt"
