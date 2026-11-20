#!/bin/bash
# K5 initial-surface test (2026-11-12): pi_i0 = ratio x p' on every row (deck pi0 mode 'ratio'), K3 constants held.
# Run with bash -lc on Esmeralda.
set -eu
cd $HOME/wp144/deck
export LD_LIBRARY_PATH=$HOME/wp144/lib
BIN=$HOME/wp144/bin/k5_0d38e57b4
for r in 1.0 1.5 2.0; do
  $HOME/wp144/venv/bin/python - $r <<'EOF'
import json, sys
r = float(sys.argv[1])
d = json.load(open("ns_params_toyoura_k3.json"))
d["name"] += f"_pi0ratio{sys.argv[1]}"
d["status"] = f"DIAGNOSTIC 2026-11-12: K3 file with pi0 mode ratio, pi_i0 = {sys.argv[1]} x p'; nothing else changed. " + d["status"]
d["pi0"]["mode"] = "ratio"
d["pi0"]["ratio"] = r
json.dump(d, open(f"ns_params_toyoura_k3_pi0r{sys.argv[1]}.json", "w"), indent=2)
EOF
done
sub() {  # leg ratio kimura
  mkdir -p runs/$1/logs
  local a="--mat norsand --sand toyoura --energy HAR --ns-params $HOME/wp144/deck/ns_params_toyoura_k3_pi0r$2.json --target 0.20 --ckpt-every 50 --kimura $3 --mesh b8 --wall 3000"
  local jid=$(sbatch --parsable --time=1:00:00 --partition=computes,new --exclude=node5 --mem=6G --job-name=ns_$1 \
    --export=ALL,LEG=$1,BIN=$BIN,ARGS="$a" leg.sbatch)
  echo "$(date -Is) $1 job $jid bin $BIN args: $a" | tee -a JOBS.txt
}
for k in V85.6 V75.1 V64.6; do for r in 1.0 1.5; do sub K5pi_kimura_${k}_b8_r$r $r $k; done; done
sub K5pi_kimura_V75.1_b8_r2.0 2.0 V75.1
