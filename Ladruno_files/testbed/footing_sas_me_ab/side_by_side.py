"""WP-138: step rows of two legs side by side at the s/B values they share.
    python side_by_side.py LEG1 LEG2"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(leg):
    return {round(float(r["s_over_B"]), 7): r for r in
            csv.DictReader(open(os.path.join(HERE, "runs", leg, "steps.csv"), newline=""))}


a, b = load(sys.argv[1]), load(sys.argv[2])
print(f"s/B       | {sys.argv[1]}: q rung it fails wall | {sys.argv[2]}: q rung it fails wall caps")
for s in sorted(b):
    if s in a:
        A, B = a[s], b[s]
        print(f"{s:.6f} | {float(A['q_kPa']):9.3f} {A['rung']} {A['iters']:>3} {A['fails_before']} "
              f"{float(A['wall_step_s']):7.1f} | {float(B['q_kPa']):9.3f} {B['rung']} "
              f"{B['iters']:>3} {B['fails_before']} {float(B['wall_step_s']):7.1f} {B['cap_step']}")
