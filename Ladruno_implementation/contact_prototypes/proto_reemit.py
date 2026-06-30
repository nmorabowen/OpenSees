#!/usr/bin/env python3
"""ADR-60 P1 oracle — the no-pass-through contract for finite-sliding NTS re-emission.

Build-free (numpy only). Demonstrates the bug and the fix on a flat master strip:

  * REFERENCE-coord broad phase (today's LadrunoContactBucketSort, built once at
    handle() from getCrds()): as the slave slides in +x across the strip, its
    candidate segment set is frozen at the reference cell -> the segment it is
    ACTUALLY over (deformed) eventually leaves the candidate set -> SILENT
    PASS-THROUGH (no adapter, zero contact force).

  * DEFORMED-coord RE-EMIT (ADR-60): rebuild the candidate set from the committed
    position each re-emit -> the active segment is ALWAYS in the candidate set
    (the no-pass-through contract, ADR-60 P1 gate (a)).

This mirrors the bucket-sort cell logic (LadrunoContactBucketSort.h): cell size =
cellFrac * median segment diagonal, 27-neighbour (here 3-neighbour in 1D-of-interest)
query. Sibling of proto_bucket_sort.py and proto_reemit_selfcheck.cpp.

Run:  python proto_reemit.py
"""
import numpy as np

CELL_FRAC = 1.0          # matches the handler default ct.cellFrac
SEG_LEN   = 1.0          # master segments are unit squares in the xy-plane
N_SEG     = 12           # strip of 12 quads along x


def master_strip(n=N_SEG, L=SEG_LEN):
    """n quad segments tiling [0,n*L] x [0,L] in z=0. Returns (n,4,3) ref coords."""
    segs = np.zeros((n, 4, 3))
    for s in range(n):
        x0 = s * L
        segs[s] = [[x0, 0, 0], [x0 + L, 0, 0], [x0 + L, L, 0], [x0, L, 0]]
    return segs


def median_diag(segs):
    diags = [max(np.linalg.norm(segs[s, i] - segs[s, j])
                 for i in range(4) for j in range(i + 1, 4))
             for s in range(len(segs))]
    return float(np.median(diags))


def candidate_segments(slave_xy, seg_centroids, cell):
    """Segments whose centroid is within the +-1 cell band of the slave (the
    bucket-sort 27-neighbour query, reduced to the relevant axes)."""
    d = np.abs(seg_centroids - slave_xy)              # (n,3)
    within = np.all(d <= 1.5 * cell, axis=1)          # 1.5*cell ~ the +-1 neighbour reach
    return set(np.nonzero(within)[0].tolist())


def active_segment(slave_xy, segs):
    """The segment the slave is geometrically over (in-bounds projection)."""
    for s in range(len(segs)):
        x0, x1 = segs[s, 0, 0], segs[s, 1, 0]
        y0, y1 = segs[s, 0, 1], segs[s, 3, 1]
        if x0 - 1e-9 <= slave_xy[0] <= x1 + 1e-9 and y0 - 1e-9 <= slave_xy[1] <= y1 + 1e-9:
            return s
    return None


def main():
    segs = master_strip()
    ref_centroids = segs.mean(axis=1)                 # (n,3)
    cell = CELL_FRAC * median_diag(segs)

    # slave starts over segment 1, slides in +x to over segment 10 (committed configs).
    y = 0.5 * SEG_LEN
    slave_path = [np.array([1.5 + t, y, 0.0]) for t in np.linspace(0.0, 8.0, 33)]

    # --- REFERENCE-coord broad phase: candidates frozen at the START position ---
    ref_anchor = slave_path[0]
    ref_cands = candidate_segments(ref_anchor, ref_centroids, cell)

    # --- DEFORMED-coord RE-EMIT: ADR-60 migration trigger (re-sort when drift > f*band) ---
    f, band = 0.5, cell
    anchor = slave_path[0]
    deformed_cands = candidate_segments(anchor, ref_centroids, cell)

    ref_miss, reemit_miss, n_resorts = 0, 0, 0
    for x in slave_path:
        act = active_segment(x, segs)
        if act is None:
            continue
        # reference path: never re-sorts
        if act not in ref_cands:
            ref_miss += 1
        # re-emit path: migration trigger
        if np.linalg.norm(x - anchor) > f * band:
            anchor = x
            deformed_cands = candidate_segments(anchor, ref_centroids, cell)
            n_resorts += 1
        if act not in deformed_cands:
            reemit_miss += 1

    print(f"cell size            : {cell:.4f}")
    print(f"reference-path misses: {ref_miss}  (slave slid off its frozen candidate set -> PASS-THROUGH)")
    print(f"re-emit-path misses  : {reemit_miss}  (re-sorts: {n_resorts})")

    ok = True
    # The bug must be reproduced: the frozen reference candidate set DOES lose the slave.
    if ref_miss == 0:
        print("[FAIL] expected the reference path to pass-through, but it never missed")
        ok = False
    else:
        print("[PASS] reference path reproduces the silent pass-through bug")
    # The fix must hold: re-emit never loses the active segment (no-pass-through).
    if reemit_miss == 0:
        print("[PASS] re-emit: active segment ALWAYS in the candidate set (no-pass-through)")
    else:
        print(f"[FAIL] re-emit missed the active segment {reemit_miss}x")
        ok = False

    print("\nALL PASS" if ok else "\nFAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
