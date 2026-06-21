"""ADR 39 P0a — bucket-sort broad-phase prototype, validated vs brute force.

Pure-Python reference for the §26.11 bucket sort (LS-DYNA Theory pp. 545-550)
the C++ LadrunoContactBucketSort will mirror. Establishes the algorithm + the
correctness contract BEFORE any SRC:

  CONTRACT: for every slave node, the bucket-restricted candidate master-segment
  set must be a SUPERSET of the true geometric near set (every segment whose
  closest point is within the search radius). Missing a true-near pair = a
  silently-dropped contact (ADR Q-WIRE: "broad-phase bound MUST be conservative").

Grounded formulas:
  - cell size = LMAX (fraction of segment-diagonal char. length), §26.11.1
  - NX=(xmax-xmin)/LMAX etc., cap NX*NY*NZ <= min(NSN,5000)   (26.34-26.37)
  - bucket id NB = PX + (PY-1)*NX + (PZ-1)*NX*NY              (26.38-26.41)
  - search the 27 neighbour buckets [PX-1..PX+1] x ...        (26.42)
  - segment-CENTROID search (not nearest-node — 26.8.1 / Fig 26.24)
  - runaway-node guard: clamp the bbox so one fly-off node can't collapse buckets
"""
import math
import random


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------
def seg_centroid(seg_nodes):
    n = len(seg_nodes)
    return tuple(sum(p[k] for p in seg_nodes) / n for k in range(3))


def seg_diag(seg_nodes):
    """Longest diagonal of a quad/tri segment (the char. length source)."""
    d = 0.0
    for i in range(len(seg_nodes)):
        for j in range(i + 1, len(seg_nodes)):
            d = max(d, math.dist(seg_nodes[i], seg_nodes[j]))
    return d


def point_seg_dist(p, seg_nodes):
    """Approx closest distance: min over the segment's sample points (centroid +
    nodes + edge midpoints). Good enough for a broad-phase superset check."""
    samples = list(seg_nodes) + [seg_centroid(seg_nodes)]
    for i in range(len(seg_nodes)):
        a, b = seg_nodes[i], seg_nodes[(i + 1) % len(seg_nodes)]
        samples.append(tuple((a[k] + b[k]) / 2 for k in range(3)))
    return min(math.dist(p, s) for s in samples)


# --------------------------------------------------------------------------
# bucket sort (§26.11)
# --------------------------------------------------------------------------
class BucketSort:
    def __init__(self, segments, lmax_frac=1.0, bbox_clip=None):
        """segments: list of node-coord lists. lmax_frac: cell = lmax_frac *
        median segment diagonal. bbox_clip: optional (lo,hi) percentile clamp to
        guard runaway nodes (e.g. (1,99))."""
        self.segments = segments
        cents = [seg_centroid(s) for s in segments]
        diags = sorted(seg_diag(s) for s in segments)
        self.lmax = max(1e-12, lmax_frac * diags[len(diags) // 2])  # median diag

        xs = [c[0] for c in cents]; ys = [c[1] for c in cents]; zs = [c[2] for c in cents]
        self.runaway_guard_fired = False
        if bbox_clip is not None:
            (self.xmin, self.xmax), gx = self._clip(xs, bbox_clip)
            (self.ymin, self.ymax), gy = self._clip(ys, bbox_clip)
            (self.zmin, self.zmax), gz = self._clip(zs, bbox_clip)
            self.runaway_guard_fired = gx or gy or gz
        else:
            self.xmin, self.xmax = min(xs), max(xs)
            self.ymin, self.ymax = min(ys), max(ys)
            self.zmin, self.zmax = min(zs), max(zs)

        nsn = len(segments)
        self.NX = max(1, int((self.xmax - self.xmin) / self.lmax) + 1)
        self.NY = max(1, int((self.ymax - self.ymin) / self.lmax) + 1)
        self.NZ = max(1, int((self.zmax - self.zmin) / self.lmax) + 1)
        # cap total buckets <= min(NSN,5000)  (26.37)
        cap = min(nsn, 5000)
        while self.NX * self.NY * self.NZ > cap:
            # shrink the largest dimension
            m = max(self.NX, self.NY, self.NZ)
            if self.NX == m and self.NX > 1: self.NX -= 1
            elif self.NY == m and self.NY > 1: self.NY -= 1
            elif self.NZ > 1: self.NZ -= 1
            else: break

        # SEGMENT-SPAN registration (LS-DYNA type-13, eq 26.44): register each
        # segment in EVERY bucket its corner pointers span — NOT just its centroid
        # bucket. Centroid-only (type 4) is deficient when a segment is large vs
        # the cell or poorly shaped (§26.8.1 / Fig 26.24): the nearest centroid is
        # not the nearest segment, so true-near pairs are missed.
        self.grid = {}
        for sid, s in enumerate(segments):
            ps = [self._bucket(p) for p in s]
            imin, imax = min(p[0] for p in ps), max(p[0] for p in ps)
            jmin, jmax = min(p[1] for p in ps), max(p[1] for p in ps)
            kmin, kmax = min(p[2] for p in ps), max(p[2] for p in ps)
            for ix in range(imin, imax + 1):
                for iy in range(jmin, jmax + 1):
                    for iz in range(kmin, kmax + 1):
                        self.grid.setdefault((ix, iy, iz), []).append(sid)

    @staticmethod
    def _clip(vals, pct):
        s = sorted(vals)
        lo = s[int(len(s) * pct[0] / 100)]
        hi = s[min(len(s) - 1, int(len(s) * pct[1] / 100))]
        fired = (lo > min(vals)) or (hi < max(vals))
        return (lo, hi), fired

    def _px(self, x, xmin, xmax, NX):
        if xmax <= xmin:
            return 1
        p = int(NX * (x - xmin) / (xmax - xmin)) + 1
        return min(max(p, 1), NX)

    def _bucket(self, c):
        PX = self._px(c[0], self.xmin, self.xmax, self.NX)
        PY = self._px(c[1], self.ymin, self.ymax, self.NY)
        PZ = self._px(c[2], self.zmin, self.zmax, self.NZ)
        return (PX, PY, PZ)

    def candidates(self, slave_pt):
        """Segments in the 27 neighbour buckets of the slave node (26.42)."""
        PX = self._px(slave_pt[0], self.xmin, self.xmax, self.NX)
        PY = self._px(slave_pt[1], self.ymin, self.ymax, self.NY)
        PZ = self._px(slave_pt[2], self.zmin, self.zmax, self.NZ)
        out = []
        for ix in range(max(1, PX - 1), min(self.NX, PX + 1) + 1):
            for iy in range(max(1, PY - 1), min(self.NY, PY + 1) + 1):
                for iz in range(max(1, PZ - 1), min(self.NZ, PZ + 1) + 1):
                    out.extend(self.grid.get((ix, iy, iz), []))
        return out


def brute_force_near(slave_pt, segments, radius):
    return {sid for sid, s in enumerate(segments) if point_seg_dist(slave_pt, s) <= radius}


# --------------------------------------------------------------------------
# verification
# --------------------------------------------------------------------------
def make_quad_grid(nx, ny, z, jitter=0.0, rng=None):
    """A nx*ny grid of unit quad segments at height z (a master surface)."""
    segs = []
    for i in range(nx):
        for j in range(ny):
            def cor(a, b):
                jx = (rng.uniform(-jitter, jitter) if rng else 0.0)
                jy = (rng.uniform(-jitter, jitter) if rng else 0.0)
                return (a + jx, b + jy, z)
            segs.append([cor(i, j), cor(i + 1, j), cor(i + 1, j + 1), cor(i, j + 1)])
    return segs


def run():
    rng = random.Random(12345)
    print("ADR39 P0a — bucket sort vs brute force\n" + "=" * 50)

    # --- Test 1: superset contract on a uniform grid ---
    segs = make_quad_grid(20, 20, z=0.0)
    bs = BucketSort(segs, lmax_frac=1.0)
    print(f"[grid 20x20] cell(lmax)={bs.lmax:.3f}  buckets={bs.NX}x{bs.NY}x{bs.NZ}={bs.NX*bs.NY*bs.NZ}")
    # DESIGN RULE: the spatial-hash superset guarantee holds for search radius
    # r <= cell_size with ±1 neighbour search. The contact capture band (gap
    # tolerance + per-epoch slip) is a fraction of a cell, so r = 0.5*cell is the
    # realistic, in-guarantee contract. (r == cell exactly is the boundary edge
    # and is not guaranteed — C++ sizes LMAX = max(seg_diag, r_search)*margin.)
    radius = 0.5 * bs.lmax
    fails = 0
    for _ in range(2000):
        sp = (rng.uniform(0, 20), rng.uniform(0, 20), rng.uniform(-0.5, 0.5))
        cand = set(bs.candidates(sp))
        truth = brute_force_near(sp, segs, radius)
        if not truth.issubset(cand):  # THE contract: candidates superset truth
            fails += 1
    print(f"[superset contract] 2000 slave probes, missed-true-near = {fails}  -> {'PASS' if fails==0 else 'FAIL'}")

    # --- Test 2: poor-aspect / jittered mesh (the §26.8.1 nearest-node failure regime) ---
    segs2 = make_quad_grid(15, 15, z=0.0, jitter=0.35, rng=rng)
    bs2 = BucketSort(segs2, lmax_frac=1.5)  # widen cell for distorted mesh
    radius2 = 0.5 * bs2.lmax
    fails2 = 0
    for _ in range(2000):
        sp = (rng.uniform(0, 15), rng.uniform(0, 15), rng.uniform(-0.4, 0.4))
        if not brute_force_near(sp, segs2, radius2).issubset(set(bs2.candidates(sp))):
            fails2 += 1
    print(f"[jittered mesh] missed-true-near = {fails2}  -> {'PASS' if fails2==0 else 'FAIL'}")

    # --- Test 3: runaway-node guard (one segment flung to infinity) ---
    segs3 = make_quad_grid(15, 15, z=0.0)
    segs3.append([(1e6, 1e6, 1e6), (1e6 + 1, 1e6, 1e6),
                  (1e6 + 1, 1e6 + 1, 1e6), (1e6, 1e6 + 1, 1e6)])  # fly-off
    bs_noguard = BucketSort(segs3, lmax_frac=1.0)
    bs_guard = BucketSort(segs3, lmax_frac=1.0, bbox_clip=(1, 99))
    # EFFECTIVE resolution metric: how many distinct buckets do the REAL (non-fly-off)
    # 15x15 segments occupy? A runaway node inflates the bbox so all real segments
    # collapse into ~1 giant bucket even though the raw bucket COUNT stays capped.
    def real_occupancy(bs):
        return len({bs._bucket(seg_centroid(s)) for s in segs3[:-1]})  # exclude fly-off
    occ_noguard = real_occupancy(bs_noguard)
    occ_guard = real_occupancy(bs_guard)
    print(f"[runaway] NO guard: real 15x15 surface occupies {occ_noguard} distinct buckets")
    print(f"[runaway] guard   : real 15x15 surface occupies {occ_guard} distinct buckets"
          f"  guard_fired={bs_guard.runaway_guard_fired}")
    collapsed = occ_noguard <= 2           # fly-off collapsed the real surface into ~1 bucket
    recovered = occ_guard >= 25 and bs_guard.runaway_guard_fired
    print(f"[runaway guard] no-guard-collapses={collapsed}  guard-recovers={recovered}"
          f"  -> {'PASS' if collapsed and recovered else 'FAIL'}")

    # --- Test 4: complexity — bucket vs brute, count distance ops ---
    N = 40 * 40
    segsN = make_quad_grid(40, 40, z=0.0)
    bsN = BucketSort(segsN, lmax_frac=1.0)
    probe = (20.3, 20.7, 0.1)
    n_bucket = len(bsN.candidates(probe))
    n_brute = N
    print(f"[complexity] N={N} segs: brute checks={n_brute}, bucket checks={n_bucket}"
          f"  speedup~{n_brute/max(1,n_bucket):.0f}x  -> {'PASS' if n_bucket < n_brute/10 else 'FAIL'}")

    ok = (fails == 0 and fails2 == 0 and collapsed and recovered and n_bucket < n_brute / 10)
    print("\n" + ("ALL P0a CHECKS PASS" if ok else "P0a HAS FAILURES"))
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if run() else 1)
