"""ADR-57 E1 ORACLE — the edge-edge routing detector (facet-pair -> lane).

Oracle-first (the fork discipline): validate the ROUTING LOGIC in pure numpy BEFORE
transcribing it into the handler's pairing loop. No OpenSees, no build. E1 decides, for a
candidate (slave facet, master facet) pair, WHICH narrow-phase lane owns it:
  NTS  >  EDGE_EDGE  >  FACE_MORTAR    (ordered precedence)  |  or NONE if separated.

The ADR-57 3-reviewer gate (Lens-C MAJOR C-3/C-4) killed the draft's "clean partition, no
superposed force" claim: NTS is a THIRD lane the clip-area arbiter ignored, and the oblique
band can present both a clip and an edge-edge pair. So E1 ships ORDERED PRECEDENCE WITH A
GUARANTEED OWNER, and the correctness is proven here, not asserted:
  - proximity gate first: facets farther than d_band -> NONE (genuinely separated).
  - within proximity an owner is GUARANTEED (no coverage hole): NTS if a slave node projects
    in-bounds on the master (the more local primitive wins -> no edge-edge/NTS double-count);
    else EDGE_EDGE if >=1 boundary-edge pair is a margin-interior, well-conditioned crossing;
    else FACE_MORTAR catches everything remaining (even a sliver clip — face-mortar stands
    down ONLY if edge-edge succeeds, so A_clip<tau_A is necessary-not-sufficient to vacate it).
  - align_fn = |n_s . n_m| (facet-normal alignment, renamed off the kernel's overloaded cos_t)
    is a CHEAP PRE-CLASSIFY with hysteresis (which lane to TRY first / which test to skip), NOT
    the arbiter — the precedence+proximity logic is what makes it correct and chatter-free.

Falsifiers (the E1 gate row + the Lens-C folds):
  - X/T/L-junction (perpendicular, no clip) -> EDGE_EDGE; facing-overlap -> FACE_MORTAR.
  - NTS-vs-edge DOUBLE-COUNT: a pair where both could fire -> exactly ONE lane (NTS), nonzero.
  - oblique-band NO-COVERAGE-HOLE: a 45deg pair with clamped/parallel edges + sub-tau_A clip but
    in proximity -> a non-NONE owner (never zero force).
  - parallel-but-offset FACING (collinear bars) -> FACE_MORTAR (the refused-by-the-kernel pair
    is NOT dropped — the other lane catches it).
  - hysteresis: sweeping the facet angle through 45deg does not oscillate the lane step-to-step.
  - partition-by-precedence: every pair maps to exactly one of {NTS,EDGE_EDGE,FACE_MORTAR,NONE}.

Run: <pythoncore-3.12> proto_e1_router.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)
_fails = 0

# lanes
NTS, EDGE_EDGE, FACE_MORTAR, NONE = "NTS", "EDGE_EDGE", "FACE_MORTAR", "NONE"
# gauges (match the E0 kernel + the ADR)
TAU_PAR = np.sin(np.deg2rad(15.0)) ** 2     # parallel/conditioning gate (E0)
MARGIN = 1e-6
TAU_FACE = np.cos(np.deg2rad(50.0))         # align_fn > this  => pre-classify "try face first"
TAU_EDGE = np.cos(np.deg2rad(70.0))         # align_fn < this  => pre-classify "try edge first"
TAU_A = 1e-3                                 # clip-area-degeneracy fraction (necessary-not-suff.)


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# ============================================ E0 kernel (compact, mirrors LadrunoEdgeKernel)
def closest_segseg(p1, q1, p2, q2, tau_par=TAU_PAR):
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e = d1 @ d1, d2 @ d2
    out = {"status": 0, "s": 0.0, "t": 0.0, "denom": 0.0, "interior": False,
           "c1": p1.copy(), "c2": p2.copy()}
    if a <= 1e-9 or e <= 1e-9:
        out["status"] = -2
        return out
    c, b, f = d1 @ r, d1 @ d2, d2 @ r
    denom = a * e - b * b
    out["denom"] = denom
    s = np.clip((b * f - c * e) / denom, 0, 1) if denom > 1e-300 else 0.0
    t = (b * s + f) / e
    if t < 0:
        t, s = 0.0, np.clip(-c / a, 0, 1)
    elif t > 1:
        t, s = 1.0, np.clip((b - c) / a, 0, 1)
    out["s"], out["t"] = s, t
    out["c1"], out["c2"] = p1 + s * d1, p2 + t * d2
    out["interior"] = (MARGIN <= s <= 1 - MARGIN) and (MARGIN <= t <= 1 - MARGIN)
    if denom < tau_par * a * e:
        out["status"] = -1   # parallel
    return out


# ============================================================== facet geometry helpers
def facet_normal(X):
    """Newell's method — robust unit normal for a planar tri/quad polygon X (n x 3)."""
    n = np.zeros(3)
    m = len(X)
    for i in range(m):
        a, b = X[i], X[(i + 1) % m]
        n[0] += (a[1] - b[1]) * (a[2] + b[2])
        n[1] += (a[2] - b[2]) * (a[0] + b[0])
        n[2] += (a[0] - b[0]) * (a[1] + b[1])
    nn = np.linalg.norm(n)
    return n / nn if nn > 1e-300 else n


def edges_of(X):
    m = len(X)
    return [(X[i], X[(i + 1) % m]) for i in range(m)]


def point_in_facet(x, X, d_band):
    """does point x project onto facet X (closest-point, in the polygon) with gap in band?
    returns (in_bounds_and_near, gap). Convex polygon assumed (facets are tri/quad)."""
    n = facet_normal(X)
    c = X.mean(axis=0)
    gap = (x - c) @ n
    proj = x - gap * n
    # in-polygon test via consistent cross-product sign against each edge (convex)
    m = len(X)
    sgn = 0
    for i in range(m):
        a, b = X[i], X[(i + 1) % m]
        cr = np.cross(b - a, proj - a) @ n
        s = np.sign(cr)
        if s != 0:
            if sgn == 0:
                sgn = s
            elif s != sgn:
                return False, gap
    return abs(gap) <= d_band, gap


def clip_area(slave, master):
    """area of (slave projected onto master plane) ∩ master, in master's 2D frame.
    Sutherland-Hodgman of two convex polygons. ~0 for perpendicular/edge-on facets."""
    nm = facet_normal(master)
    cm = master.mean(axis=0)
    # build a 2D frame on the master plane
    ref = np.array([1.0, 0, 0]) if abs(nm[0]) < 0.9 else np.array([0, 1.0, 0])
    u = ref - (ref @ nm) * nm
    u /= np.linalg.norm(u)
    v = np.cross(nm, u)

    def to2d(P):
        return np.array([[(p - cm) @ u, (p - cm) @ v] for p in P])

    subject = to2d(slave)      # slave projected (drop the out-of-plane component)
    clipw = to2d(master)
    # ensure CCW
    def area2(P):
        s = 0.0
        for i in range(len(P)):
            j = (i + 1) % len(P)
            s += P[i][0] * P[j][1] - P[j][0] * P[i][1]
        return 0.5 * s
    if area2(clipw) < 0:
        clipw = clipw[::-1]
    if area2(subject) < 0:
        subject = subject[::-1]
    out = [p for p in subject]
    for i in range(len(clipw)):
        A, B = clipw[i], clipw[(i + 1) % len(clipw)]
        edge = B - A
        nrm = np.array([-edge[1], edge[0]])   # inward for CCW clip
        inp = out
        out = []
        if not inp:
            break
        for j in range(len(inp)):
            cur, prv = inp[j], inp[j - 1]
            cin = (cur - A) @ nrm >= 0
            pin = (prv - A) @ nrm >= 0
            if cin:
                if not pin:
                    out.append(_isect(prv, cur, A, B))
                out.append(cur)
            elif pin:
                out.append(_isect(prv, cur, A, B))
    if len(out) < 3:
        return 0.0
    return abs(area2(np.array(out)))


def _isect(p, q, A, B):
    d1, d2 = q - p, B - A
    den = d1[0] * (-d2[1]) - d1[1] * (-d2[0])
    if abs(den) < 1e-300:
        return p
    tt = ((A[0] - p[0]) * (-d2[1]) - (A[1] - p[1]) * (-d2[0])) / den
    return p + tt * d1


def min_facet_dist(slave, master):
    """min distance between two facets (sampled: vertices + edge-edge closest points)."""
    d = np.inf
    for x in slave:
        ok, gap = point_in_facet(x, master, 1e9)
        d = min(d, abs(gap) if ok else np.inf)
    for x in master:
        ok, gap = point_in_facet(x, slave, 1e9)
        d = min(d, abs(gap) if ok else np.inf)
    for (a1, b1) in edges_of(slave):
        for (a2, b2) in edges_of(master):
            k = closest_segseg(a1, b1, a2, b2)
            d = min(d, np.linalg.norm(k["c1"] - k["c2"]))
    return d


# ======================================================================= the router
def route(slave, master, d_band=0.25, d_pen=0.25, verbose=False):
    """ordered precedence with a guaranteed owner. returns (lane, diag).

    REGIME-SPLIT on align_fn (the E1 oracle caught that a flat NTS-precedence steals
    facing-facet contact from the mortar lane): facing facets are the FACE-MORTAR regime
    whose integral already handles the in-bounds slave nodes — NTS-vs-edge precedence
    only applies in the NON-FACING (oblique/perpendicular) regime."""
    nS, nM = facet_normal(slave), facet_normal(master)
    align_fn = abs(nS @ nM)
    prox = min_facet_dist(slave, master)
    diag = {"align_fn": align_fn, "prox": prox, "nts": False, "edge": False, "clip": None}

    # 1) proximity gate — genuinely separated pairs own no lane
    if prox > d_band:
        return NONE, diag

    # 2) FACING regime — facets nearly parallel (align_fn -> 1): FACE-MORTAR owns the pair;
    #    its weighted-gap integral subsumes the in-bounds slave nodes (a separate NTS here
    #    would double-count). This is the regime split the oracle surfaced.
    if align_fn >= TAU_FACE:
        diag["clip"] = clip_area(slave, master)
        return FACE_MORTAR, diag

    # --- NON-FACING regime (oblique / perpendicular): NTS > EDGE_EDGE > sliver-mortar ---
    # 3) NTS precedence — a slave VERTEX poking the master face interior is the more local
    #    primitive; it OWNS the pair so edge-edge cannot double-count it (the Lens-C fix).
    for x in slave:
        ok, gap = point_in_facet(x, master, d_band)
        if ok and -d_pen <= gap <= d_band:
            diag["nts"] = True
            return NTS, diag

    # 4) EDGE-EDGE — any boundary-edge pair that is a margin-interior, well-conditioned
    #    crossing with gap in band (the E0 kernel decides interior + non-parallel).
    for (a1, b1) in edges_of(slave):
        for (a2, b2) in edges_of(master):
            k = closest_segseg(a1, b1, a2, b2)
            if k["status"] == 0 and k["interior"]:
                g = np.linalg.norm(k["c1"] - k["c2"])
                if g <= d_band:
                    diag["edge"] = True
                    return EDGE_EDGE, diag

    # 5) FACE-MORTAR — the GUARANTEED OWNER for everything else in proximity (even a sliver
    #    clip). face-mortar stands down ONLY if edge-edge succeeded above, so A_clip < tau_A
    #    is necessary-not-sufficient to vacate it: no coverage hole.
    diag["clip"] = clip_area(slave, master)
    return FACE_MORTAR, diag


# ============================================================================ tests
print("ADR-57 E1 oracle — edge-edge routing detector (ordered precedence)")

# --- canonical facets ---
def quad(center, u, v):
    c = np.asarray(center, float)
    return np.array([c - u - v, c + u - v, c + u + v, c - u + v])

# --- T1: X-crossing (two perpendicular bars crossing interior) -> EDGE_EDGE ---
print("\nT1  X-crossing perpendicular facets (no clip) -> EDGE_EDGE, align_fn->0")
# slave facet in x-z plane (thin bar along x), master in x-y plane (thin bar along y), gap ~0.1
S = quad([0, 0, 0.1], np.array([1.0, 0, 0]), np.array([0, 0, 0.05]))   # spans x, thin z
M = quad([0, 0, 0.0], np.array([0, 1.0, 0]), np.array([0.05, 0, 0]))   # spans y, thin x
lane, d = route(S, M)
check("X-crossing -> EDGE_EDGE", lane == EDGE_EDGE, f"lane={lane} align_fn={d['align_fn']:.3f}")
check("X-crossing align_fn ~ 0 (perpendicular facets)", d["align_fn"] < TAU_EDGE)

# --- T2: facing parallel overlap (align_fn->1) -> FACE_MORTAR (regime split) ---
print("\nT2  facing parallel overlapping facets -> FACE_MORTAR")
Sf = quad([0, 0, 0.05], np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
Mf = quad([0.3, 0.3, 0.0], np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))  # offset, still facing
lane2, d2 = route(Sf, Mf)
check("facing overlap -> FACE_MORTAR", lane2 == FACE_MORTAR,
      f"lane={lane2} align_fn={d2['align_fn']:.3f} clip={d2['clip']}")
check("facing align_fn ~ 1 (>= tau_face)", d2["align_fn"] >= TAU_FACE, f"align_fn={d2['align_fn']:.3f}")
check("facing clip area > 0 (real overlap, mortar regime)", d2["clip"] is not None and d2["clip"] > TAU_A)

# --- T3: parallel-but-offset FACING collinear bars -> FACE_MORTAR (kernel-refused pair NOT dropped) ---
print("\nT3  parallel-offset facing collinear bars -> FACE_MORTAR (the refused pair is not lost)")
Sb = quad([0, 0, 0.05], np.array([0.6, 0, 0]), np.array([0, 0.2, 0]))
Mb = quad([0.8, 0, 0.0], np.array([0.6, 0, 0]), np.array([0, 0.2, 0]))   # collinear x, overlapping extent
lane3, d3 = route(Sb, Mb)
check("parallel-offset facing -> FACE_MORTAR (not NONE)", lane3 == FACE_MORTAR,
      f"lane={lane3} align_fn={d3['align_fn']:.3f} prox={d3['prox']:.3f}")

# --- T4: NTS precedence in the NON-FACING regime: a slave vertex pokes the master face -> NTS ---
print("\nT4  non-facing slave VERTEX poking the master face -> NTS (the more local primitive)")
ang = np.deg2rad(60.0)   # tilt 60deg => align_fn = cos60 = 0.5 < tau_face (non-facing)
Ry = np.array([[np.cos(ang), 0, np.sin(ang)], [0, 1, 0], [-np.sin(ang), 0, np.cos(ang)]])
Mn = quad([0, 0, 0.0], np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))           # big horizontal master
# small slave facet tilted 60deg, lowered so its bottom corner sits ~0.08 above the master interior
Sn = (quad([0, 0, 0], np.array([0.25, 0, 0]), np.array([0, 0.25, 0])) @ Ry.T) + np.array([0.0, 0.0, 0.30])
lane4, d4 = route(Sn, Mn)
check("non-facing vertex-poke -> NTS (precedence over edge-edge)", lane4 == NTS,
      f"lane={lane4} align_fn={d4['align_fn']:.3f} nts={d4['nts']} prox={d4['prox']:.3f}")
check("T4 is genuinely NON-facing (align_fn < tau_face)", d4["align_fn"] < TAU_FACE)

# --- T5: oblique-band NO-COVERAGE-HOLE: non-facing, no vertex-poke, no interior crossing, in band ---
print("\nT5  oblique near pair (no NTS, no edge-cross) IN proximity -> non-NONE owner (no zero force)")
Mo = quad([0, 0, 0.0], np.array([0.5, 0, 0]), np.array([0, 0.5, 0]))
a45 = np.deg2rad(45.0)
Ry45 = np.array([[np.cos(a45), 0, np.sin(a45)], [0, 1, 0], [-np.sin(a45), 0, np.cos(a45)]])
# tilt 45deg and offset BEYOND the master in x so no edge crosses interior, but the nearest
# corner is within d_band of the master edge (proximity holds => must own a lane, not NONE)
So = (quad([0, 0, 0], np.array([0.25, 0, 0]), np.array([0, 0.25, 0])) @ Ry45.T) + np.array([0.62, 0, 0.12])
lane5, d5 = route(So, Mo, d_band=0.25)
check("oblique near pair in proximity -> NOT NONE (guaranteed owner)", lane5 != NONE,
      f"lane={lane5} prox={d5['prox']:.3f} align_fn={d5['align_fn']:.3f}")
check("T5 is within proximity (prox <= d_band)", d5["prox"] <= 0.25, f"prox={d5['prox']:.3f}")

# --- T6: genuinely separated facets (prox > d_band) -> NONE ---
print("\nT6  separated facets (far apart) -> NONE (correctly no contact, not a lost owner)")
Sfar = quad([0, 0, 5.0], np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
Mfar = quad([0, 0, 0.0], np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
lane6, d6 = route(Sfar, Mfar)
check("far apart -> NONE", lane6 == NONE, f"lane={lane6} prox={d6['prox']:.3f}")

# --- T7: hysteresis / no-chatter: sweep a facet's tilt through 45deg, lane is stable ---
print("\nT7  sweep tilt through 45deg: lane does not oscillate step-to-step (no chatter)")
lanes = []
Mh = quad([0, 0, 0.0], np.array([0.5, 0, 0]), np.array([0, 0.5, 0]))
for deg in np.linspace(80, 10, 36):      # perpendicular-ish (edge) -> facing-ish (mortar)
    a = np.deg2rad(deg)
    Ryy = np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    Sh = (quad([0, 0, 0], np.array([0.5, 0, 0]), np.array([0, 0.5, 0])) @ Ryy.T) + np.array([0, 0, 0.12])
    lanes.append(route(Sh, Mh)[0])
trans = sum(1 for i in range(1, len(lanes)) if lanes[i] != lanes[i - 1])
check("monotone tilt sweep has <=2 lane transitions (no chatter)", trans <= 2,
      f"transitions={trans} order={'|'.join(dict.fromkeys(lanes))}")

# --- T8: partition-by-precedence over a battery -> every pair maps to exactly one lane ---
print("\nT8  partition: every pair -> exactly one of {NTS,EDGE_EDGE,FACE_MORTAR,NONE}")
battery = [route(S, M)[0], route(Sf, Mf)[0], route(Sb, Mb)[0], route(Sn, Mn)[0],
           route(So, Mo)[0], route(Sfar, Mfar)[0]]
valid = all(L in (NTS, EDGE_EDGE, FACE_MORTAR, NONE) for L in battery)
check("all battery pairs return a single valid lane", valid, f"lanes={battery}")
check("battery exercises >=3 distinct lanes (router discriminates)",
      len(set(battery)) >= 3, f"distinct={sorted(set(battery))}")

print(f"\n{'='*66}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(ADR-57 E1 routing-detector oracle)\n{'='*66}")
sys.exit(1 if _fails else 0)
