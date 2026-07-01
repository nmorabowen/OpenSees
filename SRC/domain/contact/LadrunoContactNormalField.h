// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// ADR-63 P0 — LadrunoContactNormalField: header-only, OpenSees-free averaged
// nodal-normal smoothing for the NTS contact lane (ADR-47 #4a; Abaqus TG §5.1.1's
// "unit normals computed at every master node, averaging adjacent segment normals
// … a smoothly varying normal field N(X)"). NO OpenSees types — raw int/double[] +
// STL only, so the engine (LadrunoContactDomain) AND a build-free oracle
// (contact_prototypes/proto_nodal_normal_selfcheck.cpp) share ONE logic. Mirrors
// LadrunoContactProjection.h / LadrunoContactReemit.h discipline.
//
// WHY (resolves ADR-60 R3 + the NTS slice of ADR-41 Q-NORMAL): the shipped per-
// segment FACETED normal (LadrunoContactProjection::normalOriented) has a sign
// chosen PER (slave,segment) from the slave's reference side, so a slave sliding
// over a convex ridge flips the sign → silent pass-through (R3); and the normal
// JUMPS across segment junctions → chatter (Q-NORMAL). The smooth field fixes BOTH:
// a node's normal is shared by its incident segments (C0-continuous across
// junctions — NOT C1; that is the deferred Hermite #4), and the outward sense is a
// GLOBAL surface datum (a coherent winding + an -outward / slave-centroid seed),
// captured once, never per-slave → the ridge flip is structurally impossible.
//
// THREE pure stages (see Ladruno_implementation/63_ladruno_nodal_normal_smoothing_adr.md):
//   propagateOrientation() — flood-fill manifold orientation: a per-segment winding
//     sign σ_s giving globally-coherent winding from a seed. TOPOLOGICAL (depends only
//     on connectivity, not coords) ⇒ the engine caches it per R1 membership fingerprint.
//     REFUSES (status) a non-manifold edge (>2 segments), a non-orientable surface
//     (Möbius sign conflict), or a disconnected multi-shell (BLOCKER-WINDING).
//   nodalNormals() — per node, the area-weighted (Newell) sum of the consistently-wound
//     adjacent facet normals, with the GLOBAL outward sign (BLOCKER-GLOBAL-SIGN), scattered
//     back to a per-segment-per-node layout for adapter consumption.
//   smoothNormal() — interpolate the segment's nodal normals with the shape fns N[] and
//     renormalize; returns false on a degenerate query-point blend (BLOCKER-FALLBACK (b)).
//
// NO tangent math lives here: the MVP ships the symmetric frozen-field kn·BᵀB tangent
// (ADR-63 §Gate decision, RESOLVED symmetric-first); the full ∂n_smooth/∂u is the gated P3.

#ifndef LadrunoContactNormalField_h
#define LadrunoContactNormalField_h

#include <cmath>
#include <vector>
#include <map>

namespace LadrunoContactNormalField {

// ------------------------------------------------------------------ small vec3
inline double nf_dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline double nf_norm(const double a[3]) { return std::sqrt(nf_dot(a, a)); }

// ----------------------------------------------------------------- status codes
// Orientation outcome. Anything but OK ⇒ the engine refuses smoothing on this surface
// and the handler falls back to the faceted normalOriented() + -outward path (BLOCKER-WINDING).
enum Status {
    OK             = 0,
    NON_MANIFOLD   = 1,   // an edge is shared by >2 segments (refuse: orientation ill-defined)
    DISCONNECTED   = 2,   // not all segments reachable from the seed (multi-shell — refuse w/o -outward)
    NON_ORIENTABLE = 3,   // a sign conflict on a cycle (Möbius / folded — refuse)
    CLOSED         = 4    // no boundary edges (a closed shell: box/tube) — one global sign can't serve
                          // both sides (review F4); refuse, the supported target is an OPEN patch
};

// Area-weighted facet normal of one segment via Newell's method: the result's magnitude is
// 2·area and its direction is the face normal in the node-WINDING order (handles a non-planar
// quad; the linear case is exact). P = nps×3 row-major. Winding-order-dependent BY DESIGN — the
// σ_s sign from propagateOrientation() makes the per-segment windings coherent.
inline void newellAreaNormal(int nps, const double *P, double out[3]) {
    out[0] = out[1] = out[2] = 0.0;
    for (int i = 0; i < nps; i++) {
        const double *a = P + (size_t)i * 3;
        const double *b = P + (size_t)((i + 1) % nps) * 3;
        out[0] += (a[1] - b[1]) * (a[2] + b[2]);
        out[1] += (a[2] - b[2]) * (a[0] + b[0]);
        out[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    for (int d = 0; d < 3; d++) out[d] *= 0.5;   // Newell sum is 2·area·n̂; halve to area·n̂
}

// ------------------------------------------------ flood-fill manifold orientation
// Fill sigma[nSeg] (each +1/−1) with a globally-coherent winding from a seed segment, by
// propagating across shared EDGES: in a coherently oriented manifold two adjacent faces
// traverse their shared edge in OPPOSITE node order. mTags = flat nSeg*nps connectivity
// (segment s spans mTags[s*nps .. +nps), CCW/CW as declared). nps = 3 (tri) or 4 (quad).
// TOPOLOGICAL: depends only on mTags, so cache the result per membership fingerprint.
// Returns a Status; sigma is only meaningful on OK. The seed segment is index 0, sigma=+1
// (the GLOBAL outward sign is applied later by nodalNormals(), not here).
inline int propagateOrientation(const int *mTags, int nSeg, int nps, int *sigma) {
    if (mTags == 0 || sigma == 0 || nSeg <= 0 || (nps != 3 && nps != 4))
        return NON_ORIENTABLE;   // degenerate input — refuse (caller falls back to faceted)

    // Undirected edge {lo,hi} -> list of (segment, dir) where dir=+1 if this segment traverses
    // lo→hi in winding order, −1 if hi→lo. >2 entries on any edge ⇒ non-manifold.
    struct EdgeUse { int seg; int dir; };
    std::map<std::pair<int,int>, std::vector<EdgeUse> > edges;
    for (int s = 0; s < nSeg; s++)
        for (int k = 0; k < nps; k++) {
            int a = mTags[s*nps + k];
            int b = mTags[s*nps + (k + 1) % nps];
            if (a == b) return NON_ORIENTABLE;            // a collapsed edge — refuse
            int lo = a < b ? a : b, hi = a < b ? b : a;
            EdgeUse u; u.seg = s; u.dir = (a < b) ? +1 : -1;
            edges[std::make_pair(lo, hi)].push_back(u);
        }
    // Adjacency: for each interior (2-segment) edge, the relative sign that makes the pair
    // coherent. >2 ⇒ non-manifold. (1-segment edges are boundary edges of an open patch — fine.)
    struct Adj { int nbr; int relSign; };
    std::vector<std::vector<Adj> > adj((size_t)nSeg);
    int nBoundary = 0;                         // edges with exactly 1 use ⇒ open-patch boundary
    for (std::map<std::pair<int,int>, std::vector<EdgeUse> >::const_iterator it = edges.begin();
         it != edges.end(); ++it) {
        const std::vector<EdgeUse> &use = it->second;
        if (use.size() > 2) return NON_MANIFOLD;
        if (use.size() == 1) nBoundary++;
        if (use.size() == 2) {
            // coherent ⇒ opposite traversal (dir1 != dir2). relSign = −(dir1·dir2):
            //   same direction  (dir1==dir2) ⇒ relSign −1 (neighbor must flip),
            //   opposite        (dir1!=dir2) ⇒ relSign +1 (already coherent).
            int rel = -(use[0].dir * use[1].dir);
            Adj a0; a0.nbr = use[1].seg; a0.relSign = rel; adj[use[0].seg].push_back(a0);
            Adj a1; a1.nbr = use[0].seg; a1.relSign = rel; adj[use[1].seg].push_back(a1);
        }
    }
    // BFS from seg 0; σ[neighbor] = σ[cur]·relSign. A conflict on an already-visited segment
    // ⇒ non-orientable. Unvisited segments after BFS ⇒ disconnected (multi-shell).
    for (int s = 0; s < nSeg; s++) sigma[s] = 0;   // 0 = unvisited
    std::vector<int> stack;
    sigma[0] = +1; stack.push_back(0);
    int visited = 1;
    while (!stack.empty()) {
        int cur = stack.back(); stack.pop_back();
        for (size_t e = 0; e < adj[cur].size(); e++) {
            int nbr = adj[cur][e].nbr;
            int want = sigma[cur] * adj[cur][e].relSign;
            if (sigma[nbr] == 0) { sigma[nbr] = want; visited++; stack.push_back(nbr); }
            else if (sigma[nbr] != want) return NON_ORIENTABLE;   // cycle sign conflict
        }
    }
    if (visited != nSeg) return DISCONNECTED;   // a second shell the seed can't reach
    // A coherently-wound, connected, manifold surface with NO boundary edges is CLOSED (a box/tube):
    // one global outward sign cannot serve a slave that wraps around it (review F4). Refuse — the
    // supported target is an OPEN patch (which has boundary edges). nSeg==1 (a single facet) is open.
    if (nBoundary == 0 && nSeg > 1) return CLOSED;
    return OK;
}

// ------------------------------------------------------- area-weighted nodal normals
// Build the smooth nodal-normal field from CURRENT coords + the cached σ_s, with the GLOBAL
// outward sign, and SCATTER it to a per-segment-per-node layout segNodalNorm[nSeg*nps*3]
// (row-major [seg][node][xyz]) so each adapter consumes its own segment's 4 normals directly.
//
//   m_a (unnormalized) = Σ_{s ∋ a} σ_s · newellAreaNormal(s)         [area-weighted, coherent]
//   n_a                = G · normalize(m_a)
//
// G (the GLOBAL outward sign) is decided ONCE by the engine via voteSign() at the first handle and
// passed in FROZEN here (ADR-63 D2 / review F1: NEVER re-decided from a deformed config each handle —
// that could flip the whole field mid-run). G ∈ {+1,−1}. segCoords = flat nSeg*nps*3 CURRENT coords
// (the magnitudes/directions of n_a track deformation; only the sign is frozen). A node whose
// |m_a| → 0 RELATIVE to its largest incident facet contribution (a sharp fold where coherent facets
// cancel) is left ZERO ⇒ a query touching it falls back via smoothNormal()'s degeneracy gate (the
// relative test catches near-cancellation that the old absolute 1e-300 floor normalized to noise —
// review GAP-5). Returns the count of degenerate (zero) unique nodal normals.
inline int nodalNormals(const int *mTags, int nSeg, int nps, const double *segCoords,
                        const int *sigma, double G, double *segNodalNorm) {
    // accumulate per UNIQUE node tag, tracking the largest single-facet contribution per node
    std::map<int, std::vector<double> > acc;   // tag -> [x,y,z]
    std::map<int, double> accMax;              // tag -> max |σ_s·newell_s| over incident segments
    for (int s = 0; s < nSeg; s++) {
        double fn[3];
        newellAreaNormal(nps, segCoords + (size_t)s*nps*3, fn);
        double sgn = (double)sigma[s];
        double c[3] = { sgn*fn[0], sgn*fn[1], sgn*fn[2] };   // coherently-oriented area normal
        double cmag = nf_norm(c);
        for (int k = 0; k < nps; k++) {
            int tag = mTags[s*nps + k];
            std::vector<double> &v = acc[tag];
            if (v.empty()) v.assign(3, 0.0);
            for (int d = 0; d < 3; d++) v[d] += c[d];
            if (cmag > accMax[tag]) accMax[tag] = cmag;
        }
    }
    // normalize each unique node normal (apply the frozen G), tally degeneracies
    int nDeg = 0;
    std::map<int, std::vector<double> > unit;   // tag -> unit n_a (G applied); zero if degenerate
    for (std::map<int, std::vector<double> >::iterator it = acc.begin(); it != acc.end(); ++it) {
        double m[3] = { it->second[0], it->second[1], it->second[2] };
        double mn = nf_norm(m);
        std::vector<double> u(3, 0.0);
        // RELATIVE degeneracy: a coherent node sum should be O(its incident facet areas). If it has
        // collapsed to ≪ the largest single incident contribution, the facets cancel (a cusp) ⇒ leave
        // it ZERO so the query falls back, rather than normalizing numerical noise to a unit vector.
        double ref = accMax[it->first];
        if (mn > 1e-9 * (ref + 1e-300)) {
            for (int d = 0; d < 3; d++) u[d] = G * m[d] / mn;
        } else {
            nDeg++;   // leave zero ⇒ queries touching this node fall back to faceted (BLOCKER-FALLBACK)
        }
        unit[it->first] = u;
    }
    // scatter to the per-segment-per-node layout
    for (int s = 0; s < nSeg; s++)
        for (int k = 0; k < nps; k++) {
            const std::vector<double> &u = unit[mTags[s*nps + k]];
            double *dst = segNodalNorm + ((size_t)s*nps + k) * 3;
            for (int d = 0; d < 3; d++) dst[d] = u[d];
        }
    return nDeg;
}

// The coherent area-weighted vote vector Σ_s σ_s·newellAreaNormal(s) — the surface's aggregate
// oriented normal. The GLOBAL sign is sign(vote·seed); its CONFIDENCE is |cos∠(vote,seed)| ∈ [0,1]
// (≈0 ⇒ the seed is ~⟂ the field ⇒ the auto sign is a near coin-flip: a two-sided / saddle / straddle
// master — review F2/F3/F5). The engine calls this ONCE (first handle) to freeze the sign + warn if
// the confidence is low. Returns the chosen sign (+1/−1); fills *conf if non-null.
inline double voteSign(int nSeg, int nps, const double *segCoords, const int *sigma,
                       const double seed[3], double *conf) {
    double vote[3] = {0.0, 0.0, 0.0};
    for (int s = 0; s < nSeg; s++) {
        double fn[3];
        newellAreaNormal(nps, segCoords + (size_t)s*nps*3, fn);
        for (int d = 0; d < 3; d++) vote[d] += (double)sigma[s] * fn[d];
    }
    double vn = nf_norm(vote), sn = nf_norm(seed), vs = nf_dot(vote, seed);
    if (conf != 0) *conf = (vn > 1e-300 && sn > 1e-300) ? (vs >= 0 ? vs : -vs) / (vn * sn) : 0.0;
    return (vs < 0.0) ? -1.0 : +1.0;
}

// ------------------------------------------------------------- smooth normal query
// n = normalize( Σ_i N[i] · nodalNorm[i] ) over a segment's nps nodes. N[] = the master shape
// functions at the projection (LadrunoContactProjection::shape fills them); nodalNorm = the
// segment's nps nodal normals (the scattered row from nodalNormals()). Returns false on a
// degenerate blend (|Σ|→0 — antiparallel corner normals across a folded element, OR a node left
// zero by a degenerate nodal normal): the caller then falls back to the faceted normalOriented()
// for this query (BLOCKER-FALLBACK (b)).
inline bool smoothNormal(int nps, const double N[4], const double nodalNorm[4][3], double n[3]) {
    double s[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < nps; i++)
        for (int d = 0; d < 3; d++) s[d] += N[i] * nodalNorm[i][d];
    double sn = nf_norm(s);
    if (sn < 1e-12) return false;            // degenerate blend ⇒ fall back to faceted
    for (int d = 0; d < 3; d++) n[d] = s[d] / sn;
    return true;
}

} // namespace LadrunoContactNormalField

#endif
