// ADR-63 P0 build-free self-check for LadrunoContactNormalField.h
// Oracle-first gate for averaged nodal-normal smoothing (no OpenSees):
//   - flood-fill manifold orientation (coherent σ on a MIXED-winding patch)
//   - non-manifold / disconnected refusal
//   - area-weighted nodal normals + GLOBAL outward sign
//   - smoothNormal C0 continuity across a junction (the headline property)
//   - flat-master consistency (n_smooth == facet normal)
//   - degenerate-blend fallback
//
//   g++ -std=c++11 -I ../../SRC/domain/contact proto_nodal_normal_selfcheck.cpp -o nnf_check && ./nnf_check
//
// Sibling of proto_reemit_selfcheck.cpp (ADR-60 P0). See
// Ladruno_implementation/63_ladruno_nodal_normal_smoothing_adr.md.

#include "LadrunoContactNormalField.h"
#include "LadrunoContactProjection.h"   // shape() for the C0-continuity query
#include <cstdio>
#include <cmath>

using namespace LadrunoContactNormalField;

static int fails = 0;
static void check(bool ok, const char *name) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) fails++;
}
static double vlen(const double a[3]) { return std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); }

// build the field the way the engine does: freeze the global sign via voteSign, then nodalNormals(G).
static int buildField(const int *mTags, int nSeg, int nps, const double *X, const int *sigma,
                      const double seed[3], double *snn) {
    double conf;
    double G = voteSign(nSeg, nps, X, sigma, seed, &conf);
    return nodalNormals(mTags, nSeg, nps, X, sigma, G, snn);
}

int main()
{
    const double seedUp[3] = { 0, 0, 1 };   // -outward / "toward the slave" = +z

    // ---------------------------------------------------------------- flat 2-quad strip
    // A=[1,2,3,4] unit square; B=[2,5,6,3] the next square. Shared edge {2,3}. Both wound
    // CCW seen from +z ⇒ coherent (opposite traversal of the shared edge) ⇒ σ both +1.
    {
        int mTags[2*4] = { 1,2,3,4,  2,5,6,3 };
        double X[2*4*3] = {
            0,0,0, 1,0,0, 1,1,0, 0,1,0,     // A
            1,0,0, 2,0,0, 2,1,0, 1,1,0,     // B
        };
        int sigma[2];
        int st = propagateOrientation(mTags, 2, 4, sigma);
        check(st == OK, "flat strip: orientation OK");
        check(sigma[0] == 1 && sigma[1] == 1, "flat strip: coherent winding (σ +1,+1)");

        double snn[2*4*3];
        int nDeg = buildField(mTags, 2, 4, X, sigma, seedUp, snn);
        check(nDeg == 0, "flat strip: no degenerate nodal normals");
        // every nodal normal == +z
        bool allUp = true;
        for (int s = 0; s < 2; s++) for (int k = 0; k < 4; k++) {
            const double *n = snn + (size_t)(s*4+k)*3;
            if (std::fabs(n[0]) > 1e-12 || std::fabs(n[1]) > 1e-12 || std::fabs(n[2]-1) > 1e-12) allUp = false;
        }
        check(allUp, "flat strip: every nodal normal == +z (== facet normal)");
    }

    // ----------------------------------------------------- MIXED winding still coherent
    // Same flat strip but B reverse-wound [2,3,6,5] (CW ⇒ raw facet −z). propagateOrientation
    // must detect the inconsistency (σ_B = −1) and nodalNormals must STILL produce a coherent
    // +z field (σ corrects the reversed winding). This is the load-bearing winding test.
    {
        int mTags[2*4] = { 1,2,3,4,  2,3,6,5 };
        double X[2*4*3] = {
            0,0,0, 1,0,0, 1,1,0, 0,1,0,     // A  (CCW, +z)
            1,0,0, 1,1,0, 2,1,0, 2,0,0,     // B reverse-declared (CW, raw −z)
        };
        int sigma[2];
        int st = propagateOrientation(mTags, 2, 4, sigma);
        check(st == OK, "mixed winding: orientation OK");
        check(sigma[0] == 1 && sigma[1] == -1, "mixed winding: σ flips the reversed quad (+1,−1)");

        double snn[2*4*3];
        buildField(mTags, 2, 4, X, sigma, seedUp, snn);
        bool allUp = true;
        for (int s = 0; s < 2; s++) for (int k = 0; k < 4; k++) {
            const double *n = snn + (size_t)(s*4+k)*3;
            if (std::fabs(n[2]-1) > 1e-12) allUp = false;
        }
        check(allUp, "mixed winding: coherent +z field despite declared mixed winding");
    }

    // -------------------------------------------------- convex ridge (tent) C0 continuity
    // Two roof facets meeting at a ridge along y at (x=0,z=1). Ridge nodes 10,11 shared.
    //   L = [20,10,11,21]   R = [10,30,31,11]   (coherent: opposite traversal of {10,11})
    // The nodal normal at the ridge is +z by symmetry; the field tilts toward ∓x at the eaves.
    // n_smooth at the shared-edge MIDPOINT must agree evaluated from L vs from R (C0 across the junction).
    {
        int mTags[2*4] = { 20,10,11,21,  10,30,31,11 };
        double X[2*4*3] = {
            -1,0,0,  0,0,1,  0,1,1,  -1,1,0,    // L (left roof)
             0,0,1,  1,0,0,  1,1,0,   0,1,1,    // R (right roof)
        };
        int sigma[2];
        int st = propagateOrientation(mTags, 2, 4, sigma);
        check(st == OK && sigma[0] == 1 && sigma[1] == 1, "ridge: coherent tent (σ +1,+1)");

        double snn[2*4*3];
        buildField(mTags, 2, 4, X, sigma, seedUp, snn);

        // ridge nodal normal (L local node 1 = tag10, node 2 = tag11) ≈ +z (symmetric), z>0
        const double *n10 = snn + (size_t)(0*4+1)*3;   // L's node for tag 10
        const double *n11 = snn + (size_t)(0*4+2)*3;   // L's node for tag 11
        check(std::fabs(n10[0]) < 1e-9 && n10[2] > 0.5, "ridge: ridge-node normal ≈ +z (no sign flip)");
        // the field VARIES: L's eave node (tag20, local0) tilts toward −x
        const double *n20 = snn + (size_t)(0*4+0)*3;
        check(n20[0] < -1e-3, "ridge: eave normal tilts (−x) ⇒ field is non-constant (smoothing active)");

        // C0 continuity at the shared-edge midpoint, from L (ξ=+1 edge) and from R (ξ=−1 edge):
        double N[4], dN1[4], dN2[4];
        double nL[3], nR[3];
        // L: tag order [20,10,11,21] ⇒ shared edge 10(local1)-11(local2) = ξ=+1 edge; midpoint η=0
        double nodL[4][3];
        for (int k = 0; k < 4; k++) for (int d = 0; d < 3; d++) nodL[k][d] = snn[(size_t)(0*4+k)*3+d];
        LadrunoContactProjection::shape(4, 1.0, 0.0, N, dN1, dN2);
        bool okL = smoothNormal(4, N, nodL, nL);
        // R: tag order [10,30,31,11] ⇒ shared edge 10(local0)-11(local3) = ξ=−1 edge; midpoint η=0
        double nodR[4][3];
        for (int k = 0; k < 4; k++) for (int d = 0; d < 3; d++) nodR[k][d] = snn[(size_t)(1*4+k)*3+d];
        LadrunoContactProjection::shape(4, -1.0, 0.0, N, dN1, dN2);
        bool okR = smoothNormal(4, N, nodR, nR);
        double diff[3] = { nL[0]-nR[0], nL[1]-nR[1], nL[2]-nR[2] };
        check(okL && okR && vlen(diff) < 1e-12, "ridge: n_smooth C0-continuous across the junction");
    }

    // ----------------------------------------------------- flat-master n_smooth == facet
    // On a flat quad the smoothed normal must equal the facet normal to round-off (D5 consistency).
    {
        int mTags[1*4] = { 1,2,3,4 };
        double X[1*4*3] = { 0,0,0, 2,0,0, 2,3,0, 0,3,0 };   // a flat rectangle
        int sigma[1];
        propagateOrientation(mTags, 1, 4, sigma);
        double snn[1*4*3];
        buildField(mTags, 1, 4, X, sigma, seedUp, snn);
        double N[4], dN1[4], dN2[4]; LadrunoContactProjection::shape(4, 0.3, -0.4, N, dN1, dN2);
        double nod[4][3]; for (int k=0;k<4;k++) for (int d=0;d<3;d++) nod[k][d]=snn[(size_t)k*3+d];
        double ns[3]; smoothNormal(4, N, nod, ns);
        double facet[3] = { 0,0,1 };
        double diff[3] = { ns[0]-facet[0], ns[1]-facet[1], ns[2]-facet[2] };
        check(vlen(diff) < 1e-14, "flat master: n_smooth == facet normal (<1e-14)");
    }

    // ---------------------------------------------------------- global sign from the seed
    // Seeding with −z must flip the whole field (BLOCKER-GLOBAL-SIGN: a global surface datum).
    // High confidence (seed ∥ field) for a +z seed; ~0 confidence for an in-plane seed.
    {
        int mTags[1*4] = { 1,2,3,4 };
        double X[1*4*3] = { 0,0,0, 1,0,0, 1,1,0, 0,1,0 };
        int sigma[1]; propagateOrientation(mTags, 1, 4, sigma);
        double snnUp[1*4*3], snnDn[1*4*3];
        const double seedDn[3] = { 0,0,-1 };
        buildField(mTags, 1, 4, X, sigma, seedUp, snnUp);
        buildField(mTags, 1, 4, X, sigma, seedDn, snnDn);
        check(snnUp[2] > 0.5 && snnDn[2] < -0.5, "global sign: −z seed flips the field");
        double confUp, confIn; const double seedIn[3] = { 1,0,0 };  // in-plane ⇒ ⟂ the +z field
        voteSign(1, 4, X, sigma, seedUp, &confUp);
        voteSign(1, 4, X, sigma, seedIn, &confIn);
        check(confUp > 0.99, "sign confidence: ≈1 when seed ∥ field");
        check(confIn < 0.01, "sign confidence: ≈0 when seed ⟂ field (ill-conditioned ⇒ handler warns)");
    }

    // ----------------------------------------------------------------- non-manifold refuse
    // Three quads sharing one edge {2,3} ⇒ orientation ill-defined ⇒ NON_MANIFOLD.
    {
        int mTags[3*4] = { 1,2,3,4,  2,5,6,3,  2,7,8,3 };
        int sigma[3];
        int st = propagateOrientation(mTags, 3, 4, sigma);
        check(st == NON_MANIFOLD, "non-manifold (>2 seg/edge) ⇒ refuse");
    }

    // ------------------------------------------------------------------ disconnected refuse
    // Two quads sharing NO edge ⇒ a multi-shell the seed can't reach ⇒ DISCONNECTED.
    {
        int mTags[2*4] = { 1,2,3,4,  5,6,7,8 };
        int sigma[2];
        int st = propagateOrientation(mTags, 2, 4, sigma);
        check(st == DISCONNECTED, "disconnected multi-shell ⇒ refuse");
    }

    // -------------------------------------------------------------- closed surface refuse
    // A tetrahedron (4 tri facets, every edge shared by exactly 2 ⇒ NO boundary ⇒ CLOSED). One global
    // outward sign can't serve a slave wrapping around it ⇒ refuse (review F4). Winding coherent.
    {
        // outward-wound faces of tet {1,2,3,4}: [1,3,2],[1,2,4],[1,4,3],[2,3,4]
        int mTags[4*3] = { 1,3,2,  1,2,4,  1,4,3,  2,3,4 };
        int sigma[4];
        int st = propagateOrientation(mTags, 4, 3, sigma);
        check(st == CLOSED, "closed shell (tetrahedron, no boundary edge) ⇒ refuse");
    }

    // -------------------------------------------------------------- tri-3 coherent + query
    // Two triangles sharing edge {2,3}: T0=[1,2,3], T1=[2,4,3] (note 3→2 vs 2→3, opposite ⇒ coherent).
    {
        int mTags[2*3] = { 1,2,3,  2,4,3 };
        double X[2*3*3] = {
            0,0,0, 1,0,0, 0,1,0,     // T0
            1,0,0, 1,1,0, 0,1,0,     // T1
        };
        int sigma[2];
        int st = propagateOrientation(mTags, 2, 3, sigma);
        check(st == OK && sigma[0] == 1 && sigma[1] == 1, "tri-3: coherent pair (σ +1,+1)");
        double snn[2*3*3];
        int nDeg = buildField(mTags, 2, 3, X, sigma, seedUp, snn);
        check(nDeg == 0 && std::fabs(snn[2]-1) < 1e-12, "tri-3: flat field == +z");
    }

    // ------------------------------------------------------------ degenerate-blend fallback
    // Antiparallel corner normals (a folded element): the shape-fn blend at the midpoint cancels
    // ⇒ smoothNormal returns false ⇒ caller falls back to the faceted normal (BLOCKER-FALLBACK (b)).
    {
        double nod[4][3] = { {0,0,1}, {0,0,1}, {0,0,-1}, {0,0,-1} };
        double N[4], dN1[4], dN2[4];
        LadrunoContactProjection::shape(4, 0.0, 0.0, N, dN1, dN2);   // center: equal weights ⇒ Σ ≈ 0
        double n[3];
        bool ok = smoothNormal(4, N, nod, n);
        check(!ok, "degenerate blend (antiparallel corners) ⇒ smoothNormal false ⇒ fallback");
    }

    std::printf("\n%s — %d failure(s)\n", fails == 0 ? "ALL PASS" : "FAILURES", fails);
    return fails == 0 ? 0 : 1;
}
