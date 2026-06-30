// ADR-60 P0 build-free self-check for LadrunoContactReemit.h
// Oracle-first gate for the migration metric + hysteresis trigger (no OpenSees).
//
//   g++ -std=c++11 -I ../../SRC/domain/contact proto_reemit_selfcheck.cpp -o reemit_check && ./reemit_check
//
// Sibling of proto_bucket_sort.py (the ADR-39 P2.5 broad-phase oracle).

#include "LadrunoContactReemit.h"
#include <cstdio>
#include <cmath>

using namespace LadrunoContactReemit;

static int fails = 0;
static void check(bool ok, const char *name) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) fails++;
}

int main()
{
    // --- referenceBand: median segment diagonal * cellFrac ---
    // 3 quad segments (nps=4): diagonals sqrt(2), sqrt(2), 2*sqrt(2); median = sqrt(2).
    double segs[3 * 4 * 3] = {
        0,0,0,  1,0,0,  1,1,0,  0,1,0,    // unit square -> diag sqrt(2)
        2,0,0,  3,0,0,  3,1,0,  2,1,0,    // unit square -> diag sqrt(2)
        0,2,0,  2,2,0,  2,4,0,  0,4,0,    // side-2 square -> diag 2*sqrt(2)
    };
    double band = referenceBand(3, 4, segs, 1.0);
    check(std::fabs(band - std::sqrt(2.0)) < 1e-12, "referenceBand = median diag (sqrt2)");
    check(std::fabs(referenceBand(3, 4, segs, 0.5) - 0.5 * std::sqrt(2.0)) < 1e-12,
          "referenceBand scales with cellFrac");
    check(referenceBand(0, 4, segs, 1.0) == 0.0, "referenceBand degenerate -> 0");

    // --- maxMigration ---
    double anchors[2 * 3] = { 0,0,0,  5,5,5 };
    double cur[2 * 3]     = { 0.3,0,0,  5,5,5 };   // slave 0 drifted 0.3, slave 1 still
    check(std::fabs(maxMigration(2, anchors, cur) - 0.3) < 1e-12, "maxMigration picks max drift");
    double cur2[2 * 3]    = { 0,0,0,  5,5,8 };     // slave 1 drifted 3
    check(std::fabs(maxMigration(2, anchors, cur2) - 3.0) < 1e-12, "maxMigration second slave");

    // --- Trigger: hysteresis + floor (f=0.5, band=2 -> thresh=1.0, reArm 0.5, floor=3) ---
    {
        Trigger t(0.5, 0.5, 3);
        const double B = 2.0;
        check(t.update(1.5, B) == false, "no fire before floor (step1)");
        check(t.update(1.5, B) == false, "no fire before floor (step2)");
        check(t.update(1.5, B) == true,  "fire at floor when over thresh");
        check(t.update(9.9, B) == false, "disarmed: no re-fire while high");
        check(t.update(9.9, B) == false, "disarmed: still no re-fire");
        check(t.update(0.4, B) == false, "re-arm on fall-back (no fire)");
        check(t.update(1.5, B) == true,  "re-fire after re-arm + floor");
    }

    // --- Trigger: parked exactly at threshold must NOT chatter every step (gate Lens-A A5) ---
    {
        Trigger t(0.5, 0.5, 1);     // floor=1 (worst case)
        const double B = 2.0;       // thresh = 1.0
        int fires = 0;
        for (int s = 0; s < 20; s++)
            if (t.update(1.0001, B)) fires++;   // parked just over thresh, never falls back
        check(fires == 1, "parked-at-threshold fires once, not every step (hysteresis)");
    }

    // --- Trigger: degenerate band never fires ---
    {
        Trigger t(0.5, 0.5, 1);
        check(t.update(1e9, 0.0) == false, "degenerate band never fires");
    }

    // --- membershipFingerprint: order-sensitive mTags hash (R1 / BLOCKER-MEMBERSHIP) ---
    {
        int a[8] = { 10,11,21,20,  11,12,22,21 };   // 2 quads, nps=4
        int same[8] = { 10,11,21,20,  11,12,22,21 };
        int perm[8] = { 11,12,22,21,  10,11,21,20 }; // SAME tags, swapped segment order
        int chg[8]  = { 10,11,21,20,  11,12,23,21 }; // one tag differs (22 -> 23)
        unsigned long long f0 = membershipFingerprint(a, 8, 4);
        check(membershipFingerprint(same, 8, 4) == f0, "fingerprint stable under identical ordering");
        check(membershipFingerprint(perm, 8, 4) != f0, "fingerprint changes on segment reorder (permutation)");
        check(membershipFingerprint(chg,  8, 4) != f0, "fingerprint changes on a single tag change");
        check(membershipFingerprint(a, 8, 3) != f0, "fingerprint changes on nps change");
        check(membershipFingerprint(a, 4, 4) != f0, "fingerprint changes on count change");
        check(membershipFingerprint(0, 0, 4) == membershipFingerprint(0, 0, 4), "empty fingerprint well-defined");
    }

    std::printf(fails == 0 ? "\nALL PASS\n" : "\n%d FAIL\n", fails);
    return fails == 0 ? 0 : 1;
}
