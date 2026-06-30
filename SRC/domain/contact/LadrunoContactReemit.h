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

// ADR-60 P0 — LadrunoContactReemit: header-only, OpenSees-free finite-sliding
// re-emission policy for the NTS bucket-sort broad phase. NO OpenSees types —
// raw double[] + STL only, so the engine (LadrunoContactDomain) and a build-free
// oracle (contact_prototypes/proto_reemit.py / a standalone main) share ONE logic.
//
// The NTS broad phase (LadrunoContactBucketSort) builds its grid ONCE per handle()
// from REFERENCE coords, so a slave whose DEFORMED position slides past the search
// band lands on a non-candidate segment -> no adapter -> silent pass-through
// (ADR-60 §Why). This header decides WHEN to re-emit (rebuild the candidate set
// from the committed config) without any geometry/assembly knowledge:
//
//   * referenceBand()  — the search tolerance band = cellFrac * median segment
//     diagonal, computed from REFERENCE coords so it is DEFORMATION-INVARIANT and
//     cannot drift away from the grid the next re-emit builds (gate Lens-D D6).
//   * maxMigration()   — max over slaves of |x_committed - anchor|, the drift since
//     the last sort.
//   * Trigger          — hysteresis + min-steps-floor state machine that turns a
//     drift into a "re-sort now" decision (gate Lens-A A5: a bare threshold with a
//     floor of 1 re-triggers every step for a slave parked at f*band).
//
// NO slip-transfer math lives here: ADR-60 D4 re-engages friction on a fresh slot
// at a crossing (already traction-continuous) rather than rotating committed slip
// (Reviewer-C HIGH). Sustained-drag continuity is P4's patch adapter, not this file.
//
// See Ladruno_implementation/60_ladruno_finite_sliding_reemission_adr.md.

#ifndef LadrunoContactReemit_h
#define LadrunoContactReemit_h

#include <cmath>
#include <vector>
#include <algorithm>

namespace LadrunoContactReemit {

// Search-tolerance band = cellFrac * median segment diagonal, from REFERENCE coords
// (segRefCoords: flat nSeg*nps*3, row-major [seg][node][xyz] — the SAME layout the
// bucket-sort Grid consumes). Mirrors LadrunoContactBucketSort's medDiag so the band
// matches the grid cell size, but is keyed to the undeformed mesh so it stays fixed
// across re-emits while the grid itself moves to deformed coords. Returns 0 on a
// degenerate input (caller then never triggers — safe).
inline double referenceBand(int nSeg, int nps, const double *segRefCoords, double cellFrac)
{
    if (nSeg <= 0 || nps <= 0 || segRefCoords == 0) return 0.0;
    std::vector<double> diags(nSeg);
    for (int s = 0; s < nSeg; s++) {
        const double *P = segRefCoords + (size_t)s * nps * 3;
        double dmax = 0.0;
        for (int i = 0; i < nps; i++)
            for (int j = i + 1; j < nps; j++) {
                double dx = P[i*3+0] - P[j*3+0];
                double dy = P[i*3+1] - P[j*3+1];
                double dz = P[i*3+2] - P[j*3+2];
                double d = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (d > dmax) dmax = d;
            }
        diags[s] = dmax;
    }
    std::sort(diags.begin(), diags.end());
    return cellFrac * diags[nSeg / 2];
}

// Master-surface membership fingerprint (R1 / gate BLOCKER-MEMBERSHIP). A 64-bit FNV-1a hash of
// the ORDERED master node-tag sequence (+ nps + count). The NTS friction state is keyed by the
// GLOBAL segment ordinal `seg` (an index INTO this mTags vector), so if the master surface's node
// ordering changes between handles (element removal / re-mesh), ordinal N now names a DIFFERENT
// physical segment and a committed (contactTag, slaveTag, seg) friction slot would ALIAS. The
// engine fingerprints mTags each handle and, on a change, drops the contact's friction slots +
// refuses to arm re-emit that step (named error). FNV-1a mixes each tag with the running hash in
// sequence, so it is ORDER-sensitive: a pure permutation of the same tags changes the result.
// (mTags == 0 / n == 0 is well-defined — the empty-surface fingerprint.)
inline unsigned long long membershipFingerprint(const int *mTags, int n, int nps)
{
    const unsigned long long P = 1099511628211ULL;          // FNV-1a 64 prime
    unsigned long long h = 1469598103934665603ULL;          // FNV-1a 64 offset basis
    h = (h ^ (unsigned long long)(unsigned int)nps) * P;
    h = (h ^ (unsigned long long)(unsigned int)n)   * P;
    if (mTags != 0)
        for (int i = 0; i < n; i++)
            h = (h ^ (unsigned long long)(unsigned int)mTags[i]) * P;
    return h;
}

// Max migration: max over slaves of |x_current - anchor|. anchors/current are flat
// nSlave*3 (row-major [slave][xyz]); anchors = the coords fed to the grid at the last
// handle(), current = the committed coords now. O(nSlave), allocation-free.
inline double maxMigration(int nSlave, const double *anchors, const double *current)
{
    if (nSlave <= 0 || anchors == 0 || current == 0) return 0.0;
    double dmax = 0.0;
    for (int i = 0; i < nSlave; i++) {
        const double *a = anchors + (size_t)i * 3;
        const double *x = current + (size_t)i * 3;
        double dx = x[0] - a[0], dy = x[1] - a[1], dz = x[2] - a[2];
        double d = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (d > dmax) dmax = d;
    }
    return dmax;
}

// Hysteresis + min-steps re-sort trigger (one instance per contact, lives on the
// engine; ticked once per Domain::commit()). update(dmax, band) returns true exactly
// when a re-emit should fire THIS commit. Two independent fire paths:
//   * MIGRATION (the metric): fires when drift dmax > f*band, then disarms until the
//     drift falls back below reArmFrac*f*band (so a slave parked at the threshold cannot
//     re-fire every step — gate Lens-A A5) and enforces >= floorSteps commits between
//     fires (the migration floor; default 10). band <= 0 (degenerate) never fires here.
//   * FORCED CADENCE (R7 / D1 = LS-DYNA BSORT): when forceEvery > 0, fires unconditionally
//     every forceEvery commits regardless of drift/arming (the cycle override `-resortEvery`).
//     A forced fire rebuilds anchors (drift resets), so it also re-arms the migration path.
//     forceEvery == 0 (default) ⇒ no forced fires (pure migration trigger).
//
// R4: this struct PERSISTS across handle()s (the engine keys one per contactTag and refreshes
// only its CONFIG each handle via setConfig — NOT a fresh construct), so sinceLast/sinceForce/
// armed actually accumulate. Previously the engine rebuilt it every handle, making the floor +
// hysteresis vestigial (the anchor rebuild was the de-facto rate-limiter).
//
// IMPORTANT (gate Lens-A A1 / BLOCKER-GA1): the CALLER must NOT invoke update() during
// a held-load ALM augmentation sweep (Domain::isContactAugmenting()) — both the metric
// and the step counters must stand still there, else an augmentation re-commit trips the
// counter and re-handles mid-Uzawa. This struct does not know about augmentation; the
// engine gates the call.
struct Trigger {
    double f;          // drift fraction of band that arms a fire (default 0.5)
    double reArmFrac;  // re-arm when dmax < reArmFrac*f*band (default 0.5 -> 0.25*band)
    int    floorSteps; // min commits between MIGRATION fires (default 10 ~ LS-DYNA cadence)
    int    forceEvery; // R7: forced re-sort cadence in commits (LS-DYNA BSORT); 0 ⇒ off
    bool   armed;      // ready to fire (disarms on a fire, re-arms on fall-back)
    int    sinceLast;  // commits since the last MIGRATION fire
    int    sinceForce; // commits since the last FORCED fire (R7)

    Trigger(double f_ = 0.5, double reArm_ = 0.5, int floor_ = 10, int forceEvery_ = 0)
      : f(f_), reArmFrac(reArm_), floorSteps(floor_), forceEvery(forceEvery_),
        armed(true), sinceLast(0), sinceForce(0) {}

    // R4: refresh CONFIG only (keep armed/sinceLast/sinceForce) — called each handle() so a
    // re-parsed -resortFrac/-resortEvery takes effect without resetting the accumulated cadence.
    void setConfig(double f_, double reArm_, int floor_, int forceEvery_) {
        f = f_; reArmFrac = reArm_; floorSteps = floor_; forceEvery = forceEvery_;
    }

    bool update(double dmax, double band)
    {
        // R7: forced cadence (BSORT) — independent of drift/band/arming. A forced re-sort
        // rebuilds anchors, so reset the migration floor + re-arm so the two paths stay coherent.
        if (forceEvery > 0) {
            if (sinceForce < forceEvery) sinceForce++;  // saturate (no overflow on long runs)
            if (sinceForce >= forceEvery) {
                sinceForce = 0;
                sinceLast = 0;
                armed = true;
                return true;
            }
        }
        if (band <= 0.0) return false;
        if (sinceLast < floorSteps) sinceLast++;       // saturate (no overflow on long runs)
        const double thresh = f * band;
        if (!armed) {                                  // disarmed: wait for fall-back
            if (dmax < reArmFrac * thresh) armed = true;
            return false;
        }
        if (dmax > thresh && sinceLast >= floorSteps) {
            armed = false;
            sinceLast = 0;
            return true;
        }
        return false;
    }
};

} // namespace LadrunoContactReemit

#endif
