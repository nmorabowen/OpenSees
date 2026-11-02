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

// Ladruno WP-221: the state of a half-storage SOE's asymmetry guard, made
// readable by a script.
//
// `system Pardiso -matrixType 1|2` and `system ProfileSPD` store only the upper
// triangle, so a genuinely unsymmetric element tangent has its lower half
// DISCARDED and the run solves the reflected system (ADR-75 P1d). Both SOEs
// already sample the element matrices in addA() and print a one-time
// `WARNING ... UNSYMMETRIC ...` (PARDISOGenLinSOE.cpp / ProfileSPDLinSOE.cpp).
// That warning was the ONLY output: a caller that wanted to fall back to full
// storage had to scrape stderr for the word, and any rewording of the message
// would have switched its fallback off without a sound. This struct is what the
// guard knows, kept where `ladrunoSymmetryGuard` (LadrunoSymmetryGuardQuery.h)
// can hand it back, plus the `-strictSymmetry` switch that turns the latch into
// a refused solve.
//
// Header-only and dependency-free on purpose: it is a member of two upstream
// SOE classes, so it must not pull anything into their headers but Matrix.

#ifndef LadrunoSymmetryGuard_h
#define LadrunoSymmetryGuard_h

#include <Matrix.h>

struct LadrunoSymmetryGuard
{
    // Relative tolerance of the mirror compare. Unchanged from ADR-75 P1d.
    static double tolerance() { return 1.0e-8; }

    int strict;         // -strictSymmetry: check EVERY assembly (no sampling
                        // budget) and refuse to solve once latched
    int latched;        // an element matrix exceeded tolerance(); once set,
                        // stays set for the life of the SOE
    int pass;           // tangent assembly (of the current pattern) where it
                        // latched; -1 while unlatched
    double worstDev;    // max |m(i,j) - m(j,i)| of the latching matrix
    double scale;       // max |m(i,j)| of that matrix
    double maxRel;      // worst worstDev/scale seen by ANY check, latched or not
    int checks;         // element matrices compared (0 = the guard never looked)
    int checkedPasses;  // tangent assemblies in which at least one check ran
    int lastCheckedPass;
    int refused;        // solves refused under -strictSymmetry

    LadrunoSymmetryGuard()
        : strict(0), latched(0), pass(-1), worstDev(0.0), scale(0.0),
          maxRel(0.0), checks(0), checkedPasses(0), lastCheckedPass(-1),
          refused(0) {}

    // setSize() restarts the assembly count of the new pattern; the latch and
    // the totals deliberately survive (reported once per SOE, as before).
    void newPattern() { lastCheckedPass = -1; }

    // Compare the strict lower triangle of an n x n element matrix against its
    // mirror. Returns true the first time the guard latches, so the caller
    // prints its warning exactly once.
    bool check(const Matrix &m, int n, int atPass)
    {
        double dev = 0.0, sc = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < i; j++) {
                double u = m(i, j), v = m(j, i);
                double au = u < 0.0 ? -u : u;
                double av = v < 0.0 ? -v : v;
                if (au > sc) sc = au;
                if (av > sc) sc = av;
                double d = u - v;
                if (d < 0.0) d = -d;
                if (d > dev) dev = d;
            }
        checks++;
        if (atPass != lastCheckedPass) {
            checkedPasses++;
            lastCheckedPass = atPass;
        }
        if (sc > 0.0 && dev / sc > maxRel)
            maxRel = dev / sc;
        if (latched == 0 && sc > 0.0 && dev > tolerance() * sc) {
            latched = 1;
            pass = atPass;
            worstDev = dev;
            scale = sc;
            return true;
        }
        return false;
    }
};

#endif
