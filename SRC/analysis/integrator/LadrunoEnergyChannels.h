/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

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

#ifndef LadrunoEnergyChannels_h
#define LadrunoEnergyChannels_h

#include <limits>          // Ladruno — addEnergy finiteness guard
#include <OPS_Globals.h>   // Ladruno — opserr (process-once guard warning)

// Ladruno ADR-69 — named energy channels for the EnergyBalanceRecorder v2.
// See Ladruno_implementation/69_ladruno_energy_recorder_channels_adr.md.
//
// THE PROBLEM. The v1 recorder accounts four channels (KE, IE, DW, ULW) and
// folds everything else into the closure residual RES. Two classes of energy
// escape it:
//   * CLOSURE gaps — energy applied inside the integrator solve (LNVD/FLAC
//     local damping, modal damping) reaches no element/node quantity the
//     recorder reads, so it pollutes RES.
//   * ATTRIBUTION gaps — energy that IS read but lands in the wrong column.
//     The absorbing-boundary incident injection is the driving case: a
//     LysmerTriangle's stage-0 getResistingForce() returns the stale
//     internalForces member (= C*v_gnd, the incident forcing), so the
//     recorder's IE integral silently absorbs the injection work with
//     resisting-force sign. RES may then ACCIDENTALLY CLOSE while IE lies.
//
// THE SEAM. Same seam as the consistent-SMS registry
// (LadrunoMassScalingEnergy.h): the recorder holds only a Domain* and has no
// handle on the integrator or on element internals, so producers PUBLISH here
// and the recorder QUERIES at record time.
//
// SEMANTICS (differs from the SMS registry — read this):
//   * Channels store CUMULATIVE ENERGY (units of work), not rates and not
//     per-step increments. Producers integrate at their OWN cadence (per
//     commit for elements, per sub-step for integrators) and add increments
//     here. Rationale (ADR-69 F3): the recorder may sample every N steps —
//     trapezoid-integrating a published *rate* at recorder cadence aliases;
//     a running total read at any cadence does not. Sub-stepped producers
//     (Noh-Bathe LNVD fires twice per step) have no single per-step rate.
//   * Totals are MONOTONE over a MODEL lifetime and reset at wipe()
//     (Domain::clearAll — ADR-72 P3 CI find; originally process-lifetime,
//     which poisoned later models after a diverged run and absorbed small
//     increments after a huge-but-finite total; see resetOnWipe()). There is
//     deliberately no per-step clear. The recorder snapshots each channel at
//     its first record() and reports the DELTA since — so a recorder created
//     mid-run still starts from zero.
//   * declare(c) marks a channel as relevant WITHOUT adding energy (a Lysmer
//     element declares at setDomain even if no velocity load ever arrives).
//     The recorder emits a column for every declared channel; undeclared
//     channels cost nothing and print nothing, keeping the no-producer path
//     identical to v1.
//
// CHANNEL CONVENTIONS:
//   * ABSORB_LEAK stores L = sum over absorbing elements of the cumulative
//     RESISTING-SIGN injection work int R_inj^T v dt — i.e. exactly the
//     quantity that pollutes the recorder's raw IE integral. It is published
//     sign-faithful (no convention guess); the recorder derives both
//     corrections from it:  IE_display = IE_raw - L,  E_inject = -L.
//     Rebucketing L cannot change RES (it cancels); it makes IE and the new
//     E_inject column truthful.
//   * LNVD_WORK stores the cumulative FLAC local-damping dissipation
//     sum alpha*|r_i|*|v_i|*dt >= 0 (a true sink; reduces |RES|).
//   * MODAL_WORK stores the cumulative modal-damping dissipation
//     v'C_modal v integrated per commit (a true sink; reduces |RES|).
//     Published by IncrementalIntegrator::commit() (ADR-69 P2) - covers
//     integrators on the base commit path (Newmark); commit()-overriding
//     integrators (HHT family, *_TP explicit) do not publish.
//
// Single-threaded like the rest of the analysis core; header-only so the
// energy kernel and producers include it without CMake changes.

namespace Ladruno {

class EnergyChannelRegistry {
public:
    enum Channel {
        ABSORB_LEAK = 0,   // cumulative resisting-sign injection work L
        LNVD_WORK   = 1,   // cumulative FLAC local-damping dissipation
        MODAL_WORK  = 2,   // reserved (ADR-69 P2)
        NUM_CHANNELS = 3
    };

    static EnergyChannelRegistry &instance() {
        static EnergyChannelRegistry theRegistry;
        return theRegistry;
    }

    // Mark a channel relevant (column emitted) without adding energy.
    void declare(Channel c) {
        if (c >= 0 && c < NUM_CHANNELS)
            isDeclared[c] = true;
    }

    // Accumulate dE (units of work) into a channel; declares it as a side
    // effect so a producer cannot publish into an invisible channel.
    // FINITENESS GUARD (ADR-72 P3 CI find): totals are process-lifetime and
    // survive wipe() BY DESIGN (consumers baseline-delta), so a single
    // non-finite publication — e.g. an intentionally diverging explicit run
    // with LNVD active — would poison the channel for EVERY later model in
    // the process (NaN - NaN baseline => NaN RES in recorders whose model
    // never produced the channel at all). Discard non-finite dE with a
    // process-once warning; the diverged run's own energy report is already
    // non-finite through the domain sweep, so nothing honest is lost.  // Ladruno
    void addEnergy(Channel c, double dE) {
        if (c >= 0 && c < NUM_CHANNELS) {
            isDeclared[c] = true;
            if (!(dE == dE) || dE == std::numeric_limits<double>::infinity()
                            || dE == -std::numeric_limits<double>::infinity()) {
                if (!warnedNonFinite) {
                    warnedNonFinite = true;
                    opserr << "WARNING Ladruno EnergyChannelRegistry: discarded a "
                              "non-finite energy publication (diverged run?) - "
                              "channel totals stay finite so later models in this "
                              "process keep an honest balance. (printed once)\n";
                }
                return;
            }
            total[c] += dE;
        }
    }

    bool declared(Channel c) const {
        return (c >= 0 && c < NUM_CHANNELS) ? isDeclared[c] : false;
    }

    bool anyDeclared() const {
        for (int i = 0; i < NUM_CHANNELS; ++i)
            if (isDeclared[i])
                return true;
        return false;
    }

    // Cumulative total since the last model wipe; consumers subtract their
    // own first-record baseline (see header comment).
    double energy(Channel c) const {
        return (c >= 0 && c < NUM_CHANNELS) ? total[c] : 0.0;
    }

    // Called from Domain::clearAll() (the shared wipe path of every
    // interpreter). ADR-72 P3 CI find: PROCESS-lifetime totals broke in two
    // ways across models in one process — (a) a diverged run publishing
    // non-finite work poisoned every later model's RES (NaN - NaN survives
    // the baseline subtraction), and (b) a huge-but-finite total (~1e300
    // pre-overflow) ABSORBS a later model's ~1e-3 increments in double
    // precision (total + dE == total), silently zeroing its channel delta.
    // wipe() destroys every producer and consumer, so restarting the totals
    // there makes the original design promise ("a new recorder on a rebuilt
    // model starts from zero") literally true at full precision. Declared
    // flags reset too — channel columns are model-scoped.  // Ladruno
    void resetOnWipe() {
        for (int i = 0; i < NUM_CHANNELS; ++i) {
            total[i] = 0.0;
            isDeclared[i] = false;
        }
    }

private:
    EnergyChannelRegistry() {
        for (int i = 0; i < NUM_CHANNELS; ++i) {
            total[i] = 0.0;
            isDeclared[i] = false;
        }
        warnedNonFinite = false;
    }
    EnergyChannelRegistry(const EnergyChannelRegistry &);
    EnergyChannelRegistry &operator=(const EnergyChannelRegistry &);

    double total[NUM_CHANNELS];
    bool isDeclared[NUM_CHANNELS];
    bool warnedNonFinite;              // process-once guard warning  // Ladruno
};

} // namespace Ladruno

#endif // LadrunoEnergyChannels_h
