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
//   * Totals are MONOTONE over process lifetime; there is deliberately no
//     per-step clear and no reset-on-wipe hook. The recorder snapshots each
//     channel at its first record() and reports the DELTA since — so a new
//     recorder on a rebuilt model starts from zero regardless of history.
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
//   * MODAL_WORK is reserved for the modal-damping publisher (ADR-69 P2);
//     nothing declares it yet.
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
    void addEnergy(Channel c, double dE) {
        if (c >= 0 && c < NUM_CHANNELS) {
            isDeclared[c] = true;
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

    // Cumulative process-lifetime total; consumers subtract their own
    // first-record baseline (see header comment).
    double energy(Channel c) const {
        return (c >= 0 && c < NUM_CHANNELS) ? total[c] : 0.0;
    }

private:
    EnergyChannelRegistry() {
        for (int i = 0; i < NUM_CHANNELS; ++i) {
            total[i] = 0.0;
            isDeclared[i] = false;
        }
    }
    EnergyChannelRegistry(const EnergyChannelRegistry &);
    EnergyChannelRegistry &operator=(const EnergyChannelRegistry &);

    double total[NUM_CHANNELS];
    bool isDeclared[NUM_CHANNELS];
};

} // namespace Ladruno

#endif // LadrunoEnergyChannels_h
