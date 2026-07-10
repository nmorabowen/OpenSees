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

// Ladruno: modal-superposition frequency/response-domain toolkit (ADR 44).
// See Ladruno_implementation/44_ladruno_frequency_domain_adr.md.
//
//   P1a scope: modalResponseHistory — the EXACT piecewise-linear (PWL)
//   modal-superposition transient.  Once `eigen` + `modalProperties` have run,
//   each retained mode is a decoupled SDOF
//
//       q_a'' + d_a q_a' + w_a^2 q_a = f_a(t),   d_a = c_a/m_a,
//
//   integrated by a per-mode constant recurrence
//
//       [q_{n+1}; q'_{n+1}] = A_a [q_n; q'_n] + B_a [f_n; f_{n+1}]
//
//   that is EXACT for a load varying linearly within each step (Nigam-Jennings
//   / first-order hold).  The physical response is recovered by the modal sum
//   u = sum_a (phi_a * Vscale) q_a and one domain step is COMMITTED per time
//   station so ordinary recorders capture displacement / velocity /
//   acceleration / element / reaction histories exactly as in a direct run.
//
//   Excitation (P1a): uniform base acceleration u_g''(t) along a global
//   direction, relative formulation, modal load f_a = -Gamma_a u_g''(t) with
//   Gamma_a the participation factor already computed by modalProperties.
//   General load-pattern modal loads and the FRF/SSD/random post-processors
//   (P2/P3) are separate phases.
//
//   The (w, d) parametrization (rather than (w, xi)) unifies all four damping
//   branches by the discriminant Delta = w^2 - (d/2)^2 (under/critical/over)
//   AND covers a rigid-body mode carrying Rayleigh mass-proportional damping
//   (w=0, d=a0 > 0), which a pure q''=f rigid branch would silently miss.
//   Derivation + machine-precision oracle:
//   Ladruno_implementation/modal_response_p1a_spike/pwl_recurrence_oracle.py.

#ifndef LadrunoModalResponse_h
#define LadrunoModalResponse_h

#include <vector>
#include <string>

class AnalysisModel;
class Domain;
class TimeSeries;
class DomainModalProperties;

// ---------------------------------------------------------------------------
// Shared modal damping specification: how the per-mode d_a = c_a/m_a (= 2 xi w
// for w>0) is obtained.  Reused by the P2 frequency-domain post-processors; the
// P1a transient carries its own equivalent copy (kept byte-stable / untouched).
// ---------------------------------------------------------------------------
struct LadrunoModalDamping
{
    enum Kind { UNIFORM, RAYLEIGH, MODAL_LIST };
    Kind                kind = UNIFORM;
    double              xi   = 0.0;   // uniform ratio
    double              a0   = 0.0;   // Rayleigh mass factor
    double              a1   = 0.0;   // Rayleigh stiffness factor
    std::vector<double> xis;          // per absolute-mode ratios (MODAL_LIST)

    // d_a = c_a/m_a for absolute 0-based mode index at circular frequency w.
    // Rayleigh: a0 + a1 w^2 (rigid mode w=0 -> a0). Ratio channels: 2 xi w
    // (rigid mode -> 0, stiffness-proportional damping vanishing there).
    double coeff(int mode0, double w) const;
};

class LadrunoModalResponse
{
public:
    // damping specification (how the per-mode d_a = c_a/m_a is obtained)
    enum DampingKind { DAMP_UNIFORM, DAMP_RAYLEIGH, DAMP_MODAL_LIST };

    LadrunoModalResponse(AnalysisModel* theModel,
                         TimeSeries*    baseAccel,
                         int            direction,
                         double         dt,
                         int            nsteps,
                         double         t0);
    ~LadrunoModalResponse();

    // damping channels (mutually exclusive; last one set wins)
    void setDampingUniform(double xi);
    void setDampingRayleigh(double a0, double a1);
    void setDampingModalList(const std::vector<double>& xi_per_mode); // absolute mode order

    // optional 1-based mode subset (default: all modes from the last eigen run)
    void setModeSubset(const std::vector<int>& modes_1based);

    // run the whole transient, committing one domain step per time station
    int analyze();

    // exposed for the unit oracle: exact PWL recurrence blocks for one SDOF.
    //   A[2][2], B[2][2] s.t. [q1;qd1] = A[q0;qd0] + B[f0;f1].
    // w  = undamped circular frequency (>= 0), d = c/m damping coefficient (>= 0).
    static void recurrenceBlocks(double w, double d, double dt,
                                 double A[2][2], double B[2][2]);

private:
    int  check(DomainModalProperties& mp);
    // per-mode damping coefficient d_a = c_a/m_a for absolute mode index (0-based)
    double dampingCoeff(int mode0, double w) const;

private:
    AnalysisModel* m_model;
    TimeSeries*    m_baseAccel;   // ground acceleration time series (base motion)
    int            m_direction;   // 1-based global excitation direction
    double         m_dt;
    int            m_nsteps;
    double         m_t0;          // start time (base-accel sampled at t0 + n*dt)

    DampingKind         m_damp_kind;
    double              m_xi;      // uniform ratio
    double              m_a0, m_a1; // Rayleigh factors
    std::vector<double> m_xi_list; // per absolute-mode ratios

    std::vector<int>    m_modes;   // absolute 0-based mode indices to include
};

// ===========================================================================
// ADR 44 P2: modal frequency-response / steady-state dynamics.
//
//   Once `eigen` + `modalProperties` have run, the steady harmonic response of
//   a LINEAR structure to a uniform base acceleration u_g''(t) = amp e^{+iOm t}
//   (relative formulation, ADR 44 §4.3-§4.4) is the modal sum
//
//       u_hat(Om) = amp * sum_a (phi_a Vscale_a) (-Gamma_a) H_a(Om),
//       H_a(Om)   = 1 / (w_a^2 - Om^2 + i Om d_a),   d_a = c_a/m_a.
//
//   This is the frequency-domain image of the P1a transient recovery (same
//   psi_a = phi_a Vscale_a and Gamma_a = participation factor, both validated
//   there vs direct Newmark), so it inherits P1a's mass-normalization exactly.
//   Velocity FRF = i Om u_hat; relative-acceleration FRF = -Om^2 u_hat.
//
//   Sign convention is PINNED to e^{+iOm t} (the +i Om d_a term gives the
//   physical -90 deg response lag at resonance).  Verified end-to-end vs the
//   direct complex solve (K - Om^2 M + i Om C)^{-1}(-M R) in
//   Ladruno_implementation/modal_response_p2_spike/frf_oracle.py.
//
//   Frequencies are in Hz (Om = 2 pi f) throughout, in and out.
// ===========================================================================
class LadrunoFrequencyResponse
{
public:
    enum SweepKind { SWEEP_LIN, SWEEP_LOG, SWEEP_BIASED };
    enum RespKind  { RESP_DISP, RESP_VEL, RESP_ACCEL };

    explicit LadrunoFrequencyResponse(AnalysisModel* theModel);
    ~LadrunoFrequencyResponse();

    // excitation: uniform base acceleration of amplitude amp along a global dir
    void setBaseAccel(int direction, double amp) { m_direction = direction; m_amp = amp; }
    // frequency sweep in Hz
    void setSweep(double fmin, double fmax, int nf, SweepKind k)
    { m_fmin = fmin; m_fmax = fmax; m_nf = nf; m_sweep = k; }
    void setDamping(const LadrunoModalDamping& d) { m_damp = d; }
    void setModeSubset(const std::vector<int>& modes_1based);
    void setResponse(int nodeTag, int dof, RespKind r)
    { m_nodeTag = nodeTag; m_dof = dof; m_resp = r; }
    void setOutputFile(const char* fname) { m_outfile = fname ? fname : ""; }

    // Evaluate the sweep.  Fills `rows`, one per frequency:
    //   ampOnly==false -> {f, Re, Im}   (full complex FRF)
    //   ampOnly==true  -> {f, |H|}      (steady-state amplitude, SSD)
    // Also writes `rows` to the output file if one was set. Returns 0 / -1.
    int run(std::vector<std::vector<double> >& rows, bool ampOnly);

private:
    int  buildGrid(std::vector<double>& freqs) const; // Hz sample points

    AnalysisModel*      m_model;
    int                 m_direction; // 1-based global excitation dir
    double              m_amp;
    double              m_fmin, m_fmax;
    int                 m_nf;
    SweepKind           m_sweep;
    LadrunoModalDamping m_damp;
    std::vector<int>    m_modes;     // absolute 0-based mode subset (empty=all)
    int                 m_nodeTag;
    int                 m_dof;       // 1-based response DOF
    RespKind            m_resp;
    std::string         m_outfile;
};

#endif
