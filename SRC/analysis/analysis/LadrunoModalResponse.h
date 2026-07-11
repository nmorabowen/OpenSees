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
//   Excitation channels:
//     base accel (P1a): uniform u_g''(t) along a global direction, relative
//         formulation, modal load f_a = -Gamma_a u_g''(t) with Gamma_a the
//         participation factor already computed by modalProperties.
//     -load (follow-up): nodal forces P(t) = s(t) * P — P is a load pattern's
//         plain-NodalLoad shape (shared #553 assembly helper; the pattern's
//         own TimeSeries is IGNORED) and s(t) the ctor's time series sampled
//         at the stations; modal load f_a(t) = (psi_a^T P / m~_a) s(t) with
//         the #553-pinned normalization (psi = evec*Vscale, m~ =
//         generalizedMasses()).  Response is ABSOLUTE for -load (the relative
//         formulation is a base-motion concept).
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

    // switch the excitation to a nodal-force pattern P(t) = s(t)*P, where P is
    // the pattern's plain-NodalLoad shape and s(t) is the ctor's time series
    // (the `baseAccel` ctor argument doubles as s; `direction` is unused)
    void setLoadPattern(int patternTag)
    { m_useLoad = true; m_patternTag = patternTag; }

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
    TimeSeries*    m_baseAccel;   // excitation series: u_g''(t), or s(t) w/ -load
    int            m_direction;   // 1-based global excitation direction (base)
    bool           m_useLoad = false; // excitation = nodal-force pattern
    int            m_patternTag = 0;  // load-pattern tag (m_useLoad)
    double         m_dt;
    int            m_nsteps;
    double         m_t0;          // start time (series sampled at t0 + n*dt)

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
//
//   Excitation channels (ADR 44 §5.2, -load follow-up):
//     EX_BASE_ACCEL  uniform base accel, modal load  -Gamma_a  (P2, above)
//     EX_LOAD        nodal-force pattern P e^{+iOm t}: the spatial shape P is
//                    assembled from the pattern's plain NodalLoads (reference
//                    values; the pattern's own TimeSeries is IGNORED — the
//                    sweep replaces it) and the modal load is
//
//                        f_a = psi_a^T P / m~_a,   m~_a = psi_a^T M psi_a,
//
//                    where psi_a = evec_a*Vscale_a is EXACTLY the basis
//                    DomainModalProperties computes generalizedMasses() in
//                    (the -unorm scaling is applied to its local V BEFORE the
//                    GM/MPF products — DomainModalProperties.cpp ~:636-772),
//                    so m~_a == generalizedMasses()(a) directly.  -baseAccel
//                    is the special case P = -M R (then f_a == -Gamma_a).
//                    Normalization pinned vs the direct solve in
//                    modal_response_p3_spike/load_frf_oracle.py.
// ===========================================================================
class LadrunoFrequencyResponse
{
public:
    enum SweepKind { SWEEP_LIN, SWEEP_LOG, SWEEP_BIASED };
    enum RespKind  { RESP_DISP, RESP_VEL, RESP_ACCEL };
    enum ExKind    { EX_BASE_ACCEL, EX_LOAD };

    explicit LadrunoFrequencyResponse(AnalysisModel* theModel);
    ~LadrunoFrequencyResponse();

    // excitation: uniform base acceleration of amplitude amp along a global dir
    void setBaseAccel(int direction, double amp)
    { m_ex = EX_BASE_ACCEL; m_direction = direction; m_amp = amp; }
    // excitation: harmonic nodal-force pattern, amp * P * e^{+iOm t}
    void setLoadPattern(int patternTag, double amp)
    { m_ex = EX_LOAD; m_patternTag = patternTag; m_amp = amp; }
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
    ExKind              m_ex;        // excitation channel
    int                 m_direction; // 1-based global excitation dir (EX_BASE_ACCEL)
    int                 m_patternTag;// load pattern tag (EX_LOAD)
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

// ===========================================================================
// ADR 44 P3: stationary random response (PSD -> RMS).
//
//   For a stationary base-acceleration excitation given as a ONE-SIDED power
//   spectral density G(f) in Hz (sigma_ug^2 = int_0^inf G df — the engineering
//   spec convention for wind / floor-vibration / equipment curves), the
//   response auto-PSD of the requested quantity rides the P2-validated FRF:
//
//       G_xx(f) = |H_x(f)|^2 G(f),
//
//   where H_x(f) is EXACTLY what LadrunoFrequencyResponse::run() emits for the
//   same node/dof/resp (relative disp / vel / accel).  The engineering outputs
//   are then band integrals over the sweep grid [fmin, fmax] (trapezoid):
//
//       m_k  = int f^k G_xx(f) df                (spectral moments, k = 0, 2)
//       RMS  = sqrt(m_0)                         (zero-mean stationary response)
//       nu_0 = sqrt(m_2/m_0)  [Hz]               (mean zero-upcrossing rate)
//       E[peak] = g * RMS  over a duration T     (Davenport peak factor,
//           g = sqrt(2 ln(nu0 T)) + 0.5772/sqrt(2 ln(nu0 T)), needs nu0*T > 1)
//
//   CONVENTION PIN (the classic silent bug is a one-sided/two-sided or Hz /
//   rad-per-s factor): vs the textbook TWO-SIDED rad/s PSD, G(f) = 4 pi S(Om).
//   White-noise SDOF anchor in this convention: sigma_x^2 = G0/(8 xi w^3).
//   Pinned empirically (Monte-Carlo synthetic realization, 0.6% agreement) in
//   Ladruno_implementation/modal_response_p3_spike/psd_rms_oracle.py — a wrong
//   one-sided factor would read ~41% (sqrt 2), a Hz/rad mixup ~150%.
//
//   A zero-damped mode whose resonance lies INSIDE the sweep band is refused:
//   its white-noise variance integral diverges, so any finite number the grid
//   returned would be meaningless (and a -biased grid samples exactly on the
//   resonance -> inf row).  Zero-damped modes strictly outside the band keep a
//   finite in-band integrand and are legal.
//
//   -load channel (follow-up): the excitation is P(t) = s(t) * P with P the
//   pattern's nodal-load shape and s(t) a stationary scalar of one-sided PSD
//   G_s(f) [Hz] (all loads fully correlated); H_x is then the -load FRF and
//   everything downstream (G_xx = |H_x|^2 G_s, moments, stats) is unchanged.
//   SDOF anchor: white force PSD G_F -> sigma_u^2 = G_F/(8 xi w^3 m^2)
//   (modal_response_p3_spike/load_frf_oracle.py check 4).
// ===========================================================================
class LadrunoRandomResponse
{
public:
    struct Stats {
        double rms  = 0.0;   // sqrt(m0)
        double m0   = 0.0;   // int   G_xx df   over the band
        double m2   = 0.0;   // int f^2 G_xx df over the band
        double nu0  = 0.0;   // sqrt(m2/m0), Hz
        double peak = 0.0;   // Davenport expected peak (valid iff hasPeak)
        bool   hasPeak = false;
    };

    explicit LadrunoRandomResponse(AnalysisModel* theModel);
    ~LadrunoRandomResponse();

    // excitation: uniform base acceleration along a global dir, with one-sided
    // PSD G(f) [Hz] sampled from the time series (getFactor(f), f in Hz)
    void setBaseAccel(int direction)
    { m_ex = LadrunoFrequencyResponse::EX_BASE_ACCEL; m_direction = direction; }
    // excitation: nodal-force pattern shape P scaled by a stationary scalar
    // s(t) with one-sided PSD G_s(f) [Hz]: P(t) = s(t)*P (fully correlated)
    void setLoadPattern(int patternTag)
    { m_ex = LadrunoFrequencyResponse::EX_LOAD; m_patternTag = patternTag; }
    void setInputPSD(TimeSeries* psd)    { m_psd = psd; }
    void setSweep(double fmin, double fmax, int nf,
                  LadrunoFrequencyResponse::SweepKind k)
    { m_fmin = fmin; m_fmax = fmax; m_nf = nf; m_sweep = k; }
    void setDamping(const LadrunoModalDamping& d) { m_damp = d; }
    void setModeSubset(const std::vector<int>& modes_1based)
    { m_modes_1based = modes_1based; }
    void setResponse(int nodeTag, int dof, LadrunoFrequencyResponse::RespKind r)
    { m_nodeTag = nodeTag; m_dof = dof; m_resp = r; }
    void setOutputFile(const char* fname) { m_outfile = fname ? fname : ""; }
    void setDuration(double T)            { m_duration = T; } // for E[peak]

    // Evaluate: FRF sweep -> G_xx rows {f, G_in, G_xx} -> band moments/stats.
    // Also writes `rows` to the output file if one was set. Returns 0 / -1.
    int run(std::vector<std::vector<double> >& rows, Stats& out);

private:
    AnalysisModel*      m_model;
    TimeSeries*         m_psd;
    LadrunoFrequencyResponse::ExKind m_ex;
    int                 m_direction;
    int                 m_patternTag;
    double              m_fmin, m_fmax;
    int                 m_nf;
    LadrunoFrequencyResponse::SweepKind m_sweep;
    LadrunoModalDamping m_damp;
    std::vector<int>    m_modes_1based;  // pass-through to the FRF (empty=all)
    int                 m_nodeTag;
    int                 m_dof;
    LadrunoFrequencyResponse::RespKind  m_resp;
    std::string         m_outfile;
    double              m_duration;      // <=0: no peak estimate
};

#endif
