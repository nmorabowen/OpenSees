/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell (Ladruño project)
//
// References:
//   Nigam, N.C. and Jennings, P.C. (1969). Calculation of response spectra
//     from strong-motion earthquake records. Bulletin of the Seismological
//     Society of America, 59(2), 909-922. doi:10.1785/BSSA0590020909
//   Van Loan, C.F. (1978). Computing integrals involving the matrix
//     exponential. IEEE Transactions on Automatic Control, 23(3), 395-404.
//     doi:10.1109/TAC.1978.1101743
//   Chopra, A.K. (2017). Dynamics of Structures, 5th ed. Pearson
//     (Sec. 5.2, piecewise-exact recurrence).
//
// Description: modal-superposition response of a linear model on the modes
// of the last `eigen` + `modalProperties` run.
//
//   modalResponseHistory: the piecewise-linear modal-superposition
//   transient.  Each retained mode is a decoupled SDOF
//
//       q_a'' + d_a q_a' + w_a^2 q_a = f_a(t),   d_a = c_a/m_a,
//
//   integrated by a per-mode constant recurrence
//
//       [q_{n+1}; q'_{n+1}] = A_a [q_n; q'_n] + B_a [f_n; f_{n+1}]
//
//   that is exact for a load varying linearly within each step (Nigam and
//   Jennings 1969; first-order hold).  The physical response is recovered by
//   the modal sum u = sum_a (phi_a * Vscale) q_a and one domain step is
//   committed per time station, so ordinary recorders capture displacement,
//   velocity, acceleration, element and reaction histories as in a direct
//   run.
//
//   Excitation channels:
//     -baseAccel: uniform u_g''(t) along a global direction, relative
//         formulation, modal load f_a = -Gamma_a u_g''(t) with Gamma_a the
//         participation factor computed by modalProperties.
//     -load: nodal forces P(t) = s(t) * P.  P is a load pattern's
//         plain-NodalLoad shape (the pattern's own TimeSeries is ignored) and
//         s(t) the -series time series sampled at the stations; modal load
//         f_a(t) = (psi_a^T P / m~_a) s(t) with psi = evec*Vscale and
//         m~ = generalizedMasses(), independent of the eigenvector
//         normalization.  The response is absolute for -load (the relative
//         formulation is a base-motion concept).
//
//   The (w, d) parametrization (rather than (w, xi)) covers all damping
//   branches through the discriminant Delta = w^2 - (d/2)^2 (under,
//   critical, over) and a rigid-body mode carrying Rayleigh mass-proportional
//   damping (w=0, d=a0 > 0), which a pure q''=f rigid branch would miss.

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
// for w>0) is obtained.  Used by the frequency-domain commands; the transient
// carries its own equivalent copy.
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


#endif
