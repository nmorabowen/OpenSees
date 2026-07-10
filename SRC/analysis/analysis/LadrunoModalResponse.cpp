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

// Ladruno ADR 44 P1a: modalResponseHistory — exact PWL modal-superposition
// transient.  See LadrunoModalResponse.h for the formulation and
// Ladruno_implementation/44_ladruno_frequency_domain_adr.md.

#include <LadrunoModalResponse.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <DomainModalProperties.h>
#include <TimeSeries.h>
#include <Node.h>
#include <NodeIter.h>
#include <Vector.h>
#include <Matrix.h>
#include <elementAPI.h>
#include <cmath>
#include <cstring>
#include <complex>
#include <fstream>
#include <iomanip>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// shared modal-damping coefficient d_a = c_a/m_a (see header).
double
LadrunoModalDamping::coeff(int mode0, double w) const
{
    switch (kind) {
    case RAYLEIGH:
        return a0 + a1 * w * w;              // == 2 xi w for w>0; a0 for a rigid mode
    case MODAL_LIST: {
        double x = 0.0;
        if (mode0 >= 0 && mode0 < static_cast<int>(xis.size()))
            x = xis[mode0];
        return 2.0 * x * w;                  // rigid mode -> 0
    }
    case UNIFORM:
    default:
        return 2.0 * xi * w;                 // rigid mode -> 0
    }
}

// ---------------------------------------------------------------------------
// Exact PWL SDOF recurrence blocks in (w, d) form.
//   q'' + d q' + w^2 q = f(t),  f linear on the step.
//   [q1;qd1] = A [q0;qd0] + B [f0;f1],   exact for piecewise-linear f.
// Unified over all damping branches by Delta = w^2 - (d/2)^2, and covering the
// rigid-body Rayleigh case (w=0, d>0).  Verified to ~1e-14 vs a van-Loan
// first-order-hold expm oracle across 40 (w,d) cases + rigid-Rayleigh steady
// state (modal_response_p1a_spike/pwl_recurrence_oracle.py).
// ---------------------------------------------------------------------------
void
LadrunoModalResponse::recurrenceBlocks(double w, double d, double dt,
                                       double A[2][2], double B[2][2])
{
    const double CRIT_TOL = 1.0e-8; // relative closeness to the critical root
    const double half_d = 0.5 * d;

    // --- homogeneous propagator A -----------------------------------------
    double p11, p12, p21, p22;
    if (w == 0.0) {
        if (d == 0.0) {
            // double integrator
            p11 = 1.0; p12 = dt; p21 = 0.0; p22 = 1.0;
        } else {
            // roots {0, -d}: rigid mode with mass-proportional damping
            const double e = std::exp(-d * dt);
            p11 = 1.0; p12 = (1.0 - e) / d; p21 = 0.0; p22 = e;
        }
    } else {
        const double e = std::exp(-half_d * dt);
        const double Delta = w * w - half_d * half_d; // w^2 - (d/2)^2
        double c, sd;                                 // sd = sin(wd dt)/wd
        if (std::abs(Delta) <= CRIT_TOL * (w * w)) {
            c = 1.0; sd = dt;                          // critically damped limit
        } else if (Delta > 0.0) {
            const double wd = std::sqrt(Delta);        // underdamped
            c = std::cos(wd * dt); sd = std::sin(wd * dt) / wd;
        } else {
            const double wh = std::sqrt(-Delta);       // overdamped
            c = std::cosh(wh * dt); sd = std::sinh(wh * dt) / wh;
        }
        p11 = e * (c + half_d * sd);
        p12 = e * sd;
        p21 = e * (-w * w * sd);
        p22 = e * (c - half_d * sd);
    }
    A[0][0] = p11; A[0][1] = p12;
    A[1][0] = p21; A[1][1] = p22;

    // --- load block B -----------------------------------------------------
    // particular ramp solution q_p, cancelled to REST at t=0 by the
    // homogeneous part, then propagated by A to t=dt.
    if (w == 0.0 && d == 0.0) {
        B[0][0] = dt * dt / 3.0; B[0][1] = dt * dt / 6.0;
        B[1][0] = dt / 2.0;      B[1][1] = dt / 2.0;
        return;
    }
    for (int j = 0; j < 2; ++j) {
        const double f0 = (j == 0) ? 1.0 : 0.0;
        const double f1 = (j == 0) ? 0.0 : 1.0;
        const double r  = (f1 - f0) / dt;
        double qp0, qpd0, qpDt, qpdDt;
        if (w == 0.0) {
            // w=0, d>0:  v'=q', v' + d v = f0 + r t
            //   v_p(t) = (f0 + r t)/d - r/d^2
            //   q_p(t) = (f0 t + r t^2/2)/d - r t/d^2      (q_p(0)=0)
            const double invd = 1.0 / d, invd2 = invd * invd;
            qp0   = 0.0;
            qpd0  = f0 * invd - r * invd2;
            qpDt  = (f0 * dt + r * dt * dt * 0.5) * invd - r * dt * invd2;
            qpdDt = (f0 + r * dt) * invd - r * invd2;
        } else {
            // w>0, any d:  q_p(t) = (f0 + r t)/w^2 - d r / w^4
            const double w2 = w * w, invw2 = 1.0 / w2, invw4 = invw2 * invw2;
            qp0   = f0 * invw2 - d * r * invw4;
            qpd0  = r * invw2;
            qpDt  = (f0 + r * dt) * invw2 - d * r * invw4;
            qpdDt = r * invw2;
        }
        const double qh0  = -qp0;
        const double qhd0 = -qpd0;
        const double qh_dt  = p11 * qh0 + p12 * qhd0;
        const double qhd_dt = p21 * qh0 + p22 * qhd0;
        B[0][j] = qpDt  + qh_dt;
        B[1][j] = qpdDt + qhd_dt;
    }
}

// ---------------------------------------------------------------------------
LadrunoModalResponse::LadrunoModalResponse(AnalysisModel* theModel,
                                           TimeSeries*    baseAccel,
                                           int            direction,
                                           double         dt,
                                           int            nsteps,
                                           double         t0)
    : m_model(theModel)
    , m_baseAccel(baseAccel)
    , m_direction(direction)
    , m_dt(dt)
    , m_nsteps(nsteps)
    , m_t0(t0)
    , m_damp_kind(DAMP_UNIFORM)
    , m_xi(0.0)
    , m_a0(0.0), m_a1(0.0)
{
}

LadrunoModalResponse::~LadrunoModalResponse() {}

void LadrunoModalResponse::setDampingUniform(double xi)
{ m_damp_kind = DAMP_UNIFORM; m_xi = xi; }

void LadrunoModalResponse::setDampingRayleigh(double a0, double a1)
{ m_damp_kind = DAMP_RAYLEIGH; m_a0 = a0; m_a1 = a1; }

void LadrunoModalResponse::setDampingModalList(const std::vector<double>& xis)
{ m_damp_kind = DAMP_MODAL_LIST; m_xi_list = xis; }

void LadrunoModalResponse::setModeSubset(const std::vector<int>& modes_1based)
{
    m_modes.clear();
    for (int m : modes_1based) m_modes.push_back(m - 1); // to 0-based
}

double LadrunoModalResponse::dampingCoeff(int mode0, double w) const
{
    switch (m_damp_kind) {
    case DAMP_RAYLEIGH:
        // c/m = a0 + a1 w^2   (== 2 xi w for w>0; == a0 for a rigid mode)
        return m_a0 + m_a1 * w * w;
    case DAMP_MODAL_LIST: {
        double xi = 0.0;
        if (mode0 >= 0 && mode0 < static_cast<int>(m_xi_list.size()))
            xi = m_xi_list[mode0];
        return 2.0 * xi * w; // rigid mode -> 0
    }
    case DAMP_UNIFORM:
    default:
        return 2.0 * m_xi * w; // rigid mode -> 0
    }
}

int LadrunoModalResponse::check(DomainModalProperties& mp)
{
    Domain* domain = m_model->getDomainPtr();
    const Vector& ev = domain->getEigenvalues();
    if (ev.Size() < 1) {
        opserr << "LadrunoModalResponse - No eigenvalues; call 'eigen' first.\n";
        return -1;
    }
    if (ev.Size() != mp.eigenvalues().Size()) {
        opserr << "LadrunoModalResponse - modalProperties is stale vs the last "
                  "eigen run; call 'modalProperties' right after 'eigen'.\n";
        return -1;
    }
    double tol = std::max(1.0e-15, 1.0e-12 * ev.Norm());
    for (int i = 0; i < ev.Size(); ++i)
        if (std::abs(ev(i) - mp.eigenvalues()(i)) > tol) {
            opserr << "LadrunoModalResponse - eigenvalue mismatch between the "
                      "domain and DomainModalProperties (re-run modalProperties "
                      "after eigen).\n";
            return -1;
        }
    return 0;
}

int LadrunoModalResponse::analyze()
{
    if (m_model == 0 || m_model->getDomainPtr() == 0) {
        opserr << "LadrunoModalResponse - no AnalysisModel/Domain.\n";
        return -1;
    }
    Domain* domain = m_model->getDomainPtr();

    if (m_baseAccel == 0) {
        opserr << "LadrunoModalResponse - a base-acceleration time series is "
                  "required (-baseAccel).\n";
        return -1;
    }
    if (m_dt <= 0.0 || m_nsteps < 1) {
        opserr << "LadrunoModalResponse - need dt>0 and nsteps>=1.\n";
        return -1;
    }

    DomainModalProperties mp;
    if (domain->getModalProperties(mp) < 0) {
        opserr << "LadrunoModalResponse - eigen and modalProperties have not "
                  "been called.\n";
        return -1;
    }
    if (check(mp) < 0)
        return -1;

    const int num_eigen = domain->getEigenvalues().Size();
    const int ndf = mp.totalMass().Size();
    if (m_direction < 1 || m_direction > ndf) {
        opserr << "LadrunoModalResponse - direction (" << m_direction
               << ") must be in 1.." << ndf << ".\n";
        return -1;
    }
    const int exdof = m_direction - 1;

    // active mode list (default: all)
    std::vector<int> modes = m_modes;
    if (modes.empty()) {
        modes.resize(num_eigen);
        for (int i = 0; i < num_eigen; ++i) modes[i] = i;
    }
    for (int m : modes)
        if (m < 0 || m >= num_eigen) {
            opserr << "LadrunoModalResponse - requested mode " << (m + 1)
                   << " is out of range 1.." << num_eigen << ".\n";
            return -1;
        }
    // reject duplicate modes: the recovery sums per-mode contributions, so a
    // repeated index would silently double-count that mode's response.
    {
        std::vector<int> seen = modes;
        std::sort(seen.begin(), seen.end());
        if (std::adjacent_find(seen.begin(), seen.end()) != seen.end()) {
            opserr << "LadrunoModalResponse - -modes contains a duplicate index "
                      "(each mode may appear at most once).\n";
            return -1;
        }
    }
    if (m_damp_kind == DAMP_MODAL_LIST) {
        int maxmode = 0;
        for (int m : modes) maxmode = std::max(maxmode, m);
        if (static_cast<int>(m_xi_list.size()) <= maxmode) {
            opserr << "LadrunoModalResponse - -modalDamp list has "
                   << static_cast<int>(m_xi_list.size()) << " ratios but mode "
                   << (maxmode + 1) << " is used.\n";
            return -1;
        }
    }
    const int nm = static_cast<int>(modes.size());

    // negative-eigenvalue guard: a small negative lambda is numerical noise on a
    // rigid-body mode (treat as w=0); a genuinely negative one means the eigen
    // run is ill-posed and the modal transient would be meaningless — refuse
    // rather than silently reinterpret it as rigid.
    double maxAbsLambda = 0.0;
    for (int m : modes)
        maxAbsLambda = std::max(maxAbsLambda, std::abs(mp.eigenvalues()(m)));
    const double negTol = std::max(1.0e-12, 1.0e-8 * maxAbsLambda);

    // per-mode precompute
    std::vector<double> w(nm), Gamma(nm), Vscale(nm), dcoef(nm);
    std::vector<double> Amat(nm * 4), Bmat(nm * 4);
    std::vector<double> q(nm, 0.0), qd(nm, 0.0), qdd(nm, 0.0);
    for (int a = 0; a < nm; ++a) {
        const int mode0 = modes[a];
        const double lambda = mp.eigenvalues()(mode0);
        if (lambda < -negTol) {
            opserr << "LadrunoModalResponse - mode " << (mode0 + 1)
                   << " has a negative eigenvalue (" << lambda << "); the eigen "
                      "result is ill-posed for a modal transient.\n";
            return -1;
        }
        const double wa = (lambda > 0.0) ? std::sqrt(lambda) : 0.0; // clamp rigid
        w[a] = wa;
        Gamma[a] = mp.modalParticipationFactors()(mode0, exdof);
        Vscale[a] = mp.eigenVectorScaleFactors()(mode0);
        dcoef[a] = dampingCoeff(mode0, wa);
        double A[2][2], B[2][2];
        recurrenceBlocks(wa, dcoef[a], m_dt, A, B);
        Amat[a * 4 + 0] = A[0][0]; Amat[a * 4 + 1] = A[0][1];
        Amat[a * 4 + 2] = A[1][0]; Amat[a * 4 + 3] = A[1][1];
        Bmat[a * 4 + 0] = B[0][0]; Bmat[a * 4 + 1] = B[0][1];
        Bmat[a * 4 + 2] = B[1][0]; Bmat[a * 4 + 3] = B[1][1];
    }

    // --- station loop: reconstruct u,udot,uddot and COMMIT one step each ---
    // station 0 is the initial rest state at t0; then nsteps advances.
    for (int n = 0; n <= m_nsteps; ++n) {
        const double t = m_t0 + n * m_dt;
        const double ug = m_baseAccel->getFactor(t);

        // modal accelerations q'' from the ODE (for the accel recovery)
        for (int a = 0; a < nm; ++a) {
            const double f_a = -Gamma[a] * ug;
            qdd[a] = f_a - dcoef[a] * qd[a] - w[a] * w[a] * q[a];
        }

        // physical reconstruction  u = sum_a (phi_a Vscale) q_a  (+ vel/accel)
        Node* node;
        NodeIter& theNodes = domain->getNodes();
        while ((node = theNodes()) != 0) {
            const Matrix& evec = node->getEigenvectors();
            const int node_ndf = evec.noRows();
            if (node_ndf < 1) continue;

            Vector ud(node_ndf), uv(node_ndf), ua(node_ndf);
            ud.Zero(); uv.Zero(); ua.Zero();

            const int nd = std::min(node_ndf, ndf);
            for (int i = 0; i < nd; ++i) {
                // skip a pressure DOF (3D U-P node: node_ndf==4, i==3), same
                // heuristic as ResponseSpectrumAnalysis. NOTE: the 2D U-P case
                // (node_ndf==3, pressure at i==2) is ambiguous with a 2D beam
                // node (Ux,Uy,Rz) by DOF count alone, so like RSA we do not
                // special-case it. Skipped/unmodeled DOFs are left at 0 by the
                // full-vector overwrite below (correct for a from-rest modal run).
                if (ndf == 6 && node_ndf == 4 && i == 3)
                    continue;
                double sd = 0.0, sv = 0.0, sa = 0.0;
                for (int a = 0; a < nm; ++a) {
                    const double psi = evec(i, modes[a]) * Vscale[a];
                    sd += psi * q[a];
                    sv += psi * qd[a];
                    sa += psi * qdd[a];
                }
                ud(i) = sd; uv(i) = sv; ua(i) = sa;
            }
            node->setTrialDisp(ud);
            node->setTrialVel(uv);
            node->setTrialAccel(ua);
        }

        // advance domain clock, propagate to elements, commit + record
        m_model->setCurrentDomainTime(t);
        if (m_model->updateDomain() < 0) {
            opserr << "LadrunoModalResponse - updateDomain failed at step "
                   << n << ".\n";
            return -1;
        }
        if (m_model->commitDomain() < 0) {
            opserr << "LadrunoModalResponse - commitDomain failed at step "
                   << n << ".\n";
            return -1;
        }

        // advance the modal SDOFs from station n to n+1
        if (n < m_nsteps) {
            const double t_next = m_t0 + (n + 1) * m_dt;
            const double ug_next = m_baseAccel->getFactor(t_next);
            for (int a = 0; a < nm; ++a) {
                const double f0 = -Gamma[a] * ug;
                const double f1 = -Gamma[a] * ug_next;
                const double q0 = q[a], qd0 = qd[a];
                q[a]  = Amat[a * 4 + 0] * q0 + Amat[a * 4 + 1] * qd0
                      + Bmat[a * 4 + 0] * f0 + Bmat[a * 4 + 1] * f1;
                qd[a] = Amat[a * 4 + 2] * q0 + Amat[a * 4 + 3] * qd0
                      + Bmat[a * 4 + 2] * f0 + Bmat[a * 4 + 3] * f1;
            }
        }
    }

    return 0;
}

// ---------------------------------------------------------------------------
// Command entry point (openseespy + Tcl share this via the interpreter
// wrappers, mirroring OPS_ResponseSpectrumAnalysis).
//
//   modalResponseHistory -dt $dt -nsteps $n -baseAccel $tsTag -dir $dir
//       (-damp $xi | -rayleigh $a0 $a1 | -modalDamp $xi1 $xi2 ...)
//       [-modes $m1 $m2 ...] [-t0 $t0]
//
// Requires a prior `eigen N` and `modalProperties`.
// ---------------------------------------------------------------------------
int
OPS_LadrunoModalResponseHistory(void)
{
    AnalysisModel* theModel = *OPS_GetAnalysisModel();
    if (theModel == 0 || theModel->getDomainPtr() == 0) {
        opserr << "modalResponseHistory - no AnalysisModel/Domain available.\n";
        return -1;
    }

    double dt = 0.0, t0 = 0.0;
    int nsteps = 0, dir = 0;
    TimeSeries* ts = 0;
    int dampKind = -1;              // 0 uniform, 1 rayleigh, 2 modal-list
    double xi = 0.0, a0 = 0.0, a1 = 0.0;
    std::vector<double> xilist;
    std::vector<int> modes;
    int numData = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();

        if (strcmp(opt, "-dt") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&numData, &dt) < 0) {
                opserr << "modalResponseHistory - -dt requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-nsteps") == 0 || strcmp(opt, "-numSteps") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&numData, &nsteps) < 0) {
                opserr << "modalResponseHistory - -nsteps requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-t0") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&numData, &t0) < 0) {
                opserr << "modalResponseHistory - -t0 requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-dir") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&numData, &dir) < 0) {
                opserr << "modalResponseHistory - -dir requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-baseAccel") == 0) {
            int tag;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&numData, &tag) < 0) {
                opserr << "modalResponseHistory - -baseAccel requires a timeSeries tag.\n"; return -1;
            }
            ts = OPS_getTimeSeries(tag);
            if (ts == 0) {
                opserr << "modalResponseHistory - no timeSeries with tag " << tag << ".\n"; return -1;
            }
        }
        else if (strcmp(opt, "-damp") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&numData, &xi) < 0) {
                opserr << "modalResponseHistory - -damp requires a value.\n"; return -1;
            }
            // Ladruno (ADR 44 review): refuse negative damping family-wide. A
            // negative ratio is un-physical (grows the modal response) — a typo
            // trap. xi==0 stays legal (undamped modal history).
            if (xi < 0.0) {
                opserr << "modalResponseHistory - -damp xi must be >= 0 (got "
                       << xi << ").\n"; return -1;
            }
            dampKind = 0;
        }
        else if (strcmp(opt, "-rayleigh") == 0) {
            double v[2]; int nd = 2;
            if (OPS_GetNumRemainingInputArgs() < 2 || OPS_GetDoubleInput(&nd, v) < 0) {
                opserr << "modalResponseHistory - -rayleigh requires a0 a1.\n"; return -1;
            }
            a0 = v[0]; a1 = v[1]; dampKind = 1;
        }
        else if (strcmp(opt, "-modalDamp") == 0 || strcmp(opt, "-modalDamping") == 0) {
            xilist.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double item; int one = 1;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetDoubleInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1); // unconsume the non-numeric token
                    break;
                }
                xilist.push_back(item);
            }
            if (xilist.empty()) {
                opserr << "modalResponseHistory - -modalDamp requires >=1 ratio.\n"; return -1;
            }
            // Ladruno (ADR 44 review): every per-mode ratio must be >= 0.
            for (size_t j = 0; j < xilist.size(); ++j) {
                if (xilist[j] < 0.0) {
                    opserr << "modalResponseHistory - -modalDamp entry " << (int)j
                           << " = " << xilist[j] << " must be >= 0.\n"; return -1;
                }
            }
            dampKind = 2;
        }
        else if (strcmp(opt, "-modes") == 0) {
            modes.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int item; int one = 1;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                modes.push_back(item);
            }
            if (modes.empty()) {
                opserr << "modalResponseHistory - -modes requires >=1 mode index.\n"; return -1;
            }
        }
        else {
            opserr << "modalResponseHistory - unknown option '" << opt << "'.\n";
            return -1;
        }
    }

    // required arguments
    if (ts == 0) {
        opserr << "modalResponseHistory - -baseAccel <tsTag> is required (P1a).\n"; return -1;
    }
    if (dir < 1) {
        opserr << "modalResponseHistory - -dir <1..ndf> is required.\n"; return -1;
    }
    if (dt <= 0.0 || nsteps < 1) {
        opserr << "modalResponseHistory - -dt >0 and -nsteps >=1 are required.\n"; return -1;
    }
    if (dampKind < 0) {
        opserr << "modalResponseHistory - a damping channel is required "
                  "(-damp | -rayleigh | -modalDamp).\n"; return -1;
    }

    LadrunoModalResponse mrh(theModel, ts, dir, dt, nsteps, t0);
    if (dampKind == 0)      mrh.setDampingUniform(xi);
    else if (dampKind == 1) mrh.setDampingRayleigh(a0, a1);
    else                    mrh.setDampingModalList(xilist);
    if (!modes.empty())     mrh.setModeSubset(modes);

    return mrh.analyze();
}

// ===========================================================================
// ADR 44 P2: modal frequency-response / steady-state dynamics
// ===========================================================================
LadrunoFrequencyResponse::LadrunoFrequencyResponse(AnalysisModel* theModel)
    : m_model(theModel)
    , m_direction(0)
    , m_amp(1.0)
    , m_fmin(0.0), m_fmax(0.0)
    , m_nf(0)
    , m_sweep(SWEEP_LIN)
    , m_nodeTag(0)
    , m_dof(0)
    , m_resp(RESP_DISP)
{
}

LadrunoFrequencyResponse::~LadrunoFrequencyResponse() {}

void LadrunoFrequencyResponse::setModeSubset(const std::vector<int>& modes_1based)
{
    m_modes.clear();
    for (int m : modes_1based) m_modes.push_back(m - 1); // to 0-based
}

// Build the Hz sample grid.  LIN: uniform in f.  LOG: geometric.  BIASED: a
// LIN base grid PLUS a fine cluster around each in-band modal frequency f_a =
// w_a/2pi so sharp (low-damping) resonant peaks are resolved.
int
LadrunoFrequencyResponse::buildGrid(std::vector<double>& freqs) const
{
    freqs.clear();
    if (m_nf < 1 || m_fmin < 0.0)
        return -1;

    if (m_nf == 1) { freqs.push_back(m_fmin); return 0; }

    // a multi-point sweep needs a real interval; fmax==fmin with nf>1 would emit
    // nf identical rows (LIN/LOG don't de-dup). Use nf==1 for a single frequency.
    if (m_fmax <= m_fmin)
        return -1;

    if (m_sweep == SWEEP_LOG) {
        if (m_fmin <= 0.0) return -1; // log needs a positive lower bound
        const double r = std::log(m_fmax / m_fmin) / (m_nf - 1);
        for (int i = 0; i < m_nf; ++i)
            freqs.push_back(m_fmin * std::exp(r * i));
        return 0;
    }

    // LIN base grid (also the base for BIASED)
    const double df = (m_fmax - m_fmin) / (m_nf - 1);
    for (int i = 0; i < m_nf; ++i)
        freqs.push_back(m_fmin + df * i);

    if (m_sweep == SWEEP_BIASED && m_model && m_model->getDomainPtr()) {
        const Vector& ev = m_model->getDomainPtr()->getEigenvalues();
        // active modes for the cluster (default: all extracted)
        std::vector<int> modes = m_modes;
        if (modes.empty())
            for (int i = 0; i < ev.Size(); ++i) modes.push_back(i);
        const int NCLUST = 15; // points per side within the window
        for (int m : modes) {
            if (m < 0 || m >= ev.Size()) continue;
            const double lam = ev(m);
            if (lam <= 0.0) continue; // rigid mode: no resonant peak to resolve
            const double fa = std::sqrt(lam) / (2.0 * M_PI);
            if (fa < m_fmin || fa > m_fmax) continue;
            const double halfw = 0.05 * fa; // +/-5% window around the peak
            for (int k = -NCLUST; k <= NCLUST; ++k) {
                const double f = fa + halfw * (double(k) / NCLUST);
                if (f >= m_fmin && f <= m_fmax) freqs.push_back(f);
            }
        }
        std::sort(freqs.begin(), freqs.end());
        freqs.erase(std::unique(freqs.begin(), freqs.end()), freqs.end());
    }
    return 0;
}

int
LadrunoFrequencyResponse::run(std::vector<std::vector<double> >& rows, bool ampOnly)
{
    rows.clear();
    if (m_model == 0 || m_model->getDomainPtr() == 0) {
        opserr << "frequencyResponse - no AnalysisModel/Domain.\n";
        return -1;
    }
    Domain* domain = m_model->getDomainPtr();

    const Vector& ev = domain->getEigenvalues();
    if (ev.Size() < 1) {
        opserr << "frequencyResponse - no eigenvalues; call 'eigen' first.\n";
        return -1;
    }
    DomainModalProperties mp;
    if (domain->getModalProperties(mp) < 0) {
        opserr << "frequencyResponse - eigen and modalProperties have not been "
                  "called.\n";
        return -1;
    }
    if (ev.Size() != mp.eigenvalues().Size()) {
        opserr << "frequencyResponse - modalProperties is stale vs the last eigen "
                  "run; call 'modalProperties' right after 'eigen'.\n";
        return -1;
    }
    // element-wise staleness check (same as P1a): a re-run of `eigen` with the
    // SAME mode count but a changed model leaves a same-size DomainModalProperties
    // that the size test alone waves through — the FRF would then be built from
    // stale w/Gamma/Vscale while the -biased grid reads the NEW domain eigenvalues
    // (inconsistent, silent). Refuse unless the snapshot matches the domain.
    {
        const double tol = std::max(1.0e-15, 1.0e-12 * ev.Norm());
        for (int i = 0; i < ev.Size(); ++i)
            if (std::abs(ev(i) - mp.eigenvalues()(i)) > tol) {
                opserr << "frequencyResponse - eigenvalue mismatch between the domain "
                          "and DomainModalProperties (re-run modalProperties after "
                          "eigen).\n";
                return -1;
            }
    }
    const int num_eigen = ev.Size();
    const int ndf = mp.totalMass().Size();
    if (m_direction < 1 || m_direction > ndf) {
        opserr << "frequencyResponse - -dir (" << m_direction << ") must be in 1.."
               << ndf << ".\n";
        return -1;
    }
    const int exdof = m_direction - 1;

    // response node/dof
    Node* rnode = domain->getNode(m_nodeTag);
    if (rnode == 0) {
        opserr << "frequencyResponse - response node " << m_nodeTag
               << " not found.\n";
        return -1;
    }
    // probe with the non-exiting accessor first: Node::getEigenvectors() exit(0)s
    // the whole process when unset (a fully-constrained node, or one outside the
    // eigen analysis), which for a user-chosen response node would be a silent kill.
    if (rnode->getNumEigenvectors() < 1) {
        opserr << "frequencyResponse - response node " << m_nodeTag
               << " has no eigenvector (fully constrained, or not in the last "
                  "eigen analysis); pick a free node.\n";
        return -1;
    }
    const Matrix& revec = rnode->getEigenvectors();
    const int rndf = revec.noRows();
    if (m_dof < 1 || m_dof > rndf) {
        opserr << "frequencyResponse - -dof (" << m_dof << ") must be in 1.." << rndf
               << " for node " << m_nodeTag << ".\n";
        return -1;
    }
    const int rdof = m_dof - 1;

    // active mode list (default: all), 0-based, no duplicates
    std::vector<int> modes = m_modes;
    if (modes.empty()) {
        modes.resize(num_eigen);
        for (int i = 0; i < num_eigen; ++i) modes[i] = i;
    }
    for (int m : modes)
        if (m < 0 || m >= num_eigen) {
            opserr << "frequencyResponse - requested mode " << (m + 1)
                   << " is out of range 1.." << num_eigen << ".\n";
            return -1;
        }
    {
        std::vector<int> seen = modes;
        std::sort(seen.begin(), seen.end());
        if (std::adjacent_find(seen.begin(), seen.end()) != seen.end()) {
            opserr << "frequencyResponse - -modes contains a duplicate index.\n";
            return -1;
        }
    }
    if (m_damp.kind == LadrunoModalDamping::MODAL_LIST) {
        int maxmode = 0;
        for (int m : modes) maxmode = std::max(maxmode, m);
        if (static_cast<int>(m_damp.xis.size()) <= maxmode) {
            opserr << "frequencyResponse - -modalDamp list has "
                   << static_cast<int>(m_damp.xis.size()) << " ratios but mode "
                   << (maxmode + 1) << " is used.\n";
            return -1;
        }
    }
    const int nm = static_cast<int>(modes.size());

    // negative-eigenvalue guard (same policy as P1a)
    double maxAbsLambda = 0.0;
    for (int m : modes)
        maxAbsLambda = std::max(maxAbsLambda, std::abs(mp.eigenvalues()(m)));
    const double negTol = std::max(1.0e-12, 1.0e-8 * maxAbsLambda);

    // per-mode precompute: w, d, and the scalar recovery weight
    //   c_a = psi_a(node,dof) * (-Gamma_a),   psi_a = revec(rdof,mode)*Vscale
    std::vector<double> w(nm), d(nm), cw(nm);
    for (int a = 0; a < nm; ++a) {
        const int mode0 = modes[a];
        const double lambda = mp.eigenvalues()(mode0);
        if (lambda < -negTol) {
            opserr << "frequencyResponse - mode " << (mode0 + 1)
                   << " has a negative eigenvalue (" << lambda << ").\n";
            return -1;
        }
        const double wa = (lambda > 0.0) ? std::sqrt(lambda) : 0.0;
        w[a] = wa;
        d[a] = m_damp.coeff(mode0, wa);
        const double Gamma  = mp.modalParticipationFactors()(mode0, exdof);
        const double Vscale = mp.eigenVectorScaleFactors()(mode0);
        const double psi    = revec(rdof, mode0) * Vscale;
        cw[a] = psi * (-Gamma);
    }

    std::vector<double> freqs;
    if (buildGrid(freqs) < 0) {
        opserr << "frequencyResponse - bad frequency sweep "
                  "(need nf>=1, 0<=fmin<=fmax; log needs fmin>0).\n";
        return -1;
    }

    // sweep
    const std::complex<double> I(0.0, 1.0);
    for (double f : freqs) {
        const double Om = 2.0 * M_PI * f;
        std::complex<double> uc(0.0, 0.0);
        for (int a = 0; a < nm; ++a) {
            const std::complex<double> denom(w[a] * w[a] - Om * Om, Om * d[a]);
            uc += cw[a] / denom;               // psi*(-Gamma)*H_a(Om)
        }
        // response-quantity multiplier (relative), then excitation amplitude
        if (m_resp == RESP_VEL)        uc *= I * Om;
        else if (m_resp == RESP_ACCEL) uc *= -(Om * Om);
        uc *= m_amp;

        std::vector<double> row;
        row.push_back(f);
        if (ampOnly) {
            row.push_back(std::abs(uc));
        } else {
            row.push_back(uc.real());
            row.push_back(uc.imag());
        }
        rows.push_back(row);
    }

    // optional file output
    if (!m_outfile.empty()) {
        std::ofstream os(m_outfile.c_str());
        if (!os) {
            opserr << "frequencyResponse - cannot open output file '"
                   << m_outfile.c_str() << "'.\n";
            return -1;
        }
        os << std::scientific << std::setprecision(15);
        for (const std::vector<double>& r : rows) {
            for (size_t j = 0; j < r.size(); ++j)
                os << (j ? " " : "") << r[j];
            os << "\n";
        }
    }
    return 0;
}

// ---------------------------------------------------------------------------
// Shared parser for frequencyResponse / steadyStateDynamics.
//
//   frequencyResponse -freq fmin fmax nf [-lin|-log|-biased]
//       -baseAccel -dir $dir [-amp $a]
//       (-damp $xi | -rayleigh $a0 $a1 | -modalDamp $xi1 ...)
//       -node $tag -dof $dof [-resp disp|vel|accel] [-modes ...] [-out $file]
//
//   steadyStateDynamics ...same... (reports |response| instead of Re/Im)
//
// Requires a prior `eigen N` and `modalProperties`.  Returns the table to the
// interpreter via OPS_SetDoubleListsOutput (and writes -out if given).
// ---------------------------------------------------------------------------
static int
OPS_FreqResponseImpl(bool ampOnly, const char* cmd)
{
    AnalysisModel* theModel = *OPS_GetAnalysisModel();
    if (theModel == 0 || theModel->getDomainPtr() == 0) {
        opserr << cmd << " - no AnalysisModel/Domain available.\n";
        return -1;
    }

    double fmin = 0.0, fmax = 0.0, amp = 1.0;
    int nf = 0, dir = 0, nodeTag = 0, dof = 0;
    LadrunoFrequencyResponse::SweepKind sweep = LadrunoFrequencyResponse::SWEEP_LIN;
    LadrunoFrequencyResponse::RespKind  resp  = LadrunoFrequencyResponse::RESP_DISP;
    bool haveFreq = false, haveBase = false, haveDamp = false, haveNode = false;
    LadrunoModalDamping damp;
    std::vector<int> modes;
    std::string outfile;
    int one = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();

        if (strcmp(opt, "-freq") == 0) {
            double fv[2]; int nd = 2;
            if (OPS_GetNumRemainingInputArgs() < 3 || OPS_GetDoubleInput(&nd, fv) < 0
                    || OPS_GetIntInput(&one, &nf) < 0) {
                opserr << cmd << " - -freq requires fmin fmax nf.\n"; return -1;
            }
            fmin = fv[0]; fmax = fv[1]; haveFreq = true;
        }
        else if (strcmp(opt, "-lin") == 0)    sweep = LadrunoFrequencyResponse::SWEEP_LIN;
        else if (strcmp(opt, "-log") == 0)    sweep = LadrunoFrequencyResponse::SWEEP_LOG;
        else if (strcmp(opt, "-biased") == 0) sweep = LadrunoFrequencyResponse::SWEEP_BIASED;
        else if (strcmp(opt, "-baseAccel") == 0) haveBase = true;
        else if (strcmp(opt, "-dir") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &dir) < 0) {
                opserr << cmd << " - -dir requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-amp") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&one, &amp) < 0) {
                opserr << cmd << " - -amp requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-node") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &nodeTag) < 0) {
                opserr << cmd << " - -node requires a tag.\n"; return -1;
            }
            haveNode = true;
        }
        else if (strcmp(opt, "-dof") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &dof) < 0) {
                opserr << cmd << " - -dof requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-resp") == 0 || strcmp(opt, "-response") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << cmd << " - -resp requires disp|vel|accel.\n"; return -1;
            }
            const char* rname = OPS_GetString();
            if (strcmp(rname, "disp") == 0)       resp = LadrunoFrequencyResponse::RESP_DISP;
            else if (strcmp(rname, "vel") == 0)   resp = LadrunoFrequencyResponse::RESP_VEL;
            else if (strcmp(rname, "accel") == 0) resp = LadrunoFrequencyResponse::RESP_ACCEL;
            else {
                opserr << cmd << " - unknown -resp '" << rname
                       << "' (use disp|vel|accel).\n"; return -1;
            }
        }
        else if (strcmp(opt, "-damp") == 0) {
            double xi;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&one, &xi) < 0) {
                opserr << cmd << " - -damp requires a value.\n"; return -1;
            }
            // Ladruno (ADR 44 review): refuse negative damping. d_a=2 xi w_a<0
            // silently CONJUGATES the FRF (H -> H*) — identical |H|, only the phase
            // sign flips — so a typo'd '-0.05' passes every magnitude test. xi==0 is
            // legal (the honest undamped singularity, documented in the guide).
            if (xi < 0.0) {
                opserr << cmd << " - -damp xi must be >= 0 (got " << xi
                       << "); negative damping conjugates the FRF (same magnitude, "
                          "flipped phase). Use xi=0 for the undamped case.\n";
                return -1;
            }
            damp.kind = LadrunoModalDamping::UNIFORM; damp.xi = xi; haveDamp = true;
        }
        else if (strcmp(opt, "-rayleigh") == 0) {
            double v[2]; int nd = 2;
            if (OPS_GetNumRemainingInputArgs() < 2 || OPS_GetDoubleInput(&nd, v) < 0) {
                opserr << cmd << " - -rayleigh requires a0 a1.\n"; return -1;
            }
            damp.kind = LadrunoModalDamping::RAYLEIGH; damp.a0 = v[0]; damp.a1 = v[1];
            haveDamp = true;
        }
        else if (strcmp(opt, "-modalDamp") == 0 || strcmp(opt, "-modalDamping") == 0) {
            damp.xis.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double item;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetDoubleInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                damp.xis.push_back(item);
            }
            if (damp.xis.empty()) {
                opserr << cmd << " - -modalDamp requires >=1 ratio.\n"; return -1;
            }
            // Ladruno (ADR 44 review): every per-mode ratio must be >= 0 (same
            // negative-damping FRF-conjugation trap as -damp).
            for (size_t j = 0; j < damp.xis.size(); ++j) {
                if (damp.xis[j] < 0.0) {
                    opserr << cmd << " - -modalDamp entry " << (int)j << " = "
                           << damp.xis[j] << " must be >= 0 (negative damping "
                              "conjugates the FRF).\n";
                    return -1;
                }
            }
            damp.kind = LadrunoModalDamping::MODAL_LIST; haveDamp = true;
        }
        else if (strcmp(opt, "-modes") == 0) {
            modes.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int item;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                modes.push_back(item);
            }
            if (modes.empty()) {
                opserr << cmd << " - -modes requires >=1 mode index.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-out") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << cmd << " - -out requires a filename.\n"; return -1;
            }
            outfile = OPS_GetString();
        }
        else {
            opserr << cmd << " - unknown option '" << opt << "'.\n"; return -1;
        }
    }

    if (!haveFreq) { opserr << cmd << " - -freq fmin fmax nf is required.\n"; return -1; }
    if (!haveBase || dir < 1) {
        opserr << cmd << " - -baseAccel -dir <1..ndf> is required (P2).\n"; return -1;
    }
    if (!haveNode || dof < 1) {
        opserr << cmd << " - -node <tag> -dof <1..ndf> are required.\n"; return -1;
    }
    if (!haveDamp) {
        opserr << cmd << " - a damping channel is required "
                  "(-damp | -rayleigh | -modalDamp).\n"; return -1;
    }

    // Ladruno (ADR 44 review): a -biased grid places k=0 EXACTLY on each in-band
    // resonance f_a=w_a/2pi (f = fa + halfw*(k/NCLUST), k=0 => f=fa). With an
    // exactly-zero d_a the denominator w_a^2-Om^2+i Om d_a is real and hits 0
    // there => an inf/NaN FRF row. That is an honest singularity (not guarded),
    // but -biased manufactures it deliberately, so warn once — the user almost
    // never wants a zero-damping resonance-biased sweep.
    if (sweep == LadrunoFrequencyResponse::SWEEP_BIASED) {
        bool zeroDamp = false;
        if (damp.kind == LadrunoModalDamping::UNIFORM)       zeroDamp = (damp.xi == 0.0);
        else if (damp.kind == LadrunoModalDamping::RAYLEIGH) zeroDamp = (damp.a0 == 0.0 && damp.a1 == 0.0);
        else if (damp.kind == LadrunoModalDamping::MODAL_LIST) {
            for (size_t j = 0; j < damp.xis.size(); ++j)
                if (damp.xis[j] == 0.0) { zeroDamp = true; break; }
        }
        if (zeroDamp)
            opserr << "WARNING " << cmd << " - -biased with zero damping samples "
                      "EXACTLY on an undamped resonance (f_a=w_a/2pi) => inf/NaN "
                      "FRF row at that mode. Use a positive damping ratio or a "
                      "-lin/-log grid.\n";
    }

    LadrunoFrequencyResponse fr(theModel);
    fr.setBaseAccel(dir, amp);
    fr.setSweep(fmin, fmax, nf, sweep);
    fr.setDamping(damp);
    fr.setResponse(nodeTag, dof, resp);
    if (!modes.empty())  fr.setModeSubset(modes);
    if (!outfile.empty()) fr.setOutputFile(outfile.c_str());

    std::vector<std::vector<double> > rows;
    if (fr.run(rows, ampOnly) < 0)
        return -1;

    OPS_SetDoubleListsOutput(rows);
    return 0;
}

int OPS_LadrunoFrequencyResponse(void)  { return OPS_FreqResponseImpl(false, "frequencyResponse"); }
int OPS_LadrunoSteadyStateDynamics(void){ return OPS_FreqResponseImpl(true,  "steadyStateDynamics"); }
