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

// Ladruno: shared coupling-kernel for the embedded-element family. See
// LadrunoEmbeddedKernel.h and Ladruno_implementation/23_ladruno_embedded_node_adr.md.
// Every function here is a verbatim extraction of the corresponding numeric core
// from LadrunoEmbeddedRebar (ELE 33005) — same operation order ⇒ bit-identical.
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoEmbeddedKernel.h>
#include <Matrix.h>
#include <Vector.h>
#include <Node.h>
#include <math.h>

namespace LadrunoEmbedded {

double maxAbsDiagonal(const Matrix& M)
{
  double scale = 0.0;
  int n = M.noRows();
  for (int i = 0; i < n; i++) {
    double d = fabs(M(i, i));
    if (d > scale) scale = d;
  }
  return scale;
}

double effectiveCouplingStiffness(double kAxial0, double kt)
{
  return (kAxial0 > kt) ? kAxial0 : kt;
}

double massPenaltyDtcr(double kEff, double dt)
{
  return kEff * 0.25 * dt * dt;   // k_eff·(dt/2)²
}

double massPenaltyWcap(double kEff, double beta, double omega2Host)
{
  return kEff / (beta * beta * omega2Host);
}

double criticalTimeStep(double mPenalty, double kEff)
{
  if (mPenalty <= 0.0 || kEff <= 0.0) return -1.0;
  return 2.0 * sqrt(mPenalty / kEff);
}

void computeGap(Node** nodes, const Vector& N, int nHost, int ndm, Vector& g)
{
  g.resize(ndm);
  g.Zero();
  const Vector& uR = nodes[0]->getTrialDisp();
  for (int k = 0; k < ndm; k++) g(k) = uR(k);
  for (int i = 0; i < nHost; i++) {
    const Vector& uH = nodes[1 + i]->getTrialDisp();
    for (int k = 0; k < ndm; k++) g(k) -= N(i) * uH(k);
  }
}

void assembleResistingForce(Vector& P, const Vector& N, const Vector& t, int ndm)
{
  int nHost = N.Size();
  // block 0 (constrained): +t
  for (int k = 0; k < ndm; k++) P(k) = t(k);
  // block (1+i) (host i): -N_i t
  for (int i = 0; i < nHost; i++) {
    double c = -N(i);
    for (int k = 0; k < ndm; k++) P((1 + i) * ndm + k) = c * t(k);
  }
}

void assembleTangent(Matrix& K, const Vector& N, const Matrix& D, int ndm)
{
  int nHost = N.Size();
  // c = [1, -N_1, ..., -N_M]; K block(p,q) = c_p c_q D
  for (int p = 0; p < 1 + nHost; p++) {
    double cp = (p == 0) ? 1.0 : -N(p - 1);
    for (int q = 0; q < 1 + nHost; q++) {
      double cq = (q == 0) ? 1.0 : -N(q - 1);
      double cc = cp * cq;
      if (cc == 0.0) continue;
      for (int a = 0; a < ndm; a++)
        for (int b = 0; b < ndm; b++)
          K(p * ndm + a, q * ndm + b) += cc * D(a, b);
    }
  }
}

} // namespace LadrunoEmbedded
