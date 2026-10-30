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
// Description: LadrunoContactFE, see LadrunoContactFE.h.

#include "LadrunoContactFE.h"
#include <cmath>
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <Domain.h>
#include <LadrunoContactDomain.h>
#include <LadrunoContactProjection.h>
#include <LadrunoFrictionKernel.h>

// penalty normal traction tn = kn <-gap>+ (Macaulay bracket): zero for gap >= 0,
// so the contact never pulls.
static inline double
ladrunoContactTraction(double kn, double gap)
{
    return (gap < 0.0) ? kn * (-gap) : 0.0;
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes,
                                   int nps_, double kn_, const double odir[3],
                                   double kt_, double mu_, Domain *dom,
                                   int contactTag_, int segIndex_, bool consistentTan_,
                                   bool consistentNormal_, double cohesion_, double tauMax_)
  : FE_Element(tag, /*numDOF_Group=*/1 + nps_, /*ndof=*/3 * (1 + nps_)),
    resid(3 * (1 + nps_)), tang(3 * (1 + nps_), 3 * (1 + nps_)),
    mode(SEGMENT), theSlave(slaveNode), ndm(3), kn(kn_), nps(nps_),
    theDomain(dom), contactTag(contactTag_), segIndex(segIndex_),
    consistentNormal(consistentNormal_),
    kt(kt_), mu(mu_), cohesion(cohesion_), tauMax(tauMax_), consistentTan(consistentTan_)
{
    // Connectivity = slave DOF_Group + each segment-node DOF_Group; setID() fills
    // myID = [slave xyz | seg_1 xyz | ... | seg_nps xyz], the layout B assumes.
    if (slaveNode != 0) {
        DOF_Group *dg = slaveNode->getDOF_GroupPtr();
        if (dg != 0) myDOF_Groups(0) = dg->getTag();
    }
    for (int i = 0; i < nps; i++) {
        segNode[i] = segNodes[i];
        if (segNodes[i] != 0) {
            DOF_Group *dg = segNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(1 + i) = dg->getTag();
        }
    }
    for (int i = nps; i < 4; i++) segNode[i] = 0;
    // a fixed direction, not the live slave position, so the orientation stays
    // correct after the slave penetrates
    for (int d = 0; d < 3; d++) orientDir[d] = odir[d];
}

LadrunoContactFE::~LadrunoContactFE()
{
}

double
LadrunoContactFE::explicitTimeStep(Integrator *)
{
    return 0.0;   // no integrator in this build reports its explicit time step
}

int
LadrunoContactFE::setID(void)
{
    Node *nodes[8]; int nn = 0;
    switch (mode) {
    case SEGMENT:
        nodes[nn++] = theSlave;
        for (int i = 0; i < nps; i++) nodes[nn++] = segNode[i];
        break;
    default:
        return FE_Element::setID();
    }
    const int per  = ndm;
    const int ndof = myID.Size();
    if (nn != myDOF_Groups.Size() || per * nn != ndof) {
        opserr << "FATAL LadrunoContactFE::setID() - layout mismatch: " << nn << " nodes x "
               << per << " != ndof " << ndof << " (" << myDOF_Groups.Size() << " DOF_Groups)\n";
        return -3;
    }
    int current = 0;
    for (int i = 0; i < nn; i++) {
        DOF_Group *dg = (nodes[i] != 0) ? nodes[i]->getDOF_GroupPtr() : 0;
        if (dg == 0) {
            opserr << "WARNING LadrunoContactFE::setID() - 0 DOF_Group pointer\n";
            return -2;
        }
        const ID &eq = dg->getID();
        if (eq.Size() < per) {
            opserr << "FATAL LadrunoContactFE::setID() - node " << nodes[i]->getTag()
                   << " has " << eq.Size() << " equations < ndm=" << per << "\n";
            return -3;
        }
        for (int j = 0; j < per; j++) myID(current++) = eq(j);
    }
    return 0;
}

bool
LadrunoContactFE::segmentActive(double &gap, double n[3], double N[4], double *B,
                                double *gTvec) const
{
    if (mode != SEGMENT || theSlave == 0) return false;
    double Xseg[4][3], xs[3];
    const Vector &Xs = theSlave->getCrds();
    const Vector &us = theSlave->getTrialDisp();
    for (int d = 0; d < 3; d++) xs[d] = Xs(d) + us(d);
    for (int i = 0; i < nps; i++) {
        if (segNode[i] == 0) return false;
        const Vector &Xi = segNode[i]->getCrds();
        const Vector &ui = segNode[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xseg[i][d] = Xi(d) + ui(d);
    }
    if (!LadrunoContactProjection::evalSegment(nps, Xseg, xs, orientDir, gap, n, N))
        return false;
    // gap operator B (1 x ndof) over [u_s | u_1..u_nps]: [ n^T | -N_i n^T ]
    int ndof = 3 * (1 + nps);
    for (int k = 0; k < ndof; k++) B[k] = 0.0;
    for (int d = 0; d < 3; d++) B[d] = n[d];
    for (int i = 0; i < nps; i++)
        for (int d = 0; d < 3; d++) B[3 * (1 + i) + d] = -N[i] * n[d];
    // Tangential relative slip at the contact point. The closest-point
    // projection makes (x_s - xbar) parallel to n, so positions carry no
    // tangential information; the slip is the slave displacement minus the
    // master displacement interpolated at the projection, d = u_s - sum N_i u_i,
    // tangential part (the caller subtracts the engagement origin gT0). Same
    // trial configuration and same shape weights as the gap.
    if (gTvec != 0) {
        double ubar[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < nps; i++) {
            const Vector &ui = segNode[i]->getTrialDisp();
            for (int d = 0; d < 3; d++) ubar[d] += N[i] * ui(d);
        }
        double drel[3] = { us(0)-ubar[0], us(1)-ubar[1], us(2)-ubar[2] };
        LadrunoFrictionKernel::tangentPart(drel, n, gTvec);
    }
    return true;
}

void
LadrunoContactFE::addFrictionTang(double fact, const double n[3], const double N[4],
                                  double tn, const double gTeff[3], const double gpT[3],
                                  bool consistent)
{
    // K_fric = G^T K_ss G, G = [I | -N_i I]: block (a,b) of the ndof x ndof
    // tangent is w_a w_b K_ss, w = [1, -N_0, ..., -N_{nps-1}] over [slave, segment nodes]
    double Kss[3][3];
    LadrunoFrictionKernel::frictionTangentBlock(gTeff, gpT, n, tn, kn, kt, mu, consistent, Kss,
                                                cohesion, tauMax);
    double w[5];
    w[0] = 1.0;
    for (int i = 0; i < nps; i++) w[1 + i] = -N[i];
    int nn = 1 + nps;
    for (int a = 0; a < nn; a++)
        for (int b = 0; b < nn; b++) {
            double wab = fact * w[a] * w[b];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    tang(3 * a + i, 3 * b + j) += wab * Kss[i][j];
        }
}

void
LadrunoContactFE::addNormalGeomTang(double fact)
{
    if (mode != SEGMENT || theSlave == 0) return;
    double Xseg[4][3], xs[3];
    const Vector &Xs = theSlave->getCrds();
    const Vector &us = theSlave->getTrialDisp();
    for (int d = 0; d < 3; d++) xs[d] = Xs(d) + us(d);
    for (int i = 0; i < nps; i++) {
        if (segNode[i] == 0) return;
        const Vector &Xi = segNode[i]->getCrds();
        const Vector &ui = segNode[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xseg[i][d] = Xi(d) + ui(d);
    }
    double Kg[15][15];
    if (!LadrunoContactProjection::normalGeomTangent(nps, Xseg, xs, orientDir, kn, Kg))
        return;                                 // not penetrating or out of bounds
    int ndof = 3 * (1 + nps);
    for (int a = 0; a < ndof; a++)
        for (int b = 0; b < ndof; b++)
            tang(a, b) += fact * Kg[a][b];
}

const Vector &
LadrunoContactFE::getResidual(Integrator *theIntegrator)
{
    resid.Zero();
    if (mode == SEGMENT) {
        double gap, n[3], N[4], B[15], gTvec[3];   // ndof <= 3 (1 + 4) = 15
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (!segmentActive(gap, n, N, B, gTvec)) {
            // a released pair reports zero force; without this the readout kept
            // the last active value until the next handle()
            if (cd != 0 && theSlave != 0)
                cd->setNtsForce(contactTag, theSlave->getTag(), segIndex, 0.0);
        } else {
            double tn = ladrunoContactTraction(kn, gap);   // kn <-gap>+ > 0
            int ndof = 3 * (1 + nps);
            for (int k = 0; k < ndof; k++)
                resid(k) = B[k] * tn;        // r = B^T tn (slave +tn n, master -N_i tn n)

            // normal force readout (ladrunoContactForce); no effect on resid/tang
            if (cd != 0)
                cd->setNtsForce(contactTag, theSlave->getTag(), segIndex, tn);

            // --- Coulomb/Tresca friction (the tangent is in addKtToTang) ---
            // A frictionless pair returns before any state is touched.
            if (hasFriction() && cd != 0) {
                LadrunoContactDomain::FrictionState &st =
                    cd->getOrCreateFrictionState(contactTag, theSlave->getTag(), segIndex);
                // capture the tangential origin once, at the first activation, so
                // the tangential drift of a slave before contact does not become a
                // stick traction
                if (!st.engaged) {
                    for (int d = 0; d < 3; d++) st.gT0[d] = gTvec[d];
                    st.engaged = true;
                }
                double gTeff[3];
                for (int d = 0; d < 3; d++) gTeff[d] = gTvec[d] - st.gT0[d];
                double tFric[3], gpTtrial[3];
                // the cone pressure is the current normal force tn
                LadrunoFrictionKernel::frictionReturnMap(gTeff, st.gpT, tn, kt, mu,
                                                         tFric, gpTtrial, cohesion, tauMax);
                // the trial is a function of the committed state only; commit() promotes it
                for (int d = 0; d < 3; d++) st.gpTtrial[d] = gpTtrial[d];
                // as the normal block: slave += tFric, master_i += -N_i tFric
                // (tFric is the force on the slave, opposing its motion)
                for (int d = 0; d < 3; d++) resid(d) += tFric[d];
                for (int i = 0; i < nps; i++)
                    for (int d = 0; d < 3; d++)
                        resid(3 * (1 + i) + d) += -N[i] * tFric[d];
            }
        }
    }
    return resid;
}

const Matrix &
LadrunoContactFE::getTangent(Integrator *theIntegrator)
{
    // The integrator decides the K/C/M combination: an explicit integrator that
    // forms only M (and C) leaves the contact stiffness out of its matrix;
    // Newmark calls addKtToTang(c1); a static integrator addKtToTang(1).
    tang.Zero();
    if (theIntegrator != 0)
        theIntegrator->formEleTangent(this);
    return tang;
}

void
LadrunoContactFE::zeroTangent(void)
{
    tang.Zero();
}

void
LadrunoContactFE::addKtToTang(double fact)
{
    if (mode == SEGMENT) {
        // K_c = kn B^T B, exact for a flat fixed master where n is constant;
        // with consistentNormal the geometric block kn gN d^2gN/du^2 is added
        // (zero for a flat facet).
        double gap, n[3], N[4], B[15], gTvec[3];
        if (segmentActive(gap, n, N, B, gTvec)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            if (consistentNormal)
                addNormalGeomTang(fact);
            // friction tangent, from the committed gpT so that it is the
            // derivative of the residual evaluated at the same state
            if (hasFriction() && theDomain != 0) {
                LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
                if (cd != 0) {
                    LadrunoContactDomain::FrictionState &st =
                        cd->getOrCreateFrictionState(contactTag, theSlave->getTag(), segIndex);
                    double gTeff[3];
                    for (int d = 0; d < 3; d++)
                        gTeff[d] = st.engaged ? (gTvec[d] - st.gT0[d]) : 0.0;
                    double tn = ladrunoContactTraction(kn, gap);
                    addFrictionTang(fact, n, N, tn, gTeff, st.gpT, consistentTan);
                }
            }
        }
    }
}

void
LadrunoContactFE::addKiToTang(double fact)
{
    // Initial-stiffness algorithms (Newton -initial, ModifiedNewton -initial)
    // form the matrix from K_initial; the penalty K_c is used so the contact
    // stiffness is not dropped.
    if (mode == SEGMENT) {
        double gap, n[3], N[4], B[15];
        if (segmentActive(gap, n, N, B)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            // friction: the stick tangent kt G^T P_t G (SPD, so a modified or
            // initial-stiffness Newton iteration contracts); gTeff = gpT = 0
            // forces the stick branch, so no engine state is needed
            if (hasFriction()) {
                double tn = ladrunoContactTraction(kn, gap);
                double zero[3] = {0.0, 0.0, 0.0};
                addFrictionTang(fact, n, N, tn, zero, zero, false);
            }
        }
    }
}

void
LadrunoContactFE::addCtoTang(double)
{
}

void
LadrunoContactFE::addMtoTang(double)
{
}

const Vector &LadrunoContactFE::getTangForce(const Vector &, double) { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getK_Force(const Vector &, double)   { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getKi_Force(const Vector &, double)  { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getC_Force(const Vector &, double)   { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getM_Force(const Vector &, double)   { resid.Zero(); return resid; }
