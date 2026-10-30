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

// penalty normal traction tn = kn <-gap>+ (Macaulay bracket): zero for gap >= 0,
// so the contact never pulls.
static inline double
ladrunoContactTraction(double kn, double gap)
{
    return (gap < 0.0) ? kn * (-gap) : 0.0;
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes,
                                   int nps_, double kn_, const double odir[3],
                                   Domain *dom, int contactTag_, int segIndex_,
                                   bool consistentNormal_)
  : FE_Element(tag, /*numDOF_Group=*/1 + nps_, /*ndof=*/3 * (1 + nps_)),
    resid(3 * (1 + nps_)), tang(3 * (1 + nps_), 3 * (1 + nps_)),
    mode(SEGMENT), theSlave(slaveNode), ndm(3), kn(kn_), nps(nps_),
    theDomain(dom), contactTag(contactTag_), segIndex(segIndex_),
    consistentNormal(consistentNormal_)
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
LadrunoContactFE::segmentActive(double &gap, double n[3], double N[4], double *B) const
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
    return true;
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
        double gap, n[3], N[4], B[15];             // ndof <= 3 (1 + 4) = 15
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (!segmentActive(gap, n, N, B)) {
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
        double gap, n[3], N[4], B[15];
        if (segmentActive(gap, n, N, B)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            if (consistentNormal)
                addNormalGeomTang(fact);
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
