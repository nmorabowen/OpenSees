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

// Written: nmb (UANDES), 06/2026
// See ExplicitBatheSMSConsistent.h and
// Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md.

#include <ExplicitBatheSMSConsistent.h>
#include <LadrunoMassScaling.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <LinearSOE.h>
#include <DiagonalSOE.h>
#include <Channel.h>
#include <Vector.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

#include <string.h>

void *OPS_ExplicitBatheSMSConsistent(void)
{
    // Usage: integrator ExplicitBatheSMSConsistent $p $dtTarget <-maxAddedMass f>
    //                   <-lump ...> <-tangent> <-pcgTol t> <-pcgMaxIt n> <-verbose>
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING integrator ExplicitBatheSMSConsistent $p $dtTarget <options> "
                  "- needs the Noh-Bathe p and a target time step\n";
        return 0;
    }

    double pin[2] = {0.0, 0.0};
    int n2 = 2;
    if (OPS_GetDoubleInput(&n2, pin) < 0) {
        opserr << "WARNING ExplicitBatheSMSConsistent - could not read $p $dtTarget\n";
        return 0;
    }
    double p = pin[0], dtTarget = pin[1];
    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING ExplicitBatheSMSConsistent - $p must be in (0,1)\n";
        return 0;
    }
    if (dtTarget <= 0.0) {
        opserr << "WARNING ExplicitBatheSMSConsistent - $dtTarget must be positive\n";
        return 0;
    }

    double maxAddedMassFrac = 0.05;
    bool   verboseSMS = false;
    bool   cflUseTangent = false;
    int    compute_critical_timestep = 1;
    CTSLumping lumping = CTSLumping::Diagonal;
    double pcgTol = 1.0e-10;
    int    pcgMaxIt = 200;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        if (strcmp(arg, "-verbose") == 0) {
            verboseSMS = true;
        } else if (strcmp(arg, "-maxAddedMass") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) { int n1 = 1; OPS_GetDoubleInput(&n1, &maxAddedMassFrac); }
        } else if (strcmp(arg, "-tangent") == 0) {
            cflUseTangent = true;
        } else if (strcmp(arg, "-cfl") == 0 || strcmp(arg, "-criticalTimestep") == 0) {
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-pcgTol") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) { int n1 = 1; OPS_GetDoubleInput(&n1, &pcgTol); }
        } else if (strcmp(arg, "-pcgMaxIt") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) { int n1 = 1; OPS_GetIntInput(&n1, &pcgMaxIt); }
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else if (strcmp(m, "hrz") == 0)      lumping = CTSLumping::HRZ;
                else opserr << "WARNING ExplicitBatheSMSConsistent - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping diagonal)\n";
            }
        } else if (strcmp(arg, "-cflAbort") == 0 || strcmp(arg, "-recompute") == 0) {
            opserr << "WARNING ExplicitBatheSMSConsistent - " << arg << " is not supported "
                      "with mass scaling. Integrator NOT created.\n";
            return 0;
        } else {
            opserr << "WARNING ExplicitBatheSMSConsistent - unknown option " << arg << " (ignored)\n";
        }
    }

    TransientIntegrator *theIntegrator =
        new ExplicitBatheSMSConsistent(p, dtTarget, maxAddedMassFrac, verboseSMS,
                                       compute_critical_timestep, cflUseTangent, lumping,
                                       pcgTol, pcgMaxIt);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitBatheSMSConsistent integrator\n";
    return theIntegrator;
}

// Default constructor (broker)
ExplicitBatheSMSConsistent::ExplicitBatheSMSConsistent()
    : ExplicitBathe(INTEGRATOR_TAGS_ExplicitBatheSMSConsistent, 0.0, 0, false, false, 0.0,
                    false, 0, CTSLumping::Diagonal),
      dtTarget(0.0), maxAddedMassFrac(0.05), verboseSMS(false),
      lumpingSMS(CTSLumping::Diagonal), useTangentSMS(false),
      pcgTol(1.0e-10), pcgMaxIt(200), blocks(0),
      warnedLimitations(false), warnedSolver(false), warnedPCG(false)
{
}

// Main constructor
ExplicitBatheSMSConsistent::ExplicitBatheSMSConsistent(
    double p_, double dtTarget_, double maxAddedMassFrac_, bool verboseSMS_,
    int compute_critical_timestep_, bool cflUseTangent_, CTSLumping lumping_,
    double pcgTol_, int pcgMaxIt_)
    : ExplicitBathe(INTEGRATOR_TAGS_ExplicitBatheSMSConsistent, p_, compute_critical_timestep_,
                    verboseSMS_, /*cflAbort*/false, /*divergenceFactor*/0.0,
                    cflUseTangent_, /*cflRecomputeEvery*/0, lumping_),
      dtTarget(dtTarget_), maxAddedMassFrac(maxAddedMassFrac_), verboseSMS(verboseSMS_),
      lumpingSMS(lumping_), useTangentSMS(cflUseTangent_),
      pcgTol(pcgTol_), pcgMaxIt(pcgMaxIt_), blocks(0),
      warnedLimitations(false), warnedSolver(false), warnedPCG(false)
{
}

ExplicitBatheSMSConsistent::~ExplicitBatheSMSConsistent()
{
    if (blocks != 0) delete blocks;
}

int ExplicitBatheSMSConsistent::domainChanged(void)
{
    if (blocks != 0) blocks->clear();

    if (ExplicitBathe::domainChanged() < 0)
        return -1;

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "ExplicitBatheSMSConsistent::domainChanged - missing model\n";
        return -1;
    }
    if (blocks == 0) blocks = new std::vector<Ladruno::ConsistentBlock>();

    if (!warnedLimitations) {
        warnedLimitations = true;
        opserr << "ExplicitBatheSMSConsistent: consistent (Olovsson) scaling on Noh-Bathe -- "
                  "(1) elements touching an MP-constraint slave OR master node are EXCLUDED. "
                  "(2) sizing accounts for betaK damping; alphaM is not folded in. (3) both "
                  "Noh-Bathe sub-steps solve M_tilde a = r by matrix-free PCG -- REQUIRES "
                  "`system Diagonal`. (4) parallel shared/boundary nodes are not reduced.\n";
    }

    Ladruno::MassScalingReport rep =
        Ladruno::buildMassScalingConsistent(theModel, dtTarget, lumpingSMS,
                                            useTangentSMS, *blocks);

    double frac = (rep.modelMass > 0.0) ? rep.addedMass / rep.modelMass : 0.0;
    bool over = (maxAddedMassFrac > 0.0 && frac > maxAddedMassFrac);
    if (verboseSMS || over || rep.nSelfReport > 0 || rep.nMismatch > 0 || rep.nConstrained > 0) {
        opserr << "ExplicitBatheSMSConsistent: dtTarget=" << dtTarget
               << " scaled " << rep.nScaled << "/" << rep.nElems
               << " elements (Olovsson); added mass " << (100.0 * frac) << "% of model mass"
               << " (governing un-scaled dt_e=" << rep.minDtScaled << ")\n";
    }
    if (rep.nSelfReport > 0)
        opserr << "WARNING ExplicitBatheSMSConsistent: " << rep.nSelfReport
               << " throttling element(s) have a MASS-INDEPENDENT bound and were NOT scaled.\n";
    if (rep.nMismatch > 0)
        opserr << "WARNING ExplicitBatheSMSConsistent: " << rep.nMismatch
               << " element(s) had a non-node-major / DOF-mismatch mass and were SKIPPED.\n";
    if (rep.nConstrained > 0)
        opserr << "WARNING ExplicitBatheSMSConsistent: " << rep.nConstrained
               << " sub-target element(s) touch an MP-constrained node and were EXCLUDED -- "
                  "they still GOVERN at dt_e=" << rep.minDtConstrained << ".\n";
    if (over)
        opserr << "WARNING ExplicitBatheSMSConsistent: added mass " << (100.0 * frac)
               << "% exceeds -maxAddedMass cap " << (100.0 * maxAddedMassFrac) << "%\n";

    if (rep.nScaled > 0 && dynamic_cast<DiagonalSOE *>(this->getLinearSOE()) == 0)
        opserr << "WARNING ExplicitBatheSMSConsistent: " << rep.nScaled
               << " element(s) need scaling but the SOE is NOT `system Diagonal` -- the "
                  "consistent mass CANNOT be applied and dt=" << dtTarget
               << " will run UNSCALED -> expect INSTABILITY. Use `system Diagonal`.\n";

    return 0;
}

int ExplicitBatheSMSConsistent::refineAccel(Vector &a)
{
    if (blocks == 0 || blocks->empty())
        return 0;

    LinearSOE *theSOE = this->getLinearSOE();
    DiagonalSOE *theDiag = dynamic_cast<DiagonalSOE *>(theSOE);
    if (theDiag == 0) {
        if (!warnedSolver) {
            warnedSolver = true;
            opserr << "WARNING ExplicitBatheSMSConsistent::refineAccel - consistent mass "
                      "scaling REQUIRES `system Diagonal`; falling back to the un-scaled "
                      "diagonal solve (scaling NOT applied).\n";
        }
        return 0;
    }

    const double *Ainv = theDiag->getDiagonalA();   // factored: 1/mass
    int neq = theDiag->getNumEqn();
    if (Ainv == 0 || neq != a.Size())
        return 0;

    Vector r(neq);
    for (int i = 0; i < neq; ++i)
        r(i) = (Ainv[i] != 0.0) ? (a(i) / Ainv[i]) : 0.0;

    double relResid = 0.0;
    int iters = Ladruno::consistentPCG(*blocks, Ainv, neq, r, a, pcgTol, pcgMaxIt, &relResid);
    if (verboseSMS)
        opserr << "  ExplicitBatheSMSConsistent: PCG " << iters
               << " iters, rel.resid " << relResid << "\n";
    if (relResid > pcgTol && iters >= pcgMaxIt && !warnedPCG) {
        warnedPCG = true;
        opserr << "WARNING ExplicitBatheSMSConsistent::refineAccel - the mass-scaling PCG did "
                  "NOT converge (" << iters << " iters, rel.resid " << relResid
               << "); the consistent solve is incomplete -- raise -pcgMaxIt / check the mass.\n";
    }
    return 0;
}

int ExplicitBatheSMSConsistent::sendSelf(int cTag, Channel &theChannel)
{
    if (ExplicitBathe::sendSelf(cTag, theChannel) < 0)
        return -1;
    Vector data(7);
    data(0) = dtTarget;
    data(1) = maxAddedMassFrac;
    data(2) = verboseSMS ? 1.0 : 0.0;
    data(3) = (double)(int)lumpingSMS;
    data(4) = useTangentSMS ? 1.0 : 0.0;
    data(5) = pcgTol;
    data(6) = (double)pcgMaxIt;
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheSMSConsistent::sendSelf - could not send data\n";
        return -1;
    }
    return 0;
}

int ExplicitBatheSMSConsistent::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    if (ExplicitBathe::recvSelf(cTag, theChannel, theBroker) < 0)
        return -1;
    Vector data(7);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheSMSConsistent::recvSelf - could not receive data\n";
        return -1;
    }
    dtTarget         = data(0);
    maxAddedMassFrac = data(1);
    verboseSMS       = (data(2) != 0.0);
    {
        int lc = (int)data(3);
        lumpingSMS = (lc == 2) ? CTSLumping::HRZ
                   : (lc == 1) ? CTSLumping::Diagonal
                               : CTSLumping::RowSum;
    }
    useTangentSMS = (data(4) != 0.0);
    pcgTol        = data(5);
    pcgMaxIt      = (int)data(6);
    if (blocks != 0) blocks->clear();
    warnedLimitations = false;
    warnedSolver = false;
    warnedPCG = false;
    return 0;
}

void ExplicitBatheSMSConsistent::Print(OPS_Stream &s, int flag)
{
    s << "ExplicitBatheSMSConsistent (Noh-Bathe + consistent/Olovsson mass scaling, classTag 33010)\n";
    s << "  dtTarget = " << dtTarget << ", -maxAddedMass = " << maxAddedMassFrac
      << ", PCG tol = " << pcgTol << ", maxIt = " << pcgMaxIt << "\n";
    s << "  M_tilde = M_lump + sum_e beta_e[diag(m)-mm^T/M_e], solved matrix-free (PCG) at both sub-steps\n";
}
