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
// See CentralDifferenceSMSConsistent.h and
// Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md.

#include <CentralDifferenceSMSConsistent.h>
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

void *OPS_CentralDifferenceSMSConsistent(void)
{
    // Usage: integrator CentralDifferenceSMSConsistent $dtTarget <-maxAddedMass $frac>
    //                   <-verbose> <-lump rowsum|diagonal|hrz> <-tangent>
    //                   <-pcgTol $tol> <-pcgMaxIt $n>
    // NOTE: -cflAbort / -recompute are REJECTED (same MF-1 rationale as the lumped SMS).
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING integrator CentralDifferenceSMSConsistent $dtTarget <options> "
                  "- needs a target time step\n";
        return 0;
    }

    double dtTarget = 0.0;
    int nd = 1;
    if (OPS_GetDoubleInput(&nd, &dtTarget) < 0 || dtTarget <= 0.0) {
        opserr << "WARNING CentralDifferenceSMSConsistent - $dtTarget must be a positive double\n";
        return 0;
    }

    double maxAddedMassFrac = 0.05;
    bool   verboseSMS = false;
    bool   cflUseTangent = false;
    int    compute_critical_timestep = 1;
    CTSLumping lumping = CTSLumping::Diagonal;   // matches the `system Diagonal` recipe
    double pcgTol = 1.0e-10;
    int    pcgMaxIt = 200;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        if (strcmp(arg, "-verbose") == 0) {
            verboseSMS = true;
        } else if (strcmp(arg, "-maxAddedMass") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int n1 = 1;
                OPS_GetDoubleInput(&n1, &maxAddedMassFrac);
            }
        } else if (strcmp(arg, "-tangent") == 0) {
            cflUseTangent = true;
        } else if (strcmp(arg, "-cfl") == 0 || strcmp(arg, "-criticalTimestep") == 0) {
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-pcgTol") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int n1 = 1;
                OPS_GetDoubleInput(&n1, &pcgTol);
            }
        } else if (strcmp(arg, "-pcgMaxIt") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int n1 = 1;
                OPS_GetIntInput(&n1, &pcgMaxIt);
            }
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else if (strcmp(m, "hrz") == 0)      lumping = CTSLumping::HRZ;
                else opserr << "WARNING CentralDifferenceSMSConsistent - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping diagonal)\n";
            }
        } else if (strcmp(arg, "-cflAbort") == 0 || strcmp(arg, "-recompute") == 0) {
            opserr << "WARNING CentralDifferenceSMSConsistent - " << arg << " is not "
                      "supported with mass scaling (it re-runs the element-mass eigensolve, "
                      "which cannot see the scaling mass, and would spuriously abort). "
                      "Integrator NOT created.\n";
            return 0;
        } else {
            opserr << "WARNING CentralDifferenceSMSConsistent - unknown option " << arg
                   << " (ignored)\n";
        }
    }

    TransientIntegrator *theIntegrator =
        new CentralDifferenceSMSConsistent(dtTarget, maxAddedMassFrac, verboseSMS,
                                           compute_critical_timestep, cflUseTangent,
                                           lumping, pcgTol, pcgMaxIt);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating CentralDifferenceSMSConsistent integrator\n";
    return theIntegrator;
}

// Default constructor (broker)
CentralDifferenceSMSConsistent::CentralDifferenceSMSConsistent()
    : CentralDifferenceLadruno(INTEGRATOR_TAGS_CentralDifferenceSMSConsistent,
                               0, false, false, 0.0, false, 0, CTSLumping::Diagonal),
      dtTarget(0.0), maxAddedMassFrac(0.05), verboseSMS(false),
      lumpingSMS(CTSLumping::Diagonal), useTangentSMS(false),
      pcgTol(1.0e-10), pcgMaxIt(200), blocks(0),
      warnedLimitations(false), warnedSolver(false), warnedPCG(false)
{
}

// Main constructor
CentralDifferenceSMSConsistent::CentralDifferenceSMSConsistent(
    double dtTarget_, double maxAddedMassFrac_, bool verboseSMS_,
    int compute_critical_timestep_, bool cflUseTangent_, CTSLumping lumping_,
    double pcgTol_, int pcgMaxIt_)
    // cflAbort / divergenceFactor / cflRecomputeEvery hard-wired off (MF-1).
    : CentralDifferenceLadruno(INTEGRATOR_TAGS_CentralDifferenceSMSConsistent,
                               compute_critical_timestep_, verboseSMS_, /*cflAbort*/false,
                               /*divergenceFactor*/0.0, cflUseTangent_,
                               /*cflRecomputeEvery*/0, lumping_),
      dtTarget(dtTarget_), maxAddedMassFrac(maxAddedMassFrac_), verboseSMS(verboseSMS_),
      lumpingSMS(lumping_), useTangentSMS(cflUseTangent_),
      pcgTol(pcgTol_), pcgMaxIt(pcgMaxIt_), blocks(0),
      warnedLimitations(false), warnedSolver(false), warnedPCG(false)
{
}

CentralDifferenceSMSConsistent::~CentralDifferenceSMSConsistent()
{
    // No Domain mutation to restore (the consistent scaling mass lives only in `blocks`,
    // integrator-owned scratch consumed by the matrix-free PCG).
    if (blocks != 0) delete blocks;
}

int CentralDifferenceSMSConsistent::domainChanged(void)
{
    // (re)build the scaling blocks from scratch onto the ORIGINAL (un-mutated) masses.
    if (blocks != 0) blocks->clear();

    // base: (re)allocate leap-frog state + compute the pre-scaling dt_cr.
    if (CentralDifferenceLadruno::domainChanged() < 0)
        return -1;

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "CentralDifferenceSMSConsistent::domainChanged - missing model\n";
        return -1;
    }
    if (blocks == 0) blocks = new std::vector<Ladruno::ConsistentBlock>();

    if (!warnedLimitations) {
        warnedLimitations = true;
        opserr << "CentralDifferenceSMSConsistent: consistent (Olovsson) scaling -- (1) "
                  "elements touching an MP-constraint SLAVE node are EXCLUDED (same hazard "
                  "as the lumped SMS); they remain governing. (2) sizing accounts for betaK "
                  "Rayleigh damping (closed-form s=T^2+2*T*betaK/dt_e); alphaM is not folded "
                  "in. (3) the step solves M_tilde a = r by matrix-free PCG with the lumped "
                  "diagonal as preconditioner -- REQUIRES `system Diagonal`. (4) parallel "
                  "shared/boundary nodes are not reduced across ranks (sequential use).\n";
    }

    Ladruno::MassScalingReport rep =
        Ladruno::buildMassScalingConsistent(theModel, dtTarget, lumpingSMS,
                                            useTangentSMS, *blocks);

    double frac = (rep.modelMass > 0.0) ? rep.addedMass / rep.modelMass : 0.0;
    bool over = (maxAddedMassFrac > 0.0 && frac > maxAddedMassFrac);
    if (verboseSMS || over || rep.nSelfReport > 0 || rep.nMismatch > 0 || rep.nConstrained > 0) {
        opserr << "CentralDifferenceSMSConsistent: dtTarget=" << dtTarget
               << " scaled " << rep.nScaled << "/" << rep.nElems
               << " elements (Olovsson); added mass " << (100.0 * frac) << "% of model mass"
               << " (governing un-scaled dt_e=" << rep.minDtScaled << ")\n";
    }
    if (rep.nSelfReport > 0)
        opserr << "WARNING CentralDifferenceSMSConsistent: " << rep.nSelfReport
               << " throttling element(s) have a MASS-INDEPENDENT (self-reported) bound "
                  "and were NOT scaled -- lower dtTarget to satisfy them.\n";
    if (rep.nMismatch > 0)
        opserr << "WARNING CentralDifferenceSMSConsistent: " << rep.nMismatch
               << " element(s) had a non-node-major / DOF-mismatch mass and were SKIPPED.\n";
    if (rep.nConstrained > 0)
        opserr << "WARNING CentralDifferenceSMSConsistent: " << rep.nConstrained
               << " sub-target element(s) touch an MP-constrained (slave) node and were "
                  "EXCLUDED -- they still GOVERN at dt_e=" << rep.minDtConstrained
               << " < dtTarget=" << dtTarget << ".\n";
    if (over)
        opserr << "WARNING CentralDifferenceSMSConsistent: added mass " << (100.0 * frac)
               << "% exceeds -maxAddedMass cap " << (100.0 * maxAddedMassFrac)
               << "% (proceeding; Olovsson scaling preserves f1 far better than lumped)\n";

    // Up-front (not just per-step) non-Diagonal-SOE check: if scaling is actually needed
    // but the SOE is not Diagonal, refineAccel cannot apply the consistent mass and the run
    // will execute the un-scaled diagonal solve at dtTarget -> UNSTABLE (not merely slow).
    // Warn prominently here, before stepping, in addition to the per-step refineAccel guard.
    if (rep.nScaled > 0 && dynamic_cast<DiagonalSOE *>(this->getLinearSOE()) == 0)
        opserr << "WARNING CentralDifferenceSMSConsistent: " << rep.nScaled
               << " element(s) need scaling but the SOE is NOT `system Diagonal` -- the "
                  "consistent mass CANNOT be applied (the PCG preconditioner is the lumped "
                  "diagonal) and dtTarget=" << dtTarget << " will run UNSCALED -> expect "
                  "INSTABILITY. Use `system Diagonal`.\n";

    return 0;
}

int CentralDifferenceSMSConsistent::refineAccel(Vector &a)
{
    // Nothing scaled -> M_tilde == M_lump, the diagonal solve already in `a` is exact.
    if (blocks == 0 || blocks->empty())
        return 0;

    LinearSOE *theSOE = this->getLinearSOE();
    DiagonalSOE *theDiag = dynamic_cast<DiagonalSOE *>(theSOE);
    if (theDiag == 0) {
        if (!warnedSolver) {
            warnedSolver = true;
            opserr << "WARNING CentralDifferenceSMSConsistent::refineAccel - consistent "
                      "(Olovsson) mass scaling REQUIRES `system Diagonal` (the lumped "
                      "diagonal is the PCG preconditioner). The current SOE is not a "
                      "DiagonalSOE; falling back to the un-scaled diagonal solve (the "
                      "scaling mass is NOT applied -- expect the un-raised stable step).\n";
        }
        return 0;
    }

    // Post-solve the DiagonalDirectSolver stored A[i] = 1/mass_i (the factored diagonal).
    // That IS the Jacobi preconditioner; the lumped mass is 1/A[i]. The incoming `a` is
    // the diagonal solve M_lump^-1 r, so recover r = M_lump .* a = a / A.
    const double *Ainv = theDiag->getDiagonalA();
    int neq = theDiag->getNumEqn();
    if (Ainv == 0 || neq != a.Size())
        return 0;

    Vector r(neq);
    for (int i = 0; i < neq; ++i)
        r(i) = (Ainv[i] != 0.0) ? (a(i) / Ainv[i]) : 0.0;

    double relResid = 0.0;
    int iters = Ladruno::consistentPCG(*blocks, Ainv, neq, r, a,
                                       pcgTol, pcgMaxIt, &relResid);
    if (verboseSMS)
        opserr << "  CentralDifferenceSMSConsistent: PCG " << iters
               << " iters, rel.resid " << relResid << "\n";
    // M_tilde is SPD so CG must converge; a non-converged step means a near-singular
    // M_tilde (e.g. a zero-mass free DOF) -> the accel `a` is only a partial solve.
    // Surface it once so a silently-wrong run is not mistaken for success.
    if (relResid > pcgTol && iters >= pcgMaxIt && !warnedPCG) {
        warnedPCG = true;
        opserr << "WARNING CentralDifferenceSMSConsistent::refineAccel - the mass-scaling "
                  "PCG did NOT converge (" << iters << " iters, rel.resid " << relResid
               << " > tol " << pcgTol << "); the consistent mass solve is incomplete this "
                  "step. Likely a near-singular M_tilde (a zero-mass free DOF) -- raise "
                  "-pcgMaxIt or check the mass distribution.\n";
    }
    return 0;
}

int CentralDifferenceSMSConsistent::sendSelf(int cTag, Channel &theChannel)
{
    if (CentralDifferenceLadruno::sendSelf(cTag, theChannel) < 0)
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
        opserr << "CentralDifferenceSMSConsistent::sendSelf - could not send data\n";
        return -1;
    }
    return 0;
}

int CentralDifferenceSMSConsistent::recvSelf(int cTag, Channel &theChannel,
                                             FEM_ObjectBroker &theBroker)
{
    if (CentralDifferenceLadruno::recvSelf(cTag, theChannel, theBroker) < 0)
        return -1;
    Vector data(7);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "CentralDifferenceSMSConsistent::recvSelf - could not receive data\n";
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
    if (blocks != 0) { blocks->clear(); }
    warnedLimitations = false;
    warnedSolver = false;
    warnedPCG = false;
    return 0;
}

void CentralDifferenceSMSConsistent::Print(OPS_Stream &s, int flag)
{
    s << "CentralDifferenceSMSConsistent (consistent/Olovsson mass scaling, classTag 33008)\n";
    s << "  dtTarget = " << dtTarget << ", -maxAddedMass = " << maxAddedMassFrac
      << ", PCG tol = " << pcgTol << ", maxIt = " << pcgMaxIt << "\n";
    s << "  recipe: lumped diagonal precond, system Diagonal, algorithm Linear, dt <= dtTarget;\n";
    s << "          M_tilde = M_lump + sum_e beta_e[diag(m)-mm^T/M_e], solved matrix-free (PCG)\n";
}
