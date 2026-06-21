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
// See ExplicitBatheLNVDSMSConsistent.h and
// Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md.

#include <ExplicitBatheLNVDSMSConsistent.h>
#include <LadrunoMassScaling.h>
#include <LadrunoConsistentRefine.h>    // Ladruno V5: shared serial+MPI refineAccel body
#include <LadrunoMassScalingEnergy.h>   // Ladruno V4: publish M_bar blocks to the energy recorder
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

void *OPS_ExplicitBatheLNVDSMSConsistent(void)
{
    // Usage: integrator ExplicitBatheLNVDSMSConsistent $p $alpha $dtTarget <options>
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING integrator ExplicitBatheLNVDSMSConsistent $p $alpha $dtTarget "
                  "<options> - needs the Noh-Bathe p, FLAC alpha, and a target time step\n";
        return 0;
    }

    double in[3] = {0.0, 0.0, 0.0};
    int n3 = 3;
    if (OPS_GetDoubleInput(&n3, in) < 0) {
        opserr << "WARNING ExplicitBatheLNVDSMSConsistent - could not read $p $alpha $dtTarget\n";
        return 0;
    }
    double p = in[0], alpha = in[1], dtTarget = in[2];
    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING ExplicitBatheLNVDSMSConsistent - $p must be in (0,1)\n";
        return 0;
    }
    if (alpha < 0.0 || alpha >= 1.0) {   // match the ExplicitBatheLNVD base contract
        opserr << "WARNING ExplicitBatheLNVDSMSConsistent - $alpha (FLAC local damping) "
                  "must be in [0,1) (got " << alpha << ")\n";
        return 0;
    }
    if (dtTarget <= 0.0) {
        opserr << "WARNING ExplicitBatheLNVDSMSConsistent - $dtTarget must be positive\n";
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
                else opserr << "WARNING ExplicitBatheLNVDSMSConsistent - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping diagonal)\n";
            }
        } else if (strcmp(arg, "-cflAbort") == 0 || strcmp(arg, "-recompute") == 0) {
            opserr << "WARNING ExplicitBatheLNVDSMSConsistent - " << arg << " is not supported "
                      "with mass scaling. Integrator NOT created.\n";
            return 0;
        } else {
            opserr << "WARNING ExplicitBatheLNVDSMSConsistent - unknown option " << arg << " (ignored)\n";
        }
    }

    TransientIntegrator *theIntegrator =
        new ExplicitBatheLNVDSMSConsistent(p, alpha, dtTarget, maxAddedMassFrac, verboseSMS,
                                           compute_critical_timestep, cflUseTangent, lumping,
                                           pcgTol, pcgMaxIt);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitBatheLNVDSMSConsistent integrator\n";
    return theIntegrator;
}

// Default constructor (broker)
ExplicitBatheLNVDSMSConsistent::ExplicitBatheLNVDSMSConsistent()
    : ExplicitBatheLNVD(INTEGRATOR_TAGS_ExplicitBatheLNVDSMSConsistent, 0.0, 0.0, 0, false,
                        false, 0.0, false, 0, CTSLumping::Diagonal),
      dtTarget(0.0), maxAddedMassFrac(0.05), verboseSMS(false),
      lumpingSMS(CTSLumping::Diagonal), useTangentSMS(false),
      pcgTol(1.0e-10), pcgMaxIt(200), blocks(0),
      warnedLimitations(false), warnedSolver(false), warnedPCG(false)
{
}

// Main constructor
ExplicitBatheLNVDSMSConsistent::ExplicitBatheLNVDSMSConsistent(
    double p_, double alpha_, double dtTarget_, double maxAddedMassFrac_, bool verboseSMS_,
    int compute_critical_timestep_, bool cflUseTangent_, CTSLumping lumping_,
    double pcgTol_, int pcgMaxIt_)
    : ExplicitBatheLNVD(INTEGRATOR_TAGS_ExplicitBatheLNVDSMSConsistent, p_, alpha_,
                        compute_critical_timestep_, verboseSMS_, /*cflAbort*/false,
                        /*divergenceFactor*/0.0, cflUseTangent_, /*cflRecomputeEvery*/0,
                        lumping_),
      dtTarget(dtTarget_), maxAddedMassFrac(maxAddedMassFrac_), verboseSMS(verboseSMS_),
      lumpingSMS(lumping_), useTangentSMS(cflUseTangent_),
      pcgTol(pcgTol_), pcgMaxIt(pcgMaxIt_), blocks(0),
      warnedLimitations(false), warnedSolver(false), warnedPCG(false)
{
}

ExplicitBatheLNVDSMSConsistent::~ExplicitBatheLNVDSMSConsistent()
{
    Ladruno::MassScalingEnergyRegistry::instance().clear(this);   // Ladruno V4 (owner-guarded)
    if (blocks != 0) delete blocks;
}

int ExplicitBatheLNVDSMSConsistent::domainChanged(void)
{
    if (blocks != 0) blocks->clear();

    if (ExplicitBatheLNVD::domainChanged() < 0)
        return -1;

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "ExplicitBatheLNVDSMSConsistent::domainChanged - missing model\n";
        return -1;
    }
    if (blocks == 0) blocks = new std::vector<Ladruno::ConsistentBlock>();

    if (!warnedLimitations) {
        warnedLimitations = true;
        opserr << "ExplicitBatheLNVDSMSConsistent: consistent (Olovsson) scaling on Noh-Bathe+"
                  "LNVD -- (1) elements touching an MP-constraint slave OR master node are "
                  "EXCLUDED. (2) sizing accounts for betaK damping; alphaM and the FLAC alpha "
                  "are separate. (3) both sub-steps solve M_tilde a = r by matrix-free PCG -- "
                  "REQUIRES `system Diagonal` (serial) or `system MPIDiagonal` (OpenSeesMP). "
                  "(4) in PARALLEL the distributed PCG reduces shared/boundary-node coupling "
                  "across ranks -- use `system MPIDiagonal`.\n";
    }

    Ladruno::MassScalingReport rep =
        Ladruno::buildMassScalingConsistent(theModel, dtTarget, lumpingSMS,
                                            useTangentSMS, *blocks);

    double frac = (rep.modelMass > 0.0) ? rep.addedMass / rep.modelMass : 0.0;
    bool over = (maxAddedMassFrac > 0.0 && frac > maxAddedMassFrac);
    if (verboseSMS || over || rep.nSelfReport > 0 || rep.nMismatch > 0 || rep.nConstrained > 0) {
        opserr << "ExplicitBatheLNVDSMSConsistent: dtTarget=" << dtTarget
               << " scaled " << rep.nScaled << "/" << rep.nElems
               << " elements (Olovsson); added mass " << (100.0 * frac) << "% of model mass"
               << " (governing un-scaled dt_e=" << rep.minDtScaled << ")\n";
    }
    if (rep.nSelfReport > 0)
        opserr << "WARNING ExplicitBatheLNVDSMSConsistent: " << rep.nSelfReport
               << " throttling element(s) have a MASS-INDEPENDENT bound and were NOT scaled.\n";
    if (rep.nMismatch > 0)
        opserr << "WARNING ExplicitBatheLNVDSMSConsistent: " << rep.nMismatch
               << " element(s) had a non-node-major / DOF-mismatch mass and were SKIPPED.\n";
    if (rep.nConstrained > 0)
        opserr << "WARNING ExplicitBatheLNVDSMSConsistent: " << rep.nConstrained
               << " sub-target element(s) touch an MP-constrained node and were EXCLUDED -- "
                  "they still GOVERN at dt_e=" << rep.minDtConstrained << ".\n";
    if (over)
        opserr << "WARNING ExplicitBatheLNVDSMSConsistent: added mass " << (100.0 * frac)
               << "% exceeds -maxAddedMass cap " << (100.0 * maxAddedMassFrac) << "%\n";

    {
        LinearSOE *soe = this->getLinearSOE();
        bool okSOE = (dynamic_cast<DiagonalSOE *>(soe) != 0) ||
                     (soe != 0 && soe->isDistributedDiagonal());
        if (rep.nScaled > 0 && !okSOE)
            opserr << "WARNING ExplicitBatheLNVDSMSConsistent: " << rep.nScaled
                   << " element(s) need scaling but the SOE is NOT `system Diagonal`/"
                      "`system MPIDiagonal` -- the consistent mass CANNOT be applied and dt="
                   << dtTarget << " will run UNSCALED -> expect INSTABILITY.\n";
    }

    // Ladruno V4: publish the node-major M_bar blocks (keyed by element tag) for the
    // EnergyBalanceRecorder. Empty when nothing is scaled -> registry stays inactive.
    {
        std::map<int, Matrix> ebBlocks;
        for (size_t i = 0; i < blocks->size(); ++i)
            ebBlocks[(*blocks)[i].eleTag] = (*blocks)[i].Mbar;
        Ladruno::MassScalingEnergyRegistry::instance().publish(this, ebBlocks);
    }

    return 0;
}

int ExplicitBatheLNVDSMSConsistent::refineAccel(Vector &a)
{
    if (blocks == 0)
        return 0;
    return Ladruno::applyConsistentRefine(*blocks, this->getLinearSOE(), a,
                                          pcgTol, pcgMaxIt, verboseSMS,
                                          warnedSolver, warnedPCG,
                                          "ExplicitBatheLNVDSMSConsistent");
}

int ExplicitBatheLNVDSMSConsistent::sendSelf(int cTag, Channel &theChannel)
{
    if (ExplicitBatheLNVD::sendSelf(cTag, theChannel) < 0)
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
        opserr << "ExplicitBatheLNVDSMSConsistent::sendSelf - could not send data\n";
        return -1;
    }
    return 0;
}

int ExplicitBatheLNVDSMSConsistent::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    if (ExplicitBatheLNVD::recvSelf(cTag, theChannel, theBroker) < 0)
        return -1;
    Vector data(7);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheLNVDSMSConsistent::recvSelf - could not receive data\n";
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
    Ladruno::MassScalingEnergyRegistry::instance().clear(this);   // Ladruno V4
    warnedLimitations = false;
    warnedSolver = false;
    warnedPCG = false;
    return 0;
}

void ExplicitBatheLNVDSMSConsistent::Print(OPS_Stream &s, int flag)
{
    s << "ExplicitBatheLNVDSMSConsistent (Noh-Bathe + FLAC LNVD + consistent/Olovsson mass "
         "scaling, classTag 33012)\n";
    s << "  dtTarget = " << dtTarget << ", -maxAddedMass = " << maxAddedMassFrac
      << ", PCG tol = " << pcgTol << ", maxIt = " << pcgMaxIt << "\n";
}
