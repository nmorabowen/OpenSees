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

void *OPS_CentralDifferenceSMSConsistent(void)
{
    // Usage: integrator CentralDifferenceSMSConsistent $dtTarget <-maxAddedMass $frac>
    //                   <-verbose> <-lump rowsum|diagonal|hrz> <-tangent>
    //                   <-pcgTol $tol> <-pcgMaxIt $n>
    // NOTE: -cflAbort and -recompute are DOWNGRADED to report-only with SMS (ADR-36 MF-1,
    // ADR-52 W1-E3a): their inherited path re-runs the element-mass eigensolve, which cannot
    // see the scaling mass; rather than reject the run we keep the integrator and report the
    // pre-scaling dt_cr instead.
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
            // W1-E3a (ADR-52): do NOT refuse the run. Under SMS these flags cannot do
            // their job -- the element-mass eigensolve can't see the scaling mass, so an
            // abort/recompute on the un-augmented pencil would be wrong (MF-1). Instead of
            // rejecting, DOWNGRADE to report-only: keep the integrator and force the
            // pre-scaling dt_cr to be reported so the user can still sanity-check the
            // stability margin.
            opserr << "NOTE CentralDifferenceSMSConsistent - " << arg << " is downgraded to "
                      "REPORT-ONLY under mass scaling (no abort/recompute on the "
                      "un-augmented element pencil, MF-1); the pre-scaling dt_cr will be "
                      "reported each domainChanged.\n";
            verboseSMS = true;
            if (strcmp(arg, "-recompute") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
                // consume the trailing N so it does not fall into the unknown-option
                // branch (peek-as-string idiom: under openseespy a numeric arg reads
                // as garbage text, but never as something starting with '-').
                const char *peek = OPS_GetString();
                if (peek != 0 && peek[0] == '-')
                    OPS_ResetCurrentInputArg(-1);   // next flag, not our N — un-read it
            }
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
    // Ladruno V4: retire our blocks from the energy-recorder conduit (owner-guarded:
    // a no-op if a newer integrator has since published).
    Ladruno::MassScalingEnergyRegistry::instance().clear(this);
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
                  "diagonal as preconditioner -- REQUIRES `system Diagonal` (serial) or "
                  "`system MPIDiagonal` (OpenSeesMP). (4) in PARALLEL the distributed PCG "
                  "reduces shared/boundary-node coupling across ranks (global inner products + "
                  "shared-DOF M_bar assembly) -- use `system MPIDiagonal`.\n";
    }

    Ladruno::MassScalingReport rep =
        Ladruno::buildMassScalingConsistent(theModel, dtTarget, lumpingSMS,
                                            useTangentSMS, *blocks);

    // Push the POST-SCALING effective stable step to the base (see the lumped
    // sibling): dtTarget capped by any still-governing excluded/self-reported step.
    {
        double effLimit = dtTarget;
        if (rep.minDtConstrained > 0.0 && rep.minDtConstrained < effLimit)
            effLimit = rep.minDtConstrained;
        if (rep.minDtSelfReport > 0.0 && rep.minDtSelfReport < effLimit)
            effLimit = rep.minDtSelfReport;
        this->setSMSEffectiveLimit(effLimit);
    }

    double frac = (rep.modelMass > 0.0) ? rep.addedMass / rep.modelMass : 0.0;
    bool over = (maxAddedMassFrac > 0.0 && frac > maxAddedMassFrac);
    if (verboseSMS || over || rep.nSelfReport > 0 || rep.nMismatch > 0 || rep.nConstrained > 0) {
        opserr << "CentralDifferenceSMSConsistent: dtTarget=" << dtTarget
               << " scaled " << rep.nScaled << "/" << rep.nElems
               << " elements (Olovsson); added mass " << (100.0 * frac)
               << "% of total element (translational) mass"
               << "; PRE-SCALING dt_cr estimate=" << rep.minDtScaled
               << " (governing un-scaled element step; AFTER scaling the run is stable "
                  "at dt <= dtTarget=" << dtTarget << ")\n";
        // Ladruno (ADR-73 P3b, panel battery-critic m-9): the line above is NOT
        // honored for overlay-owned cells — the Olovsson centroid-preserving
        // blocks under-scale the undrained coupling mode (measured; the loud
        // builder warning owns the details). Qualify rather than contradict.
        if (rep.nOverlayAugScaled > 0)
            opserr << "  EXCEPT " << rep.nOverlayAugScaled
                   << " porous-overlay element(s): consistent scaling does NOT "
                      "deliver dtTarget for them (see the under-delivery "
                      "warning; ADR-73 §12 P3b item 5)\n";
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
               << "% of total element (translational) mass exceeds -maxAddedMass cap "
               << (100.0 * maxAddedMassFrac)
               << "% (proceeding; Olovsson scaling preserves f1 far better than lumped)\n";

    // Up-front (not just per-step) non-Diagonal-SOE check: if scaling is actually needed
    // but the SOE is not Diagonal, refineAccel cannot apply the consistent mass and the run
    // will execute the un-scaled diagonal solve at dtTarget -> UNSTABLE (not merely slow).
    // Warn prominently here, before stepping, in addition to the per-step refineAccel guard.
    {
        LinearSOE *soe = this->getLinearSOE();
        bool okSOE = (dynamic_cast<DiagonalSOE *>(soe) != 0) ||
                     (soe != 0 && soe->isDistributedDiagonal());
        if (rep.nScaled > 0 && !okSOE)
            opserr << "WARNING CentralDifferenceSMSConsistent: " << rep.nScaled
                   << " element(s) need scaling but the SOE is NOT `system Diagonal`/"
                      "`system MPIDiagonal` -- the consistent mass CANNOT be applied (the PCG "
                      "preconditioner is the lumped diagonal) and dtTarget=" << dtTarget
                   << " will run UNSCALED -> expect INSTABILITY.\n";
    }

    // Ladruno V4: publish the node-major M_bar blocks (keyed by element tag) so the
    // EnergyBalanceRecorder can add the cross-node 1/2 v^T M_bar v its Node/Element
    // getMass() cannot see. Empty when nothing is scaled -> registry stays inactive.
    {
        std::map<int, Matrix> ebBlocks;
        for (size_t i = 0; i < blocks->size(); ++i)
            ebBlocks[(*blocks)[i].eleTag] = (*blocks)[i].Mbar;
        Ladruno::MassScalingEnergyRegistry::instance().publish(this, ebBlocks);
    }

    return 0;
}

int CentralDifferenceSMSConsistent::refineAccel(Vector &a)
{
    // Nothing built -> M_tilde == M_lump, the diagonal solve already in `a` is exact.
    // (In a parallel run every rank's domainChanged allocated `blocks`, so this is symmetric.)
    if (blocks == 0)
        return 0;
    // Serial OR distributed (MPIDiagonalSOE) consistent solve -> a := M_tilde^-1 r.
    // The helper picks the path at runtime via the LinearSOE base virtuals.
    return Ladruno::applyConsistentRefine(*blocks, this->getLinearSOE(), a,
                                          pcgTol, pcgMaxIt, verboseSMS,
                                          warnedSolver, warnedPCG,
                                          "CentralDifferenceSMSConsistent");
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
    Ladruno::MassScalingEnergyRegistry::instance().clear(this);   // Ladruno V4
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
