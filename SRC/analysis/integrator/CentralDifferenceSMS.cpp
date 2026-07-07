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
// See CentralDifferenceSMS.h and Ladruno_implementation/36_ladruno_selective_mass_scaling_adr.md.

#include <CentralDifferenceSMS.h>
#include <LadrunoMassScaling.h>
#include <LadrunoMassScalingEnergy.h>   // Ladruno ADR-69 P3: publish injected ΔM (KE_ms advisory)
#include <AnalysisModel.h>
#include <Domain.h>
#include <Channel.h>
#include <Vector.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

#include <string.h>

void *OPS_CentralDifferenceSMS(void)
{
    // Usage: integrator CentralDifferenceSMS $dtTarget <-maxAddedMass $frac>
    //                                        <-verbose> <-lump rowsum|diagonal|hrz>
    //                                        <-tangent> <-cfl>
    // NOTE: -cflAbort and -recompute are DOWNGRADED to report-only with SMS (ADR-36
    // MF-1, ADR-52 W1-E3a): their inherited path re-runs the element-mass eigensolve,
    // which cannot see the nodal augmentation; rather than reject the run we keep the
    // integrator and report the pre-scaling dt_cr instead.
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING integrator CentralDifferenceSMS $dtTarget <options> "
                  "- needs a target time step\n";
        return 0;
    }

    double dtTarget = 0.0;
    int nd = 1;
    if (OPS_GetDoubleInput(&nd, &dtTarget) < 0 || dtTarget <= 0.0) {
        opserr << "WARNING CentralDifferenceSMS - $dtTarget must be a positive double\n";
        return 0;
    }

    double maxAddedMassFrac = 0.05;   // default: warn past 5% added mass
    bool   verboseSMS = false;
    bool   cflUseTangent = false;
    int    compute_critical_timestep = 1;   // compute the (pre-scaling) dt_cr once
    CTSLumping lumping = CTSLumping::Diagonal;

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
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else if (strcmp(m, "hrz") == 0)      lumping = CTSLumping::HRZ;
                else opserr << "WARNING CentralDifferenceSMS - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping diagonal)\n";
            }
        } else if (strcmp(arg, "-cflAbort") == 0 || strcmp(arg, "-recompute") == 0) {
            // W1-E3a (ADR-52): do NOT refuse the run. Under SMS these flags cannot do
            // their job -- the element-mass eigensolve can't see the nodal augmentation,
            // so an abort/recompute on the un-augmented pencil would be wrong (MF-1).
            // Instead of rejecting, DOWNGRADE to report-only: keep the integrator and
            // force the pre-scaling dt_cr to be reported so the user can still
            // sanity-check the stability margin.
            opserr << "NOTE CentralDifferenceSMS - " << arg << " is downgraded to "
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
            opserr << "WARNING CentralDifferenceSMS - unknown option " << arg
                   << " (ignored)\n";
        }
    }

    TransientIntegrator *theIntegrator =
        new CentralDifferenceSMS(dtTarget, maxAddedMassFrac, verboseSMS,
                                 compute_critical_timestep, cflUseTangent, lumping);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating CentralDifferenceSMS integrator\n";
    return theIntegrator;
}

// Default constructor (broker)
CentralDifferenceSMS::CentralDifferenceSMS()
    : CentralDifferenceLadruno(INTEGRATOR_TAGS_CentralDifferenceSMS,
                               0, false, false, 0.0, false, 0, CTSLumping::Diagonal),
      dtTarget(0.0), maxAddedMassFrac(0.05), verboseSMS(false),
      lumpingSMS(CTSLumping::Diagonal), useTangentSMS(false), scaled(false),
      appliedDomain(0), warnedLimitations(false)
{
}

// Main constructor
CentralDifferenceSMS::CentralDifferenceSMS(double dtTarget_, double maxAddedMassFrac_,
                                           bool verboseSMS_,
                                           int compute_critical_timestep_,
                                           bool cflUseTangent_, CTSLumping lumping_)
    // cflAbort=false, divergenceFactor=0, cflRecomputeEvery=0 are HARD-WIRED off:
    // the parser ACCEPTS -cflAbort/-recompute but DOWNGRADES them to report-only
    // (W1-E3a), so the base never aborts/recomputes on the un-augmented element
    // pencil (MF-1).
    : CentralDifferenceLadruno(INTEGRATOR_TAGS_CentralDifferenceSMS,
                               compute_critical_timestep_, verboseSMS_, /*cflAbort*/false,
                               /*divergenceFactor*/0.0, cflUseTangent_,
                               /*cflRecomputeEvery*/0, lumping_),
      dtTarget(dtTarget_), maxAddedMassFrac(maxAddedMassFrac_), verboseSMS(verboseSMS_),
      lumpingSMS(lumping_), useTangentSMS(cflUseTangent_), scaled(false),
      appliedDomain(0), warnedLimitations(false)
{
}

CentralDifferenceSMS::~CentralDifferenceSMS()
{
    // SMS-NO-RESTORE: the injected fictitious mass is a mutation of persistent Domain
    // node state. Restore it on teardown so a LATER stage / integrator on the same
    // model (the normal gravity->modal->SMS->post workflow) sees the ORIGINAL masses.
    // We restore against the Domain we actually committed to (appliedDomain), which is
    // alive whenever the integrator is swapped out mid-session; at full program exit
    // node-mass correctness is moot.
    if (scaled && appliedDomain != 0 && !injected.empty())
        Ladruno::applyMassScaling(appliedDomain, injected, -1.0);
    // Ladruno ADR-69 P3: retire the published ΔM (owner-guarded).
    Ladruno::MassScalingEnergyRegistry::instance().clearNodal(this);
}

void CentralDifferenceSMS::removeScaling(void)
{
    if (scaled && appliedDomain != 0 && !injected.empty())
        Ladruno::applyMassScaling(appliedDomain, injected, -1.0);  // subtract baseline
    injected.clear();
    scaled = false;
    appliedDomain = 0;
    // Ladruno ADR-69 P3: the registry must never outlive the Domain injection.
    Ladruno::MassScalingEnergyRegistry::instance().clearNodal(this);
}

int CentralDifferenceSMS::domainChanged(void)
{
    // 1) re-baseline: undo any prior injection so we add onto the ORIGINAL masses
    //    (re-entrant domainChanged must not compound).
    removeScaling();

    // 2) base: (re)allocate leap-frog state + compute the pre-scaling dt_cr.
    if (CentralDifferenceLadruno::domainChanged() < 0)
        return -1;

    // 3) scaling pass: size per-element fictitious mass and inject it into nodes.
    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = (theModel != 0) ? theModel->getDomainPtr() : 0;
    if (theModel == 0 || theDomain == 0) {
        opserr << "CentralDifferenceSMS::domainChanged - missing model/domain\n";
        return -1;
    }

    // one-time v1-limitation warnings (review SMS-CONSTRAINTS / -MASS-DAMP / -BETAK / -MPI).
    // These are honest disclosures of what v1 does NOT handle, mirroring the -cflAbort
    // hard-reject: the feature still runs, but the user must heed them.
    if (!warnedLimitations) {
        warnedLimitations = true;
        opserr << "CentralDifferenceSMS: v1 limitations -- (1) elements touching an "
                  "MP-constraint SLAVE node (equalDOF/rigidDiaphragm/rigidLink) are EXCLUDED "
                  "from scaling (their injected mass would not land through the constraint); "
                  "they remain governing, so dtTarget is NOT delivered for them -- lower dt or "
                  "remove the constraint (count reported below). SP/fix nodes are fine; RBE2/"
                  "RBE3 couplings are elements handled by the self-report skip, not here. "
                  "(2) sizing now accounts for stiffness-proportional Rayleigh damping -- ALL "
                  "three slots betaK+betaKinit+betaKcomm summed -- (it shrinks the explicit "
                  "step) via the closed-form s=T^2+2*T*betaK/dt_e; "
                  "mass-proportional (alphaM) damping is NOT folded in (it does not reduce the "
                  "high-frequency stable step). (3) in PARALLEL the injected lumped mass on a "
                  "shared/boundary node IS summed across ranks by a distributed/MPI diagonal "
                  "solver (OpenSeesMP `system MPIDiagonal`, OpenSeesSP `system Diagonal`) -- the "
                  "explicit M^-1 solve reads that reduced diagonal (validated bit-identical to "
                  "serial). The CONSISTENT (Olovsson) variant is ALSO parallel-safe (ADR-38 V5): "
                  "its matrix-free PCG uses GLOBAL inner products (globalReduceSum) and shared-DOF "
                  "assembly under `system MPIDiagonal` (see LadrunoConsistentRefine.h).\n";
    }

    Ladruno::MassScalingReport rep =
        Ladruno::buildMassScaling(theModel, dtTarget, lumpingSMS, useTangentSMS, injected);
    Ladruno::applyMassScaling(theDomain, injected, +1.0);
    scaled = true;
    appliedDomain = theDomain;   // remember where we committed, for restore (M5)

    // Ladruno ADR-69 P3: publish the committed ΔM so the -v2 EnergyBalanceRecorder
    // can report the fictitious added-mass KE share (KE_ms advisory column). The
    // injection is already inside Node::getMass() — the recorder MEASURES, never
    // adds. Empty map when nothing scaled -> registry side stays inactive.
    Ladruno::MassScalingEnergyRegistry::instance().publishNodal(this, injected);

    // Push the POST-SCALING effective stable step to the base so its step-1 dt_cr
    // report stops warning "expect INSTABILITY" against the pre-scaling pencil:
    // dtTarget, capped by any excluded (constrained) or self-reported element that
    // scaling could not lift and therefore still governs.
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
        opserr << "CentralDifferenceSMS: dtTarget=" << dtTarget
               << " scaled " << rep.nScaled << "/" << rep.nElems
               << " elements; added mass " << (100.0 * frac)
               << "% of total element (translational) mass"
               << "; PRE-SCALING dt_cr estimate=" << rep.minDtScaled
               << " (governing un-scaled element step; AFTER scaling the run is stable "
                  "at dt <= dtTarget=" << dtTarget << ")\n";
    }
    if (rep.nSelfReport > 0)
        opserr << "WARNING CentralDifferenceSMS: " << rep.nSelfReport
               << " throttling element(s) have a MASS-INDEPENDENT (self-reported) "
                  "stability bound and were NOT scaled -- lower dtTarget to satisfy them.\n";
    if (rep.nMismatch > 0)
        opserr << "WARNING CentralDifferenceSMS: " << rep.nMismatch
               << " element(s) had a non-node-major / DOF-mismatch mass and were SKIPPED "
                  "(not scaled) -- their dt is not recovered.\n";
    if (rep.nConstrained > 0)
        opserr << "WARNING CentralDifferenceSMS: " << rep.nConstrained
               << " sub-target element(s) touch an MP-constrained (slave) node and were "
                  "EXCLUDED from scaling -- they still GOVERN at dt_e=" << rep.minDtConstrained
               << " < dtTarget=" << dtTarget << "; lower dt to dt_e or remove the constraint.\n";
    if (over)
        opserr << "WARNING CentralDifferenceSMS: added mass " << (100.0 * frac)
               << "% of total element (translational) mass exceeds -maxAddedMass cap "
               << (100.0 * maxAddedMassFrac)
               << "% (proceeding; the scaled inertia shifts global frequencies)\n";

    return 0;
}

int CentralDifferenceSMS::sendSelf(int cTag, Channel &theChannel)
{
    // chain: base sends its Vector(7), then SMS its own params. Correct for the
    // streaming (MPI) channel; the FileDatastore dbTag-collision case is a v1
    // limitation (mass scaling is sequential-first).
    if (CentralDifferenceLadruno::sendSelf(cTag, theChannel) < 0)
        return -1;
    Vector data(5);
    data(0) = dtTarget;
    data(1) = maxAddedMassFrac;
    data(2) = verboseSMS ? 1.0 : 0.0;
    data(3) = (double)(int)lumpingSMS;
    data(4) = useTangentSMS ? 1.0 : 0.0;
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "CentralDifferenceSMS::sendSelf - could not send data\n";
        return -1;
    }
    return 0;
}

int CentralDifferenceSMS::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    if (CentralDifferenceLadruno::recvSelf(cTag, theChannel, theBroker) < 0)
        return -1;
    Vector data(5);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "CentralDifferenceSMS::recvSelf - could not receive data\n";
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
    // If THIS object is live and already scaled (recvSelf on a reused object, not
    // the usual fresh broker construction), restore the injected nodal mass FIRST —
    // clearing `injected` below would otherwise leak the fictitious ΔM into the
    // Domain permanently (review 2026-07-01).
    if (scaled && appliedDomain != 0 && !injected.empty())
        Ladruno::applyMassScaling(appliedDomain, injected, -1.0);
    injected.clear();
    scaled = false;
    appliedDomain = 0;
    Ladruno::MassScalingEnergyRegistry::instance().clearNodal(this);   // Ladruno ADR-69 P3
    warnedLimitations = false;
    return 0;
}

void CentralDifferenceSMS::Print(OPS_Stream &s, int flag)
{
    s << "CentralDifferenceSMS (selective mass scaling, classTag 33007)\n";
    s << "  dtTarget = " << dtTarget << ", -maxAddedMass = " << maxAddedMassFrac << "\n";
    s << "  recipe: lumped/diagonal mass, system Diagonal, algorithm Linear, dt <= dtTarget\n";
}
