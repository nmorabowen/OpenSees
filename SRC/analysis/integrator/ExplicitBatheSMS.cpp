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
// See ExplicitBatheSMS.h and Ladruno_implementation/36_ladruno_selective_mass_scaling_adr.md.

#include <ExplicitBatheSMS.h>
#include <LadrunoMassScaling.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <Channel.h>
#include <Vector.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

#include <string.h>

void *OPS_ExplicitBatheSMS(void)
{
    // Usage: integrator ExplicitBatheSMS $p $dtTarget <-maxAddedMass f>
    //                   <-lump rowsum|diagonal|hrz> <-tangent> <-verbose>
    // NOTE: -cflAbort and -recompute are DOWNGRADED to report-only with SMS (ADR-36
    // MF-1, ADR-52 W1-E3a): their inherited path re-runs the element-mass eigensolve,
    // which cannot see the nodal augmentation; rather than reject the run we keep the
    // integrator and report the pre-scaling dt_cr instead.
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING integrator ExplicitBatheSMS $p $dtTarget <options> "
                  "- needs the Noh-Bathe p and a target time step\n";
        return 0;
    }

    double pin[2] = {0.0, 0.0};
    int n2 = 2;
    if (OPS_GetDoubleInput(&n2, pin) < 0) {
        opserr << "WARNING ExplicitBatheSMS - could not read $p $dtTarget\n";
        return 0;
    }
    double p = pin[0], dtTarget = pin[1];
    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING ExplicitBatheSMS - $p must be in (0,1) (Noh-Bathe sub-step)\n";
        return 0;
    }
    if (dtTarget <= 0.0) {
        opserr << "WARNING ExplicitBatheSMS - $dtTarget must be a positive double\n";
        return 0;
    }

    double maxAddedMassFrac = 0.05;
    bool   verboseSMS = false;
    bool   cflUseTangent = false;
    int    compute_critical_timestep = 1;
    CTSLumping lumping = CTSLumping::Diagonal;   // matches the system Diagonal run

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
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else if (strcmp(m, "hrz") == 0)      lumping = CTSLumping::HRZ;
                else opserr << "WARNING ExplicitBatheSMS - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping diagonal)\n";
            }
        } else if (strcmp(arg, "-cflAbort") == 0 || strcmp(arg, "-recompute") == 0) {
            // W1-E3a (ADR-52): do NOT refuse the run. Under SMS these flags cannot do
            // their job -- the element-mass eigensolve can't see the nodal augmentation,
            // so an abort/recompute on the un-augmented pencil would be wrong (MF-1).
            // Instead of rejecting, DOWNGRADE to report-only: keep the integrator and
            // force the pre-scaling dt_cr to be reported so the user can still
            // sanity-check the stability margin.
            opserr << "NOTE ExplicitBatheSMS - " << arg << " is downgraded to "
                      "REPORT-ONLY under mass scaling (no abort/recompute on the "
                      "un-augmented element pencil, MF-1); the pre-scaling dt_cr will be "
                      "reported each domainChanged.\n";
            verboseSMS = true;
        } else {
            opserr << "WARNING ExplicitBatheSMS - unknown option " << arg << " (ignored)\n";
        }
    }

    TransientIntegrator *theIntegrator =
        new ExplicitBatheSMS(p, dtTarget, maxAddedMassFrac, verboseSMS,
                             compute_critical_timestep, cflUseTangent, lumping);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitBatheSMS integrator\n";
    return theIntegrator;
}

// Default constructor (broker)
ExplicitBatheSMS::ExplicitBatheSMS()
    : ExplicitBathe(INTEGRATOR_TAGS_ExplicitBatheSMS, 0.0, 0, false, false, 0.0,
                    false, 0, CTSLumping::Diagonal),
      dtTarget(0.0), maxAddedMassFrac(0.05), verboseSMS(false),
      lumpingSMS(CTSLumping::Diagonal), useTangentSMS(false), scaled(false),
      appliedDomain(0), warnedLimitations(false)
{
}

// Main constructor
ExplicitBatheSMS::ExplicitBatheSMS(double p_, double dtTarget_, double maxAddedMassFrac_,
                                   bool verboseSMS_, int compute_critical_timestep_,
                                   bool cflUseTangent_, CTSLumping lumping_)
    : ExplicitBathe(INTEGRATOR_TAGS_ExplicitBatheSMS, p_, compute_critical_timestep_,
                    verboseSMS_, /*cflAbort*/false, /*divergenceFactor*/0.0,
                    cflUseTangent_, /*cflRecomputeEvery*/0, lumping_),
      dtTarget(dtTarget_), maxAddedMassFrac(maxAddedMassFrac_), verboseSMS(verboseSMS_),
      lumpingSMS(lumping_), useTangentSMS(cflUseTangent_), scaled(false),
      appliedDomain(0), warnedLimitations(false)
{
}

ExplicitBatheSMS::~ExplicitBatheSMS()
{
    if (scaled && appliedDomain != 0 && !injected.empty())
        Ladruno::applyMassScaling(appliedDomain, injected, -1.0);
}

void ExplicitBatheSMS::removeScaling(void)
{
    if (scaled && appliedDomain != 0 && !injected.empty())
        Ladruno::applyMassScaling(appliedDomain, injected, -1.0);
    injected.clear();
    scaled = false;
    appliedDomain = 0;
}

int ExplicitBatheSMS::domainChanged(void)
{
    removeScaling();   // re-baseline onto the ORIGINAL masses

    if (ExplicitBathe::domainChanged() < 0)
        return -1;

    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = (theModel != 0) ? theModel->getDomainPtr() : 0;
    if (theModel == 0 || theDomain == 0) {
        opserr << "ExplicitBatheSMS::domainChanged - missing model/domain\n";
        return -1;
    }

    if (!warnedLimitations) {
        warnedLimitations = true;
        opserr << "ExplicitBatheSMS: v1 limitations -- (1) elements touching an MP-constraint "
                  "SLAVE node are EXCLUDED from scaling (they remain governing). (2) sizing "
                  "accounts for betaK Rayleigh damping (closed-form s=T^2+2*T*betaK/dt_e); alphaM "
                  "is not folded in. (3) in PARALLEL the injected lumped mass on a shared node IS "
                  "summed across ranks by a distributed/MPI diagonal solver (`system MPIDiagonal` "
                  "/ OpenSeesSP `system Diagonal`); the CONSISTENT (Olovsson) variant is ALSO "
                  "parallel-safe (ADR-38 V5): its matrix-free PCG uses GLOBAL inner products "
                  "(globalReduceSum) and shared-DOF assembly under `system MPIDiagonal` "
                  "(see LadrunoConsistentRefine.h).\n";
    }

    Ladruno::MassScalingReport rep =
        Ladruno::buildMassScaling(theModel, dtTarget, lumpingSMS, useTangentSMS, injected);
    Ladruno::applyMassScaling(theDomain, injected, +1.0);
    scaled = true;
    appliedDomain = theDomain;

    double frac = (rep.modelMass > 0.0) ? rep.addedMass / rep.modelMass : 0.0;
    bool over = (maxAddedMassFrac > 0.0 && frac > maxAddedMassFrac);
    if (verboseSMS || over || rep.nSelfReport > 0 || rep.nMismatch > 0 || rep.nConstrained > 0) {
        opserr << "ExplicitBatheSMS: dtTarget=" << dtTarget
               << " scaled " << rep.nScaled << "/" << rep.nElems
               << " elements; added mass " << (100.0 * frac) << "% of model mass"
               << "; PRE-SCALING dt_cr estimate=" << rep.minDtScaled
               << " (governing un-scaled element step; AFTER scaling the run is stable "
                  "at dt <= dtTarget=" << dtTarget << ")\n";
    }
    if (rep.nSelfReport > 0)
        opserr << "WARNING ExplicitBatheSMS: " << rep.nSelfReport
               << " throttling element(s) have a MASS-INDEPENDENT (self-reported) bound and "
                  "were NOT scaled -- lower dtTarget to satisfy them.\n";
    if (rep.nMismatch > 0)
        opserr << "WARNING ExplicitBatheSMS: " << rep.nMismatch
               << " element(s) had a non-node-major / DOF-mismatch mass and were SKIPPED.\n";
    if (rep.nConstrained > 0)
        opserr << "WARNING ExplicitBatheSMS: " << rep.nConstrained
               << " sub-target element(s) touch an MP-constrained (slave) node and were "
                  "EXCLUDED -- they still GOVERN at dt_e=" << rep.minDtConstrained
               << " < dtTarget=" << dtTarget << ".\n";
    if (over)
        opserr << "WARNING ExplicitBatheSMS: added mass " << (100.0 * frac)
               << "% exceeds -maxAddedMass cap " << (100.0 * maxAddedMassFrac)
               << "% (proceeding; the scaled inertia shifts global frequencies)\n";

    return 0;
}

int ExplicitBatheSMS::sendSelf(int cTag, Channel &theChannel)
{
    if (ExplicitBathe::sendSelf(cTag, theChannel) < 0)
        return -1;
    Vector data(5);
    data(0) = dtTarget;
    data(1) = maxAddedMassFrac;
    data(2) = verboseSMS ? 1.0 : 0.0;
    data(3) = (double)(int)lumpingSMS;
    data(4) = useTangentSMS ? 1.0 : 0.0;
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheSMS::sendSelf - could not send data\n";
        return -1;
    }
    return 0;
}

int ExplicitBatheSMS::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    if (ExplicitBathe::recvSelf(cTag, theChannel, theBroker) < 0)
        return -1;
    Vector data(5);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheSMS::recvSelf - could not receive data\n";
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
    injected.clear();
    scaled = false;
    appliedDomain = 0;
    warnedLimitations = false;
    return 0;
}

void ExplicitBatheSMS::Print(OPS_Stream &s, int flag)
{
    s << "ExplicitBatheSMS (Noh-Bathe + lumped selective mass scaling, classTag 33009)\n";
    s << "  dtTarget = " << dtTarget << ", -maxAddedMass = " << maxAddedMassFrac << "\n";
    s << "  recipe: lumped/diagonal mass, system Diagonal, algorithm Linear, dt <= dtTarget\n";
}
