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
// See ExplicitBatheLNVDSMS.h and Ladruno_implementation/36_ladruno_selective_mass_scaling_adr.md.

#include <ExplicitBatheLNVDSMS.h>
#include <LadrunoMassScaling.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <Channel.h>
#include <Vector.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

#include <string.h>

void *OPS_ExplicitBatheLNVDSMS(void)
{
    // Usage: integrator ExplicitBatheLNVDSMS $p $alpha $dtTarget <-maxAddedMass f>
    //                   <-lump rowsum|diagonal|hrz> <-tangent> <-verbose>
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING integrator ExplicitBatheLNVDSMS $p $alpha $dtTarget <options> "
                  "- needs the Noh-Bathe p, FLAC alpha, and a target time step\n";
        return 0;
    }

    double in[3] = {0.0, 0.0, 0.0};
    int n3 = 3;
    if (OPS_GetDoubleInput(&n3, in) < 0) {
        opserr << "WARNING ExplicitBatheLNVDSMS - could not read $p $alpha $dtTarget\n";
        return 0;
    }
    double p = in[0], alpha = in[1], dtTarget = in[2];
    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING ExplicitBatheLNVDSMS - $p must be in (0,1)\n";
        return 0;
    }
    if (dtTarget <= 0.0) {
        opserr << "WARNING ExplicitBatheLNVDSMS - $dtTarget must be a positive double\n";
        return 0;
    }

    double maxAddedMassFrac = 0.05;
    bool   verboseSMS = false;
    bool   cflUseTangent = false;
    int    compute_critical_timestep = 1;
    CTSLumping lumping = CTSLumping::Diagonal;

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
                else opserr << "WARNING ExplicitBatheLNVDSMS - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping diagonal)\n";
            }
        } else if (strcmp(arg, "-cflAbort") == 0 || strcmp(arg, "-recompute") == 0) {
            opserr << "WARNING ExplicitBatheLNVDSMS - " << arg << " is not supported with mass "
                      "scaling (it re-runs the un-augmented element eigensolve). NOT created.\n";
            return 0;
        } else {
            opserr << "WARNING ExplicitBatheLNVDSMS - unknown option " << arg << " (ignored)\n";
        }
    }

    TransientIntegrator *theIntegrator =
        new ExplicitBatheLNVDSMS(p, alpha, dtTarget, maxAddedMassFrac, verboseSMS,
                                 compute_critical_timestep, cflUseTangent, lumping);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitBatheLNVDSMS integrator\n";
    return theIntegrator;
}

// Default constructor (broker)
ExplicitBatheLNVDSMS::ExplicitBatheLNVDSMS()
    : ExplicitBatheLNVD(INTEGRATOR_TAGS_ExplicitBatheLNVDSMS, 0.0, 0.0, 0, false, false, 0.0,
                        false, 0, CTSLumping::Diagonal),
      dtTarget(0.0), maxAddedMassFrac(0.05), verboseSMS(false),
      lumpingSMS(CTSLumping::Diagonal), useTangentSMS(false), scaled(false),
      appliedDomain(0), warnedLimitations(false)
{
}

// Main constructor
ExplicitBatheLNVDSMS::ExplicitBatheLNVDSMS(double p_, double alpha_, double dtTarget_,
                                           double maxAddedMassFrac_, bool verboseSMS_,
                                           int compute_critical_timestep_,
                                           bool cflUseTangent_, CTSLumping lumping_)
    : ExplicitBatheLNVD(INTEGRATOR_TAGS_ExplicitBatheLNVDSMS, p_, alpha_,
                        compute_critical_timestep_, verboseSMS_, /*cflAbort*/false,
                        /*divergenceFactor*/0.0, cflUseTangent_, /*cflRecomputeEvery*/0,
                        lumping_),
      dtTarget(dtTarget_), maxAddedMassFrac(maxAddedMassFrac_), verboseSMS(verboseSMS_),
      lumpingSMS(lumping_), useTangentSMS(cflUseTangent_), scaled(false),
      appliedDomain(0), warnedLimitations(false)
{
}

ExplicitBatheLNVDSMS::~ExplicitBatheLNVDSMS()
{
    if (scaled && appliedDomain != 0 && !injected.empty())
        Ladruno::applyMassScaling(appliedDomain, injected, -1.0);
}

void ExplicitBatheLNVDSMS::removeScaling(void)
{
    if (scaled && appliedDomain != 0 && !injected.empty())
        Ladruno::applyMassScaling(appliedDomain, injected, -1.0);
    injected.clear();
    scaled = false;
    appliedDomain = 0;
}

int ExplicitBatheLNVDSMS::domainChanged(void)
{
    removeScaling();

    if (ExplicitBatheLNVD::domainChanged() < 0)
        return -1;

    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = (theModel != 0) ? theModel->getDomainPtr() : 0;
    if (theModel == 0 || theDomain == 0) {
        opserr << "ExplicitBatheLNVDSMS::domainChanged - missing model/domain\n";
        return -1;
    }

    if (!warnedLimitations) {
        warnedLimitations = true;
        opserr << "ExplicitBatheLNVDSMS: v1 limitations -- (1) elements touching an "
                  "MP-constraint SLAVE node are EXCLUDED. (2) sizing accounts for betaK "
                  "damping (closed-form s=T^2+2*T*betaK/dt_e); alphaM not folded in. NOTE the "
                  "FLAC local damping (alpha) is separate from sizing. (3) parallel "
                  "shared/boundary nodes are not reduced across ranks.\n";
    }

    Ladruno::MassScalingReport rep =
        Ladruno::buildMassScaling(theModel, dtTarget, lumpingSMS, useTangentSMS, injected);
    Ladruno::applyMassScaling(theDomain, injected, +1.0);
    scaled = true;
    appliedDomain = theDomain;

    double frac = (rep.modelMass > 0.0) ? rep.addedMass / rep.modelMass : 0.0;
    bool over = (maxAddedMassFrac > 0.0 && frac > maxAddedMassFrac);
    if (verboseSMS || over || rep.nSelfReport > 0 || rep.nMismatch > 0 || rep.nConstrained > 0) {
        opserr << "ExplicitBatheLNVDSMS: dtTarget=" << dtTarget
               << " scaled " << rep.nScaled << "/" << rep.nElems
               << " elements; added mass " << (100.0 * frac) << "% of model mass"
               << " (governing un-scaled dt_e=" << rep.minDtScaled << ")\n";
    }
    if (rep.nSelfReport > 0)
        opserr << "WARNING ExplicitBatheLNVDSMS: " << rep.nSelfReport
               << " throttling element(s) have a MASS-INDEPENDENT bound and were NOT scaled.\n";
    if (rep.nMismatch > 0)
        opserr << "WARNING ExplicitBatheLNVDSMS: " << rep.nMismatch
               << " element(s) had a non-node-major / DOF-mismatch mass and were SKIPPED.\n";
    if (rep.nConstrained > 0)
        opserr << "WARNING ExplicitBatheLNVDSMS: " << rep.nConstrained
               << " sub-target element(s) touch an MP-constrained (slave) node and were "
                  "EXCLUDED -- they still GOVERN at dt_e=" << rep.minDtConstrained << ".\n";
    if (over)
        opserr << "WARNING ExplicitBatheLNVDSMS: added mass " << (100.0 * frac)
               << "% exceeds -maxAddedMass cap " << (100.0 * maxAddedMassFrac) << "%\n";

    return 0;
}

int ExplicitBatheLNVDSMS::sendSelf(int cTag, Channel &theChannel)
{
    if (ExplicitBatheLNVD::sendSelf(cTag, theChannel) < 0)
        return -1;
    Vector data(5);
    data(0) = dtTarget;
    data(1) = maxAddedMassFrac;
    data(2) = verboseSMS ? 1.0 : 0.0;
    data(3) = (double)(int)lumpingSMS;
    data(4) = useTangentSMS ? 1.0 : 0.0;
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheLNVDSMS::sendSelf - could not send data\n";
        return -1;
    }
    return 0;
}

int ExplicitBatheLNVDSMS::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    if (ExplicitBatheLNVD::recvSelf(cTag, theChannel, theBroker) < 0)
        return -1;
    Vector data(5);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheLNVDSMS::recvSelf - could not receive data\n";
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

void ExplicitBatheLNVDSMS::Print(OPS_Stream &s, int flag)
{
    s << "ExplicitBatheLNVDSMS (Noh-Bathe + FLAC LNVD + lumped mass scaling, classTag 33011)\n";
    s << "  dtTarget = " << dtTarget << ", -maxAddedMass = " << maxAddedMassFrac << "\n";
}
