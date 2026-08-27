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

// Ladruno: LadrunoSANISAND -- Manzari-Dafalias (2004) with the two low-stress
// constants made settable, carried on the wire, and echoed at construction.
//
// `ManzariDafalias::initialize()` assigns, unconditionally and with no way to
// override:
//
//     m_Pmin      = 1.0e-4 * m_P_atm;   //  0.0101 kPa at P_atm = 101
//     m_Presidual = 1.0e-2 * m_P_atm;   //  1.01   kPa at P_atm = 101
//
// `m_Presidual` is threaded through ~30 mean-stress evaluations as
// `p = one3*GetTrace(stress) + m_Presidual`, including the one that feeds
// `GetPSI()`, so it displaces the state parameter psi and with it M^b, M^d and
// the dilatancy. It acts as an apparent cohesion c = p_r*tan(phi) ~ 0.95 kPa:
// negligible at 100 kPa, +20 % on the model's own bounding-surface identity at
// 1 kPa. It is also NOT on the vanilla wire (`Vector(97)` carries `m_Pmin` at
// data(96) and nothing for `m_Presidual`), so an OpenSeesMP worker or a
// database-restored material runs at p_r = 0 while the serial process beside it
// runs at 1.01 kPa.
//
// This class changes exactly two constants and nothing else. Defaults follow
// NTUASand02 (Gorini): p_residual = 0 (a cohesionless sand has no cohesion) and
// p_min = 1.0e-3 * P_atm.
//
// DESIGN NOTE (do not "clean up"):
//   `m_Presidual` and `m_Pmin` are PROTECTED DATA read at run time by every base
//   integrator, so this subclass never has to override an integrator -- it only
//   has to win the LAST WRITE. Every method it must control for that is already
//   virtual through NDMaterial/MovableObject: revertToStart, sendSelf, recvSelf,
//   both getCopy forms, Print.
//
//   `initialize()` is NOT virtual and MUST NOT BECOME VIRTUAL. The declaration
//   below is deliberately a *shadow*: the base constructor's `initialize()` call
//   still binds to the base version (correct -- the base must run its own init),
//   and the derived constructor re-applies the two constants afterwards.
//   `ManzariDafaliasRO` (UWmaterials/ManzariDafaliasRO.h:93-97) shadows
//   `initialize()` and two `GetElasticModuli` overloads with identical
//   signatures; adding `virtual` to any of those in the base would silently turn
//   those shadows into overrides and make every existing RO deck run
//   Ramberg-Osgood elasticity from inside the base integrators.
//
// Vanilla `ManzariDafalias` is not edited.
// See Ladruno_implementation/86_ladruno_sanisand_adr.md.
// classTags 33019 (base) / 33020 (3D) / 33021 (PlaneStrain).
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoSANISAND_h
#define LadrunoSANISAND_h

// Quote-include so it resolves against this file's own directory first
// (SRC/material/nD/UWmaterials/) and against SRC/material/nD/ second -- the
// same spelling FEM_ObjectBrokerAllClasses.cpp uses.
#include "UWmaterials/ManzariDafalias.h"

class LadrunoSANISAND : public ManzariDafalias
{
  public:

    // full constructor, explicit classTag -- used by the parser and by the
    // LadrunoSANISAND3D / LadrunoSANISANDPlaneStrain wrappers.
    // Defaults of the five optional integration args match the base's
    // equivalent (int tag, int classTag, ...) form: 1, 0, 1, 1e-7, 1e-7.
    LadrunoSANISAND(int tag, int classTag, double G0, double nu, double e_init, double Mc, double c, double lambda_c,
                    double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0,
                    double nd, double z_max, double cz, double mDen,
                    int integrationScheme = 1, int tangentType = 0, int JacoType = 1,
                    double TolF = 1.0e-7, double TolR = 1.0e-7,
                    double Presidual = 0.0, double Pmin = -1.0, int honorTolR = 0);   // Ladruno

    // full constructor, classTag defaults to ND_TAG_LadrunoSANISAND.
    // Defaults of the five optional integration args match the base's
    // equivalent (int tag, ...) form: 2, 2, 1, 1e-7, 1e-7.
    LadrunoSANISAND(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c,
                    double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0,
                    double nd, double z_max, double cz, double mDen,
                    int integrationScheme = 2, int tangentType = 2, int JacoType = 1,
                    double TolF = 1.0e-7, double TolR = 1.0e-7,
                    double Presidual = 0.0, double Pmin = -1.0, int honorTolR = 0);   // Ladruno

    // specific-type null constructor (used by the wrappers' null constructors)
    LadrunoSANISAND(int classTag);

    // null constructor (FEM_ObjectBroker)
    LadrunoSANISAND();

    // destructor
    ~LadrunoSANISAND();

    // --- the six overrides. See ADR 86 section 4.4; none is optional. -------

    // Replicates the base ops_InitialStateAnalysis guard but re-initialises
    // through THIS class's initialize(), so p_residual is not restored to
    // 1.0e-2*P_atm mid-analysis.
    int revertToStart(void);

    // The base getCopy(const char*) hardcodes `new ManzariDafalias3D(...)` /
    // `new ManzariDafaliasPlaneStrain(...)`, which would strip these settings at
    // every Gauss point.
    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *type);

    // Base Vector(97) wire format unchanged; one extra Vector(4) follows it.
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

  protected:

    double mPresidualInput;   // Ladruno: p_residual as given by the user (>= 0)
    double mPminInput;        // Ladruno: p_min as given; < 0 is the SENTINEL for
                              //          "resolve to 1.0e-3 * P_atm" (P_atm is not
                              //          known until the base ctor has run)
    int    mHonorTolR;        // Ladruno: PR-2 seam placeholder. Always 0 in PR-1;
                              //          any other value is rejected at parse time.

    // SHADOW of the non-virtual ManzariDafalias::initialize(). Same signature on
    // purpose -- see the DESIGN NOTE above. DO NOT add `virtual` here or in the
    // base.
    void initialize();

    // The single "win the last write" helper. Called from every constructor
    // (after the base ctor has returned), from revertToStart via initialize(),
    // and from recvSelf.
    void applyLadrunoConstants(void);   // Ladruno

    // Per-construction echo (NOT latched -- see ADR 86 section 4.4).
    void echoLadrunoConstants(void);    // Ladruno
};

#endif
