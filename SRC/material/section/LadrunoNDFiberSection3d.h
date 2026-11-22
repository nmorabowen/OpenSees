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

// Ladruno WP-245 (program #1131, rcbf-3): an NDFiber 3D section with a chosen
// shear-strain profile across the section.
//
//   section LadrunoNDFiber $tag -shearProfile uniform   [-alpha $a] [-noCentroid] { fibres }
//   section LadrunoNDFiber $tag -shearProfile parabolic -shearDepth $hy $hz
//                               [-shearResultant modelC|conjugate] [-noCentroid] { fibres }
//
// It is a thin subclass of the stock NDFiberSection3d. Fibre storage, addFiber,
// the centroid and the fibre/patch/layer builders are the stock ones; only the
// fibre kinematics, the resultants and the tangent differ. The shear map is in
// LadrunoFiberShearKernel.h; read its header for the two resultant rules (and
// why neither gives Timoshenko's 5/6).
//
// Fail-closed:
//   * -shearProfile is required: no silent default.
//   * a fibre outside the parabolic band (|y| > hy/2 or |z| > hz/2 from the
//     centroid) is refused at the first state determination.
//   * stock DDM sensitivities use the stock sqrt(alpha) kinematics, so they are
//     refused under the parabolic profile.
// The modelC resultant gives an unsymmetric tangent: use an unsymmetric solver
// (plain `system Pardiso`, UmfPack, Mumps); not Pardiso -symmetric / -spd.

#ifndef LadrunoNDFiberSection3d_h
#define LadrunoNDFiberSection3d_h

#include <NDFiberSection3d.h>
#include "LadrunoFiberShearKernel.h"

class LadrunoNDFiberSection3d : public NDFiberSection3d
{
  public:
    LadrunoNDFiberSection3d();
    LadrunoNDFiberSection3d(int tag, int sizeHint, const ladruno_shear::Spec& spec, bool compCentroid);
    ~LadrunoNDFiberSection3d() {}

    const char *getClassType(void) const {return "LadrunoNDFiberSection3d";};

    // Reads the options after the tag from the OPS input stream (consumes them
    // all). Returns 0, or -1 after a named refusal. Shared by every interpreter.
    static int parseOptions(int tag, ladruno_shear::Spec& spec, bool& compCentroid);

    int   setTrialSectionDeformation(const Vector &deforms);
    const Matrix &getInitialTangent(void);
    int   revertToLastCommit(void);
    int   revertToStart(void);

    SectionForceDeformation *getCopy(void);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

    int setParameter(const char **argv, int argc, Parameter &param);
    const Vector& getStressResultantSensitivity(int gradIndex, bool conditional);
    const Matrix& getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    const ladruno_shear::Spec& getShearSpec(void) const { return spec; }

  private:
    // Rebuild s and ks from the fibres' current stress and tangent (after a revert).
    void assembleFromMaterials(void);
    bool checkGeometry(void);

    ladruno_shear::Spec spec;
    int geomState;   // 0 unchecked, 1 ok, -1 refused
};

#endif
