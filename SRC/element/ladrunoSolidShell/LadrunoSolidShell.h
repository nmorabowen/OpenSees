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

// LadrunoSolidShell — 8-node ANS/EAS solid-shell element (Ladruno fork, ADR 66).
//
// A 3-translational-DOF-per-node solid-shell carrying the full triaxial stress
// state (sigma_33 != 0) through shell-like geometry: the narrow punching /
// bearing / 3D-crush specialist the director shell (ASDShellQ4 +
// LayeredShellFiberSection) structurally cannot be. NOT a co-equal flexural
// host (ADR 66 D1) — walls/slabs in bending stay on the director shell.
//
//   element('LadrunoSolidShell', tag, n1..n8, matTag
//           [, '-nz', n]                   # through-thickness points (default 2)
//           [, '-quadz', 'gauss'|'lobatto']# z-rule (default: gauss for nz=2,
//                                          #          lobatto for nz>=3)
//           [, '-formulation', 'ans'|'std']# default ans; std = plain
//                                          # displacement brick (locking A/B)
//           [, '-geom', 'linear'])         # corot/finite reserved (P5.4)
//
// Nodes 1-4 = bottom face (zeta=-1, counterclockwise), 5-8 = top face; the
// thickness direction is DEFINED by the node ordering (the 1-5 .. 4-8 fibers).
//
// Formulation (covariant, small strain; ADR 66 §5) — 'ans' applies, in order:
//  1. ANS transverse shear (Dvorkin-Bathe / MITC4 tying): covariant gamma_13
//     sampled at (0,-+1,zeta) and interpolated linearly in eta; gamma_23
//     sampled at (-+1,0,zeta), linear in xi. Cures transverse-shear locking.
//  2. ANS thickness strain (Betsch-Stein): covariant E_33 bilinearly
//     interpolated from the four midsurface corners (+-1,+-1,0). Cures
//     curvature-thickness (trapezoidal) locking.
//  3. EAS on E_33 (Simo-Rifai, adapted from LadrunoBrick::formEAStrue): one
//     enhanced mode Etilde_33 = zeta*alpha mapped through the CENTER transform,
//     eps_enh = (j0/j(xi)) * zeta * T0(:,3) * alpha. Cures Poisson-thickness
//     locking. alpha is solved each form pass by an inner Newton on
//     int G^T sigma dV = 0 against the LIVE material stress, statically
//     condensed (K* = Kdd - Kda Kaa^-1 Kad), committed (alphaCommit) and
//     serialized — state-dependent EAS from day one (ADR 66 §3); the
//     softening-stability gate (G6) is the P5.2 deliverable.
//
// The covariant strain E_ij = sym(g_i . u_,j) is transformed to Cartesian
// engineering Voigt {xx,yy,zz,xy,yz,zx} per GP via the contravariant-basis
// 6x6 map T (eps_ab = E_ij g^i_a g^j_b), so ANY 3D NDMaterial drops in
// unchanged (LadrunoConcrete3D is the headline host).
//
// Integration: 2x2 in-plane Gauss x selectable nz through thickness
// (-quadz gauss nz in [2,5] | lobatto nz in [3,7]; Lobatto puts points on the
// faces for surface crack initiation). One material state per GP (4*nz).
//
// getCharacteristicLength() returns the IN-PLANE projected size
// sqrt(midsurface area) — never the thickness — so crack-band materials
// regularize over the dominant in-plane band (ADR 66 D6). Through-thickness
// bands regularize over t/n_elem_z: mesh guidance, not code, in v1.
//
// v1 scope: -geom linear only (the SolidTransformation seam is wired at P5.4
// together with the cond(S)-guarded corot; parser rejects corot/finite).
// No element loads (surface pressure = P5.3); no per-element damping.  // Ladruno

#ifndef LadrunoSolidShell_h
#define LadrunoSolidShell_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <LadrunoMassCache.h>   // Ladruno (ADR-77 G2 ext)
#include <Node.h>
#include <NDMaterial.h>

class Response;

class LadrunoSolidShell : public Element {

 public:
     // Ladruno (ADR-77 G2 ext): escape = -noMassCache
     void setMassCache(bool s) { massCache.setEnabled(s); }


  // Formulation selector. Ordinals are serialized — order is load-bearing.
  //   ANS — the full solid-shell recipe (ANS shear + ANS E33 + EAS-on-E33).
  //   STD — plain covariant displacement brick at the same quadrature: the
  //         locking NEGATIVE CONTROL for the G2/G4 gates (no ANS, no EAS).
  enum class Formulation { ANS, STD };

  // Through-thickness quadrature rule. Ordinals serialized.
  enum class QuadZ { GAUSS, LOBATTO };

  static const int NEAS = 1;      // enhanced E33 modes (slots reserved for the
                                  // bilinear-in-plane extension, ADR 66 O1)
  static const int NZ_MIN = 2;
  static const int NZ_MAX = 7;    // gauss caps at 5; lobatto at 7 (zRule)

  // through-thickness rule tables (public: the OPS parser validates the
  // (rule, nz) combination up front with the same tables the element uses)
  static int zRuleStatic(QuadZ q, int n, double *pts, double *wts);

  // null constructor (broker)
  LadrunoSolidShell();

  // full constructor
  LadrunoSolidShell(int tag,
                    int node1, int node2, int node3, int node4,
                    int node5, int node6, int node7, int node8,
                    NDMaterial &theMaterial,
                    int nzPts = 2,
                    QuadZ zQuad = QuadZ::GAUSS,
                    Formulation formulation = Formulation::ANS,
                    int matype = 0);

  virtual ~LadrunoSolidShell();

  const char *getClassType(void) const { return "LadrunoSolidShell"; };

  // domain
  void setDomain(Domain *theDomain);

  // topology
  int getNumExternalNodes(void) const;
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  int getNumDOF(void);

  // In-plane projected characteristic length sqrt(midsurface area) — see the
  // file header (ADR 66 D6). Never the thickness.  // Ladruno
  double getCharacteristicLength(void);

  // state
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  void Print(OPS_Stream &s, int flag);

  // matrices / vectors
  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);

  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);

  // parallel / database
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  // output
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);

 private:

  // -------- attributes --------
  ID connectedExternalNodes;        // 8 node tags
  Node *nodePointers[8];
  LadrunoMassCache massCache;   // Ladruno (ADR-77 G2 ext): per-instance mass cache,
                  // guard-checked (rho/thickness/coords); -noMassCache escape

  NDMaterial **materialPointers;    // 4*nz material states (allocated in ctor/recvSelf)

  Formulation formulation;
  QuadZ quadz;
  int nz;                           // through-thickness point count
  int numGP;                        // = 4*nz
  int massType;                     // 0 consistent (default), 1 row-sum lumped
                                    // (== HRZ on the trilinear shape; the
                                    // explicit-path mass, ADR 66 G9)  // Ladruno

  Vector alpha;                     // NEAS enhanced parameters (trial)
  Vector alphaCommit;               // committed enhanced parameters (serialized)

  // center enhanced-mode map (Simo-Rifai): 3rd column of the center covariant->
  // Cartesian transform T0, and the center Jacobian determinant j0. Geometry-only
  // — rebuilt in setDomain (nothing extra travels in sendSelf besides alphaCommit).
  double t0col3[6];
  double j0det;
  bool   easBuilt;

  double X[8][3];                   // reference nodal coordinates (setDomain)

  Vector *load;
  Matrix *Ki;

  // -------- static workspace --------
  static Matrix stiff;
  static Vector resid;
  static Matrix mass;

  static const double natCoord[8][3];   // node natural coordinates (brick order)

  // -------- private methods --------
  int  nGP(void) const { return numGP; }

  // through-thickness rule: fills pts/wts (size nz), returns 0 or -1 (bad combo)
  int zRule(double *pts, double *wts) const;

  // trilinear shape functions + natural gradients at xi=(xi,eta,zeta)
  static void shapeFcn(const double xi[3], double N[8], double dN[8][3]);
  // covariant base J[i][j] = dx_i/dxi_j from reference coords; returns detJ
  static double jacobian(const double X8[8][3], const double dN[8][3], double J[3][3]);
  static void invert3(const double J[3][3], double det, double Jinv[3][3]);
  // 6x6 covariant->Cartesian strain map (engineering Voigt both sides)
  static void buildT(const double Jinv[3][3], double T[6][6]);

  // one covariant B row (engineering shear) at natural point xi: out[24]
  void covRow(int row, const double xi[3], double out[24]) const;
  // full 6x24 covariant B + J data at xi (no ANS substitution)
  void covB(const double xi[3], double Bcov[6][24], double J[3][3],
            double &detJ) const;
  // Cartesian 6x24 B at a GP: covariant B, ANS-substituted rows (formulation==
  // ANS), transformed by T(Jinv). Also returns detJ.
  void cartB(const double xi[3], Matrix &B, double &detJ) const;

  void buildEAS(void);              // cache t0col3 / j0det (center map)
  void formANS(int tang_flag, bool useInitialTangent);  // inner Newton + condensation
  void formInertiaTerms(int tangFlag);
};

#endif
