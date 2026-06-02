/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// SolidTransformationCorot — corotational geometry method (v2). See the header
// for the full theory and the v2.0 tangent scope.

#include <SolidTransformationCorot.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
#include <OPS_Globals.h>   // opserr

// =========================================================================== //
//  File-local 3x3 linear-algebra on row-major double[9] (M(i,j) = M[3*i+j]).
//  Kept allocation-free and local; the corot wrapper is the only consumer.
// =========================================================================== //
namespace {

inline void mat3_zero(double A[9]) { for (int i = 0; i < 9; i++) A[i] = 0.0; }
inline void mat3_iden(double A[9]) { mat3_zero(A); A[0] = A[4] = A[8] = 1.0; }

inline double det3(const double A[9])
{
  return A[0]*(A[4]*A[8] - A[5]*A[7])
       - A[1]*(A[3]*A[8] - A[5]*A[6])
       + A[2]*(A[3]*A[7] - A[4]*A[6]);
}

// inverse of a general 3x3; returns det (0 ⇒ singular, Ainv untouched).
inline double inv3(const double A[9], double Ainv[9])
{
  double d = det3(A);
  if (d == 0.0) return 0.0;
  double di = 1.0 / d;
  Ainv[0] =  (A[4]*A[8] - A[5]*A[7]) * di;
  Ainv[1] = -(A[1]*A[8] - A[2]*A[7]) * di;
  Ainv[2] =  (A[1]*A[5] - A[2]*A[4]) * di;
  Ainv[3] = -(A[3]*A[8] - A[5]*A[6]) * di;
  Ainv[4] =  (A[0]*A[8] - A[2]*A[6]) * di;
  Ainv[5] = -(A[0]*A[5] - A[2]*A[3]) * di;
  Ainv[6] =  (A[3]*A[7] - A[4]*A[6]) * di;
  Ainv[7] = -(A[0]*A[7] - A[1]*A[6]) * di;
  Ainv[8] =  (A[0]*A[4] - A[1]*A[3]) * di;
  return d;
}

// C = A * B
inline void matmul3(const double A[9], const double B[9], double C[9])
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      C[3*i+j] = A[3*i]*B[j] + A[3*i+1]*B[3+j] + A[3*i+2]*B[6+j];
}

// C = Aᵀ * B
inline void matTmul3(const double A[9], const double B[9], double C[9])
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      C[3*i+j] = A[i]*B[j] + A[3+i]*B[3+j] + A[6+i]*B[6+j];
}

inline void transpose3(const double A[9], double At[9])
{
  At[0]=A[0]; At[1]=A[3]; At[2]=A[6];
  At[3]=A[1]; At[4]=A[4]; At[5]=A[7];
  At[6]=A[2]; At[7]=A[5]; At[8]=A[8];
}

// y = A * x  (3-vector)
inline void matvec3(const double A[9], const double x[3], double y[3])
{
  y[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2];
  y[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2];
  y[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2];
}

// y = Aᵀ * x  (3-vector)
inline void matTvec3(const double A[9], const double x[3], double y[3])
{
  y[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
  y[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
  y[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
}

// axial vector of skew(M) = ½(M − Mᵀ):  v = ½[M21−M12, M02−M20, M10−M01]
inline void axialSkew(const double M[9], double v[3])
{
  v[0] = 0.5 * (M[7] - M[5]);
  v[1] = 0.5 * (M[2] - M[6]);
  v[2] = 0.5 * (M[3] - M[1]);
}

// skew (cross-product) matrix of a vector: spin(p) w = p × w
inline void spin3(const double p[3], double S[9])
{
  S[0] = 0.0;   S[1] = -p[2]; S[2] =  p[1];
  S[3] = p[2];  S[4] = 0.0;   S[5] = -p[0];
  S[6] = -p[1]; S[7] = p[0];  S[8] = 0.0;
}

inline double frob3(const double A[9])
{
  double s = 0.0; for (int i = 0; i < 9; i++) s += A[i]*A[i]; return sqrt(s);
}

// Polar rotation R = polar(H) by scaled Higham Newton iteration:
//   X_{k+1} = ½ ( γ_k X_k + γ_k⁻¹ X_k^{−ᵀ} ),  γ_k = ( ‖X⁻¹‖_F / ‖X‖_F )^{½}.
// Converges quadratically to the orthogonal polar factor for any det(H) > 0.
// Returns 0 on success; -1 if H is singular / det ≤ 0 (collapsed/inverted).
inline int polarRotation(const double H[9], double R[9])
{
  double X[9]; for (int i = 0; i < 9; i++) X[i] = H[i];
  if (det3(X) <= 0.0) return -1;            // reflection or collapse — reject

  const int    maxIter = 50;
  const double tol      = 1.0e-15;
  double Xinv[9], XinvT[9], Xnew[9];

  for (int it = 0; it < maxIter; it++) {
    if (inv3(X, Xinv) == 0.0) return -1;
    transpose3(Xinv, XinvT);
    double g = sqrt(frob3(Xinv) / frob3(X));
    if (!(g > 0.0)) g = 1.0;
    double diff = 0.0;
    for (int i = 0; i < 9; i++) {
      Xnew[i] = 0.5 * (g * X[i] + (1.0/g) * XinvT[i]);
      double d = Xnew[i] - X[i];
      diff += d*d;
    }
    for (int i = 0; i < 9; i++) X[i] = Xnew[i];
    if (sqrt(diff) < tol) break;
  }
  for (int i = 0; i < 9; i++) R[i] = X[i];
  return (det3(R) > 0.0) ? 0 : -1;
}

} // anonymous namespace

// =========================================================================== //

SolidTransformationCorot::SolidTransformationCorot()
  : SolidTransformation(), numNodes(0)
{
  mat3_iden(Rmat);
  for (int a = 0; a < MAX_NODES; a++)
    for (int d = 0; d < 3; d++) { xrel[a][d] = 0.0; Xrel[a][d] = 0.0; }
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < MAX_DOF; j++) Gamma[i][j] = 0.0;
}

SolidTransformationCorot::~SolidTransformationCorot()
{
}

SolidTransformation *
SolidTransformationCorot::getCopy() const
{
  return new SolidTransformationCorot();
}

// --------------------------------------------------------------------------- //
// (seam 0) refresh the corotational frame: extract R = polar(H), the centred
// coordinates, and the spin map Γ = ∂θ/∂u. Everything the const seam methods
// need is computed and cached here, once per residual/tangent evaluation.
// --------------------------------------------------------------------------- //
int
SolidTransformationCorot::update(int nn, const Matrix &refCrds,
                                 const Matrix &curCrds)
{
  if (nn < 1 || nn > MAX_NODES) {
    opserr << "SolidTransformationCorot::update - numNodes " << nn
           << " out of range [1, " << (int)MAX_NODES << "]\n";
    return -1;
  }
  numNodes = nn;

  // centroids
  double Xc[3] = {0,0,0}, xc[3] = {0,0,0};
  for (int a = 0; a < nn; a++)
    for (int d = 0; d < 3; d++) {
      Xc[d] += refCrds(a, d);
      xc[d] += curCrds(a, d);
    }
  for (int d = 0; d < 3; d++) { Xc[d] /= nn; xc[d] /= nn; }

  // centred coordinates and Procrustes cross-covariance H = Σ x_a⁰ (X_a⁰)ᵀ
  double H[9]; mat3_zero(H);
  for (int a = 0; a < nn; a++) {
    for (int d = 0; d < 3; d++) {
      Xrel[a][d] = refCrds(a, d) - Xc[d];
      xrel[a][d] = curCrds(a, d) - xc[d];
    }
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        H[3*i+j] += xrel[a][i] * Xrel[a][j];
  }

  // R = polar(H)
  if (polarRotation(H, Rmat) != 0) {
    opserr << "SolidTransformationCorot::update - polar decomposition failed "
              "(collapsed / inverted element: det(H) <= 0)\n";
    // leave a SAFE state so a consumer that ignores this return value (e.g.
    // computeLocalDisp) degrades to an identity-rotation response rather than
    // running on stale/garbage R, Γ.  // Ladruno (review M2)
    mat3_iden(Rmat);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < MAX_DOF; j++) Gamma[i][j] = 0.0;
    return -1;
  }

  // stretch S = sym(Rᵀ H), and A_S = tr(S) I − S, then A_S⁻¹
  double S[9]; matTmul3(Rmat, H, S);
  for (int i = 0; i < 3; i++)
    for (int j = i+1; j < 3; j++) {
      double avg = 0.5*(S[3*i+j] + S[3*j+i]);
      S[3*i+j] = S[3*j+i] = avg;
    }
  double trS = S[0] + S[4] + S[8];
  double AS[9];
  for (int i = 0; i < 9; i++) AS[i] = -S[i];
  AS[0] += trS; AS[4] += trS; AS[8] += trS;
  double ASinv[9];
  if (inv3(AS, ASinv) == 0.0) {
    opserr << "SolidTransformationCorot::update - (tr(S) I - S) singular\n";
    // degrade to a MATCHED identity gradient/Hessian pair (R=I, Γ=0): with
    // Γ=0 the force loses its −Γᵀm term, so to keep f and K a consistent pair
    // we also drop R (mirrors the polar-failure branch). Practically
    // unreachable once det(H)>0 passed, but keeps Newton from chasing a
    // non-gradient force.  // Ladruno (sweep #13)
    mat3_iden(Rmat);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < MAX_DOF; j++) Gamma[i][j] = 0.0;
    return -1;
  }

  // Spin map Γ[:,3b+k] = A_S⁻¹ · ( 2 · axial( skew( Rᵀ δH ) ) ),
  //   δH = e_k ⊗ X_b⁰  ⇒  (Rᵀ δH)(i,j) = Rᵀ(i,k) X_b⁰(j) = R[3k+i] X_b⁰(j).
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < MAX_DOF; j++) Gamma[i][j] = 0.0;

  for (int b = 0; b < nn; b++) {
    for (int k = 0; k < 3; k++) {
      double M[9];
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          M[3*i+j] = Rmat[3*k+i] * Xrel[b][j];   // Rᵀ(i,k) X_b⁰(j)
      double v[3]; axialSkew(M, v);
      double rhs[3] = {2.0*v[0], 2.0*v[1], 2.0*v[2]};
      double w[3]; matvec3(ASinv, rhs, w);        // BODY spin ω = A_S⁻¹(2 axial(skew))
      // map to the SPATIAL spin used by the force/stiffness seams
      // (δR = spin(δθ) R ⇒ δθ = R ω). Without this R· the geometric terms are
      // in the wrong frame (caught by the standalone Γ-vs-FD gate).  // Ladruno
      double ws[3]; matvec3(Rmat, w, ws);
      for (int i = 0; i < 3; i++)
        Gamma[i][3*b + k] = ws[i];
    }
  }
  return 0;
}

// --------------------------------------------------------------------------- //
// (seam 2) deformational displacement u_d,a = Rᵀ(x_a − x_c) − (X_a − X_c).
// uGlobal is unused (the centred coordinates were cached in update()).
// --------------------------------------------------------------------------- //
int
SolidTransformationCorot::localizeDisp(const Vector & /*uGlobal*/,
                                       Vector &uCore) const
{
  for (int a = 0; a < numNodes; a++) {
    double ud[3]; matTvec3(Rmat, xrel[a], ud);     // Rᵀ x_a⁰
    for (int d = 0; d < 3; d++)
      uCore(3*a + d) = ud[d] - Xrel[a][d];
  }
  return 0;
}

// --------------------------------------------------------------------------- //
// (seam 3) internal force to global (exact gradient of the corotated energy):
//   f_global,b = R f_d,b − (1/n) R Σ_a f_d,a − Γ_{:,b}ᵀ m ,
//   m = Σ_a (x_a − x_c) × (R f_d,a).
// The −(1/n)RΣf_d centroid back-reaction is the gradient of u_d's centroid
// dependence; it vanishes for the pure internal force (Σ B ᵀσ = 0, self-
// equilibrated) but is NON-zero once the brick folds a body/applied load into
// fCore before this call — so keep it (sweep #9).
// fGlobal MAY alias fCore — read the core force into a local buffer first.
// --------------------------------------------------------------------------- //
int
SolidTransformationCorot::globalizeForce(const Vector &fCore,
                                         Vector &fGlobal) const
{
  const int n = numNodes;
  if (n == 0) {                       // queried before first update(): identity
    if (&fGlobal != &fCore) fGlobal = fCore;
    return 0;
  }

  // local copy of the core nodal forces, the rotated forces R f_d,a, the
  // moment m, and the net rotated force fSum (for the centroid term)
  double Rf[MAX_NODES][3];
  double m[3] = {0,0,0};
  double fSum[3] = {0,0,0};
  for (int a = 0; a < n; a++) {
    double fd[3] = { fCore(3*a), fCore(3*a+1), fCore(3*a+2) };
    matvec3(Rmat, fd, Rf[a]);                       // R f_d,a
    for (int d = 0; d < 3; d++) fSum[d] += Rf[a][d];
    // m += x_a⁰ × (R f_d,a)
    m[0] += xrel[a][1]*Rf[a][2] - xrel[a][2]*Rf[a][1];
    m[1] += xrel[a][2]*Rf[a][0] - xrel[a][0]*Rf[a][2];
    m[2] += xrel[a][0]*Rf[a][1] - xrel[a][1]*Rf[a][0];
  }

  // f_global,b = R f_d,b − (1/n) Σ R f_d,a − Γ_{:,b}ᵀ m
  const double invn = 1.0 / n;
  for (int b = 0; b < n; b++)
    for (int d = 0; d < 3; d++) {
      double gtm = Gamma[0][3*b+d]*m[0] + Gamma[1][3*b+d]*m[1]
                 + Gamma[2][3*b+d]*m[2];
      fGlobal(3*b + d) = Rf[b][d] - invn*fSum[d] - gtm;
    }
  return 0;
}

// --------------------------------------------------------------------------- //
// (seam 3 + K_geo) stiffness to global. v2.0:
//   K = R̄ k_d R̄ᵀ  +  K_geo ,
//   R̄ = blockdiag(R) (rotate every 3x3 block: K(a,b) = R k_d(a,b) Rᵀ),
//   K_geo[3b+i, c] = − [ spin(R f_d,b) Γ ]_{i,c}  (dominant load/spin term),
//   then symmetrise K. The polar-Hessian (∂Γ/∂u·m) and lever-arm terms are
//   deferred (see header / ADR §5). kGlobal MAY alias kCore.
// --------------------------------------------------------------------------- //
int
SolidTransformationCorot::globalizeStiff(const Matrix &kCore,
                                         const Vector &fCore,
                                         Matrix &kGlobal) const
{
  const int n    = numNodes;
  if (n == 0) {                       // queried before first update(): identity
    if (&kGlobal != &kCore) kGlobal = kCore;
    return 0;
  }
  const int ndof = 3 * n;

  // scratch (aliasing-safe). Stack-local so the seam is reentrant / thread-safe
  // (ndof ≤ 30 — 8-node hex = 24, 10-node tet = 30; allocation is cheap); the
  // ctor zero-initialises.
  Matrix K(ndof, ndof);

  // (1) rotated material stiffness  K(a,b) = R k_d(a,b) Rᵀ , block by block.
  double Rt[9]; transpose3(Rmat, Rt);
  for (int a = 0; a < n; a++)
    for (int b = 0; b < n; b++) {
      double blk[9], tmp[9], out[9];
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          blk[3*i+j] = kCore(3*a+i, 3*b+j);
      matmul3(Rmat, blk, tmp);          // R · blk
      matmul3(tmp, Rt, out);            // (R · blk) · Rᵀ
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          K(3*a+i, 3*b+j) = out[3*i+j];
    }

  // (2) geometric / load-stiffness. The single load term G1[3b+i, c] =
  //     − [spin(R f_d,b) Γ]_{i,c} is NOT symmetric by itself; the consistent
  //     dominant geometric stiffness is its symmetric pair G1 + G1ᵀ (the G2
  //     lever-arm term equals G1ᵀ at equilibrium). Assemble BOTH halves so the
  //     full magnitude survives — assembling only G1 and then averaging would
  //     halve it (review finding C1). TWO higher-order pieces remain deferred:
  //     the polar-Hessian ∂Γ/∂u·m AND the material-frame coupling
  //     R k_d (∂Rᵀ/∂u)x⁰; both are O(strain·force) and are what the FD-tangent
  //     gate tolerates (add both for quadratic Newton in v2.1).  // Ladruno (sweep #10)
  for (int b = 0; b < n; b++) {
    double fd[3] = { fCore(3*b), fCore(3*b+1), fCore(3*b+2) };
    double Rfb[3]; matvec3(Rmat, fd, Rfb);
    double Sp[9]; spin3(Rfb, Sp);                  // spin(R f_d,b)
    for (int i = 0; i < 3; i++)
      for (int c = 0; c < ndof; c++) {
        double v = Sp[3*i]*Gamma[0][c] + Sp[3*i+1]*Gamma[1][c]
                 + Sp[3*i+2]*Gamma[2][c];
        K(3*b + i, c) += -v;   // G1 into rows
        K(c, 3*b + i) += -v;   // G1ᵀ into columns ⇒ symmetric G1 + G1ᵀ
      }
  }

  // (3) symmetrise — now a rounding-level no-op (material R k_d Rᵀ and the
  //     G1+G1ᵀ geometric part are both symmetric by construction).
  for (int i = 0; i < ndof; i++)
    for (int j = i+1; j < ndof; j++) {
      double avg = 0.5 * (K(i, j) + K(j, i));
      K(i, j) = K(j, i) = avg;
    }

  kGlobal = K;
  return 0;
}

// --------------------------------------------------------------------------- //
void
SolidTransformationCorot::getCurrentFrame(Matrix &R) const
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      R(i, j) = Rmat[3*i + j];
}
