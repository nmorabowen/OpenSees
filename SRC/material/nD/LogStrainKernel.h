/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// LogStrainKernel — the PURE numerical core of the logarithmic (Hencky)
// strain-space finite-strain adaptor. Plain doubles + <math.h>, NO OpenSees
// dependency, so it can be unit-tested standalone against the numpy oracle
// (tests/logstrain_reference.py) WITHOUT building OpenSees.
//
// Tensors: 3×3 symmetric stored row-major in double[9]; 4th-order [3][3][3][3];
// 6×6 stored row-major in double[36]; 6-vectors in double[6] with Voigt order
// {00,11,22,01,12,20} and ENGINEERING shear (γ_ij = 2 ε_ij).
//
// References: de Souza Neto, Perić & Owen (2008): Box 14.3 (MATISU), eq. 14.30
// (½ln Bᵉ), §14.5 (14.99–14.102) spatial tangent, Appendix A.5/A.6 (the
// isotropic-tensor-function derivative with the three eigenvalue-degeneracy
// branches — A.52/A.53).

#ifndef LogStrainKernel_h
#define LogStrainKernel_h

#include <math.h>

namespace logstrain_kernel {

// Voigt {00,11,22,01,12,20} representatives and (i,j)→Voigt map.
const int REP_I[6] = {0, 1, 2, 0, 1, 2};
const int REP_J[6] = {0, 1, 2, 1, 2, 0};

inline int VI(int i, int j) {
  if (i == j) return i;
  if ((i == 0 && j == 1) || (i == 1 && j == 0)) return 3;
  if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 4;
  return 5;
}

inline double halfLog(double x)      { return 0.5 * log(x); }
inline double halfLogPrime(double x) { return 0.5 / x; }

// ---- 3×3 helpers (row-major double[9]) ---------------------------------- //
inline void mat3_mul(const double A[9], const double B[9], double C[9]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += A[3*i+k] * B[3*k+j];
      C[3*i+j] = s;
    }
}
inline void mat3_transpose(const double A[9], double At[9]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) At[3*i+j] = A[3*j+i];
}
inline double mat3_det(const double A[9]) {
  return A[0]*(A[4]*A[8]-A[5]*A[7]) - A[1]*(A[3]*A[8]-A[5]*A[6])
       + A[2]*(A[3]*A[7]-A[4]*A[6]);
}
inline double mat3_inv(const double A[9], double Ai[9]) {
  double det = mat3_det(A);
  double id = (det != 0.0) ? 1.0/det : 0.0;
  Ai[0] =  (A[4]*A[8]-A[5]*A[7])*id;  Ai[1] = -(A[1]*A[8]-A[2]*A[7])*id;
  Ai[2] =  (A[1]*A[5]-A[2]*A[4])*id;  Ai[3] = -(A[3]*A[8]-A[5]*A[6])*id;
  Ai[4] =  (A[0]*A[8]-A[2]*A[6])*id;  Ai[5] = -(A[0]*A[5]-A[2]*A[3])*id;
  Ai[6] =  (A[3]*A[7]-A[4]*A[6])*id;  Ai[7] = -(A[0]*A[7]-A[1]*A[6])*id;
  Ai[8] =  (A[0]*A[4]-A[1]*A[3])*id;
  return det;
}

// ---- symmetric 3×3 eigensolver (cyclic Jacobi; robust at degeneracy) ---- //
// eval[a] eigenvalues; evec[3*i + a] = component i of eigenvector a.
inline void jacobi3(const double Ain[9], double eval[3], double evec[9]) {
  double a[3][3], v[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) { a[i][j] = Ain[3*i+j]; v[i][j] = (i==j)?1.0:0.0; }
  for (int sweep = 0; sweep < 50; sweep++) {
    double off = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
    if (off < 1.0e-300) break;
    const int pq[3][2] = {{0,1},{0,2},{1,2}};
    for (int r = 0; r < 3; r++) {
      int p = pq[r][0], q = pq[r][1];
      if (fabs(a[p][q]) < 1.0e-300) continue;
      double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
      double t = (theta >= 0.0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta*theta + 1.0));
      double c = 1.0 / sqrt(t*t + 1.0), s = t * c;
      double app = a[p][p], aqq = a[q][q], apq = a[p][q];
      a[p][p] = c*c*app - 2.0*s*c*apq + s*s*aqq;
      a[q][q] = s*s*app + 2.0*s*c*apq + c*c*aqq;
      a[p][q] = a[q][p] = 0.0;
      for (int k = 0; k < 3; k++) {
        if (k != p && k != q) {
          double akp = a[k][p], akq = a[k][q];
          a[k][p] = a[p][k] = c*akp - s*akq;
          a[k][q] = a[q][k] = s*akp + c*akq;
        }
        double vkp = v[k][p], vkq = v[k][q];
        v[k][p] = c*vkp - s*vkq;
        v[k][q] = s*vkp + c*vkq;
      }
    }
  }
  for (int A = 0; A < 3; A++) {
    eval[A] = a[A][A];
    for (int i = 0; i < 3; i++) evec[3*i + A] = v[i][A];
  }
}

// ---- isotropic tensor function Y = Σ y(xᵢ) Eᵢ  (A.51) ------------------- //
inline void isoFunction(const double X[9], double (*y)(double), double Y[9]) {
  double val[3], vec[9];
  jacobi3(X, val, vec);
  for (int i = 0; i < 9; i++) Y[i] = 0.0;
  for (int a = 0; a < 3; a++) {
    double ya = y(val[a]);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        Y[3*i+j] += ya * vec[3*i+a] * vec[3*j+a];
  }
}

// ---- derivative dY/dX (A.52/A.53) with the three multiplicity branches -- //
inline void isoFunctionDeriv(const double X[9], double (*y)(double),
                             double (*yp)(double), double D[3][3][3][3]) {
  double val[3], vec[9];
  jacobi3(X, val, vec);
  double E[3][3][3];
  for (int a = 0; a < 3; a++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) E[a][i][j] = vec[3*i+a]*vec[3*j+a];

  double IS[3][3][3][3], dX2[3][3][3][3];
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) {
    double dik=(i==k), djl=(j==l), dil=(i==l), djk=(j==k);
    IS[i][j][k][l]  = 0.5*(dik*djl + dil*djk);
    dX2[i][j][k][l] = 0.5*(dik*X[3*l+j] + dil*X[3*k+j] + djl*X[3*i+k] + djk*X[3*i+l]);
  }
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
    D[i][j][k][l] = 0.0;

  double scale = fabs(val[0]);
  if (fabs(val[1]) > scale) scale = fabs(val[1]);
  if (fabs(val[2]) > scale) scale = fabs(val[2]);
  if (scale < 1.0) scale = 1.0;
  const double rtol = 1.0e-9;
  bool e01 = fabs(val[0]-val[1]) <= rtol*scale;
  bool e12 = fabs(val[1]-val[2]) <= rtol*scale;
  bool e20 = fabs(val[2]-val[0]) <= rtol*scale;

  if (e01 && e12 && e20) {                     // triple
    double yp0 = yp(val[0]);
    for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
      D[i][j][k][l] = yp0 * IS[i][j][k][l];
    return;
  }
  if (e01 || e12 || e20) {                     // two equal (A.52 second / A.53)
    int a; if (e01) a = 2; else if (e12) a = 0; else a = 1;
    int c = (a+1)%3;
    double xa = val[a], xc = val[c];
    double ya = y(xa), yc = y(xc), ypa = yp(xa), ypc = yp(xc);
    double d = xa - xc;
    double s1 = (ya - yc)/(d*d) - ypc/d;
    double s2 = 2.0*xc*(ya - yc)/(d*d) - (xa + xc)/d*ypc;
    double s3 = 2.0*(ya - yc)/(d*d*d) - (ypa + ypc)/(d*d);
    double s4 = xc*s3, s5 = xc*s3, s6 = xc*xc*s3;
    for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) {
      double dij=(i==j), dkl=(k==l);
      D[i][j][k][l] = s1*dX2[i][j][k][l] - s2*IS[i][j][k][l]
                    - s3*X[3*i+j]*X[3*k+l] + s4*X[3*i+j]*dkl
                    + s5*dij*X[3*k+l]      - s6*dij*dkl;
    }
    return;
  }
  for (int a = 0; a < 3; a++) {                // all distinct (A.52 first)
    int b = (a+1)%3, c = (a+2)%3;
    double xa=val[a], xb=val[b], xc=val[c];
    double coef = y(xa) / ((xa-xb)*(xa-xc));
    double f1 = (xa-xb) + (xa-xc), f2 = (xb-xc), ypa = yp(xa);
    for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) {
      double EaEa = E[a][i][j]*E[a][k][l];
      double EbEb = E[b][i][j]*E[b][k][l];
      double EcEc = E[c][i][j]*E[c][k][l];
      D[i][j][k][l] += coef*( dX2[i][j][k][l] - (xb+xc)*IS[i][j][k][l]
                              - f1*EaEa - f2*(EbEb - EcEc) ) + ypa*EaEa;
    }
  }
}

// ===================================================================== //
//  High-level steps (used by both LogStrainNDMaterial and the unit test) //
// ===================================================================== //

// (i)+(ii) Bᵉᵗʳ = F_Δ Bᵉ_n F_Δᵀ,  F_Δ = F F_n⁻¹  (symmetrised).
inline void trial_Be(const double F[9], const double Fn[9],
                     const double Be_n[9], double Betr[9]) {
  double Fninv[9]; mat3_inv(Fn, Fninv);
  double Fd[9];    mat3_mul(F, Fninv, Fd);
  double tmp[9];   mat3_mul(Fd, Be_n, tmp);
  double FdT[9];   mat3_transpose(Fd, FdT);
  mat3_mul(tmp, FdT, Betr);
  Betr[1] = Betr[3] = 0.5*(Betr[1]+Betr[3]);
  Betr[2] = Betr[6] = 0.5*(Betr[2]+Betr[6]);
  Betr[5] = Betr[7] = 0.5*(Betr[5]+Betr[7]);
}

// εᵉᵗʳ = ½ ln Bᵉᵗʳ → engineering-shear Voigt 6-vector.
inline void hencky_voigt(const double Betr[9], double eps6[6]) {
  double e[9]; isoFunction(Betr, halfLog, e);
  eps6[0]=e[0]; eps6[1]=e[4]; eps6[2]=e[8];
  eps6[3]=2.0*e[1]; eps6[4]=2.0*e[5]; eps6[5]=2.0*e[2];
}

// σ = J⁻¹ τ  and material spatial tangent c = (1/2J)[D : L : B]  (14.99,
// geometric term excluded — that is the element's K_geo). D6/tau6 come from
// the inner small-strain law (engineering-Voigt).
inline void assemble_material(const double Betr[9], const double D6[36],
                              const double tau6[6], double J,
                              double sigma6[6], double c6[36]) {
  double invJ = (J != 0.0) ? 1.0/J : 0.0;
  for (int a = 0; a < 6; a++) sigma6[a] = tau6[a] * invJ;

  double L[3][3][3][3];
  isoFunctionDeriv(Betr, halfLog, halfLogPrime, L);   // ∂(½ln)/∂B
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
    L[i][j][k][l] *= 2.0;                              // ∂lnB/∂B

  double Bt[3][3][3][3], Dt[3][3][3][3];
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) {
    double dik=(i==k), djk=(j==k);
    Bt[i][j][k][l] = dik*Betr[3*j+l] + djk*Betr[3*i+l];
    Dt[i][j][k][l] = D6[6*VI(i,j) + VI(k,l)];
  }
  double DL[3][3][3][3];
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int r=0;r<3;r++) for(int s=0;s<3;s++) {
    double acc = 0.0;
    for (int p=0;p<3;p++) for(int q=0;q<3;q++) acc += Dt[i][j][p][q]*L[p][q][r][s];
    DL[i][j][r][s] = acc;
  }
  for (int I = 0; I < 6; I++) {
    int i = REP_I[I], j = REP_J[I];
    for (int Jn = 0; Jn < 6; Jn++) {
      int k = REP_I[Jn], l = REP_J[Jn];
      double acc = 0.0;
      for (int r=0;r<3;r++) for(int s=0;s<3;s++) acc += DL[i][j][r][s]*Bt[r][s][k][l];
      c6[6*I + Jn] = invJ * 0.5 * acc;
    }
  }
}

// FULL 4th-order spatial constitutive modulus c_ijkl = (1/2J)[D:L:B] for ALL 81
// components (NO Voigt collapse). The (k,l) pair of c is NOT minor-symmetric (B,
// eq.14.102, is not k<->l symmetric), so the 6×6 form above is LOSSY in (k,l); a
// finite-strain ELEMENT that needs the consistent tangent must use this full
// tensor (it then forms a_ijkl = c_ijkl − σ_il δ_jk and contracts with the full
// nodal gradients). geometric term still excluded — it is the element's K_geo.
inline void spatial_tangent_full(const double Betr[9], const double D6[36],
                                 double J, double c[3][3][3][3]) {
  double invJ = (J != 0.0) ? 1.0/J : 0.0;

  double L[3][3][3][3];
  isoFunctionDeriv(Betr, halfLog, halfLogPrime, L);   // ∂(½ln)/∂B
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
    L[i][j][k][l] *= 2.0;                              // ∂lnB/∂B

  double Bt[3][3][3][3], Dt[3][3][3][3];
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) {
    double dik=(i==k), djk=(j==k);
    Bt[i][j][k][l] = dik*Betr[3*j+l] + djk*Betr[3*i+l];
    Dt[i][j][k][l] = D6[6*VI(i,j) + VI(k,l)];
  }
  double DL[3][3][3][3];
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int r=0;r<3;r++) for(int s=0;s<3;s++) {
    double acc = 0.0;
    for (int p=0;p<3;p++) for(int q=0;q<3;q++) acc += Dt[i][j][p][q]*L[p][q][r][s];
    DL[i][j][r][s] = acc;
  }
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) {
    double acc = 0.0;
    for (int r=0;r<3;r++) for(int s=0;s<3;s++) acc += DL[i][j][r][s]*Bt[r][s][k][l];
    c[i][j][k][l] = invJ * 0.5 * acc;
  }
}

// Reconstruct the elastic left Cauchy–Green tensor from an elastic Hencky strain
// given in engineering-shear Voigt:  bᵉ = exp[2 εᵉ]  (eq. 14.93). Used by the
// plastic-inner state protocol to update the committed bᵉ from the recovered
// εᵉ_{n+1} = Cᵉ:τ.
inline double _expd(double x) { return exp(x); }
inline void be_from_hencky_voigt(const double e6[6], double Be[9]) {
  double two_eps[9];
  two_eps[0] = 2.0*e6[0];      two_eps[4] = 2.0*e6[1];      two_eps[8] = 2.0*e6[2];
  two_eps[1] = two_eps[3] = e6[3];   // 2·(½ γ) = γ = e6[3]   (engineering→tensor·2)
  two_eps[5] = two_eps[7] = e6[4];
  two_eps[2] = two_eps[6] = e6[5];
  isoFunction(two_eps, _expd, Be);
}

} // namespace logstrain_kernel

#endif
