// Standalone cross-check of the CHANNEL-B spatial-tangent block of
// SRC/material/nD/LadrunoJ2Finite.cpp::getSpatialTangentTensor — the one code path
// the rest of the test suite does not reach numerically (the adversarial review's
// one medium finding). It mirrors, line for line, the channel-B math: the R-only
// perturbation (cauchyFromR), the central-FD dSdF, and the velocity-gradient push
// cmatB_ijkl = Σ_m dSdF[i][j][k][m]·F_lm. Channel B is PURELY constitutive (∂σ/∂α̃),
// with NO geometric/initial-stress term, so the assembled cmatB can be compared
// directly — with zero OpenSees commit semantics — to the numpy oracle's
// channelB_dsigma_dF mapped through the SAME Σ_m·F_lm transform.
//
//   g++ -O2 -std=c++17 -I ../SRC/material/nD ladrunoj2_finite_native_tangentB_check.cpp -o lj2fintB

#include <cstdio>
#include <cmath>
#include <LadrunoJ2Kernel.h>
#include <LogStrainKernel.h>

using namespace logstrain_kernel;
using ladruno_j2_kernel::Params;
static const double ZERO6[6] = {0,0,0,0,0,0};

static void polarRotation(const double Fd[9], double R[9]) {
  double FdT[9]; mat3_transpose(Fd, FdT);
  double C[9];   mat3_mul(FdT, Fd, C);
  double val[3], vec[9]; jacobi3(C, val, vec);
  double vmax = fmax(val[0], fmax(val[1], val[2]));
  double vfloor = 1.0e-30 * (vmax > 0.0 ? vmax : 1.0);
  double Uinv[9]; for (int i=0;i<9;i++) Uinv[i]=0.0;
  for (int a=0;a<3;a++){ double lam=(val[a]>vfloor)?val[a]:vfloor; double s=1.0/sqrt(lam);
    for (int i=0;i<3;i++) for(int j=0;j<3;j++) Uinv[3*i+j]+=s*vec[3*i+a]*vec[3*j+a]; }
  mat3_mul(Fd, Uinv, R);
}
static void rotateSym6(const double R[9], const double a6[6], double o[6]) {
  double A[9]; A[0]=a6[0];A[4]=a6[1];A[8]=a6[2];
  A[1]=A[3]=a6[3]; A[5]=A[7]=a6[4]; A[2]=A[6]=a6[5];
  double RA[9]; mat3_mul(R,A,RA); double RT[9]; mat3_transpose(R,RT);
  double RAR[9]; mat3_mul(RA,RT,RAR);
  o[0]=RAR[0];o[1]=RAR[4];o[2]=RAR[8];
  o[3]=0.5*(RAR[1]+RAR[3]); o[4]=0.5*(RAR[5]+RAR[7]); o[5]=0.5*(RAR[2]+RAR[6]);
}

struct Mat {
  Params p; int nB;
  double Fn[9], Be_n[9], ebarP_n, alpha_n[8][6];
  void init(const Params& pp){ p=pp; nB=pp.nBack;
    for(int i=0;i<9;i++){Fn[i]=0;Be_n[i]=0;} Fn[0]=Fn[4]=Fn[8]=1; Be_n[0]=Be_n[4]=Be_n[8]=1;
    ebarP_n=0; for(int k=0;k<8;k++) for(int i=0;i<6;i++) alpha_n[k][i]=0; }

  void step(const double F[9]) {  // setTrialF + commit (build the committed base)
    double Fninv[9]; mat3_inv(Fn,Fninv);
    double Fd[9]; mat3_mul(F,Fninv,Fd);
    double Betr[9]; trial_Be(F,Fn,Be_n,Betr);
    double eng[6]; hencky_voigt(Betr,eng);
    double st[6]={eng[0],eng[1],eng[2],0.5*eng[3],0.5*eng[4],0.5*eng[5]};
    double R[9]; polarRotation(Fd,R);
    double at[8][6]; for(int k=0;k<nB;k++) rotateSym6(R,alpha_n[k],at[k]);
    double s6[6],D[6][6],epsP[6],eb,al[8][6],dG;
    ladruno_j2_kernel::returnMap(p,st,ZERO6,ebarP_n,at,s6,D,epsP,eb,al,dG);
    double ee[6]={eng[0]-epsP[0],eng[1]-epsP[1],eng[2]-epsP[2],
                  eng[3]-2*epsP[3],eng[4]-2*epsP[4],eng[5]-2*epsP[5]};
    double Beu[9]; be_from_hencky_voigt(ee,Beu);
    for(int i=0;i<9;i++){Fn[i]=F[i];Be_n[i]=Beu[i];} ebarP_n=eb;
    for(int k=0;k<nB;k++) for(int i=0;i<6;i++) alpha_n[k][i]=al[k][i];
  }

  // Cauchy sigma (3x3) at trial F with the committed base, perturbing only R.
  void cauchyFromR(const double st[6], const double R[9], double invJ, double sig[3][3]) {
    double at[8][6]; for(int k=0;k<nB;k++) rotateSym6(R,alpha_n[k],at[k]);
    double s6[6],D[6][6],epsP[6],eb,al[8][6],dG;
    ladruno_j2_kernel::returnMap(p,st,ZERO6,ebarP_n,at,s6,D,epsP,eb,al,dG);
    sig[0][0]=s6[0]*invJ; sig[1][1]=s6[1]*invJ; sig[2][2]=s6[2]*invJ;
    sig[0][1]=sig[1][0]=s6[3]*invJ; sig[1][2]=sig[2][1]=s6[4]*invJ; sig[0][2]=sig[2][0]=s6[5]*invJ;
  }

  int trialPlastic(const double F[9]) {  // 1 if the eval step yields
    double Fninv[9]; mat3_inv(Fn,Fninv); double Fd[9]; mat3_mul(F,Fninv,Fd);
    double Betr[9]; trial_Be(F,Fn,Be_n,Betr); double eng[6]; hencky_voigt(Betr,eng);
    double st[6]={eng[0],eng[1],eng[2],0.5*eng[3],0.5*eng[4],0.5*eng[5]};
    double R[9]; polarRotation(Fd,R); double at[8][6]; for(int k=0;k<nB;k++) rotateSym6(R,alpha_n[k],at[k]);
    double s6[6],D[6][6],epsP[6],eb,al[8][6],dG;
    ladruno_j2_kernel::returnMap(p,st,ZERO6,ebarP_n,at,s6,D,epsP,eb,al,dG);
    return dG>0.0 ? 1 : 0;
  }

  // The channel-B block mirror: cB[i][j][k][l] = sum_m dSdF[i][j][k][m] * F[l][m].
  void channelB(const double F[9], double cB[3][3][3][3]) {
    double Fninv[9]; mat3_inv(Fn,Fninv);
    double Betr[9]; trial_Be(F,Fn,Be_n,Betr);
    double eng[6]; hencky_voigt(Betr,eng);
    double st[6]={eng[0],eng[1],eng[2],0.5*eng[3],0.5*eng[4],0.5*eng[5]};
    double J=mat3_det(F), invJ=1.0/J;
    const double h=1.0e-7;
    double dSdF[3][3][3][3];
    for(int k=0;k<3;k++) for(int m=0;m<3;m++){
      double Fp[9],Fm[9]; for(int q=0;q<9;q++){Fp[q]=F[q];Fm[q]=F[q];}
      Fp[3*k+m]+=h; Fm[3*k+m]-=h;
      double Fdp[9]; mat3_mul(Fp,Fninv,Fdp); double Rp[9]; polarRotation(Fdp,Rp);
      double Fdm[9]; mat3_mul(Fm,Fninv,Fdm); double Rm[9]; polarRotation(Fdm,Rm);
      double sp[3][3],sm[3][3];
      cauchyFromR(st,Rp,invJ,sp); cauchyFromR(st,Rm,invJ,sm);
      for(int i=0;i<3;i++) for(int j=0;j<3;j++) dSdF[i][j][k][m]=(sp[i][j]-sm[i][j])/(2.0*h);
    }
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++){
      double acc=0.0; for(int m=0;m<3;m++) acc+=dSdF[i][j][k][m]*F[3*l+m];
      cB[i][j][k][l]=acc;
    }
  }
};

static void rot(double ax,double ay,double az,double th,double R[9]){
  double n=sqrt(ax*ax+ay*ay+az*az); ax/=n;ay/=n;az/=n;
  double c=cos(th),s=sin(th),C=1-c;
  R[0]=c+ax*ax*C; R[1]=ax*ay*C-az*s; R[2]=ax*az*C+ay*s;
  R[3]=ay*ax*C+az*s; R[4]=c+ay*ay*C; R[5]=ay*az*C-ax*s;
  R[6]=az*ax*C-ay*s; R[7]=az*ay*C+ax*s; R[8]=c+az*az*C;
}
static void diagStretch(double a, double F[9]){
  for(int i=0;i<9;i++)F[i]=0; F[0]=exp(a); F[4]=exp(-a/2); F[8]=exp(-a/2);
}

int main(){
  Params p; p.K=1.5e3; p.G=7.0e2; p.sig0=10.0; p.Qinf=6.0; p.bIso=25.0; p.Hiso=40.0;
  p.nBack=3; double C[3]={600,350,120}, g[3]={120,60,8};
  for(int k=0;k<8;k++){p.C[k]=(k<3)?C[k]:0; p.gam[k]=(k<3)?g[k]:0;}
  Mat m; m.init(p);
  // committed plastic base: 3 symmetric stretch steps (alpha built, R=I)
  for(int s=1;s<=3;s++){ double F[9]; diagStretch(0.05*s, F); m.step(F); }
  // rotated + extra stretch eval ⇒ f_D carries rotation ⇒ channel B active
  double U[9]; diagStretch(0.20, U); double R[9]; rot(0.3,1.0,0.5,0.35,R);
  double Fe[9]; mat3_mul(R,U,Fe);
  printf("PLASTIC %d\n", m.trialPlastic(Fe));
  double cB[3][3][3][3]; m.channelB(Fe, cB);
  printf("FEVAL"); for(int q=0;q<9;q++) printf(" %.15e", Fe[q]); printf("\n");
  printf("CB");
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
    printf(" %.15e", cB[i][j][k][l]);
  printf("\n");
  return 0;
}
