// Standalone cross-check of the finite-strain-native combined-hardening J2
// composition (polar co-rotation of the backstress + the shared return map),
// using ONLY the OpenSees-free kernels — NO OpenSees build required. It mirrors,
// line for line, the math of SRC/material/nD/LadrunoJ2Finite.cpp::setTrialF (the
// engineering<->tensor Voigt bridging is the part the numpy oracle cannot catch,
// since the oracle works in 3x3 tensors throughout). The python driver
// tests/test_ladrunoJ2_finite_native_cpp.py compiles this with g++, runs the SAME
// path through tests/ladrunoj2_finite_native_reference.py, and compares.
//
//   g++ -O2 -std=c++17 -I ../SRC/material/nD ladrunoj2_finite_native_check.cpp -o lj2fin

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
  double Uinv[9]; for (int i=0;i<9;i++) Uinv[i]=0.0;
  for (int a=0;a<3;a++){ double s=(val[a]>0.0)?1.0/sqrt(val[a]):0.0;
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
  // setTrialF + commit; writes Cauchy sig6, ebarP, total backstress tot6
  void step(const double F[9], double sig6[6], double& ebarP, double tot6[6]) {
    double Fninv[9]; mat3_inv(Fn,Fninv);
    double Fd[9]; mat3_mul(F,Fninv,Fd);
    double Betr[9]; trial_Be(F,Fn,Be_n,Betr);
    double eng[6]; hencky_voigt(Betr,eng);
    double st[6]={eng[0],eng[1],eng[2],0.5*eng[3],0.5*eng[4],0.5*eng[5]};
    double R[9]; polarRotation(Fd,R);
    double at[8][6]; for(int k=0;k<nB;k++) rotateSym6(R,alpha_n[k],at[k]);
    double s6[6],D[6][6],epsP[6],eb,al[8][6],dG;
    ladruno_j2_kernel::returnMap(p,st,ZERO6,ebarP_n,at,s6,D,epsP,eb,al,dG);
    double J=mat3_det(F), invJ=1.0/J;
    for(int i=0;i<6;i++) sig6[i]=s6[i]*invJ;
    double ee[6]={eng[0]-epsP[0],eng[1]-epsP[1],eng[2]-epsP[2],
                  eng[3]-2*epsP[3],eng[4]-2*epsP[4],eng[5]-2*epsP[5]};
    double Beu[9]; be_from_hencky_voigt(ee,Beu);
    // commit
    for(int i=0;i<9;i++){Fn[i]=F[i];Be_n[i]=Beu[i];} ebarP_n=eb;
    for(int k=0;k<nB;k++) for(int i=0;i<6;i++) alpha_n[k][i]=al[k][i];
    ebarP=ebarP_n;
    for(int i=0;i<6;i++){ tot6[i]=0; for(int k=0;k<nB;k++) tot6[i]+=alpha_n[k][i]; }
  }
};

static void rot(double ax,double ay,double az,double th,double R[9]){
  double n=sqrt(ax*ax+ay*ay+az*az); ax/=n;ay/=n;az/=n;
  double c=cos(th),s=sin(th),C=1-c;
  R[0]=c+ax*ax*C; R[1]=ax*ay*C-az*s; R[2]=ax*az*C+ay*s;
  R[3]=ay*ax*C+az*s; R[4]=c+ay*ay*C; R[5]=ay*az*C-ax*s;
  R[6]=az*ax*C-ay*s; R[7]=az*ay*C+ax*s; R[8]=c+az*az*C;
}

int main(){
  Params p; p.K=1.5e3; p.G=7.0e2; p.sig0=10.0; p.Qinf=6.0; p.bIso=25.0; p.Hiso=40.0;
  p.nBack=3; double C[3]={600,350,120}, g[3]={120,60,8};
  for(int k=0;k<8;k++){p.C[k]=(k<3)?C[k]:0; p.gam[k]=(k<3)?g[k]:0;}
  Mat m; m.init(p);
  // rotating + stretching plastic path (mirror of the oracle test)
  for(int k=1;k<=9;k++){
    double s=1.0+0.03*k, gg=0.015*k;
    double U[9]={s,gg,0, 0,1.0/sqrt(s),0, 0,0,1.0/sqrt(s)};
    double Rr[9]; rot(0,0,1, 0.07*k, Rr);
    double F[9]; mat3_mul(Rr,U,F);
    double sig[6],eb,tot[6]; m.step(F,sig,eb,tot);
    printf("%d", k);
    for(int i=0;i<6;i++) printf(" %.15e", sig[i]);
    printf(" %.15e", eb);
    for(int i=0;i<6;i++) printf(" %.15e", tot[i]);
    printf("\n");
  }
  return 0;
}
