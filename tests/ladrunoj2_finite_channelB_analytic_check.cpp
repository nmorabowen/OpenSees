// Standalone cross-check of the ANALYTIC channel-B block of
// SRC/material/nD/LadrunoJ2Finite.cpp::getSpatialTangentTensor (ADR-16 P3), using ONLY
// the OpenSees-free kernels. It mirrors, line for line, the material's analytic chain:
// the polar-rotation derivative polarRotationDeriv (Sylvester axial solve), the
// return-map backstress sensitivity ∂τ/∂α̃, and the velocity-gradient push
// cmatB_ijkl = Σ_m dSdF[i][j][k][m]·F_lm. Channel B is purely constitutive (no
// geometric term), so the assembled cmatB compares directly — zero commit semantics —
// to the numpy oracle's channelB_dsigma_dF mapped through the SAME Σ_m·F_lm transform.
// The python driver tests/test_ladrunoJ2_finite_channelB_cpp.py compiles this and
// compares to tests/ladrunoj2_finite_native_reference.py.
//
//   g++ -O2 -std=c++17 -I ../SRC/material/nD ladrunoj2_finite_channelB_analytic_check.cpp -o lj2cba

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
static inline double ddot3(const double A[9], const double B[9])
{ double s=0; for(int i=0;i<9;i++) s+=A[i]*B[i]; return s; }

// dR for f = R U (polar), perturbation df: (tr U·I − U) ω = axial(A−Aᵀ), A = Rᵀ df.
static void polarRotationDeriv(const double f[9], const double df[9],
                               const double R[9], double dR[9]) {
  double RT[9]; mat3_transpose(R, RT);
  double U[9];  mat3_mul(RT, f, U);
  for(int i=0;i<3;i++) for(int j=i+1;j<3;j++){ double a=0.5*(U[3*i+j]+U[3*j+i]); U[3*i+j]=U[3*j+i]=a; }
  double A[9]; mat3_mul(RT, df, A);
  double w[3] = { A[7]-A[5], A[2]-A[6], A[3]-A[1] };
  double trU = U[0]+U[4]+U[8];
  double K[9]; for(int i=0;i<9;i++) K[i]=-U[i]; K[0]+=trU; K[4]+=trU; K[8]+=trU;
  double Kinv[9]; mat3_inv(K, Kinv);
  double om[3]; for(int i=0;i<3;i++) om[i]=Kinv[3*i]*w[0]+Kinv[3*i+1]*w[1]+Kinv[3*i+2]*w[2];
  double Om[9] = {0,-om[2],om[1], om[2],0,-om[0], -om[1],om[0],0};
  mat3_mul(R, Om, dR);
}

struct Mat {
  Params p; int nB;
  double Fn[9], Be_n[9], ebarP_n, alpha_n[8][6];
  void init(const Params& pp){ p=pp; nB=pp.nBack;
    for(int i=0;i<9;i++){Fn[i]=0;Be_n[i]=0;} Fn[0]=Fn[4]=Fn[8]=1; Be_n[0]=Be_n[4]=Be_n[8]=1;
    ebarP_n=0; for(int k=0;k<8;k++) for(int i=0;i<6;i++) alpha_n[k][i]=0; }

  void step(const double F[9]) {
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

  int trialPlastic(const double F[9]) {
    double Fninv[9]; mat3_inv(Fn,Fninv); double Fd[9]; mat3_mul(F,Fninv,Fd);
    double Betr[9]; trial_Be(F,Fn,Be_n,Betr); double eng[6]; hencky_voigt(Betr,eng);
    double st[6]={eng[0],eng[1],eng[2],0.5*eng[3],0.5*eng[4],0.5*eng[5]};
    double R[9]; polarRotation(Fd,R); double at[8][6]; for(int k=0;k<nB;k++) rotateSym6(R,alpha_n[k],at[k]);
    double s6[6],D[6][6],epsP[6],eb,al[8][6],dG;
    ladruno_j2_kernel::returnMap(p,st,ZERO6,ebarP_n,at,s6,D,epsP,eb,al,dG);
    return dG>0.0 ? 1 : 0;
  }

  // ANALYTIC channel B — mirror of LadrunoJ2Finite::channelBAnalytic.
  void channelB(const double F[9], double cB[3][3][3][3]) {
    const double G=p.G, root23=sqrt(2.0/3.0);
    double Fninv[9]; mat3_inv(Fn,Fninv);
    double Fd[9]; mat3_mul(F,Fninv,Fd);
    double R[9]; polarRotation(Fd,R);
    double Betr[9]; trial_Be(F,Fn,Be_n,Betr);
    double eng[6]; hencky_voigt(Betr,eng);
    double st[6]={eng[0],eng[1],eng[2],0.5*eng[3],0.5*eng[4],0.5*eng[5]};
    double J=mat3_det(F), invJ=1.0/J;
    // converged dG at the probe (recompute)
    double at[8][6]; for(int k=0;k<nB;k++) rotateSym6(R,alpha_n[k],at[k]);
    double s6[6],D[6][6],epsP[6],eb,al[8][6],dG;
    ladruno_j2_kernel::returnMap(p,st,ZERO6,ebarP_n,at,s6,D,epsP,eb,al,dG);

    // aRef (committed α as 3×3) + aux (3×3)
    double aRef[8][9];
    for(int q=0;q<nB;q++){ aRef[q][0]=alpha_n[q][0];aRef[q][4]=alpha_n[q][1];aRef[q][8]=alpha_n[q][2];
      aRef[q][1]=aRef[q][3]=alpha_n[q][3]; aRef[q][5]=aRef[q][7]=alpha_n[q][4]; aRef[q][2]=aRef[q][6]=alpha_n[q][5]; }
    double aChi[8][9];
    for(int q=0;q<nB;q++){ double RA[9]; mat3_mul(R,aRef[q],RA); double RT[9]; mat3_transpose(R,RT); mat3_mul(RA,RT,aChi[q]); }
    double eps[9]={st[0],st[3],st[5], st[3],st[1],st[4], st[5],st[4],st[2]};
    double trE=eps[0]+eps[4]+eps[8];
    double sTr[9]; for(int i=0;i<9;i++) sTr[i]=2.0*G*eps[i]; sTr[0]-=2.0*G*trE/3.0; sTr[4]-=2.0*G*trE/3.0; sTr[8]-=2.0*G*trE/3.0;
    double Dk[8]; double M[9]; for(int i=0;i<9;i++)M[i]=sTr[i]; double Mp[9]={0,0,0,0,0,0,0,0,0};
    for(int q=0;q<nB;q++){ Dk[q]=1.0+root23*p.gam[q]*dG;
      for(int i=0;i<9;i++){ M[i]-=aChi[q][i]/Dk[q]; Mp[i]+=aChi[q][i]*root23*p.gam[q]/(Dk[q]*Dk[q]); } }
    double normM=sqrt(ddot3(M,M)); double n[9]; for(int i=0;i<9;i++)n[i]=M[i]/normM;
    double dtheta=2.0*G; for(int q=0;q<nB;q++) dtheta+=(2.0/3.0)*p.C[q]/(Dk[q]*Dk[q]);
    double pbar=ebarP_n+root23*dG;
    double h=dtheta+(2.0/3.0)*ladruno_j2_kernel::yieldSlope(p,pbar)-ddot3(n,Mp);

    double dSdF[3][3][3][3];
    for(int k=0;k<3;k++) for(int mIdx=0;mIdx<3;mIdx++){
      double dFd[9]={0,0,0,0,0,0,0,0,0};
      dFd[3*k+0]=Fninv[3*mIdx+0]; dFd[3*k+1]=Fninv[3*mIdx+1]; dFd[3*k+2]=Fninv[3*mIdx+2];
      double dR[9]; polarRotationDeriv(Fd,dFd,R,dR);
      double dRT[9]; mat3_transpose(dR,dRT); double RT[9]; mat3_transpose(R,RT);
      double ds[9]={0,0,0,0,0,0,0,0,0};
      for(int q=0;q<nB;q++){
        double t1[9],dA[9],t2[9];
        mat3_mul(dR,aRef[q],t1); mat3_mul(t1,RT,dA);
        mat3_mul(R,aRef[q],t1);  mat3_mul(t1,dRT,t2);
        for(int i=0;i<9;i++) dA[i]+=t2[i];
        double dGp=-ddot3(n,dA)/(h*Dk[q]);
        double dM[9]; for(int i=0;i<9;i++) dM[i]=-dA[i]/Dk[q]+Mp[i]*dGp;
        double ndM=ddot3(n,dM);
        for(int i=0;i<9;i++){ double dn=(dM[i]-n[i]*ndM)/normM; ds[i]+=-2.0*G*(dGp*n[i]+dG*dn); }
      }
      for(int i=0;i<3;i++) for(int j=0;j<3;j++) dSdF[i][j][k][mIdx]=ds[3*i+j]*invJ;
    }
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++){
      double acc=0.0; for(int mm=0;mm<3;mm++) acc+=dSdF[i][j][k][mm]*F[3*l+mm];
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
static void diagStretch(double a, double F[9]){ for(int i=0;i<9;i++)F[i]=0; F[0]=exp(a); F[4]=exp(-a/2); F[8]=exp(-a/2); }

int main(){
  Params p; p.K=1.5e3; p.G=7.0e2; p.sig0=10.0; p.Qinf=6.0; p.bIso=25.0; p.Hiso=40.0;
  p.nBack=3; double C[3]={600,350,120}, g[3]={120,60,8};
  for(int k=0;k<8;k++){p.C[k]=(k<3)?C[k]:0; p.gam[k]=(k<3)?g[k]:0;}
  Mat m; m.init(p);
  for(int s=1;s<=3;s++){ double F[9]; diagStretch(0.05*s, F); m.step(F); }
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
