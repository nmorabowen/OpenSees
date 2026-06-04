// Standalone validator for the Lemaitre damage layer (returnMapDamaged) in
// LadrunoJ2Kernel.h. Pure C++ + <math.h>, NO OpenSees build required.
//
//   * V0  damage-OFF == returnMap byte-identical (the regression tripwire)
//   * grows D monotonically on a uniaxial plastic path past the threshold
//   * FD-checks the CONSISTENT tangent dsigma/deps (analytic vs central FD) at
//     several damaged states, for uniaxial and 3D-mixed (shear-inclusive) loading
//   * triaxiality: a hydrostatically-confined state damages faster than uniaxial
//     at the same accumulated plastic strain (R_v signature)
//
// Build: g++ -O2 -std=c++17 -I SRC/material/nD tests/ladruno_damage_check.cpp -o ldc
// Exit code 0 = all checks pass.

#include <LadrunoJ2Kernel.h>
#include <cstdio>
#include <cmath>

using namespace ladruno_j2_kernel;

struct State { double epsP[6]; double ebarP; double alpha[MAXBACK][6]; double D; };

static void zero(State& st) {
  for (int i=0;i<6;i++) st.epsP[i]=0;
  st.ebarP=0; st.D=0;
  for (int k=0;k<MAXBACK;k++) for(int i=0;i<6;i++) st.alpha[k][i]=0;
}

// one damaged trial at total tensor-strain eps from committed state st (no commit)
static int trial(const Params& p, const double eps[6], const State& st,
                 double stress[6], double Dtan[6][6], double& Dout) {
  double epsP[6], ebarP, alpha[MAXBACK][6], dG;
  return returnMapDamaged(p, eps, st.epsP, st.ebarP, st.alpha, st.D,
                          stress, Dtan, epsP, ebarP, alpha, dG, Dout);
}

static void commit(const Params& p, const double eps[6], State& st) {
  double stress[6], Dtan[6][6], epsP[6], ebarP, alpha[MAXBACK][6], dG, Dout;
  returnMapDamaged(p, eps, st.epsP, st.ebarP, st.alpha, st.D,
                   stress, Dtan, epsP, ebarP, alpha, dG, Dout);
  for (int i=0;i<6;i++) st.epsP[i]=epsP[i];
  st.ebarP=ebarP; st.D=Dout;
  for (int k=0;k<p.nBack;k++) for(int i=0;i<6;i++) st.alpha[k][i]=alpha[k][i];
}

static int fails = 0;
static void check(bool ok, const char* msg, double got=0, double ref=0) {
  if (!ok) { printf("  FAIL: %s (got %.6g ref %.6g)\n", msg, got, ref); fails++; }
  else     printf("  ok:   %s\n", msg);
}

// central-difference tangent at committed state st, trial strain eps0
static void fdTangent(const Params& p, const double eps0[6], const State& st,
                      double Dfd[6][6]) {
  const double w[6] = {1,1,1,0.5,0.5,0.5};   // eng->tensor strain factor
  for (int j=0;j<6;j++) {
    double h = 1.0e-7;
    double ep[6], em[6];
    for (int i=0;i<6;i++) { ep[i]=eps0[i]; em[i]=eps0[i]; }
    ep[j] += h*w[j];   // +h on ENGINEERING strain comp j => +h*w on tensor comp
    em[j] -= h*w[j];
    double sp[6], sm[6], Dt[6][6], Dd;
    trial(p, ep, st, sp, Dt, Dd);
    trial(p, em, st, sm, Dt, Dd);
    for (int I=0;I<6;I++) Dfd[I][j] = (sp[I]-sm[I])/(2.0*h);
  }
}

static double maxTangErr(const Params& p, const double eps[6], const State& st) {
  double Da[6][6], Dfd[6][6], stress[6], Dout;
  trial(p, eps, st, stress, Da, Dout);
  fdTangent(p, eps, st, Dfd);
  double scale=0, err=0;
  for (int I=0;I<6;I++) for(int J=0;J<6;J++) {
    scale = fmax(scale, fabs(Da[I][J]));
    err   = fmax(err, fabs(Da[I][J]-Dfd[I][J]));
  }
  return err/fmax(scale,1.0);
}

int main()
{
  const double K=1500.0, G=700.0;

  // ---- V0: damage OFF must equal returnMap byte-identical ----
  printf("V0 damage-off bit-identical:\n");
  {
    Params p{K,G,8.0,3.0,20.0,5.0,3,{500,300,150},{100,50,10}};  // dmg off by default
    double eps[6]={0.5,-0.2,-0.3,0.4,0.1,0.2};
    for (int i=3;i<6;i++) eps[i]*=0.5; // tensor comps
    State st; zero(st);
    double s1[6],D1[6][6],eP[6],eb,al[MAXBACK][6],dg;
    int st1 = returnMap(p, eps, st.epsP, st.ebarP, st.alpha, s1,D1,eP,eb,al,dg);
    double s2[6],D2[6][6],Dout; int st2 = trial(p, eps, st, s2, D2, Dout);
    bool same=true;
    for (int i=0;i<6;i++) if (s1[i]!=s2[i]) same=false;
    for (int I=0;I<6;I++) for(int J=0;J<6;J++) if (D1[I][J]!=D2[I][J]) same=false;
    check(same && st1==st2 && Dout==0.0, "stress+tangent identical, D=0", Dout, 0.0);
  }

  // ---- damaging uniaxial path: D grows monotonically past threshold ----
  printf("damage growth (uniaxial, pD=0.01):\n");
  {
    Params p{K,G,5.0,4.0,30.0,10.0,1,{400},{80}};
    p.dmg.on=true; p.dmg.r=2.0; p.dmg.s=1.0; p.dmg.pD=0.01; p.dmg.Dc=0.99;
    State st; zero(st);
    double mult[8]={0.005,0.01,0.02,0.03,0.05,0.07,0.10,0.14};
    double dir[6]={1,0,0,0,0,0};
    double lastD=-1; bool mono=true, zeroPre=true, grewPost=false;
    for (int s=0;s<8;s++){
      double eps[6]; for(int i=0;i<6;i++) eps[i]=mult[s]*dir[i];
      commit(p, eps, st);
      if (st.ebarP < 0.01 && st.D>1e-14) zeroPre=false;     // no damage below threshold
      if (st.D + 1e-15 < lastD) mono=false;                  // never decreases
      if (st.D > 1e-6) grewPost=true;
      lastD=st.D;
    }
    check(zeroPre, "no damage below pD");
    check(mono,    "D monotonic non-decreasing");
    check(grewPost && st.D<0.99, "D grew into (0,Dc)", st.D);
  }

  // ---- consistent tangent FD: uniaxial + 3D-mixed, several damaged states ----
  printf("consistent-tangent FD (rel err < 1e-5):\n");
  {
    Params p{K,G,5.0,3.0,25.0,8.0,2,{500,250},{60,12}};
    p.dmg.on=true; p.dmg.r=3.0; p.dmg.s=1.0; p.dmg.pD=0.0; p.dmg.Dc=0.99;

    // uniaxial: build damage over a few committed steps, FD-check each trial
    State su; zero(su);
    double umul[5]={0.01,0.02,0.03,0.05,0.07};
    double worstU=0;
    for (int s=0;s<5;s++){
      double eps[6]={umul[s],0,0,0,0,0};
      double e=maxTangErr(p,eps,su); worstU=fmax(worstU,e);
      commit(p,eps,su);
    }
    check(worstU<1e-5, "uniaxial damaged tangent", worstU, 1e-5);

    // 3D mixed incl. shear (exercises the deviatoric + shear-convention terms)
    State sm; zero(sm);
    double base[6]={0.5,-0.2,-0.3,0.4,0.1,0.2};
    double mmul[5]={0.01,0.02,0.03,0.045,0.06};
    double worstM=0;
    for (int s=0;s<5;s++){
      double eps[6]; for(int i=0;i<6;i++){ eps[i]=mmul[s]*base[i]; if(i>=3) eps[i]*=0.5; }
      double e=maxTangErr(p,eps,sm); worstM=fmax(worstM,e);
      commit(p,eps,sm);
    }
    check(worstM<1e-5, "3D-mixed damaged tangent", worstM, 1e-5);

    // s != 1 exponent (nonlinear damage rate) tangent
    Params p2 = p; p2.dmg.s=2.0; p2.dmg.r=0.5;
    State s2; zero(s2);
    double worst2=0;
    for (int s=0;s<5;s++){
      double eps[6]={umul[s],0,0,0,0,0};
      double e=maxTangErr(p2,eps,s2); worst2=fmax(worst2,e);
      commit(p2,eps,s2);
    }
    check(worst2<1e-5, "uniaxial tangent, s=2 exponent", worst2, 1e-5);
  }

  // ---- lch regularization: lchScale linearly scales the damage increment ----
  printf("lch regularization (lchScale scales dD, tangent still consistent):\n");
  {
    Params p1{K,G,5.0,3.0,25.0,8.0,1,{400},{80}};
    p1.dmg.on=true; p1.dmg.r=3.0; p1.dmg.s=1.0; p1.dmg.pD=0.0; p1.dmg.Dc=0.99;
    Params p2 = p1; p2.dmg.lchScale = 2.0;        // double the band scale
    State a; zero(a); State b; zero(b);
    double eps[6]={0.03,0,0,0,0,0};
    // single step from virgin state: dD(lchScale=2) == 2*dD(lchScale=1)
    double s1[6],D1[6][6],eP[6],eb,al[MAXBACK][6],dg,Da,Db;
    returnMapDamaged(p1, eps, a.epsP, a.ebarP, a.alpha, a.D, s1,D1,eP,eb,al,dg,Da);
    returnMapDamaged(p2, eps, b.epsP, b.ebarP, b.alpha, b.D, s1,D1,eP,eb,al,dg,Db);
    check(fabs(Db - 2.0*Da) <= 1e-12*(Da+1.0), "lchScale=2 doubles dD", Db, 2.0*Da);
    // tangent stays FD-consistent under regularization
    double e2 = maxTangErr(p2, eps, b);
    check(e2 < 1e-5, "regularized tangent FD-consistent", e2, 1e-5);
  }

  // ---- triaxiality (R_v) signature: confined state damages faster ----
  printf("triaxiality R_v signature:\n");
  {
    Params p{K,G,5.0,0.0,0.0,20.0,0,{0},{0}};   // pure linear iso, no kin (clean)
    p.dmg.on=true; p.dmg.r=5.0; p.dmg.s=1.0; p.dmg.pD=0.0; p.dmg.Dc=0.99;

    // uniaxial strain path
    State su; zero(su);
    // confined: add hydrostatic tension on top of the same deviatoric drive so the
    // triaxiality sigmaH/sigma_eq is higher at comparable equivalent plastic strain
    State sc; zero(sc);
    for (int s=1;s<=6;s++){
      double a = 0.01*s;
      double epsU[6]={a,0,0,0,0,0};
      double epsC[6]={a, 0.5*a, 0.5*a, 0,0,0};  // adds hydrostatic part => higher triaxiality
      commit(p, epsU, su);
      commit(p, epsC, sc);
    }
    // compare damage per unit accumulated plastic strain
    double rateU = su.D/fmax(su.ebarP,1e-30);
    double rateC = sc.D/fmax(sc.ebarP,1e-30);
    printf("  uniaxial D/p=%.6g  confined D/p=%.6g\n", rateU, rateC);
    check(rateC > rateU, "higher triaxiality => faster damage", rateC, rateU);
  }

  // ---- IMPL-EX: frozen-D~ degrade gives stress/tangent = (1-D~)*effective, SPD,
  //      and D_out is the IMPLICIT value independent of the override ----
  printf("IMPL-EX (dScaleOverride): SPD (1-D~)*effective, implicit D returned:\n");
  {
    Params pon{K,G,5.0,3.0,25.0,8.0,1,{400},{80}};
    pon.dmg.on=true; pon.dmg.r=3.0; pon.dmg.s=1.0; pon.dmg.pD=0.0; pon.dmg.Dc=0.99;
    Params poff = pon; poff.dmg.on=false;       // effective (D-independent) reference
    State st; zero(st);
    double eps[6]={0.04,0,0,0,0,0};
    // effective stress + tangent (damage off)
    double se[6],De[6][6],eP[6],eb,al[MAXBACK][6],dg,Doff;
    returnMapDamaged(poff, eps, st.epsP, st.ebarP, st.alpha, st.D, se,De,eP,eb,al,dg,Doff);
    // implicit-damage run (to read the implicit D)
    double si[6],Di[6][6],Dimpl;
    returnMapDamaged(pon, eps, st.epsP, st.ebarP, st.alpha, st.D, si,Di,eP,eb,al,dg,Dimpl,0,-1.0);
    // IMPL-EX run with a frozen extrapolated D~ != Dimpl
    double Dtilde = 0.5*Dimpl;                   // arbitrary frozen extrapolation
    double sx[6],Dx[6][6],Dout_ix;
    returnMapDamaged(pon, eps, st.epsP, st.ebarP, st.alpha, st.D, sx,Dx,eP,eb,al,dg,Dout_ix,0,Dtilde);
    // D_out must be the IMPLICIT value, independent of the override
    check(fabs(Dout_ix - Dimpl) <= 1e-14, "IMPL-EX returns implicit D", Dout_ix, Dimpl);
    // stress == (1-D~)*effective ; tangent == (1-D~)*effective (SPD, no rank-one term)
    double errS=0, errT=0, scaleT=0;
    for (int i=0;i<6;i++) errS=fmax(errS, fabs(sx[i]-(1.0-Dtilde)*se[i]));
    for (int I=0;I<6;I++) for(int J=0;J<6;J++){
      scaleT=fmax(scaleT, fabs(De[I][J]));
      errT=fmax(errT, fabs(Dx[I][J]-(1.0-Dtilde)*De[I][J]));
    }
    check(errS <= 1e-9*(fabs(se[0])+1.0), "IMPL-EX stress = (1-D~)*effective", errS, 0.0);
    check(errT <= 1e-9*scaleT, "IMPL-EX tangent = (1-D~)*Dtan_eff (SPD)", errT, 0.0);
  }

  printf("\n%s (%d failure%s)\n", fails? "FAILED":"ALL PASS", fails, fails==1?"":"s");
  return fails ? 1 : 0;
}
