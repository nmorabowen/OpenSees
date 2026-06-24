// Standalone C++ check of LadrunoContactProjection::normalGeomTangent (B3): builds the
// FULL frictionless NTS tangent K = kn(B⊗B) + K_geom from the header, FD-checks it
// against the assembled residual, and confirms the flat-master slave block is 0. No
// OpenSees, no build — g++ -I SRC/domain/contact this file. Mirrors proto_b3_normal_tangent.py.
#include <cstdio>
#include <cmath>
#include "LadrunoContactProjection.h"
namespace P = LadrunoContactProjection;

static int NPS;
static void gatherXxs(const double Xin[4][3], const double xs[3], double q[15]) {
    for (int d = 0; d < 3; d++) q[d] = xs[d];
    for (int i = 0; i < NPS; i++) for (int d = 0; d < 3; d++) q[3*(1+i)+d] = Xin[i][d];
}
static void splitq(const double q[15], double X[4][3], double xs[3]) {
    for (int d = 0; d < 3; d++) xs[d] = q[d];
    for (int i = 0; i < NPS; i++) for (int d = 0; d < 3; d++) X[i][d] = q[3*(1+i)+d];
}

// residual r = B·tn (penalty), B = [n | −N_i n], tn = kn·<−gN>+ ; 0 if not active.
static void residual(int nps, const double X[4][3], const double xs[3],
                     const double refDir[3], double kn, double r[15]) {
    int ndof = 3*(1+nps);
    for (int k = 0; k < ndof; k++) r[k] = 0.0;
    double gap, n[3], N[4];
    if (!P::evalSegment(nps, X, xs, refDir, gap, n, N)) return;   // not penetrating/oob
    double tn = kn * (-gap);
    for (int d = 0; d < 3; d++) r[d] = n[d]*tn;
    for (int i = 0; i < nps; i++) for (int d = 0; d < 3; d++) r[3*(1+i)+d] = -N[i]*n[d]*tn;
}

// full analytic K = kn·B⊗B + normalGeomTangent
static void tangentAnalytic(int nps, const double X[4][3], const double xs[3],
                            const double refDir[3], double kn, double K[15][15]) {
    int ndof = 3*(1+nps);
    for (int a = 0; a < ndof; a++) for (int b = 0; b < ndof; b++) K[a][b] = 0.0;
    double gap, n[3], N[4];
    if (!P::evalSegment(nps, X, xs, refDir, gap, n, N)) return;
    double B[15] = {0};
    for (int d = 0; d < 3; d++) B[d] = n[d];
    for (int i = 0; i < nps; i++) for (int d = 0; d < 3; d++) B[3*(1+i)+d] = -N[i]*n[d];
    for (int a = 0; a < ndof; a++) for (int b = 0; b < ndof; b++) K[a][b] += kn*B[a]*B[b];
    double Kg[15][15];
    if (P::normalGeomTangent(nps, X, xs, refDir, kn, Kg))
        for (int a = 0; a < ndof; a++) for (int b = 0; b < ndof; b++) K[a][b] += Kg[a][b];
}

static double check(int nps, const double X[4][3], const double xs[3],
                    const double refDir[3], double kn, const char *label) {
    NPS = nps;
    int ndof = 3*(1+nps);
    // FD tangent K = −∂r/∂q
    double q0[15]; gatherXxs(X, xs, q0);
    double Kfd[15][15];
    double h = 1e-7;
    for (int j = 0; j < ndof; j++) {
        double qp[15], qm[15];
        for (int k = 0; k < ndof; k++) { qp[k]=q0[k]; qm[k]=q0[k]; }
        qp[j]+=h; qm[j]-=h;
        double Xp[4][3], xsp[3], Xm[4][3], xsm[3], rp[15], rm[15];
        splitq(qp, Xp, xsp); splitq(qm, Xm, xsm);
        residual(nps, Xp, xsp, refDir, kn, rp);
        residual(nps, Xm, xsm, refDir, kn, rm);
        for (int i = 0; i < ndof; i++) Kfd[i][j] = -(rp[i]-rm[i])/(2*h);
    }
    double Kan[15][15]; tangentAnalytic(nps, X, xs, refDir, kn, Kan);
    double num=0, den=0, sym=0, kref=0;
    for (int a = 0; a < ndof; a++) for (int b = 0; b < ndof; b++) {
        double e = Kfd[a][b]-Kan[a][b]; num+=e*e; den+=Kfd[a][b]*Kfd[a][b];
        double s = Kan[a][b]-Kan[b][a]; sym+=s*s; kref+=Kan[a][b]*Kan[a][b];
    }
    double rel = std::sqrt(num/(den+1e-300));
    double asym = std::sqrt(sym/(kref+1e-300));
    // slave-block geom norm
    double Kg[15][15]; double hss=0;
    bool act = P::normalGeomTangent(nps, X, xs, refDir, kn, Kg);
    if (act) for (int a=0;a<3;a++) for (int b=0;b<3;b++) hss+=Kg[a][b]*Kg[a][b];
    printf("[%s] analytic-vs-FD=%.2e  sym=%.2e  |Kgeom_ss|=%.3e\n",
           label, rel, asym, std::sqrt(hss));
    return rel;
}

int main() {
    double refdir[3] = {0,0,1};
    double kn = 1e6;
    int bad = 0;

    // flat quad, fixed-equivalent slave-block must be 0
    double Xf[4][3] = {{-.5,-.5,0},{.5,-.5,0},{.5,.5,0},{-.5,.5,0}};
    double xsf[3] = {0.07,-0.04,-0.02};
    if (check(4, Xf, xsf, refdir, kn, "flat quad") > 1e-5) bad++;

    // warped (curved) quad: nodes 0,2 +w ; 1,3 −w
    double w=0.15;
    double Xc[4][3] = {{-.5,-.5,+w},{.5,-.5,-w},{.5,.5,+w},{-.5,.5,-w}};
    double xsc[3] = {0.05,-0.03,-0.02};
    if (check(4, Xc, xsc, refdir, kn, "curved quad") > 1e-5) bad++;

    // tri-3
    double Xt[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,0}};
    double xst[3] = {0.3,0.3,-0.01};
    if (check(3, Xt, xst, refdir, kn, "tri-3") > 1e-5) bad++;

    printf(bad==0 ? "B3 CPP STANDALONE: PASS\n" : "B3 CPP STANDALONE: FAIL (%d)\n", bad);
    return bad;
}
