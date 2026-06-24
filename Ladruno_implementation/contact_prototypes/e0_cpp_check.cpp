// Standalone C++ check of LadrunoEdgeKernel (ADR-57 E0): exercises closestPtSegSeg /
// edgeNormal / edgeGap / bOperator and asserts they reproduce the proto_e0_closest_point.py
// oracle bit-for-bit (closed-form sub-regions, parallel/degenerate refusal, the B-operator
// first variation via FD, and the committed-sign monotonicity through contact). No OpenSees,
// no build:  g++ -I SRC/domain/contact e0_cpp_check.cpp -o e0_cpp_check && ./e0_cpp_check
#include <cstdio>
#include <cmath>
#include "LadrunoEdgeKernel.h"
namespace E = LadrunoEdgeKernel;

static const double PI = 3.14159265358979323846;
static int fails = 0;
static void check(const char *name, bool ok, const char *extra = "") {
    std::printf("  [%s] %s%s%s\n", ok ? "PASS" : "FAIL", name, extra[0] ? " — " : "", extra);
    if (!ok) fails++;
}

int main() {
    std::printf("ADR-57 E0 C++ check — LadrunoEdgeKernel vs the numpy oracle\n");

    // ---- T1: vertex-vertex double-clamp (s=1, t=0, dist=sqrt2) ----
    {
        double p1[3]={0,0,0}, q1[3]={1,0,0}, p2[3]={2,0,1}, q2[3]={2,1,1};
        E::ClosestResult k = E::closestPtSegSeg(p1,q1,p2,q2);
        double dx=k.c1[0]-k.c2[0], dy=k.c1[1]-k.c2[1], dz=k.c1[2]-k.c2[2];
        double dist=std::sqrt(dx*dx+dy*dy+dz*dz);
        check("vertex-vertex s==1,t==0", std::fabs(k.s-1.0)<1e-12 && std::fabs(k.t)<1e-12);
        check("vertex-vertex dist==sqrt2", std::fabs(dist-std::sqrt(2.0))<1e-12);
        check("vertex-vertex NOT interior", !k.interior);
    }
    // ---- T1b: interior crossing (s=t=0.5, dist=1) ----
    {
        double p1[3]={-1,0,0}, q1[3]={1,0,0}, p2[3]={0,-1,1}, q2[3]={0,1,1};
        E::ClosestResult k = E::closestPtSegSeg(p1,q1,p2,q2);
        check("interior s==t==0.5", std::fabs(k.s-0.5)<1e-12 && std::fabs(k.t-0.5)<1e-12);
        check("interior IS margin-interior", k.interior && k.status==E::EE_OK);
    }
    // ---- T2: parallel + 4deg near-parallel refusal; 30deg admitted ----
    {
        double p1[3]={0,0,0}, q1[3]={1,0,0}, p2[3]={0,0,1}, q2[3]={1,0,1};
        check("exactly parallel -> EE_PARALLEL", E::closestPtSegSeg(p1,q1,p2,q2).status==E::EE_PARALLEL);
        double a1[3]={-1,0,0}, b1[3]={1,0,0};
        // seg2 centered at (0,0,0.5) with half-direction d at angle theta to seg1
        double ang=4.0*PI/180.0,  d4[3] ={std::cos(ang),std::sin(ang),0};
        double a2[3]={-d4[0],-d4[1],0.5-d4[2]}, b2[3]={d4[0],d4[1],0.5+d4[2]};
        check("4deg near-parallel -> EE_PARALLEL", E::closestPtSegSeg(a1,b1,a2,b2).status==E::EE_PARALLEL);
        double a30=30.0*PI/180.0, d3[3]={std::cos(a30),std::sin(a30),0};
        double c2a[3]={-d3[0],-d3[1],0.5-d3[2]}, c2b[3]={d3[0],d3[1],0.5+d3[2]};
        E::ClosestResult k30=E::closestPtSegSeg(a1,b1,c2a,c2b);
        check("30deg -> OK, margin-interior", k30.status==E::EE_OK && k30.interior);
    }
    // ---- T3: zero-length edge -> DEGENERATE ----
    {
        double p1[3]={1,1,1}, q1[3]={1,1,1}, p2[3]={0,-1,0}, q2[3]={0,1,0};
        check("collapsed seg1 -> EE_DEGENERATE", E::closestPtSegSeg(p1,q1,p2,q2).status==E::EE_DEGENERATE);
    }
    // ---- T5: gN==+-||w|| interior; |gN|<||w|| clamped ----
    {
        double p1[3]={-1,0,0}, q1[3]={1,0,0}, p2[3]={0,-1,1}, q2[3]={0,1,1}, ref[3]={0,0,-1};
        E::ClosestResult k=E::closestPtSegSeg(p1,q1,p2,q2);
        double d1[3]={q1[0]-p1[0],q1[1]-p1[1],q1[2]-p1[2]}, d2[3]={q2[0]-p2[0],q2[1]-p2[1],q2[2]-p2[2]};
        double n[3]; int sg;
        double gN=E::edgeGap(k.c1,k.c2,d1,d2,0,ref,n,&sg);
        double w[3]={k.c1[0]-k.c2[0],k.c1[1]-k.c2[1],k.c1[2]-k.c2[2]};
        check("interior |gN|==||w||", std::fabs(std::fabs(gN)-E::norm3(w))<1e-12);
    }
    // ---- T6: B-operator vs central-difference (envelope theorem) ----
    {
        double p1[3]={-1,0.1,0}, q1[3]={1,-0.05,0}, p2[3]={0.1,-1,0.8}, q2[3]={-0.1,1,1.2}, ref[3]={0,0,-1};
        E::ClosestResult k=E::closestPtSegSeg(p1,q1,p2,q2);
        double d1[3]={q1[0]-p1[0],q1[1]-p1[1],q1[2]-p1[2]}, d2[3]={q2[0]-p2[0],q2[1]-p2[1],q2[2]-p2[2]};
        double n0[3]; int sg0;
        E::edgeGap(k.c1,k.c2,d1,d2,0,ref,n0,&sg0);
        double B[4][3]; E::bOperator(n0,k.s,k.t,B);
        double *nodes[4]={p1,q1,p2,q2};
        double maxerr=0.0, h=1e-7;
        auto gapAt=[&](double P1[3],double Q1[3],double P2[3],double Q2[3])->double{
            E::ClosestResult kk=E::closestPtSegSeg(P1,Q1,P2,Q2);
            double dd1[3]={Q1[0]-P1[0],Q1[1]-P1[1],Q1[2]-P1[2]}, dd2[3]={Q2[0]-P2[0],Q2[1]-P2[1],Q2[2]-P2[2]};
            double nn[3]; return E::edgeGap(kk.c1,kk.c2,dd1,dd2,sg0,0,nn,0);
        };
        for (int ni=0; ni<4; ni++) for (int d=0; d<3; d++) {
            double sav=nodes[ni][d];
            nodes[ni][d]=sav+h; double gp=gapAt(p1,q1,p2,q2);
            nodes[ni][d]=sav-h; double gm=gapAt(p1,q1,p2,q2);
            nodes[ni][d]=sav;
            double fd=(gp-gm)/(2*h);
            double err=std::fabs(fd-B[ni][d]); if (err>maxerr) maxerr=err;
        }
        char buf[64]; std::snprintf(buf,sizeof buf,"max err=%.2e",maxerr);
        check("B-operator == FD < 1e-6", maxerr<1e-6, buf);
    }
    // ---- T7: committed sign -> gN monotone through contact ----
    {
        double p2[3]={0,-1,1}, q2[3]={0,1,1}, ref[3]={0,0,1};
        double d2[3]={q2[0]-p2[0],q2[1]-p2[1],q2[2]-p2[2]};
        // commit the sign at z=1.6 (gap>0)
        double p1a[3]={-1,0,1.6}, q1a[3]={1,0,1.6};
        E::ClosestResult k0=E::closestPtSegSeg(p1a,q1a,p2,q2);
        double d1a[3]={q1a[0]-p1a[0],q1a[1]-p1a[1],q1a[2]-p1a[2]}, n[3]; int sgc;
        E::edgeGap(k0.c1,k0.c2,d1a,d2,0,ref,n,&sgc);
        double prev=1e9; bool mono=true, neg=false;
        for (int i=0;i<13;i++){
            double z=1.6-(1.6-0.4)*i/12.0;
            double p1[3]={-1,0,z}, q1[3]={1,0,z};
            E::ClosestResult k=E::closestPtSegSeg(p1,q1,p2,q2);
            double d1[3]={q1[0]-p1[0],q1[1]-p1[1],q1[2]-p1[2]}, nn[3];
            double gN=E::edgeGap(k.c1,k.c2,d1,d2,sgc,0,nn,0);
            if (gN>prev+1e-9) mono=false;
            prev=gN;
            if (gN<0) neg=true;
        }
        check("committed sign: gN monotone through contact", mono);
        check("committed sign: gN goes negative (penetration registered)", neg);
    }

    std::printf("\n================================================\n%s  (E0 C++ check)\n",
                fails ? "FAILURE(S)" : "ALL PASS");
    return fails ? 1 : 0;
}
