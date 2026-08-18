// PRE-T2 EXPERIMENT (ADR-85 How/4): reuse-vs-scalar for the 2D friction return map.
// Path A = the SHIPPED 3D LadrunoFrictionKernel evaluated on plane-constrained configs.
// Path B = an independent scalar (1D tangent space) return map.
// Dumps both at %.17g so a numpy referee can diff at the bit level.
#include <cstdio>
#include <cmath>
#include <cstdint>
#include "LadrunoFrictionKernel.h"

// ---- Path B: independent scalar return map (1D tangent space) --------------
// s      = signed tangential slip measure along t_hat
// sp     = committed signed plastic slip
// returns signed friction force applied to the slave (already negated), and
// the trial plastic slip + the scalar tangent.
static bool returnMap1D(double s, double sp, double N, double kt, double mu,
                        double &tFric, double &spTrial, double cohesion, double tauMax)
{
    double tTtr = kt * (s - sp);
    double capC = (N > 0.0) ? (mu * N + cohesion) : 0.0;
    double cap  = (tauMax > 0.0 && tauMax < capC) ? tauMax : capC;
    double nrm  = std::fabs(tTtr);                 // <-- exact, no sqrt
    if (cap <= 0.0) { tFric = 0.0; spTrial = s; return true; }
    if (nrm <= cap) { tFric = -tTtr; spTrial = sp; return false; }
    double sgn  = (tTtr > 0.0) ? 1.0 : -1.0;       // <-- exact, no normalize
    double dlam = (nrm - cap) / kt;
    tFric   = -cap * sgn;
    spTrial = sp + dlam * sgn;
    return true;
}
static double tangent1D(double s, double sp, double N, double kn, double kt,
                        double mu, bool consistent, double cohesion, double tauMax)
{
    double tTtr = kt * (s - sp);
    double capC = (N > 0.0) ? (mu * N + cohesion) : 0.0;
    bool capped = (tauMax > 0.0 && tauMax < capC);
    double cap  = capped ? tauMax : capC;
    double nrm  = std::fabs(tTtr);
    if (cap <= 0.0) return 0.0;
    if (nrm <= cap) return kt;                     // stick
    return 0.0;                                    // slip: 1D tangent space => EXACTLY 0
    // (the consistent d_TN coupling is a tangential-normal cross term, reported separately)
}

static uint64_t bits(double x){ uint64_t u; __builtin_memcpy(&u,&x,8); return u; }

int main(){
    // deterministic LCG so the artifact is reproducible
    uint64_t st = 88172645463325252ULL;
    auto rnd=[&](double lo,double hi){ st^=st<<13; st^=st>>7; st^=st<<17;
        return lo+(hi-lo)*((double)(st>>11)/9007199254740992.0); };

    printf("# case theta kt kn N mu coh tmax s sp | A_t A_sp A_Ktt A_leak_z A_leak_n A_slip | B_t B_sp B_Ktt B_slip\n");
    int nThresh=0, nBranchDiff=0, ncase=0;
    for (int i=0;i<200000;i++){
        double theta = rnd(-3.14159265358979,3.14159265358979);
        double kt=std::pow(10.0,rnd(1,7)), kn=std::pow(10.0,rnd(1,7));
        double N = (rnd(0,1)<0.15)?0.0:std::pow(10.0,rnd(-2,5));
        double mu=rnd(0,1.2);
        double coh=(rnd(0,1)<0.3)?std::pow(10.0,rnd(-2,3)):0.0;
        double tmax=(rnd(0,1)<0.3)?std::pow(10.0,rnd(-2,4)):0.0;
        double sp=rnd(-1,1)*std::pow(10.0,rnd(-6,1));
        double s;
        // 40%: park the trial traction right AT the cone radius (threshold probing)
        double capC=(N>0.0)?(mu*N+coh):0.0;
        double cap=(tmax>0.0&&tmax<capC)?tmax:capC;
        if (rnd(0,1)<0.4 && cap>0.0){
            double target=cap/kt*((rnd(0,1)<0.5)?1.0:-1.0);
            // perturb by a few ulp to sit astride the stick/slip boundary
            double pert=1.0+ (rnd(-1,1))*4.0*2.220446049250313e-16;
            s = sp + target*pert; nThresh++;
        } else {
            s = sp + rnd(-1,1)*std::pow(10.0,rnd(-8,2));
        }

        double c=std::cos(theta), sn=std::sin(theta);
        double n[3]={c,sn,0.0};
        double th[3]={-sn,c,0.0};                       // t_hat = perp(n)
        double gT[3]={s*th[0], s*th[1], 0.0};
        double gp[3]={sp*th[0],sp*th[1],0.0};

        // ---- Path A: shipped 3D kernel, plane-constrained ----
        double tF[3],gpt[3],Kss[3][3];
        bool slipA=LadrunoFrictionKernel::frictionReturnMap(gT,gp,N,kt,mu,tF,gpt,coh,tmax);
        LadrunoFrictionKernel::frictionTangentBlock(gT,gp,n,N,kn,kt,mu,false,Kss,coh,tmax);
        double A_t  = tF[0]*th[0]+tF[1]*th[1];          // project back onto t_hat
        double A_sp = gpt[0]*th[0]+gpt[1]*th[1];
        double A_Ktt= th[0]*(Kss[0][0]*th[0]+Kss[0][1]*th[1]+Kss[0][2]*th[2])
                    + th[1]*(Kss[1][0]*th[0]+Kss[1][1]*th[1]+Kss[1][2]*th[2])
                    + th[2]*(Kss[2][0]*th[0]+Kss[2][1]*th[1]+Kss[2][2]*th[2]);
        double A_lz = std::fabs(tF[2])+std::fabs(gpt[2]);   // out-of-plane leakage
        double A_ln = std::fabs(tF[0]*n[0]+tF[1]*n[1]);     // normal-direction leakage

        // ---- Path B: independent scalar ----
        double B_t,B_sp; bool slipB=returnMap1D(s,sp,N,kt,mu,B_t,B_sp,coh,tmax);
        double B_Ktt=tangent1D(s,sp,N,kn,kt,mu,false,coh,tmax);

        if (slipA!=slipB) nBranchDiff++;
        ncase++;
        if (i<4000)
          printf("%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g | %.17g %.17g %.17g %.17g %.17g %d | %.17g %.17g %.17g %d\n",
             i,theta,kt,kn,N,mu,coh,tmax,s,sp,
             A_t,A_sp,A_Ktt,A_lz,A_ln,(int)slipA, B_t,B_sp,B_Ktt,(int)slipB);
    }
    fprintf(stderr,"cases=%d threshold_probes=%d STICK_SLIP_BRANCH_DISAGREEMENTS=%d\n",
            ncase,nThresh,nBranchDiff);
    return 0;
}
