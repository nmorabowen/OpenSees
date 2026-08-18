// BYTE-IDENTITY REFUTER (ADR-85 T2 pre-merge condition).
// Adversarial: TRY TO PROVE the shipped 3D friction kernel changed in this PR.
// Compiles the PRE-PR (HEAD~1) and POST-PR headers into two separate namespaces and
// compares frictionReturnMap / frictionTangentBlock BIT-FOR-BIT over randomized
// inputs, including the branch-boundary cases (cap<=0 free slip, exact stick/slip
// threshold, tauMax binding, separated N<=0). Any single differing bit refutes the
// "pure append => 3D untouched" claim.
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace OLDK { namespace LadrunoFrictionKernel {
#include "old_body.inc"
}}
namespace NEWK { namespace LadrunoFrictionKernel {
#include "new_body.inc"
}}

static uint64_t bits(double x){ uint64_t u; std::memcpy(&u,&x,8); return u; }

int main(){
    uint64_t st = 0x9E3779B97F4A7C15ULL;
    auto rnd=[&](double lo,double hi){ st^=st<<13; st^=st>>7; st^=st<<17;
        return lo+(hi-lo)*((double)(st>>11)/9007199254740992.0); };
    long long n=0, diffs=0;
    for (int i=0;i<400000;i++){
        double gT[3],gp[3],nn[3];
        for(int d=0;d<3;d++){ gT[d]=rnd(-1,1)*std::pow(10.0,rnd(-8,2));
                              gp[d]=rnd(-1,1)*std::pow(10.0,rnd(-8,2)); }
        double th=rnd(-3.14159,3.14159), ph=rnd(-1.5707,1.5707);
        nn[0]=std::cos(ph)*std::cos(th); nn[1]=std::cos(ph)*std::sin(th); nn[2]=std::sin(ph);
        double N   = (rnd(0,1)<0.2)?((rnd(0,1)<0.5)?0.0:-std::pow(10.0,rnd(-2,3)))
                                   : std::pow(10.0,rnd(-2,5));
        double kn  = std::pow(10.0,rnd(1,9));
        double kt  = std::pow(10.0,rnd(1,9));
        double mu  = (rnd(0,1)<0.15)?0.0:rnd(0,1.5);
        double coh = (rnd(0,1)<0.35)?std::pow(10.0,rnd(-3,3)):0.0;
        double tmx = (rnd(0,1)<0.35)?std::pow(10.0,rnd(-3,4)):0.0;
        bool   cons= (rnd(0,1)<0.5);

        double tO[3],gO[3],tN[3],gN[3];
        bool sO = OLDK::LadrunoFrictionKernel::frictionReturnMap(gT,gp,N,kt,mu,tO,gO,coh,tmx);
        bool sN = NEWK::LadrunoFrictionKernel::frictionReturnMap(gT,gp,N,kt,mu,tN,gN,coh,tmx);
        if (sO!=sN) diffs++;
        for(int d=0;d<3;d++){ if(bits(tO[d])!=bits(tN[d])) diffs++;
                              if(bits(gO[d])!=bits(gN[d])) diffs++; }

        double KO[3][3],KN2[3][3];
        OLDK::LadrunoFrictionKernel::frictionTangentBlock(gT,gp,nn,N,kn,kt,mu,cons,KO,coh,tmx);
        NEWK::LadrunoFrictionKernel::frictionTangentBlock(gT,gp,nn,N,kn,kt,mu,cons,KN2,coh,tmx);
        for(int a=0;a<3;a++) for(int b=0;b<3;b++)
            if(bits(KO[a][b])!=bits(KN2[a][b])) diffs++;

        double cO = OLDK::LadrunoFrictionKernel::frictionCap(N,mu,coh,tmx);
        double cN = NEWK::LadrunoFrictionKernel::frictionCap(N,mu,coh,tmx);
        if(bits(cO)!=bits(cN)) diffs++;
        n++;
    }
    printf("configs=%lld  DIFFERING BITS/BRANCHES=%lld\n", n, diffs);
    printf("%s\n", diffs==0
        ? "REFUTATION FAILED: the 3D friction kernel is BIT-IDENTICAL across this PR."
        : "*** REFUTED: the 3D friction kernel CHANGED -- FREEZE ***");
    return diffs==0?0:1;
}
