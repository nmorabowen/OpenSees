// Standalone g++ gate for Phase-3b crack-band regularization in LadrunoRCKernel.h.
// Mirrors tests/_testbed/rc_shell_ref.py run_R1: the regularized specific fracture energy must
// scale as G_f0*(lch_ref/lch) so the PHYSICAL energy g*lch is mesh-objective (constant).
//   g++ -std=c++17 -I SRC/material/nD rc_reg_gpp.cpp -o rc_reg && ./rc_reg
#include "LadrunoRCKernel.h"
#include <cstdio>
#include <cmath>
using namespace ladruno_rc_kernel;

int main()
{
    const double E = 30000.0;
    // softening tension backbone: ft=3 @ 1e-4, down to 0.5 @ 1e-3, 0.0 @ 4e-3 (d via 1-y/q)
    double Te[4] = { 0.0, 0.0001, 0.0010, 0.0040 };
    double Ts[4] = { 0.0, 3.0,    0.5,    0.0 };
    double Td[4] = { 0.0, 0.0,    1.0 - 0.5/5.0, 0.999 };
    Backbone b0; buildBackbone(b0, E, Te, Ts, Td, 4);

    int bd, p1, p2;
    double g0 = fractureEnergy(b0, E, &bd, &p1, &p2);
    double lch_ref = 50.0;
    double phys_ref = g0 * lch_ref;
    printf("G_f0(specific)=%.6f  physical=g*lch=%.5f (must stay constant)  bounded=%d\n", g0, phys_ref, bd);
    if (!bd) { printf("REG FAIL: backbone not bounded\n"); return 1; }

    int nbad = 0;
    double lchs[4] = { 50.0, 25.0, 100.0, 200.0 };
    for (int j = 0; j < 4; ++j) {
        double lch = lchs[j];
        Backbone b = b0;                              // fresh copy of the un-regularized backbone
        regularize(b, lch, lch_ref, E);
        int bd2, q1, q2;
        double g = fractureEnergy(b, E, &bd2, &q1, &q2);
        double want = g0 * (lch_ref / lch);
        double phys = g * lch;
        double e_spec = fabs(g - want) / (want + 1e-30);
        double e_phys = fabs(phys - phys_ref) / (phys_ref + 1e-30);
        printf("  lch=%6.1f  g_reg=%.6f  want=%.6f  g*lch=%.5f  |phys-ref|/ref=%.2e\n",
               lch, g, want, phys, e_phys);
        if (e_spec > 2e-3 || e_phys > 2e-3) nbad++;
    }
    // reduce check: lch == lch_ref must be a no-op (lch_scale==1 quick return)
    {
        Backbone b = b0; regularize(b, lch_ref, lch_ref, E);
        double dmax = 0.0;
        for (int i = 0; i < b0.n; ++i) { double d = fabs(b.x[i]-b0.x[i]) + fabs(b.q[i]-b0.q[i]); if (d>dmax) dmax=d; }
        printf("  lch==lch_ref no-op: max|delta|=%.2e\n", dmax);
        if (dmax > 1e-12) nbad++;
    }
    if (nbad == 0) { printf("GPP REG GATE: PASS\n"); return 0; }
    printf("GPP REG GATE: FAIL (%d)\n", nbad);
    return 1;
}
