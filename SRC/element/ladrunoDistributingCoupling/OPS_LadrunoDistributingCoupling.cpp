/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno: OPS parser for the RBE3 / distributing-coupling element (ADR 28).
//
//   element LadrunoDistributingCoupling tag refNode N i1 ... iN
//           [-w w1 ... wN]                 # weights (default equal); -area = deferred
//           [-k {Kt | auto}] [-kAlpha a] [-host eleTag]   # auto needs a representative -host
//           [-kr Kr]                       # rotational penalty (default derived K_t*ell^2)
//           [-enforce {penalty | al}]
//           [-bipenalty {-dtcr dt | -wcap beta}]
//           [-absolute]                    # opt out of initial-gap (offset) capture
//
//   The reference node (6-DOF in 3D, 3-DOF in 2D) is interpolated from the
//   independent node set (translations); a load at the reference distributes to
//   the set without stiffening it. -host names ONE representative element among
//   the independents' parents, used ONLY to scale -k auto / -wcap (RBE3 has no
//   single host element). Default penalty K_t = 1e12; K_r = K_t * (Sum w_i|r_i|^2/W)
//   unless -kr is given. -bipenalty (with -dtcr) bounds the explicit step of the
//   massless reference node. See ADR 28 for the formulation.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoDistributingCoupling.h>
#include <ID.h>
#include <Vector.h>
#include <Domain.h>
#include <Element.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>

void* OPS_LadrunoDistributingCoupling(void)
{
  int ndm = OPS_GetNDM();
  if (ndm != 2 && ndm != 3) {
    opserr << "WARNING LadrunoDistributingCoupling: model ndm must be 2 or 3\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING insufficient args\n"
           << "Want: element LadrunoDistributingCoupling tag refNode N i1..iN "
           << "[-w w1..wN] [-k {Kt|auto}] [-kAlpha a] [-host eleTag] [-kr Kr] "
           << "[-enforce {penalty|al}] [-bipenalty {-dtcr dt | -wcap beta}] [-absolute]\n";
    return 0;
  }

  int idata[3];                    // tag, refNode, N
  int n = 3;
  if (OPS_GetIntInput(&n, idata) < 0) {
    opserr << "WARNING LadrunoDistributingCoupling: invalid tag/refNode/N\n";
    return 0;
  }
  int tag = idata[0], refNode = idata[1], N = idata[2];
  if (N < 1) {
    opserr << "WARNING LadrunoDistributingCoupling: N (number of independent nodes) "
           << "must be >= 1; got " << N << "\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < N) {
    opserr << "WARNING LadrunoDistributingCoupling: need " << N
           << " independent node tags\n";
    return 0;
  }
  ID indep(N);
  for (int i = 0; i < N; i++) {
    int h; n = 1;
    if (OPS_GetIntInput(&n, &h) < 0) {
      opserr << "WARNING LadrunoDistributingCoupling: invalid independent node " << i << "\n";
      return 0;
    }
    indep(i) = h;
  }

  Vector weights(N);
  for (int i = 0; i < N; i++) weights(i) = 1.0;   // default equal weights

  double Kt = 1.0e12;
  bool ktAuto = false;
  double kAlpha = 1.0e3;
  int hostEleTag = -1;
  double Kr = 0.0;
  bool krUser = false;
  int enforce = 0;
  bool bipenalty = false;
  int bpMode = 0;
  double bpDt = 0.0, bpBeta = 0.0;
  bool bpBudgetSet = false;
  bool initGapCapture = true;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    if (strcmp(opt, "-w") == 0) {
      n = N;
      if (OPS_GetDoubleInput(&n, &weights(0)) < 0) {
        opserr << "WARNING LadrunoDistributingCoupling: -w wants " << N << " weights\n";
        return 0;
      }
      for (int i = 0; i < N; i++)
        if (weights(i) < 0.0) {
          opserr << "WARNING LadrunoDistributingCoupling: -w weight " << i
                 << " is negative (" << weights(i) << "); weights must be >= 0\n";
          return 0;
        }
    }
    else if (strcmp(opt, "-area") == 0) {
      opserr << "INFO LadrunoDistributingCoupling " << tag
             << ": -area (tributary-area weights) not implemented in v1; using -w "
             << "(default equal). Supply -w explicitly for area weights.\n";
    }
    else if (strcmp(opt, "-k") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoDistributingCoupling: -k wants a value or 'auto'\n";
        return 0;
      }
      char kTok[64];
      OPS_GetStringFromAll(kTok, sizeof(kTok));
      if (strcmp(kTok, "auto") == 0) ktAuto = true;
      else {
        char* endp = 0;
        double kv = strtod(kTok, &endp);
        if (endp == kTok || *endp != '\0') {
          opserr << "WARNING LadrunoDistributingCoupling: -k wants a number or 'auto', got '"
                 << kTok << "'\n";
          return 0;
        }
        Kt = kv;
      }
    }
    else if (strcmp(opt, "-kAlpha") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &kAlpha) < 0) {
        opserr << "WARNING LadrunoDistributingCoupling: -kAlpha wants a value\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-host") == 0) {
      int eleTag; n = 1;
      if (OPS_GetIntInput(&n, &eleTag) < 0) {
        opserr << "WARNING LadrunoDistributingCoupling: -host wants an element tag\n";
        return 0;
      }
      hostEleTag = eleTag;
    }
    else if (strcmp(opt, "-kr") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &Kr) < 0) {
        opserr << "WARNING LadrunoDistributingCoupling: -kr wants a value\n";
        return 0;
      }
      krUser = true;
    }
    else if (strcmp(opt, "-enforce") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoDistributingCoupling: -enforce wants penalty|al\n";
        return 0;
      }
      const char* mode = OPS_GetString();
      if (strcmp(mode, "penalty") == 0)  enforce = 0;
      else if (strcmp(mode, "al") == 0)  enforce = 1;
      else {
        opserr << "WARNING LadrunoDistributingCoupling: unknown -enforce '" << mode
               << "' (want penalty|al)\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-bipenalty") == 0) {
      bipenalty = true;
    }
    else if (strcmp(opt, "-dtcr") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &bpDt) < 0) {
        opserr << "WARNING LadrunoDistributingCoupling: -dtcr wants a target step\n";
        return 0;
      }
      bipenalty = true; bpMode = 0; bpBudgetSet = true;
    }
    else if (strcmp(opt, "-wcap") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &bpBeta) < 0) {
        opserr << "WARNING LadrunoDistributingCoupling: -wcap wants a frequency ratio beta\n";
        return 0;
      }
      bipenalty = true; bpMode = 1; bpBudgetSet = true;
    }
    else if (strcmp(opt, "-absolute") == 0 || strcmp(opt, "-noInitGap") == 0) {
      initGapCapture = false;
    }
    else {
      opserr << "WARNING LadrunoDistributingCoupling: unknown option '" << opt << "'\n";
      return 0;
    }
  }

  if (ktAuto && hostEleTag < 0) {
    opserr << "WARNING LadrunoDistributingCoupling: -k auto requires a representative "
           << "-host element (RBE3 has no single host); supply -host eleTag or a numeric -k\n";
    return 0;
  }
  if (bipenalty && !bpBudgetSet) {
    opserr << "WARNING LadrunoDistributingCoupling: -bipenalty needs a budget — "
           << "-dtcr <dt> (explicit) or -wcap <beta> (with -host)\n";
    return 0;
  }
  if (bipenalty && bpMode == 1 && hostEleTag < 0) {
    opserr << "WARNING LadrunoDistributingCoupling: -wcap requires a representative -host "
           << "element for omega_host; use -dtcr instead\n";
    return 0;
  }
  if (bipenalty && enforce == 1) {
    opserr << "WARNING LadrunoDistributingCoupling: -bipenalty ignored with -enforce al "
           << "(augmented Lagrangian needs no mass penalty)\n";
    bipenalty = false; bpBudgetSet = false;
  }

  Element* e = new LadrunoDistributingCoupling(tag, ndm, refNode, indep, weights, Kt,
                                               Kr, krUser, enforce,
                                               bipenalty, bpMode, bpDt, bpBeta,
                                               kAlpha, hostEleTag, ktAuto, initGapCapture);
  if (e == 0) {
    opserr << "WARNING LadrunoDistributingCoupling: could not create element\n";
    return 0;
  }
  return e;
}
