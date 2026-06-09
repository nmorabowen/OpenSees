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

// Ladruno: OPS parser for the RBE2 / kinematic-coupling element (ADR 29).
//
//   element LadrunoKinematicCoupling tag refNode N s1 ... sN
//           [-dof c1 ... cK]               # dependent components on each slave
//                                          #   1..ndm = translations, ndm+1..ndm+nrot = rotations
//                                          #   default = every DOF the slave possesses
//           [-k {Kt | auto}] [-kAlpha a] [-host eleTag]   # auto needs a representative -host
//           [-kr Kr]                       # rotational-tie penalty (default derived K_t*ell^2)
//           [-enforce {penalty | al}]
//           [-bipenalty {-dtcr dt | -wcap beta}]   # OFF by default (R is often massed)
//           [-absolute]                    # opt out of initial-gap (offset) capture
//
//   The reference node R RIGIDLY drives the slave set: u_i = u_R + theta_R x (x_i-x_R),
//   theta_i = theta_R where tied. Adds rigid penalty stiffness to the slaves. -host names
//   ONE representative element among the slaves' parents, used ONLY to scale -k auto /
//   -wcap. Default penalty K_t = 1e12; K_r = K_t * ell^2 unless -kr. See ADR 29.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoKinematicCoupling.h>
#include <ID.h>
#include <Domain.h>
#include <Element.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>

void* OPS_LadrunoKinematicCoupling(void)
{
  int ndm = OPS_GetNDM();
  if (ndm != 2 && ndm != 3) {
    opserr << "WARNING LadrunoKinematicCoupling: model ndm must be 2 or 3\n";
    return 0;
  }
  int nrot = (ndm == 3) ? 3 : 1;
  int cmax = ndm + nrot;

  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING insufficient args\n"
           << "Want: element LadrunoKinematicCoupling tag refNode N s1..sN "
           << "[-dof c1..cK] [-k {Kt|auto}] [-kAlpha a] [-host eleTag] [-kr Kr] "
           << "[-enforce {penalty|al}] [-bipenalty {-dtcr dt | -wcap beta}] [-absolute]\n";
    return 0;
  }

  int idata[3];                    // tag, refNode, N
  int n = 3;
  if (OPS_GetIntInput(&n, idata) < 0) {
    opserr << "WARNING LadrunoKinematicCoupling: invalid tag/refNode/N\n";
    return 0;
  }
  int tag = idata[0], refNode = idata[1], N = idata[2];
  if (N < 1) {
    opserr << "WARNING LadrunoKinematicCoupling: N (number of slave nodes) must be >= 1; got "
           << N << "\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < N) {
    opserr << "WARNING LadrunoKinematicCoupling: need " << N << " slave node tags\n";
    return 0;
  }
  ID slaves(N);
  for (int i = 0; i < N; i++) {
    int h; n = 1;
    if (OPS_GetIntInput(&n, &h) < 0) {
      opserr << "WARNING LadrunoKinematicCoupling: invalid slave node " << i << "\n";
      return 0;
    }
    slaves(i) = h;
  }

  ID dofSel(0);                    // empty = default all available
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
    if (strcmp(opt, "-dof") == 0) {
      // greedily read following integer tokens as components until a non-integer/flag.
      int nsel = 0;
      while (OPS_GetNumRemainingInputArgs() > 0) {
        char tok[64];
        OPS_GetStringFromAll(tok, sizeof(tok));
        char* endp = 0;
        long cv = strtol(tok, &endp, 10);
        if (endp == tok || *endp != '\0') {       // not an integer — it's the next flag
          OPS_ResetCurrentInputArg(-1);
          break;
        }
        if (cv < 1 || cv > cmax) {
          opserr << "WARNING LadrunoKinematicCoupling: -dof component " << cv
                 << " out of range for " << ndm << "D; valid 1.." << cmax << "\n";
          return 0;
        }
        dofSel[nsel++] = (int)cv;
      }
      if (nsel == 0) {
        opserr << "WARNING LadrunoKinematicCoupling: -dof needs at least one component\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-k") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoKinematicCoupling: -k wants a value or 'auto'\n";
        return 0;
      }
      char kTok[64];
      OPS_GetStringFromAll(kTok, sizeof(kTok));
      if (strcmp(kTok, "auto") == 0) ktAuto = true;
      else {
        char* endp = 0;
        double kv = strtod(kTok, &endp);
        if (endp == kTok || *endp != '\0') {
          opserr << "WARNING LadrunoKinematicCoupling: -k wants a number or 'auto', got '"
                 << kTok << "'\n";
          return 0;
        }
        Kt = kv;
      }
    }
    else if (strcmp(opt, "-kAlpha") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &kAlpha) < 0) {
        opserr << "WARNING LadrunoKinematicCoupling: -kAlpha wants a value\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-host") == 0) {
      int eleTag; n = 1;
      if (OPS_GetIntInput(&n, &eleTag) < 0) {
        opserr << "WARNING LadrunoKinematicCoupling: -host wants an element tag\n";
        return 0;
      }
      hostEleTag = eleTag;
    }
    else if (strcmp(opt, "-kr") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &Kr) < 0) {
        opserr << "WARNING LadrunoKinematicCoupling: -kr wants a value\n";
        return 0;
      }
      krUser = true;
    }
    else if (strcmp(opt, "-enforce") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoKinematicCoupling: -enforce wants penalty|al\n";
        return 0;
      }
      const char* mode = OPS_GetString();
      if (strcmp(mode, "penalty") == 0)  enforce = 0;
      else if (strcmp(mode, "al") == 0)  enforce = 1;
      else {
        opserr << "WARNING LadrunoKinematicCoupling: unknown -enforce '" << mode
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
        opserr << "WARNING LadrunoKinematicCoupling: -dtcr wants a target step\n";
        return 0;
      }
      bipenalty = true; bpMode = 0; bpBudgetSet = true;
    }
    else if (strcmp(opt, "-wcap") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &bpBeta) < 0) {
        opserr << "WARNING LadrunoKinematicCoupling: -wcap wants a frequency ratio beta\n";
        return 0;
      }
      bipenalty = true; bpMode = 1; bpBudgetSet = true;
    }
    else if (strcmp(opt, "-absolute") == 0 || strcmp(opt, "-noInitGap") == 0) {
      initGapCapture = false;
    }
    else {
      opserr << "WARNING LadrunoKinematicCoupling: unknown option '" << opt << "'\n";
      return 0;
    }
  }

  if (ktAuto && hostEleTag < 0) {
    opserr << "WARNING LadrunoKinematicCoupling: -k auto requires a representative -host "
           << "element (no single host for a node set); supply -host eleTag or a numeric -k\n";
    return 0;
  }
  if (bipenalty && !bpBudgetSet) {
    opserr << "WARNING LadrunoKinematicCoupling: -bipenalty needs a budget — "
           << "-dtcr <dt> (explicit) or -wcap <beta> (with -host)\n";
    return 0;
  }
  if (bipenalty && bpMode == 1 && hostEleTag < 0) {
    opserr << "WARNING LadrunoKinematicCoupling: -wcap requires a representative -host "
           << "element for omega_host; use -dtcr instead\n";
    return 0;
  }
  if (bipenalty && enforce == 1) {
    opserr << "WARNING LadrunoKinematicCoupling: -bipenalty ignored with -enforce al "
           << "(augmented Lagrangian needs no mass penalty)\n";
    bipenalty = false; bpBudgetSet = false;
  }

  Element* e = new LadrunoKinematicCoupling(tag, ndm, refNode, slaves, dofSel, Kt,
                                            Kr, krUser, enforce,
                                            bipenalty, bpMode, bpDt, bpBeta,
                                            kAlpha, hostEleTag, ktAuto, initGapCapture);
  if (e == 0) {
    opserr << "WARNING LadrunoKinematicCoupling: could not create element\n";
    return 0;
  }
  return e;
}
