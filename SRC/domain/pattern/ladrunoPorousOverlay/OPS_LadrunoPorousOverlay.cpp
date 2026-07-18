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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 07/2026
//
// Parser/factory for LadrunoPorousOverlay — the persistent-fluid staggered u-p
// overlay LoadPattern (ADR 73, PATTERN_TAG 33022). The pattern-command
// dispatcher (OpenSeesPatternCommands.cpp) has already consumed the type token
// "LadrunoPorousOverlay", so the arg cursor sits on $tag here (the H5DRM wiring).
//
//   pattern LadrunoPorousOverlay $tag \
//       -region  $e1 $e2 ...                 ;# solid cells (flat list / {*}$list)
//       -Kf $Kf -rhoF $rhof -poro $n \
//       -perm $k1 $k2 <$k3>                  ;# k̄ = k_hydraulic/gammaW, ndm values
//       -moduli $E $nu                       ;# REQUIRED (v1): L / α-stab moduli
//       -drained $nd1 $nd2 ...               ;# p-fixed boundary (>=1 region node)
//       <-alpha $biotAlpha> <-Ks $Ks> <-thick $t> \
//       <-layer $e1 ... -perm $k1 $k2 <$k3> <-poro $n> <-moduli $E $nu>> ...
//       <-stab auto <$a0> | off | $alpha> \
//       <-fsL classic | oedometric | zero | $scale> \  ;# zero = explicit lane (P3)
//       <-staticMode hold | steady> \
//       <-onRemoval keep | drain $kF> \
//       <-fluidUpdate implicit | explicit> \  ;# explicit = fatal-NYI (P3b)
//       <-subcycle auto | $N> \
//       <-pInit steady | hydrostatic $gw <$zw> | $nd $val ...> \
//       <-fluidBody $b1 $b2 <$b3>> <-dynSeepage on|off> \
//       <-record $file <$everyN>>
//
// UNKNOWN FLAGS ARE FATAL (ADR 71/73 house rule — a mistyped porous flag
// silently changes physics). The pattern's TimeSeries/load-factor is IGNORED by
// design (the overlay owns its own force amplitudes; one notice at first use).

#include "LadrunoPorousOverlay.h"
#include <elementAPI.h>
#include <classTags.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <string>

// OPS_GetStringFromAll portability shim (see OPS_LadrunoUP.cpp): the classic-Tcl
// backend never writes into `buf`, so always use the RETURN value.
static const char* ovGetTok(char* buf, int len)
{
  const char* s = OPS_GetStringFromAll(buf, len);
  if (s == 0) return 0;
  if (s != buf) {
    strncpy(buf, s, (size_t)len - 1);
    buf[len - 1] = '\0';
  }
  return buf;
}

// greedily read integer tokens into `out` until a non-integer token (a flag),
// which is pushed back. Returns the count read.
static int ovReadIntList(std::vector<int>& out)
{
  int n = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    char tok[64];
    ovGetTok(tok, sizeof(tok));
    char* endp = 0;
    long v = strtol(tok, &endp, 10);
    if (endp == tok || *endp != '\0') {   // not a pure integer → next flag
      OPS_ResetCurrentInputArg(-1);
      break;
    }
    out.push_back((int)v);
    n++;
  }
  return n;
}

void* OPS_LadrunoPorousOverlay()
{
  int ndm = OPS_GetNDM();
  if (ndm != 2 && ndm != 3) {
    opserr << "ERROR LadrunoPorousOverlay -- model ndm must be 2 or 3, got "
           << ndm << "\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "ERROR LadrunoPorousOverlay -- missing $tag\n";
    return 0;
  }
  int tag;
  int one = 1;
  if (OPS_GetIntInput(&one, &tag) < 0) {
    opserr << "ERROR LadrunoPorousOverlay -- invalid $tag\n";
    return 0;
  }

  // ---- defaults -------------------------------------------------------------
  std::vector<int> region, drained;
  double Kf = 0.0;        bool KfGiven = false;
  double rhoF = 0.0;      bool rhoFGiven = false;
  double poro = 0.0;      bool poroGiven = false;
  double E = 0.0, nu = 0.0;   bool moduliGiven = false;
  double permG[3] = {0.0, 0.0, 0.0};   bool permGiven = false;
  double biotAlpha = 1.0;
  double Ks = -1.0;                    // <=0 => infinite grains
  double thick = 1.0;
  int    stabMode = LadrunoPorousOverlay::STAB_AUTO;
  double stabValue = 0.25;
  int    fsLMode = LadrunoPorousOverlay::FSL_CLASSIC;
  double fsLScale = 1.0;
  int    staticMode = LadrunoPorousOverlay::SM_MARCH;
  int    onRemoval = LadrunoPorousOverlay::RM_KEEP;
  double onRemovalKf = 1.0;
  int    fluidUpdate = LadrunoPorousOverlay::FU_IMPLICIT;
  int    subcycleMode = LadrunoPorousOverlay::SC_FIXED;
  int    subcycleN = 1;
  int    dynSeepage = 0;
  int    pInitMode = LadrunoPorousOverlay::PI_ZERO;
  double pInitGw = 0.0, pInitZw = 0.0;
  std::vector<int>    pInitNodes;
  std::vector<double> pInitVals;
  double fluidBody[3] = {0.0, 0.0, 0.0};
  std::vector<LadrunoPorousOverlay::Layer> layers;
  std::string recordFile;
  int recordEveryN = 1;

  int num;

  // ---- flag loop (UNKNOWN FLAGS FATAL) --------------------------------------
  while (OPS_GetNumRemainingInputArgs() > 0) {
    char opt[64];
    ovGetTok(opt, sizeof(opt));

    if (strcmp(opt, "-region") == 0) {
      region.clear();
      if (ovReadIntList(region) < 1) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -region needs >=1 element tag\n";
        return 0;
      }

    } else if (strcmp(opt, "-drained") == 0) {
      drained.clear();
      if (ovReadIntList(drained) < 1) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -drained needs >=1 node tag\n";
        return 0;
      }

    } else if (strcmp(opt, "-Kf") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &Kf) < 0 || Kf <= 0.0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag << " -- bad -Kf (>0)\n";
        return 0;
      }
      KfGiven = true;

    } else if (strcmp(opt, "-rhoF") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &rhoF) < 0 || rhoF < 0.0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag << " -- bad -rhoF (>=0)\n";
        return 0;
      }
      rhoFGiven = true;

    } else if (strcmp(opt, "-poro") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &poro) < 0 || poro <= 0.0 || poro > 1.0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- bad -poro (0 < n <= 1)\n";
        return 0;
      }
      poroGiven = true;

    } else if (strcmp(opt, "-perm") == 0) {
      double kk[3];
      num = ndm;
      if (OPS_GetDoubleInput(&num, kk) < 0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag << " -- -perm wants "
               << ndm << " values (k_hydraulic/gammaW per axis)\n";
        return 0;
      }
      for (int i = 0; i < ndm; i++) {
        if (kk[i] < 0.0) {
          opserr << "ERROR LadrunoPorousOverlay " << tag
                 << " -- -perm values must be >= 0\n";
          return 0;
        }
        permG[i] = kk[i];
      }
      permGiven = true;

    } else if (strcmp(opt, "-moduli") == 0) {
      double m[2];
      num = 2;
      if (OPS_GetDoubleInput(&num, m) < 0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -moduli wants $E $nu\n";
        return 0;
      }
      E = m[0]; nu = m[1];
      if (E <= 0.0 || nu <= -1.0 || nu >= 0.5) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- bad -moduli (want E>0, -1<nu<0.5); got E=" << E
               << " nu=" << nu << "\n";
        return 0;
      }
      moduliGiven = true;

    } else if (strcmp(opt, "-alpha") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &biotAlpha) < 0 ||
          biotAlpha <= 0.0 || biotAlpha > 1.0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- bad -alpha (0 < biotAlpha <= 1)\n";
        return 0;
      }

    } else if (strcmp(opt, "-Ks") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &Ks) < 0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- bad -Ks (grain bulk; <=0 = incompressible grains)\n";
        return 0;
      }

    } else if (strcmp(opt, "-thick") == 0 || strcmp(opt, "-thickness") == 0) {
      if (ndm == 3) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -thick is 2D only; remove it in 3D\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &thick) < 0 || thick <= 0.0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag << " -- bad -thick (>0)\n";
        return 0;
      }

    } else if (strcmp(opt, "-layer") == 0) {
      LadrunoPorousOverlay::Layer Ly;
      if (ovReadIntList(Ly.eleTags) < 1) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -layer needs >=1 element tag before its -perm/-poro/-moduli\n";
        return 0;
      }
      bool layerPerm = false;
      // consume the layer-scoped sub-options (perm/poro/moduli) that follow
      while (OPS_GetNumRemainingInputArgs() > 0) {
        char sub[64];
        ovGetTok(sub, sizeof(sub));
        if (strcmp(sub, "-perm") == 0) {
          double kk[3]; num = ndm;
          if (OPS_GetDoubleInput(&num, kk) < 0) {
            opserr << "ERROR LadrunoPorousOverlay " << tag
                   << " -- -layer -perm wants " << ndm << " values\n";
            return 0;
          }
          for (int i = 0; i < ndm; i++) {
            if (kk[i] < 0.0) {
              opserr << "ERROR LadrunoPorousOverlay " << tag
                     << " -- -layer -perm values must be >= 0\n";
              return 0;
            }
            Ly.perm[i] = kk[i];
          }
          layerPerm = true;
        } else if (strcmp(sub, "-poro") == 0) {
          num = 1;
          if (OPS_GetDoubleInput(&num, &Ly.poro) < 0 ||
              Ly.poro <= 0.0 || Ly.poro > 1.0) {
            opserr << "ERROR LadrunoPorousOverlay " << tag
                   << " -- bad -layer -poro (0 < n <= 1)\n";
            return 0;
          }
        } else if (strcmp(sub, "-moduli") == 0) {
          double m[2]; num = 2;
          if (OPS_GetDoubleInput(&num, m) < 0) {
            opserr << "ERROR LadrunoPorousOverlay " << tag
                   << " -- -layer -moduli wants $E $nu\n";
            return 0;
          }
          Ly.E = m[0]; Ly.nu = m[1];
          if (Ly.E <= 0.0 || Ly.nu <= -1.0 || Ly.nu >= 0.5) {
            opserr << "ERROR LadrunoPorousOverlay " << tag
                   << " -- bad -layer -moduli (E>0, -1<nu<0.5)\n";
            return 0;
          }
        } else {
          OPS_ResetCurrentInputArg(-1);   // not a layer sub-option → back to outer
          break;
        }
      }
      if (!layerPerm) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- every -layer must set its own -perm\n";
        return 0;
      }
      layers.push_back(Ly);

    } else if (strcmp(opt, "-stab") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -stab wants 'auto <a0>' | 'off' | a value\n";
        return 0;
      }
      char s[64];
      ovGetTok(s, sizeof(s));
      if (strcmp(s, "off") == 0) {
        stabMode = LadrunoPorousOverlay::STAB_OFF; stabValue = 0.0;
      } else if (strcmp(s, "auto") == 0) {
        stabMode = LadrunoPorousOverlay::STAB_AUTO; stabValue = 0.25;
        if (OPS_GetNumRemainingInputArgs() > 0) {
          char a[64]; ovGetTok(a, sizeof(a));
          char* endp = 0; double a0 = strtod(a, &endp);
          if (endp == a || *endp != '\0') {
            OPS_ResetCurrentInputArg(-1);
          } else if (a0 <= 0.0) {
            opserr << "ERROR LadrunoPorousOverlay " << tag
                   << " -- -stab auto a0 must be > 0\n";
            return 0;
          } else {
            stabValue = a0;
          }
        }
      } else {
        char* endp = 0; double a = strtod(s, &endp);
        if (endp == s || *endp != '\0' || a < 0.0) {
          opserr << "ERROR LadrunoPorousOverlay " << tag
                 << " -- bad -stab '" << s << "' (want auto|off|a value >= 0)\n";
          return 0;
        }
        stabMode = LadrunoPorousOverlay::STAB_MANUAL; stabValue = a;
      }

    } else if (strcmp(opt, "-fsL") == 0) {
      char s[64];
      ovGetTok(s, sizeof(s));
      if (strcmp(s, "classic") == 0) {
        fsLMode = LadrunoPorousOverlay::FSL_CLASSIC;
      } else if (strcmp(s, "oedometric") == 0 || strcmp(s, "oed") == 0) {
        fsLMode = LadrunoPorousOverlay::FSL_OEDOMETRIC;
      } else if (strcmp(s, "zero") == 0) {
        // Ladruno (ADR-73 P3 §3.1): the EXPLICIT-lane setting. Bypasses the
        // oedometric floor (the floor guards the fs1/iterated lanes).
        fsLMode = LadrunoPorousOverlay::FSL_ZERO;
        static bool fsLZeroAdvisoryShown = false;      // one-time, process-wide
        if (!fsLZeroAdvisoryShown) {
          fsLZeroAdvisoryShown = true;
          opserr << "LadrunoPorousOverlay: -fsL zero is the EXPLICIT-lane "
                    "setting (ADR-73 P3): run under an explicit integrator "
                    "(e.g. CentralDifferenceLadruno) at dt <= 0.5x the "
                    "UNDRAINED critical-step pencil (criticalTimeStep is "
                    "overlay-aware). A quasi-static implicit fs1 march with "
                    "L = 0 is the naive drained split and diverges in ~4 "
                    "steps at soil coupling (measured) -- it will fail "
                    "loudly, not wrongly. (printed once)\n";
        }
      } else {
        char* endp = 0; double v = strtod(s, &endp);
        if (endp == s || *endp != '\0' || v <= 0.0) {
          opserr << "ERROR LadrunoPorousOverlay " << tag
                 << " -- bad -fsL '" << s
                 << "' (want classic|oedometric|zero|a positive scale)\n";
          return 0;
        }
        fsLMode = LadrunoPorousOverlay::FSL_MANUAL; fsLScale = v;
      }

    } else if (strcmp(opt, "-staticMode") == 0) {
      char s[64];
      ovGetTok(s, sizeof(s));
      if (strcmp(s, "hold") == 0)        staticMode = LadrunoPorousOverlay::SM_HOLD;
      else if (strcmp(s, "steady") == 0) staticMode = LadrunoPorousOverlay::SM_STEADY;
      else {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- bad -staticMode '" << s << "' (want hold|steady)\n";
        return 0;
      }

    } else if (strcmp(opt, "-onRemoval") == 0) {
      char s[64];
      ovGetTok(s, sizeof(s));
      if (strcmp(s, "keep") == 0) {
        onRemoval = LadrunoPorousOverlay::RM_KEEP;
      } else if (strcmp(s, "drain") == 0) {
        onRemoval = LadrunoPorousOverlay::RM_DRAIN;
        num = 1;
        if (OPS_GetDoubleInput(&num, &onRemovalKf) < 0 || onRemovalKf <= 0.0) {
          opserr << "ERROR LadrunoPorousOverlay " << tag
                 << " -- -onRemoval drain wants a positive $kFactor\n";
          return 0;
        }
      } else {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- bad -onRemoval '" << s << "' (want keep|drain $kF)\n";
        return 0;
      }

    } else if (strcmp(opt, "-fluidUpdate") == 0) {
      char s[64];
      ovGetTok(s, sizeof(s));
      if (strcmp(s, "implicit") == 0) {
        fluidUpdate = LadrunoPorousOverlay::FU_IMPLICIT;
      } else if (strcmp(s, "explicit") == 0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -fluidUpdate explicit is NYI until P3b (ADR §3.4 item 1)\n";
        return 0;
      } else {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- bad -fluidUpdate '" << s << "' (want implicit)\n";
        return 0;
      }

    } else if (strcmp(opt, "-subcycle") == 0) {
      char s[64];
      ovGetTok(s, sizeof(s));
      if (strcmp(s, "auto") == 0) {
        subcycleMode = LadrunoPorousOverlay::SC_AUTO; subcycleN = 1;
      } else {
        char* endp = 0; long v = strtol(s, &endp, 10);
        if (endp == s || *endp != '\0' || v < 1) {
          opserr << "ERROR LadrunoPorousOverlay " << tag
                 << " -- bad -subcycle '" << s << "' (want auto | a positive N)\n";
          return 0;
        }
        subcycleMode = LadrunoPorousOverlay::SC_FIXED; subcycleN = (int)v;
      }

    } else if (strcmp(opt, "-dynSeepage") == 0) {
      char s[64];
      ovGetTok(s, sizeof(s));
      if (strcmp(s, "on") == 0)       dynSeepage = 1;
      else if (strcmp(s, "off") == 0) dynSeepage = 0;
      else {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- bad -dynSeepage '" << s << "' (want on|off)\n";
        return 0;
      }

    } else if (strcmp(opt, "-fluidBody") == 0) {
      double ff[3]; num = ndm;
      if (OPS_GetDoubleInput(&num, ff) < 0) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -fluidBody wants " << ndm << " values\n";
        return 0;
      }
      for (int i = 0; i < ndm; i++) fluidBody[i] = ff[i];

    } else if (strcmp(opt, "-pInit") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -pInit wants steady|hydrostatic|a list\n";
        return 0;
      }
      char s[64];
      ovGetTok(s, sizeof(s));
      if (strcmp(s, "steady") == 0) {
        pInitMode = LadrunoPorousOverlay::PI_STEADY;
      } else if (strcmp(s, "hydrostatic") == 0) {
        pInitMode = LadrunoPorousOverlay::PI_HYDRO;
        num = 1;
        if (OPS_GetDoubleInput(&num, &pInitGw) < 0 || pInitGw <= 0.0) {
          opserr << "ERROR LadrunoPorousOverlay " << tag
                 << " -- -pInit hydrostatic wants $gw (>0) <$zw>\n";
          return 0;
        }
        pInitZw = 0.0;                      // optional phreatic level, default 0
        if (OPS_GetNumRemainingInputArgs() > 0) {
          char z[64]; ovGetTok(z, sizeof(z));
          char* endp = 0; double zv = strtod(z, &endp);
          if (endp == z || *endp != '\0') OPS_ResetCurrentInputArg(-1);
          else pInitZw = zv;
        }
      } else {
        // explicit list: $nd $val [$nd $val ...] — s holds the first node token
        OPS_ResetCurrentInputArg(-1);
        pInitMode = LadrunoPorousOverlay::PI_LIST;
        while (OPS_GetNumRemainingInputArgs() > 0) {
          char nt[64]; ovGetTok(nt, sizeof(nt));
          char* endp = 0; long nd = strtol(nt, &endp, 10);
          if (endp == nt || *endp != '\0') { OPS_ResetCurrentInputArg(-1); break; }
          if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "ERROR LadrunoPorousOverlay " << tag
                   << " -- -pInit list: node " << nd << " has no value\n";
            return 0;
          }
          double val; num = 1;
          if (OPS_GetDoubleInput(&num, &val) < 0) {
            opserr << "ERROR LadrunoPorousOverlay " << tag
                   << " -- -pInit list: bad value for node " << nd << "\n";
            return 0;
          }
          pInitNodes.push_back((int)nd);
          pInitVals.push_back(val);
        }
        if (pInitNodes.empty()) {
          opserr << "ERROR LadrunoPorousOverlay " << tag
                 << " -- -pInit list is empty (want $nd $val ... or "
                    "steady|hydrostatic)\n";
          return 0;
        }
      }

    } else if (strcmp(opt, "-record") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR LadrunoPorousOverlay " << tag
               << " -- -record wants a $file <$everyN>\n";
        return 0;
      }
      char f[512]; ovGetTok(f, sizeof(f));
      recordFile = f;
      if (OPS_GetNumRemainingInputArgs() > 0) {   // optional everyN
        char e[64]; ovGetTok(e, sizeof(e));
        char* endp = 0; long ev = strtol(e, &endp, 10);
        if (endp == e || *endp != '\0') OPS_ResetCurrentInputArg(-1);
        else if (ev >= 1) recordEveryN = (int)ev;
        else {
          opserr << "ERROR LadrunoPorousOverlay " << tag
                 << " -- -record everyN must be >= 1\n";
          return 0;
        }
      }

    } else if (strcmp(opt, "-substep") == 0) {
      opserr << "ERROR LadrunoPorousOverlay " << tag
             << " -- -substep was REMOVED post-E7.3b (the accuracy claim was "
                "refuted; ADR §12 P0 pin 2). Use -subcycle for the explicit "
                "lane\n";
      return 0;

    } else {
      opserr << "ERROR LadrunoPorousOverlay " << tag
             << " -- unknown flag '" << opt
             << "' (unknown flags are fatal for this pattern)\n";
      return 0;
    }
  }

  // ---- required-argument audit ----------------------------------------------
  if (region.empty() || !KfGiven || !rhoFGiven || !poroGiven ||
      !permGiven || !moduliGiven || drained.empty()) {
    opserr << "ERROR LadrunoPorousOverlay " << tag << " -- missing required:";
    if (region.empty())   opserr << " -region";
    if (!KfGiven)         opserr << " -Kf";
    if (!rhoFGiven)       opserr << " -rhoF";
    if (!poroGiven)       opserr << " -poro";
    if (!permGiven)       opserr << " -perm";
    if (!moduliGiven)     opserr << " -moduli";
    if (drained.empty())  opserr << " -drained";
    opserr << "\n";
    return 0;
  }

  // storage-domain police (kernel oneOverQbar contract): n <= biotAlpha.
  if (poro > biotAlpha) {
    opserr << "ERROR LadrunoPorousOverlay " << tag << " -- porosity n=" << poro
           << " exceeds Biot alpha=" << biotAlpha
           << " (storage 1/Qbar needs 0 < n <= alpha <= 1)\n";
    return 0;
  }

  LadrunoPorousOverlay* thePattern = new LadrunoPorousOverlay(
      tag, region, drained, Kf, rhoF, biotAlpha, Ks, poro, E, nu, permG, thick,
      stabMode, stabValue, fsLMode, fsLScale, staticMode, onRemoval, onRemovalKf,
      fluidUpdate, subcycleMode, subcycleN, dynSeepage, pInitMode, pInitGw, pInitZw,
      pInitNodes, pInitVals, fluidBody, layers,
      recordFile.empty() ? 0 : recordFile.c_str(), recordEveryN);

  if (thePattern == 0) {
    opserr << "ERROR LadrunoPorousOverlay " << tag << " -- could not create pattern\n";
    return 0;
  }
  return thePattern;
}
