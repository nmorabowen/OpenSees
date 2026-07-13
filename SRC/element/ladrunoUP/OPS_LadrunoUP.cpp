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
// Parser/factory for LadrunoUP — the unified Biot u-p (solid displacement +
// pore pressure) saturated-porous continuum element (ADR 71, ELE 33017).
//
//   element LadrunoUP $tag $n1 .. $nk $matTag
//       <-thick $t>                          ;# 2D only, default 1.0
//       -Kf $Kf -poro $n -rhoF $rhof         ;# REQUIRED
//       -perm $k1 $k2 <$k3>                  ;# REQUIRED (ndm values), k_hydraulic/gammaW
//       <-permH $k1 $k2 <$k3> -gammaW $gw>   ;# sugar: conductivity / gammaW internally
//       <-alpha $biotAlpha> <-Ks $Ks>        ;# defaults 1.0 / <=0 (infinite grains)
//       <-body $b1 $b2 <$b3>>                ;# accelerations, default 0
//       <-fluidBody $f1 $f2 <$f3>>           ;# defaults to -body
//       <-formulation std|bbar> <-pOrder equal|linear> <-lumped>
//       <-stab auto <$alpha0> | off | $alpha>   ;# equal-order default: auto (a0=0.25)
//       <-dynSeepage on|off>                 ;# default OFF (P4 B5 adjudication)
//       <-geom linear>                       ;# only accepted value (axis reserved)
//
// (ndm, k) selects the shape provider: (2,3) T3 · (2,4) Q4 · (2,6) Bézier T6 ·
// (3,8) H8 · (3,10) Bézier Tet10. The quadratic Bézier shapes (2,6)/(3,10) are
// Taylor–Hood only and REQUIRE -pOrder linear (quadratic u, linear vertex p);
// -pOrder equal (or omitted) on a quadratic shape is fatal — equal-order
// quadratic is a reserved axis with no theory/gate (ADR §3.3, §4.1 legality
// matrix). Linear shapes accept equal (linear is a documented synonym).
//
// UNKNOWN FLAGS ARE FATAL — a deliberate break from the family's
// warn-and-continue idiom (ADR 71 §4.1): a mistyped u-p flag silently changes
// the physics, so the parser refuses instead of ignoring.
//
// Honest-p contract consequence (ADR 71 §3.2): the assembled effective tangent
// is UNSYMMETRIC — a general solver is mandatory (UmfPack / SuperLU /
// FullGeneral / BandGeneral serial; MUMPS SYM=0 in the MPI targets). The
// no-`system`-command default ProfileSPD silently DROPS one coupling block and
// returns plausible wrong pore pressures — hence the one-line solver notice
// printed at element creation (once per process).

#include <LadrunoUP.h>
#include <NDMaterial.h>
#include <elementAPI.h>
#include <ID.h>
#include <Vector.h>
#include <string.h>
#include <stdlib.h>

// OPS_GetStringFromAll portability shim — the two interpreter backends
// disagree on the buffer contract: the modern-interpreter version
// (OpenSeesCommands.cpp:1185) copies the token into `buf` AND returns it,
// but the classic-Tcl version (elementAPI_TCL.cpp:494) just returns
// OPS_GetString() and NEVER touches `buf`. Reading `buf` after the call —
// the family idiom — therefore reads garbage under classic Tcl. Always go
// through the RETURN VALUE, normalized into the caller's buffer here.
static const char *upGetTok(char *buf, int len)
{
  const char *s = OPS_GetStringFromAll(buf, len);
  if (s == 0)
    return 0;
  if (s != buf) {
    strncpy(buf, s, (size_t)len - 1);
    buf[len - 1] = '\0';
  }
  return buf;
}

void *OPS_LadrunoUP()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  // ndm gate only, here. The model-builder NDF check is DEFERRED until the
  // node count is known (below): the Taylor–Hood modeling pattern necessarily
  // toggles the builder between ndf=ndm+1 (vertex carriers) and ndf=ndm
  // (mid-edge nodes), so at element-command time the builder can legally sit
  // in either state — setDomain does the strict per-node validation. Linear
  // (equal-order) shapes keep the early builder gate for a clear P1-style
  // error. (3.B integration fix, WP3.A pin 2 rationale.)
  if (ndm != 2 && ndm != 3) {
    opserr << "ERROR LadrunoUP -- model ndm must be 2 or 3, got " << ndm << "\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "ERROR LadrunoUP -- insufficient arguments\n"
              "Want: element LadrunoUP tag n1 .. nk matTag -Kf Kf -poro n "
              "-rhoF rhof -perm k1 k2 <k3> <options>\n";
    return 0;
  }

  // ---- leading integers: tag, n1..nk, matTag (k = 3|4|8 at P1) -------------
  // Greedy peek via OPS_GetStringFromAll (NOT OPS_GetString — openseespy
  // numeric args do not survive OPS_GetString) + strtol, pushing the first
  // non-integer token back with OPS_ResetCurrentInputArg(-1) (family idiom,
  // cf. OPS_LadrunoKinematicCoupling.cpp).
  const int maxLead = 13;   // tag + 10 nodes + matTag + 1 slack for the count check
  int lead[maxLead];
  int nLead = 0;
  while (OPS_GetNumRemainingInputArgs() > 0 && nLead < maxLead) {
    char tok[64];
    upGetTok(tok, sizeof(tok));
    char *endp = 0;
    long v = strtol(tok, &endp, 10);
    if (endp == tok || *endp != '\0') {     // not an integer — first flag
      OPS_ResetCurrentInputArg(-1);
      break;
    }
    lead[nLead++] = (int)v;
  }

  if (nLead < 5) {
    opserr << "ERROR LadrunoUP -- need tag + >=3 node tags + matTag before the "
              "flags (got " << nLead << " leading integers)\n";
    return 0;
  }

  int tag = lead[0];
  int matTag = lead[nLead - 1];
  int nNodes = nLead - 2;

  // ---- (ndm, nNodes) provider legality (ADR 71 §4.1) -----------------------
  // P3 opens the quadratic Bézier shapes (Taylor–Hood). isQuadratic gates the
  // -pOrder handling below (they demand -pOrder linear; equal-order quadratic
  // is a reserved axis).
  bool isQuadratic = (ndm == 2 && nNodes == 6) || (ndm == 3 && nNodes == 10);
  if (!((ndm == 2 && (nNodes == 3 || nNodes == 4 || nNodes == 6)) ||
        (ndm == 3 && (nNodes == 8 || nNodes == 10)))) {
    opserr << "ERROR LadrunoUP " << tag << " -- no shape provider for ndm="
           << ndm << " with " << nNodes
           << " nodes (want (2,3) T3, (2,4) Q4, (2,6) Bézier T6, (3,8) H8, "
              "or (3,10) Bézier Tet10)\n";
    return 0;
  }

  // Deferred builder-NDF gate (see the ndm-gate comment above): equal-order
  // linear shapes demand ndf = ndm+1 at the builder; Taylor–Hood quadratic
  // shapes skip it (mixed-ndf pattern) — setDomain validates per node.
  if (!isQuadratic && ndf != ndm + 1) {
    opserr << "ERROR LadrunoUP " << tag << " -- model ndm/ndf must be 2/3 "
              "(2D u-p) or 3/4 (3D u-p) for the equal-order shapes; got "
           << ndm << "/" << ndf << ". The pressure DOF needs ndf = ndm+1.\n";
    return 0;
  }

  ID nodeTags(nNodes);
  for (int i = 0; i < nNodes; i++)
    nodeTags(i) = lead[1 + i];

  NDMaterial *mat = OPS_getNDMaterial(matTag);
  if (mat == 0) {
    opserr << "ERROR LadrunoUP " << tag << " -- NDMaterial " << matTag
           << " not found\n";
    return 0;
  }

  // ---- defaults -------------------------------------------------------------
  double thick = 1.0;      bool thickGiven = false;
  double Kf = 0.0;         bool KfGiven = false;
  double poro = 0.0;       bool poroGiven = false;
  double rhoF = 0.0;       bool rhoFGiven = false;
  Vector perm(ndm);        bool permGiven = false;   // k_hydraulic/gammaW
  Vector permH(ndm);       bool permHGiven = false;  // hydraulic conductivity
  double gammaW = 0.0;     bool gammaWGiven = false;
  double biotAlpha = 1.0;
  double Ks = -1.0;                                  // <=0 => infinite grains
  Vector body(ndm);        // zero-initialized by Vector
  Vector fluidBody(ndm);   bool fluidBodyGiven = false;
  int formulation = 0;     // 0 = std, 1 = bbar
  int pOrder = 0;          // 0 = equal-order; 1 = Taylor–Hood (linear vertex p)
  bool pOrderGiven = false;
  bool lumped = false;
  int stabMode = 1;        // 0 = off, 1 = auto, 2 = manual — default auto
  double stabValue = 0.25; // auto: alpha0; manual: alpha
  bool stabGiven = false;
  bool dynSeepage = false;  // P4 B5 adjudication: default OFF. The u-p
                            // dynamic-seepage drive (b - u_ddot) diverges under
                            // quasi-static dt-refinement (P1, ZS84 sweep) AND
                            // misbehaves in genuine dynamics (P4 B5: wandering
                            // pressure plateau, unbounded shallow-station
                            // growth) - trial-accel noise feeds f_seep. Opt-in
                            // '-dynSeepage on' retained (research/SWANDYNE
                            // parity; the +G residual term stays FD-gated).

  int num;

  // ---- flag loop (UNKNOWN FLAGS FATAL) --------------------------------------
  while (OPS_GetNumRemainingInputArgs() > 0) {
    char opt[64];
    upGetTok(opt, sizeof(opt));

    if (strcmp(opt, "-thick") == 0 || strcmp(opt, "-thickness") == 0) {
      if (ndm == 3) {
        opserr << "ERROR LadrunoUP " << tag
               << " -- -thick is 2D only (plane strain); remove it in 3D\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &thick) < 0 || thick <= 0.0) {
        opserr << "ERROR LadrunoUP " << tag << " -- bad -thick (want t > 0)\n";
        return 0;
      }
      thickGiven = true;

    } else if (strcmp(opt, "-Kf") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &Kf) < 0 || Kf <= 0.0) {
        opserr << "ERROR LadrunoUP " << tag << " -- bad -Kf (want Kf > 0)\n";
        return 0;
      }
      KfGiven = true;

    } else if (strcmp(opt, "-poro") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &poro) < 0 || poro <= 0.0 || poro > 1.0) {
        opserr << "ERROR LadrunoUP " << tag
               << " -- bad -poro (want 0 < n <= 1)\n";
        return 0;
      }
      poroGiven = true;

    } else if (strcmp(opt, "-rhoF") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &rhoF) < 0 || rhoF < 0.0) {
        opserr << "ERROR LadrunoUP " << tag << " -- bad -rhoF (want rhof >= 0)\n";
        return 0;
      }
      rhoFGiven = true;

    } else if (strcmp(opt, "-perm") == 0) {
      if (permHGiven) {
        opserr << "ERROR LadrunoUP " << tag
               << " -- give -perm OR -permH/-gammaW, not both\n";
        return 0;
      }
      double kk[3];
      num = ndm;
      if (OPS_GetDoubleInput(&num, kk) < 0) {
        opserr << "ERROR LadrunoUP " << tag << " -- -perm wants " << ndm
               << " values (k_hydraulic/gammaW per axis)\n";
        return 0;
      }
      for (int i = 0; i < ndm; i++) {
        if (kk[i] < 0.0) {
          opserr << "ERROR LadrunoUP " << tag
                 << " -- -perm values must be >= 0\n";
          return 0;
        }
        perm(i) = kk[i];
      }
      permGiven = true;

    } else if (strcmp(opt, "-permH") == 0) {
      if (permGiven) {
        opserr << "ERROR LadrunoUP " << tag
               << " -- give -perm OR -permH/-gammaW, not both\n";
        return 0;
      }
      double kk[3];
      num = ndm;
      if (OPS_GetDoubleInput(&num, kk) < 0) {
        opserr << "ERROR LadrunoUP " << tag << " -- -permH wants " << ndm
               << " values (hydraulic conductivity per axis)\n";
        return 0;
      }
      for (int i = 0; i < ndm; i++) {
        if (kk[i] < 0.0) {
          opserr << "ERROR LadrunoUP " << tag
                 << " -- -permH values must be >= 0\n";
          return 0;
        }
        permH(i) = kk[i];
      }
      permHGiven = true;

    } else if (strcmp(opt, "-gammaW") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &gammaW) < 0 || gammaW <= 0.0) {
        opserr << "ERROR LadrunoUP " << tag
               << " -- bad -gammaW (want gw > 0)\n";
        return 0;
      }
      gammaWGiven = true;

    } else if (strcmp(opt, "-alpha") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &biotAlpha) < 0 ||
          biotAlpha <= 0.0 || biotAlpha > 1.0) {
        opserr << "ERROR LadrunoUP " << tag
               << " -- bad -alpha (want 0 < biotAlpha <= 1)\n";
        return 0;
      }

    } else if (strcmp(opt, "-Ks") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &Ks) < 0) {
        opserr << "ERROR LadrunoUP " << tag
               << " -- bad -Ks (grain bulk modulus; <= 0 = incompressible "
                  "grains)\n";
        return 0;
      }

    } else if (strcmp(opt, "-body") == 0) {
      double bb[3];
      num = ndm;
      if (OPS_GetDoubleInput(&num, bb) < 0) {
        opserr << "ERROR LadrunoUP " << tag << " -- -body wants " << ndm
               << " values (accelerations)\n";
        return 0;
      }
      for (int i = 0; i < ndm; i++)
        body(i) = bb[i];

    } else if (strcmp(opt, "-fluidBody") == 0) {
      double ff[3];
      num = ndm;
      if (OPS_GetDoubleInput(&num, ff) < 0) {
        opserr << "ERROR LadrunoUP " << tag << " -- -fluidBody wants " << ndm
               << " values (accelerations)\n";
        return 0;
      }
      for (int i = 0; i < ndm; i++)
        fluidBody(i) = ff[i];
      fluidBodyGiven = true;

    } else if (strcmp(opt, "-formulation") == 0 || strcmp(opt, "-form") == 0) {
      char f[32];
      upGetTok(f, sizeof(f));
      if (strcmp(f, "std") == 0)       formulation = 0;
      else if (strcmp(f, "bbar") == 0) formulation = 1;
      else {
        opserr << "ERROR LadrunoUP " << tag << " -- unknown -formulation '" << f
               << "' (want std|bbar; ssp/eas/uri are not u-p axes, ADR 71 "
                  "§3.4)\n";
        return 0;
      }

    } else if (strcmp(opt, "-pOrder") == 0) {
      char p[32];
      upGetTok(p, sizeof(p));
      if (strcmp(p, "equal") == 0) {
        pOrder = 0;
      } else if (strcmp(p, "linear") == 0) {
        // On a linear shape every node is a vertex — 'linear' is a documented
        // synonym of 'equal' (identical interpolation). On a quadratic Bézier
        // shape 'linear' is the Taylor–Hood cure: quadratic u, linear vertex p.
        pOrder = isQuadratic ? 1 : 0;
      } else {
        opserr << "ERROR LadrunoUP " << tag << " -- unknown -pOrder '" << p
               << "' (want equal|linear)\n";
        return 0;
      }
      pOrderGiven = true;

    } else if (strcmp(opt, "-lumped") == 0) {
      lumped = true;

    } else if (strcmp(opt, "-stab") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR LadrunoUP " << tag
               << " -- -stab wants 'auto <alpha0>' | 'off' | a value\n";
        return 0;
      }
      char sTok[64];
      upGetTok(sTok, sizeof(sTok));
      if (strcmp(sTok, "off") == 0) {
        stabMode = 0;
        stabValue = 0.0;
      } else if (strcmp(sTok, "auto") == 0) {
        stabMode = 1;
        stabValue = 0.25;                      // papers' pinned default alpha0
        if (OPS_GetNumRemainingInputArgs() > 0) {
          char aTok[64];                       // optional alpha0 — peek
          upGetTok(aTok, sizeof(aTok));
          char *endp = 0;
          double a0 = strtod(aTok, &endp);
          if (endp == aTok || *endp != '\0') {
            OPS_ResetCurrentInputArg(-1);      // not a number — next flag
          } else {
            if (a0 <= 0.0) {
              opserr << "ERROR LadrunoUP " << tag
                     << " -- -stab auto alpha0 must be > 0 (got " << a0
                     << "; a non-positive alpha0 DEstabilizes)\n";
              return 0;
            }
            if (a0 < 0.1 || a0 > 0.5) {
              // once per process, like the solver/stab notices — per-element
              // repetition spams large meshes (1.F-ii finding)
              static bool a0RangeWarned = false;
              if (!a0RangeWarned) {
                a0RangeWarned = true;
                opserr << "WARNING LadrunoUP " << tag << " -- -stab auto alpha0="
                       << a0 << " is outside the papers' supported range "
                          "[0.1, 0.5] (McGann 2012/2015) (printed once)\n";
              }
            }
            stabValue = a0;
          }
        }
      } else {
        char *endp = 0;
        double a = strtod(sTok, &endp);
        if (endp == sTok || *endp != '\0') {
          opserr << "ERROR LadrunoUP " << tag << " -- bad -stab value '" << sTok
                 << "' (want 'auto <alpha0>' | 'off' | a numeric alpha)\n";
          return 0;
        }
        if (a < 0.0) {
          opserr << "ERROR LadrunoUP " << tag
                 << " -- manual -stab alpha must be >= 0 (a negative alpha "
                    "anti-stabilizes the pressure block)\n";
          return 0;
        }
        stabMode = 2;
        stabValue = a;
      }
      stabGiven = true;

    } else if (strcmp(opt, "-dynSeepage") == 0) {
      char d[32];
      upGetTok(d, sizeof(d));
      if (strcmp(d, "on") == 0)       dynSeepage = true;
      else if (strcmp(d, "off") == 0) dynSeepage = false;
      else {
        opserr << "ERROR LadrunoUP " << tag << " -- bad -dynSeepage '" << d
               << "' (want on|off)\n";
        return 0;
      }

    } else if (strcmp(opt, "-geom") == 0) {
      char g[32];
      upGetTok(g, sizeof(g));
      if (strcmp(g, "linear") != 0) {
        opserr << "ERROR LadrunoUP " << tag << " -- -geom '" << g
               << "' not supported: only 'linear' is accepted (the axis is "
                  "reserved; finite-strain u-p is a research phase, ADR 71 "
                  "§2.4)\n";
        return 0;
      }

    } else {
      // DELIBERATE break from the family's warn-and-continue: a mistyped u-p
      // flag silently changes physics (ADR 71 §4.1).
      opserr << "ERROR LadrunoUP " << tag << " -- unknown flag '" << opt
             << "' (unknown flags are fatal for this element)\n";
      return 0;
    }
  }

  // ---- required-argument audit ----------------------------------------------
  if (!KfGiven || !poroGiven || !rhoFGiven || !(permGiven || permHGiven)) {
    opserr << "ERROR LadrunoUP " << tag << " -- missing required argument(s):";
    if (!KfGiven)   opserr << " -Kf";
    if (!poroGiven) opserr << " -poro";
    if (!rhoFGiven) opserr << " -rhoF";
    if (!(permGiven || permHGiven)) opserr << " -perm (or -permH/-gammaW)";
    opserr << "\n";
    return 0;
  }

  // -permH/-gammaW sugar: both or neither.
  if (permHGiven != gammaWGiven) {
    opserr << "ERROR LadrunoUP " << tag
           << " -- -permH and -gammaW come together (conductivity / gammaW "
              "= the -perm value)\n";
    return 0;
  }
  if (permHGiven) {
    for (int i = 0; i < ndm; i++)
      perm(i) = permH(i) / gammaW;
  }

  // Storage-domain police (kernel oneOverQbar contract): n <= biotAlpha.
  if (poro > biotAlpha) {
    opserr << "ERROR LadrunoUP " << tag << " -- porosity n=" << poro
           << " exceeds Biot alpha=" << biotAlpha
           << " (storage 1/Qbar needs 0 < n <= alpha <= 1)\n";
    return 0;
  }

  // Quadratic Bézier shapes are Taylor–Hood ONLY: they REQUIRE -pOrder linear.
  // Omitted or -pOrder equal on a quadratic shape is fatal — equal-order
  // quadratic is a reserved axis with no theory or gate (ADR §3.3, §4.1).
  if (isQuadratic && pOrder != 1) {
    opserr << "ERROR LadrunoUP " << tag << " -- the " << nNodes
           << "-node Bézier shape is Taylor–Hood only and requires '-pOrder "
              "linear' (quadratic u, linear vertex p)"
           << (pOrderGiven ? "; '-pOrder equal'"
                           : "; -pOrder was omitted, which")
           << " selects equal-order quadratic — a reserved axis with no theory "
              "or gate\n";
    return 0;
  }

  // -stab on a Taylor–Hood pair is fatal (TH is inf-sup stable without it).
  if (pOrder == 1 && stabGiven) {
    opserr << "ERROR LadrunoUP " << tag
           << " -- -stab is not accepted on Taylor–Hood pairs (they are "
              "inf-sup stable without it)\n";
    return 0;
  }

  // -fluidBody defaults to -body when omitted.
  if (!fluidBodyGiven)
    fluidBody = body;

  (void)thickGiven;   // 2D default 1.0 is fine silently; 3D always 1.0

  // Equal-order pair with -stab omitted => AUTO (alpha0 = 0.25), one-line
  // notice (once per process — B4's '-stab off' leg is always an explicit
  // opt-out, ADR 71 §4.1).
  if (!stabGiven && pOrder == 0) {
    static bool stabNoticeShown = false;
    if (!stabNoticeShown) {
      stabNoticeShown = true;
      opserr << "LadrunoUP: equal-order u-p pair — '-stab' omitted, defaulting "
                "to 'auto' (alpha0=0.25, McGann); use '-stab off' to opt out\n";
    }
  }

  // Mandatory one-line solver notice (once per process): the honest-p tangent
  // is unsymmetric, and the no-system default ProfileSPD silently drops a
  // coupling block (ADR 71 §3.2 ⟨FW-F1⟩).
  {
    static bool solverNoticeShown = false;
    if (!solverNoticeShown) {
      solverNoticeShown = true;
      opserr << "LadrunoUP: unsymmetric tangent — use a general solver "
                "(UmfPack/SuperLU/FullGeneral/BandGeneral; MUMPS SYM=0 in "
                "MPI); the ProfileSPD default silently drops a coupling "
                "block\n";
    }
  }

  // ---- build (pinned WP1.A ctor, verbatim) -----------------------------------
  return new LadrunoUP(tag, ndm, nodeTags, *mat,
                       thick, Kf, poro, rhoF, perm,
                       biotAlpha, Ks, body, fluidBody,
                       formulation, pOrder, lumped,
                       stabMode, stabValue, dynSeepage);
}
