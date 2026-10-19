/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
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
                                                                        
// $Revision$
// $Date$
// $URL$

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for LadrunoDispBeamColumn3d.

#include <LadrunoDispBeamColumn3d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>   // Ladruno (ADR 33) Tier-2: embedded cohesive hinge
#include <LadrunoCohesiveHinge.h>
#include <Bilin.h>
#include <ModIMKPeakOriented.h>
#include <NDMaterial.h>         // Ladruno (ADR 34): coupled biaxial cohesive hinge (order-2)
#include <LadrunoCohesiveHingeBiaxial.h>
#include <CrdTransf.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <stdlib.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>
#include <Parameter.h>
#include <math.h>
#include <cmath>
#include <elementAPI.h>
#include <string>
#include <LadrunoResponseTokens.h>
#include <LadrunoThreads.h>

// Ladruno (ADR 32): active element for crack-band materials (see 2D sibling).
extern thread_local Element *ops_TheActiveElement;   // Ladruno WP-107: thread_local, see Element.cpp

static thread_local Matrix K(12,12);
static thread_local Vector P(12);
static thread_local double workArea[200];

static bool ldbc3dGetHingeImplexError(UniaxialMaterial* mat, double& errOut)
{
  if (mat == 0) return false;
  if (LadrunoCohesiveHinge* ch = dynamic_cast<LadrunoCohesiveHinge*>(mat)) {
    errOut = ch->getImplexError();
    return true;
  }
  if (Bilin* b = dynamic_cast<Bilin*>(mat)) {
    errOut = b->getImplexError();
    return true;
  }
  if (ModIMKPeakOriented* m = dynamic_cast<ModIMKPeakOriented*>(mat)) {
    errOut = m->getImplexError();
    return true;
  }
  return false;
}

static bool ldbc3dHingeIsImplex(UniaxialMaterial* mat)
{
  if (mat == 0) return false;
  if (LadrunoCohesiveHinge* ch = dynamic_cast<LadrunoCohesiveHinge*>(mat))
    return ch->isImplex();
  if (Bilin* b = dynamic_cast<Bilin*>(mat))
    return b->isImplex();
  if (ModIMKPeakOriented* m = dynamic_cast<ModIMKPeakOriented*>(mat))
    return m->isImplex();
  return false;
}

static bool ldbc3dHingeBiaxialIsImplex(NDMaterial* mat)
{
  if (mat == 0) return false;
  if (LadrunoCohesiveHingeBiaxial* ch = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(mat))
    return ch->isImplex();
  return false;
}

static bool
ladrunoPeekNextInt3d(int &outTag)
{
    if (OPS_GetNumRemainingInputArgs() < 1)
        return false;
    char tok[128];
    const char *s = OPS_GetStringFromAll(tok, sizeof(tok));
    OPS_ResetCurrentInputArg(-1);
    if (s == 0 || *s == '\0')
        return false;
    char *endp = 0;
    long v = strtol(s, &endp, 10);
    if (endp == s || *endp != '\0')
        return false;
    int one = 1;
    if (OPS_GetIntInput(&one, &outTag) < 0) {
        outTag = (int)v;
    }
    return true;
}

static bool
ladrunoPeekNextDouble3d(double &outVal)
{
    if (OPS_GetNumRemainingInputArgs() < 1)
        return false;
    char tok[128];
    const char *s = OPS_GetStringFromAll(tok, sizeof(tok));
    OPS_ResetCurrentInputArg(-1);
    if (s == 0 || *s == '\0')
        return false;
    char *endp = 0;
    double v = strtod(s, &endp);
    if (endp == s || *endp != '\0')
        return false;
    int one = 1;
    if (OPS_GetDoubleInput(&one, &outVal) < 0) {
        outVal = v;
    }
    return true;
}

static inline double ladrunoDisp3dScalePM(double P_axial, double P0, double alphaPM)
{
  if (P0 <= 0.0)
    return 1.0;
  double n = P_axial / P0;
  double f = 1.0 - alphaPM * n - (1.0 - alphaPM) * n * n;
  if (f < 0.05) f = 0.05;
  if (f > 2.5)  f = 2.5;
  return f;
}

void* OPS_LadrunoDispBeamColumn3d()
{
    if (OPS_GetNDM() != 3 || OPS_GetNDF() != 6) {
        opserr << "WARNING LadrunoDispBeamColumn3d -- model must be ndm 3/ndf 6\n";
        return 0;
    }
    int dampingTag = 0;
    Damping* theDamping = 0;
    if(OPS_GetNumRemainingInputArgs() < 5) {
	opserr<<"insufficient arguments:eleTag,iNode,jNode,transfTag,integrationTag <-mass mass> <-cmass>\n";
	return 0;
    }

    // inputs: 
    int iData[5];
    int numData = 5;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    // options
    double mass = 0.0;
    int cmass = 0;
    int lchMode = 0;       // Ladruno (ADR 32): 0=ip (default), 1=element, 2=user
    double userLch = 0.0;
    int nlGeom = 0;        // Ladruno (ADR 32): 0=linear basic strain, 1=½(θy²+θz²) bowing (-nl)
    UniaxialMaterial *hingeMatZ = 0;   // Ladruno (ADR 33) Tier-2: embedded strong-axis (Mz) cohesive hinge (end I / legacy)
    UniaxialMaterial *hingeMatY = 0;   // Ladruno (ADR 33 PR-3b): embedded weak-axis (My) cohesive hinge (end I)
    NDMaterial *hingeMatC = 0;         // Ladruno (ADR 34): coupled biaxial (Mz-My) cohesive hinge (end I)
    UniaxialMaterial *hingeMatZJ = 0;  // Ladruno (WP-200): end-J strong-axis (Mz) cohesive hinge
    UniaxialMaterial *hingeMatYJ = 0;  // Ladruno (WP-200): end-J weak-axis (My) cohesive hinge
    NDMaterial *hingeMatCJ = 0;        // Ladruno (WP-200): end-J coupled biaxial (Mz-My) cohesive hinge
    int dualHingeMode = 0;
    double xiI = 0.0;
    double xiJ = 1.0;
    bool eleImplex3d = false;
    double elePmP0_3d = 0.0;
    double elePmAlpha_3d = 0.5;
    numData = 1;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type,"-cMass") == 0 || strcmp(type,"-cmass") == 0) {
	    cmass = 1;
	} else if(strcmp(type,"-mass") == 0) {
	    if(OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"WARNING LadrunoDispBeamColumn3d: -mass needs a value\n";
		return 0;
	    }
	    if(OPS_GetDoubleInput(&numData,&mass) < 0) {
		opserr<<"WARNING: invalid mass\n";
		return 0;
	    }
	}
    else if (strcmp(type, "-damp") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -damp needs a damping tag\n";
            return 0;
        }
        if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
        theDamping = OPS_getDamping(dampingTag);
        if (theDamping == 0) {
            opserr << "damping not found\n";
            return 0;
        }
    }
    // Ladruno (ADR 32): -lch {ip|element|<value>} regularization length
    else if (strcmp(type, "-lch") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -lch needs an argument {ip|element|<value>}\n";
            return 0;
        }
        const char* lchArg = OPS_GetString();
        if (strcmp(lchArg, "ip") == 0) {
            lchMode = 0;
        } else if (strcmp(lchArg, "element") == 0) {
            lchMode = 1;
            opserr << "WARNING LadrunoDispBeamColumn3d: -lch element regularizes softening over "
                   << "the WHOLE element length (multi-IP energy double-counting); A/B debugging only.\n";
        } else {
            userLch = atof(lchArg);
            if (!std::isfinite(userLch) || userLch <= 0.0) {
                opserr << "WARNING LadrunoDispBeamColumn3d: -lch value must be ip|element|<finite positive value>, got '"
                       << lchArg << "'\n";
                return 0;
            }
            lchMode = 2;
        }
    }
    // Ladruno (ADR 32): -nl enables the ½(θy²+θz²) biaxial bowing (nonlinear) basic strain.
    else if (strcmp(type, "-nl") == 0) {
        nlGeom = 1;
    }
    // Ladruno (ADR 33) Tier-2 PR-3a: -hinge $matTag embeds the strong-axis (Mz) cohesive
    // rotation-jump hinge carrying any UniaxialMaterial (LadrunoCohesiveHinge default).
    else if (strcmp(type, "-hinge") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hinge needs a uniaxialMaterial tag\n";
            return 0;
        }
        int hingeTag;
        if (OPS_GetIntInput(&numData, &hingeTag) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid -hinge material tag\n";
            return 0;
        }
        hingeMatZ = OPS_getUniaxialMaterial(hingeTag);
        if (hingeMatZ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hinge uniaxialMaterial " << hingeTag
                   << " not found\n";
            return 0;
        }
    }
    else if (strcmp(type, "-hingeI") == 0 || strcmp(type, "-hingeZI") == 0 || strcmp(type, "-hinge1") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " needs a uniaxialMaterial tag\n";
            return 0;
        }
        int hingeTag;
        if (OPS_GetIntInput(&numData, &hingeTag) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid " << type << " material tag\n";
            return 0;
        }
        hingeMatZ = OPS_getUniaxialMaterial(hingeTag);
        if (hingeMatZ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " uniaxialMaterial " << hingeTag
                   << " not found\n";
            return 0;
        }
        dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeJ") == 0 || strcmp(type, "-hingeZJ") == 0 || strcmp(type, "-hinge2") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " needs a uniaxialMaterial tag\n";
            return 0;
        }
        int hingeTagJ;
        if (OPS_GetIntInput(&numData, &hingeTagJ) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid " << type << " material tag\n";
            return 0;
        }
        hingeMatZJ = OPS_getUniaxialMaterial(hingeTagJ);
        if (hingeMatZJ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " uniaxialMaterial " << hingeTagJ
                   << " not found\n";
            return 0;
        }
        dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeIJ") == 0 || strcmp(type, "-hingeZIJ") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " needs at least one uniaxialMaterial tag\n";
            return 0;
        }
        int tagI = 0;
        if (OPS_GetIntInput(&numData, &tagI) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid " << type << " material tag\n";
            return 0;
        }
        int tagJ = tagI;
        int optJ = 0;
        if (ladrunoPeekNextInt3d(optJ))
            tagJ = optJ;
        hingeMatZ  = OPS_getUniaxialMaterial(tagI);
        hingeMatZJ = OPS_getUniaxialMaterial(tagJ);
        if (hingeMatZ == 0 || hingeMatZJ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " uniaxialMaterial not found\n";
            return 0;
        }
        dualHingeMode = 1;
    }
    // Ladruno (ADR 33) Tier-2 PR-3b: -hingeY $matTag adds the weak-axis (My) cohesive
    // rotation-jump, making the hinge biaxial (true coupled 2x2 condensation). Requires -hinge.
    else if (strcmp(type, "-hingeY") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hingeY needs a uniaxialMaterial tag\n";
            return 0;
        }
        int hingeTagY;
        if (OPS_GetIntInput(&numData, &hingeTagY) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid -hingeY material tag\n";
            return 0;
        }
        hingeMatY = OPS_getUniaxialMaterial(hingeTagY);
        if (hingeMatY == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hingeY uniaxialMaterial " << hingeTagY
                   << " not found\n";
            return 0;
        }
    }
    else if (strcmp(type, "-hingeYI") == 0 || strcmp(type, "-hingeY1") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " needs a uniaxialMaterial tag\n";
            return 0;
        }
        int hingeTagY;
        if (OPS_GetIntInput(&numData, &hingeTagY) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid " << type << " material tag\n";
            return 0;
        }
        hingeMatY = OPS_getUniaxialMaterial(hingeTagY);
        if (hingeMatY == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " uniaxialMaterial " << hingeTagY
                   << " not found\n";
            return 0;
        }
        dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeYJ") == 0 || strcmp(type, "-hingeY2") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " needs a uniaxialMaterial tag\n";
            return 0;
        }
        int hingeTagYJ;
        if (OPS_GetIntInput(&numData, &hingeTagYJ) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid " << type << " material tag\n";
            return 0;
        }
        hingeMatYJ = OPS_getUniaxialMaterial(hingeTagYJ);
        if (hingeMatYJ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " uniaxialMaterial " << hingeTagYJ
                   << " not found\n";
            return 0;
        }
        dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeYIJ") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hingeYIJ needs at least one uniaxialMaterial tag\n";
            return 0;
        }
        int tagYI = 0;
        if (OPS_GetIntInput(&numData, &tagYI) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid -hingeYIJ material tag\n";
            return 0;
        }
        int tagYJ = tagYI;
        int optYJ = 0;
        if (ladrunoPeekNextInt3d(optYJ))
            tagYJ = optYJ;
        hingeMatY  = OPS_getUniaxialMaterial(tagYI);
        hingeMatYJ = OPS_getUniaxialMaterial(tagYJ);
        if (hingeMatY == 0 || hingeMatYJ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hingeYIJ uniaxialMaterial not found\n";
            return 0;
        }
        dualHingeMode = 1;
    }
    // Ladruno (ADR 34): -hingeBiaxial $matTag embeds the COUPLED biaxial cohesive law (an order-2
    // NDMaterial, LadrunoCohesiveHingeBiaxial) carrying the Mz-My interaction surface. Exclusive
    // with the block-diagonal -hinge/-hingeY (and -nl).
    else if (strcmp(type, "-hingeBiaxial") == 0 || strcmp(type, "-hingeC") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hingeBiaxial needs an nDMaterial tag\n";
            return 0;
        }
        int hingeTagC;
        if (OPS_GetIntInput(&numData, &hingeTagC) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid -hingeBiaxial material tag\n";
            return 0;
        }
        hingeMatC = OPS_getNDMaterial(hingeTagC);
        if (hingeMatC == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hingeBiaxial nDMaterial " << hingeTagC
                   << " not found\n";
            return 0;
        }
    }
    else if (strcmp(type, "-hingeBiaxialI") == 0 || strcmp(type, "-hingeBiaxial1") == 0 || strcmp(type, "-hingeCI") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " needs an nDMaterial tag\n";
            return 0;
        }
        int hingeTagC;
        if (OPS_GetIntInput(&numData, &hingeTagC) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid " << type << " material tag\n";
            return 0;
        }
        hingeMatC = OPS_getNDMaterial(hingeTagC);
        if (hingeMatC == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " nDMaterial " << hingeTagC
                   << " not found\n";
            return 0;
        }
        dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeBiaxialJ") == 0 || strcmp(type, "-hingeBiaxial2") == 0 || strcmp(type, "-hingeCJ") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " needs an nDMaterial tag\n";
            return 0;
        }
        int hingeTagCJ;
        if (OPS_GetIntInput(&numData, &hingeTagCJ) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid " << type << " material tag\n";
            return 0;
        }
        hingeMatCJ = OPS_getNDMaterial(hingeTagCJ);
        if (hingeMatCJ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " nDMaterial " << hingeTagCJ
                   << " not found\n";
            return 0;
        }
        dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeBiaxialIJ") == 0 || strcmp(type, "-hingeCIJ") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " needs at least one nDMaterial tag\n";
            return 0;
        }
        int tagCI = 0;
        if (OPS_GetIntInput(&numData, &tagCI) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: invalid " << type << " material tag\n";
            return 0;
        }
        int tagCJ = tagCI;
        int optCJ = 0;
        if (ladrunoPeekNextInt3d(optCJ))
            tagCJ = optCJ;
        hingeMatC  = OPS_getNDMaterial(tagCI);
        hingeMatCJ = OPS_getNDMaterial(tagCJ);
        if (hingeMatC == 0 || hingeMatCJ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: " << type << " nDMaterial not found\n";
            return 0;
        }
        dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeLoc") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hingeLoc needs xiI [xiJ]\n";
            return 0;
        }
        if (OPS_GetDoubleInput(&numData, &xiI) < 0 || !std::isfinite(xiI) || xiI < 0.0 || xiI > 1.0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -hingeLoc xiI must be in [0, 1]\n";
            return 0;
        }
        double optXiJ = 0.0;
        if (ladrunoPeekNextDouble3d(optXiJ)) {
            if (!std::isfinite(optXiJ) || optXiJ < 0.0 || optXiJ > 1.0) {
                opserr << "WARNING LadrunoDispBeamColumn3d: -hingeLoc xiJ must be in [0, 1]\n";
                return 0;
            }
            xiJ = optXiJ;
        } else {
            xiJ = 1.0 - xiI;
        }
        dualHingeMode = 1;
    }
    else if (strcmp(type, "-implex") == 0 || strcmp(type, "-IMPLEX") == 0) {
        eleImplex3d = true;
    }
    else if (strcmp(type, "-pm") == 0 || strcmp(type, "-PM") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &elePmP0_3d) < 0 || elePmP0_3d <= 0.0) {
            opserr << "WARNING LadrunoDispBeamColumn3d: -pm needs P0 > 0\n";
            return 0;
        }
        double optAlpha = elePmAlpha_3d;
        if (ladrunoPeekNextDouble3d(optAlpha)) {
            elePmAlpha_3d = optAlpha;
        }
    }
    else {
        opserr << "WARNING LadrunoDispBeamColumn3d: unknown option '" << type << "'\n";
        return 0;
    }
    }

    // Ladruno (ADR 34): the coupled biaxial cohesive replaces the two block-diagonal scalar laws —
    // reject mixing it with -hinge/-hingeY.
    if ((hingeMatC != 0 || hingeMatCJ != 0) &&
        (hingeMatZ != 0 || hingeMatY != 0 || hingeMatZJ != 0 || hingeMatYJ != 0)) {
        opserr << "WARNING LadrunoDispBeamColumn3d: -hingeBiaxial (coupled) is mutually exclusive "
                  "with -hinge/-hingeY (block-diagonal); use one or the other.\n";
        return 0;
    }

    // Ladruno (ADR 33) PR-3b: the weak-axis jump is an ADD-ON to the strong-axis hinge — the
    // coupled biaxial 2x2 condenses both bending jumps together, so -hingeY requires -hinge.
    if ((hingeMatY != 0 && hingeMatZ == 0) || (hingeMatYJ != 0 && hingeMatZJ == 0)) {
        opserr << "WARNING LadrunoDispBeamColumn3d: -hingeY requires -hinge (the biaxial hinge "
                  "condenses the strong- and weak-axis jumps together).\n";
        return 0;
    }

    // Ladruno (ADR 33) Tier-2: the embedded hinge and the -nl bowing strain are mutually
    // exclusive in v1 — the ½(θz²+θy²) membrane term couples the rotation jump into the axial
    // channel. Reject the combination loudly rather than silently mis-condense.
    if ((hingeMatZ != 0 || hingeMatY != 0 || hingeMatC != 0 ||
         hingeMatZJ != 0 || hingeMatYJ != 0 || hingeMatCJ != 0) && nlGeom != 0) {
        opserr << "WARNING LadrunoDispBeamColumn3d: -hinge and -nl are mutually exclusive in v1 "
                  "(the bowing strain couples the rotation jump into the axial channel); "
                  "use one or the other.\n";
        return 0;
    }

    // check transf
    CrdTransf* theTransf = OPS_getCrdTransf(iData[3]);
    if(theTransf == 0) {
	opserr<<"coord transfomration not found\n";
	return 0;
    }

    // check beam integrataion
    BeamIntegrationRule* theRule = OPS_getBeamIntegrationRule(iData[4]);
    if(theRule == 0) {
	opserr<<"beam integration not found\n";
	return 0;
    }
    BeamIntegration* bi = theRule->getBeamIntegration();
    if(bi == 0) {
	opserr<<"beam integration is null\n";
	return 0;
    }

    // check sections
    const ID& secTags = theRule->getSectionTags();
    if (secTags.Size() <= 0 || secTags.Size() > 20) {
	opserr << "WARNING LadrunoDispBeamColumn3d: number of sections (" << secTags.Size()
	       << ") must be in 1..20\n";
	return 0;
    }
    SectionForceDeformation** sections = new SectionForceDeformation *[secTags.Size()];
    for(int i=0; i<secTags.Size(); i++) {
	sections[i] = OPS_getSectionForceDeformation(secTags(i));
	if(sections[i] == 0) {
	    opserr<<"section "<<secTags(i)<<"not found\n";
		delete [] sections;
	    return 0;
	}
    }

    UniaxialMaterial *hZ_tmp = 0, *hY_tmp = 0, *hZJ_tmp = 0, *hYJ_tmp = 0;
    NDMaterial *hC_tmp = 0, *hCJ_tmp = 0;
    if (eleImplex3d || elePmP0_3d > 0.0) {
        auto cfgUni = [&](UniaxialMaterial *&hPtr, UniaxialMaterial *&tmpPtr) {
            if (hPtr != 0) {
                tmpPtr = hPtr->getCopy();
                if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge*>(tmpPtr)) {
                    if (eleImplex3d) ch->setImplex(true);
                    if (elePmP0_3d > 0.0) ch->setPM(elePmP0_3d, elePmAlpha_3d);
                } else if (Bilin *b = dynamic_cast<Bilin*>(tmpPtr)) {
                    if (eleImplex3d) b->setImplex(true);
                } else if (ModIMKPeakOriented *m = dynamic_cast<ModIMKPeakOriented*>(tmpPtr)) {
                    if (eleImplex3d) m->setImplex(true);
                }
                hPtr = tmpPtr;
            }
        };
        auto cfgND = [&](NDMaterial *&hPtr, NDMaterial *&tmpPtr) {
            if (hPtr != 0) {
                tmpPtr = hPtr->getCopy();
                if (LadrunoCohesiveHingeBiaxial *ch = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(tmpPtr)) {
                    if (eleImplex3d) ch->setImplex(true);
                    if (elePmP0_3d > 0.0) ch->setPM(elePmP0_3d, elePmAlpha_3d);
                }
                hPtr = tmpPtr;
            }
        };
        cfgUni(hingeMatZ, hZ_tmp);
        cfgUni(hingeMatY, hY_tmp);
        cfgUni(hingeMatZJ, hZJ_tmp);
        cfgUni(hingeMatYJ, hYJ_tmp);
        cfgND(hingeMatC, hC_tmp);
        cfgND(hingeMatCJ, hCJ_tmp);
    }
    
    Element *theEle =  new LadrunoDispBeamColumn3d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					    *bi,*theTransf,mass,cmass, theDamping,lchMode,userLch,nlGeom,
					    hingeMatZ,hingeMatY,hingeMatC,hingeMatZJ,hingeMatYJ,hingeMatCJ,
					    dualHingeMode,xiI,xiJ,elePmP0_3d,elePmAlpha_3d);
    delete hZ_tmp;
    delete hY_tmp;
    delete hZJ_tmp;
    delete hYJ_tmp;
    delete hC_tmp;
    delete hCJ_tmp;
    delete [] sections;
    return theEle;
}

LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration &bi,
				   CrdTransf &coordTransf, double r, int cm,
				   Damping *damping, int lchM, double userL, int nlG,
				   UniaxialMaterial *hingeMatZ, UniaxialMaterial *hingeMatY,
				   NDMaterial *hingeMatC,
				   UniaxialMaterial *hingeMatZJ, UniaxialMaterial *hingeMatYJ,
				   NDMaterial *hingeMatCJ,
				   int dualMode, double xiI, double xiJ,
				   double pmP0_, double pmAlpha_)
:Element (tag, ELE_TAG_LadrunoDispBeamColumn3d),
numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2),
Q(12), q(6), rho(r), cMass(cm), parameterID(0), theDamping(0),
current_section_lch(0.0), lchMode(lchM), userLch(userL), nlGeom(nlG),
hingeOn(0), theHingeZ(0), hingeJumpZ(0.0), hingeJumpCommitZ(0.0), hingeKaaZ(0.0), hingeKaaCommitZ(0.0),
hingeMscaleZ(0.0), hingeMscaleCommitZ(0.0),
theHingeY(0), hingeJumpY(0.0), hingeJumpCommitY(0.0), hingeMscaleY(0.0), hingeMscaleCommitY(0.0),
hingeKaaInvZZ(0.0), hingeKaaInvZY(0.0), hingeKaaInvYY(0.0),
hingeKaaInvCommitZZ(0.0), hingeKaaInvCommitZY(0.0), hingeKaaInvCommitYY(0.0),
theHingeC(0),
dualHingeMode(dualMode), hingeXiI(xiI), hingeXiJ(xiJ),
theHingeZJ(0), theHingeYJ(0), theHingeCJ(0),
hingeJumpZJ(0.0), hingeJumpCommitZJ(0.0), hingeJumpYJ(0.0), hingeJumpCommitYJ(0.0),
hingeMscaleZJ(0.0), hingeMscaleCommitZJ(0.0), hingeMscaleYJ(0.0), hingeMscaleCommitYJ(0.0),
pmP0(pmP0_), pmAlpha(pmAlpha_),
postRestoreSync(false), vRestoreValid(false)
{
  for (int k = 0; k < 6; k++) vRestore[k] = 0.0;
  hingeKvZ[0] = hingeKvZ[1] = hingeKvZ[2] = hingeKvZ[3] = hingeKvZ[4] = hingeKvZ[5] = 0.0;
  hingeKvCommitZ[0] = hingeKvCommitZ[1] = hingeKvCommitZ[2] = hingeKvCommitZ[3] = hingeKvCommitZ[4] = hingeKvCommitZ[5] = 0.0;
  hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
  hingeKvCommitY[0] = hingeKvCommitY[1] = hingeKvCommitY[2] = hingeKvCommitY[3] = hingeKvCommitY[4] = hingeKvCommitY[5] = 0.0;
  for (int a = 0; a < 6; a++)
    for (int m = 0; m < 4; m++)
      hingeDualKv[a][m] = hingeDualKvCommit[a][m] = 0.0;
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < 4; n++)
      hingeDualKaaInv[m][n] = hingeDualKaaInvCommit[m][n] = 0.0;
  if (numSections < 0 || numSections > maxNumSections) {
    opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - numSections ("
           << numSections << ") exceeds maxNumSections (" << maxNumSections << ")\n";
    numSections = (numSections < 0) ? 0 : maxNumSections;
  }
  // Ladruno (ADR 33) Tier-2: own a copy of the strong-axis cohesive hinge material
  if (hingeMatZ != 0) {
    theHingeZ = hingeMatZ->getCopy();
    if (theHingeZ == 0) {
      opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to copy hinge material\n";
      exit(-1);
    }
    hingeOn = 1;
  }
  // Ladruno (ADR 33 PR-3b): own a copy of the weak-axis cohesive hinge material (biaxial add-on)
  if (hingeMatY != 0) {
    theHingeY = hingeMatY->getCopy();
    if (theHingeY == 0) {
      opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to copy weak-axis hinge material\n";
      exit(-1);
    }
    hingeOn = 1;
  }
  // Ladruno (ADR 34): own a copy of the coupled biaxial cohesive law (order-2 NDMaterial)
  if (hingeMatC != 0) {
    theHingeC = hingeMatC->getCopy();
    if (theHingeC == 0) {
      opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to copy coupled biaxial hinge material\n";
      exit(-1);
    }
    hingeOn = 1;
  }
  // Ladruno (WP-200): end-J hinge materials
  if (hingeMatZJ != 0) {
    theHingeZJ = hingeMatZJ->getCopy();
    if (theHingeZJ == 0) {
      opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to copy end-J hinge material\n";
      exit(-1);
    }
    hingeOn = 1;
    dualHingeMode = 1;
  }
  if (hingeMatYJ != 0) {
    theHingeYJ = hingeMatYJ->getCopy();
    if (theHingeYJ == 0) {
      opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to copy end-J weak-axis hinge material\n";
      exit(-1);
    }
    hingeOn = 1;
    dualHingeMode = 1;
  }
  if (hingeMatCJ != 0) {
    theHingeCJ = hingeMatCJ->getCopy();
    if (theHingeCJ == 0) {
      opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to copy end-J coupled biaxial hinge material\n";
      exit(-1);
    }
    hingeOn = 1;
    dualHingeMode = 1;
  }

  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
  
  if (theSections == 0) {
    opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to allocate section model pointer\n";
    exit(-1);
  }
  
  for (int i = 0; i < numSections; i++) {
    
    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();
    
    // Check allocation
    if (theSections[i] == 0) {
      opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d -- failed to get a copy of section model\n";
      exit(-1);
    }
  }
  
  beamInt = bi.getCopy();
  
  if (beamInt == 0) {
    opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to copy beam integration\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy3d();
  
  if (crdTransf == 0) {
    opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to copy coordinate transformation\n";
    exit(-1);
  }
  
  if (damping)
  {
    theDamping =(*damping).getCopy();
    
    if (!theDamping) {
      opserr << "LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d - failed to copy damping\n";
      exit(-1);
    }
  }
  
  // Set connected external node IDs
  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;


  theNodes[0] = 0;
  theNodes[1] = 0;

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
}

LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d()
:Element (0, ELE_TAG_LadrunoDispBeamColumn3d),
numSections(0), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2),
Q(12), q(6), rho(0.0), cMass(0), parameterID(0),
theDamping(0),
current_section_lch(0.0), lchMode(0), userLch(0.0), nlGeom(0),
hingeOn(0), theHingeZ(0), hingeJumpZ(0.0), hingeJumpCommitZ(0.0), hingeKaaZ(0.0), hingeKaaCommitZ(0.0),
hingeMscaleZ(0.0), hingeMscaleCommitZ(0.0),
theHingeY(0), hingeJumpY(0.0), hingeJumpCommitY(0.0), hingeMscaleY(0.0), hingeMscaleCommitY(0.0),
hingeKaaInvZZ(0.0), hingeKaaInvZY(0.0), hingeKaaInvYY(0.0),
hingeKaaInvCommitZZ(0.0), hingeKaaInvCommitZY(0.0), hingeKaaInvCommitYY(0.0),
theHingeC(0),
dualHingeMode(0), hingeXiI(0.0), hingeXiJ(1.0),
theHingeZJ(0), theHingeYJ(0), theHingeCJ(0),
hingeJumpZJ(0.0), hingeJumpCommitZJ(0.0), hingeJumpYJ(0.0), hingeJumpCommitYJ(0.0),
hingeMscaleZJ(0.0), hingeMscaleCommitZJ(0.0), hingeMscaleYJ(0.0), hingeMscaleCommitYJ(0.0),
pmP0(0.0), pmAlpha(0.5),
postRestoreSync(false), vRestoreValid(false)
{
  for (int k = 0; k < 6; k++) vRestore[k] = 0.0;
  hingeKvZ[0] = hingeKvZ[1] = hingeKvZ[2] = hingeKvZ[3] = hingeKvZ[4] = hingeKvZ[5] = 0.0;
  hingeKvCommitZ[0] = hingeKvCommitZ[1] = hingeKvCommitZ[2] = hingeKvCommitZ[3] = hingeKvCommitZ[4] = hingeKvCommitZ[5] = 0.0;
  hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
  hingeKvCommitY[0] = hingeKvCommitY[1] = hingeKvCommitY[2] = hingeKvCommitY[3] = hingeKvCommitY[4] = hingeKvCommitY[5] = 0.0;
  for (int a = 0; a < 6; a++)
    for (int m = 0; m < 4; m++)
      hingeDualKv[a][m] = hingeDualKvCommit[a][m] = 0.0;
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < 4; n++)
      hingeDualKaaInv[m][n] = hingeDualKaaInvCommit[m][n] = 0.0;
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  theNodes[0] = 0;
  theNodes[1] = 0;
}

LadrunoDispBeamColumn3d::~LadrunoDispBeamColumn3d()
{    
  for (int i = 0; i < numSections; i++) {
    if (theSections[i])
      delete theSections[i];
  }
  
  // Delete the array of pointers to SectionForceDeformation pointer arrays
  if (theSections)
    delete [] theSections;
  
  if (crdTransf)
    delete crdTransf;

  if (beamInt != 0)
    delete beamInt;

	if (theDamping) delete theDamping;

  if (theHingeZ != 0)   // Ladruno (ADR 33) Tier-2
    delete theHingeZ;

  if (theHingeY != 0)   // Ladruno (ADR 33 PR-3b)
    delete theHingeY;

  if (theHingeC != 0)   // Ladruno (ADR 34)
    delete theHingeC;

  if (theHingeZJ != 0)  // Ladruno (WP-200)
    delete theHingeZJ;

  if (theHingeYJ != 0)
    delete theHingeYJ;

  if (theHingeCJ != 0)
    delete theHingeCJ;
}

int
LadrunoDispBeamColumn3d::getNumExternalNodes() const
{
    return 2;
}

const ID&
LadrunoDispBeamColumn3d::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
LadrunoDispBeamColumn3d::getNodePtrs()
{

    return theNodes;
}

int
LadrunoDispBeamColumn3d::getNumDOF()
{
    return 12;
}

void
LadrunoDispBeamColumn3d::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    if (theNodes[0] == 0 || theNodes[1] == 0) {
	opserr << "FATAL ERROR LadrunoDispBeamColumn3d (tag: " << this->getTag()
	       << "), node not found in domain" << endln;
	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    if (dofNd1 != 6 || dofNd2 != 6) {
	opserr << "FATAL ERROR LadrunoDispBeamColumn3d (tag: " << this->getTag()
	       << "), has differing number of DOFs at its nodes" << endln;
	return;
    }

	if (crdTransf->initialize(theNodes[0], theNodes[1])) {
		// Add some error check
	}

  // initialize the damping
  if (theDamping && theDamping->setDomain(theDomain, 6)) {
    opserr << "LadrunoDispBeamColumn3d::setDomain(): Error initializing damping";  
    exit(0);
  }

	double L = crdTransf->getInitialLength();

	if (L == 0.0) {
		// Add some error check
	}

    this->DomainComponent::setDomain(theDomain);

	this->update();
	if (!theDamping)
	  this->getResistingForce();
}

int
LadrunoDispBeamColumn3d::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping)
  {
    if (theDamping) delete theDamping;

    theDamping =(*damping).getCopy();
    
    if (!theDamping) {
      opserr << "LadrunoDispBeamColumn3d::setDamping -- failed to get copy of damping\n";
      return -1;
    }
    if (theDamping->setDomain(theDomain, 6)) {
      opserr << "LadrunoDispBeamColumn3d::setDamping -- Error initializing damping\n";
      return -2;
    }
  }
  
  return 0;
}

int
LadrunoDispBeamColumn3d::commitState()
{
    int retVal = 0;
    postRestoreSync = false;
    vRestoreValid = false;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "LadrunoDispBeamColumn3d::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    if (theDamping) retVal += theDamping->commitState();

    // Ladruno (ADR 33) Tier-2: lock in the converged strong-axis jump + the cohesive law's
    // irreversible history. Commit is the ONLY place the jump advances (no resurrection on revert).
    if (hingeOn) {
      hingeJumpCommitZ = hingeJumpZ;
      hingeKaaCommitZ = hingeKaaZ;
      for (int m = 0; m < 6; m++) hingeKvCommitZ[m] = hingeKvZ[m];
      hingeMscaleCommitZ = hingeMscaleZ;
      if (theHingeZ) retVal += theHingeZ->commitState();
      if (theHingeY) {                       // PR-3b biaxial
        hingeJumpCommitY = hingeJumpY;
        for (int m = 0; m < 6; m++) hingeKvCommitY[m] = hingeKvY[m];
        hingeMscaleCommitY = hingeMscaleY;
        hingeKaaInvCommitZZ = hingeKaaInvZZ;
        hingeKaaInvCommitZY = hingeKaaInvZY;
        hingeKaaInvCommitYY = hingeKaaInvYY;
        retVal += theHingeY->commitState();
      }
      if (theHingeC) {                       // ADR 34 coupled (jumps committed, plus the law's history)
        hingeJumpCommitY = hingeJumpY;
        for (int m = 0; m < 6; m++) hingeKvCommitY[m] = hingeKvY[m];
        hingeMscaleCommitY = hingeMscaleY;
        hingeKaaInvCommitZZ = hingeKaaInvZZ;
        hingeKaaInvCommitZY = hingeKaaInvZY;
        hingeKaaInvCommitYY = hingeKaaInvYY;
        retVal += theHingeC->commitState();
      }
      if (dualHingeMode) {
        hingeJumpCommitZJ = hingeJumpZJ;
        hingeJumpCommitYJ = hingeJumpYJ;
        hingeMscaleCommitZJ = hingeMscaleZJ;
        hingeMscaleCommitYJ = hingeMscaleYJ;
        for (int a = 0; a < 6; a++)
          for (int m = 0; m < 4; m++)
            hingeDualKvCommit[a][m] = hingeDualKv[a][m];
        for (int m = 0; m < 4; m++)
          for (int n = 0; n < 4; n++)
            hingeDualKaaInvCommit[m][n] = hingeDualKaaInv[m][n];
        if (theHingeZJ) retVal += theHingeZJ->commitState();
        if (theHingeYJ) retVal += theHingeYJ->commitState();
        if (theHingeCJ) retVal += theHingeCJ->commitState();
      }
    }

    return retVal;
}

int
LadrunoDispBeamColumn3d::revertToLastCommit()
{
    int retVal = 0;
    postRestoreSync = true;
    vRestoreValid = false;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    if (theDamping) retVal += theDamping->revertToLastCommit();

    // Ladruno (ADR 33) Tier-2: restore the jump (and the cohesive law) to the committed state
    // so a rejected step does not leave a stale/open hinge.
    if (hingeOn) {
      hingeJumpZ = hingeJumpCommitZ;
      hingeKaaZ = hingeKaaCommitZ;
      for (int m = 0; m < 6; m++) hingeKvZ[m] = hingeKvCommitZ[m];
      hingeMscaleZ = hingeMscaleCommitZ;
      if (theHingeZ && !ldbc3dHingeIsImplex(theHingeZ)) retVal += theHingeZ->revertToLastCommit();
      if (theHingeY) {                       // PR-3b biaxial
        hingeJumpY = hingeJumpCommitY;
        for (int m = 0; m < 6; m++) hingeKvY[m] = hingeKvCommitY[m];
        hingeMscaleY = hingeMscaleCommitY;
        hingeKaaInvZZ = hingeKaaInvCommitZZ;
        hingeKaaInvZY = hingeKaaInvCommitZY;
        hingeKaaInvYY = hingeKaaInvCommitYY;
        if (!ldbc3dHingeIsImplex(theHingeY)) retVal += theHingeY->revertToLastCommit();
      }
      if (theHingeC) {                       // ADR 34 coupled
        hingeJumpY = hingeJumpCommitY;
        for (int m = 0; m < 6; m++) hingeKvY[m] = hingeKvCommitY[m];
        hingeMscaleY = hingeMscaleCommitY;
        hingeKaaInvZZ = hingeKaaInvCommitZZ;
        hingeKaaInvZY = hingeKaaInvCommitZY;
        hingeKaaInvYY = hingeKaaInvCommitYY;
        if (!ldbc3dHingeBiaxialIsImplex(theHingeC)) retVal += theHingeC->revertToLastCommit();
      }
      if (dualHingeMode) {
        hingeJumpZJ = hingeJumpCommitZJ;
        hingeJumpYJ = hingeJumpCommitYJ;
        hingeMscaleZJ = hingeMscaleCommitZJ;
        hingeMscaleYJ = hingeMscaleCommitYJ;
        for (int a = 0; a < 6; a++)
          for (int m = 0; m < 4; m++)
            hingeDualKv[a][m] = hingeDualKvCommit[a][m];
        for (int m = 0; m < 4; m++)
          for (int n = 0; n < 4; n++)
            hingeDualKaaInv[m][n] = hingeDualKaaInvCommit[m][n];
        if (theHingeZJ && !ldbc3dHingeIsImplex(theHingeZJ)) retVal += theHingeZJ->revertToLastCommit();
        if (theHingeYJ && !ldbc3dHingeIsImplex(theHingeYJ)) retVal += theHingeYJ->revertToLastCommit();
        if (theHingeCJ && !ldbc3dHingeBiaxialIsImplex(theHingeCJ)) retVal += theHingeCJ->revertToLastCommit();
      }
    }

    if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
      retVal += this->update();
      if (!theDamping) this->getResistingForce();
    }

    return retVal;
}

int
LadrunoDispBeamColumn3d::revertToStart()
{
    int retVal = 0;
    postRestoreSync = false;
    vRestoreValid = false;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    if (theDamping) retVal += theDamping->revertToStart();

    // Ladruno (ADR 33) Tier-2: reset the embedded hinge to a closed, undamaged state
    if (hingeOn) {
      hingeJumpZ = 0.0;
      hingeJumpCommitZ = 0.0;
      hingeKaaZ = 0.0;
      hingeKaaCommitZ = 0.0;
      hingeKvZ[0] = hingeKvZ[1] = hingeKvZ[2] = hingeKvZ[3] = hingeKvZ[4] = hingeKvZ[5] = 0.0;
      hingeKvCommitZ[0] = hingeKvCommitZ[1] = hingeKvCommitZ[2] = hingeKvCommitZ[3] = hingeKvCommitZ[4] = hingeKvCommitZ[5] = 0.0;
      hingeMscaleZ = 0.0;
      hingeMscaleCommitZ = 0.0;
      if (theHingeZ) retVal += theHingeZ->revertToStart();
      if (theHingeY) {                       // PR-3b biaxial
        hingeJumpY = 0.0;
        hingeJumpCommitY = 0.0;
        hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
        hingeKvCommitY[0] = hingeKvCommitY[1] = hingeKvCommitY[2] = hingeKvCommitY[3] = hingeKvCommitY[4] = hingeKvCommitY[5] = 0.0;
        hingeMscaleY = 0.0;
        hingeMscaleCommitY = 0.0;
        hingeKaaInvZZ = hingeKaaInvZY = hingeKaaInvYY = 0.0;
        hingeKaaInvCommitZZ = hingeKaaInvCommitZY = hingeKaaInvCommitYY = 0.0;
        retVal += theHingeY->revertToStart();
      }
      if (theHingeC) {                       // ADR 34 coupled
        hingeJumpY = 0.0;
        hingeJumpCommitY = 0.0;
        hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
        hingeKvCommitY[0] = hingeKvCommitY[1] = hingeKvCommitY[2] = hingeKvCommitY[3] = hingeKvCommitY[4] = hingeKvCommitY[5] = 0.0;
        hingeMscaleY = 0.0;
        hingeMscaleCommitY = 0.0;
        hingeKaaInvZZ = hingeKaaInvZY = hingeKaaInvYY = 0.0;
        hingeKaaInvCommitZZ = hingeKaaInvCommitZY = hingeKaaInvCommitYY = 0.0;
        retVal += theHingeC->revertToStart();
      }
      if (dualHingeMode) {
        hingeJumpZJ = 0.0;
        hingeJumpCommitZJ = 0.0;
        hingeJumpYJ = 0.0;
        hingeJumpCommitYJ = 0.0;
        hingeMscaleZJ = 0.0;
        hingeMscaleCommitZJ = 0.0;
        hingeMscaleYJ = 0.0;
        hingeMscaleCommitYJ = 0.0;
        for (int a = 0; a < 6; a++)
          for (int m = 0; m < 4; m++)
            hingeDualKv[a][m] = hingeDualKvCommit[a][m] = 0.0;
        for (int m = 0; m < 4; m++)
          for (int n = 0; n < 4; n++)
            hingeDualKaaInv[m][n] = hingeDualKaaInvCommit[m][n] = 0.0;
        if (theHingeZJ) retVal += theHingeZJ->revertToStart();
        if (theHingeYJ) retVal += theHingeYJ->revertToStart();
        if (theHingeCJ) retVal += theHingeCJ->revertToStart();
      }
    }

    if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
      retVal += this->update();
      if (!theDamping) this->getResistingForce();
    }

    return retVal;
}

int
LadrunoDispBeamColumn3d::update(void)
{
  int err = 0;

  // Ladruno (ADR 32): defensively re-assert the active element (see 2D sibling).
  ops_TheActiveElement = this;

  // Update the transformation
  crdTransf->update();

  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();
  if (postRestoreSync) {
    if (!vRestoreValid) {
      for (int k = 0; k < 6; k++) vRestore[k] = v(k);
      vRestoreValid = true;
    } else if (v(0) != vRestore[0] || v(1) != vRestore[1] || v(2) != vRestore[2] ||
               v(3) != vRestore[3] || v(4) != vRestore[4] || v(5) != vRestore[5]) {
      postRestoreSync = false;
    }
  }

  double L = crdTransf->getInitialLength();
  if (L <= 0.0)
    return -1;
  double oneOverL = 1.0/L;

  // Ladruno (ADR 33 / WP-200): when embedded hinge(s) are active, dispatch to dual-hinge
  // solver (when dualHingeMode is set) or the legacy single-hinge solvers.
  if (hingeOn) {
    if (dualHingeMode)
      return solveDualHingeJump3d(v, L);
    if (theHingeY != 0 || theHingeC != 0)   // PR-3b block-diagonal OR ADR-34 coupled -> biaxial solve
      return solveHingeJumpBiaxial(v, L);
    return solveHingeJump(v, L);
  }

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  // Ladruno (ADR 32): section weights for the per-IP tributary length wt[i]*L.
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Vector e(workArea, order);

    double xi6 = 6.0*xi[i];
    // Ladruno (ADR 32): biaxial transverse slopes for the ½(θz²+θy²) bowing term
    double zeta = xi[i];
    double thetaz = (3.0*zeta*zeta - 4.0*zeta + 1.0)*v(1) + (3.0*zeta*zeta - 2.0*zeta)*v(2);
    double thetay = (3.0*zeta*zeta - 4.0*zeta + 1.0)*v(3) + (3.0*zeta*zeta - 2.0*zeta)*v(4);

    int j;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	// Ladruno (ADR 32): nlGeom adds the ½(θz²+θy²) membrane/bowing strain (DispBeamColumnNL3d)
	e(j) = oneOverL*v(0) + (nlGeom ? 0.5*thetaz*thetaz + 0.5*thetay*thetay : 0.0);
	break;
      case SECTION_RESPONSE_MZ:
	e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2));
	break;
      case SECTION_RESPONSE_MY:
	e(j) = oneOverL*((xi6-4.0)*v(3) + (xi6-2.0)*v(4));
	break;
      case SECTION_RESPONSE_T:
	e(j) = oneOverL*v(5);
	break;
      default:
	e(j) = 0.0;
	break;
      }
    }

    // Ladruno (ADR 32) Tier-1: per-IP characteristic length (reference), set immediately
    // before setTrialSectionDeformation so the crack-band material latches the right band.
    if (lchMode == 1)        current_section_lch = L;
    else if (lchMode == 2)   current_section_lch = userLch;
    else                     current_section_lch = wt[i]*L;

    // Set the section deformations
    err += theSections[i]->setTrialSectionDeformation(e);
  }

  if (err != 0) {
    LadrunoOutputGuard og;
    opserr << "LadrunoDispBeamColumn3d::update() - failed setTrialSectionDeformations()\n";
    return err;
  }
  return 0;
}

// Ladruno (ADR 33) Tier-2 PR-3a: inner Newton on the scalar strong-axis rotation jump alpha_z.
//
// 3D analog of the 2D LadrunoDispBeamColumn2d::solveHingeJump, restricted to the Mz axis on the
// 6-DOF basic system. The bulk section sees kappa_z_bulk = B_z*v + Gbar_z*alpha_z (Gbar_z = -1/L),
// the My/T/axial channels stay linear/untouched. Enhancement equilibrium and tangent:
//   h_z(alpha_z) = -sum_k wt_k*M_sec,k(MZ) + M_coh_z(alpha_z) = 0
//   K_aa,z       = (1/L) sum_k wt_k*ks(MZ,MZ) + dM_coh_z/dalpha_z   (guarded; indefinite on softening)
// The condensation 6-vector K_v-alpha_z includes the CROSS-tangent rows (ks(P,MZ), ks(MY,MZ),
// ks(T,MZ)) so the condensed off-diagonals of the 6x6 basic kb are right. The basic FORCE needs no
// correction (at h_z=0 the sections hold kappa_z_bulk, so q = sum wt*B^T*s IS the condensed force).
int
LadrunoDispBeamColumn3d::solveHingeJump(const Vector &v, double L)
{
  if (L <= 0.0)
    return -1;
  double oneOverL = 1.0/L;
  const double Gbar = -oneOverL;   // bounded incompatible-mode enhancement on the Mz channel

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  double alpha = hingeJumpCommitZ;   // warm start from the committed jump (path-correct)
  LadrunoCohesiveHinge *chZ = dynamic_cast<LadrunoCohesiveHinge*>(theHingeZ);
  const bool freezeZ = postRestoreSync && ldbc3dHingeIsImplex(theHingeZ);

  const int maxIter = 50;
  double Kva[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double Kaa = 0.0;
  double damp = 1.0;
  double hPrev = 1.0e300;

  for (int iter = 0; iter < maxIter; iter++) {

    int err = 0;
    double sumM = 0.0;     // sum_k wt_k * M_sec,k(MZ)
    double sumP = 0.0;     // sum_k wt_k * N_sec,k(P)
    double sumEI = 0.0;    // sum_k wt_k * ks(MZ,MZ)
    for (int a = 0; a < 6; a++) Kva[a] = 0.0;

    // --- set every IP's bulk-enhanced section strain at the current alpha_z ---
    for (int i = 0; i < numSections; i++) {

      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      double xi6 = 6.0*xi[i];

      Vector e(workArea, order);
      for (int j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          e(j) = oneOverL*v(0); break;
        case SECTION_RESPONSE_MZ:
          // strong-axis bulk curvature = B_z*v + Gbar_z*alpha_z  (the -alpha_z/L offset)
          e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)) + Gbar*alpha; break;
        case SECTION_RESPONSE_MY:
          e(j) = oneOverL*((xi6-4.0)*v(3) + (xi6-2.0)*v(4)); break;
        case SECTION_RESPONSE_T:
          e(j) = oneOverL*v(5); break;
        default:
          e(j) = 0.0; break;
        }
      }

      // Tier-1 lch latch (identical to the no-hinge update loop)
      if (lchMode == 1)       current_section_lch = L;
      else if (lchMode == 2)  current_section_lch = userLch;
      else                    current_section_lch = wt[i]*L;

      err += theSections[i]->setTrialSectionDeformation(e);

      // accumulate residual / condensation terms from this IP (scan for the channel indices)
      const Vector &s  = theSections[i]->getStressResultant();
      const Matrix &ks = theSections[i]->getSectionTangent();
      int idxP = -1, idxMZ = -1, idxMY = -1, idxT = -1;
      for (int j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:  idxP  = j; break;
        case SECTION_RESPONSE_MZ: idxMZ = j; break;
        case SECTION_RESPONSE_MY: idxMY = j; break;
        case SECTION_RESPONSE_T:  idxT  = j; break;
        default: break;
        }
      }
      if (idxP >= 0) sumP += wt[i]*s(idxP);
      if (idxMZ >= 0) {
        double EI = ks(idxMZ, idxMZ);
        sumM  += wt[i]*s(idxMZ);
        sumEI += wt[i]*EI;
        // K_v-alpha_z column: basic-DOF coupling of alpha_z (driver column ks(*,MZ))
        Kva[1] += -oneOverL*wt[i]*(xi6-4.0)*EI;
        Kva[2] += -oneOverL*wt[i]*(xi6-2.0)*EI;
        if (idxP  >= 0) Kva[0] += -oneOverL*wt[i]*ks(idxP,  idxMZ);  // axial-flexural coupling
        if (idxMY >= 0) {                                            // cross-axis (My-Mz) coupling
          Kva[3] += -oneOverL*wt[i]*(xi6-4.0)*ks(idxMY, idxMZ);
          Kva[4] += -oneOverL*wt[i]*(xi6-2.0)*ks(idxMY, idxMZ);
        }
        if (idxT  >= 0) Kva[5] += -oneOverL*wt[i]*ks(idxT,  idxMZ);  // torsion-flexural coupling
      }
    }

    if (err != 0) {
      LadrunoOutputGuard og;
      opserr << "LadrunoDispBeamColumn3d::solveHingeJump() - failed setTrialSectionDeformations()\n";
      return err;
    }

    sumP += q0[0];
    sumM += -(1.0 - hingeXiI)*q0[1] + hingeXiI*q0[2];

    if (chZ != 0) chZ->setAxialForce(sumP);
    double fPMZ = (chZ == 0 && pmP0 > 0.0) ? ladrunoDisp3dScalePM(sumP, pmP0, pmAlpha) : 1.0;

    // --- cohesive law at the current jump: M_coh_z(alpha_z), dM_coh_z/dalpha_z ---
    if (!freezeZ)
      theHingeZ->setTrialStrain(alpha);
    double Mcoh = fPMZ * theHingeZ->getStress();
    double Kcoh = fPMZ * theHingeZ->getTangent();

    double h = freezeZ ? 0.0 : (-sumM + Mcoh); // enhancement residual
    Kaa      = oneOverL*sumEI + Kcoh;        // (1/L) sum wt*EI + dM_coh/dalpha

    // GUARDED reciprocal: floor |K_aa| against the positive bulk term (sign-discontinuous at
    // activation / indefinite on the softening branch).
    double KaaFloor = 1.0e-8*(fabs(oneOverL*sumEI) + fabs(Kcoh)) + 1.0e-300;
    double KaaSolve = Kaa;
    if (fabs(KaaSolve) < KaaFloor)
      KaaSolve = (KaaSolve < 0.0) ? -KaaFloor : KaaFloor;

    // STABLE moment scale (does not collapse when a fully-broken LINEAR hinge carries M -> 0)
    double Mhere = fabs(Mcoh) + fabs(sumM);
    if (Mhere > hingeMscaleZ) hingeMscaleZ = Mhere;
    if (fabs(h) <= 1.0e-11*hingeMscaleZ + 1.0e-12) {
      hingeKaaZ = KaaSolve;
      for (int a = 0; a < 6; a++) hingeKvZ[a] = Kva[a];
      hingeJumpZ = alpha;
      if (postRestoreSync) {
        hingeKaaCommitZ = hingeKaaZ;
        for (int a = 0; a < 6; a++) hingeKvCommitZ[a] = hingeKvZ[a];
      }
      return 0;
    }

    double hnorm = fabs(h) / (hingeMscaleZ + 1.0e-12);
    if (hnorm > hPrev) damp *= 0.5;
    else               damp = (damp * 1.5 < 1.0) ? damp * 1.5 : 1.0;
    hPrev = hnorm;

    alpha -= damp * (h/KaaSolve);            // guarded Newton step on the scalar jump
  }

  hingeKaaZ = (fabs(Kaa) < 1.0e-300) ? 1.0e-300 : Kaa;
  for (int a = 0; a < 6; a++) hingeKvZ[a] = Kva[a];
  hingeJumpZ = alpha;
  {
    LadrunoOutputGuard og;
    opserr << "WARNING LadrunoDispBeamColumn3d (tag " << this->getTag()
           << "): embedded-hinge inner Newton did not converge in " << maxIter << " iters\n";
  }
  return 0;
}

// Ladruno (ADR 33) PR-3b: eigenvalue-floored inverse of a symmetric 2x2 [[a,b],[b,d]].
static void
ladrunoFlooredInv2x2(double a, double b, double d, double floor,
                     double &iZZ, double &iZY, double &iYY)
{
  double half = 0.5*(a + d);
  double diff = 0.5*(a - d);
  double rad  = sqrt(diff*diff + b*b);
  double lam1 = half + rad;   // q1: eigenvector (lam1 - d, b) (robust; for symmetric 2x2)
  double lam2 = half - rad;   // q2: the orthonormal complement of q1

  double q1x, q1y;
  if (fabs(b) > 1.0e-300) {
    q1x = lam1 - d;  q1y = b;
  } else {
    // already diagonal: principal axes ARE the coordinate axes
    if (a >= d) { q1x = 1.0; q1y = 0.0; }
    else        { q1x = 0.0; q1y = 1.0; }
  }
  double n1 = sqrt(q1x*q1x + q1y*q1y);
  if (n1 < 1.0e-300) { q1x = 1.0; q1y = 0.0; n1 = 1.0; }
  q1x /= n1; q1y /= n1;
  double q2x = -q1y, q2y = q1x;     // orthonormal complement

  double f1 = (fabs(lam1) < floor) ? ((lam1 < 0.0) ? -floor : floor) : lam1;
  double f2 = (fabs(lam2) < floor) ? ((lam2 < 0.0) ? -floor : floor) : lam2;
  double i1 = 1.0/f1, i2 = 1.0/f2;

  iZZ = i1*q1x*q1x + i2*q2x*q2x;
  iZY = i1*q1x*q1y + i2*q2x*q2y;
  iYY = i1*q1y*q1y + i2*q2y*q2y;
}

// Ladruno (WP-200): eigenvalue-floored inverse of the active submatrix of a symmetric 4x4 matrix
// using Jacobi eigendecomposition. Inactive DOFs (active[m] == false) have their inverse rows/cols
// set to 0.0. Active eigenvalues have their magnitudes floored to `floor` while preserving sign.
static void
ladrunoFlooredInv4x4(const double Kaa[4][4], const bool active[4], double floor,
                     double invKaa[4][4])
{
  double A[4][4];
  double V[4][4];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      V[i][j] = (i == j) ? 1.0 : 0.0;
      if (active[i] && active[j])
        A[i][j] = 0.5 * (Kaa[i][j] + Kaa[j][i]);
      else
        A[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  for (int sweep = 0; sweep < 15; sweep++) {
    double off = 0.0;
    for (int p = 0; p < 3; p++)
      for (int q = p + 1; q < 4; q++)
        off += fabs(A[p][q]);
    if (off < 1.0e-30)
      break;
    for (int p = 0; p < 3; p++) {
      for (int q = p + 1; q < 4; q++) {
        double apq = A[p][q];
        if (fabs(apq) < 1.0e-300)
          continue;
        double diff = A[q][q] - A[p][p];
        double t;
        if (fabs(apq) < 1.0e-18 * fabs(diff)) {
          t = apq / diff;
        } else {
          double phi = 0.5 * diff / apq;
          t = 1.0 / (fabs(phi) + sqrt(1.0 + phi * phi));
          if (phi < 0.0) t = -t;
        }
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;
        double tau = s / (1.0 + c);
        A[p][q] = 0.0;
        A[q][p] = 0.0;
        A[p][p] -= t * apq;
        A[q][q] += t * apq;
        for (int r = 0; r < 4; r++) {
          if (r != p && r != q) {
            double arp = A[r][p];
            double arq = A[r][q];
            A[r][p] = A[p][r] = arp - s * (arq + tau * arp);
            A[r][q] = A[q][r] = arq + s * (arp - tau * arq);
          }
          double vrp = V[r][p];
          double vrq = V[r][q];
          V[r][p] = vrp - s * (vrq + tau * vrp);
          V[r][q] = vrq + s * (vrp - tau * vrq);
        }
      }
    }
  }

  double invLam[4];
  for (int k = 0; k < 4; k++) {
    double lam = A[k][k];
    double f = (fabs(lam) < floor) ? ((lam < 0.0) ? -floor : floor) : lam;
    invLam[k] = 1.0 / f;
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (!active[i] || !active[j]) {
        invKaa[i][j] = 0.0;
      } else {
        double sum = 0.0;
        for (int k = 0; k < 4; k++)
          sum += V[i][k] * invLam[k] * V[j][k];
        invKaa[i][j] = sum;
      }
    }
  }
}

int
LadrunoDispBeamColumn3d::solveHingeJumpBiaxial(const Vector &v, double L)
{
  if (L <= 0.0)
    return -1;
  double oneOverL = 1.0/L;
  const double Gbar = -oneOverL;   // bounded incompatible-mode enhancement (per bending axis)

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  double az = hingeJumpCommitZ;    // warm start from the committed jumps (path-correct)
  double ay = hingeJumpCommitY;
  LadrunoCohesiveHinge *chZ = dynamic_cast<LadrunoCohesiveHinge*>(theHingeZ);
  LadrunoCohesiveHinge *chY = dynamic_cast<LadrunoCohesiveHinge*>(theHingeY);
  LadrunoCohesiveHingeBiaxial *chC = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeC);

  const int maxIter = 50;
  double KvZ[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double KvY[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double iZZ = 0.0, iZY = 0.0, iYY = 0.0;

  double damp = 1.0;
  double hnormPrev = 1.0e300;

  for (int iter = 0; iter < maxIter; iter++) {

    int err = 0;
    double sumMz = 0.0, sumMy = 0.0;       // sum wt*M_sec(MZ/MY)
    double sumP  = 0.0;                     // sum wt*N_sec(P)
    double sumEIz = 0.0, sumEIy = 0.0;     // sum wt*ks(MZ,MZ)/ks(MY,MY)
    double sumCpl = 0.0;                    // sum wt*½(ks(MZ,MY)+ks(MY,MZ))
    for (int a = 0; a < 6; a++) { KvZ[a] = 0.0; KvY[a] = 0.0; }

    // --- set every IP's bulk-enhanced section strain at the current (alpha_z, alpha_y) ---
    for (int i = 0; i < numSections; i++) {

      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      double xi6 = 6.0*xi[i];

      Vector e(workArea, order);
      for (int j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          e(j) = oneOverL*v(0); break;
        case SECTION_RESPONSE_MZ:
          e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)) + Gbar*az; break;
        case SECTION_RESPONSE_MY:
          e(j) = oneOverL*((xi6-4.0)*v(3) + (xi6-2.0)*v(4)) + Gbar*ay; break;
        case SECTION_RESPONSE_T:
          e(j) = oneOverL*v(5); break;
        default:
          e(j) = 0.0; break;
        }
      }

      // Tier-1 lch latch (identical to the no-hinge update loop)
      if (lchMode == 1)       current_section_lch = L;
      else if (lchMode == 2)  current_section_lch = userLch;
      else                    current_section_lch = wt[i]*L;

      err += theSections[i]->setTrialSectionDeformation(e);

      const Vector &s  = theSections[i]->getStressResultant();
      const Matrix &ks = theSections[i]->getSectionTangent();
      int idxP = -1, idxMZ = -1, idxMY = -1, idxT = -1;
      for (int j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:  idxP  = j; break;
        case SECTION_RESPONSE_MZ: idxMZ = j; break;
        case SECTION_RESPONSE_MY: idxMY = j; break;
        case SECTION_RESPONSE_T:  idxT  = j; break;
        default: break;
        }
      }
      if (idxP >= 0) sumP += wt[i]*s(idxP);
      // strong-axis (alpha_z) column of K_v-alpha + the Mz residual/diagonal terms
      if (idxMZ >= 0) {
        double EIz = ks(idxMZ, idxMZ);
        sumMz  += wt[i]*s(idxMZ);
        sumEIz += wt[i]*EIz;
        KvZ[1] += -oneOverL*wt[i]*(xi6-4.0)*EIz;
        KvZ[2] += -oneOverL*wt[i]*(xi6-2.0)*EIz;
        if (idxP  >= 0) KvZ[0] += -oneOverL*wt[i]*ks(idxP,  idxMZ);
        if (idxMY >= 0) {
          KvZ[3] += -oneOverL*wt[i]*(xi6-4.0)*ks(idxMY, idxMZ);
          KvZ[4] += -oneOverL*wt[i]*(xi6-2.0)*ks(idxMY, idxMZ);
        }
        if (idxT  >= 0) KvZ[5] += -oneOverL*wt[i]*ks(idxT,  idxMZ);
      }
      // weak-axis (alpha_y) column of K_v-alpha + the My residual/diagonal terms
      if (idxMY >= 0) {
        double EIy = ks(idxMY, idxMY);
        sumMy  += wt[i]*s(idxMY);
        sumEIy += wt[i]*EIy;
        KvY[3] += -oneOverL*wt[i]*(xi6-4.0)*EIy;
        KvY[4] += -oneOverL*wt[i]*(xi6-2.0)*EIy;
        if (idxP  >= 0) KvY[0] += -oneOverL*wt[i]*ks(idxP,  idxMY);
        if (idxMZ >= 0) {
          KvY[1] += -oneOverL*wt[i]*(xi6-4.0)*ks(idxMZ, idxMY);
          KvY[2] += -oneOverL*wt[i]*(xi6-2.0)*ks(idxMZ, idxMY);
        }
        if (idxT  >= 0) KvY[5] += -oneOverL*wt[i]*ks(idxT,  idxMY);
      }
      // symmetric bending-bending coupling K_aa(z,y) = (1/L)sum wt*½(ks(MZ,MY)+ks(MY,MZ))
      if (idxMZ >= 0 && idxMY >= 0)
        sumCpl += wt[i]*0.5*(ks(idxMZ, idxMY) + ks(idxMY, idxMZ));
    }

    if (err != 0) {
      LadrunoOutputGuard og;
      opserr << "LadrunoDispBeamColumn3d::solveHingeJumpBiaxial() - failed setTrialSectionDeformations()\n";
      return err;
    }

    sumP  += q0[0];
    sumMz += -(1.0 - hingeXiI)*q0[1] + hingeXiI*q0[2];
    sumMy += -(1.0 - hingeXiI)*q0[3] + hingeXiI*q0[4];

    // --- cohesive moments + cohesive tangent at the current jumps ---
    double Mz, My, Czz, Cyy, Czy;
    if (theHingeC != 0) {
      if (chC != 0) chC->setAxialForce(sumP);
      if (!(postRestoreSync && ldbc3dHingeBiaxialIsImplex(theHingeC))) {
        static thread_local Vector jumpVec(2);
        jumpVec(0) = az; jumpVec(1) = ay;
        theHingeC->setTrialStrain(jumpVec);
      }
      const Vector &Mvec = theHingeC->getStress();
      const Matrix &Kc   = theHingeC->getTangent();
      Mz = Mvec(0); My = Mvec(1);
      Czz = Kc(0,0); Cyy = Kc(1,1);
      Czy = 0.5*(Kc(0,1) + Kc(1,0));   // symmetrize the cohesive coupling
    } else {
      if (chZ != 0) chZ->setAxialForce(sumP);
      if (chY != 0) chY->setAxialForce(sumP);
      double fPMZ = (chZ == 0 && pmP0 > 0.0) ? ladrunoDisp3dScalePM(sumP, pmP0, pmAlpha) : 1.0;
      double fPMY = (chY == 0 && pmP0 > 0.0) ? ladrunoDisp3dScalePM(sumP, pmP0, pmAlpha) : 1.0;
      if (!(postRestoreSync && ldbc3dHingeIsImplex(theHingeZ)))
        theHingeZ->setTrialStrain(az);
      Mz = fPMZ * theHingeZ->getStress();  Czz = fPMZ * theHingeZ->getTangent();
      if (!(postRestoreSync && ldbc3dHingeIsImplex(theHingeY)))
        theHingeY->setTrialStrain(ay);
      My = fPMY * theHingeY->getStress();  Cyy = fPMY * theHingeY->getTangent();
      Czy = 0.0;
    }

    double hz = -sumMz + Mz;
    double hy = -sumMy + My;

    double Kaa_zz = oneOverL*sumEIz + Czz;
    double Kaa_yy = oneOverL*sumEIy + Cyy;
    double Kaa_zy = oneOverL*sumCpl + Czy;

    double floor = 1.0e-8*(fabs(oneOverL*sumEIz) + fabs(oneOverL*sumEIy) + fabs(Czz) + fabs(Cyy)) + 1.0e-300;
    ladrunoFlooredInv2x2(Kaa_zz, Kaa_zy, Kaa_yy, floor, iZZ, iZY, iYY);

    double MhereZ = fabs(Mz) + fabs(sumMz);
    double MhereY = fabs(My) + fabs(sumMy);
    if (MhereZ > hingeMscaleZ) hingeMscaleZ = MhereZ;
    if (MhereY > hingeMscaleY) hingeMscaleY = MhereY;

    if (postRestoreSync && (ldbc3dHingeBiaxialIsImplex(theHingeC) ||
                            ldbc3dHingeIsImplex(theHingeZ) ||
                            ldbc3dHingeIsImplex(theHingeY))) {
      hingeJumpZ = az; hingeJumpY = ay;
      for (int a = 0; a < 6; a++) { hingeKvZ[a] = KvZ[a]; hingeKvY[a] = KvY[a]; }
      hingeKaaInvZZ = iZZ; hingeKaaInvZY = iZY; hingeKaaInvYY = iYY;
      return 0;
    }

    bool okZ = fabs(hz) <= 1.0e-11*hingeMscaleZ + 1.0e-12;
    bool okY = fabs(hy) <= 1.0e-11*hingeMscaleY + 1.0e-12;
    if (okZ && okY) {
      hingeJumpZ = az; hingeJumpY = ay;
      for (int a = 0; a < 6; a++) { hingeKvZ[a] = KvZ[a]; hingeKvY[a] = KvY[a]; }
      hingeKaaInvZZ = iZZ; hingeKaaInvZY = iZY; hingeKaaInvYY = iYY;
      return 0;
    }

    double hnorm = sqrt((hz/(hingeMscaleZ+1.0e-12))*(hz/(hingeMscaleZ+1.0e-12))
                      + (hy/(hingeMscaleY+1.0e-12))*(hy/(hingeMscaleY+1.0e-12)));
    if (hnorm > hnormPrev) damp *= 0.5;
    else                   damp = (damp*1.5 < 1.0) ? damp*1.5 : 1.0;
    hnormPrev = hnorm;

    // guarded coupled Newton step: dalpha = -damp * K_aa_floored^{-1} h
    az -= damp*(iZZ*hz + iZY*hy);
    ay -= damp*(iZY*hz + iYY*hy);
  }

  hingeJumpZ = az; hingeJumpY = ay;
  for (int a = 0; a < 6; a++) { hingeKvZ[a] = KvZ[a]; hingeKvY[a] = KvY[a]; }
  hingeKaaInvZZ = iZZ; hingeKaaInvZY = iZY; hingeKaaInvYY = iYY;
  {
    LadrunoOutputGuard og;
    opserr << "WARNING LadrunoDispBeamColumn3d (tag " << this->getTag()
           << "): biaxial embedded-hinge inner Newton did not converge in " << maxIter << " iters\n";
  }
  return 0;
}

// Ladruno (WP-200): 4-DOF dual-hinge inner Newton for [alpha_zI, alpha_yI, alpha_zJ, alpha_yJ]
// with Hermite end-rotation modes at normalized locations xi_I, xi_J in [0, 1].
int
LadrunoDispBeamColumn3d::solveDualHingeJump3d(const Vector &v, double L)
{
  if (L <= 0.0)
    return -1;
  double oneOverL = 1.0 / L;

  const double betaI  = 1.0 - hingeXiI;
  const double gammaI = hingeXiI;
  const double betaJ  = 1.0 - hingeXiJ;
  const double gammaJ = hingeXiJ;

  // W is the 6x4 projection matrix -d(v_tilde)/d(alpha)
  double W[6][4] = {{0.0}};
  W[1][0] = betaI;  W[2][0] = gammaI; // alpha_zI -> v1, v2
  W[3][1] = betaI;  W[4][1] = gammaI; // alpha_yI -> v3, v4
  W[1][2] = betaJ;  W[2][2] = gammaJ; // alpha_zJ -> v1, v2
  W[3][3] = betaJ;  W[4][3] = gammaJ; // alpha_yJ -> v3, v4

  bool active[4];
  active[0] = (theHingeZ  != 0 || theHingeC  != 0);
  active[1] = (theHingeY  != 0 || theHingeC  != 0);
  active[2] = (theHingeZJ != 0 || theHingeCJ != 0);
  active[3] = (theHingeYJ != 0 || theHingeCJ != 0);

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  double alpha[4];
  alpha[0] = active[0] ? hingeJumpCommitZ  : 0.0;
  alpha[1] = active[1] ? hingeJumpCommitY  : 0.0;
  alpha[2] = active[2] ? hingeJumpCommitZJ : 0.0;
  alpha[3] = active[3] ? hingeJumpCommitYJ : 0.0;

  LadrunoCohesiveHinge *chZI = dynamic_cast<LadrunoCohesiveHinge*>(theHingeZ);
  LadrunoCohesiveHinge *chYI = dynamic_cast<LadrunoCohesiveHinge*>(theHingeY);
  LadrunoCohesiveHingeBiaxial *chCI = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeC);
  LadrunoCohesiveHinge *chZJ = dynamic_cast<LadrunoCohesiveHinge*>(theHingeZJ);
  LadrunoCohesiveHinge *chYJ = dynamic_cast<LadrunoCohesiveHinge*>(theHingeYJ);
  LadrunoCohesiveHingeBiaxial *chCJ = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeCJ);

  const bool anyImplex = ldbc3dHingeBiaxialIsImplex(theHingeC) ||
                         ldbc3dHingeBiaxialIsImplex(theHingeCJ) ||
                         ldbc3dHingeIsImplex(theHingeZ) ||
                         ldbc3dHingeIsImplex(theHingeY) ||
                         ldbc3dHingeIsImplex(theHingeZJ) ||
                         ldbc3dHingeIsImplex(theHingeYJ);

  const int maxIter = 50;
  double Kva[6][4] = {{0.0}};
  double invKaa[4][4] = {{0.0}};
  double damp = 1.0;
  double hnormPrev = 1.0e300;

  for (int iter = 0; iter < maxIter; iter++) {
    int err = 0;
    double qSec[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double Kbulk[6][6] = {{0.0}};

    double v1t = v(1) - betaI * alpha[0] - betaJ * alpha[2];
    double v2t = v(2) - gammaI * alpha[0] - gammaJ * alpha[2];
    double v3t = v(3) - betaI * alpha[1] - betaJ * alpha[3];
    double v4t = v(4) - gammaI * alpha[1] - gammaJ * alpha[3];

    for (int i = 0; i < numSections; i++) {
      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      double xi6 = 6.0 * xi[i];
      double b1 = xi6 - 4.0;
      double b2 = xi6 - 2.0;

      Vector e(workArea, order);
      for (int j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:  e(j) = oneOverL * v(0); break;
        case SECTION_RESPONSE_MZ: e(j) = oneOverL * (b1 * v1t + b2 * v2t); break;
        case SECTION_RESPONSE_MY: e(j) = oneOverL * (b1 * v3t + b2 * v4t); break;
        case SECTION_RESPONSE_T:  e(j) = oneOverL * v(5); break;
        default:                  e(j) = 0.0; break;
        }
      }

      if (lchMode == 1)       current_section_lch = L;
      else if (lchMode == 2)  current_section_lch = userLch;
      else                    current_section_lch = wt[i] * L;

      err += theSections[i]->setTrialSectionDeformation(e);

      const Vector &s  = theSections[i]->getStressResultant();
      const Matrix &ks = theSections[i]->getSectionTangent();
      double wti = wt[i] * oneOverL;

      Matrix ka(workArea, order, 6);
      ka.Zero();
      for (int j = 0; j < order; j++) {
        double si = s(j) * wt[i];
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          qSec[0] += si;
          for (int k = 0; k < order; k++) ka(k, 0) += ks(k, j) * wti;
          break;
        case SECTION_RESPONSE_MZ:
          qSec[1] += b1 * si;
          qSec[2] += b2 * si;
          for (int k = 0; k < order; k++) {
            double tmp = ks(k, j) * wti;
            ka(k, 1) += b1 * tmp;
            ka(k, 2) += b2 * tmp;
          }
          break;
        case SECTION_RESPONSE_MY:
          qSec[3] += b1 * si;
          qSec[4] += b2 * si;
          for (int k = 0; k < order; k++) {
            double tmp = ks(k, j) * wti;
            ka(k, 3) += b1 * tmp;
            ka(k, 4) += b2 * tmp;
          }
          break;
        case SECTION_RESPONSE_T:
          qSec[5] += si;
          for (int k = 0; k < order; k++) ka(k, 5) += ks(k, j) * wti;
          break;
        default: break;
        }
      }
      for (int j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          for (int k = 0; k < 6; k++) Kbulk[0][k] += ka(j, k);
          break;
        case SECTION_RESPONSE_MZ:
          for (int k = 0; k < 6; k++) {
            double tmp = ka(j, k);
            Kbulk[1][k] += b1 * tmp;
            Kbulk[2][k] += b2 * tmp;
          }
          break;
        case SECTION_RESPONSE_MY:
          for (int k = 0; k < 6; k++) {
            double tmp = ka(j, k);
            Kbulk[3][k] += b1 * tmp;
            Kbulk[4][k] += b2 * tmp;
          }
          break;
        case SECTION_RESPONSE_T:
          for (int k = 0; k < 6; k++) Kbulk[5][k] += ka(j, k);
          break;
        default: break;
        }
      }
    }

    if (err != 0) {
      LadrunoOutputGuard og;
      opserr << "LadrunoDispBeamColumn3d::solveDualHingeJump3d() - failed setTrialSectionDeformations()\n";
      return err;
    }

    for (int m = 0; m < 5; m++)
      qSec[m] += q0[m];

    double sumP = qSec[0];
    double Mbulk[4] = {0.0, 0.0, 0.0, 0.0};
    for (int m = 0; m < 4; m++) {
      for (int a = 0; a < 6; a++) {
        Mbulk[m] += W[a][m] * qSec[a];
        double kv = 0.0;
        if (active[m]) {
          for (int b = 0; b < 6; b++)
            kv -= Kbulk[a][b] * W[b][m];
        }
        Kva[a][m] = kv;
      }
    }

    double Kaa[4][4] = {{0.0}};
    for (int m = 0; m < 4; m++) {
      for (int n = 0; n < 4; n++) {
        if (active[m] && active[n]) {
          double sum = 0.0;
          for (int a = 0; a < 6; a++)
            sum -= W[a][m] * Kva[a][n];
          Kaa[m][n] = sum;
        }
      }
    }

    double Mcoh[4] = {0.0, 0.0, 0.0, 0.0};
    double CdiagSum = 0.0;
    // End I cohesive
    if (theHingeC != 0) {
      if (chCI != 0) chCI->setAxialForce(sumP);
      if (!(postRestoreSync && ldbc3dHingeBiaxialIsImplex(theHingeC))) {
        static thread_local Vector jumpI(2);
        jumpI(0) = alpha[0]; jumpI(1) = alpha[1];
        theHingeC->setTrialStrain(jumpI);
      }
      const Vector &MI = theHingeC->getStress();
      const Matrix &KcI = theHingeC->getTangent();
      Mcoh[0] = MI(0); Mcoh[1] = MI(1);
      double c00 = KcI(0, 0), c11 = KcI(1, 1), c01 = 0.5 * (KcI(0, 1) + KcI(1, 0));
      Kaa[0][0] += c00; Kaa[1][1] += c11;
      Kaa[0][1] += c01; Kaa[1][0] += c01;
      CdiagSum += fabs(c00) + fabs(c11);
    } else {
      if (theHingeZ != 0) {
        if (chZI != 0) chZI->setAxialForce(sumP);
        double fPMZI = (chZI == 0 && pmP0 > 0.0) ? ladrunoDisp3dScalePM(sumP, pmP0, pmAlpha) : 1.0;
        if (!(postRestoreSync && ldbc3dHingeIsImplex(theHingeZ)))
          theHingeZ->setTrialStrain(alpha[0]);
        Mcoh[0] = fPMZI * theHingeZ->getStress();
        double c00 = fPMZI * theHingeZ->getTangent();
        Kaa[0][0] += c00;
        CdiagSum += fabs(c00);
      }
      if (theHingeY != 0) {
        if (chYI != 0) chYI->setAxialForce(sumP);
        double fPMYI = (chYI == 0 && pmP0 > 0.0) ? ladrunoDisp3dScalePM(sumP, pmP0, pmAlpha) : 1.0;
        if (!(postRestoreSync && ldbc3dHingeIsImplex(theHingeY)))
          theHingeY->setTrialStrain(alpha[1]);
        Mcoh[1] = fPMYI * theHingeY->getStress();
        double c11 = fPMYI * theHingeY->getTangent();
        Kaa[1][1] += c11;
        CdiagSum += fabs(c11);
      }
    }
    // End J cohesive
    if (theHingeCJ != 0) {
      if (chCJ != 0) chCJ->setAxialForce(sumP);
      if (!(postRestoreSync && ldbc3dHingeBiaxialIsImplex(theHingeCJ))) {
        static thread_local Vector jumpJ(2);
        jumpJ(0) = alpha[2]; jumpJ(1) = alpha[3];
        theHingeCJ->setTrialStrain(jumpJ);
      }
      const Vector &MJ = theHingeCJ->getStress();
      const Matrix &KcJ = theHingeCJ->getTangent();
      Mcoh[2] = MJ(0); Mcoh[3] = MJ(1);
      double c22 = KcJ(0, 0), c33 = KcJ(1, 1), c23 = 0.5 * (KcJ(0, 1) + KcJ(1, 0));
      Kaa[2][2] += c22; Kaa[3][3] += c33;
      Kaa[2][3] += c23; Kaa[3][2] += c23;
      CdiagSum += fabs(c22) + fabs(c33);
    } else {
      if (theHingeZJ != 0) {
        if (chZJ != 0) chZJ->setAxialForce(sumP);
        double fPMZJ = (chZJ == 0 && pmP0 > 0.0) ? ladrunoDisp3dScalePM(sumP, pmP0, pmAlpha) : 1.0;
        if (!(postRestoreSync && ldbc3dHingeIsImplex(theHingeZJ)))
          theHingeZJ->setTrialStrain(alpha[2]);
        Mcoh[2] = fPMZJ * theHingeZJ->getStress();
        double c22 = fPMZJ * theHingeZJ->getTangent();
        Kaa[2][2] += c22;
        CdiagSum += fabs(c22);
      }
      if (theHingeYJ != 0) {
        if (chYJ != 0) chYJ->setAxialForce(sumP);
        double fPMYJ = (chYJ == 0 && pmP0 > 0.0) ? ladrunoDisp3dScalePM(sumP, pmP0, pmAlpha) : 1.0;
        if (!(postRestoreSync && ldbc3dHingeIsImplex(theHingeYJ)))
          theHingeYJ->setTrialStrain(alpha[3]);
        Mcoh[3] = fPMYJ * theHingeYJ->getStress();
        double c33 = fPMYJ * theHingeYJ->getTangent();
        Kaa[3][3] += c33;
        CdiagSum += fabs(c33);
      }
    }

    double h[4] = {0.0, 0.0, 0.0, 0.0};
    for (int m = 0; m < 4; m++)
      if (active[m]) h[m] = -Mbulk[m] + Mcoh[m];

    double bulkDiagSum = 0.0;
    for (int a = 1; a <= 4; a++) bulkDiagSum += fabs(Kbulk[a][a]);
    double floor = 1.0e-8 * (bulkDiagSum + CdiagSum) + 1.0e-300;
    ladrunoFlooredInv4x4(Kaa, active, floor, invKaa);

    double *mscales[4] = {&hingeMscaleZ, &hingeMscaleY, &hingeMscaleZJ, &hingeMscaleYJ};
    bool allOk = true;
    double hnormSq = 0.0;
    for (int m = 0; m < 4; m++) {
      if (active[m]) {
        double Mhere = fabs(Mcoh[m]) + fabs(Mbulk[m]);
        if (Mhere > *mscales[m]) *mscales[m] = Mhere;
        if (fabs(h[m]) > 1.0e-11 * (*mscales[m]) + 1.0e-12)
          allOk = false;
        double hn = h[m] / (*mscales[m] + 1.0e-12);
        hnormSq += hn * hn;
      }
    }

    if (allOk || (postRestoreSync && anyImplex)) {
      hingeJumpZ  = alpha[0]; hingeJumpY  = alpha[1];
      hingeJumpZJ = alpha[2]; hingeJumpYJ = alpha[3];
      for (int a = 0; a < 6; a++)
        for (int m = 0; m < 4; m++)
          hingeDualKv[a][m] = Kva[a][m];
      for (int m = 0; m < 4; m++)
        for (int n = 0; n < 4; n++)
          hingeDualKaaInv[m][n] = invKaa[m][n];
      return 0;
    }

    double hnorm = sqrt(hnormSq);
    if (hnorm > hnormPrev) damp *= 0.5;
    else                   damp = (damp * 1.5 < 1.0) ? damp * 1.5 : 1.0;
    hnormPrev = hnorm;

    for (int m = 0; m < 4; m++) {
      if (active[m]) {
        double dalpha = 0.0;
        for (int n = 0; n < 4; n++)
          dalpha += invKaa[m][n] * h[n];
        alpha[m] -= damp * dalpha;
      }
    }
  }

  hingeJumpZ  = alpha[0]; hingeJumpY  = alpha[1];
  hingeJumpZJ = alpha[2]; hingeJumpYJ = alpha[3];
  for (int a = 0; a < 6; a++)
    for (int m = 0; m < 4; m++)
      hingeDualKv[a][m] = Kva[a][m];
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < 4; n++)
      hingeDualKaaInv[m][n] = invKaa[m][n];
  {
    LadrunoOutputGuard og;
    opserr << "WARNING LadrunoDispBeamColumn3d (tag " << this->getTag()
           << "): dual embedded-hinge inner Newton did not converge in " << maxIter << " iters\n";
  }
  return 0;
}

double
LadrunoDispBeamColumn3d::getCharacteristicLength(void)
{
  // Ladruno (ADR 32) Tier-1: report the localizing IP's tributary length (set in
  // update()) so crack-band materials regularize over the correct band. Mirrors
  // ForceBeamColumn / the 2D sibling.
  if (current_section_lch > 0.0)
    return current_section_lch;
  return Element::getCharacteristicLength();
}

const Matrix&
LadrunoDispBeamColumn3d::getTangentStiff()
{
  static thread_local Matrix kb(6,6);
  
  // Zero for integral
  q.Zero();
  
  double L = crdTransf->getInitialLength();
  if (L <= 0.0) {
    K.Zero();
    return K;
  }

  this->getBasicStiff(kb);

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points to assemble basic force q
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    double xi6 = 6.0*xi[i];

    // Get the section stress resultant
    const Vector &s = theSections[i]->getStressResultant();
        
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));
    double si;
    for (int j = 0; j < order; j++) {
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si;
	break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si;
	break;
      case SECTION_RESPONSE_MY:
	q(3) += (xi6-4.0)*si; q(4) += (xi6-2.0)*si;
	break;
      case SECTION_RESPONSE_T:
	q(5) += si;
	break;
      default:
	break;
      }
    }

  }

  // Ladruno (ADR 32) Stage-1: ½(θz²+θy²) nonlinear geometry. Add the biaxial
  // bowing contribution to the basic force q (kb was already built by getBasicStiff).
  if (nlGeom) {
    const Vector &v = crdTransf->getBasicTrialDisp();
    for (int i = 0; i < numSections; i++) {
      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      double zeta = xi[i];
      double c1nl = 3.0*zeta*zeta - 4.0*zeta + 1.0;
      double c2nl = 3.0*zeta*zeta - 2.0*zeta;
      double thetaz = c1nl*v(1) + c2nl*v(2);
      double thetay = c1nl*v(3) + c2nl*v(4);
      const Vector &s = theSections[i]->getStressResultant();
      for (int j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_MZ) {
          for (int jj = 0; jj < order; jj++)
            if (code(jj) == SECTION_RESPONSE_P) {
              q(1) += c1nl*thetaz*s(jj)*wt[i]*L;
              q(2) += c2nl*thetaz*s(jj)*wt[i]*L;
            }
        } else if (code(j) == SECTION_RESPONSE_MY) {
          for (int jj = 0; jj < order; jj++)
            if (code(jj) == SECTION_RESPONSE_P) {
              q(3) += c1nl*thetay*s(jj)*wt[i]*L;
              q(4) += c2nl*thetay*s(jj)*wt[i]*L;
            }
        }
      }
    }
  }

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  // Transform to global stiffness
  K = crdTransf->getGlobalStiffMatrix(kb, q);
  //   opserr << this->getTag() << " " << K;
  return K;
}

void
LadrunoDispBeamColumn3d::getBasicStiff(Matrix &kb, int initial)
{
  // Zero for integral
  kb.Zero();
  
  double L = crdTransf->getInitialLength();
  if (L <= 0.0)
    return;
  double oneOverL = 1.0/L;

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    
    double xi6 = 6.0*xi[i];
    int j, k;
    double tmp;

    if (nlGeom && !initial) {
      // Ladruno (ADR 32) Stage-1: ½(θz²+θy²) biaxial nonlinear-geometry tangent
      // (verbatim from DispBeamColumnNL3d::getBasicStiff).
      const Vector &v = crdTransf->getBasicTrialDisp();
      double zeta = xi[i];
      double c1 = 3.0*zeta*zeta - 4.0*zeta + 1.0;
      double c2 = 3.0*zeta*zeta - 2.0*zeta;
      double thetaz = c1*v(1) + c2*v(2);
      double thetay = c1*v(3) + c2*v(4);

      const Matrix &ks = theSections[i]->getSectionTangent();
      const Vector &s  = theSections[i]->getStressResultant();
      double wti = wt[i]*oneOverL;

      // geometric (axial-force) term: int_0^L N * C'*C dx (both bending axes)
      for (j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_P) {
          tmp = s(j)*wt[i]*L;
          kb(1,1) += tmp*c1*c1; kb(2,1) += tmp*c2*c1; kb(1,2) += tmp*c1*c2; kb(2,2) += tmp*c2*c2;
          kb(3,3) += tmp*c1*c1; kb(4,3) += tmp*c2*c1; kb(3,4) += tmp*c1*c2; kb(4,4) += tmp*c2*c2;
        }
      }

      Matrix B(order,6);
      Matrix Cz(order,6);
      Matrix Cy(order,6);
      static thread_local Matrix Cz1(1,6);
      static thread_local Matrix Cy1(1,6);
      Cz1.Zero();
      Cy1.Zero();
      for (j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          B(j,0) = 1.0;
          Cz(j,1) = c1; Cz(j,2) = c2; Cz1(0,1) = c1; Cz1(0,2) = c2;
          Cy(j,3) = c1; Cy(j,4) = c2; Cy1(0,3) = c1; Cy1(0,4) = c2; break;
        case SECTION_RESPONSE_T:  B(j,5) = 1.0; break;
        case SECTION_RESPONSE_MZ: B(j,1) = xi6-4.0; B(j,2) = xi6-2.0; break;
        case SECTION_RESPONSE_MY: B(j,3) = xi6-4.0; B(j,4) = xi6-2.0; break;
        default: break;
        }
      }

      kb.addMatrixTripleProduct(1.0, B, ks, wti);
      Matrix kC(order,6);
      kC.addMatrixProduct(0.0, ks, Cz, 1.0);
      kb.addMatrixTransposeProduct(1.0, B, kC, thetaz*wt[i]);
      kC.addMatrixProduct(0.0, ks, Cy, 1.0);
      kb.addMatrixTransposeProduct(1.0, B, kC, thetay*wt[i]);

      Matrix ks1(1,order);
      static thread_local Matrix ksB(1,6);
      for (j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_P) {
          for (int jj = 0; jj < order; jj++) ks1(0,jj) = ks(j,jj);
          ksB.addMatrixProduct(0.0, ks1, B, 1.0);
          kb.addMatrixTransposeProduct(1.0, Cz1, ksB, thetaz*wt[i]);
          ksB.addMatrixProduct(0.0, ks1, Cz, 1.0);
          kb.addMatrixTransposeProduct(1.0, Cz1, ksB, thetaz*thetaz*wt[i]*L);
          ksB.addMatrixProduct(0.0, ks1, Cz, 1.0);
          kb.addMatrixTransposeProduct(1.0, Cy1, ksB, thetaz*thetay*wt[i]*L);
          ksB.addMatrixProduct(0.0, ks1, B, 1.0);
          kb.addMatrixTransposeProduct(1.0, Cy1, ksB, thetay*wt[i]);
          ksB.addMatrixProduct(0.0, ks1, Cy, 1.0);
          kb.addMatrixTransposeProduct(1.0, Cy1, ksB, thetay*thetay*wt[i]*L);
          ksB.addMatrixProduct(0.0, ks1, Cy, 1.0);
          kb.addMatrixTransposeProduct(1.0, Cz1, ksB, thetay*thetaz*wt[i]*L);
        }
      }
      continue;
    }

    // ----- linear basic strain (stock), also used for initial -----
    Matrix ka(workArea, order, 6);
    ka.Zero();

    // Get the section tangent stiffness
    const Matrix &ks = (initial) ? theSections[i]->getInitialTangent() : theSections[i]->getSectionTangent();

    // Perform numerical integration
    double wti = wt[i]*oneOverL;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < order; k++)
	  ka(k,0) += ks(k,j)*wti;
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < order; k++) {
	  tmp = ks(k,j)*wti;
	  ka(k,1) += (xi6-4.0)*tmp;
	  ka(k,2) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_MY:
	for (k = 0; k < order; k++) {
	  tmp = ks(k,j)*wti;
	  ka(k,3) += (xi6-4.0)*tmp;
	  ka(k,4) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (k = 0; k < order; k++)
	  ka(k,5) += ks(k,j)*wti;
	break;
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 6; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 6; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_MY:
	for (k = 0; k < 6; k++) {
	  tmp = ka(j,k);
	  kb(3,k) += (xi6-4.0)*tmp;
	  kb(4,k) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (k = 0; k < 6; k++)
	  kb(5,k) += ka(j,k);
	break;
      default:
	break;
      }
    }

  }

  // Ladruno (ADR 33 / WP-200): statically condense the rotation jump(s) out of the 6x6
  // basic tangent stiffness when !initial.
  if (!initial && hingeOn) {
    if (dualHingeMode) {
      for (int a = 0; a < 6; a++)
        for (int b = 0; b < 6; b++) {
          double sub = 0.0;
          for (int m = 0; m < 4; m++) {
            if (hingeDualKv[a][m] == 0.0) continue;
            for (int n = 0; n < 4; n++)
              sub += hingeDualKv[a][m] * hingeDualKaaInv[m][n] * hingeDualKv[b][n];
          }
          kb(a, b) -= sub;
        }
    } else if (theHingeY != 0 || theHingeC != 0) {
      double iZZ = hingeKaaInvZZ, iZY = hingeKaaInvZY, iYY = hingeKaaInvYY;
      for (int a = 0; a < 6; a++)
        for (int b = 0; b < 6; b++)
          kb(a,b) -= hingeKvZ[a]*(iZZ*hingeKvZ[b] + iZY*hingeKvY[b])
                   + hingeKvY[a]*(iZY*hingeKvZ[b] + iYY*hingeKvY[b]);
    } else if (hingeKaaZ != 0.0) {
      double invKaa = 1.0/hingeKaaZ;
      for (int a = 0; a < 6; a++)
        for (int b = 0; b < 6; b++)
          kb(a,b) -= hingeKvZ[a]*hingeKvZ[b]*invKaa;
    }
  }
}

const Matrix&
LadrunoDispBeamColumn3d::getInitialStiff()
{
  static thread_local Matrix kb(6,6);

  this->getBasicStiff(kb, 1);

  if (hingeOn) {
    double L = crdTransf->getInitialLength();
    if (L > 0.0) {
      if (dualHingeMode) {
        const double betaI  = 1.0 - hingeXiI;
        const double gammaI = hingeXiI;
        const double betaJ  = 1.0 - hingeXiJ;
        const double gammaJ = hingeXiJ;
        double W[6][4] = {{0.0}};
        W[1][0] = betaI;  W[2][0] = gammaI;
        W[3][1] = betaI;  W[4][1] = gammaI;
        W[1][2] = betaJ;  W[2][2] = gammaJ;
        W[3][3] = betaJ;  W[4][3] = gammaJ;
        bool active[4];
        active[0] = (theHingeZ  != 0 || theHingeC  != 0);
        active[1] = (theHingeY  != 0 || theHingeC  != 0);
        active[2] = (theHingeZJ != 0 || theHingeCJ != 0);
        active[3] = (theHingeYJ != 0 || theHingeCJ != 0);
        double Kva0[6][4] = {{0.0}};
        double Kaa0[4][4] = {{0.0}};
        for (int m = 0; m < 4; m++) {
          if (active[m]) {
            for (int a = 0; a < 6; a++) {
              double kv = 0.0;
              for (int b = 0; b < 6; b++)
                kv -= kb(a, b) * W[b][m];
              Kva0[a][m] = kv;
            }
          }
        }
        for (int m = 0; m < 4; m++)
          for (int n = 0; n < 4; n++)
            if (active[m] && active[n]) {
              double sum = 0.0;
              for (int a = 0; a < 6; a++)
                sum -= W[a][m] * Kva0[a][n];
              Kaa0[m][n] = sum;
            }
        double CdiagSum = 0.0;
        bool anyFinite = false;
        if (theHingeC != 0) {
          const Matrix &KcI0 = theHingeC->getInitialTangent();
          double c00 = KcI0(0,0), c11 = KcI0(1,1), c01 = 0.5*(KcI0(0,1)+KcI0(1,0));
          if (std::isfinite(c00) && std::isfinite(c11) && c00 > 0.0 && c11 > 0.0) {
            Kaa0[0][0] += c00; Kaa0[1][1] += c11;
            Kaa0[0][1] += c01; Kaa0[1][0] += c01;
            CdiagSum += fabs(c00) + fabs(c11);
            if (c00 < 1.0e20 || c11 < 1.0e20) anyFinite = true;
          } else { active[0] = active[1] = false; }
        } else {
          if (theHingeZ != 0) {
            double c00 = theHingeZ->getInitialTangent();
            if (std::isfinite(c00) && c00 > 0.0 && c00 < 1.0e20) {
              Kaa0[0][0] += c00; CdiagSum += fabs(c00); anyFinite = true;
            } else { active[0] = false; }
          }
          if (theHingeY != 0) {
            double c11 = theHingeY->getInitialTangent();
            if (std::isfinite(c11) && c11 > 0.0 && c11 < 1.0e20) {
              Kaa0[1][1] += c11; CdiagSum += fabs(c11); anyFinite = true;
            } else { active[1] = false; }
          }
        }
        if (theHingeCJ != 0) {
          const Matrix &KcJ0 = theHingeCJ->getInitialTangent();
          double c22 = KcJ0(0,0), c33 = KcJ0(1,1), c23 = 0.5*(KcJ0(0,1)+KcJ0(1,0));
          if (std::isfinite(c22) && std::isfinite(c33) && c22 > 0.0 && c33 > 0.0) {
            Kaa0[2][2] += c22; Kaa0[3][3] += c33;
            Kaa0[2][3] += c23; Kaa0[3][2] += c23;
            CdiagSum += fabs(c22) + fabs(c33);
            if (c22 < 1.0e20 || c33 < 1.0e20) anyFinite = true;
          } else { active[2] = active[3] = false; }
        } else {
          if (theHingeZJ != 0) {
            double c22 = theHingeZJ->getInitialTangent();
            if (std::isfinite(c22) && c22 > 0.0 && c22 < 1.0e20) {
              Kaa0[2][2] += c22; CdiagSum += fabs(c22); anyFinite = true;
            } else { active[2] = false; }
          }
          if (theHingeYJ != 0) {
            double c33 = theHingeYJ->getInitialTangent();
            if (std::isfinite(c33) && c33 > 0.0 && c33 < 1.0e20) {
              Kaa0[3][3] += c33; CdiagSum += fabs(c33); anyFinite = true;
            } else { active[3] = false; }
          }
        }
        if (anyFinite) {
          double bulkDiagSum = fabs(kb(1,1)) + fabs(kb(2,2)) + fabs(kb(3,3)) + fabs(kb(4,4));
          double floor = 1.0e-8 * (bulkDiagSum + CdiagSum) + 1.0e-300;
          double invKaa0[4][4] = {{0.0}};
          ladrunoFlooredInv4x4(Kaa0, active, floor, invKaa0);
          for (int a = 0; a < 6; a++)
            for (int b = 0; b < 6; b++) {
              double sub = 0.0;
              for (int m = 0; m < 4; m++) {
                if (!active[m] || Kva0[a][m] == 0.0) continue;
                for (int n = 0; n < 4; n++)
                  if (active[n]) sub += Kva0[a][m] * invKaa0[m][n] * Kva0[b][n];
              }
              kb(a, b) -= sub;
            }
        }
      } else {
        double oneOverL = 1.0/L;
        double xi[maxNumSections];
        double wt[maxNumSections];
        beamInt->getSectionLocations(numSections, L, xi);
        beamInt->getSectionWeights(numSections, L, wt);
        double sumEIz0 = 0.0, sumEIy0 = 0.0, sumCpl0 = 0.0;
        double KvZ0[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double KvY0[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (int i = 0; i < numSections; i++) {
          int order = theSections[i]->getOrder();
          const ID &code = theSections[i]->getType();
          const Matrix &ks0 = theSections[i]->getInitialTangent();
          double xi6 = 6.0*xi[i];
          int idxP = -1, idxMZ = -1, idxMY = -1, idxT = -1;
          for (int j = 0; j < order; j++) {
            switch (code(j)) {
            case SECTION_RESPONSE_P:  idxP  = j; break;
            case SECTION_RESPONSE_MZ: idxMZ = j; break;
            case SECTION_RESPONSE_MY: idxMY = j; break;
            case SECTION_RESPONSE_T:  idxT  = j; break;
            default: break;
            }
          }
          if (idxMZ >= 0) {
            double EIz = ks0(idxMZ, idxMZ);
            sumEIz0 += wt[i]*EIz;
            KvZ0[1] += -oneOverL*wt[i]*(xi6-4.0)*EIz;
            KvZ0[2] += -oneOverL*wt[i]*(xi6-2.0)*EIz;
            if (idxP  >= 0) KvZ0[0] += -oneOverL*wt[i]*ks0(idxP,  idxMZ);
            if (idxMY >= 0) {
              KvZ0[3] += -oneOverL*wt[i]*(xi6-4.0)*ks0(idxMY, idxMZ);
              KvZ0[4] += -oneOverL*wt[i]*(xi6-2.0)*ks0(idxMY, idxMZ);
            }
            if (idxT  >= 0) KvZ0[5] += -oneOverL*wt[i]*ks0(idxT,  idxMZ);
          }
          if (idxMY >= 0) {
            double EIy = ks0(idxMY, idxMY);
            sumEIy0 += wt[i]*EIy;
            KvY0[3] += -oneOverL*wt[i]*(xi6-4.0)*EIy;
            KvY0[4] += -oneOverL*wt[i]*(xi6-2.0)*EIy;
            if (idxP  >= 0) KvY0[0] += -oneOverL*wt[i]*ks0(idxP,  idxMY);
            if (idxMZ >= 0) {
              KvY0[1] += -oneOverL*wt[i]*(xi6-4.0)*ks0(idxMZ, idxMY);
              KvY0[2] += -oneOverL*wt[i]*(xi6-2.0)*ks0(idxMZ, idxMY);
            }
            if (idxT  >= 0) KvY0[5] += -oneOverL*wt[i]*ks0(idxT,  idxMY);
          }
          if (idxMZ >= 0 && idxMY >= 0)
            sumCpl0 += wt[i]*0.5*(ks0(idxMZ, idxMY) + ks0(idxMY, idxMZ));
        }
        if (theHingeY != 0 || theHingeC != 0) {
          double Czz0 = 0.0, Cyy0 = 0.0, Czy0 = 0.0;
          if (theHingeC != 0) {
            const Matrix &Kc0 = theHingeC->getInitialTangent();
            Czz0 = Kc0(0,0); Cyy0 = Kc0(1,1);
            Czy0 = 0.5*(Kc0(0,1) + Kc0(1,0));
          } else {
            if (theHingeZ != 0) Czz0 = theHingeZ->getInitialTangent();
            if (theHingeY != 0) Cyy0 = theHingeY->getInitialTangent();
          }
          if (std::isfinite(Czz0) && std::isfinite(Cyy0) && Czz0 > 0.0 && Cyy0 > 0.0 &&
              (Czz0 < 1.0e20 || Cyy0 < 1.0e20)) {
            double Kaa_zz = oneOverL*sumEIz0 + Czz0;
            double Kaa_yy = oneOverL*sumEIy0 + Cyy0;
            double Kaa_zy = oneOverL*sumCpl0 + Czy0;
            double floor = 1.0e-8*(fabs(oneOverL*sumEIz0) + fabs(oneOverL*sumEIy0) + fabs(Czz0) + fabs(Cyy0)) + 1.0e-300;
            double iZZ0 = 0.0, iZY0 = 0.0, iYY0 = 0.0;
            ladrunoFlooredInv2x2(Kaa_zz, Kaa_zy, Kaa_yy, floor, iZZ0, iZY0, iYY0);
            for (int a = 0; a < 6; a++)
              for (int b = 0; b < 6; b++)
                kb(a,b) -= KvZ0[a]*(iZZ0*KvZ0[b] + iZY0*KvY0[b])
                         + KvY0[a]*(iZY0*KvZ0[b] + iYY0*KvY0[b]);
          }
        } else if (theHingeZ != 0) {
          double Kcoh0 = theHingeZ->getInitialTangent();
          if (std::isfinite(Kcoh0) && Kcoh0 > 0.0 && Kcoh0 < 1.0e20) {
            double Kaa0 = oneOverL*sumEIz0 + Kcoh0;
            double KaaFloor = 1.0e-8*(fabs(oneOverL*sumEIz0) + fabs(Kcoh0)) + 1.0e-300;
            if (fabs(Kaa0) < KaaFloor)
              Kaa0 = (Kaa0 < 0.0) ? -KaaFloor : KaaFloor;
            double invKaa0 = 1.0/Kaa0;
            for (int a = 0; a < 6; a++)
              for (int b = 0; b < 6; b++)
                kb(a,b) -= KvZ0[a]*KvZ0[b]*invKaa0;
          }
        }
      }
    }
  }

  if(theDamping) kb *= theDamping->getStiffnessMultiplier();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);

  return K;
}

const Matrix&
LadrunoDispBeamColumn3d::getMass()
{
  K.Zero();
  
  if (rho == 0.0)
    return K;
  
  double L = crdTransf->getInitialLength();
  if (cMass == 0)  {
    // lumped mass matrix
    double m = 0.5*rho*L;
    K(0,0) = K(1,1) = K(2,2) = K(6,6) = K(7,7) = K(8,8) = m;
  } else  {
    // consistent mass matrix
    static thread_local Matrix ml(12,12);
    ml.Zero();
    double m = rho*L/420.0;
    ml(0,0) = ml(6,6) = m*140.0;
    ml(0,6) = ml(6,0) = m*70.0;
    //ml(3,3) = ml(9,9) = m*(Jx/A)*140.0;  // CURRENTLY NO TORSIONAL MASS 
    //ml(3,9) = ml(9,3) = m*(Jx/A)*70.0;   // CURRENTLY NO TORSIONAL MASS
    
    ml(2,2) = ml(8,8) = m*156.0;
    ml(2,8) = ml(8,2) = m*54.0;
    ml(4,4) = ml(10,10) = m*4.0*L*L;
    ml(4,10) = ml(10,4) = -m*3.0*L*L;
    ml(2,4) = ml(4,2) = -m*22.0*L;
    ml(8,10) = ml(10,8) = -ml(2,4);
    ml(2,10) = ml(10,2) = m*13.0*L;
    ml(4,8) = ml(8,4) = -ml(2,10);
    
    ml(1,1) = ml(7,7) = m*156.0;
    ml(1,7) = ml(7,1) = m*54.0;
    ml(5,5) = ml(11,11) = m*4.0*L*L;
    ml(5,11) = ml(11,5) = -m*3.0*L*L;
    ml(1,5) = ml(5,1) = m*22.0*L;
    ml(7,11) = ml(11,7) = -ml(1,5);
    ml(1,11) = ml(11,1) = -m*13.0*L;
    ml(5,7) = ml(7,5) = -ml(1,11);
    
    // transform local mass matrix to global system
    K = crdTransf->getGlobalMatrixFromLocal(ml);
  }
  
  return K;
}

void
LadrunoDispBeamColumn3d::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  postRestoreSync = false;

  return;
}

int 
LadrunoDispBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  postRestoreSync = false;
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = crdTransf->getInitialLength();

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wz = data(1)*loadFactor;  // Transverse
    double wx = data(2)*loadFactor;  // Axial (+ve from node I to J)

    double Vy = 0.5*wy*L;
    double Mz = Vy*L/6.0; // wy*L*L/12
    double Vz = 0.5*wz*L;
    double My = Vz*L/6.0; // wz*L*L/12
    double P = wx*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= Vy;
    p0[2] -= Vy;
    p0[3] -= Vz;
    p0[4] -= Vz;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= Mz;
    q0[2] += Mz;
    q0[3] += My;
    q0[4] -= My;
  }
  else if (type == LOAD_TAG_BeamUniformMoment) {
    double my = data(1)*loadFactor;  // About y
    double mz = data(2)*loadFactor;  // About z

    // Reactions in basic system
    p0[1] += mz;
    p0[2] -= mz;
    p0[3] -= my;
    p0[4] += my;
  }
  else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
    // Ladruno: trapezoidal 3D partial load support
    double wya = data(0) * loadFactor;  // Transverse Y at start
    double wza = data(1) * loadFactor;  // Transverse Z at start
    double waa = data(2) * loadFactor;  // Axial at start
    double a   = data(3) * L;
    double b   = data(4) * L;
    double wyb = data(5) * loadFactor;  // Transverse Y at end
    double wzb = data(6) * loadFactor;  // Transverse Z at end
    double wab = data(7) * loadFactor;  // Axial at end

    double ell = b - a;
    if (ell <= 0.0 || a < 0.0 || b > L)
      return 0;

    // Axial
    double P  = 0.5 * ell * (waa + wab);
    double PJ = (a * P + (ell * ell / 6.0) * (waa + 2.0 * wab)) / L;

    // Transverse Y (bending about z)
    double J0y = 0.5 * ell * (wya + wyb);
    double J1y = (ell * ell / 6.0) * (wya + 2.0 * wyb);
    double J2y = (ell * ell * ell / 12.0) * (wya + 3.0 * wyb);
    double J3y = (ell * ell * ell * ell / 20.0) * (wya + 4.0 * wyb);
    double I0y = J0y;
    double I1y = a * J0y + J1y;
    double I2y = a * a * J0y + 2.0 * a * J1y + J2y;
    double I3y = a * a * a * J0y + 3.0 * a * a * J1y + 3.0 * a * J2y + J3y;
    double Vy2 = I1y / L;
    double Vy1 = I0y - Vy2;
    double M1y = -I1y + (2.0 / L) * I2y - (1.0 / (L * L)) * I3y;
    double M2y = (1.0 / L) * I2y - (1.0 / (L * L)) * I3y;

    // Transverse Z (bending about y)
    double J0z = 0.5 * ell * (wza + wzb);
    double J1z = (ell * ell / 6.0) * (wza + 2.0 * wzb);
    double J2z = (ell * ell * ell / 12.0) * (wza + 3.0 * wzb);
    double J3z = (ell * ell * ell * ell / 20.0) * (wza + 4.0 * wzb);
    double I0z = J0z;
    double I1z = a * J0z + J1z;
    double I2z = a * a * J0z + 2.0 * a * J1z + J2z;
    double I3z = a * a * a * J0z + 3.0 * a * a * J1z + 3.0 * a * J2z + J3z;
    double Vz2 = I1z / L;
    double Vz1 = I0z - Vz2;
    double M1z = -I1z + (2.0 / L) * I2z - (1.0 / (L * L)) * I3z;
    double M2z = (1.0 / L) * I2z - (1.0 / (L * L)) * I3z;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= Vy1;
    p0[2] -= Vy2;
    p0[3] -= Vz1;
    p0[4] -= Vz2;

    // Fixed end forces in basic system
    q0[0] -= PJ;
    q0[1] += M1y;
    q0[2] += M2y;
    q0[3] -= M1z;
    q0[4] -= M2z;
  }
  else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0)*loadFactor;
    double Pz = data(1)*loadFactor;
    double N  = data(2)*loadFactor;
    double aOverL = data(3);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1, V2;
    V1 = Py*(1.0-aOverL);
    V2 = Py*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;
    V1 = Pz*(1.0-aOverL);
    V2 = Pz*aOverL;
    p0[3] -= V1;
    p0[4] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1, M2;
    M1 = -a * b2 * Py * L2;
    M2 = a2 * b * Py * L2;
    q0[1] += M1;
    q0[2] += M2;
    M1 = -a * b2 * Pz * L2;
    M2 = a2 * b * Pz * L2;
    q0[3] -= M1;
    q0[4] -= M2;
  }
  else {
    opserr << "LadrunoDispBeamColumn3d::addLoad() -- load type unknown for element with tag: " << 
      this->getTag() << endln;
    return -1;
  }

  return 0;
}

int 
LadrunoDispBeamColumn3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0) 
    return 0;
  
  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "LadrunoDispBeamColumn3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }
  
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    Q(0) -= m*Raccel1(0);
    Q(1) -= m*Raccel1(1);
    Q(2) -= m*Raccel1(2);
    Q(6) -= m*Raccel2(0);
    Q(7) -= m*Raccel2(1);
    Q(8) -= m*Raccel2(2);

  } else  {
    // use matrix vector multip. for consistent mass matrix
    static thread_local Vector Raccel(12);
    for (int i=0; i<6; i++)  {
      Raccel(i)   = Raccel1(i);
      Raccel(i+6) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }
  
  return 0;
}

const Vector&
LadrunoDispBeamColumn3d::getResistingForce()
{
  double L = crdTransf->getInitialLength();

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Zero for integration
  q.Zero();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    double xi6 = 6.0*xi[i];

    // Get section stress resultant
    const Vector &s = theSections[i]->getStressResultant();

    // Ladruno (ADR 32) Stage-1: biaxial slopes for the ½(θz²+θy²) bowing force
    const Vector &vb = crdTransf->getBasicTrialDisp();
    double zeta = xi[i];
    double c1nl = 3.0*zeta*zeta - 4.0*zeta + 1.0;
    double c2nl = 3.0*zeta*zeta - 2.0*zeta;
    double thetaz = c1nl*vb(1) + c2nl*vb(2);
    double thetay = c1nl*vb(3) + c2nl*vb(4);

    // Perform numerical integration on internal force
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));

    double si;
    for (int j = 0; j < order; j++) {
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si;
	break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si;
	if (nlGeom) {
	  for (int jj = 0; jj < order; jj++)
	    if (code(jj) == SECTION_RESPONSE_P) {
	      q(1) += c1nl*thetaz*s(jj)*wt[i]*L;
	      q(2) += c2nl*thetaz*s(jj)*wt[i]*L;
	    }
	}
	break;
      case SECTION_RESPONSE_MY:
	q(3) += (xi6-4.0)*si; q(4) += (xi6-2.0)*si;
	if (nlGeom) {
	  for (int jj = 0; jj < order; jj++)
	    if (code(jj) == SECTION_RESPONSE_P) {
	      q(3) += c1nl*thetay*s(jj)*wt[i]*L;
	      q(4) += c2nl*thetay*s(jj)*wt[i]*L;
	    }
	}
	break;
      case SECTION_RESPONSE_T:
	q(5) += si;
	break;
      default:
	break;
      }
    }

  }
  
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  if (theDamping) theDamping->update(q);

  // Transform forces
  Vector p0Vec(p0, 5);
  P = crdTransf->getGlobalResistingForce(q, p0Vec);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (rho != 0)
    P.addVector(1.0, Q, -1.0);
  
  return P;
}

const Vector &
LadrunoDispBeamColumn3d::getDampingForce(void)
{
  crdTransf->update();

  return crdTransf->getGlobalResistingForce(theDamping->getDampingForce(), Vector(5));
}

const Vector&
LadrunoDispBeamColumn3d::getResistingForceIncInertia()
{
  // SNAPSHOT into a function-local buffer, not the shared static P: betaK
  // Rayleigh re-enters getTangentStiff(), so accumulating into a buffer that
  // path could refill would silently drop inertia (LEDGER_quirks "MUST snapshot
  // the shared static `resid`"). Returning `res` also keeps the result distinct
  // from the P that getResistingForce() returns. Same operations, same order,
  // so results are bit-identical (WP-116).  // Ladruno
  static thread_local Vector res(12);
  res = this->getResistingForce();
  
  if (theDamping) res += this->getDampingForce();
  
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
  
    res(0) += m*accel1(0);
    res(1) += m*accel1(1);
    res(2) += m*accel1(2);
    res(6) += m*accel2(0);
    res(7) += m*accel2(1);
    res(8) += m*accel2(2);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static thread_local Vector accel(12);
    for (int i=0; i<6; i++)  {
      accel(i)   = accel1(i);
      accel(i+6) = accel2(i);
    }
    res.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }
    
    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      res.addVector(1.0, this->getRayleighDampingForces(), 1.0);

  } else {

    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      res.addVector(1.0, this->getRayleighDampingForces(), 1.0);
  }
  
  return res;
}

int
LadrunoDispBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;
  
  static thread_local Vector data(57);  // Ladruno (WP-200: +17 dual-hinge + 2 PM + 10 q0/p0)
  data.Zero();
  data(0) = this->getTag();
  data(1) = connectedExternalNodes(0);
  data(2) = connectedExternalNodes(1);
  data(3) = numSections;
  data(4) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
      crdTransf->setDbTag(crdTransfDbTag);
  }
  data(5) = crdTransfDbTag;
  data(6) = beamInt->getClassTag();
  int beamIntDbTag  = beamInt->getDbTag();
  if (beamIntDbTag  == 0) {
    beamIntDbTag = theChannel.getDbTag();
    if (beamIntDbTag  != 0) 
      beamInt->setDbTag(beamIntDbTag);
  }
  data(7) = beamIntDbTag;
  data(8) = rho;
  data(9) = cMass;
  data(10) = alphaM;
  data(11) = betaK;
  data(12) = betaK0;
  data(13) = betaKc;
  
  data(14) = 0;
  data(15) = 0;
  if (theDamping) {
    data(14) = theDamping->getClassTag();
    int dmpDbTag = theDamping->getDbTag();
    if (dmpDbTag == 0) {
      dmpDbTag = theChannel.getDbTag();
      if (dmpDbTag != 0)
	      theDamping->setDbTag(dmpDbTag);
	  }
    data(15) = dmpDbTag;
  }

  // Ladruno (ADR 32): regularization-length config + nonlinear-geometry flag
  data(16) = lchMode;
  data(17) = userLch;
  data(18) = nlGeom;

  // Ladruno (ADR 33) Tier-2: embedded-hinge flag + committed strong-axis jump
  data(19) = hingeOn;
  data(20) = hingeJumpCommitZ;
  // Ladruno (ADR 33 PR-3b): weak-axis (biaxial) hinge presence flag + committed jump
  data(21) = (theHingeY != 0) ? 1 : 0;
  data(22) = hingeJumpCommitY;
  // Ladruno (ADR 34): coupled biaxial cohesive (NDMaterial) presence + classTag/dbTag
  data(23) = (theHingeC != 0) ? 1 : 0;
  data(24) = 0;
  data(25) = 0;
  if (theHingeC != 0) {
    data(24) = theHingeC->getClassTag();
    int cDbTag = theHingeC->getDbTag();
    if (cDbTag == 0) {
      cDbTag = theChannel.getDbTag();
      if (cDbTag != 0) theHingeC->setDbTag(cDbTag);
    }
    data(25) = theHingeC->getDbTag();
  }
  data(26) = hingeMscaleZ;
  data(27) = hingeMscaleY;

  // Ladruno (WP-200): dual-hinge fields (folded into data Vector to avoid FE_Datastore ID size collisions)
  data(28) = dualHingeMode;
  data(29) = hingeXiI;
  data(30) = hingeXiJ;
  data(31) = hingeJumpCommitZJ;
  data(32) = hingeJumpCommitYJ;
  data(33) = hingeMscaleZJ;
  data(34) = hingeMscaleYJ;
  if (theHingeZ != 0) {
    data(35) = theHingeZ->getClassTag();
    int hDb = theHingeZ->getDbTag();
    if (hDb == 0) { hDb = theChannel.getDbTag(); if (hDb != 0) theHingeZ->setDbTag(hDb); }
    data(36) = theHingeZ->getDbTag();
  }
  if (theHingeY != 0) {
    data(37) = theHingeY->getClassTag();
    int hDb = theHingeY->getDbTag();
    if (hDb == 0) { hDb = theChannel.getDbTag(); if (hDb != 0) theHingeY->setDbTag(hDb); }
    data(38) = theHingeY->getDbTag();
  }
  if (theHingeZJ != 0) {
    data(39) = theHingeZJ->getClassTag();
    int hDb = theHingeZJ->getDbTag();
    if (hDb == 0) { hDb = theChannel.getDbTag(); if (hDb != 0) theHingeZJ->setDbTag(hDb); }
    data(40) = theHingeZJ->getDbTag();
  }
  if (theHingeYJ != 0) {
    data(41) = theHingeYJ->getClassTag();
    int hDb = theHingeYJ->getDbTag();
    if (hDb == 0) { hDb = theChannel.getDbTag(); if (hDb != 0) theHingeYJ->setDbTag(hDb); }
    data(42) = theHingeYJ->getDbTag();
  }
  if (theHingeCJ != 0) {
    data(43) = theHingeCJ->getClassTag();
    int hDb = theHingeCJ->getDbTag();
    if (hDb == 0) { hDb = theChannel.getDbTag(); if (hDb != 0) theHingeCJ->setDbTag(hDb); }
    data(44) = theHingeCJ->getDbTag();
  }
  data(45) = pmP0;
  data(46) = pmAlpha;
  for (int m = 0; m < 5; m++) {
    data(47 + m) = q0[m];
    data(52 + m) = p0[m];
  }

  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send data Vector\n";
     return -1;
  }
  
  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send crdTranf\n";
     return -1;
  }      

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send beamInt\n";
    return -1;
  }      
  
  ID idSections(2*numSections);
  loc = 0;
  for (i = 0; i<numSections; i++) {
    int sectClassTag = theSections[i]->getClassTag();
    int sectDbTag = theSections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      theSections[i]->setDbTag(sectDbTag);
    }

    idSections(loc) = sectClassTag;
    idSections(loc+1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0)  {
    opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - section " << j << "failed to send itself\n";
      return -1;
    }
  }

  // Ask the Damping to send itself
  if (theDamping && theDamping->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf -- could not send Damping\n";
      return -1;
  }

  if (hingeOn) {
    if (theHingeZ != 0 && theHingeZ->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send hinge material\n";
      return -1;
    }
    if (theHingeY != 0 && theHingeY->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send weak-axis hinge material\n";
      return -1;
    }
    if (theHingeC != 0 && theHingeC->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send coupled biaxial hinge material\n";
      return -1;
    }
    if (theHingeZJ != 0 && theHingeZJ->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send end-J hinge material\n";
      return -1;
    }
    if (theHingeYJ != 0 && theHingeYJ->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send end-J weak-axis hinge material\n";
      return -1;
    }
    if (theHingeCJ != 0 && theHingeCJ->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send end-J coupled biaxial hinge material\n";
      return -1;
    }
  }

  return 0;
}

int
LadrunoDispBeamColumn3d::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  int i;
  
  static thread_local Vector data(57);

  if (theChannel.recvVector(dbTag, commitTag, data) < 0)  {
    opserr << "LadrunoDispBeamColumn3d::recvSelf() - failed to recv data Vector\n";
    return -1;
  }
  
  this->setTag((int)data(0));
  connectedExternalNodes(0) = (int)data(1);
  connectedExternalNodes(1) = (int)data(2);
  int nSect = (int)data(3);
  if (nSect <= 0 || nSect > maxNumSections) {
    opserr << "LadrunoDispBeamColumn3d::recvSelf() - invalid nSect (" << nSect
            << "); must be in [1, " << maxNumSections << "]\n";
    return -1;
  }
  int crdTransfClassTag = (int)data(4);
  int crdTransfDbTag = (int)data(5);
  
  int beamIntClassTag = (int)data(6);
  int beamIntDbTag = (int)data(7);
  
  rho = data(8);
  cMass = (int)data(9);
  
  alphaM = data(10);
  betaK = data(11);
  betaK0 = data(12);
  betaKc = data(13);

  lchMode = (int)data(16);
  userLch = data(17);
  nlGeom = (int)data(18);
  current_section_lch = 0.0;

  hingeOn = (int)data(19);
  hingeJumpCommitZ = data(20);
  hingeJumpZ = hingeJumpCommitZ;
  hingeKaaZ = 0.0;
  hingeKaaCommitZ = 0.0;
  hingeKvZ[0] = hingeKvZ[1] = hingeKvZ[2] = hingeKvZ[3] = hingeKvZ[4] = hingeKvZ[5] = 0.0;
  hingeKvCommitZ[0] = hingeKvCommitZ[1] = hingeKvCommitZ[2] = hingeKvCommitZ[3] = hingeKvCommitZ[4] = hingeKvCommitZ[5] = 0.0;
  hingeMscaleZ = data(26);
  hingeMscaleCommitZ = hingeMscaleZ;

  int hasHingeY = (int)data(21);
  hingeJumpCommitY = data(22);
  hingeJumpY = hingeJumpCommitY;
  hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
  hingeKvCommitY[0] = hingeKvCommitY[1] = hingeKvCommitY[2] = hingeKvCommitY[3] = hingeKvCommitY[4] = hingeKvCommitY[5] = 0.0;
  hingeMscaleY = data(27);
  hingeMscaleCommitY = hingeMscaleY;
  hingeKaaInvZZ = hingeKaaInvZY = hingeKaaInvYY = 0.0;
  hingeKaaInvCommitZZ = hingeKaaInvCommitZY = hingeKaaInvCommitYY = 0.0;

  int hasHingeC = (int)data(23);
  int cClassTag = (int)data(24);
  int cDbTag    = (int)data(25);

  dualHingeMode = (int)data(28);
  hingeXiI = data(29);
  hingeXiJ = data(30);
  hingeJumpCommitZJ = data(31);
  hingeJumpZJ = hingeJumpCommitZJ;
  hingeJumpCommitYJ = data(32);
  hingeJumpYJ = hingeJumpCommitYJ;
  hingeMscaleZJ = data(33);
  hingeMscaleCommitZJ = hingeMscaleZJ;
  hingeMscaleYJ = data(34);
  hingeMscaleCommitYJ = hingeMscaleYJ;
  for (int a = 0; a < 6; a++)
    for (int m = 0; m < 4; m++)
      hingeDualKv[a][m] = hingeDualKvCommit[a][m] = 0.0;
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < 4; n++)
      hingeDualKaaInv[m][n] = hingeDualKaaInvCommit[m][n] = 0.0;

  int zIClassTag = (int)data(35), zIDbTag = (int)data(36);
  int yIClassTag = (int)data(37), yIDbTag = (int)data(38);
  int zJClassTag = (int)data(39), zJDbTag = (int)data(40);
  int yJClassTag = (int)data(41), yJDbTag = (int)data(42);
  int cJClassTag = (int)data(43), cJDbTag = (int)data(44);
  pmP0 = data(45);
  pmAlpha = data(46);
  for (int m = 0; m < 5; m++) {
    q0[m] = data(47 + m);
    p0[m] = data(52 + m);
  }

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "LadrunoDispBeamColumn3d::recvSelf() - " <<
	  "failed to obtain a CrdTrans object with classTag" <<
	  crdTransfClassTag << endln;
	return -2;	  
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }      

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
      if (beamInt != 0)
	  delete beamInt;

      beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

      if (beamInt == 0) {
	opserr << "LadrunoDispBeamColumn3d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	  beamIntClassTag << endln;
	exit(-1);
      }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to recv beam integration\n";
     return -3;
  }      
  
  ID idSections(2*nSect);
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "LadrunoDispBeamColumn3d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  if (numSections != nSect) {

    if (numSections != 0) {
      for (int i=0; i<numSections; i++)
	delete theSections[i];
      delete [] theSections;
    }

    theSections = new SectionForceDeformation *[nSect];
    if (theSections == 0) {
      opserr << "LadrunoDispBeamColumn3d::recvSelf() - out of memory creating sections array of size" <<
	nSect << endln;
      exit(-1);
    }    

    numSections = nSect;
    loc = 0;
    
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      theSections[i] = theBroker.getNewSection(sectClassTag);
      if (theSections[i] == 0) {
	opserr << "LadrunoDispBeamColumn3d::recvSelf() - Broker could not create Section of class type" <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "LadrunoDispBeamColumn3d::recvSelf() - section " <<
	  i << "failed to recv itself\n";
	return -1;
      }     
    }

  } else {

    loc = 0;
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;

      if (theSections[i]->getClassTag() !=  sectClassTag) {
	delete theSections[i];
	theSections[i] = theBroker.getNewSection(sectClassTag);
	if (theSections[i] == 0) {
	  opserr << "LadrunoDispBeamColumn3d::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  exit(-1);
	}
      }

      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "LadrunoDispBeamColumn3d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
	return -1;
      }     
    }
  }


  // Check if the Damping is null; if so, get a new one
  int dmpTag = (int)data(14);
  if (dmpTag) {
    if (theDamping == 0) {
      theDamping = theBroker.getNewDamping(dmpTag);
      if (theDamping == 0) {
        opserr << "LadrunoDispBeamColumn3d::recvSelf -- could not get a Damping\n";
        exit(-1);
      }
    }
  
    if (theDamping->getClassTag() != dmpTag) {
      delete theDamping;
      theDamping = theBroker.getNewDamping(dmpTag);
      if (theDamping == 0) {
        opserr << "LadrunoDispBeamColumn3d::recvSelf -- could not get a Damping\n";
        exit(-1);
      }
    }
  
    theDamping->setDbTag((int)data(15));
    if (theDamping->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "LadrunoDispBeamColumn3d::recvSelf -- could not receive Damping\n";
      exit(-1);
    }
  }
  else {
    if (theDamping) {
      delete theDamping;
      theDamping = 0;
    }
  }

  if (hingeOn) {
    if (zIClassTag != 0) {
      if (theHingeZ == 0 || theHingeZ->getClassTag() != zIClassTag) {
        if (theHingeZ) delete theHingeZ;
        theHingeZ = theBroker.getNewUniaxialMaterial(zIClassTag);
        if (theHingeZ == 0) return -1;
      }
      theHingeZ->setDbTag(zIDbTag);
      if (theHingeZ->recvSelf(commitTag, theChannel, theBroker) < 0) return -1;
    } else if (theHingeZ) { delete theHingeZ; theHingeZ = 0; }

    if (yIClassTag != 0) {
      if (theHingeY == 0 || theHingeY->getClassTag() != yIClassTag) {
        if (theHingeY) delete theHingeY;
        theHingeY = theBroker.getNewUniaxialMaterial(yIClassTag);
        if (theHingeY == 0) return -1;
      }
      theHingeY->setDbTag(yIDbTag);
      if (theHingeY->recvSelf(commitTag, theChannel, theBroker) < 0) return -1;
    } else if (theHingeY) { delete theHingeY; theHingeY = 0; }

    if (hasHingeC && cClassTag != 0) {
      if (theHingeC == 0 || theHingeC->getClassTag() != cClassTag) {
        if (theHingeC) delete theHingeC;
        theHingeC = theBroker.getNewNDMaterial(cClassTag);
        if (theHingeC == 0) return -1;
      }
      theHingeC->setDbTag(cDbTag);
      if (theHingeC->recvSelf(commitTag, theChannel, theBroker) < 0) return -1;
    } else if (theHingeC) { delete theHingeC; theHingeC = 0; }

    if (zJClassTag != 0) {
      if (theHingeZJ == 0 || theHingeZJ->getClassTag() != zJClassTag) {
        if (theHingeZJ) delete theHingeZJ;
        theHingeZJ = theBroker.getNewUniaxialMaterial(zJClassTag);
        if (theHingeZJ == 0) return -1;
      }
      theHingeZJ->setDbTag(zJDbTag);
      if (theHingeZJ->recvSelf(commitTag, theChannel, theBroker) < 0) return -1;
    } else if (theHingeZJ) { delete theHingeZJ; theHingeZJ = 0; }

    if (yJClassTag != 0) {
      if (theHingeYJ == 0 || theHingeYJ->getClassTag() != yJClassTag) {
        if (theHingeYJ) delete theHingeYJ;
        theHingeYJ = theBroker.getNewUniaxialMaterial(yJClassTag);
        if (theHingeYJ == 0) return -1;
      }
      theHingeYJ->setDbTag(yJDbTag);
      if (theHingeYJ->recvSelf(commitTag, theChannel, theBroker) < 0) return -1;
    } else if (theHingeYJ) { delete theHingeYJ; theHingeYJ = 0; }

    if (cJClassTag != 0) {
      if (theHingeCJ == 0 || theHingeCJ->getClassTag() != cJClassTag) {
        if (theHingeCJ) delete theHingeCJ;
        theHingeCJ = theBroker.getNewNDMaterial(cJClassTag);
        if (theHingeCJ == 0) return -1;
      }
      theHingeCJ->setDbTag(cJDbTag);
      if (theHingeCJ->recvSelf(commitTag, theChannel, theBroker) < 0) return -1;
    } else if (theHingeCJ) { delete theHingeCJ; theHingeCJ = 0; }
  }
  else {
    if (theHingeZ)  { delete theHingeZ;  theHingeZ  = 0; }
    if (theHingeY)  { delete theHingeY;  theHingeY  = 0; }
    if (theHingeC)  { delete theHingeC;  theHingeC  = 0; }
    if (theHingeZJ) { delete theHingeZJ; theHingeZJ = 0; }
    if (theHingeYJ) { delete theHingeYJ; theHingeYJ = 0; }
    if (theHingeCJ) { delete theHingeCJ; theHingeCJ = 0; }
  }

  postRestoreSync = true;
  vRestoreValid = false;
  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
    this->update();
    this->getResistingForce();
  }

  return 0;
}

void
LadrunoDispBeamColumn3d::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "\nLadrunoDispBeamColumn3d, element id:  " << this->getTag() << endln;
		s << "\tConnected external nodes:  " << connectedExternalNodes;
		s << "\tCoordTransf: " << crdTransf->getTag() << endln;
		s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;

		double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
		double L = crdTransf->getInitialLength();
		double oneOverL = 1.0 / L;

		N = q(0);
		Mz1 = q(1);
		Mz2 = q(2);
		Vy = (Mz1 + Mz2)*oneOverL;
		My1 = q(3);
		My2 = q(4);
		Vz = -(My1 + My2)*oneOverL;
		T = q(5);

		s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
			<< -N + p0[0] << ' ' << Mz1 << ' ' << Vy + p0[1] << ' ' << My1 << ' ' << Vz + p0[3] << ' ' << -T << endln;
		s << "\tEnd 2 Forces (P Mz Vy My Vz T): "
			<< N << ' ' << Mz2 << ' ' << -Vy + p0[2] << ' ' << My2 << ' ' << -Vz + p0[4] << ' ' << T << endln;
		s << "Number of sections: " << numSections << endln;
		beamInt->Print(s, flag);

		for (int i = 0; i < numSections; i++) {
		  //opserr << "Section Type: " << theSections[i]->getClassTag() << endln;
		  theSections[i]->Print(s,flag);
		}
		//  if (rho != 0)
		//    opserr << "Mass: \n" << this->getMass();
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"LadrunoDispBeamColumn3d\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
		s << "\"sections\": [";
		for (int i = 0; i < numSections - 1; i++)
			s << "\"" << theSections[i]->getTag() << "\", ";
		s << "\"" << theSections[numSections - 1]->getTag() << "\"], ";
		s << "\"integration\": ";
		beamInt->Print(s, flag);
		s << ", \"massperlength\": " << rho << ", ";
		s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
	}
}


int
LadrunoDispBeamColumn3d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numModes)
{
    static thread_local Vector v1(3);
    static thread_local Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response*
LadrunoDispBeamColumn3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","LadrunoDispBeamColumn3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    // global force - 
    if (LadrunoResp::isGlobalForce(argv[0])) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");


      theResponse = new ElementResponse(this, 1, P);

    // local force -
    }  else if (LadrunoResp::isLocalForce(argv[0])) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Vy_2");
      output.tag("ResponseType","Vz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");

      theResponse = new ElementResponse(this, 2, P);
    }
    else if (LadrunoResp::isBasicForce(argv[0])) {
      output.tag("ResponseType","N");
      output.tag("ResponseType","M1");
      output.tag("ResponseType","M2");

      theResponse = new ElementResponse(this, 9, Vector(6));
    }
    else if (strcmp(argv[0],"basicStiffness") == 0) {
      output.tag("ResponseType","N");
      output.tag("ResponseType","M1");
      output.tag("ResponseType","M2");

      theResponse = new ElementResponse(this, 19, Matrix(6,6));
    // global damping force - 
    } else if (theDamping && (strcmp(argv[0],"globalDampingForce") == 0 || strcmp(argv[0],"globalDampingForces") == 0)) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");


      theResponse = new ElementResponse(this, 21, P);

    // local damping force -
    } else if (theDamping && (strcmp(argv[0],"localDampingForce") == 0 || strcmp(argv[0],"localDampingForces") == 0)) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Vy_2");
      output.tag("ResponseType","Vz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");

      theResponse = new ElementResponse(this, 22, P);

    } else if (theDamping && (strcmp(argv[0],"basicDampingForce") == 0 || strcmp(argv[0],"basicDampingForces") == 0)) {

      output.tag("ResponseType","N");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Mz_2");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","T");
    
      theResponse = new ElementResponse(this, 23, Vector(6));

    // chord rotation -
    }  else if (LadrunoResp::isBasicDef(argv[0])) {

      output.tag("ResponseType","eps");
      output.tag("ResponseType","thetaZ_1");
      output.tag("ResponseType","thetaZ_2");
      output.tag("ResponseType","thetaY_1");
      output.tag("ResponseType","thetaY_2");
      output.tag("ResponseType","thetaX");

      theResponse = new ElementResponse(this, 3, Vector(6));

    // plastic rotation -
    } else if (LadrunoResp::isPlasticDef(argv[0])) {

    output.tag("ResponseType","epsP");
    output.tag("ResponseType","thetaZP_1");
    output.tag("ResponseType","thetaZP_2");
    output.tag("ResponseType","thetaYP_1");
    output.tag("ResponseType","thetaYP_2");
    output.tag("ResponseType","thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(6));
  

  } else if (LadrunoResp::isRayleighForce(argv[0])) {

    theResponse =  new ElementResponse(this, 12, P);

  }   
    else if (strcmp(argv[0],"integrationPoints") == 0)
      theResponse = new ElementResponse(this, 10, Vector(numSections));

    else if (strcmp(argv[0],"integrationWeights") == 0)
      theResponse = new ElementResponse(this, 11, Vector(numSections));

    else if (strcmp(argv[0],"sectionTags") == 0)
      theResponse = new ElementResponse(this, 110, ID(numSections));

  // section response -
  else if (strcmp(argv[0],"sectionX") == 0) {
      if (argc > 2) {
	float sectionLoc = atof(argv[1]);
	
	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamInt->getSectionLocations(numSections, L, xi);
	
	sectionLoc /= L;
	
	float minDistance = fabs(xi[0]-sectionLoc);
	int sectionNum = 0;
	for (int i = 1; i < numSections; i++) {
	  if (fabs(xi[i]-sectionLoc) < minDistance) {
	    minDistance = fabs(xi[i]-sectionLoc);
	    sectionNum = i;
	  }
	}
	
	output.tag("GaussPointOutput");
	output.attr("number",sectionNum+1);
	output.attr("eta",xi[sectionNum]*L);
	
	theResponse = theSections[sectionNum]->setResponse(&argv[2], argc-2, output);
      }
    }
    
    else if (strcmp(argv[0],"section") == 0) { 
      if (argc > 1) {
	
	int sectionNum = atoi(argv[1]);

	if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {
	  
	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamInt->getSectionLocations(numSections, L, xi);
	  
	  output.tag("GaussPointOutput");
	  output.attr("number",sectionNum);
	  output.attr("eta",xi[sectionNum-1]*L);
	  
	  theResponse =  theSections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
	  
	  output.endTag();
	} else if (sectionNum == 0) { // argv[1] was not an int, we want all sections, 
	
	  CompositeResponse *theCResponse = new CompositeResponse();
	  int numResponse = 0;
	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamInt->getSectionLocations(numSections, L, xi);
	  
	  for (int i=0; i<numSections; i++) {
	    
	    output.tag("GaussPointOutput");
	    output.attr("number",i+1);
	    output.attr("eta",xi[i]*L);
	    
	    Response *theSectionResponse = theSections[i]->setResponse(&argv[1], argc-1, output);
	    
	    output.endTag();	  
	    
	    if (theSectionResponse != 0) {
	      numResponse = theCResponse->addResponse(theSectionResponse);
	    }
	  }
	  
	  if (numResponse == 0) // no valid responses found
	    delete theCResponse;
	  else
	    theResponse = theCResponse;
	}
      }
    }
    // by SAJalali
    else if (LadrunoResp::isEnergy(argv[0])) {
      theResponse = new ElementResponse(this, 13, 0.0);
    }
    else if (LadrunoResp::isLocalAxes(argv[0])) {
      theResponse = new ElementResponse(this, 30, Vector(9));
    }
    // Ladruno (ADR 33 / WP-200): forward hinge responses for ends I and J
    else if (hingeOn && argc > 1 &&
             (strcmp(argv[0],"hinge") == 0 || strcmp(argv[0],"hingeI") == 0 ||
              strcmp(argv[0],"hingeZI") == 0 || strcmp(argv[0],"hinge1") == 0)) {
      UniaxialMaterial *hm = (theHingeZ != 0) ? theHingeZ : theHingeZJ;
      if (hm != 0)
        theResponse = hm->setResponse(&argv[1], argc-1, output);
    }
    else if (hingeOn && argc > 1 &&
             (strcmp(argv[0],"hingeJ") == 0 || strcmp(argv[0],"hingeZJ") == 0 || strcmp(argv[0],"hinge2") == 0)) {
      UniaxialMaterial *hm = (theHingeZJ != 0) ? theHingeZJ : theHingeZ;
      if (hm != 0)
        theResponse = hm->setResponse(&argv[1], argc-1, output);
    }
    else if (hingeOn && argc > 1 &&
             (strcmp(argv[0],"hingeY") == 0 || strcmp(argv[0],"hingeYI") == 0 || strcmp(argv[0],"hingeY1") == 0)) {
      UniaxialMaterial *hm = (theHingeY != 0) ? theHingeY : theHingeYJ;
      if (hm != 0)
        theResponse = hm->setResponse(&argv[1], argc-1, output);
    }
    else if (hingeOn && argc > 1 &&
             (strcmp(argv[0],"hingeYJ") == 0 || strcmp(argv[0],"hingeY2") == 0)) {
      UniaxialMaterial *hm = (theHingeYJ != 0) ? theHingeYJ : theHingeY;
      if (hm != 0)
        theResponse = hm->setResponse(&argv[1], argc-1, output);
    }
    else if (hingeOn && argc > 1 &&
             (strcmp(argv[0],"hingeBiaxial") == 0 || strcmp(argv[0],"hingeBiaxialI") == 0 ||
              strcmp(argv[0],"hingeBiaxial1") == 0 || strcmp(argv[0],"hingeC") == 0 || strcmp(argv[0],"hingeCI") == 0)) {
      NDMaterial *hm = (theHingeC != 0) ? theHingeC : theHingeCJ;
      if (hm != 0)
        theResponse = hm->setResponse(&argv[1], argc-1, output);
    }
    else if (hingeOn && argc > 1 &&
             (strcmp(argv[0],"hingeBiaxialJ") == 0 || strcmp(argv[0],"hingeBiaxial2") == 0 || strcmp(argv[0],"hingeCJ") == 0)) {
      NDMaterial *hm = (theHingeCJ != 0) ? theHingeCJ : theHingeC;
      if (hm != 0)
        theResponse = hm->setResponse(&argv[1], argc-1, output);
    }
    else if (hingeOn && LadrunoResp::is(argv[0], "hingeSummary")) {
      output.tag("ResponseType", "thetaZ_i");
      output.tag("ResponseType", "Mz_i");
      output.tag("ResponseType", "thetaZ_j");
      output.tag("ResponseType", "Mz_j");
      output.tag("ResponseType", "thetaY_i");
      output.tag("ResponseType", "My_i");
      output.tag("ResponseType", "thetaY_j");
      output.tag("ResponseType", "My_j");
      theResponse = new ElementResponse(this, 31, Vector(8));
    }
    else if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
      output.tag("ResponseType", "implexError");
      theResponse = new ElementResponse(this, 32, 0.0);
    }
    else if (strcmp(argv[0], "avgImplexError") == 0 || strcmp(argv[0], "AvgImplexError") == 0) {
      output.tag("ResponseType", "avgImplexError");
      theResponse = new ElementResponse(this, 33, 0.0);
    }
    else if (strcmp(argv[0], "hingeJump") == 0 || strcmp(argv[0], "hingeJumpI") == 0) {
      output.tag("ResponseType", "thetaZ_i");
      output.tag("ResponseType", "thetaY_i");
      theResponse = new ElementResponse(this, 34, Vector(2));
    }
    else if (strcmp(argv[0], "hingeJumpJ") == 0) {
      output.tag("ResponseType", "thetaZ_j");
      output.tag("ResponseType", "thetaY_j");
      theResponse = new ElementResponse(this, 35, Vector(2));
    }

    if (theResponse == 0)
      theResponse = crdTransf->setResponse(argv, argc, output);

  output.endTag();

  if (theResponse == 0)
    return this->Element::setResponse(argv, argc, output);
  return theResponse;
}

int 
LadrunoDispBeamColumn3d::getResponse(int responseID, Information &eleInfo)
{
  double N, V, M1, M2, T;
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  // Ladruno: local axes (vx,vy,vz dir cosines) from the CrdTransf
  else if (responseID == 30) {
    static thread_local Vector la(9);
    static thread_local Vector vx(3), vy(3), vz(3);
    crdTransf->getLocalAxes(vx, vy, vz);
    for (int i = 0; i < 3; i++) { la(i) = vx(i); la(i + 3) = vy(i); la(i + 6) = vz(i); }
    return eleInfo.setVector(la);
  }

  else if (responseID == 31) {
    static thread_local Vector hs(8);
    hs.Zero();
    if (hingeOn) {
      hs(0) = hingeJumpZ;
      hs(2) = hingeJumpZJ;
      hs(4) = hingeJumpY;
      hs(6) = hingeJumpYJ;
      double fPM = (pmP0 > 0.0) ? ladrunoDisp3dScalePM(q(0), pmP0, pmAlpha) : 1.0;
      if (theHingeC != 0) {
        const Vector &s = theHingeC->getStress();
        if (s.Size() >= 2) {
          hs(1) = s(0);
          hs(5) = s(1);
        }
      } else {
        if (theHingeZ != 0) {
          double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeZ) == 0) ? fPM : 1.0;
          hs(1) = f * theHingeZ->getStress();
        }
        if (theHingeY != 0) {
          double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeY) == 0) ? fPM : 1.0;
          hs(5) = f * theHingeY->getStress();
        }
      }
      if (theHingeCJ != 0) {
        const Vector &sJ = theHingeCJ->getStress();
        if (sJ.Size() >= 2) {
          hs(3) = sJ(0);
          hs(7) = sJ(1);
        }
      } else {
        if (theHingeZJ != 0) {
          double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeZJ) == 0) ? fPM : 1.0;
          hs(3) = f * theHingeZJ->getStress();
        }
        if (theHingeYJ != 0) {
          double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeYJ) == 0) ? fPM : 1.0;
          hs(7) = f * theHingeYJ->getStress();
        }
      }
    }
    return eleInfo.setVector(hs);
  }

  else if (responseID == 32 || responseID == 33) {
    double maxErr = 0.0, sumErr = 0.0;
    int cnt = 0;
    UniaxialMaterial *uniList[4] = {theHingeZ, theHingeY, theHingeZJ, theHingeYJ};
    for (int k = 0; k < 4; k++) {
      double e = 0.0;
      if (ldbc3dGetHingeImplexError(uniList[k], e)) {
        if (e > maxErr) maxErr = e;
        sumErr += e;
        cnt++;
      }
    }
    NDMaterial *ndList[2] = {theHingeC, theHingeCJ};
    for (int k = 0; k < 2; k++) {
      if (LadrunoCohesiveHingeBiaxial *ch = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(ndList[k])) {
        double e = ch->getImplexError();
        if (e > maxErr) maxErr = e;
        sumErr += e;
        cnt++;
      }
    }
    double val = (responseID == 32) ? maxErr : ((cnt > 0) ? (sumErr / cnt) : 0.0);
    return eleInfo.setDouble(val);
  }

  else if (responseID == 34) {
    static thread_local Vector hjI(2);
    hjI(0) = hingeOn ? hingeJumpZ : 0.0;
    hjI(1) = hingeOn ? hingeJumpY : 0.0;
    return eleInfo.setVector(hjI);
  }

  else if (responseID == 35) {
    static thread_local Vector hjJ(2);
    hjJ(0) = hingeOn ? hingeJumpZJ : 0.0;
    hjJ(1) = hingeOn ? hingeJumpYJ : 0.0;
    return eleInfo.setVector(hjJ);
  }

  else if (responseID == 12) {
    static thread_local Vector damp(12);
    damp.Zero();
    damp.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    return eleInfo.setVector(damp);
  }
    
  else if (responseID == 2) {
    // Axial
    N = q(0);
    P(6) =  N;
    P(0) = -N+p0[0];
    
    // Torsion
    T = q(5);
    P(9) =  T;
    P(3) = -T;
    
    // Moments about z and shears along y
    M1 = q(1);
    M2 = q(2);
    P(5)  = M1;
    P(11) = M2;
    V = (M1+M2)*oneOverL;
    P(1) =  V+p0[1];
    P(7) = -V+p0[2];
    
    // Moments about y and shears along z
    M1 = q(3);
    M2 = q(4);
    P(4)  = M1;
    P(10) = M2;
    V = (M1+M2)*oneOverL;
    P(2) = -V+p0[3];
    P(8) =  V+p0[4];

    return eleInfo.setVector(P);
  }

  else if (responseID == 21)
    return eleInfo.setVector(this->getDampingForce());

  else if (responseID == 22) {
    Vector Sd(6);
    Sd = theDamping->getDampingForce();
    // Axial
    N = Sd(0);
    P(6) =  N;
    P(0) = -N;
    
    // Torsion
    T = Sd(5);
    P(9) =  T;
    P(3) = -T;
    
    // Moments about z and shears along y
    M1 = Sd(1);
    M2 = Sd(2);
    P(5)  = M1;
    P(11) = M2;
    V = (M1+M2)*oneOverL;
    P(1) =  V;
    P(7) = -V;
    
    // Moments about y and shears along z
    M1 = Sd(3);
    M2 = Sd(4);
    P(4)  = M1;
    P(10) = M2;
    V = (M1+M2)*oneOverL;
    P(2) = -V;
    P(8) =  V;
    return eleInfo.setVector(P);
  }

  else if (responseID == 23)
    return eleInfo.setVector(theDamping->getDampingForce());

  else if (responseID == 9) {
    return eleInfo.setVector(q);
  }

  else if (responseID == 19) {
    static thread_local Matrix kb(6,6);
    this->getBasicStiff(kb);
    return eleInfo.setMatrix(kb);
  }  
  // Chord rotation
  else if (responseID == 3)
    return eleInfo.setVector(crdTransf->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static thread_local Vector vp(6);
    static thread_local Vector ve(6);
    static thread_local Vector qSec(6);
    static thread_local Matrix kb(6,6);
    this->getBasicStiff(kb,1);
    for (int m = 0; m < 5; m++)
      qSec(m) = q(m) - q0[m];
    qSec(5) = q(5);
    double L = crdTransf->getInitialLength();
    if (hingeOn && !dualHingeMode && L > 0.0) {
      double oneOverL = 1.0/L;
      double xi[maxNumSections];
      double wt[maxNumSections];
      beamInt->getSectionLocations(numSections, L, xi);
      beamInt->getSectionWeights(numSections, L, wt);
      double KvZ0[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      double KvY0[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      for (int i = 0; i < numSections; i++) {
        int order = theSections[i]->getOrder();
        const ID &code = theSections[i]->getType();
        const Matrix &ks0 = theSections[i]->getInitialTangent();
        double xi6 = 6.0*xi[i];
        int idxP = -1, idxMZ = -1, idxMY = -1, idxT = -1;
        for (int j = 0; j < order; j++) {
          switch (code(j)) {
          case SECTION_RESPONSE_P:  idxP  = j; break;
          case SECTION_RESPONSE_MZ: idxMZ = j; break;
          case SECTION_RESPONSE_MY: idxMY = j; break;
          case SECTION_RESPONSE_T:  idxT  = j; break;
          default: break;
          }
        }
        if (idxMZ >= 0) {
          double EIz = ks0(idxMZ, idxMZ);
          KvZ0[1] += -oneOverL*wt[i]*(xi6-4.0)*EIz;
          KvZ0[2] += -oneOverL*wt[i]*(xi6-2.0)*EIz;
          if (idxP  >= 0) KvZ0[0] += -oneOverL*wt[i]*ks0(idxP,  idxMZ);
          if (idxMY >= 0) {
            KvZ0[3] += -oneOverL*wt[i]*(xi6-4.0)*ks0(idxMY, idxMZ);
            KvZ0[4] += -oneOverL*wt[i]*(xi6-2.0)*ks0(idxMY, idxMZ);
          }
          if (idxT  >= 0) KvZ0[5] += -oneOverL*wt[i]*ks0(idxT,  idxMZ);
        }
        if (idxMY >= 0) {
          double EIy = ks0(idxMY, idxMY);
          KvY0[3] += -oneOverL*wt[i]*(xi6-4.0)*EIy;
          KvY0[4] += -oneOverL*wt[i]*(xi6-2.0)*EIy;
          if (idxP  >= 0) KvY0[0] += -oneOverL*wt[i]*ks0(idxP,  idxMY);
          if (idxMZ >= 0) {
            KvY0[1] += -oneOverL*wt[i]*(xi6-4.0)*ks0(idxMZ, idxMY);
            KvY0[2] += -oneOverL*wt[i]*(xi6-2.0)*ks0(idxMZ, idxMY);
          }
          if (idxT  >= 0) KvY0[5] += -oneOverL*wt[i]*ks0(idxT,  idxMY);
        }
      }
      for (int a = 0; a < 6; a++)
        qSec(a) -= KvZ0[a]*hingeJumpZ + KvY0[a]*hingeJumpY;
      kb.Solve(qSec, ve);
      vp = crdTransf->getBasicTrialDisp();
      vp -= ve;
      vp(1) += (1.0 - hingeXiI)*hingeJumpZ;
      vp(2) += hingeXiI*hingeJumpZ;
      vp(3) += (1.0 - hingeXiI)*hingeJumpY;
      vp(4) += hingeXiI*hingeJumpY;
    } else {
      // Note: in dualHingeMode, bulk sections are evaluated at v_tilde = v - W*alpha,
      // so qSec = kb0 * (v - W*alpha - vp_bulk) and v - kb0^-1*qSec = vp_bulk + W*alpha
      // already includes the Hermite-projected hinge jumps W*alpha.
      kb.Solve(qSec, ve);
      vp = crdTransf->getBasicTrialDisp();
      vp -= ve;
    }
    return eleInfo.setVector(vp);
  }

  else if (responseID == 10) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamInt->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i]*L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 11) {
    double L = crdTransf->getInitialLength();
    double wts[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i]*L;
    return eleInfo.setVector(weights);
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = theSections[i]->getTag();
    return eleInfo.setID(tags);
  }
  
  //by SAJalali
  else if (responseID == 13) {
	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamInt->getSectionWeights(numSections, L, xi);
	  double energy = 0;
	  for (int i = 0; i < numSections; i++) {
		  energy += theSections[i]->getEnergy()*xi[i] * L;
	  }
	  if (hingeOn) {
	    double fPM = (pmP0 > 0.0) ? ladrunoDisp3dScalePM(q(0), pmP0, pmAlpha) : 1.0;
	    if (theHingeZ != 0) {
	      double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeZ) == 0) ? fPM : 1.0;
	      energy += f * theHingeZ->getEnergy();
	    }
	    if (theHingeY != 0) {
	      double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeY) == 0) ? fPM : 1.0;
	      energy += f * theHingeY->getEnergy();
	    }
	    if (theHingeC != 0) {
	      LadrunoCohesiveHingeBiaxial *hc = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeC);
	      if (hc != 0)
	        energy += hc->getEnergy();
	    }
	    if (theHingeZJ != 0) {
	      double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeZJ) == 0) ? fPM : 1.0;
	      energy += f * theHingeZJ->getEnergy();
	    }
	    if (theHingeYJ != 0) {
	      double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeYJ) == 0) ? fPM : 1.0;
	      energy += f * theHingeYJ->getEnergy();
	    }
	    if (theHingeCJ != 0) {
	      LadrunoCohesiveHingeBiaxial *hcJ = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeCJ);
	      if (hcJ != 0)
	        energy += hcJ->getEnergy();
	    }
	  }
	  return eleInfo.setDouble(energy);
  }

  else
    return this->Element::getResponse(responseID, eleInfo);
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
LadrunoDispBeamColumn3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // don't do anything if MaterialStageParameter calls this element
  if (strcmp(argv[0],"updateMaterialStage") == 0) {
      return -1;
  }
  
  // If the parameter belongs to the element itself
  if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(1, this);
  }

  if (hingeOn && argc >= 2) {
    if (strcmp(argv[0], "hingeI") == 0 || strcmp(argv[0], "hingeZI") == 0 || strcmp(argv[0], "hinge1") == 0) {
      if (theHingeZ != 0) return theHingeZ->setParameter(&argv[1], argc - 1, param);
      if (theHingeC != 0) return theHingeC->setParameter(&argv[1], argc - 1, param);
      return -1;
    }
    if (strcmp(argv[0], "hingeJ") == 0 || strcmp(argv[0], "hingeZJ") == 0 || strcmp(argv[0], "hinge2") == 0) {
      if (theHingeZJ != 0) return theHingeZJ->setParameter(&argv[1], argc - 1, param);
      if (theHingeCJ != 0) return theHingeCJ->setParameter(&argv[1], argc - 1, param);
      return -1;
    }
    if (strcmp(argv[0], "hinge") == 0 || strcmp(argv[0], "hingeZ") == 0 ||
        strcmp(argv[0], "hingeIJ") == 0 || strcmp(argv[0], "hingeZIJ") == 0) {
      int res = -1;
      if (theHingeZ != 0) {
        int ok = theHingeZ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeZJ != 0) {
        int ok = theHingeZJ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeC != 0) {
        int ok = theHingeC->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeCJ != 0) {
        int ok = theHingeCJ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      return res;
    }
    if (strcmp(argv[0], "hingeY") == 0 || strcmp(argv[0], "hingeYIJ") == 0) {
      int res = -1;
      if (theHingeY != 0) {
        int ok = theHingeY->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeYJ != 0) {
        int ok = theHingeYJ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      return res;
    }
    if (strcmp(argv[0], "hingeYI") == 0 || strcmp(argv[0], "hingeY1") == 0) {
      if (theHingeY != 0) return theHingeY->setParameter(&argv[1], argc - 1, param);
      return -1;
    }
    if (strcmp(argv[0], "hingeYJ") == 0 || strcmp(argv[0], "hingeY2") == 0) {
      if (theHingeYJ != 0) return theHingeYJ->setParameter(&argv[1], argc - 1, param);
      return -1;
    }
    if (strcmp(argv[0], "hingeBiaxial") == 0 || strcmp(argv[0], "hingeC") == 0 ||
        strcmp(argv[0], "hingeBiaxialIJ") == 0 || strcmp(argv[0], "hingeCIJ") == 0) {
      int res = -1;
      if (theHingeC != 0) {
        int ok = theHingeC->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeCJ != 0) {
        int ok = theHingeCJ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      return res;
    }
    if (strcmp(argv[0], "hingeBiaxialI") == 0 || strcmp(argv[0], "hingeBiaxial1") == 0 ||
        strcmp(argv[0], "hingeCI") == 0) {
      if (theHingeC != 0) return theHingeC->setParameter(&argv[1], argc - 1, param);
      return -1;
    }
    if (strcmp(argv[0], "hingeBiaxialJ") == 0 || strcmp(argv[0], "hingeBiaxial2") == 0 ||
        strcmp(argv[0], "hingeCJ") == 0) {
      if (theHingeCJ != 0) return theHingeCJ->setParameter(&argv[1], argc - 1, param);
      return -1;
    }
  }

  // damping
  if (strstr(argv[0], "damp") != 0) {

    if (argc < 2 || !theDamping)
      return -1;

    return theDamping->setParameter(&argv[1], argc-1, param);
  }

  if (strstr(argv[0],"sectionX") != 0) {
    if (argc < 3)
		return -1;
      
	float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamInt->getSectionLocations(numSections, L, xi);
      
      sectionLoc /= L;

      float minDistance = fabs(xi[0]-sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
	if (fabs(xi[i]-sectionLoc) < minDistance) {
	  minDistance = fabs(xi[i]-sectionLoc);
	  sectionNum = i;
	}
	  }  
	return theSections[sectionNum]->setParameter(&argv[2], argc-2, param);
  }
  // If the parameter belongs to a section or lower
  if (strstr(argv[0],"section") != 0) {
    
    if (argc < 3)
      return -1;
    
    // Get section number
    int sectionNum = atoi(argv[1]);
    
    if (sectionNum > 0 && sectionNum <= numSections)
      return theSections[sectionNum-1]->setParameter(&argv[2], argc-2, param);
    else
      return -1;
  }
  
  if (strstr(argv[0],"integration") != 0) {
    
    if (argc < 2)
      return -1;

    return beamInt->setParameter(&argv[1], argc-1, param);
  }

  // Default, send to every object
  int ok = 0;
  int result = -1;

  for (int i = 0; i < numSections; i++) {
    ok = theSections[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  
  ok = beamInt->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  if (hingeOn) {
    UniaxialMaterial *uniList[4] = {theHingeZ, theHingeY, theHingeZJ, theHingeYJ};
    for (int k = 0; k < 4; k++) {
      if (uniList[k] != 0) {
        ok = uniList[k]->setParameter(argv, argc, param);
        if (ok != -1) result = ok;
      }
    }
    NDMaterial *ndList[2] = {theHingeC, theHingeCJ};
    for (int k = 0; k < 2; k++) {
      if (ndList[k] != 0) {
        ok = ndList[k]->setParameter(argv, argc, param);
        if (ok != -1) result = ok;
      }
    }
  }

  return result;
}

int
LadrunoDispBeamColumn3d::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1;  
}




int
LadrunoDispBeamColumn3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}

const Matrix &
LadrunoDispBeamColumn3d::getInitialStiffSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Matrix &
LadrunoDispBeamColumn3d::getMassSensitivity(int gradNumber)
{
  K.Zero();
  
  if (rho == 0.0 || parameterID != 1)
    return K;
  
  double L = crdTransf->getInitialLength();
  if (cMass == 0)  {
    // lumped mass matrix
    //double m = 0.5*rho*L;
    double m = 0.5*L;
    K(0,0) = K(1,1) = K(2,2) = K(6,6) = K(7,7) = K(8,8) = m;
  } else  {
    // consistent mass matrix
    static thread_local Matrix ml(12,12);
    //double m = rho*L/420.0;
    double m = L/420.0;
    ml(0,0) = ml(6,6) = m*140.0;
    ml(0,6) = ml(6,0) = m*70.0;
    
    ml(2,2) = ml(8,8) = m*156.0;
    ml(2,8) = ml(8,2) = m*54.0;
    ml(4,4) = ml(10,10) = m*4.0*L*L;
    ml(4,10) = ml(10,4) = -m*3.0*L*L;
    ml(2,4) = ml(4,2) = -m*22.0*L;
    ml(8,10) = ml(10,8) = -ml(2,4);
    ml(2,10) = ml(10,2) = m*13.0*L;
    ml(4,8) = ml(8,4) = -ml(2,10);
    
    ml(1,1) = ml(7,7) = m*156.0;
    ml(1,7) = ml(7,1) = m*54.0;
    ml(5,5) = ml(11,11) = m*4.0*L*L;
    ml(5,11) = ml(11,5) = -m*3.0*L*L;
    ml(1,5) = ml(5,1) = m*22.0*L;
    ml(7,11) = ml(11,7) = -ml(1,5);
    ml(1,11) = ml(11,1) = -m*13.0*L;
    ml(5,7) = ml(7,5) = -ml(1,11);
    
    // transform local mass matrix to global system
    K = crdTransf->getGlobalMatrixFromLocal(ml);
  }
  
  return K;
}



const Vector &
LadrunoDispBeamColumn3d::getResistingForceSensitivity(int gradNumber)
{
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Zero for integration
  static thread_local Vector dqdh(6);
  dqdh.Zero();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    
    double xi6 = 6.0*xi[i];
    double wti = wt[i];
    
    // Get section stress resultant gradient
    const Vector &dsdh = theSections[i]->getStressResultantSensitivity(gradNumber,true);
    
    // Perform numerical integration on internal force gradient
    double sensi;
    for (int j = 0; j < order; j++) {
      sensi = dsdh(j)*wti;
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	dqdh(0) += sensi; 
	break;
      case SECTION_RESPONSE_MZ:
	dqdh(1) += (xi6-4.0)*sensi; 
	dqdh(2) += (xi6-2.0)*sensi; 
	break;
      case SECTION_RESPONSE_MY:
	dqdh(3) += (xi6-4.0)*sensi; 
	dqdh(4) += (xi6-2.0)*sensi; 
	break;
      case SECTION_RESPONSE_T:
	dqdh(5) += sensi; 
	break;
      default:
	break;
      }
    }
  }
  
  // Transform forces
  static thread_local Vector dp0dh(6);		// No distributed loads

  P.Zero();

  //////////////////////////////////////////////////////////////

  if (crdTransf->isShapeSensitivity()) {
    
    // Perform numerical integration to obtain basic stiffness matrix
    // Some extra declarations
    static thread_local Matrix kbmine(6,6);
    kbmine.Zero();
    q.Zero();
    
    double tmp;
    
    int j, k;
    
    for (int i = 0; i < numSections; i++) {
      
      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      
      double xi6 = 6.0*xi[i];
      double wti = wt[i];
      
      const Vector &s = theSections[i]->getStressResultant();
      const Matrix &ks = theSections[i]->getSectionTangent();
      
      Matrix ka(workArea, order, 6);
      ka.Zero();
      
      double si;
      for (j = 0; j < order; j++) {
	si = s(j)*wti;
	switch(code(j)) {
	case SECTION_RESPONSE_P:
	  q(0) += si;
	  for (k = 0; k < order; k++) {
	    ka(k,0) += ks(k,j)*wti;
	  }
	  break;
	case SECTION_RESPONSE_MZ:
	  q(1) += (xi6-4.0)*si; 
	  q(2) += (xi6-2.0)*si;
	  for (k = 0; k < order; k++) {
	    tmp = ks(k,j)*wti;
	    ka(k,1) += (xi6-4.0)*tmp;
	    ka(k,2) += (xi6-2.0)*tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  q(3) += (xi6-4.0)*si; 
	  q(4) += (xi6-2.0)*si;
	  for (k = 0; k < order; k++) {
	    tmp = ks(k,j)*wti;
	    ka(k,3) += (xi6-4.0)*tmp;
	    ka(k,4) += (xi6-2.0)*tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  q(5) += si;
	  for (k = 0; k < order; k++) {
	    ka(k,5) += ks(k,j)*wti;
	  }
	  break;
	default:
	  break;
	}
      }
      for (j = 0; j < order; j++) {
	switch (code(j)) {
	case SECTION_RESPONSE_P:
	  for (k = 0; k < 6; k++) {
	    kbmine(0,k) += ka(j,k);
	  }
	  break;
	case SECTION_RESPONSE_MZ:
	  for (k = 0; k < 6; k++) {
	    tmp = ka(j,k);
	    kbmine(1,k) += (xi6-4.0)*tmp;
	    kbmine(2,k) += (xi6-2.0)*tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  for (k = 0; k < 6; k++) {
	    tmp = ka(j,k);
	    kbmine(3,k) += (xi6-4.0)*tmp;
	    kbmine(4,k) += (xi6-2.0)*tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  for (k = 0; k < 6; k++) {
	    kbmine(5,k) += ka(j,k);
	  }
	  break;
	default:
	  break;
	}
      }
    }      
    
    const Vector &A_u = crdTransf->getBasicTrialDisp();
    double dLdh = crdTransf->getdLdh();
    double d1overLdh = -dLdh/(L*L);
    // a^T k_s dadh v
    dqdh.addMatrixVector(1.0, kbmine, A_u, d1overLdh);

    // k dAdh u
    const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity();
    dqdh.addMatrixVector(1.0, kbmine, dAdh_u, oneOverL);

    // dAdh^T q
    P += crdTransf->getGlobalResistingForceShapeSensitivity(q, dp0dh, gradNumber);
  }
  
  // A^T (dqdh + k dAdh u)
  P += crdTransf->getGlobalResistingForce(dqdh, dp0dh);
  
  return P;
}



// NEW METHOD
int
LadrunoDispBeamColumn3d::commitSensitivity(int gradNumber, int numGrads)
{
  // Get basic deformation and sensitivities
  const Vector &v = crdTransf->getBasicTrialDisp();
  
  static thread_local Vector dvdh(6);
  dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);

  // Some extra declarations
  double d1oLdh = crdTransf->getd1overLdh();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    
    Vector e(workArea, order);
    
    double xi6 = 6.0*xi[i];
    
    for (int j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	e(j) = oneOverL*dvdh(0)
	  + d1oLdh*v(0); 
	break;
      case SECTION_RESPONSE_MZ:
	e(j) = oneOverL*((xi6-4.0)*dvdh(1) + (xi6-2.0)*dvdh(2))
	  + d1oLdh*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)); 
	break;
      case SECTION_RESPONSE_MY:
	e(j) = oneOverL*((xi6-4.0)*dvdh(3) + (xi6-2.0)*dvdh(4))
	  + d1oLdh*((xi6-4.0)*v(3) + (xi6-2.0)*v(4)); 
	break;
      case SECTION_RESPONSE_T:
	e(j) = oneOverL*dvdh(5)
	  + d1oLdh*v(5); 
	break;
      default:
	e(j) = 0.0; 
	break;
      }
    }
    
    // Set the section deformations
    theSections[i]->commitSensitivity(e,gradNumber,numGrads);
  }
  
  return 0;
}


// AddingSensitivity:END /////////////////////////////////////////////

bool
LadrunoDispBeamColumn3d::ladrunoThreadSafeUpdate(void) const
{
  if (theDamping != 0)
    return false;
  for (int i = 0; i < numSections; i++) {
    if (theSections[i] == 0 || !theSections[i]->ladrunoThreadSafeUpdate())
      return false;
  }
  if (hingeOn) {
    if (theHingeZ != 0 && !theHingeZ->ladrunoThreadSafeUpdate())
      return false;
    if (theHingeY != 0 && !theHingeY->ladrunoThreadSafeUpdate())
      return false;
    if (theHingeC != 0 && !theHingeC->ladrunoThreadSafeUpdate())
      return false;
    if (theHingeZJ != 0 && !theHingeZJ->ladrunoThreadSafeUpdate())
      return false;
    if (theHingeYJ != 0 && !theHingeYJ->ladrunoThreadSafeUpdate())
      return false;
    if (theHingeCJ != 0 && !theHingeCJ->ladrunoThreadSafeUpdate())
      return false;
  }
  return true;
}

bool
LadrunoDispBeamColumn3d::ladrunoThreadSafeForm(void) const
{
  if (theDamping != 0)
    return false;
  for (int i = 0; i < numSections; i++) {
    if (theSections[i] == 0 || !theSections[i]->ladrunoThreadSafeForm())
      return false;
  }
  if (hingeOn) {
    if (theHingeZ != 0 && !theHingeZ->ladrunoThreadSafeForm())
      return false;
    if (theHingeY != 0 && !theHingeY->ladrunoThreadSafeForm())
      return false;
    if (theHingeC != 0 && !theHingeC->ladrunoThreadSafeForm())
      return false;
    if (theHingeZJ != 0 && !theHingeZJ->ladrunoThreadSafeForm())
      return false;
    if (theHingeYJ != 0 && !theHingeYJ->ladrunoThreadSafeForm())
      return false;
    if (theHingeCJ != 0 && !theHingeCJ->ladrunoThreadSafeForm())
      return false;
  }
  return true;
}


