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
// Description: This file contains the class definition for LadrunoDispBeamColumn2d.

#include <LadrunoDispBeamColumn2d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>   // Ladruno (ADR 32) Tier-2: embedded cohesive hinge
#include <LadrunoCohesiveHinge.h>
#include <Bilin.h>
#include <ModIMKPeakOriented.h>
#include <CrdTransf.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <stdlib.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <LadrunoResponseTokens.h>
#include <math.h>
#include <cmath>
#include <ElementalLoad.h>
#include <elementAPI.h>
#include <string>
#include <string.h>
#include <map>
#include <ElementIter.h>
#include <LadrunoThreads.h>

// Ladruno (ADR 32): the element currently being updated, set by the framework
// (Domain::update) before update() and read by crack-band materials inside
// setTrialStrain via getCharacteristicLength(). Declared here for the defensive
// re-assignment at the top of update().
extern thread_local Element *ops_TheActiveElement;   // Ladruno WP-107: thread_local, see Element.cpp

static thread_local Matrix K(6,6);
static thread_local Vector P(6);
static thread_local double workArea[100];

static bool ldbc2dTryOptionalDouble(double& out)
{
  if (OPS_GetNumRemainingInputArgs() < 1) return false;
  char tok[64];
  const char* s = OPS_GetStringFromAll(tok, sizeof(tok));
  OPS_ResetCurrentInputArg(-1);
  if (s == 0 || s[0] == '\0') return false;
  char* endp = 0;
  (void)strtod(s, &endp);
  if (endp == s || *endp != '\0') return false;
  int one = 1;
  return (OPS_GetDoubleInput(&one, &out) == 0);
}

static bool ldbc2dTryOptionalInt(int& out)
{
  if (OPS_GetNumRemainingInputArgs() < 1) return false;
  char tok[64];
  const char* s = OPS_GetStringFromAll(tok, sizeof(tok));
  OPS_ResetCurrentInputArg(-1);
  if (s == 0 || s[0] == '\0') return false;
  char* endp = 0;
  long v = strtol(s, &endp, 10);
  if (endp == s || *endp != '\0') return false;
  int one = 1;
  if (OPS_GetIntInput(&one, &out) < 0) return false;
  (void)v;
  return true;
}

static bool ldbc2dGetHingeImplexError(UniaxialMaterial* mat, double& errOut)
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

static bool ldbc2dHingeIsImplex(UniaxialMaterial* mat)
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

static void invertSym2x2Guarded2d(double A, double B, double C, double (&Hinv)[2][2])
{
  double avg = 0.5*(A + B);
  double hdiff = 0.5*(A - B);
  double disc = sqrt(hdiff*hdiff + C*C);
  double lam1 = avg + disc;
  double lam2 = avg - disc;
  double v1z, v1y, v2z, v2y;
  if (fabs(C) > 1.0e-300) {
    double n1 = sqrt(C*C + (lam1 - A)*(lam1 - A));
    v1z = C/n1;           v1y = (lam1 - A)/n1;
    v2z = -v1y;           v2y =  v1z;
  } else {
    if (A >= B) { v1z = 1.0; v1y = 0.0; v2z = 0.0; v2y = 1.0; }
    else        { v1z = 0.0; v1y = 1.0; v2z = -1.0; v2y = 0.0; }
  }
  double scale = fabs(A) + fabs(B) + 2.0*fabs(C);
  double lamFloor = 1.0e-8*scale + 1.0e-300;
  double l1s = (fabs(lam1) < lamFloor) ? ((lam1 < 0.0) ? -lamFloor : lamFloor) : lam1;
  double l2s = (fabs(lam2) < lamFloor) ? ((lam2 < 0.0) ? -lamFloor : lamFloor) : lam2;
  double il1 = 1.0/l1s, il2 = 1.0/l2s;
  Hinv[0][0] = il1*v1z*v1z + il2*v2z*v2z;
  Hinv[0][1] = il1*v1z*v1y + il2*v2z*v2y;
  Hinv[1][0] = Hinv[0][1];
  Hinv[1][1] = il1*v1y*v1y + il2*v2y*v2y;
}

static inline double ladrunoDisp2dScalePM(double P_axial, double P0, double alphaPM)
{
  if (P0 <= 0.0)
    return 1.0;
  double n = P_axial / P0;
  double f = 1.0 - alphaPM * n - (1.0 - alphaPM) * n * n;
  if (f < 0.05) f = 0.05;
  if (f > 2.5)  f = 2.5;
  return f;
}

void* OPS_LadrunoDispBeamColumn2d()
{
    if (OPS_GetNDM() != 2 || OPS_GetNDF() != 3) {
        opserr << "WARNING LadrunoDispBeamColumn2d -- model must be ndm 2/ndf 3\n";
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
    int nlGeom = 0;        // Ladruno (ADR 32): 0=linear basic strain, 1=½θ² bowing (-nl)
    UniaxialMaterial *hingeMat = 0;   // Ladruno (ADR 32) Tier-2: embedded cohesive hinge (-hinge $matTag)
    UniaxialMaterial *hingeMatJ = 0;  // Ladruno: second end-J embedded cohesive hinge (-hingeJ $matTagJ)
    int dualHingeMode = 0;            // 0=legacy constant G=-1/L, 1=Hermite end-rotation modes G_I(xi), G_J(xi)
    double hingeXiI = 0.0;
    double hingeXiJ = 1.0;
    bool eleImplex = false;
    double elePmP0 = 0.0;
    double elePmAlpha = 0.5;
    int spanLoadMode = 1;
    numData = 1;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type, "-cMass") == 0 || strcmp(type, "-cmass") == 0) {
	    cmass = 1;
	} else if(strcmp(type,"-mass") == 0) {
	    if(OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"WARNING LadrunoDispBeamColumn2d: -mass needs a value\n";
		return 0;
	    }
	    if(OPS_GetDoubleInput(&numData,&mass) < 0) {
		opserr<<"WARNING: invalid mass\n";
		return 0;
	    }
	}
    else if (strcmp(type, "-damp") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn2d: -damp needs a damping tag\n";
            return 0;
        }
        if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
        theDamping = OPS_getDamping(dampingTag);
        if (theDamping == 0) {
            opserr << "damping not found\n";
            return 0;
        }
    }
    // Ladruno (ADR 32): -lch {ip|element|<value>} selects the regularization
    // (characteristic) length reported to crack-band materials.
    else if (strcmp(type, "-lch") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn2d: -lch needs an argument {ip|element|<value>}\n";
            return 0;
        }
        char lchBuf[64];
        const char* lchArg = OPS_GetStringFromAll(lchBuf, sizeof(lchBuf));
        if (lchArg == 0) {
            opserr << "WARNING LadrunoDispBeamColumn2d: failed to read -lch argument\n";
            return 0;
        }
        if (strcmp(lchArg, "ip") == 0 || strcmp(lchArg, "IP") == 0 || strcmp(lchArg, "wtL") == 0) {
            lchMode = 0;
            userLch = 0.0;
        } else if (strcmp(lchArg, "element") == 0 || strcmp(lchArg, "L") == 0) {
            lchMode = 1;
            userLch = 0.0;
            opserr << "WARNING LadrunoDispBeamColumn2d: -lch element regularizes softening over "
                   << "the WHOLE element length; this re-enables multi-integration-point energy "
                   << "double-counting and is for A/B debugging only. Use -lch ip (default).\n";
        } else {
            // treat the argument as a numeric crack-band width
            char* endptr = 0;
            double val = strtod(lchArg, &endptr);
            // reject nan/inf/infinity, trailing garbage, and <= 0;
            // a non-finite band would silently poison the material's regularization
            if (endptr == lchArg || *endptr != '\0' || !std::isfinite(val) || val <= 0.0) {
                opserr << "WARNING LadrunoDispBeamColumn2d: -lch value must be ip|element|<finite positive value>, got '"
                       << lchArg << "'\n";
                return 0;
            }
            userLch = val;
            lchMode = 2;
        }
    }
    // Ladruno (ADR 32): -nl enables the ½θ² bowing (nonlinear) basic strain.
    else if (strcmp(type, "-nl") == 0) {
        nlGeom = 1;
    }
    else if (strcmp(type, "-spanLoad") == 0 || strcmp(type, "-spanLoadEnhancement") == 0) {
        int optFlag = 1;
        if (ldbc2dTryOptionalInt(optFlag))
            spanLoadMode = (optFlag != 0) ? 1 : 0;
        else
            spanLoadMode = 1;
    }
    else if (strcmp(type, "-noSpanLoad") == 0) {
        spanLoadMode = 0;
    }
    // Ladruno (ADR 32) Tier-2: -hinge $matTag embeds a discrete cohesive rotation-jump
    // hinge carrying any UniaxialMaterial (LadrunoCohesiveHinge is the intended default).
    else if (strcmp(type, "-hinge") == 0 || strcmp(type, "-hingeZ") == 0 ||
             strcmp(type, "-hingeI") == 0 || strcmp(type, "-hingeZI") == 0 ||
             strcmp(type, "-hinge1") == 0) {
        const bool allowPair = (strcmp(type, "-hinge") == 0 || strcmp(type, "-hingeZ") == 0);
        if (!allowPair)
            dualHingeMode = 1;
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn2d: " << type << " needs a uniaxialMaterial tag\n";
            return 0;
        }
        int hingeTag;
        if (OPS_GetIntInput(&numData, &hingeTag) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn2d: invalid " << type << " material tag\n";
            return 0;
        }
        hingeMat = OPS_getUniaxialMaterial(hingeTag);
        if (hingeMat == 0) {
            opserr << "WARNING LadrunoDispBeamColumn2d: " << type << " uniaxialMaterial " << hingeTag
                   << " not found\n";
            return 0;
        }
        int optJ = 0;
        if (allowPair && ldbc2dTryOptionalInt(optJ)) {
            dualHingeMode = 1;
            hingeMatJ = OPS_getUniaxialMaterial(optJ);
            if (hingeMatJ == 0) {
                opserr << "WARNING LadrunoDispBeamColumn2d: " << type << " uniaxialMaterial " << optJ
                       << " not found\n";
                return 0;
            }
        }
    }
    else if (strcmp(type, "-hingeJ") == 0 || strcmp(type, "-hingeZJ") == 0 || strcmp(type, "-hinge2") == 0) {
        dualHingeMode = 1;
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn2d: " << type << " needs a uniaxialMaterial tag\n";
            return 0;
        }
        int hingeTagJ;
        if (OPS_GetIntInput(&numData, &hingeTagJ) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn2d: invalid " << type << " material tag\n";
            return 0;
        }
        hingeMatJ = OPS_getUniaxialMaterial(hingeTagJ);
        if (hingeMatJ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn2d: " << type << " uniaxialMaterial " << hingeTagJ
                   << " not found\n";
            return 0;
        }
    }
    else if (strcmp(type, "-hingeIJ") == 0 || strcmp(type, "-hingeZIJ") == 0) {
        dualHingeMode = 1;
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoDispBeamColumn2d: " << type << " needs at least one uniaxialMaterial tag\n";
            return 0;
        }
        int tagI = 0;
        if (OPS_GetIntInput(&numData, &tagI) < 0) {
            opserr << "WARNING LadrunoDispBeamColumn2d: invalid " << type << " material tag\n";
            return 0;
        }
        int tagJ = tagI;
        int optJ = 0;
        if (ldbc2dTryOptionalInt(optJ))
            tagJ = optJ;
        hingeMat = OPS_getUniaxialMaterial(tagI);
        hingeMatJ = OPS_getUniaxialMaterial(tagJ);
        if (hingeMat == 0 || hingeMatJ == 0) {
            opserr << "WARNING LadrunoDispBeamColumn2d: " << type << " uniaxialMaterial not found\n";
            return 0;
        }
    }
    else if (strcmp(type, "-hingeLoc") == 0) {
        dualHingeMode = 1;
        if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &hingeXiI) < 0 ||
            !std::isfinite(hingeXiI) || hingeXiI < 0.0 || hingeXiI > 1.0) {
            opserr << "WARNING LadrunoDispBeamColumn2d: -hingeLoc needs xiI in [0, 1]\n";
            return 0;
        }
        double optXiJ = 0.0;
        if (ldbc2dTryOptionalDouble(optXiJ)) {
            if (!std::isfinite(optXiJ) || optXiJ < 0.0 || optXiJ > 1.0) {
                opserr << "WARNING LadrunoDispBeamColumn2d: -hingeLoc xiJ must be in [0, 1]\n";
                return 0;
            }
            hingeXiJ = optXiJ;
        } else {
            hingeXiJ = 1.0 - hingeXiI;
        }
    }
    else if (strcmp(type, "-implex") == 0 || strcmp(type, "-IMPLEX") == 0) {
        eleImplex = true;
    }
    else if (strcmp(type, "-pm") == 0 || strcmp(type, "-PM") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &elePmP0) < 0 || elePmP0 <= 0.0) {
            opserr << "WARNING LadrunoDispBeamColumn2d: -pm needs P0 > 0\n";
            return 0;
        }
        double optAlpha = elePmAlpha;
        if (ldbc2dTryOptionalDouble(optAlpha)) {
            elePmAlpha = optAlpha;
        }
    }
    else {
        opserr << "WARNING LadrunoDispBeamColumn2d: unknown option '" << type << "'\n";
        return 0;
    }
    }

    if (hingeMatJ != 0)
        dualHingeMode = 1;

    // Ladruno (ADR 32) Tier-2: the embedded hinge and the -nl bowing strain are mutually
    // exclusive in v1 — the ½θ² membrane term couples the rotation jump into the axial
    // channel (the "work-conjugate to the basic moment" story breaks), so the cross-term
    // algebra is deferred. Reject the combination loudly rather than silently mis-condense.
    if ((hingeMat != 0 || hingeMatJ != 0) && nlGeom != 0) {
        opserr << "WARNING LadrunoDispBeamColumn2d: -hinge and -nl are mutually exclusive in v1 "
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
	opserr << "WARNING LadrunoDispBeamColumn2d: number of sections (" << secTags.Size()
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

    UniaxialMaterial *hI_tmp = 0, *hJ_tmp = 0;
    if (eleImplex || elePmP0 > 0.0) {
        auto cfgHinge = [&](UniaxialMaterial *&hPtr, UniaxialMaterial *&tmpPtr) {
            if (hPtr != 0) {
                tmpPtr = hPtr->getCopy();
                if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge*>(tmpPtr)) {
                    if (eleImplex) ch->setImplex(true);
                    if (elePmP0 > 0.0) ch->setPM(elePmP0, elePmAlpha);
                } else if (Bilin *b = dynamic_cast<Bilin*>(tmpPtr)) {
                    if (eleImplex) b->setImplex(true);
                } else if (ModIMKPeakOriented *m = dynamic_cast<ModIMKPeakOriented*>(tmpPtr)) {
                    if (eleImplex) m->setImplex(true);
                }
                hPtr = tmpPtr;
            }
        };
        cfgHinge(hingeMat, hI_tmp);
        cfgHinge(hingeMatJ, hJ_tmp);
    }
    
    Element *theEle =  new LadrunoDispBeamColumn2d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					    *bi,*theTransf,mass,cmass,theDamping,lchMode,userLch,nlGeom,
					    hingeMat,hingeMatJ,dualHingeMode,hingeXiI,hingeXiJ,elePmP0,elePmAlpha,spanLoadMode);
    delete hI_tmp;
    delete hJ_tmp;
    delete [] sections;
    return theEle;
}

LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration& bi,
				   CrdTransf &coordTransf, double r, int cm,
				   Damping *damping, int lchM, double userL, int nlG,
				   UniaxialMaterial *hingeMat,
				   UniaxialMaterial *hingeMatJ,
				   int dualMode, double xiI, double xiJ,
				   double pmP0_, double pmAlpha_, int spanMode)
:Element (tag, ELE_TAG_LadrunoDispBeamColumn2d),
 numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
  connectedExternalNodes(2),
  Q(6), q(3), rho(r), cMass(cm), parameterID(0), theDamping(0),
  current_section_lch(0.0), lchMode(lchM), userLch(userL), nlGeom(nlG),
  hingeOn(0), dualHingeMode(dualMode), hingeXiI(xiI), hingeXiJ(xiJ),
  theHinge(0), theHingeJ(0),
  hingeJump(0.0), hingeJumpCommit(0.0),
  hingeJumpJ(0.0), hingeJumpJCommit(0.0),
  hingeKaa(0.0), hingeKaaCommit(0.0),
  hingeMscale(0.0), hingeMscaleCommit(0.0),
  pmP0(pmP0_), pmAlpha(pmAlpha_),
  postRestoreSync(false), vRestoreValid(false),
  spanLoadMode(spanMode), hasSpanLoad(false), spanLoadDirty(false),
  spanBetaP(0.0), spanBetaPCommit(0.0),
  spanBetaZ(0.0), spanBetaZCommit(0.0)
{
  vRestore[0] = vRestore[1] = vRestore[2] = 0.0;
  hingeKv[0] = hingeKv[1] = hingeKv[2] = 0.0;
  hingeKvCommit[0] = hingeKvCommit[1] = hingeKvCommit[2] = 0.0;
  q0Enh[0] = q0Enh[1] = q0Enh[2] = 0.0;
  for (int i = 0; i < maxNumSections; i++) {
    s0Fixed[i][0] = 0.0;
    s0Fixed[i][1] = 0.0;
  }
  for (int m = 0; m < 2; m++) {
    for (int n = 0; n < 2; n++) {
      hingeHinv[m][n] = 0.0;
      hingeHinvCommit[m][n] = 0.0;
      spanKbbInv[m][n] = 0.0;
      spanKbbInvCommit[m][n] = 0.0;
    }
    for (int a = 0; a < 3; a++) {
      hingeKva[a][m] = 0.0;
      hingeKvaCommit[a][m] = 0.0;
      spanKvb[a][m] = 0.0;
      spanKvbCommit[a][m] = 0.0;
    }
  }
  if (numSections < 0 || numSections > maxNumSections) {
    opserr << "LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d - numSections ("
           << numSections << ") exceeds maxNumSections (" << maxNumSections << ")\n";
    numSections = (numSections < 0) ? 0 : maxNumSections;
  }
  // Ladruno (ADR 32) Tier-2: own a copy of the cohesive hinge material (like the sections)
  if (hingeMat != 0) {
    theHinge = hingeMat->getCopy();
    if (theHinge == 0) {
      opserr << "LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d - failed to copy hinge material\n";
      exit(-1);
    }
    hingeOn = 1;
  }
  if (hingeMatJ != 0) {
    theHingeJ = hingeMatJ->getCopy();
    if (theHingeJ == 0) {
      opserr << "LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d - failed to copy hingeJ material\n";
      exit(-1);
    }
    hingeOn = 1;
    dualHingeMode = 1;
  }
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
    
  if (theSections == 0) {
    opserr << "LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d - failed to allocate section model pointer\n";
    exit(-1);
  }

  for (int i = 0; i < numSections; i++) {
    
    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();
    
    // Check allocation
    if (theSections[i] == 0) {
      opserr << "LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d -- failed to get a copy of section model\n";
      exit(-1);
    }
  }
  
  beamInt = bi.getCopy();
  
  if (beamInt == 0) {
    opserr << "LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d - failed to copy beam integration\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy2d();
  
  if (crdTransf == 0) {
    opserr << "LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d - failed to copy coordinate transformation\n";
    exit(-1);
  }
  
  if (damping)
  {
    theDamping =(*damping).getCopy();
    
    if (!theDamping) {
      opserr << "LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d - failed to copy damping\n";
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
  
  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
}

LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d()
:Element (0, ELE_TAG_LadrunoDispBeamColumn2d),
 numSections(0), theSections(0), crdTransf(0), beamInt(0),
 connectedExternalNodes(2),
  Q(6), q(3), rho(0.0), cMass(0), parameterID(0),
  theDamping(0),
  current_section_lch(0.0), lchMode(0), userLch(0.0), nlGeom(0),
  hingeOn(0), dualHingeMode(0), hingeXiI(0.0), hingeXiJ(1.0),
  theHinge(0), theHingeJ(0),
  hingeJump(0.0), hingeJumpCommit(0.0),
  hingeJumpJ(0.0), hingeJumpJCommit(0.0),
  hingeKaa(0.0), hingeKaaCommit(0.0),
  hingeMscale(0.0), hingeMscaleCommit(0.0),
  pmP0(0.0), pmAlpha(0.5),
  postRestoreSync(false), vRestoreValid(false),
  spanLoadMode(1), hasSpanLoad(false), spanLoadDirty(false),
  spanBetaP(0.0), spanBetaPCommit(0.0),
  spanBetaZ(0.0), spanBetaZCommit(0.0)
{
    vRestore[0] = vRestore[1] = vRestore[2] = 0.0;
    hingeKv[0] = hingeKv[1] = hingeKv[2] = 0.0;
    hingeKvCommit[0] = hingeKvCommit[1] = hingeKvCommit[2] = 0.0;
    q0Enh[0] = q0Enh[1] = q0Enh[2] = 0.0;
    for (int i = 0; i < maxNumSections; i++) {
      s0Fixed[i][0] = 0.0;
      s0Fixed[i][1] = 0.0;
    }
    for (int m = 0; m < 2; m++) {
      for (int n = 0; n < 2; n++) {
        hingeHinv[m][n] = 0.0;
        hingeHinvCommit[m][n] = 0.0;
        spanKbbInv[m][n] = 0.0;
        spanKbbInvCommit[m][n] = 0.0;
      }
      for (int a = 0; a < 3; a++) {
        hingeKva[a][m] = 0.0;
        hingeKvaCommit[a][m] = 0.0;
        spanKvb[a][m] = 0.0;
        spanKvbCommit[a][m] = 0.0;
      }
    }
    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;

    theNodes[0] = 0;
    theNodes[1] = 0;
}

LadrunoDispBeamColumn2d::~LadrunoDispBeamColumn2d()
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

  if (theHinge != 0)   // Ladruno (ADR 32) Tier-2
    delete theHinge;
  if (theHingeJ != 0)
    delete theHingeJ;
}

int
LadrunoDispBeamColumn2d::getNumExternalNodes() const
{
    return 2;
}

const ID&
LadrunoDispBeamColumn2d::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
LadrunoDispBeamColumn2d::getNodePtrs()
{
    return theNodes;
}

int
LadrunoDispBeamColumn2d::getNumDOF()
{
    return 6;
}

void
LadrunoDispBeamColumn2d::setDomain(Domain *theDomain)
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
      opserr << "WARNING LadrunoDispBeamColumn2d (tag: " << this->getTag() << "), node not found in domain" << endln;
      return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    if (dofNd1 != 3 || dofNd2 != 3) {
	opserr << "FATAL ERROR LadrunoDispBeamColumn2d (tag: " << this->getTag()
	       << "), has differing number of DOFs at its nodes" << endln;
	return;
    }

	if (crdTransf->initialize(theNodes[0], theNodes[1])) {
		// Add some error check
	}

  // initialize the damping
  if (theDamping && theDamping->setDomain(theDomain, 3)) {
    opserr << "LadrunoDispBeamColumn2d::setDomain(): Error initializing damping";  
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
LadrunoDispBeamColumn2d::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping)
  {
    if (theDamping) delete theDamping;

    theDamping =(*damping).getCopy();
    
    if (!theDamping) {
      opserr << "LadrunoDispBeamColumn2d::setDamping -- failed to get copy of damping\n";
      return -1;
    }
    if (theDamping->setDomain(theDomain, 3)) {
      opserr << "LadrunoDispBeamColumn2d::setDamping -- Error initializing damping\n";
      return -2;
    }
  }
  
  return 0;
}

int
LadrunoDispBeamColumn2d::commitState()
{
    int retVal = 0;
    postRestoreSync = false;
    vRestoreValid = false;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "LadrunoDispBeamColumn2d::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    if (theDamping) retVal += theDamping->commitState();

    // Ladruno (ADR 32) Tier-2: lock in the converged jump + the cohesive law's irreversible
    // history (kappaMax). The commit is the ONLY place the jump advances, so a rejected step
    // (revertToLastCommit) cannot resurrect/lose a hinge.
    if (hingeOn) {
      hingeJumpCommit = hingeJump;
      hingeJumpJCommit = hingeJumpJ;
      hingeKaaCommit = hingeKaa;
      hingeKvCommit[0] = hingeKv[0];
      hingeKvCommit[1] = hingeKv[1];
      hingeKvCommit[2] = hingeKv[2];
      for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 2; n++) hingeHinvCommit[m][n] = hingeHinv[m][n];
        for (int a = 0; a < 3; a++) hingeKvaCommit[a][m] = hingeKva[a][m];
      }
      hingeMscaleCommit = hingeMscale;
      if (theHinge) retVal += theHinge->commitState();
      if (theHingeJ) retVal += theHingeJ->commitState();
    }
    spanBetaPCommit = spanBetaP;
    spanBetaZCommit = spanBetaZ;
    for (int m = 0; m < 2; m++) {
      for (int n = 0; n < 2; n++) spanKbbInvCommit[m][n] = spanKbbInv[m][n];
      for (int a = 0; a < 3; a++) spanKvbCommit[a][m] = spanKvb[a][m];
    }

    return retVal;
}

int
LadrunoDispBeamColumn2d::revertToLastCommit()
{
    int retVal = 0;
    postRestoreSync = true;
    vRestoreValid = false;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    if (theDamping) retVal += theDamping->revertToLastCommit();

    // Ladruno (ADR 32) Tier-2: restore the jump (and the cohesive law) to the committed
    // state so a rejected step does not leave a stale/open hinge.
    if (hingeOn) {
      hingeJump = hingeJumpCommit;
      hingeJumpJ = hingeJumpJCommit;
      hingeKaa = hingeKaaCommit;
      hingeKv[0] = hingeKvCommit[0];
      hingeKv[1] = hingeKvCommit[1];
      hingeKv[2] = hingeKvCommit[2];
      for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 2; n++) hingeHinv[m][n] = hingeHinvCommit[m][n];
        for (int a = 0; a < 3; a++) hingeKva[a][m] = hingeKvaCommit[a][m];
      }
      hingeMscale = hingeMscaleCommit;
      if (theHinge && !ldbc2dHingeIsImplex(theHinge)) retVal += theHinge->revertToLastCommit();
      if (theHingeJ && !ldbc2dHingeIsImplex(theHingeJ)) retVal += theHingeJ->revertToLastCommit();
    }
    spanBetaP = spanBetaPCommit;
    spanBetaZ = spanBetaZCommit;
    for (int m = 0; m < 2; m++) {
      for (int n = 0; n < 2; n++) spanKbbInv[m][n] = spanKbbInvCommit[m][n];
      for (int a = 0; a < 3; a++) spanKvb[a][m] = spanKvbCommit[a][m];
    }

    if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
      retVal += this->update();
      if (!theDamping)
        this->getResistingForce();
    }

    return retVal;
}

int
LadrunoDispBeamColumn2d::revertToStart()
{
    int retVal = 0;
    postRestoreSync = false;
    vRestoreValid = false;
    spanLoadDirty = false;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    if (theDamping) retVal += theDamping->revertToStart();

    // Ladruno (ADR 32) Tier-2: reset the embedded hinge to a closed, undamaged state
    if (hingeOn) {
      hingeJump = 0.0;
      hingeJumpCommit = 0.0;
      hingeJumpJ = 0.0;
      hingeJumpJCommit = 0.0;
      hingeKaa = 0.0;
      hingeKaaCommit = 0.0;
      hingeKv[0] = hingeKv[1] = hingeKv[2] = 0.0;
      hingeKvCommit[0] = hingeKvCommit[1] = hingeKvCommit[2] = 0.0;
      for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 2; n++) {
          hingeHinv[m][n] = 0.0;
          hingeHinvCommit[m][n] = 0.0;
        }
        for (int a = 0; a < 3; a++) {
          hingeKva[a][m] = 0.0;
          hingeKvaCommit[a][m] = 0.0;
        }
      }
      hingeMscale = 0.0;
      hingeMscaleCommit = 0.0;
      if (theHinge) retVal += theHinge->revertToStart();
      if (theHingeJ) retVal += theHingeJ->revertToStart();
    }
    spanBetaP = 0.0;
    spanBetaPCommit = 0.0;
    spanBetaZ = 0.0;
    spanBetaZCommit = 0.0;
    for (int m = 0; m < 2; m++) {
      for (int n = 0; n < 2; n++) spanKbbInv[m][n] = spanKbbInvCommit[m][n] = 0.0;
      for (int a = 0; a < 3; a++) spanKvb[a][m] = spanKvbCommit[a][m] = 0.0;
    }

    if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
      retVal += this->update();
      if (!theDamping)
        this->getResistingForce();
    }

    return retVal;
}

int
LadrunoDispBeamColumn2d::update(void)
{
  int err = 0;
  spanLoadDirty = false;

  // Ladruno (ADR 32): defensively re-assert the active element so a crack-band
  // material's once-only getCharacteristicLength() read latches THIS element's
  // per-IP length even if some path strained the section outside Domain::update().
  ops_TheActiveElement = this;

  // Update the transformation
  crdTransf->update();

  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();
  if (postRestoreSync) {
    if (!vRestoreValid) {
      vRestore[0] = v(0); vRestore[1] = v(1); vRestore[2] = v(2);
      vRestoreValid = true;
    } else if (v(0) != vRestore[0] || v(1) != vRestore[1] || v(2) != vRestore[2]) {
      postRestoreSync = false;
    }
  }

  double L = crdTransf->getInitialLength();
  if (L <= 0.0)
    return -1;
  double oneOverL = 1.0/L;

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Ladruno (ADR 32) Tier-2: when an embedded hinge is active, the bulk section strains
  // include the -alpha/L enhancement and alpha is converged by an inner Newton (which
  // sets the section deformations itself). GATED so the no-hinge path below is untouched.
  if (hingeOn) {
    if (dualHingeMode || theHingeJ != 0)
      err = solveDualHingeJump(v, L);
    else
      err = solveHingeJump(v, L);
    if (err != 0)
      return err;
  } else if (spanLoadMode && hasSpanLoad && !nlGeom) {
    err = solveSpanLoadBubble(v, L);
    if (err != 0)
      return err;
  } else {
    spanBetaZ = 0.0;
    spanBetaP = 0.0;
    for (int m = 0; m < 2; m++)
      for (int n = 0; n < 2; n++)
        spanKbbInv[m][n] = 0.0;
    for (int a = 0; a < 3; a++) {
      spanKvb[a][0] = spanKvb[a][1] = 0.0;
    }

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();

      Vector e(workArea, order);

      double xi6 = 6.0*xi[i];
      // Ladruno (ADR 32) Stage-1: transverse slope theta(zeta) for the ½θ² bowing term
      double zeta = xi[i];
      double theta = (3.0*zeta*zeta - 4.0*zeta + 1.0)*v(1) + (3.0*zeta*zeta - 2.0*zeta)*v(2);

      int j;
      for (j = 0; j < order; j++) {
        switch(code(j)) {
        case SECTION_RESPONSE_P:
          // Ladruno (ADR 32): nlGeom adds the ½θ² membrane/bowing strain (DispBeamColumnNL2d)
          e(j) = oneOverL*v(0) + (nlGeom ? 0.5*theta*theta : 0.0); break;
        case SECTION_RESPONSE_MZ:
          e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)); break;
        default:
          e(j) = 0.0; break;
        }
      }

      if (lchMode == 1)        current_section_lch = L;          // -lch element (debug)
      else if (lchMode == 2)   current_section_lch = userLch;    // -lch <value>
      else                     current_section_lch = wt[i]*L;    // -lch ip (default)

      // Set the section deformations
      err += theSections[i]->setTrialSectionDeformation(e);
    }

    if (err != 0) {
      LadrunoOutputGuard og;
      opserr << "LadrunoDispBeamColumn2d::update() - failed setTrialSectionDeformations()\n";
      return err;
    }
  }

  // Keep basic force q synchronized at the end of update() across all branches
  q.Zero();
  for (int i = 0; i < numSections; i++) {
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    double xi6 = 6.0 * xi[i];
    const Vector &s = theSections[i]->getStressResultant();
    double zeta = xi[i];
    double c1nl = 3.0 * zeta * zeta - 4.0 * zeta + 1.0;
    double c2nl = 3.0 * zeta * zeta - 2.0 * zeta;
    double thetaNL = c1nl * v(1) + c2nl * v(2);
    for (int j = 0; j < order; j++) {
      double si = s(j) * wt[i];
      switch (code(j)) {
      case SECTION_RESPONSE_P:
        q(0) += si;
        break;
      case SECTION_RESPONSE_MZ:
        q(1) += (xi6 - 4.0) * si;
        q(2) += (xi6 - 2.0) * si;
        if (nlGeom) {
          for (int jj = 0; jj < order; jj++) {
            if (code(jj) == SECTION_RESPONSE_P) {
              q(1) += c1nl * thetaNL * s(jj) * wt[i] * L;
              q(2) += c2nl * thetaNL * s(jj) * wt[i] * L;
            }
          }
        }
        break;
      default:
        break;
      }
    }
  }
  if (!hingeOn && !nlGeom && spanLoadMode && hasSpanLoad) {
    q(0) += q0[0] - q0Enh[0];
    q(1) += q0[1] - q0Enh[1];
    q(2) += q0[2] - q0Enh[2];
  } else {
    q(0) += q0[0];
    q(1) += q0[1];
    q(2) += q0[2];
  }

  return 0;
}

int
LadrunoDispBeamColumn2d::solveSpanLoadBubble(const Vector &v, double L)
{
  if (L <= 0.0)
    return -1;
  double oneOverL = 1.0 / L;

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  double phi[maxNumSections];
  double G11 = 0.0, G12 = 0.0, G22 = 0.0, g1 = 0.0, g2 = 0.0;
  for (int i = 0; i < numSections; i++) {
    double b1 = 6.0 * xi[i] - 4.0;
    double b2 = 6.0 * xi[i] - 2.0;
    double p0i = 1.0 - 6.0 * xi[i] + 6.0 * xi[i] * xi[i];
    phi[i] = p0i;
    G11 += wt[i] * b1 * b1;
    G12 += wt[i] * b1 * b2;
    G22 += wt[i] * b2 * b2;
    g1  += wt[i] * b1 * p0i;
    g2  += wt[i] * b2 * p0i;
  }
  double detG = G11 * G22 - G12 * G12;
  double phiNorm2 = 0.0;
  if (numSections >= 3 && detG > 1.0e-12) {
    double c1 = ( G22 * g1 - G12 * g2) / detG;
    double c2 = (-G12 * g1 + G11 * g2) / detG;
    for (int i = 0; i < numSections; i++) {
      phi[i] -= c1 * (6.0 * xi[i] - 4.0) + c2 * (6.0 * xi[i] - 2.0);
      phiNorm2 += wt[i] * phi[i] * phi[i];
    }
  } else {
    for (int i = 0; i < numSections; i++)
      phi[i] = 0.0;
  }

  double e0P[maxNumSections], e0Mz[maxNumSections];
  double mScale = 0.0, pScale = 0.0;
  bool hasAnyP = false, hasAnyMZ = false;
  double Kbb0_PP = 0.0, Kbb0_ZZ = 0.0;
  for (int i = 0; i < numSections; i++) {
    e0P[i] = 0.0;
    e0Mz[i] = 0.0;
    pScale += fabs(s0Fixed[i][0]);
    mScale += fabs(s0Fixed[i][1]);
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    const Matrix &ks0 = theSections[i]->getInitialTangent();
    int idxP = -1, idxMZ = -1;
    for (int j = 0; j < order; j++) {
      if (code(j) == SECTION_RESPONSE_P)  { idxP  = j; hasAnyP  = true; }
      if (code(j) == SECTION_RESPONSE_MZ) { idxMZ = j; hasAnyMZ = true; }
    }
    double kPP = (idxP >= 0) ? ks0(idxP, idxP) : 0.0;
    double kMM = (idxMZ >= 0) ? ks0(idxMZ, idxMZ) : 0.0;
    double kPM = (idxP >= 0 && idxMZ >= 0) ? 0.5 * (ks0(idxP, idxMZ) + ks0(idxMZ, idxP)) : 0.0;
    if (idxP >= 0)  Kbb0_PP += oneOverL * wt[i] * phi[i] * phi[i] * fabs(kPP);
    if (idxMZ >= 0) Kbb0_ZZ += oneOverL * wt[i] * phi[i] * phi[i] * fabs(kMM);
    double det0 = kPP * kMM - kPM * kPM;
    if (idxP >= 0 && idxMZ >= 0 && det0 > 1.0e-12 * fabs(kPP * kMM) && det0 > 1.0e-20) {
      e0P[i]  = ( kMM * s0Fixed[i][0] - kPM * s0Fixed[i][1]) / det0;
      e0Mz[i] = (-kPM * s0Fixed[i][0] + kPP * s0Fixed[i][1]) / det0;
    } else {
      if (idxP >= 0 && fabs(kPP) > 1.0e-20)  e0P[i]  = s0Fixed[i][0] / kPP;
      if (idxMZ >= 0 && fabs(kMM) > 1.0e-20) e0Mz[i] = s0Fixed[i][1] / kMM;
    }
  }
  if (mScale * oneOverL > pScale)
    pScale = mScale * oneOverL;

  const bool useBubble = (phiNorm2 > 1.0e-8) && (hasAnyP || hasAnyMZ);
  double bp = (useBubble && hasAnyP)  ? spanBetaPCommit : 0.0;
  double bz = (useBubble && hasAnyMZ) ? spanBetaZCommit : 0.0;
  const int maxIter = useBubble ? 25 : 1;
  double damp = 1.0;
  double rPrev = 1.0e300;
  double Kvb[3][2] = {{0.0}};
  double invKbb[2][2] = {{0.0}};

  for (int iter = 0; iter < maxIter; iter++) {
    int err = 0;
    double rp = 0.0, rz = 0.0;
    double Kbb_PP = 0.0, Kbb_ZZ = 0.0, Kbb_PZ = 0.0;
    for (int a = 0; a < 3; a++) {
      Kvb[a][0] = 0.0;
      Kvb[a][1] = 0.0;
    }

    for (int i = 0; i < numSections; i++) {
      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      double xi6 = 6.0 * xi[i];

      Vector e(workArea, order);
      for (int j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          e(j) = oneOverL * (v(0) + phi[i] * bp) + e0P[i];
          break;
        case SECTION_RESPONSE_MZ:
          e(j) = oneOverL * ((xi6 - 4.0) * v(1) + (xi6 - 2.0) * v(2) + phi[i] * bz) + e0Mz[i];
          break;
        default:
          e(j) = 0.0;
          break;
        }
      }

      if (lchMode == 1)      current_section_lch = L;
      else if (lchMode == 2) current_section_lch = userLch;
      else                   current_section_lch = wt[i] * L;

      err += theSections[i]->setTrialSectionDeformation(e);

      if (useBubble) {
        const Vector &s = theSections[i]->getStressResultant();
        const Matrix &ks = theSections[i]->getSectionTangent();
        int idxP = -1, idxMZ = -1;
        for (int j = 0; j < order; j++) {
          if (code(j) == SECTION_RESPONSE_P)  idxP  = j;
          if (code(j) == SECTION_RESPONSE_MZ) idxMZ = j;
        }
        double wPhi = oneOverL * wt[i] * phi[i];
        if (idxP >= 0) {
          double sP = s(idxP);
          if (fabs(sP) > pScale) pScale = fabs(sP);
          rp += wt[i] * phi[i] * (sP - s0Fixed[i][0]);
          double EA = ks(idxP, idxP);
          Kbb_PP += wPhi * phi[i] * EA;
          Kvb[0][0] += wPhi * EA;
          if (idxMZ >= 0) {
            double kPM = 0.5 * (ks(idxP, idxMZ) + ks(idxMZ, idxP));
            Kvb[1][0] += wPhi * (xi6 - 4.0) * kPM;
            Kvb[2][0] += wPhi * (xi6 - 2.0) * kPM;
          }
        }
        if (idxMZ >= 0) {
          double sMz = s(idxMZ);
          if (fabs(sMz) > mScale) mScale = fabs(sMz);
          if (mScale * oneOverL > pScale) pScale = mScale * oneOverL;
          rz += wt[i] * phi[i] * (sMz - s0Fixed[i][1]);
          double EI = ks(idxMZ, idxMZ);
          Kbb_ZZ += wPhi * phi[i] * EI;
          Kvb[1][1] += wPhi * (xi6 - 4.0) * EI;
          Kvb[2][1] += wPhi * (xi6 - 2.0) * EI;
          if (idxP >= 0) {
            double kPM = 0.5 * (ks(idxP, idxMZ) + ks(idxMZ, idxP));
            Kvb[0][1] += wPhi * kPM;
            Kbb_PZ += wPhi * phi[i] * kPM;
          }
        }
      }
    }

    if (err != 0) {
      LadrunoOutputGuard og;
      opserr << "LadrunoDispBeamColumn2d::solveSpanLoadBubble() - failed setTrialSectionDeformations()\n";
      return err;
    }

    if (!useBubble) {
      spanBetaP = 0.0;
      spanBetaZ = 0.0;
      for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 2; n++) spanKbbInv[m][n] = 0.0;
        for (int a = 0; a < 3; a++) spanKvb[a][m] = 0.0;
      }
      return 0;
    }

    // Diagonal-equilibrated 2x2 guarded inverse so axial and flexural units never cross-floor
    double dP = (hasAnyP  && Kbb0_PP > 1.0e-30) ? sqrt(Kbb0_PP) : 1.0;
    double dZ = (hasAnyMZ && Kbb0_ZZ > 1.0e-30) ? sqrt(Kbb0_ZZ) : 1.0;
    double sPP = hasAnyP  ? (Kbb_PP / (dP * dP)) : 1.0;
    double sZZ = hasAnyMZ ? (Kbb_ZZ / (dZ * dZ)) : 1.0;
    double sPZ = (hasAnyP && hasAnyMZ) ? (Kbb_PZ / (dP * dZ)) : 0.0;
    double sInv[2][2] = {{0.0}};
    invertSym2x2Guarded2d(sPP, sZZ, sPZ, sInv);
    invKbb[0][0] = hasAnyP  ? (sInv[0][0] / (dP * dP)) : 0.0;
    invKbb[1][1] = hasAnyMZ ? (sInv[1][1] / (dZ * dZ)) : 0.0;
    invKbb[0][1] = (hasAnyP && hasAnyMZ) ? (sInv[0][1] / (dP * dZ)) : 0.0;
    invKbb[1][0] = invKbb[0][1];

    bool okP = !hasAnyP  || (fabs(rp) <= 1.0e-11 * pScale + 1.0e-12);
    bool okZ = !hasAnyMZ || (fabs(rz) <= 1.0e-11 * mScale + 1.0e-12);
    if (okP && okZ) {
      spanBetaP = bp;
      spanBetaZ = bz;
      for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 2; n++) spanKbbInv[m][n] = invKbb[m][n];
        for (int a = 0; a < 3; a++) spanKvb[a][m] = Kvb[a][m];
      }
      if (postRestoreSync) {
        for (int m = 0; m < 2; m++) {
          for (int n = 0; n < 2; n++) spanKbbInvCommit[m][n] = invKbb[m][n];
          for (int a = 0; a < 3; a++) spanKvbCommit[a][m] = Kvb[a][m];
        }
      }
      return 0;
    }

    double rnP = hasAnyP  ? (rp / (pScale + 1.0e-12)) : 0.0;
    double rnZ = hasAnyMZ ? (rz / (mScale + 1.0e-12)) : 0.0;
    double rNorm = sqrt(rnP * rnP + rnZ * rnZ);
    if (rNorm > rPrev) damp *= 0.5;
    else               damp = (damp * 1.5 < 1.0) ? damp * 1.5 : 1.0;
    rPrev = rNorm;

    if (hasAnyP)  bp -= damp * (invKbb[0][0] * rp + invKbb[0][1] * rz);
    if (hasAnyMZ) bz -= damp * (invKbb[1][0] * rp + invKbb[1][1] * rz);
  }

  spanBetaP = bp;
  spanBetaZ = bz;
  for (int m = 0; m < 2; m++) {
    for (int n = 0; n < 2; n++) spanKbbInv[m][n] = invKbb[m][n];
    for (int a = 0; a < 3; a++) spanKvb[a][m] = Kvb[a][m];
  }
  return 0;
}

int
LadrunoDispBeamColumn2d::solveHingeJump(const Vector &v, double L)
{
  if (L <= 0.0)
    return -1;
  double oneOverL = 1.0/L;
  const double Gbar = -oneOverL;   // bounded incompatible-mode enhancement

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  double alpha = hingeJumpCommit;   // warm start from the committed jump (path-correct)
  LadrunoCohesiveHinge *cohHinge = dynamic_cast<LadrunoCohesiveHinge*>(theHinge);
  const bool freezeImplex = postRestoreSync && ldbc2dHingeIsImplex(theHinge);

  const int maxIter = 50;
  double Kva[3] = {0.0, 0.0, 0.0};
  double Kaa = 0.0;
  double damp = 1.0;
  double hnormPrev = 1.0e300;

  for (int iter = 0; iter < maxIter; iter++) {

    int err = 0;
    double sumM = 0.0;     // sum_k wt_k * M_sec,k(MZ)
    double sumEI = 0.0;    // sum_k wt_k * EI_k(MZ,MZ)
    double sumP = 0.0;     // sum_k wt_k * N_sec,k(P)
    Kva[0] = Kva[1] = Kva[2] = 0.0;

    // --- set every IP's bulk-enhanced section strain at the current alpha ---
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
          // bulk curvature = B*v + Gbar*alpha  (the -alpha/L unloading offset)
          e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)) + Gbar*alpha; break;
        default:
          e(j) = 0.0; break;
        }
      }

      // Tier-1 lch latch (identical to the no-hinge update loop)
      if (lchMode == 1)       current_section_lch = L;
      else if (lchMode == 2)  current_section_lch = userLch;
      else                    current_section_lch = wt[i]*L;

      err += theSections[i]->setTrialSectionDeformation(e);

      // accumulate residual / condensation terms from this IP
      const Vector &s  = theSections[i]->getStressResultant();
      const Matrix &ks = theSections[i]->getSectionTangent();
      int idxMZ = -1, idxP = -1;
      for (int j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_MZ) idxMZ = j;
        if (code(j) == SECTION_RESPONSE_P)  idxP  = j;
      }
      if (idxP >= 0)
        sumP += wt[i]*s(idxP);
      if (idxMZ >= 0) {
        double EI = ks(idxMZ, idxMZ);
        sumM  += wt[i]*s(idxMZ);
        sumEI += wt[i]*EI;
        Kva[1] += -oneOverL*wt[i]*(xi6-4.0)*EI;
        Kva[2] += -oneOverL*wt[i]*(xi6-2.0)*EI;
        if (idxP >= 0)
          Kva[0] += -oneOverL*wt[i]*ks(idxP, idxMZ);   // axial-flexural coupling
      }
    }

    if (err != 0) {
      LadrunoOutputGuard og;
      opserr << "LadrunoDispBeamColumn2d::solveHingeJump() - failed setTrialSectionDeformations()\n";
      return err;
    }

    sumP += q0[0];
    q(0) = sumP;
    sumM += -(1.0 - hingeXiI)*q0[1] + hingeXiI*q0[2];

    // --- cohesive law at the current jump: M_coh(alpha), dM_coh/dalpha ---
    if (cohHinge != 0)
      cohHinge->setAxialForce(sumP);
    double fPM = (cohHinge == 0 && pmP0 > 0.0) ? ladrunoDisp2dScalePM(sumP, pmP0, pmAlpha) : 1.0;
    if (!freezeImplex)
      theHinge->setTrialStrain(alpha);
    double Mcoh = fPM * theHinge->getStress();
    double Kcoh = fPM * theHinge->getTangent();

    double h = freezeImplex ? 0.0 : (-sumM + Mcoh); // enhancement residual
    Kaa      = oneOverL*sumEI + Kcoh;        // (1/L) sum wt*EI + dM_coh/dalpha

    // GUARDED reciprocal: floor |K_aa| against the positive bulk term so the sign-
    // discontinuity at activation / the indefinite softening branch cannot divide by ~0.
    double KaaFloor = 1.0e-8*(fabs(oneOverL*sumEI) + fabs(Kcoh)) + 1.0e-300;
    double KaaSolve = Kaa;
    if (fabs(KaaSolve) < KaaFloor)
      KaaSolve = (KaaSolve < 0.0) ? -KaaFloor : KaaFloor;

    double Mhere = fabs(Mcoh) + fabs(sumM);
    if (Mhere > hingeMscale) hingeMscale = Mhere;
    if (fabs(h) <= 1.0e-11*hingeMscale + 1.0e-12) {
      hingeKaa  = KaaSolve;                  // guarded, cached for the condensed tangent
      hingeKv[0] = Kva[0]; hingeKv[1] = Kva[1]; hingeKv[2] = Kva[2];
      hingeJump = alpha;
      if (postRestoreSync) {
        hingeKaaCommit = hingeKaa;
        for (int m = 0; m < 3; m++) hingeKvCommit[m] = hingeKv[m];
      }
      return 0;
    }

    double hnorm = fabs(h) / (hingeMscale + 1.0e-12);
    if (hnorm > hnormPrev) damp *= 0.5;
    else                   damp = (damp*1.5 < 1.0) ? damp*1.5 : 1.0;
    hnormPrev = hnorm;

    alpha -= damp*(h/KaaSolve);              // guarded Newton step on the scalar jump
  }

  // not converged within maxIter: cache the last state and let the global solver cut the
  // step. PR-2a's monotonic single-element gates essentially never reach this.
  hingeKaa  = (fabs(Kaa) < 1.0e-300) ? 1.0e-300 : Kaa;
  for (int m = 0; m < 3; m++) hingeKv[m] = Kva[m];
  hingeJump = alpha;
  {
    LadrunoOutputGuard og;
    opserr << "WARNING LadrunoDispBeamColumn2d (tag " << this->getTag()
           << "): embedded-hinge inner Newton did not converge in " << maxIter << " iters\n";
  }
  return 0;
}

// Ladruno: Dual end-hinge (or localized Hermite mode) solver with 2x2 static condensation.
// Uses Hermite end-rotation modes:
//   G_k(xi) = -[(1 - xi_k)*(6*xi - 4) + xi_k*(6*xi - 2)] / L
// At xi_I = 0, G_I(xi) = -(6*xi - 4)/L = -b_1(xi) (work-conjugate to basic end-I moment q_1);
// at xi_J = 1, G_J(xi) = -(6*xi - 2)/L = -b_2(xi) (work-conjugate to basic end-J moment q_2).
// Under double-curvature seismic bending, end I and end J carry independent rotation jumps
// alpha_I and alpha_J simultaneously, solved via a coupled 2x2 inner Newton and condensed
// onto the 3x3 basic stiffness matrix.
int
LadrunoDispBeamColumn2d::solveDualHingeJump(const Vector &v, double L)
{
  if (L <= 0.0)
    return -1;
  double oneOverL = 1.0/L;

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  const bool hasI = (theHinge != 0);
  const bool hasJ = (theHingeJ != 0);
  const bool freezeI = postRestoreSync && hasI && ldbc2dHingeIsImplex(theHinge);
  const bool freezeJ = postRestoreSync && hasJ && ldbc2dHingeIsImplex(theHingeJ);
  const double betaI  = 1.0 - hingeXiI;
  const double gammaI = hingeXiI;
  const double betaJ  = 1.0 - hingeXiJ;
  const double gammaJ = hingeXiJ;

  double alphaI = hasI ? hingeJumpCommit  : 0.0;
  double alphaJ = hasJ ? hingeJumpJCommit : 0.0;

  LadrunoCohesiveHinge *cohI = hasI ? dynamic_cast<LadrunoCohesiveHinge*>(theHinge)  : 0;
  LadrunoCohesiveHinge *cohJ = hasJ ? dynamic_cast<LadrunoCohesiveHinge*>(theHingeJ) : 0;

  const int maxIter = 50;
  double damp = 1.0;
  double hnormPrev = 1.0e300;
  double Hinv[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double Kva[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

  for (int iter = 0; iter < maxIter; iter++) {
    int err = 0;
    double qsec[3] = {0.0, 0.0, 0.0};
    double Kbulk[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    double v1tilde = v(1) - (hasI ? betaI*alphaI : 0.0) - (hasJ ? betaJ*alphaJ : 0.0);
    double v2tilde = v(2) - (hasI ? gammaI*alphaI : 0.0) - (hasJ ? gammaJ*alphaJ : 0.0);

    for (int i = 0; i < numSections; i++) {
      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      double xi6 = 6.0*xi[i];
      double b1 = xi6 - 4.0;
      double b2 = xi6 - 2.0;

      Vector e(workArea, order);
      for (int j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          e(j) = oneOverL*v(0); break;
        case SECTION_RESPONSE_MZ:
          e(j) = oneOverL*(b1*v1tilde + b2*v2tilde); break;
        default:
          e(j) = 0.0; break;
        }
      }

      if (lchMode == 1)       current_section_lch = L;
      else if (lchMode == 2)  current_section_lch = userLch;
      else                    current_section_lch = wt[i]*L;

      err += theSections[i]->setTrialSectionDeformation(e);

      const Vector &s  = theSections[i]->getStressResultant();
      const Matrix &ks = theSections[i]->getSectionTangent();
      double wti = wt[i]*oneOverL;

      for (int j = 0; j < order; j++) {
        double si = s(j)*wt[i];
        if (code(j) == SECTION_RESPONSE_P) {
          qsec[0] += si;
        } else if (code(j) == SECTION_RESPONSE_MZ) {
          qsec[1] += b1*si;
          qsec[2] += b2*si;
        }
        for (int k = 0; k < order; k++) {
          double kjk = ks(j, k)*wti;
          int rJ = (code(j) == SECTION_RESPONSE_P) ? 0 : ((code(j) == SECTION_RESPONSE_MZ) ? 1 : -1);
          int rK = (code(k) == SECTION_RESPONSE_P) ? 0 : ((code(k) == SECTION_RESPONSE_MZ) ? 1 : -1);
          if (rJ == 0 && rK == 0) {
            Kbulk[0][0] += kjk;
          } else if (rJ == 0 && rK == 1) {
            Kbulk[0][1] += kjk*b1;
            Kbulk[0][2] += kjk*b2;
          } else if (rJ == 1 && rK == 0) {
            Kbulk[1][0] += b1*kjk;
            Kbulk[2][0] += b2*kjk;
          } else if (rJ == 1 && rK == 1) {
            Kbulk[1][1] += b1*kjk*b1;
            Kbulk[1][2] += b1*kjk*b2;
            Kbulk[2][1] += b2*kjk*b1;
            Kbulk[2][2] += b2*kjk*b2;
          }
        }
      }
    }

    if (err != 0) {
      LadrunoOutputGuard og;
      opserr << "LadrunoDispBeamColumn2d::solveDualHingeJump() - failed setTrialSectionDeformations()\n";
      return err;
    }

    qsec[0] += q0[0];
    qsec[1] += q0[1];
    qsec[2] += q0[2];
    q(0) = qsec[0];

    if (cohI != 0) cohI->setAxialForce(qsec[0]);
    if (cohJ != 0) cohJ->setAxialForce(qsec[0]);
    double fPMI = (cohI == 0 && pmP0 > 0.0) ? ladrunoDisp2dScalePM(qsec[0], pmP0, pmAlpha) : 1.0;
    double fPMJ = (cohJ == 0 && pmP0 > 0.0) ? ladrunoDisp2dScalePM(qsec[0], pmP0, pmAlpha) : 1.0;

    double McohI = 0.0, KcohI = 0.0, hI = 0.0;
    if (hasI) {
      if (!freezeI)
        theHinge->setTrialStrain(alphaI);
      McohI = fPMI * theHinge->getStress();
      KcohI = fPMI * theHinge->getTangent();
      double qBulkI = betaI*qsec[1] + gammaI*qsec[2];
      hI = freezeI ? 0.0 : (-qBulkI + McohI);
    }

    double McohJ = 0.0, KcohJ = 0.0, hJ = 0.0;
    if (hasJ) {
      if (!freezeJ)
        theHingeJ->setTrialStrain(alphaJ);
      McohJ = fPMJ * theHingeJ->getStress();
      KcohJ = fPMJ * theHingeJ->getTangent();
      double qBulkJ = betaJ*qsec[1] + gammaJ*qsec[2];
      hJ = freezeJ ? 0.0 : (-qBulkJ + McohJ);
    }

    // Build Kva = Kbulk * W (3x2) where W_10 = -betaI, W_20 = -gammaI, W_11 = -betaJ, W_21 = -gammaJ
    for (int a = 0; a < 3; a++) {
      Kva[a][0] = hasI ? -(Kbulk[a][1]*betaI + Kbulk[a][2]*gammaI) : 0.0;
      Kva[a][1] = hasJ ? -(Kbulk[a][1]*betaJ + Kbulk[a][2]*gammaJ) : 0.0;
    }

    // Build 2x2 H = W^T * Kbulk * W + diag(KcohI, KcohJ)
    double H00 = hasI ? (-(betaI*Kva[1][0] + gammaI*Kva[2][0]) + KcohI) : 1.0;
    double H11 = hasJ ? (-(betaJ*Kva[1][1] + gammaJ*Kva[2][1]) + KcohJ) : 1.0;
    double H01 = (hasI && hasJ) ? -0.5*((betaI*Kva[1][1] + gammaI*Kva[2][1]) +
                                        (betaJ*Kva[1][0] + gammaJ*Kva[2][0])) : 0.0;

    invertSym2x2Guarded2d(H00, H11, H01, Hinv);
    if (!hasI) { Hinv[0][0] = 0.0; Hinv[0][1] = 0.0; Hinv[1][0] = 0.0; }
    if (!hasJ) { Hinv[1][1] = 0.0; Hinv[0][1] = 0.0; Hinv[1][0] = 0.0; }

    double Mhere = fabs(McohI) + fabs(McohJ) + fabs(qsec[1]) + fabs(qsec[2]);
    if (Mhere > hingeMscale) hingeMscale = Mhere;

    double hMag = fabs(hI) + fabs(hJ);
    if (hMag <= 1.0e-11*hingeMscale + 1.0e-12) {
      for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 2; n++) hingeHinv[m][n] = Hinv[m][n];
        for (int a = 0; a < 3; a++) hingeKva[a][m] = Kva[a][m];
      }
      hingeJump  = alphaI;
      hingeJumpJ = alphaJ;
      if (postRestoreSync) {
        for (int m = 0; m < 2; m++) {
          for (int n = 0; n < 2; n++) hingeHinvCommit[m][n] = hingeHinv[m][n];
          for (int a = 0; a < 3; a++) hingeKvaCommit[a][m] = hingeKva[a][m];
        }
      }
      return 0;
    }

    double hnorm = hMag / (hingeMscale + 1.0e-12);
    if (hnorm > hnormPrev) damp *= 0.5;
    else                   damp = (damp*1.5 < 1.0) ? damp*1.5 : 1.0;
    hnormPrev = hnorm;

    if (hasI && !freezeI) alphaI -= damp*(Hinv[0][0]*hI + Hinv[0][1]*hJ);
    if (hasJ && !freezeJ) alphaJ -= damp*(Hinv[1][0]*hI + Hinv[1][1]*hJ);
  }

  for (int m = 0; m < 2; m++) {
    for (int n = 0; n < 2; n++) hingeHinv[m][n] = Hinv[m][n];
    for (int a = 0; a < 3; a++) hingeKva[a][m] = Kva[a][m];
  }
  hingeJump  = alphaI;
  hingeJumpJ = alphaJ;
  {
    LadrunoOutputGuard og;
    opserr << "WARNING LadrunoDispBeamColumn2d (tag " << this->getTag()
           << "): dual embedded-hinge inner Newton did not converge in " << maxIter << " iters\n";
  }
  return 0;
}

double
LadrunoDispBeamColumn2d::getCharacteristicLength(void)
{
  if (lchMode == 1)
    return crdTransf ? crdTransf->getInitialLength() : Element::getCharacteristicLength();
  if (lchMode == 2 && userLch > 0.0)
    return userLch;
  if (current_section_lch > 0.0)
    return current_section_lch;
  return Element::getCharacteristicLength();
}

void
LadrunoDispBeamColumn2d::getBasicStiff(Matrix &kb, int initial)
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

  // Ladruno (ADR 32) Stage-1: basic trial disp needed for the NL bowing terms.
  const Vector &v = crdTransf->getBasicTrialDisp();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    double xi6 = 6.0*xi[i];
    int j, k;
    double tmp;

    if (nlGeom && !initial) {
      // Ladruno (ADR 32) Stage-1: ½θ² nonlinear-geometry tangent (verbatim from
      // DispBeamColumnNL2d::getBasicStiff): material + geometric (axial-force) terms.
      double zeta = xi[i];
      double c1 = 3.0*zeta*zeta - 4.0*zeta + 1.0;
      double c2 = 3.0*zeta*zeta - 2.0*zeta;
      double theta = c1*v(1) + c2*v(2);

      const Matrix &ks = theSections[i]->getSectionTangent();
      const Vector &s  = theSections[i]->getStressResultant();
      double wti = wt[i]*oneOverL;

      // geometric (axial-force) term: int_0^L N * C'*C dx
      for (j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_P) {
          tmp = s(j)*wt[i]*L;
          kb(1,1) += tmp*c1*c1; kb(2,1) += tmp*c2*c1;
          kb(1,2) += tmp*c1*c2; kb(2,2) += tmp*c2*c2;
        }
      }

      Matrix B(order,3);
      Matrix C(order,3);
      static thread_local Matrix C1(1,3);
      C1.Zero();
      for (j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          B(j,0) = 1.0; C(j,1) = c1; C(j,2) = c2; C1(0,1) = c1; C1(0,2) = c2; break;
        case SECTION_RESPONSE_MZ:
          B(j,1) = xi6-4.0; B(j,2) = xi6-2.0; break;
        default: break;
        }
      }

      // B'*ks*B  +  B'*ks*C*theta  +  C1'*ks*B*theta  +  C1'*ks*C*theta²*L
      kb.addMatrixTripleProduct(1.0, B, ks, wti);
      Matrix kC(order,3);
      kC.addMatrixProduct(0.0, ks, C, 1.0);
      kb.addMatrixTransposeProduct(1.0, B, kC, theta*wt[i]);

      Matrix ks1(1,order);
      static thread_local Matrix ksB(1,3);
      for (j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_P) {
          for (int jj = 0; jj < order; jj++) ks1(0,jj) = ks(j,jj);
          ksB.addMatrixProduct(0.0, ks1, B, 1.0);
          kb.addMatrixTransposeProduct(1.0, C1, ksB, theta*wt[i]);
          ksB.addMatrixProduct(0.0, ks1, C, 1.0);
          kb.addMatrixTransposeProduct(1.0, C1, ksB, theta*theta*wt[i]*L);
        }
      }
      continue;
    }

    // ----- linear basic strain (stock DispBeamColumn2d), also used for initial -----
    Matrix ka(workArea, order, 3);
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
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 3; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 3; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }
  }

  // Ladruno (ADR 32) Tier-2: statically condense the embedded jump(s) out of the basic
  // stiffness BEFORE crdTransf (and for basicStiffness response)
  if (!initial && hingeOn) {
    if (dualHingeMode || theHingeJ != 0) {
      for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
          double c = 0.0;
          for (int m = 0; m < 2; m++)
            for (int n = 0; n < 2; n++)
              c += hingeKva[a][m]*hingeHinv[m][n]*hingeKva[b][n];
          kb(a,b) -= c;
        }
      }
    } else if (hingeKaa != 0.0) {
      double invKaa = 1.0/hingeKaa;
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
          kb(a,b) -= hingeKv[a]*hingeKv[b]*invKaa;
    }
  } else if (!initial && !hingeOn && !nlGeom && spanLoadMode && hasSpanLoad) {
    // Ladruno (WP-200): static condensation of 2-DOF quartic span-load bubble onto kb
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        double c = 0.0;
        for (int m = 0; m < 2; m++)
          for (int n = 0; n < 2; n++)
            c += spanKvb[a][m] * spanKbbInv[m][n] * spanKvb[b][n];
        kb(a, b) -= c;
      }
    }
  }
}

const Matrix&
LadrunoDispBeamColumn2d::getTangentStiff()
{
  if (spanLoadDirty && theNodes[0] != 0 && theNodes[1] != 0) {
    Element *savedActive = ops_TheActiveElement;
    ops_TheActiveElement = this;
    this->update();
    ops_TheActiveElement = savedActive;
  }

  static thread_local Matrix kb(3,3);

  this->getBasicStiff(kb);

  // Zero for integral
  q.Zero();

  double L = crdTransf->getInitialLength();
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

    // Get the section tangent stiffness and stress resultant
    const Vector &s = theSections[i]->getStressResultant();

    // Ladruno (ADR 32) Stage-1: ½θ² bowing contribution to the basic force
    const Vector &vb = crdTransf->getBasicTrialDisp();
    double zeta = xi[i];
    double c1nl = 3.0*zeta*zeta - 4.0*zeta + 1.0;
    double c2nl = 3.0*zeta*zeta - 2.0*zeta;
    double thetaNL = c1nl*vb(1) + c2nl*vb(2);

    double si;
    for (int j = 0; j < order; j++) {
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si; break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si;
	if (nlGeom) {
	  for (int jj = 0; jj < order; jj++) {
	    if (code(jj) == SECTION_RESPONSE_P) {
	      q(1) += c1nl*thetaNL*s(jj)*wt[i]*L;
	      q(2) += c2nl*thetaNL*s(jj)*wt[i]*L;
	    }
	  }
	}
	break;
      default:
	break;
      }
    }

  }

  // Add effects of element loads, q = q(v) + q0 (or quadrature-truncation correction q0 - q0Enh)
  if (!hingeOn && !nlGeom && spanLoadMode && hasSpanLoad) {
    q(0) += q0[0] - q0Enh[0];
    q(1) += q0[1] - q0Enh[1];
    q(2) += q0[2] - q0Enh[2];
  } else {
    q(0) += q0[0];
    q(1) += q0[1];
    q(2) += q0[2];
  }

  // Transform to global stiffness
  K = crdTransf->getGlobalStiffMatrix(kb, q);

  return K;
}

const Matrix&
LadrunoDispBeamColumn2d::getInitialStiff()
{
  static thread_local Matrix kb(3,3);
  this->getBasicStiff(kb, 1);

  if (hingeOn && (dualHingeMode || theHingeJ != 0)) {
    const bool hasI = (theHinge != 0);
    const bool hasJ = (theHingeJ != 0);
    double KcohI0 = hasI ? theHinge->getInitialTangent()  : 0.0;
    double KcohJ0 = hasJ ? theHingeJ->getInitialTangent() : 0.0;
    bool okI = !hasI || (std::isfinite(KcohI0) && KcohI0 > 0.0 && KcohI0 < 1.0e20);
    bool okJ = !hasJ || (std::isfinite(KcohJ0) && KcohJ0 > 0.0 && KcohJ0 < 1.0e20);
    if (okI && okJ && (hasI || hasJ)) {
      const double betaI  = 1.0 - hingeXiI;
      const double gammaI = hingeXiI;
      const double betaJ  = 1.0 - hingeXiJ;
      const double gammaJ = hingeXiJ;
      double Kva0[3][2];
      for (int a = 0; a < 3; a++) {
        Kva0[a][0] = hasI ? -(kb(a,1)*betaI + kb(a,2)*gammaI) : 0.0;
        Kva0[a][1] = hasJ ? -(kb(a,1)*betaJ + kb(a,2)*gammaJ) : 0.0;
      }
      double H00 = hasI ? (-(betaI*Kva0[1][0] + gammaI*Kva0[2][0]) + KcohI0) : 1.0;
      double H11 = hasJ ? (-(betaJ*Kva0[1][1] + gammaJ*Kva0[2][1]) + KcohJ0) : 1.0;
      double H01 = (hasI && hasJ) ? -0.5*((betaI*Kva0[1][1] + gammaI*Kva0[2][1]) +
                                          (betaJ*Kva0[1][0] + gammaJ*Kva0[2][0])) : 0.0;
      double Hinv0[2][2];
      invertSym2x2Guarded2d(H00, H11, H01, Hinv0);
      if (!hasI) { Hinv0[0][0] = 0.0; Hinv0[0][1] = 0.0; Hinv0[1][0] = 0.0; }
      if (!hasJ) { Hinv0[1][1] = 0.0; Hinv0[0][1] = 0.0; Hinv0[1][0] = 0.0; }
      for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
          double c = 0.0;
          for (int m = 0; m < 2; m++)
            for (int n = 0; n < 2; n++)
              c += Kva0[a][m]*Hinv0[m][n]*Kva0[b][n];
          kb(a,b) -= c;
        }
      }
    }
  } else if (hingeOn && theHinge != 0) {
    double Kcoh0 = theHinge->getInitialTangent();
    if (std::isfinite(Kcoh0) && Kcoh0 > 0.0 && Kcoh0 < 1.0e20) {
      double L = crdTransf->getInitialLength();
      if (L > 0.0) {
        double oneOverL = 1.0/L;
        double xi[maxNumSections];
        double wt[maxNumSections];
        beamInt->getSectionLocations(numSections, L, xi);
        beamInt->getSectionWeights(numSections, L, wt);
        double sumEI0 = 0.0;
        double Kva0[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < numSections; i++) {
          int order = theSections[i]->getOrder();
          const ID &code = theSections[i]->getType();
          const Matrix &ks0 = theSections[i]->getInitialTangent();
          double xi6 = 6.0*xi[i];
          int idxMZ = -1, idxP = -1;
          for (int j = 0; j < order; j++) {
            if (code(j) == SECTION_RESPONSE_MZ) idxMZ = j;
            if (code(j) == SECTION_RESPONSE_P)  idxP  = j;
          }
          if (idxMZ >= 0) {
            double EI0 = ks0(idxMZ, idxMZ);
            sumEI0 += wt[i]*EI0;
            Kva0[1] += -oneOverL*wt[i]*(xi6-4.0)*EI0;
            Kva0[2] += -oneOverL*wt[i]*(xi6-2.0)*EI0;
            if (idxP >= 0)
              Kva0[0] += -oneOverL*wt[i]*ks0(idxP, idxMZ);
          }
        }
        double Kaa0 = oneOverL*sumEI0 + Kcoh0;
        double KaaFloor = 1.0e-8*(fabs(oneOverL*sumEI0) + fabs(Kcoh0)) + 1.0e-300;
        if (fabs(Kaa0) < KaaFloor)
          Kaa0 = (Kaa0 < 0.0) ? -KaaFloor : KaaFloor;
        double invKaa0 = 1.0/Kaa0;
        for (int a = 0; a < 3; a++)
          for (int b = 0; b < 3; b++)
            kb(a,b) -= Kva0[a]*Kva0[b]*invKaa0;
      }
    }
  }

  if(theDamping) kb *= theDamping->getStiffnessMultiplier();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);

  return K;
}

const Matrix&
LadrunoDispBeamColumn2d::getMass()
{
  K.Zero();

  if (rho == 0.0)
    return K;
  
  double L = crdTransf->getInitialLength();
  if (cMass == 0)  {
    // lumped mass matrix
    double m = 0.5*rho*L;
    K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;
  } else  {
    // consistent mass matrix
    static thread_local Matrix ml(6,6);
    ml.Zero();
    double m = rho*L/420.0;
    ml(0,0) = ml(3,3) = m*140.0;
    ml(0,3) = ml(3,0) = m*70.0;
    
    ml(1,1) = ml(4,4) = m*156.0;
    ml(1,4) = ml(4,1) = m*54.0;
    ml(2,2) = ml(5,5) = m*4.0*L*L;
    ml(2,5) = ml(5,2) = -m*3.0*L*L;
    ml(1,2) = ml(2,1) = m*22.0*L;
    ml(4,5) = ml(5,4) = -ml(1,2);
    ml(1,5) = ml(5,1) = -m*13.0*L;
    ml(2,4) = ml(4,2) = -ml(1,5);
    
    // transform local mass matrix to global system
    K = crdTransf->getGlobalMatrixFromLocal(ml);
  }
  
  return K;
}

void
LadrunoDispBeamColumn2d::zeroLoad(void)
{
  Q.Zero();
  if (q0[0] != 0.0 || q0[1] != 0.0 || q0[2] != 0.0 || hasSpanLoad != 0) {
    postRestoreSync = false;
    spanLoadDirty = true;
  }

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  
  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  hasSpanLoad = 0;
  q0Enh[0] = 0.0;
  q0Enh[1] = 0.0;
  q0Enh[2] = 0.0;
  for (int i = 0; i < maxNumSections; i++) {
    s0Fixed[i][0] = 0.0;
    s0Fixed[i][1] = 0.0;
  }
  
  return;
}

int 
LadrunoDispBeamColumn2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = crdTransf->getInitialLength();
  postRestoreSync = false;

  double xi[maxNumSections];
  double wt[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  beamInt->getSectionWeights(numSections, L, wt);
  bool addedSpanFixed = false;
  
  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wtLoad = data(0)*loadFactor;  // Transverse (+ve upward)
    double waLoad = data(1)*loadFactor;  // Axial (+ve from node I to J)

    double V = 0.5*wtLoad*L;
    double M = V*L/6.0; // wtLoad*L*L/12
    double P = waLoad*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= V;
    p0[2] -= V;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= M;
    q0[2] += M;

    for (int i = 0; i < numSections; i++) {
      double x = xi[i];
      s0Fixed[i][0] += P * (0.5 - x);
      s0Fixed[i][1] += M * (1.0 - 6.0 * x + 6.0 * x * x);
    }
    addedSpanFixed = true;
  }
  else if (type == LOAD_TAG_BeamUniformMoment) {
    double mz = data(2)*loadFactor;  // About z

    // Reactions in basic system
    p0[1] += mz;
    p0[2] -= mz;

    // Fixed end forces in basic system
    //q0[1] -= 0.0;
    //q0[2] += 0.0;
    spanLoadDirty = true;
  }
  else if (type == LOAD_TAG_Beam2dPartialUniformLoad) {
    // Ladruno: cancellation-free trapezoidal partial load support
    double wya = data(0) * loadFactor;  // Transverse at start
    double wyb = data(1) * loadFactor;  // Transverse at end
    double waa = data(2) * loadFactor;  // Axial at start
    double wab = data(3) * loadFactor;  // Axial at end
    double a = data(4) * L;
    double b = data(5) * L;
    double ell = b - a;
    if (ell <= 0.0 || a < 0.0 || b > L)
      return 0;

    double J0 = 0.5 * ell * (wya + wyb);
    double J1 = (ell * ell / 6.0) * (wya + 2.0 * wyb);
    double J2 = (ell * ell * ell / 12.0) * (wya + 3.0 * wyb);
    double J3 = (ell * ell * ell * ell / 20.0) * (wya + 4.0 * wyb);
    double I0 = J0;
    double I1 = a * J0 + J1;
    double I2 = a * a * J0 + 2.0 * a * J1 + J2;
    double I3 = a * a * a * J0 + 3.0 * a * a * J1 + 3.0 * a * J2 + J3;

    double V2 = I1 / L;
    double V1 = I0 - V2;
    double M1 = -I1 + (2.0 / L) * I2 - (1.0 / (L * L)) * I3;
    double M2 = (1.0 / L) * I2 - (1.0 / (L * L)) * I3;

    double P = 0.5 * ell * (waa + wab);
    double PJ = (a * P + (ell * ell / 6.0) * (waa + 2.0 * wab)) / L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= V1;
    p0[2] -= V2;

    // Fixed end forces in basic system
    q0[0] -= PJ;
    q0[1] += M1;
    q0[2] += M2;

    for (int i = 0; i < numSections; i++) {
      double x = xi[i] * L;
      double Ncant = 0.0;
      double Mcant = 0.0;
      if (x > a && x < b) {
        double u = x - a;
        double s = u / ell;
        Ncant = u * (waa + 0.5 * s * (wab - waa));
        Mcant = u * u * (0.5 * wya + (1.0 / 6.0) * s * (wyb - wya));
      } else if (x >= b) {
        Ncant = P;
        Mcant = (x - a) * J0 - J1;
      }
      s0Fixed[i][0] += (P - PJ) - Ncant;
      s0Fixed[i][1] += (xi[i] - 1.0) * M1 + xi[i] * M2 - V1 * x + Mcant;
    }
    addedSpanFixed = true;
  }
  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1 = -a * b2 * P * L2;
    double M2 = a2 * b * P * L2;
    q0[1] += M1;
    q0[2] += M2;

    for (int i = 0; i < numSections; i++) {
      double x = xi[i] * L;
      double Ncant = (x > a) ? N : 0.0;
      double Mcant = (x > a) ? P * (x - a) : 0.0;
      s0Fixed[i][0] += N * (1.0 - aOverL) - Ncant;
      s0Fixed[i][1] += (xi[i] - 1.0) * M1 + xi[i] * M2 - V1 * x + Mcant;
    }
    addedSpanFixed = true;
  }
  else {
    opserr << "LadrunoDispBeamColumn2d::LadrunoDispBeamColumn2d -- load type unknown for element with tag: "
	   << this->getTag() << "LadrunoDispBeamColumn2d::addLoad()\n"; 
			    
    return -1;
  }

  if (addedSpanFixed) {
    hasSpanLoad = 1;
    spanLoadDirty = true;
    q0Enh[0] = 0.0;
    q0Enh[1] = 0.0;
    q0Enh[2] = 0.0;
    for (int i = 0; i < numSections; i++) {
      double xi6 = 6.0 * xi[i];
      q0Enh[0] += wt[i] * s0Fixed[i][0];
      q0Enh[1] += wt[i] * (xi6 - 4.0) * s0Fixed[i][1];
      q0Enh[2] += wt[i] * (xi6 - 2.0) * s0Fixed[i][1];
    }
  }

  return 0;
}

int 
LadrunoDispBeamColumn2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0) 
    return 0;
  
  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
    opserr << "LadrunoDispBeamColumn2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }
  
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    Q(0) -= m*Raccel1(0);
    Q(1) -= m*Raccel1(1);
    Q(3) -= m*Raccel2(0);
    Q(4) -= m*Raccel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static thread_local Vector Raccel(6);
    for (int i=0; i<3; i++)  {
      Raccel(i)   = Raccel1(i);
      Raccel(i+3) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }
    
    return 0;
}

const Vector&
LadrunoDispBeamColumn2d::getResistingForce()
{
  if (spanLoadDirty && theNodes[0] != 0 && theNodes[1] != 0) {
    Element *savedActive = ops_TheActiveElement;
    ops_TheActiveElement = this;
    this->update();
    ops_TheActiveElement = savedActive;
  }
  double L = crdTransf->getInitialLength();
  
  double oneOverL = 1.0/L;
  
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
  
    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    
    // Get section stress resultant
    const Vector &s = theSections[i]->getStressResultant();

    // Ladruno (ADR 32) Stage-1: ½θ² bowing contribution to the basic force
    const Vector &vb = crdTransf->getBasicTrialDisp();
    double zeta = xi[i];
    double c1nl = 3.0*zeta*zeta - 4.0*zeta + 1.0;
    double c2nl = 3.0*zeta*zeta - 2.0*zeta;
    double thetaNL = c1nl*vb(1) + c2nl*vb(2);

    // Perform numerical integration on internal force
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));

    double si;
    for (int j = 0; j < order; j++) {
      //si = s(j)*wts(i);
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si; break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si;
	if (nlGeom) {
	  for (int jj = 0; jj < order; jj++) {
	    if (code(jj) == SECTION_RESPONSE_P) {
	      q(1) += c1nl*thetaNL*s(jj)*wt[i]*L;
	      q(2) += c2nl*thetaNL*s(jj)*wt[i]*L;
	    }
	  }
	}
	break;
      default:
	break;
      }
    }
  }
  
  // Add effects of element loads, q = q(v) + q0 (or quadrature-truncation correction q0 - q0Enh)
  if (!hingeOn && !nlGeom && spanLoadMode && hasSpanLoad) {
    q(0) += q0[0] - q0Enh[0];
    q(1) += q0[1] - q0Enh[1];
    q(2) += q0[2] - q0Enh[2];
  } else {
    q(0) += q0[0];
    q(1) += q0[1];
    q(2) += q0[2];
  }

  if (theDamping) theDamping->update(q);

  // Vector for reactions in basic system
  Vector p0Vec(p0, 3);

  P = crdTransf->getGlobalResistingForce(q, p0Vec);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (rho != 0)
    P.addVector(1.0, Q, -1.0);
  
  return P;
}

const Vector &
LadrunoDispBeamColumn2d::getDampingForce(void)
{
  crdTransf->update();

  return crdTransf->getGlobalResistingForce(theDamping->getDampingForce(), Vector(3));
}

const Vector&
LadrunoDispBeamColumn2d::getResistingForceIncInertia()
{
  // SNAPSHOT into a function-local buffer, not the shared static P: betaK
  // Rayleigh re-enters getTangentStiff(), so accumulating into a buffer that
  // path could refill would silently drop inertia (LEDGER_quirks "MUST snapshot
  // the shared static `resid`"). Returning `res` also keeps the result distinct
  // from the P that getResistingForce() returns. Same operations, same order,
  // so results are bit-identical (WP-116).  // Ladruno
  static thread_local Vector res(6);
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
    res(3) += m*accel2(0);
    res(4) += m*accel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static thread_local Vector accel(6);
    for (int i=0; i<3; i++)  {
      accel(i)   = accel1(i);
      accel(i+3) = accel2(i);
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
LadrunoDispBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;
  
  static thread_local Vector data(87);  // Ladruno (ADR 32 + dual end hinges + PM + q0/p0 + spanLoad)
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
    int dbTag = theDamping->getDbTag();
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	      theDamping->setDbTag(dbTag);
	  }
    data(15) = dbTag;
  }

  // Ladruno (ADR 32): regularization-length config + nonlinear-geometry flag
  data(16) = lchMode;
  data(17) = userLch;
  data(18) = nlGeom;

  // Ladruno (ADR 32) Tier-2: embedded-hinge flags + committed jumps + moment scale
  data(19) = hingeOn;
  data(20) = hingeJumpCommit;
  data(21) = hingeMscale;
  data(22) = dualHingeMode;
  data(23) = hingeXiI;
  data(24) = hingeXiJ;
  data(25) = (theHinge != 0) ? 1.0 : 0.0;
  data(26) = (theHingeJ != 0) ? 1.0 : 0.0;
  data(27) = hingeJumpJCommit;
  data(28) = 0.0;
  data(29) = 0.0;
  if (hingeOn && theHinge) {
    data(28) = theHinge->getClassTag();
    int hDbTag = theHinge->getDbTag();
    if (hDbTag == 0) {
      hDbTag = theChannel.getDbTag();
      if (hDbTag != 0) theHinge->setDbTag(hDbTag);
    }
    data(29) = hDbTag;
  }
  data(30) = 0.0;
  data(31) = 0.0;
  if (hingeOn && theHingeJ) {
    data(30) = theHingeJ->getClassTag();
    int hDbTagJ = theHingeJ->getDbTag();
    if (hDbTagJ == 0) {
      hDbTagJ = theChannel.getDbTag();
      if (hDbTagJ != 0) theHingeJ->setDbTag(hDbTagJ);
    }
    data(31) = hDbTagJ;
  }
  data(32) = pmP0;
  data(33) = pmAlpha;
  for (int m = 0; m < 3; m++) {
    data(34 + m) = q0[m];
    data(37 + m) = p0[m];
  }
  data(40) = spanLoadMode;
  data(41) = hasSpanLoad;
  data(42) = spanBetaZCommit;
  for (int m = 0; m < 3; m++)
    data(43 + m) = q0Enh[m];
  for (int s = 0; s < maxNumSections; s++) {
    data(46 + 2 * s)     = s0Fixed[s][0];
    data(46 + 2 * s + 1) = s0Fixed[s][1];
  }
  data(86) = spanBetaPCommit;

  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "LadrunoDispBeamColumn2d::sendSelf() - failed to send data Vector\n";
     return -1;
  }
  
  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "LadrunoDispBeamColumn2d::sendSelf() - failed to send crdTranf\n";
     return -1;
  }      

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoDispBeamColumn2d::sendSelf() - failed to send beamInt\n";
    return -1;
  }      

  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //
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
    opserr << "LadrunoDispBeamColumn2d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn2d::sendSelf() - section " << 
	j << "failed to send itself\n";
      return -1;
    }
  }

  // Ask the Damping to send itself
  if (theDamping && theDamping->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn2d::sendSelf -- could not send Damping\n";
      return -1;
  }

  // Ladruno (ADR 32) Tier-2: send the embedded hinge materials (classTag+dbTag already in data Vector)
  if (hingeOn && theHinge) {
    if (theHinge->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn2d::sendSelf() - failed to send hinge material\n";
      return -1;
    }
  }
  if (hingeOn && theHingeJ) {
    if (theHingeJ->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn2d::sendSelf() - failed to send hingeJ material\n";
      return -1;
    }
  }

  return 0;
}

int
LadrunoDispBeamColumn2d::recvSelf(int commitTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;
  
  static thread_local Vector data(87);  // Ladruno (ADR 32 + dual end hinges + PM + q0/p0 + spanLoad)

  if (theChannel.recvVector(dbTag, commitTag, data) < 0)  {
    opserr << "LadrunoDispBeamColumn2d::recvSelf() - failed to recv data Vector\n";
    return -1;
  }

  this->setTag((int)data(0));
  connectedExternalNodes(0) = (int)data(1);
  connectedExternalNodes(1) = (int)data(2);
  int nSect = (int)data(3);
  if (nSect <= 0 || nSect > maxNumSections) {
    opserr << "LadrunoDispBeamColumn2d::recvSelf() - invalid numSections (" << nSect << ")\n";
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

  // Ladruno (ADR 32): regularization-length config (current_section_lch is transient,
  // recomputed in update(), so it is not serialized).
  lchMode = (int)data(16);
  userLch = data(17);
  nlGeom = (int)data(18);
  current_section_lch = 0.0;

  // Ladruno (ADR 32) Tier-2: embedded-hinge flags + committed jumps + moment scale
  hingeOn = (int)data(19);
  hingeJumpCommit = data(20);
  hingeMscale = data(21);
  dualHingeMode = (int)data(22);
  hingeXiI = data(23);
  hingeXiJ = data(24);
  int hasHingeI = (data(25) != 0.0);
  int hasHingeJ = (data(26) != 0.0);
  hingeJumpJCommit = data(27);
  int hClassTag = (int)data(28);
  int hDbTag = (int)data(29);
  int hClassTagJ = (int)data(30);
  int hDbTagJ = (int)data(31);
  pmP0 = data(32);
  pmAlpha = data(33);
  for (int m = 0; m < 3; m++) {
    q0[m] = data(34 + m);
    p0[m] = data(37 + m);
  }
  spanLoadMode = (int)data(40);
  hasSpanLoad = (int)data(41);
  spanBetaZCommit = data(42);
  spanBetaZ = spanBetaZCommit;
  spanBetaPCommit = data(86);
  spanBetaP = spanBetaPCommit;
  for (int m = 0; m < 3; m++)
    q0Enh[m] = data(43 + m);
  for (int m = 0; m < 2; m++) {
    for (int n = 0; n < 2; n++) spanKbbInv[m][n] = spanKbbInvCommit[m][n] = 0.0;
    for (int a = 0; a < 3; a++) spanKvb[a][m] = spanKvbCommit[a][m] = 0.0;
  }
  for (int s = 0; s < maxNumSections; s++) {
    s0Fixed[s][0] = data(46 + 2 * s);
    s0Fixed[s][1] = data(46 + 2 * s + 1);
  }
  postRestoreSync = true;
  spanLoadDirty = false;
  vRestoreValid = false;
  hingeMscaleCommit = hingeMscale;
  hingeJump = hingeJumpCommit;
  hingeJumpJ = hingeJumpJCommit;
  hingeKaa = 0.0;
  hingeKaaCommit = 0.0;
  hingeKv[0] = hingeKv[1] = hingeKv[2] = 0.0;
  hingeKvCommit[0] = hingeKvCommit[1] = hingeKvCommit[2] = 0.0;
  for (int m = 0; m < 2; m++) {
    for (int n = 0; n < 2; n++) {
      hingeHinv[m][n] = 0.0;
      hingeHinvCommit[m][n] = 0.0;
    }
    for (int a = 0; a < 3; a++) {
      hingeKva[a][m] = 0.0;
      hingeKvaCommit[a][m] = 0.0;
    }
  }

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "LadrunoDispBeamColumn2d::recvSelf() - failed to obtain a CrdTrans object with classTag " <<
	  crdTransfClassTag << endln;
	  return -2;	  
      }
  }
  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "LadrunoDispBeamColumn2d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }      

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
      if (beamInt != 0)
	  delete beamInt;

      beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

      if (beamInt == 0) {
	opserr << "LadrunoDispBeamColumn2d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	  beamIntClassTag << endln;
	exit(-1);
      }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "LadrunoDispBeamColumn2d::sendSelf() - failed to recv beam integration\n";
     return -3;
  }      

  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*nSect);
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "LadrunoDispBeamColumn2d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  //
  // now receive the sections
  //
  
  if (numSections != nSect) {

    // delete the old
    if (numSections != 0) {
      for (int i=0; i<numSections; i++)
	delete theSections[i];
      delete [] theSections;
    }

    // create a new array to hold pointers
    theSections = new SectionForceDeformation *[nSect];
    if (theSections == 0) {
opserr << "LadrunoDispBeamColumn2d::recvSelf() - out of memory creating sections array of size " <<
  nSect << endln;
      return -1;
    }    

    // create a section and recvSelf on it
    numSections = nSect;
    loc = 0;
    
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      theSections[i] = theBroker.getNewSection(sectClassTag);
      if (theSections[i] == 0) {
	opserr << "LadrunoDispBeamColumn2d::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "LadrunoDispBeamColumn2d::recvSelf() - section " << i << " failed to recv itself\n";
	return -1;
      }     
    }

  } else {

    loc = 0;
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;

      // check of correct type
      if (theSections[i]->getClassTag() !=  sectClassTag) {
	// delete the old section[i] and create a new one
	delete theSections[i];
	theSections[i] = theBroker.getNewSection(sectClassTag);
	if (theSections[i] == 0) {
	opserr << "LadrunoDispBeamColumn2d::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "LadrunoDispBeamColumn2d::recvSelf() - section " << i << " failed to recv itself\n";
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
        opserr << "LadrunoDispBeamColumn2d::recvSelf -- could not get a Damping\n";
        exit(-1);
      }
    }
  
    if (theDamping->getClassTag() != dmpTag) {
      delete theDamping;
      theDamping = theBroker.getNewDamping(dmpTag);
      if (theDamping == 0) {
        opserr << "LadrunoDispBeamColumn2d::recvSelf -- could not get a Damping\n";
        exit(-1);
      }
    }
  
    theDamping->setDbTag((int)data(15));
    if (theDamping->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "LadrunoDispBeamColumn2d::recvSelf -- could not receive Damping\n";
      exit(-1);
    }
  }
  else {
    if (theDamping) {
      delete theDamping;
      theDamping = 0;
    }
  }

  // Ladruno (ADR 32) Tier-2: reconstruct the embedded hinge materials (classTag+dbTag in data Vector).
  if (hingeOn && hasHingeI) {
    if (theHinge == 0 || theHinge->getClassTag() != hClassTag) {
      if (theHinge) delete theHinge;
      theHinge = theBroker.getNewUniaxialMaterial(hClassTag);
      if (theHinge == 0) {
        opserr << "LadrunoDispBeamColumn2d::recvSelf() - Broker could not create hinge material of class "
               << hClassTag << endln;
        return -1;
      }
    }
    theHinge->setDbTag(hDbTag);
    if (theHinge->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "LadrunoDispBeamColumn2d::recvSelf() - hinge material failed to recv itself\n";
      return -1;
    }
  }
  else {
    if (theHinge) { delete theHinge; theHinge = 0; }
  }

  if (hingeOn && hasHingeJ) {
    if (theHingeJ == 0 || theHingeJ->getClassTag() != hClassTagJ) {
      if (theHingeJ) delete theHingeJ;
      theHingeJ = theBroker.getNewUniaxialMaterial(hClassTagJ);
      if (theHingeJ == 0) {
        opserr << "LadrunoDispBeamColumn2d::recvSelf() - Broker could not create hingeJ material of class "
               << hClassTagJ << endln;
        return -1;
      }
    }
    theHingeJ->setDbTag(hDbTagJ);
    if (theHingeJ->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "LadrunoDispBeamColumn2d::recvSelf() - hingeJ material failed to recv itself\n";
      return -1;
    }
  }
  else {
    if (theHingeJ) { delete theHingeJ; theHingeJ = 0; }
  }

  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
    this->update();
    if (!theDamping) this->getResistingForce();
  }

  return 0;
}

void
LadrunoDispBeamColumn2d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nLadrunoDispBeamColumn2d, element id:  " << this->getTag() << endln;
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tCoordTransf: " << crdTransf->getTag() << endln;
    s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;
    
    double L = crdTransf->getInitialLength();
    double P  = q(0);
    double M1 = q(1);
    double M2 = q(2);
    double V = (M1+M2)/L;
    
    s << "\tEnd 1 Forces (P V M): " << -P+p0[0]
      << " " << V+p0[1] << " " << M1 << endln;
    s << "\tEnd 2 Forces (P V M): " << P
      << " " << -V+p0[2] << " " << M2 << endln;
    
    beamInt->Print(s, flag);
    
    for (int i = 0; i < numSections; i++)
      theSections[i]->Print(s,flag);
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": " << this->getTag() << ", ";
	s << "\"type\": \"LadrunoDispBeamColumn2d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
    s << "\"sections\": [" ;
    for (int i = 0; i < numSections-1; i++)
      s << "\"" << theSections[i]->getTag() << "\", ";
    s << "\"" << theSections[numSections-1]->getTag() << "\"], ";
    s << "\"integration\": ";
    beamInt->Print(s, flag);
	s << ", \"massperlength\": " << rho << ", ";
    s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
  }
}


int
LadrunoDispBeamColumn2d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes)
{
    static thread_local Vector v1(3);
    static thread_local Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response*
LadrunoDispBeamColumn2d::setResponse(const char **argv, int argc,
			      OPS_Stream &output)
{
  if (spanLoadDirty && theNodes[0] != 0 && theNodes[1] != 0) {
    Element *savedActive = ops_TheActiveElement;
    ops_TheActiveElement = this;
    this->update();
    ops_TheActiveElement = savedActive;
  }

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","LadrunoDispBeamColumn2d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

  // global force - 
  if (LadrunoResp::isGlobalForce(argv[0])) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 1, P);
  
  
  // local force -
  } else if (LadrunoResp::isLocalForce(argv[0])) {

    output.tag("ResponseType","N1");
    output.tag("ResponseType","V1");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","N2");
    output.tag("ResponseType","V2");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 2, P);
  

  // basic force -
  } else if (LadrunoResp::isBasicForce(argv[0])) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 9, Vector(3));

  // global damping force - 
  } else if (theDamping && (strcmp(argv[0],"globalDampingForce") == 0 || strcmp(argv[0],"globalDampingForces") == 0)) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 21, P);
  
  
  // local damping force -
  } else if (theDamping && (strcmp(argv[0],"localDampingForce") == 0 || strcmp(argv[0],"localDampingForces") == 0)) {

    output.tag("ResponseType","N1");
    output.tag("ResponseType","V1");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","N2");
    output.tag("ResponseType","V2");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 22, P);
  

  // basic damping force -
  } else if (theDamping && (strcmp(argv[0],"basicDampingForce") == 0 || strcmp(argv[0],"basicDampingForces") == 0)) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 23, Vector(3));

  // basic stiffness -
  } else if (strcmp(argv[0],"basicStiffness") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 19, Matrix(3,3));

  // chord rotation -
  } else if (LadrunoResp::isBasicDef(argv[0])) {

    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta1");
    output.tag("ResponseType","theta2");

    theResponse =  new ElementResponse(this, 3, Vector(3));
  
  // plastic rotation -
  } else if (LadrunoResp::isPlasticDef(argv[0])) {

    output.tag("ResponseType","epsP");
    output.tag("ResponseType","theta1P");
    output.tag("ResponseType","theta2P");

    theResponse =  new ElementResponse(this, 4, Vector(3));

  } else if (LadrunoResp::isRayleighForce(argv[0])) {

    theResponse =  new ElementResponse(this, 12, P);
  }

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
      output.endTag();
    }
  }
  else if (strcmp(argv[0],"section") == 0) {

    if (argc > 1) {
      
      int sectionNum = atoi(argv[1]);
      
      if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

	output.tag("GaussPointOutput");
	output.attr("number",sectionNum);
	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamInt->getSectionLocations(numSections, L, xi);
	output.attr("eta",xi[sectionNum-1]*L);

	if (strcmp(argv[2],"dsdh") != 0) {
	  theResponse = theSections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
	} else {
	  int order = theSections[sectionNum-1]->getOrder();
	  theResponse = new ElementResponse(this, 76, Vector(order));
	  Information &info = theResponse->getInformation();
	  info.theInt = sectionNum;
	}
		
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
  
  // curvature sensitivity along element length
  else if (strcmp(argv[0],"dcurvdh") == 0)
    theResponse = new ElementResponse(this, 5, Vector(numSections));
  
  // basic deformation sensitivity
  else if (strcmp(argv[0],"dvdh") == 0)
    theResponse = new ElementResponse(this, 6, Vector(3));
  
  else if (strcmp(argv[0],"integrationPoints") == 0)
    theResponse = new ElementResponse(this, 7, Vector(numSections));
  
  else if (strcmp(argv[0],"integrationWeights") == 0)
    theResponse = new ElementResponse(this, 8, Vector(numSections));

  else if (strcmp(argv[0],"sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(numSections));
  
  else if (LadrunoResp::isEnergy(argv[0])) {
    theResponse = new ElementResponse(this, 10, 0.0);
  }
  // Ladruno: expose the element local frame (from the CrdTransf) as 9 packed
  // direction cosines so the Ladruno recorder can record MODEL/LOCAL_AXES instead of
  // falling back to a silent identity quaternion (apeGmsh beam-orientation gap).
  else if (LadrunoResp::isLocalAxes(argv[0])) {
    theResponse = new ElementResponse(this, 30, Vector(9));
  }
  // Ladruno (ADR 32) Tier-2: forward "hinge <resp>" / "hingeZ <resp>" (end-agnostic fallback)
  // and explicit "hingeI" / "hingeZI" / "hinge1" or "hingeJ" / "hingeZJ" / "hinge2" (no cross-end fallback).
  else if (hingeOn && argc > 1 &&
           (strcmp(argv[0],"hinge") == 0 || strcmp(argv[0],"hingeZ") == 0)) {
    UniaxialMaterial *hm = (theHinge != 0) ? theHinge : theHingeJ;
    if (hm != 0)
      theResponse = hm->setResponse(&argv[1], argc-1, output);
  }
  else if (hingeOn && theHinge != 0 && argc > 1 &&
           (strcmp(argv[0],"hingeI") == 0 || strcmp(argv[0],"hingeZI") == 0 || strcmp(argv[0],"hinge1") == 0)) {
    theResponse = theHinge->setResponse(&argv[1], argc-1, output);
  }
  else if (hingeOn && theHingeJ != 0 && argc > 1 &&
           (strcmp(argv[0],"hingeJ") == 0 || strcmp(argv[0],"hingeZJ") == 0 || strcmp(argv[0],"hinge2") == 0)) {
    theResponse = theHingeJ->setResponse(&argv[1], argc-1, output);
  }
  else if (hingeOn && (theHinge != 0 || theHingeJ != 0) && LadrunoResp::is(argv[0], "hingeSummary")) {
    output.tag("ResponseType", "theta_i");
    output.tag("ResponseType", "M_i");
    output.tag("ResponseType", "theta_j");
    output.tag("ResponseType", "M_j");
    theResponse = new ElementResponse(this, 31, Vector(4));
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
    output.tag("ResponseType", "hingeJumpI");
    theResponse = new ElementResponse(this, 34, 0.0);
  }
  else if (strcmp(argv[0], "hingeJumpJ") == 0) {
    output.tag("ResponseType", "hingeJumpJ");
    theResponse = new ElementResponse(this, 35, 0.0);
  }

  if (theResponse == 0)
    theResponse = crdTransf->setResponse(argv, argc, output);
  
  output.endTag();

  if (theResponse == 0)
    return Element::setResponse(argv, argc, output);
  else
    return theResponse;
}

int 
LadrunoDispBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
  if (spanLoadDirty && theNodes[0] != 0 && theNodes[1] != 0) {
    Element *savedActive = ops_TheActiveElement;
    ops_TheActiveElement = this;
    this->update();
    ops_TheActiveElement = savedActive;
  }

  double V;
  double L = crdTransf->getInitialLength();

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 12) {
    // own buffer, not the shared static P -- same reason as in
    // getResistingForceIncInertia (WP-116).  // Ladruno
    static thread_local Vector damp(6);
    damp.Zero();
    damp.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    return eleInfo.setVector(damp);

  } else if (responseID == 2) {
      P(3) =  q(0);
      P(0) = -q(0)+p0[0];
      P(2) = q(1);
      P(5) = q(2);
      V = (q(1)+q(2))/L;
      P(1) =  V+p0[1];
      P(4) = -V+p0[2];
      return eleInfo.setVector(P);
  }

  // Ladruno: local axes (vx,vy,vz dir cosines) from the CrdTransf
  else if (responseID == 30) {
    static thread_local Vector la(9);
    static thread_local Vector vx(3), vy(3), vz(3);
    crdTransf->getLocalAxes(vx, vy, vz);
    for (int i = 0; i < 3; i++) { la(i) = vx(i); la(i + 3) = vy(i); la(i + 6) = vz(i); }
    return eleInfo.setVector(la);
  }

  else if (responseID == 31) {
    static thread_local Vector hs(4);
    hs.Zero();
    if (hingeOn && theHinge != 0) {
      double fPMI = (dynamic_cast<LadrunoCohesiveHinge*>(theHinge) == 0 && pmP0 > 0.0)
                      ? ladrunoDisp2dScalePM(q(0), pmP0, pmAlpha) : 1.0;
      hs(0) = hingeJump;
      hs(1) = fPMI * theHinge->getStress();
    }
    if (hingeOn && theHingeJ != 0) {
      double fPMJ = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeJ) == 0 && pmP0 > 0.0)
                      ? ladrunoDisp2dScalePM(q(0), pmP0, pmAlpha) : 1.0;
      hs(2) = hingeJumpJ;
      hs(3) = fPMJ * theHingeJ->getStress();
    }
    return eleInfo.setVector(hs);
  }

  else if (responseID == 32 || responseID == 33) {
    double e1 = 0.0, e2 = 0.0;
    int cnt = 0;
    if (ldbc2dGetHingeImplexError(theHinge, e1))
      cnt++;
    if (ldbc2dGetHingeImplexError(theHingeJ, e2))
      cnt++;
    double val = (responseID == 32) ? ((e1 > e2) ? e1 : e2) : ((cnt > 0) ? (e1 + e2) / cnt : 0.0);
    return eleInfo.setDouble(val);
  }

  else if (responseID == 34) {
    return eleInfo.setDouble(hingeOn ? hingeJump : 0.0);
  }

  else if (responseID == 35) {
    return eleInfo.setDouble(hingeOn ? hingeJumpJ : 0.0);
  }

  else if (responseID == 9) {
    return eleInfo.setVector(q);
  }

  else if (responseID == 19) {
    static thread_local Matrix kb(3,3);
    this->getBasicStiff(kb);
    return eleInfo.setMatrix(kb);
  }

  else if (responseID == 21)
    return eleInfo.setVector(this->getDampingForce());

  else if (responseID == 22) {
    Vector Sd(3);
    Sd = theDamping->getDampingForce();
    P(3) =  Sd(0);
    P(0) = -Sd(0);
    P(2) = Sd(1);
    P(5) = Sd(2);
    V = (Sd(1)+Sd(2))/L;
    P(1) =  V;
    P(4) = -V;
    return eleInfo.setVector(P);
  }

  else if (responseID == 23)
    return eleInfo.setVector(theDamping->getDampingForce());

  // Chord rotation
  else if (responseID == 3) {
    return eleInfo.setVector(crdTransf->getBasicTrialDisp());
  }

  // Plastic rotation
  else if (responseID == 4) {
    static thread_local Vector vp(3);
    static thread_local Vector ve(3);
    static thread_local Vector qSec(3);
    static thread_local Matrix kb(3,3);
    this->getBasicStiff(kb, 1);
    qSec(0) = q(0) - q0[0];
    qSec(1) = q(1) - q0[1];
    qSec(2) = q(2) - q0[2];
    double alpha_eI = 0.0, alpha_eJ = 0.0;
    if (hingeOn) {
      if (theHinge != 0 && dynamic_cast<LadrunoCohesiveHinge*>(theHinge) == 0) {
        double k0 = theHinge->getInitialTangent();
        if (k0 > 0.0)
          alpha_eI = theHinge->getStress() / k0;
      }
      if (theHingeJ != 0 && dynamic_cast<LadrunoCohesiveHinge*>(theHingeJ) == 0) {
        double k0J = theHingeJ->getInitialTangent();
        if (k0J > 0.0)
          alpha_eJ = theHingeJ->getStress() / k0J;
      }
    }
    if (hingeOn && !dualHingeMode && theHingeJ == 0 && theHinge != 0 && L > 0.0) {
      double oneOverL = 1.0/L;
      double xi[maxNumSections];
      double wt[maxNumSections];
      beamInt->getSectionLocations(numSections, L, xi);
      beamInt->getSectionWeights(numSections, L, wt);
      double Kva0[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < numSections; i++) {
        int order = theSections[i]->getOrder();
        const ID &code = theSections[i]->getType();
        const Matrix &ks0 = theSections[i]->getInitialTangent();
        double xi6 = 6.0*xi[i];
        int idxMZ = -1, idxP = -1;
        for (int j = 0; j < order; j++) {
          if (code(j) == SECTION_RESPONSE_MZ) idxMZ = j;
          if (code(j) == SECTION_RESPONSE_P)  idxP  = j;
        }
        if (idxMZ >= 0) {
          double EI0 = ks0(idxMZ, idxMZ);
          Kva0[1] += -oneOverL*wt[i]*(xi6-4.0)*EI0;
          Kva0[2] += -oneOverL*wt[i]*(xi6-2.0)*EI0;
          if (idxP >= 0)
            Kva0[0] += -oneOverL*wt[i]*ks0(idxP, idxMZ);
        }
      }
      qSec(0) -= Kva0[0]*hingeJump;
      qSec(1) -= Kva0[1]*hingeJump;
      qSec(2) -= Kva0[2]*hingeJump;
      kb.Solve(qSec, ve);
      vp = crdTransf->getBasicTrialDisp();
      vp -= ve;
      double alpha_p = hingeJump - alpha_eI;
      vp(1) += (1.0 - hingeXiI)*alpha_p;
      vp(2) += hingeXiI*alpha_p;
    } else {
      // Note: in dualHingeMode, bulk sections are evaluated at v_tilde = v - W*alpha,
      // so qSec = kb0 * (v - W*alpha - vp_bulk) and v - kb0^-1*qSec = vp_bulk + W*alpha
      // includes the Hermite-projected hinge jumps W*alpha; subtract the elastic hinge
      // component W*alpha_e so only plastic hinge rotation W*alpha_p remains.
      kb.Solve(qSec, ve);
      vp = crdTransf->getBasicTrialDisp();
      vp -= ve;
      if (hingeOn) {
        double w1I = 1.0 - hingeXiI;
        double w2I = hingeXiI;
        double w1J = 1.0 - hingeXiJ;
        double w2J = hingeXiJ;
        vp(1) -= w1I * alpha_eI + w1J * alpha_eJ;
        vp(2) -= w2I * alpha_eI + w2J * alpha_eJ;
      }
    }
    return eleInfo.setVector(vp);
  }

  // Curvature sensitivity
  else if (responseID == 5) {
    Vector curv(numSections);
    return eleInfo.setVector(curv);
  }

  // Basic deformation sensitivity
  else if (responseID == 6) {  
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(1);
    return eleInfo.setVector(dvdh);
  }

  else if (responseID == 7) {
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = xi[i]*L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 8) {
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wt[i]*L;
    return eleInfo.setVector(weights);
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = theSections[i]->getTag();
    return eleInfo.setID(tags);
  }
  
  //by SAJalali
  else if (responseID == 10) {
	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamInt->getSectionWeights(numSections, L, xi);
	  double energy = 0;
	  double P_axial = q0[0];
	  for (int i = 0; i < numSections; i++) {
		  energy += theSections[i]->getEnergy()*xi[i] * L;
		  if (hingeOn && pmP0 > 0.0) {
		    int order = theSections[i]->getOrder();
		    const ID &code = theSections[i]->getType();
		    const Vector &s = theSections[i]->getStressResultant();
		    for (int j = 0; j < order; j++) {
		      if (code(j) == SECTION_RESPONSE_P)
		        P_axial += xi[i] * s(j);
		    }
		  }
	  }
	  if (hingeOn && theHinge != 0) {
		  double fPMI = (dynamic_cast<LadrunoCohesiveHinge*>(theHinge) == 0 && pmP0 > 0.0)
		                  ? ladrunoDisp2dScalePM(P_axial, pmP0, pmAlpha) : 1.0;
		  energy += fPMI * theHinge->getEnergy();
	  }
	  if (hingeOn && theHingeJ != 0) {
		  double fPMJ = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeJ) == 0 && pmP0 > 0.0)
		                  ? ladrunoDisp2dScalePM(P_axial, pmP0, pmAlpha) : 1.0;
		  energy += fPMJ * theHingeJ->getEnergy();
	  }
	  return eleInfo.setDouble(energy);
  }

  else
    return Element::getResponse(responseID, eleInfo);
}

int 
LadrunoDispBeamColumn2d::getResponseSensitivity(int responseID, int gradNumber,
					 Information &eleInfo)
{
  // Basic deformation sensitivity
  if (responseID == 3) {  
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    return eleInfo.setVector(dvdh);
  }

  // Basic force sensitivity
  else if (responseID == 9) {
    static thread_local Vector dqdh(3);

    dqdh.Zero();

    return eleInfo.setVector(dqdh);
  }

  // dsdh
  else if (responseID == 76) {

    int sectionNum = eleInfo.theInt;
    int order = theSections[sectionNum-1]->getOrder();
    const ID &code = theSections[sectionNum-1]->getType();

    Vector dsdh(order);
    dsdh = theSections[sectionNum-1]->getStressResultantSensitivity(gradNumber, true);

    const Vector &v = crdTransf->getBasicTrialDisp();
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0/L;

    const Matrix &ks = theSections[sectionNum-1]->getSectionTangent();

    Vector dedh(order);

    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);

    double x = xi[sectionNum-1];
    double xi6 = 6.0*x;

    int j;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	dedh(j) = oneOverL*dvdh(0); break;
      case SECTION_RESPONSE_MZ:
	dedh(j) = oneOverL*((xi6-4.0)*dvdh(1) + (xi6-2.0)*dvdh(2)); break;
      default:
	dedh(j) = 0.0; break;
      }
    }

    dsdh.addMatrixVector(1.0, ks, dedh, 1.0);

    return eleInfo.setVector(dsdh);
  }

  else
    return -1;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
LadrunoDispBeamColumn2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // If the parameter belongs to the element itself
  if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0], "lch") == 0 || strcmp(argv[0], "Lch") == 0 || strcmp(argv[0], "characteristicLength") == 0) {
    param.setValue(this->getCharacteristicLength());
    return param.addObject(2, this);
  }

  if ((strcmp(argv[0], "hinge") == 0 || strcmp(argv[0], "hingeZ") == 0 ||
       strcmp(argv[0], "hingeIJ") == 0 || strcmp(argv[0], "hingeZIJ") == 0) && argc >= 2) {
    int res = -1;
    if (theHinge) {
      int r = theHinge->setParameter(&argv[1], argc - 1, param);
      if (r != -1) res = r;
    }
    if (theHingeJ) {
      int r = theHingeJ->setParameter(&argv[1], argc - 1, param);
      if (r != -1) res = r;
    }
    if (res != -1) param.addObject(3, this);
    return res;
  }
  if ((strcmp(argv[0], "hingeI") == 0 || strcmp(argv[0], "hingeZI") == 0 || strcmp(argv[0], "hinge1") == 0) && argc >= 2) {
    if (!theHinge) return -1;
    int res = theHinge->setParameter(&argv[1], argc - 1, param);
    if (res != -1) param.addObject(3, this);
    return res;
  }
  if ((strcmp(argv[0], "hingeJ") == 0 || strcmp(argv[0], "hingeZJ") == 0 || strcmp(argv[0], "hinge2") == 0) && argc >= 2) {
    if (!theHingeJ) return -1;
    int res = theHingeJ->setParameter(&argv[1], argc - 1, param);
    if (res != -1) param.addObject(3, this);
    return res;
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
	int res = theSections[sectionNum]->setParameter(&argv[2], argc-2, param);
	if (res != -1) param.addObject(3, this);
	return res;
  }
  if (strcmp(argv[0], "allSections") == 0) {
    if (argc < 2)
      return -1;
    int res = -1;
    for (int i = 0; i < numSections; i++) {
      int ok = theSections[i]->setParameter(&argv[1], argc - 1, param);
      if (ok != -1) res = ok;
    }
    if (res != -1) param.addObject(3, this);
    return res;
  }
  // If the parameter belongs to a section or lower
  if (strstr(argv[0],"section") != 0) {
    
    if (argc < 3)
      return -1;
    
    // Get section number: 1...Np
    int sectionNum = atoi(argv[1]);
    
    if (sectionNum > 0 && sectionNum <= numSections) {
      int res = theSections[sectionNum-1]->setParameter(&argv[2], argc-2, param);
      if (res != -1) param.addObject(3, this);
      return res;
    }
    
    else
      return -1;
  }
  
  if (strstr(argv[0],"integration") != 0) {
    
    if (argc < 2)
      return -1;

    int res = beamInt->setParameter(&argv[1], argc-1, param);
    if (res != -1) param.addObject(3, this);
    return res;
  }
  
  int result =-1;
  // Default, send to every object
  int ok = 0;
  for (int i = 0; i < numSections; i++) {
    ok = theSections[i]->setParameter(argv, argc, param);
    if (ok != -1) 
      result = ok;
  }
  ok = beamInt->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;
  if (theHinge) {
    ok = theHinge->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  if (theHingeJ) {
    ok = theHingeJ->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  if (result != -1)
    param.addObject(3, this);
  return result;
}

int
LadrunoDispBeamColumn2d::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  if (parameterID == 2 || parameterID == 3) {
    if (parameterID == 2) {
      if (info.theDouble <= 0.0)
        return -1;
      lchMode = 2;
      userLch = info.theDouble;
      current_section_lch = userLch;
    }
    if (theNodes[0] != 0 && theNodes[1] != 0) {
      ops_TheActiveElement = this;
      int rc = this->update();
      ops_TheActiveElement = 0;
      return rc;
    }
    return 0;
  }
  return -1;  
}




int
LadrunoDispBeamColumn2d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}



const Matrix &
LadrunoDispBeamColumn2d::getInitialStiffSensitivity(int gradNumber)
{
  static thread_local Matrix kb(3,3);

  // Zero for integral
  kb.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Matrix ka(workArea, order, 3);
    ka.Zero();

    double xi6 = 6.0*xi[i];

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangentSensitivity(gradNumber);
        
    // Perform numerical integration
    double wti = wt[i]*oneOverL;
    double tmp;
    int j, k;
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
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 3; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 3; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }
  }

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);
  
  return K;
}

const Matrix &
LadrunoDispBeamColumn2d::getMassSensitivity(int gradNumber)
{
  K.Zero();

  if (rho == 0.0 || parameterID != 1)
    return K;
  
  double L = crdTransf->getInitialLength();
  if (cMass == 0)  {
    // lumped mass matrix
    double m = 0.5*L;
    K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;
  } else  {
    // consistent mass matrix
    static thread_local Matrix ml(6,6);
    ml.Zero();
    double m = L/420.0;
    ml(0,0) = ml(3,3) = m*140.0;
    ml(0,3) = ml(3,0) = m*70.0;
    
    ml(1,1) = ml(4,4) = m*156.0;
    ml(1,4) = ml(4,1) = m*54.0;
    ml(2,2) = ml(5,5) = m*4.0*L*L;
    ml(2,5) = ml(5,2) = -m*3.0*L*L;
    ml(1,2) = ml(2,1) = m*22.0*L;
    ml(4,5) = ml(5,4) = -ml(1,2);
    ml(1,5) = ml(5,1) = -m*13.0*L;
    ml(2,4) = ml(4,2) = -ml(1,5);
    
    // transform local mass matrix to global system
    K = crdTransf->getGlobalMatrixFromLocal(ml);
  }
  
  return K;
}



const Vector &
LadrunoDispBeamColumn2d::getResistingForceSensitivity(int gradNumber)
{
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  double dLdh = crdTransf->getdLdh();
  double d1oLdh = crdTransf->getd1overLdh();

  double dptsdh[maxNumSections];
  beamInt->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double dwtsdh[maxNumSections];
  beamInt->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

  // Zero for integration
  static thread_local Vector dqdh(3);
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
      default:
	break;
      }
    }

    const Vector &s = theSections[i]->getStressResultant();

    double dxi6dh = 6.0*dptsdh[i];
    double dwtLdh = wt[i]*dLdh + dwtsdh[i]*L;

    // Perform numerical integration on internal force gradient
    for (int j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	break;
      case SECTION_RESPONSE_MZ:
	break;
      default:
	break;
      }
    }

  }
  
  // Transform forces
  static thread_local Vector dp0dh(3);		// No distributed loads

  P.Zero();

  //////////////////////////////////////////////////////////////

  if (crdTransf->isShapeSensitivity()) {
    
    // Perform numerical integration to obtain basic stiffness matrix
    // Some extra declarations
    static thread_local Matrix kbmine(3,3);
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
      
      Matrix ka(workArea, order, 3);
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
	default:
	  break;
	}
      }
      for (j = 0; j < order; j++) {
	switch (code(j)) {
	case SECTION_RESPONSE_P:
	  for (k = 0; k < 3; k++) {
	    kbmine(0,k) += ka(j,k);
	  }
	  break;
	case SECTION_RESPONSE_MZ:
	  for (k = 0; k < 3; k++) {
	    tmp = ka(j,k);
	    kbmine(1,k) += (xi6-4.0)*tmp;
	    kbmine(2,k) += (xi6-2.0)*tmp;
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
LadrunoDispBeamColumn2d::commitSensitivity(int gradNumber, int numGrads)
{
  // Get basic deformation and sensitivities
  const Vector &v = crdTransf->getBasicTrialDisp();
  
  static thread_local Vector dvdh(3);
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
LadrunoDispBeamColumn2d::ladrunoThreadSafeUpdate(void) const
{
  if (theDamping != 0 || crdTransf == 0 || beamInt == 0 || theSections == 0)
    return false;
  for (int i = 0; i < numSections; i++) {
    if (theSections[i] == 0 || !theSections[i]->ladrunoThreadSafeUpdate())
      return false;
  }
  if (hingeOn) {
    if (theHinge == 0 && theHingeJ == 0) return false;
    if (theHinge != 0 && !theHinge->ladrunoThreadSafeUpdate()) return false;
    if (theHingeJ != 0 && !theHingeJ->ladrunoThreadSafeUpdate()) return false;
  }
  return true;
}

bool
LadrunoDispBeamColumn2d::ladrunoThreadSafeForm(void) const
{
  if (theDamping != 0 || crdTransf == 0 || beamInt == 0 || theSections == 0)
    return false;
  for (int i = 0; i < numSections; i++) {
    if (theSections[i] == 0 || !theSections[i]->ladrunoThreadSafeForm())
      return false;
  }
  if (hingeOn) {
    if (theHinge == 0 && theHingeJ == 0) return false;
    if (theHinge != 0 && !theHinge->ladrunoThreadSafeForm()) return false;
    if (theHingeJ != 0 && !theHingeJ->ladrunoThreadSafeForm()) return false;
  }
  return true;
}



