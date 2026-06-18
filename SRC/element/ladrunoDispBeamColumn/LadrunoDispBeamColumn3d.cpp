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
#include <NDMaterial.h>         // Ladruno (ADR 34): coupled biaxial cohesive hinge (order-2)
#include <CrdTransf.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
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

// Ladruno (ADR 32): active element for crack-band materials (see 2D sibling).
extern Element *ops_TheActiveElement;

Matrix LadrunoDispBeamColumn3d::K(12,12);
Vector LadrunoDispBeamColumn3d::P(12);
double LadrunoDispBeamColumn3d::workArea[200];

void* OPS_LadrunoDispBeamColumn3d()
{
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
    UniaxialMaterial *hingeMatZ = 0;  // Ladruno (ADR 33) Tier-2: embedded strong-axis (Mz) cohesive hinge
    UniaxialMaterial *hingeMatY = 0;  // Ladruno (ADR 33 PR-3b): embedded weak-axis (My) cohesive hinge (biaxial)
    NDMaterial *hingeMatC = 0;        // Ladruno (ADR 34): coupled biaxial (Mz-My) cohesive hinge (order-2 NDMaterial)
    numData = 1;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type,"-cMass") == 0) {
	    cmass = 1;
	} else if(strcmp(type,"-mass") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&mass) < 0) {
		    opserr<<"WARNING: invalid mass\n";
		    return 0;
		}
	    }
	}
    else if (strcmp(type, "-damp") == 0) {

        if (OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
            theDamping = OPS_getDamping(dampingTag);
            if (theDamping == 0) {
                opserr << "damping not found\n";
                return 0;
            }
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
    }

    // Ladruno (ADR 34): the coupled biaxial cohesive replaces the two block-diagonal scalar laws —
    // reject mixing it with -hinge/-hingeY.
    if (hingeMatC != 0 && (hingeMatZ != 0 || hingeMatY != 0)) {
        opserr << "WARNING LadrunoDispBeamColumn3d: -hingeBiaxial (coupled) is mutually exclusive "
                  "with -hinge/-hingeY (block-diagonal); use one or the other.\n";
        return 0;
    }

    // Ladruno (ADR 33) PR-3b: the weak-axis jump is an ADD-ON to the strong-axis hinge — the
    // coupled biaxial 2x2 condenses both bending jumps together, so -hingeY requires -hinge.
    if (hingeMatY != 0 && hingeMatZ == 0) {
        opserr << "WARNING LadrunoDispBeamColumn3d: -hingeY requires -hinge (the biaxial hinge "
                  "condenses the strong- and weak-axis jumps together).\n";
        return 0;
    }

    // Ladruno (ADR 33) Tier-2: the embedded hinge and the -nl bowing strain are mutually
    // exclusive in v1 — the ½(θz²+θy²) membrane term couples the rotation jump into the axial
    // channel. Reject the combination loudly rather than silently mis-condense.
    if ((hingeMatZ != 0 || hingeMatY != 0 || hingeMatC != 0) && nlGeom != 0) {
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
    SectionForceDeformation** sections = new SectionForceDeformation *[secTags.Size()];
    for(int i=0; i<secTags.Size(); i++) {
	sections[i] = OPS_getSectionForceDeformation(secTags(i));
	if(sections[i] == 0) {
	    opserr<<"section "<<secTags(i)<<"not found\n";
		delete [] sections;
	    return 0;
	}
    }
    
    Element *theEle =  new LadrunoDispBeamColumn3d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					    *bi,*theTransf,mass,cmass, theDamping,lchMode,userLch,nlGeom,hingeMatZ,hingeMatY,hingeMatC);
    delete [] sections;
    return theEle;
}

LadrunoDispBeamColumn3d::LadrunoDispBeamColumn3d(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration &bi,
				   CrdTransf &coordTransf, double r, int cm,
				   Damping *damping, int lchM, double userL, int nlG,
				   UniaxialMaterial *hingeMatZ, UniaxialMaterial *hingeMatY,
				   NDMaterial *hingeMatC)
:Element (tag, ELE_TAG_LadrunoDispBeamColumn3d),
numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2),
Q(12), q(6), rho(r), cMass(cm), parameterID(0), theDamping(0),
current_section_lch(0.0), lchMode(lchM), userLch(userL), nlGeom(nlG),
hingeOn(0), theHingeZ(0), hingeJumpZ(0.0), hingeJumpCommitZ(0.0), hingeKaaZ(0.0),
hingeMscaleZ(0.0),
theHingeY(0), hingeJumpY(0.0), hingeJumpCommitY(0.0), hingeMscaleY(0.0),
hingeKaaInvZZ(0.0), hingeKaaInvZY(0.0), hingeKaaInvYY(0.0),
theHingeC(0)
{
  hingeKvZ[0] = hingeKvZ[1] = hingeKvZ[2] = hingeKvZ[3] = hingeKvZ[4] = hingeKvZ[5] = 0.0;
  hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
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
hingeOn(0), theHingeZ(0), hingeJumpZ(0.0), hingeJumpCommitZ(0.0), hingeKaaZ(0.0),
hingeMscaleZ(0.0),
theHingeY(0), hingeJumpY(0.0), hingeJumpCommitY(0.0), hingeMscaleY(0.0),
hingeKaaInvZZ(0.0), hingeKaaInvZY(0.0), hingeKaaInvYY(0.0),
theHingeC(0)
{
  hingeKvZ[0] = hingeKvZ[1] = hingeKvZ[2] = hingeKvZ[3] = hingeKvZ[4] = hingeKvZ[5] = 0.0;
  hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
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
	//opserr << "FATAL ERROR LadrunoDispBeamColumn3d (tag: %d), node not found in domain",
	//	this->getTag());
	
	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    if (dofNd1 != 6 || dofNd2 != 6) {
	//opserr << "FATAL ERROR LadrunoDispBeamColumn3d (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());
	
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
      if (theHingeZ) retVal += theHingeZ->commitState();
      if (theHingeY) {                       // PR-3b biaxial
        hingeJumpCommitY = hingeJumpY;
        retVal += theHingeY->commitState();
      }
      if (theHingeC) {                       // ADR 34 coupled (jumps committed, plus the law's history)
        hingeJumpCommitY = hingeJumpY;
        retVal += theHingeC->commitState();
      }
    }

    return retVal;
}

int
LadrunoDispBeamColumn3d::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    if (theDamping) retVal += theDamping->revertToLastCommit();

    // Ladruno (ADR 33) Tier-2: restore the jump (and the cohesive law) to the committed state
    // so a rejected step does not leave a stale/open hinge.
    if (hingeOn) {
      hingeJumpZ = hingeJumpCommitZ;
      if (theHingeZ) retVal += theHingeZ->revertToLastCommit();
      if (theHingeY) {                       // PR-3b biaxial
        hingeJumpY = hingeJumpCommitY;
        retVal += theHingeY->revertToLastCommit();
      }
      if (theHingeC) {                       // ADR 34 coupled
        hingeJumpY = hingeJumpCommitY;
        retVal += theHingeC->revertToLastCommit();
      }
    }

    return retVal;
}

int
LadrunoDispBeamColumn3d::revertToStart()
{
    int retVal = 0;

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
      hingeKvZ[0] = hingeKvZ[1] = hingeKvZ[2] = hingeKvZ[3] = hingeKvZ[4] = hingeKvZ[5] = 0.0;
      hingeMscaleZ = 0.0;
      if (theHingeZ) retVal += theHingeZ->revertToStart();
      if (theHingeY) {                       // PR-3b biaxial
        hingeJumpY = 0.0;
        hingeJumpCommitY = 0.0;
        hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
        hingeMscaleY = 0.0;
        hingeKaaInvZZ = hingeKaaInvZY = hingeKaaInvYY = 0.0;
        retVal += theHingeY->revertToStart();
      }
      if (theHingeC) {                       // ADR 34 coupled
        hingeJumpY = 0.0;
        hingeJumpCommitY = 0.0;
        hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
        hingeMscaleY = 0.0;
        hingeKaaInvZZ = hingeKaaInvZY = hingeKaaInvYY = 0.0;
        retVal += theHingeC->revertToStart();
      }
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

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  // Ladruno (ADR 33) Tier-2: when the embedded hinge is active, the bulk Mz (and, biaxial, My)
  // section strain includes the -alpha_*/L enhancement and the jump(s) are converged by an inner
  // Newton (which sets the section deformations itself). GATED so the no-hinge path below is
  // untouched. PR-3b: theHingeY != 0 routes to the coupled biaxial solve; Z-only stays PR-3a.
  if (hingeOn) {
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
  double oneOverL = 1.0/L;
  const double Gbar = -oneOverL;   // bounded incompatible-mode enhancement on the Mz channel

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  double alpha = hingeJumpCommitZ;   // warm start from the committed jump (path-correct)

  const int maxIter = 30;
  double Kva[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double Kaa = 0.0;

  for (int iter = 0; iter < maxIter; iter++) {

    int err = 0;
    double sumM = 0.0;     // sum_k wt_k * M_sec,k(MZ)
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
      opserr << "LadrunoDispBeamColumn3d::solveHingeJump() - failed setTrialSectionDeformations()\n";
      return err;
    }

    // --- cohesive law at the current jump: M_coh_z(alpha_z), dM_coh_z/dalpha_z ---
    theHingeZ->setTrialStrain(alpha);
    double Mcoh = theHingeZ->getStress();
    double Kcoh = theHingeZ->getTangent();

    double h = -sumM + Mcoh;                 // enhancement residual
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
      return 0;
    }

    alpha -= h/KaaSolve;                     // guarded Newton step on the scalar jump
  }

  hingeKaaZ = (fabs(Kaa) < 1.0e-300) ? 1.0e-300 : Kaa;
  for (int a = 0; a < 6; a++) hingeKvZ[a] = Kva[a];
  hingeJumpZ = alpha;
  opserr << "WARNING LadrunoDispBeamColumn3d (tag " << this->getTag()
         << "): embedded-hinge inner Newton did not converge in " << maxIter << " iters\n";
  return 0;
}

// Ladruno (ADR 33) PR-3b: eigenvalue-floored inverse of a symmetric 2x2 [[a,b],[b,d]].
// Each eigenvalue's MAGNITUDE is floored to `floor` (its SIGN is preserved — a softening-induced
// negative eigenvalue is physical), then the inverse is reconstructed spectrally
// A^{-1} = sum_i (1/lam_i) q_i q_i^T. This bounds EVERY mode, unlike a det-floored adjugate whose
// numerator blows up ~1e8x along the near-null eigenvector at staggered activation (ADR 33 risk).
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

// Ladruno (ADR 33) Tier-2 PR-3b: unified inner Newton on the BIAXIAL rotation jump
// alpha = [alpha_z, alpha_y] on the Mz/My rows of the 6-DOF basic system. ONE
// setTrialSectionDeformation per IP sets BOTH bulk curvatures (a second pass would overwrite the
// strain vector and corrupt the section state from the first). The bulk sections see
//   kappa_z_bulk = (1/L)((6xi-4)v1+(6xi-2)v2) - alpha_z/L,   kappa_y_bulk = (1/L)((6xi-4)v3+(6xi-2)v4) - alpha_y/L,
// the axial/torsion channels stay linear. Enhancement equilibrium and the TRUE coupled tangent:
//   h_z = -sum wt*s(MZ) + M_coh_z(alpha_z) = 0,   h_y = -sum wt*s(MY) + M_coh_y(alpha_y) = 0
//   K_aa = [[ (1/L)sum wt*ks(MZ,MZ) + dMz/dalphaz ,  (1/L)sum wt*½(ks(MZ,MY)+ks(MY,MZ)) ],
//           [  (sym)                              ,  (1/L)sum wt*ks(MY,MY) + dMy/dalphay ]]
// solved with an EIGENVALUE-FLOORED 2x2 inverse. The 6x2 K_v-alpha (hingeKvZ | hingeKvY) carries
// the cross-tangent rows so the condensed off-diagonals of the 6x6 basic kb are right. The basic
// FORCE needs no correction (at h=0 the sections hold kappa_bulk, so q = sum wt*B^T*s IS condensed).
int
LadrunoDispBeamColumn3d::solveHingeJumpBiaxial(const Vector &v, double L)
{
  double oneOverL = 1.0/L;
  const double Gbar = -oneOverL;   // bounded incompatible-mode enhancement (per bending axis)

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  double az = hingeJumpCommitZ;    // warm start from the committed jumps (path-correct)
  double ay = hingeJumpCommitY;

  const int maxIter = 50;
  double KvZ[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double KvY[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double iZZ = 0.0, iZY = 0.0, iYY = 0.0;

  // Ladruno (ADR 34): adaptive step damping — the coupled cohesive (theHingeC) has a
  // sign-discontinuous tangent at the elliptical onset (r=1) and a frozen-mix tangent that is
  // only exact on radial jump paths, so a full Newton step can overshoot the activation kink and
  // oscillate. Halving the step whenever the residual grows tames it; the damping is a no-op when
  // the residual decreases monotonically, so the block-diagonal (PR-3b) path is unaffected.
  double damp = 1.0;
  double hnormPrev = 1.0e300;

  for (int iter = 0; iter < maxIter; iter++) {

    int err = 0;
    double sumMz = 0.0, sumMy = 0.0;       // sum wt*M_sec(MZ/MY)
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
      opserr << "LadrunoDispBeamColumn3d::solveHingeJumpBiaxial() - failed setTrialSectionDeformations()\n";
      return err;
    }

    // --- cohesive moments + cohesive tangent at the current jumps ---
    // theHingeC (ADR 34): ONE coupled law -> moment vector + a FULL 2x2 cohesive tangent (the
    // off-diagonal Czy enters K_aa). Else (PR-3b): two block-diagonal scalar laws (Czy = 0).
    double Mz, My, Czz, Cyy, Czy;
    if (theHingeC != 0) {
      static Vector jumpVec(2);
      jumpVec(0) = az; jumpVec(1) = ay;
      theHingeC->setTrialStrain(jumpVec);
      const Vector &Mvec = theHingeC->getStress();
      const Matrix &Kc   = theHingeC->getTangent();
      Mz = Mvec(0); My = Mvec(1);
      Czz = Kc(0,0); Cyy = Kc(1,1);
      Czy = 0.5*(Kc(0,1) + Kc(1,0));   // symmetrize the cohesive coupling
    } else {
      theHingeZ->setTrialStrain(az);
      Mz = theHingeZ->getStress();  Czz = theHingeZ->getTangent();
      theHingeY->setTrialStrain(ay);
      My = theHingeY->getStress();  Cyy = theHingeY->getTangent();
      Czy = 0.0;
    }

    double hz = -sumMz + Mz;
    double hy = -sumMy + My;

    // TRUE coupled 2x2 K_aa (bulk + cohesive on the diagonal; off-diagonal = bulk ks(MZ,MY) +
    // the cohesive coupling Czy — nonzero only for the ADR-34 coupled interaction surface)
    double Kaa_zz = oneOverL*sumEIz + Czz;
    double Kaa_yy = oneOverL*sumEIy + Cyy;
    double Kaa_zy = oneOverL*sumCpl + Czy;

    // EIGENVALUE floor: bounds every mode against the positive bulk terms + the cohesive tangents
    double floor = 1.0e-8*(fabs(oneOverL*sumEIz) + fabs(oneOverL*sumEIy) + fabs(Czz) + fabs(Cyy)) + 1.0e-300;
    ladrunoFlooredInv2x2(Kaa_zz, Kaa_zy, Kaa_yy, floor, iZZ, iZY, iYY);

    // per-axis STABLE moment scales (do not collapse when a fully-broken LINEAR hinge carries M->0)
    double MhereZ = fabs(Mz) + fabs(sumMz);
    double MhereY = fabs(My) + fabs(sumMy);
    if (MhereZ > hingeMscaleZ) hingeMscaleZ = MhereZ;
    if (MhereY > hingeMscaleY) hingeMscaleY = MhereY;

    bool okZ = fabs(hz) <= 1.0e-11*hingeMscaleZ + 1.0e-12;
    bool okY = fabs(hy) <= 1.0e-11*hingeMscaleY + 1.0e-12;
    if (okZ && okY) {
      hingeJumpZ = az; hingeJumpY = ay;
      for (int a = 0; a < 6; a++) { hingeKvZ[a] = KvZ[a]; hingeKvY[a] = KvY[a]; }
      hingeKaaInvZZ = iZZ; hingeKaaInvZY = iZY; hingeKaaInvYY = iYY;
      return 0;
    }

    // adaptive damping: shrink the step when the residual grew (overshot the activation kink),
    // relax back toward the full step otherwise. No-op while the residual decreases monotonically.
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
  opserr << "WARNING LadrunoDispBeamColumn3d (tag " << this->getTag()
         << "): biaxial embedded-hinge inner Newton did not converge in " << maxIter << " iters\n";
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
  static Matrix kb(6,6);
  
  // Zero for integral
  kb.Zero();
  q.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Matrix ka(workArea, order, 6);
    ka.Zero();

    double xi6 = 6.0*xi[i];

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getSectionTangent();
    const Vector &s = theSections[i]->getStressResultant();
        
    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
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
    
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));
    double si;
    for (j = 0; j < order; j++) {
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

  // Ladruno (ADR 32) Stage-1: ½(θz²+θy²) nonlinear geometry. The loop above filled
  // the LINEAR kb and the linear part of q; overwrite kb with the NL geometric
  // tangent and add the biaxial bowing contribution to the basic force q.
  if (nlGeom) {
    this->getBasicStiff(kb);
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

  // Ladruno (ADR 33) Tier-2: statically condense the strong-axis jump alpha_z out of the 6x6
  // basic stiffness BEFORE crdTransf (the PINNED INVARIANT — CorotCrdTransf3d builds its entire
  // geometric stiffness from q/pb and reads kb only through Tp^T kb Tp, so it exposes no seam for
  // element-internal DOFs): a guarded rank-1 update kb -= hingeKvZ hingeKvZ^T / hingeKaaZ. The kb
  // here is the INLINE-built linear basic stiffness (the nlGeom path is dead under -hinge). The
  // basic FORCE q needs NO correction (sections already hold the converged kappa_z_bulk).
  if (hingeOn) {
    if (theHingeY != 0 || theHingeC != 0) {
      // PR-3b biaxial / ADR-34 coupled: rank-2 kb -= K_v-alpha * K_aa_floored^{-1} * K_v-alpha^T, with the 6x2
      // K_v-alpha = [hingeKvZ | hingeKvY] and the symmetric eigenvalue-floored 2x2 inverse.
      double iZZ = hingeKaaInvZZ, iZY = hingeKaaInvZY, iYY = hingeKaaInvYY;
      for (int a = 0; a < 6; a++)
        for (int b = 0; b < 6; b++)
          kb(a,b) -= hingeKvZ[a]*(iZZ*hingeKvZ[b] + iZY*hingeKvY[b])
                   + hingeKvY[a]*(iZY*hingeKvZ[b] + iYY*hingeKvY[b]);
    } else if (hingeKaaZ != 0.0) {
      // PR-3a strong-axis only: guarded rank-1 kb -= hingeKvZ hingeKvZ^T / hingeKaaZ.
      double invKaa = 1.0/hingeKaaZ;
      for (int a = 0; a < 6; a++)
        for (int b = 0; b < 6; b++)
          kb(a,b) -= hingeKvZ[a]*hingeKvZ[b]*invKaa;
    }
  }

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
  double oneOverL = 1.0/L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
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
      static Matrix Cz1(1,6);
      static Matrix Cy1(1,6);
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
      static Matrix ksB(1,6);
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
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
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
  // Ladruno (ADR 32): the damping stiffness-multiplier is applied by getInitialStiff
  // (mirroring the 2D sibling), NOT here. Stock 3D applied it in getBasicStiff, but that
  // made getTangentStiff's -nl path (which calls getBasicStiff) pick up the multiplier on
  // the tangent while the inlined linear path did not — a branch-dependent inconsistency
  // under -damp. Applying it only in getInitialStiff keeps both tangent branches consistent.
}

const Matrix&
LadrunoDispBeamColumn3d::getInitialStiff()
{
  static Matrix kb(6,6);

  this->getBasicStiff(kb, 1);
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
    static Matrix ml(12,12);
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

  return;
}

int 
LadrunoDispBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
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
    static Vector Raccel(12);
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
  P = this->getResistingForce();
  
  if (theDamping) P += this->getDampingForce();
  
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
  
    P(0) += m*accel1(0);
    P(1) += m*accel1(1);
    P(2) += m*accel1(2);
    P(6) += m*accel2(0);
    P(7) += m*accel2(1);
    P(8) += m*accel2(2);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(12);
    for (int i=0; i<6; i++)  {
      accel(i)   = accel1(i);
      accel(i+6) = accel2(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }
    
    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

  } else {

    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
  }
  
  return P;
}

int
LadrunoDispBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;
  
  static Vector data(26);  // Ladruno (ADR 32 +3; ADR 33 +2; PR-3b +2; ADR 34 +3 hingeOnC/classTag/dbTag)
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

  // Ladruno (ADR 33) Tier-2: embedded-hinge flag + committed strong-axis jump (the hinge
  // material itself is sent after the sections, mirroring the section dbTag/classTag pattern).
  data(19) = hingeOn;
  data(20) = hingeJumpCommitZ;
  // Ladruno (ADR 33 PR-3b): weak-axis (biaxial) hinge presence flag + committed jump
  data(21) = (theHingeY != 0) ? 1 : 0;
  data(22) = hingeJumpCommitY;
  // Ladruno (ADR 34): coupled biaxial cohesive (NDMaterial) presence + classTag/dbTag (folded into
  // data to avoid any FE_Datastore ID-size collision; the material itself sendSelf's below).
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
    opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  //
  // send the sections
  //
  
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

  // Ladruno (ADR 33) Tier-2: send the embedded hinge material metadata as ONE ID — size 2
  // (strong-axis only) or size 4 (biaxial). Two separate same-size IDs sharing the element
  // dbTag+commitTag would COLLIDE in the FE_Datastore (key = dbTag/commitTag/size), so the
  // second silently overwrites the first; the combined ID (size 2 vs 4) keys distinctly. Each
  // material then sends itself under its own (distinct) material dbTag. Z-only keeps the PR-3a stream.
  if (hingeOn && theHingeZ) {
    int nHingeID = (theHingeY != 0) ? 4 : 2;
    ID hingeID(nHingeID);
    hingeID(0) = theHingeZ->getClassTag();
    int hDbTag = theHingeZ->getDbTag();
    if (hDbTag == 0) {
      hDbTag = theChannel.getDbTag();
      if (hDbTag != 0) theHingeZ->setDbTag(hDbTag);
    }
    hingeID(1) = hDbTag;
    if (theHingeY != 0) {
      hingeID(2) = theHingeY->getClassTag();
      int hDbTagY = theHingeY->getDbTag();
      if (hDbTagY == 0) {
        hDbTagY = theChannel.getDbTag();
        if (hDbTagY != 0) theHingeY->setDbTag(hDbTagY);
      }
      hingeID(3) = hDbTagY;
    }
    if (theChannel.sendID(dbTag, commitTag, hingeID) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send hinge ID\n";
      return -1;
    }
    if (theHingeZ->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send hinge material\n";
      return -1;
    }
    if (theHingeY != 0 && theHingeY->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send weak-axis hinge material\n";
      return -1;
    }
  }

  // Ladruno (ADR 34): send the coupled biaxial cohesive material (mutually exclusive with Z/Y, so
  // no ID block above; its classTag/dbTag rode in the data Vector).
  if (theHingeC != 0 && theHingeC->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoDispBeamColumn3d::sendSelf() - failed to send coupled biaxial hinge material\n";
    return -1;
  }

  return 0;
}

int
LadrunoDispBeamColumn3d::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;
  
  static Vector data(26);  // Ladruno (ADR 32 +3; ADR 33 +2; PR-3b +2; ADR 34 +3 hingeOnC/classTag/dbTag)

  if (theChannel.recvVector(dbTag, commitTag, data) < 0)  {
    opserr << "LadrunoDispBeamColumn3d::recvSelf() - failed to recv data Vector\n";
    return -1;
  }
  
  this->setTag((int)data(0));
  connectedExternalNodes(0) = (int)data(1);
  connectedExternalNodes(1) = (int)data(2);
  int nSect = (int)data(3);
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

  // Ladruno (ADR 32): regularization-length config + nonlinear-geometry flag
  lchMode = (int)data(16);
  userLch = data(17);
  nlGeom = (int)data(18);
  current_section_lch = 0.0;

  // Ladruno (ADR 33) Tier-2: embedded-hinge flag + committed strong-axis jump (the material is
  // reconstructed after the sections; hingeKaaZ/hingeKvZ are transient, rebuilt in update()).
  hingeOn = (int)data(19);
  hingeJumpCommitZ = data(20);
  hingeJumpZ = hingeJumpCommitZ;
  hingeKaaZ = 0.0;
  hingeKvZ[0] = hingeKvZ[1] = hingeKvZ[2] = hingeKvZ[3] = hingeKvZ[4] = hingeKvZ[5] = 0.0;
  hingeMscaleZ = 0.0;

  // Ladruno (ADR 33 PR-3b): weak-axis (biaxial) hinge — presence flag + committed jump (the
  // material is reconstructed after the strong-axis one; the 2x2 inverse cache is transient).
  int hasHingeY = (int)data(21);
  hingeJumpCommitY = data(22);
  hingeJumpY = hingeJumpCommitY;
  hingeKvY[0] = hingeKvY[1] = hingeKvY[2] = hingeKvY[3] = hingeKvY[4] = hingeKvY[5] = 0.0;
  hingeMscaleY = 0.0;
  hingeKaaInvZZ = hingeKaaInvZY = hingeKaaInvYY = 0.0;

  // Ladruno (ADR 34): coupled biaxial cohesive (NDMaterial) presence + classTag/dbTag
  int hasHingeC = (int)data(23);
  int cClassTag = (int)data(24);
  int cDbTag    = (int)data(25);

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
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*nSect);
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "LadrunoDispBeamColumn3d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  //
  // now receive the sections
  //
  
  if (numSections != nSect) {

    //
    // we do not have correct number of sections, must delete the old and create
    // new ones before can recvSelf on the sections
    //

    // delete the old
    if (numSections != 0) {
      for (int i=0; i<numSections; i++)
	delete theSections[i];
      delete [] theSections;
    }

    // create a new array to hold pointers
    theSections = new SectionForceDeformation *[nSect];
    if (theSections == 0) {
      opserr << "LadrunoDispBeamColumn3d::recvSelf() - out of memory creating sections array of size" <<
	nSect << endln;
      exit(-1);
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

    // 
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //
    
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
	  opserr << "LadrunoDispBeamColumn3d::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvSelf on it
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
  
    // Check that the Damping is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (theDamping->getClassTag() != dmpTag) {
      delete theDamping;
      theDamping = theBroker.getNewDamping(dmpTag);
      if (theDamping == 0) {
        opserr << "LadrunoDispBeamColumn3d::recvSelf -- could not get a Damping\n";
        exit(-1);
      }
    }
  
    // Now, receive the Damping
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

  // Ladruno (ADR 33) Tier-2: reconstruct the block-diagonal hinge material(s) from the combined ID
  // (size 2 strong-axis only, size 4 biaxial). Skipped entirely for the coupled law (ADR 34), which
  // is mutually exclusive with -hinge/-hingeY and sent NO Z/Y ID block (only theHingeC, below).
  if (hingeOn && hasHingeC == 0) {
    int nHingeID = hasHingeY ? 4 : 2;
    ID hingeID(nHingeID);
    if (theChannel.recvID(dbTag, commitTag, hingeID) < 0) {
      opserr << "LadrunoDispBeamColumn3d::recvSelf() - failed to recv hinge ID\n";
      return -1;
    }
    int hClassTag = hingeID(0);
    int hDbTag = hingeID(1);
    if (theHingeZ == 0 || theHingeZ->getClassTag() != hClassTag) {
      if (theHingeZ) delete theHingeZ;
      theHingeZ = theBroker.getNewUniaxialMaterial(hClassTag);
      if (theHingeZ == 0) {
        opserr << "LadrunoDispBeamColumn3d::recvSelf() - Broker could not create hinge material of class "
               << hClassTag << endln;
        return -1;
      }
    }
    theHingeZ->setDbTag(hDbTag);
    if (theHingeZ->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "LadrunoDispBeamColumn3d::recvSelf() - hinge material failed to recv itself\n";
      return -1;
    }

    // Ladruno (ADR 33 PR-3b): the weak-axis material (sent after the strong-axis one)
    if (hasHingeY) {
      int hClassTagY = hingeID(2);
      int hDbTagY = hingeID(3);
      if (theHingeY == 0 || theHingeY->getClassTag() != hClassTagY) {
        if (theHingeY) delete theHingeY;
        theHingeY = theBroker.getNewUniaxialMaterial(hClassTagY);
        if (theHingeY == 0) {
          opserr << "LadrunoDispBeamColumn3d::recvSelf() - Broker could not create weak-axis hinge material of class "
                 << hClassTagY << endln;
          return -1;
        }
      }
      theHingeY->setDbTag(hDbTagY);
      if (theHingeY->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "LadrunoDispBeamColumn3d::recvSelf() - weak-axis hinge material failed to recv itself\n";
        return -1;
      }
    }
    else {
      if (theHingeY) { delete theHingeY; theHingeY = 0; }
    }
  }
  else {
    if (theHingeZ) { delete theHingeZ; theHingeZ = 0; }
    if (theHingeY) { delete theHingeY; theHingeY = 0; }
  }

  // Ladruno (ADR 34): reconstruct the coupled biaxial cohesive material (NDMaterial broker).
  if (hingeOn && hasHingeC) {
    if (theHingeC == 0 || theHingeC->getClassTag() != cClassTag) {
      if (theHingeC) delete theHingeC;
      theHingeC = theBroker.getNewNDMaterial(cClassTag);
      if (theHingeC == 0) {
        opserr << "LadrunoDispBeamColumn3d::recvSelf() - Broker could not create coupled hinge nDMaterial of class "
               << cClassTag << endln;
        return -1;
      }
    }
    theHingeC->setDbTag(cDbTag);
    if (theHingeC->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "LadrunoDispBeamColumn3d::recvSelf() - coupled hinge material failed to recv itself\n";
      return -1;
    }
  }
  else {
    if (theHingeC) { delete theHingeC; theHingeC = 0; }
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
    static Vector v1(3);
    static Vector v2(3);

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

    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

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
    }  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

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
    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {
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
    }  else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
	      || strcmp(argv[0],"basicDeformation") == 0) {

      output.tag("ResponseType","eps");
      output.tag("ResponseType","thetaZ_1");
      output.tag("ResponseType","thetaZ_2");
      output.tag("ResponseType","thetaY_1");
      output.tag("ResponseType","thetaY_2");
      output.tag("ResponseType","thetaX");

      theResponse = new ElementResponse(this, 3, Vector(6));

    // plastic rotation -
    } else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

    output.tag("ResponseType","epsP");
    output.tag("ResponseType","thetaZP_1");
    output.tag("ResponseType","thetaZP_2");
    output.tag("ResponseType","thetaYP_1");
    output.tag("ResponseType","thetaYP_2");
    output.tag("ResponseType","thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(6));
  

  } else if (strcmp(argv[0],"RayleighForces") == 0 || strcmp(argv[0],"rayleighForces") == 0) {

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
    else if (strcmp(argv[0], "energy") == 0) {
      theResponse = new ElementResponse(this, 13, 0.0);
    }
    // Ladruno: expose the element local frame (from the CrdTransf) as 9 packed
    // direction cosines so the Ladruno recorder can record MODEL/LOCAL_AXES instead of
    // falling back to a silent identity quaternion (apeGmsh beam-orientation gap).
    else if (strcmp(argv[0],"localAxes") == 0) {
      theResponse = new ElementResponse(this, 30, Vector(9));
    }
    // Ladruno (ADR 33) Tier-2: forward "hinge <resp>" to the embedded strong-axis cohesive
    // material (hinge stress = M_coh_z, hinge strain = the jump [[theta_z]], hinge energy =
    // int Mz d[[theta_z]], hinge kappaMax/damage/...). Lets the gates read the hinge state.
    else if (hingeOn && theHingeZ != 0 && argc > 1 && strcmp(argv[0],"hinge") == 0) {
      theResponse = theHingeZ->setResponse(&argv[1], argc-1, output);
    }
    // Ladruno (ADR 33 PR-3b): forward "hingeY <resp>" to the embedded weak-axis cohesive material.
    else if (hingeOn && theHingeY != 0 && argc > 1 && strcmp(argv[0],"hingeY") == 0) {
      theResponse = theHingeY->setResponse(&argv[1], argc-1, output);
    }
    // Ladruno (ADR 34): forward "hingeBiaxial <resp>" to the coupled biaxial cohesive material.
    else if (hingeOn && theHingeC != 0 && argc > 1 &&
             (strcmp(argv[0],"hingeBiaxial") == 0 || strcmp(argv[0],"hingeC") == 0)) {
      theResponse = theHingeC->setResponse(&argv[1], argc-1, output);
    }

    if (theResponse == 0)
      theResponse = crdTransf->setResponse(argv, argc, output);
    
  output.endTag();
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
    static Vector la(9);
    static Vector vx(3), vy(3), vz(3);
    crdTransf->getLocalAxes(vx, vy, vz);
    for (int i = 0; i < 3; i++) { la(i) = vx(i); la(i + 3) = vy(i); la(i + 6) = vz(i); }
    return eleInfo.setVector(la);
  }

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());
    
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
    static Matrix kb(6,6);
    this->getBasicStiff(kb);
    return eleInfo.setMatrix(kb);
  }  
  // Chord rotation
  else if (responseID == 3)
    return eleInfo.setVector(crdTransf->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(6);
    static Vector ve(6);
    static Matrix kb(6,6);
    this->getBasicStiff(kb,1);
    kb.Solve(q, ve);
    vp = crdTransf->getBasicTrialDisp();
    vp -= ve;
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
	  return eleInfo.setDouble(energy);
  }

  else
    return -1;
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
    static Matrix ml(12,12);
    //double m = rho*L/420.0;
    double m = L/420.0;
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



const Vector &
LadrunoDispBeamColumn3d::getResistingForceSensitivity(int gradNumber)
{
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Zero for integration
  static Vector dqdh(6);
  dqdh.Zero();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    
    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    //double wti = wts(i);
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
  static Vector dp0dh(6);		// No distributed loads

  P.Zero();

  //////////////////////////////////////////////////////////////

  if (crdTransf->isShapeSensitivity()) {
    
    // Perform numerical integration to obtain basic stiffness matrix
    // Some extra declarations
    static Matrix kbmine(6,6);
    kbmine.Zero();
    q.Zero();
    
    double tmp;
    
    int j, k;
    
    for (int i = 0; i < numSections; i++) {
      
      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      
      //double xi6 = 6.0*pts(i,0);
      double xi6 = 6.0*xi[i];
      //double wti = wts(i);
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
  
  static Vector dvdh(6);
  dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);

  // Some extra declarations
  double d1oLdh = crdTransf->getd1overLdh();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    
    Vector e(workArea, order);
    
    //double xi6 = 6.0*pts(i,0);
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

