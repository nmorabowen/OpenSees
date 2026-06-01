/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// Authors: Nicolas Mora Bowen, Patricio Palacios, Jose Abell, Guppi (Ladruño)
// Created: 05/2026
//
// Factory function for the BezierTet10 element.
// Parses user input from both Tcl and Python interfaces.
//
// Usage:
//   element BezierTet10 $tag $nd1 $nd2 $nd3 $nd4 $nd5 $nd6 $nd7 $nd8 $nd9 $nd10
//                       $matTag
//                       <-bbar> <-cMass> <-rho $r>
//                       <-bodyForce $b1 $b2 $b3> <-pressure $p>
//
// Required arguments:
//   $tag        - unique element tag
//   $nd1-$nd10  - node tags (4 corners + 6 mid-edge), TenNodeTetrahedron order:
//                 v1 v2 v3 v4, e12 e23 e13 e14 e34 e24
//   $matTag     - 3D NDMaterial tag
//
// Optional arguments:
//   -bbar               - use B-bar formulation (near-incompressibility)
//   -cMass              - consistent mass (default is lumped ρVe/10, Eq. 57)
//   -rho $r             - mass density (else taken from the material)
//   -bodyForce $b1..$b3 - body force per unit volume (rampable via SelfWeight)
//   -pressure $p        - volume "pressure" hack acting in +z (as BezierTri6)
//
// Example (Python):
//   ops.nDMaterial('ElasticIsotropic', 1, 1000.0, 0.3)
//   ops.element('BezierTet10', 1, 1,2,3,4,5,6,7,8,9,10, 1, '-bbar')

#include "BezierTet10.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>

#include <string.h>

// Built into OpenSees statically (registered in OpenSeesElementCommands.cpp via
// `void* OPS_BezierTet10();`), so use plain C++ linkage — NOT the extern "C"
// OPS_Export DLL style, which would not resolve against that declaration.
void *OPS_BezierTet10()
{
    // Required: tag nd1..nd10 matTag = 12 remaining args (eleType consumed)
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();
    if (numRemainingArgs < 12) {
        opserr << "WARNING insufficient arguments\n"
               << "Want: element BezierTet10 tag nd1 nd2 nd3 nd4 nd5 nd6 nd7 "
                  "nd8 nd9 nd10 matTag <-bbar> <-cMass> <-rho r> "
                  "<-bodyForce b1 b2 b3> <-pressure p>\n";
        return 0;
    }

    // ─── Parse integer data: tag + 10 node tags ───────────────
    int iData[11];
    int numData = 11;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data for BezierTet10\n";
        return 0;
    }

    int tag = iData[0];

    // ─── Parse material tag ───────────────────────────────────
    int matTag;
    numData = 1;
    if (OPS_GetIntInput(&numData, &matTag) != 0) {
        opserr << "WARNING invalid matTag for BezierTet10 " << tag << "\n";
        return 0;
    }

    NDMaterial *theMat = OPS_getNDMaterial(matTag);
    if (theMat == 0) {
        opserr << "WARNING NDMaterial " << matTag
               << " not found for BezierTet10 " << tag << "\n";
        return 0;
    }

    // ─── Parse optional arguments ─────────────────────────────
    double rho = 0.0;
    double pressure = 0.0;
    double b1 = 0.0, b2 = 0.0, b3 = 0.0;
    bool useBbar = false;
    bool cMass = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *option = OPS_GetString();

        if (strcmp(option, "-bbar") == 0 || strcmp(option, "-Bbar") == 0) {
            useBbar = true;
        }
        else if (strcmp(option, "-cMass") == 0 || strcmp(option, "-cmass") == 0) {
            cMass = true;
        }
        else if (strcmp(option, "-rho") == 0 || strcmp(option, "-mass") == 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &rho) != 0) {
                opserr << "WARNING invalid rho for BezierTet10 " << tag << "\n";
                return 0;
            }
        }
        else if (strcmp(option, "-pressure") == 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &pressure) != 0) {
                opserr << "WARNING invalid pressure for BezierTet10 " << tag << "\n";
                return 0;
            }
        }
        else if (strcmp(option, "-bodyForce") == 0 || strcmp(option, "-body") == 0) {
            double bData[3];
            numData = 3;
            if (OPS_GetDoubleInput(&numData, bData) != 0) {
                opserr << "WARNING invalid body force for BezierTet10 " << tag << "\n";
                return 0;
            }
            b1 = bData[0];
            b2 = bData[1];
            b3 = bData[2];
        }
        else {
            opserr << "WARNING unknown option '" << option
                   << "' for BezierTet10 " << tag << "\n";
        }
    }

    // ─── Create element ───────────────────────────────────────
    Element *theEle = new BezierTet10(tag,
                                      iData[1], iData[2], iData[3], iData[4],
                                      iData[5], iData[6], iData[7], iData[8],
                                      iData[9], iData[10],
                                      *theMat, rho, b1, b2, b3,
                                      useBbar, cMass, pressure);

    if (theEle == 0) {
        opserr << "WARNING ran out of memory creating BezierTet10 " << tag << "\n";
        return 0;
    }

    return theEle;
}
