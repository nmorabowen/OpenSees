/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// Authors: Nicolas Mora Bowen, Patricio Palacios, Jose Abell, Guppi (Ladruño)
// Created: 03/2026
//
// Factory function for the BezierTri6 element.
// Parses user input from both Tcl and Python interfaces.
//
// Usage:
//   element BezierTri6 $tag $nd1 $nd2 $nd3 $nd4 $nd5 $nd6
//                      $thick $type $matTag
//                      <-bbar>
//                      <-pressure $p>
//                      <-bodyForce $b1 $b2>
//
// Required arguments:
//   $tag      - unique element tag
//   $nd1-$nd6 - node tags (3 corners + 3 mid-edge)
//               Ordering: corner1, corner2, corner3, mid12, mid23, mid31
//   $thick    - element thickness
//   $type     - material type: "PlaneStress" or "PlaneStrain"
//   $matTag   - NDMaterial tag
//
// Optional arguments:
//   -bbar               - use B-bar formulation (PlaneStrain/3D only)
//   -cMass              - consistent mass (default is lumped ρVe/6, Eq. 56)
//   -pressure $p        - surface pressure
//   -rho $r             - mass density (else taken from the material)
//   -bodyForce $b1 $b2  - body force per unit volume (rampable via SelfWeight)
//
// Example (Tcl):
//   nDMaterial ElasticIsotropic 1 1000.0 0.3
//   element BezierTri6 1  1 2 3 4 5 6  0.1 PlaneStrain 1 -bbar
//
// Example (Python):
//   ops.nDMaterial('ElasticIsotropic', 1, 1000.0, 0.3)
//   ops.element('BezierTri6', 1, 1, 2, 3, 4, 5, 6,
//               0.1, 'PlaneStrain', 1, '-bbar')

#include "BezierTri6.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>

#include <string.h>

// Built into OpenSees statically (registered in OpenSeesElementCommands.cpp via
// `void* OPS_BezierTri6();`), so use plain C++ linkage — NOT the extern "C"
// OPS_Export DLL style, which would not resolve against that declaration.
void *OPS_BezierTri6()
{
    // ─── Check minimum argument count ─────────────────────────
    // Required: eleType tag nd1 nd2 nd3 nd4 nd5 nd6 thick type matTag
    //          = 10 remaining args (eleType is already consumed)
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();
    if (numRemainingArgs < 10) {
        opserr << "WARNING insufficient arguments\n"
               << "Want: element BezierTri6 tag nd1 nd2 nd3 nd4 nd5 nd6 "
               << "thick type matTag <-bbar> <-cMass> <-pressure p> "
                  "<-rho r> <-bodyForce b1 b2>\n";
        return 0;
    }

    // ─── Parse integer data: tag + 6 node tags ────────────────
    int iData[7];
    int numData = 7;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data for BezierTri6\n";
        return 0;
    }

    int tag = iData[0];
    int nd1 = iData[1], nd2 = iData[2], nd3 = iData[3];
    int nd4 = iData[4], nd5 = iData[5], nd6 = iData[6];

    // ─── Parse thickness ──────────────────────────────────────
    double thick;
    numData = 1;
    if (OPS_GetDoubleInput(&numData, &thick) != 0) {
        opserr << "WARNING invalid thickness for BezierTri6 " << tag << "\n";
        return 0;
    }

    // ─── Parse type string ────────────────────────────────────
    const char *type = OPS_GetString();
    if (type == 0) {
        opserr << "WARNING missing type for BezierTri6 " << tag << "\n";
        return 0;
    }

    // Validate type
    if (strcmp(type, "PlaneStress") != 0 &&
        strcmp(type, "PlaneStrain") != 0) {
        opserr << "WARNING invalid type '" << type
               << "' for BezierTri6 " << tag
               << " (use PlaneStress or PlaneStrain)\n";
        return 0;
    }

    // ─── Parse material tag ───────────────────────────────────
    int matTag;
    numData = 1;
    if (OPS_GetIntInput(&numData, &matTag) != 0) {
        opserr << "WARNING invalid matTag for BezierTri6 " << tag << "\n";
        return 0;
    }

    // Look up the material
    NDMaterial *theMat = OPS_getNDMaterial(matTag);
    if (theMat == 0) {
        opserr << "WARNING NDMaterial " << matTag
               << " not found for BezierTri6 " << tag << "\n";
        return 0;
    }

    // ─── Parse optional arguments ─────────────────────────────
    double pressure = 0.0;
    double rho = 0.0;
    double b1 = 0.0, b2 = 0.0;
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
        else if (strcmp(option, "-pressure") == 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &pressure) != 0) {
                opserr << "WARNING invalid pressure for BezierTri6 " << tag << "\n";
                return 0;
            }
        }
        else if (strcmp(option, "-rho") == 0 || strcmp(option, "-mass") == 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &rho) != 0) {
                opserr << "WARNING invalid rho for BezierTri6 " << tag << "\n";
                return 0;
            }
        }
        else if (strcmp(option, "-bodyForce") == 0 || strcmp(option, "-body") == 0) {
            double bData[2];
            numData = 2;
            if (OPS_GetDoubleInput(&numData, bData) != 0) {
                opserr << "WARNING invalid body force for BezierTri6 " << tag << "\n";
                return 0;
            }
            b1 = bData[0];
            b2 = bData[1];
        }
        else {
            opserr << "WARNING unknown option '" << option
                   << "' for BezierTri6 " << tag << "\n";
        }
    }

    // ─── Guard: B-bar is for volumetric locking (3D / plane strain /
    //     axisymmetric). In plane stress ε_zz is free, there is no
    //     volumetric locking, and the (1/3) deviatoric split (Kadapa
    //     Eq. 45) does not apply. Disable B-bar with a warning rather
    //     than silently producing a wrong formulation.
    if (useBbar && strcmp(type, "PlaneStress") == 0) {
        opserr << "WARNING BezierTri6 " << tag
               << " - '-bbar' is not valid for PlaneStress (no volumetric "
                  "locking); ignoring -bbar\n";
        useBbar = false;
    }

    // ─── Create element ───────────────────────────────────────
    Element *theEle = new BezierTri6(tag,
                                      nd1, nd2, nd3, nd4, nd5, nd6,
                                      *theMat, type, thick,
                                      pressure, rho, b1, b2,
                                      useBbar, cMass);

    if (theEle == 0) {
        opserr << "WARNING ran out of memory creating BezierTri6 " << tag << "\n";
        return 0;
    }

    return theEle;
}
