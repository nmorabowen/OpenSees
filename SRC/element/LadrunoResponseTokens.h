/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno — one canonical vocabulary for element recorder tokens.
//
// Upstream every element hand-rolls its own strcmp chain, so the accepted
// spelling of the SAME quantity drifts between elements: vanilla Brick takes
// "stresses" but not "stress", FourNodeQuad takes both. The fork inherited that
// drift into every new element, which makes `recorder Element ... stress`
// silently produce an empty file on one element and data on its neighbour.
//
// This header holds the single alias table. An element asks
//
//     if (LadrunoResp::is(argv[0], "stress")) { ... }
//
// and every registered spelling of "stress" resolves to that branch.
//
// TWO RULES, both load-bearing:
//
//  1. Matching is SYMMETRIC — is(a, b) is true whenever a and b canonicalize
//     to the same row. So is(argv[0], "stresses") behaves identically to
//     is(argv[0], "stress"); an element cannot pick the "wrong" spelling.
//
//  2. An alias NEVER changes what the element emits. The ResponseType labels,
//     the response ID and the vector width stay keyed to the branch, not to
//     what the user typed. Otherwise a recorder's column headers would depend
//     on the spelling in the script, and MPCO/STKO readers would break.
//
// A token that is not in the table degrades to an exact strcmp, so registering
// a row is optional and forgetting one is never worse than upstream behaviour.
//
// Matching is case-sensitive, as everywhere else in OpenSees; case variants
// that already exist in the wild (rayleighForces) are registered explicitly.

#ifndef LadrunoResponseTokens_h
#define LadrunoResponseTokens_h

#include <string.h>

namespace LadrunoResp {

struct Row {
  const char *canonical;
  const char *alias[6];   // NULL-terminated
};

// Canonical spelling of s, or s itself when unregistered.
//
// The table is a flat list on purpose: it is walked once per setResponse call
// (i.e. once per recorder, at setup), never in the analysis loop.
//
// Rows are grouped by element family only for reading. Canonicalization is
// GLOBAL — two families must not claim the same spelling for different
// quantities. Check the table before adding a token to a new element.
inline const char *canon(const char *s)
{
  static const Row table[] = {

    // ---- forces (element-wide) -------------------------------------------
    {"force",              {"forces", "couplingForce", 0}},
    {"globalForce",        {"globalForces", 0}},
    {"localForce",         {"localForces", "frameForce", "interfaceForce", 0}},
    {"basicForce",         {"basicForces", 0}},
    {"globalDampingForce", {"globalDampingForces", 0}},
    {"localDampingForce",  {"localDampingForces", 0}},
    {"basicDampingForce",  {"basicDampingForces", 0}},
    {"RayleighForces",     {"rayleighForces", 0}},

    // ---- continuum / plane / shell ---------------------------------------
    {"stress",             {"stresses", 0}},
    {"strain",             {"strains", 0}},
    {"stressPlaneStrain",  {"stressesPlaneStrain", 0}},
    {"stressTotal",        {"stressesTotal", 0}},
    {"stiff",              {"stiffness", "tangent", "tangentStiff",
                            "tangentStiffness", 0}},
    {"stiffInitial",       {"initialStiffness", "initialStiff", "stiffInit", 0}},
    {"basicStiffness",     {"basicStiff", 0}},
    {"charLength",         {"characteristicLength", "lch", 0}},
    {"material",           {"integrPoint", "integrationPoint", 0}},
    {"hourglassEnergy",    {"hgEnergy", "hourglassWork", "hourglassDissipation",
                            "hgDissipation", 0}},
    {"Jbar",               {"jbar", 0}},
    {"alpha",              {"eas", 0}},
    {"flux",               {"darcyFlux", 0}},

    // ---- beams ------------------------------------------------------------
    // chordRotation / basicDeformation / deformation are the SAME quantity
    // (the basic trial displacement v); upstream DispBeamColumn already groups
    // the first two, LadrunoIMKBeam used the third.
    {"chordRotation",      {"chordDeformation", "basicDeformation",
                            "deformation", "deformations", 0}},
    {"plasticRotation",    {"plasticDeformation", "hingeRotation",
                            "hingeRotations", 0}},
    {"hingeMoment",        {"hingeMoments", 0}},

    // ---- penalty / tie / coupling family ----------------------------------
    {"penalty",            {"kt", "k", "penaltyStiffness", 0}},
    {"rotPenalty",         {"kr", 0}},
    {"massPenalty",        {"mpenalty", 0}},
    {"lambda",             {"augLambda", 0}},
    {"lambdaP",            {"augLambdaP", "plambda", 0}},
    {"lambdaR",            {"augLambdaR", "rlambda", 0}},
    {"dtCritical",         {"dtcr", 0}},
    {"pressureGap",        {"pgap", 0}},
    {"pressureForce",      {"pforce", 0}},
    {"rotationGap",        {"rgap", 0}},
    {"rotationForce",      {"rforce", "moment", 0}},
    {"frameGap",           {"localGap", 0}},
    {"initGap",            {"offset", 0}},
    {"barAxis",            {"dir", 0}},
    {"bondEnergy",         {"bondDissipation", 0}},
    {"tiedDOFs",           {"nGap", 0}},
    {"rotationAxes",       {"nKept", 0}},

    // ---- rigid body -------------------------------------------------------
    {"comDisp",            {"centroidDisp", 0}},
    {"orientation",        {"quaternion", 0}},
    {"angularVel",         {"omega", 0}},
    {"angularMom",         {"L", 0}},
    {"omegaBody",          {"angularVelBody", 0}},

    {0, {0}}
  };

  if (s == 0)
    return 0;
  for (const Row *r = table; r->canonical != 0; ++r) {
    if (strcmp(s, r->canonical) == 0)
      return r->canonical;
    for (int i = 0; r->alias[i] != 0; ++i)
      if (strcmp(s, r->alias[i]) == 0)
        return r->canonical;
  }
  return s;
}

// True when a and b name the same quantity. Unregistered tokens fall back to
// an exact strcmp, i.e. exactly the upstream behaviour.
inline bool is(const char *a, const char *b)
{
  if (a == 0 || b == 0)
    return false;
  return strcmp(canon(a), canon(b)) == 0;
}

} // namespace LadrunoResp

#endif
