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

// Ladruno (relativeMassDamping): see RelativeMassDamping.h.

#include <RelativeMassDamping.h>
#include <TimeSeries.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <MeshRegion.h>
#include <Matrix.h>
#include <ID.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

#include <map>
#include <vector>
#include <cstring>

RelativeMassDamping::RelativeMassDamping(int t, int n, const int *d, TimeSeries **s)
  : tag(t), numRef(n), dofs(0), series(0), factors(0),
    lastTime(0.0), haveCache(false), result(0)
{
  if (numRef > 0) {
    dofs    = new int[numRef];
    series  = new TimeSeries *[numRef];
    factors = new double[numRef];
    for (int i = 0; i < numRef; i++) {
      dofs[i]    = d[i];
      series[i]  = (s[i] != 0) ? s[i]->getCopy() : 0;
      factors[i] = 0.0;
    }
  }
}

RelativeMassDamping::~RelativeMassDamping()
{
  for (int i = 0; i < numRef; i++)
    if (series[i] != 0)
      delete series[i];
  if (dofs != 0)    delete [] dofs;
  if (series != 0)  delete [] series;
  if (factors != 0) delete [] factors;
}

const Vector &
RelativeMassDamping::getVelocity(double time, int numDOF)
{
  if (!haveCache || time != lastTime) {
    for (int i = 0; i < numRef; i++)
      factors[i] = (series[i] != 0) ? series[i]->getFactor(time) : 0.0;
    lastTime = time;
    haveCache = true;
  }
  if (result.Size() != numDOF)
    result.resize(numDOF);
  result.Zero();
  for (int i = 0; i < numRef; i++)
    if (dofs[i] >= 0 && dofs[i] < numDOF)
      result(dofs[i]) += factors[i];
  return result;
}

// ------------------------------------------------------------------ name

// The command name, chosen by the author 2026-10-10; a rename is this one line.
const char *OPS_RelativeMassDampingCmdName = "relativeMassDamping";

// ------------------------------------------------------------------ registry

static std::map<int, RelativeMassDamping *> theRelativeMassDampings;

int
OPS_addRelativeMassDamping(RelativeMassDamping *theRef)
{
  if (theRef == 0)
    return -1;
  if (theRelativeMassDampings.count(theRef->getTag()) != 0)
    return -1;
  theRelativeMassDampings[theRef->getTag()] = theRef;
  return 0;
}

RelativeMassDamping *
OPS_getRelativeMassDamping(int tag)
{
  std::map<int, RelativeMassDamping *>::iterator it = theRelativeMassDampings.find(tag);
  return (it == theRelativeMassDampings.end()) ? 0 : it->second;
}

int
OPS_numRelativeMassDampings(void)
{
  return (int)theRelativeMassDampings.size();
}

void
OPS_clearAllRelativeMassDampings(void)
{
  for (std::map<int, RelativeMassDamping *>::iterator it = theRelativeMassDampings.begin();
       it != theRelativeMassDampings.end(); ++it)
    delete it->second;
  theRelativeMassDampings.clear();
}

// ------------------------------------------------------------------ command

int
OPS_RelativeMassDamping(void)
{
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0) {
    opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " - no domain\n";
    return -1;
  }
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " tag -vel dof tsTag ... (-region r ... | -node n ... | -allMassNodes)\n";
    return -1;
  }

  int numData = 1;
  int tag;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " - invalid tag\n";
    return -1;
  }
  if (OPS_getRelativeMassDamping(tag) != 0) {
    opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " - tag " << tag << " already exists\n";
    return -1;
  }

  std::vector<int> dofs;
  std::vector<TimeSeries *> series;
  std::vector<int> regions, nodes;
  bool allMassNodes = false;
  enum {NONE, REG, NOD} mode = NONE;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-vel") == 0) {
      mode = NONE;
      int dv[2];
      numData = 2;
      if (OPS_GetNumRemainingInputArgs() < 2 || OPS_GetIntInput(&numData, dv) < 0) {
        opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " -vel dof tsTag - invalid input\n";
        return -1;
      }
      TimeSeries *ts = OPS_getTimeSeries(dv[1]);
      if (ts == 0) {
        opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " - time series " << dv[1] << " not found\n";
        return -1;
      }
      if (dv[0] < 1) {
        opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " - dof must be >= 1\n";
        return -1;
      }
      dofs.push_back(dv[0] - 1);
      series.push_back(ts);
    } else if (strcmp(opt, "-region") == 0) {
      mode = REG;
    } else if (strcmp(opt, "-node") == 0) {
      mode = NOD;
    } else if (strcmp(opt, "-allMassNodes") == 0) {
      mode = NONE;
      allMassNodes = true;
    } else if (mode != NONE) {
      OPS_ResetCurrentInputArg(-1);
      int v;
      numData = 1;
      if (OPS_GetIntInput(&numData, &v) < 0) {
        opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " - invalid token " << opt << "\n";
        return -1;
      }
      if (mode == REG) regions.push_back(v); else nodes.push_back(v);
    } else {
      opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " - unknown option " << opt << "\n";
      return -1;
    }
  }

  if (dofs.empty()) {
    opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " - no -vel given\n";
    return -1;
  }
  int nsel = (allMassNodes ? 1 : 0) + (regions.empty() ? 0 : 1) + (nodes.empty() ? 0 : 1);
  if (nsel != 1) {
    opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " - give exactly one of -region, -node, -allMassNodes\n";
    return -1;
  }

  // collect candidate node tags
  std::vector<int> cand;
  if (allMassNodes) {
    NodeIter &it = theDomain->getNodes();
    Node *nd;
    while ((nd = it()) != 0)
      cand.push_back(nd->getTag());
  } else if (!regions.empty()) {
    for (size_t r = 0; r < regions.size(); r++) {
      MeshRegion *reg = theDomain->getRegion(regions[r]);
      if (reg == 0) {
        opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " - region " << regions[r] << " not found\n";
        return -1;
      }
      const ID &rn = reg->getNodes();
      for (int i = 0; i < rn.Size(); i++)
        cand.push_back(rn(i));
    }
  } else {
    cand = nodes;
  }

  RelativeMassDamping *theRef = new RelativeMassDamping(tag, (int)dofs.size(), &dofs[0], &series[0]);
  if (OPS_addRelativeMassDamping(theRef) < 0) {
    delete theRef;
    opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " - could not add\n";
    return -1;
  }

  // attach to nodes with nodal mass (each node once even if listed twice)
  std::map<int, bool> seen;
  int nAttached = 0, nAlphaZero = 0;
  std::vector<double> sumAlphaM(dofs.size(), 0.0);
  for (size_t k = 0; k < cand.size(); k++) {
    if (seen.count(cand[k])) continue;
    seen[cand[k]] = true;
    Node *nd = theDomain->getNode(cand[k]);
    if (nd == 0) continue;               // not on this process
    const Matrix &M = nd->getMass();
    int n = nd->getNumberDOF();
    bool hasMass = false;
    for (int i = 0; i < M.noRows() && !hasMass; i++)
      for (int j = 0; j < M.noCols(); j++)
        if (M(i, j) != 0.0) {hasMass = true; break;}
    if (!hasMass) continue;
    if (nd->getMassDampingReference() != 0 && nd->getMassDampingReference() != theRef) {
      opserr << "WARNING " << OPS_RelativeMassDampingCmdName << " " << tag << " - node " << cand[k]
             << " already carries reference " << nd->getMassDampingReference()->getTag() << "; replaced\n";
    }
    nd->setMassDampingReference(theRef);
    nAttached++;
    double aM = nd->getRayleighDampingFactor();
    if (aM == 0.0) nAlphaZero++;
    for (size_t d = 0; d < dofs.size(); d++)
      if (dofs[d] < n)
        sumAlphaM[d] += aM * M(dofs[d], dofs[d]);
  }

  std::vector<double> out;
  out.push_back((double)nAttached);
  out.push_back((double)nAlphaZero);
  for (size_t d = 0; d < dofs.size(); d++)
    out.push_back(sumAlphaM[d]);
  numData = (int)out.size();
  OPS_SetDoubleOutput(&numData, &out[0], false);
  return 0;
}
