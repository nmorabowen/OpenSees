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

// Ladruno WP-192 -- see LadrunoThreadedAssembly.h for the strategy.

#include <LadrunoThreadedAssembly.h>
#include <LadrunoThreads.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Integrator.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <OPS_Globals.h>
#include <profiler/Profiler.h>
#include <vector>
#include <string.h>

namespace {

// One parallel compute pass over a segment fes[0, n).
struct LadrunoFECtx {
  FE_Element **fes;
  Integrator *integ;
  int kind;
  Element *active;                        // the master's ops_TheActiveElement
  std::vector<std::vector<double> > *slots;
  std::vector<int> *dims;
};

int
ladrunoFEBody(int k, void *c)
{
  LadrunoFECtx *ctx = (LadrunoFECtx *)c;
  // The serial tangent/residual loops never touch ops_TheActiveElement, so a
  // material reading it there sees whatever Domain::update left behind; give
  // every worker that same value (it is thread_local, WP-107).
  ops_TheActiveElement = ctx->active;
  FE_Element *fe = ctx->fes[k];
  std::vector<double> &slot = (*ctx->slots)[k];
  if (ctx->kind == LADRUNO_FE_TANGENT) {
    const Matrix &K = fe->getTangent(ctx->integ);
    const int nr = K.noRows(), nc = K.noCols();
    slot.resize((size_t)nr * (size_t)nc);
    for (int j = 0; j < nc; j++)
      for (int i = 0; i < nr; i++)
        slot[(size_t)j * nr + i] = K(i, j);      // column-major, as Matrix
    (*ctx->dims)[2 * k] = nr;
    (*ctx->dims)[2 * k + 1] = nc;
  } else {
    const Vector &R = fe->getResidual(ctx->integ);
    const int n = R.Size();
    slot.resize((size_t)n);
    for (int i = 0; i < n; i++)
      slot[i] = R(i);
    (*ctx->dims)[2 * k] = n;
    (*ctx->dims)[2 * k + 1] = 1;
  }
  return 0;
}

// The announcement latch: once per (FE count, safe count, threads, kind).
struct LadrunoFEAnnounce {
  int nFE, nSafe, nThreads;
  LadrunoFEAnnounce() : nFE(-1), nSafe(-1), nThreads(-1) {}
};
LadrunoFEAnnounce ladrunoFEAnnounced[2];

} // namespace

bool
ladrunoThreadedFELoop(AnalysisModel *theModel, LinearSOE *theSOE,
                      Integrator *theIntegrator, int kind,
                      int failCode, const char *who, int &result)
{
  const char *why = 0;
  const int nThreads = ladruno_loopThreads(&why);
  if (nThreads <= 1)
    return false;                  // Domain::update already announced why
  if (!ladruno_loopEnabled(kind == LADRUNO_FE_TANGENT ? LADRUNO_LOOP_TANGENT : LADRUNO_LOOP_RESIDUAL))
    return false;                  // WP-192 diagnostic mask (LADRUNO_THREADS_LOOPS)
  if (theModel == 0 || theSOE == 0)
    return false;
  if (ops_profiler::theProfiler().deep())
    return false;                  // H6, announced by Domain::update

  // H3: snapshot the shared mutable FE cursor.
  std::vector<FE_Element *> fes;
  {
    FE_EleIter &it = theModel->getFEs();
    FE_Element *fe;
    while ((fe = it()) != 0)
      fes.push_back(fe);
  }
  const int nFE = (int)fes.size();
  if (nFE < 2)
    return false;

  // Per-FE opt-in. ladrunoThreadSafeFE() also allocates the element's Rayleigh
  // pool slot (serially), so it is evaluated here, before any thread starts.
  std::vector<char> safe(nFE, 0);
  int nSafe = 0;
  int firstUnsafe = -1;
  for (int i = 0; i < nFE; i++) {
    if (fes[i]->ladrunoThreadSafeFE()) {
      safe[i] = 1;
      nSafe++;
    } else if (firstUnsafe < 0)
      firstUnsafe = i;
  }
  if (nSafe < 2)
    return false;

  // H4: Node's trial-state getters allocate lazily on first touch.
  Domain *theDomain = theModel->getDomainPtr();
  if (theDomain != 0) {
    NodeIter &nIt = theDomain->getNodes();
    Node *nd;
    while ((nd = nIt()) != 0) {
      nd->getTrialDisp();
      nd->getTrialVel();
      nd->getTrialAccel();
    }
  }

  LadrunoFEAnnounce &a = ladrunoFEAnnounced[kind == LADRUNO_FE_TANGENT ? 0 : 1];
  if (a.nFE != nFE || a.nSafe != nSafe || a.nThreads != nThreads) {
    a.nFE = nFE; a.nSafe = nSafe; a.nThreads = nThreads;
    opserr << "ladrunoThreads: element " << (kind == LADRUNO_FE_TANGENT ? "tangent" : "residual")
           << " loop threaded on " << nThreads << " threads: " << nSafe << " of " << nFE
           << " FE_Elements computed in parallel, scattered in serial order";
    if (firstUnsafe >= 0) {
      Element *e = fes[firstUnsafe]->getElement();
      opserr << "; " << (nFE - nSafe) << " run serially in place (first: FE "
             << fes[firstUnsafe]->getTag();
      if (e != 0)
        opserr << ", element " << e->getTag() << " classTag " << e->getClassTag();
      else
        opserr << ", a constraint FE";
      opserr << ")";
    }
    opserr << "\n";
  }

  // Chunk: large enough to amortise the fork/join, small enough to bound the
  // buffer (C x idSize^2 doubles: 512 x 24^2 x 8 B = 2.4 MB for shells).
  const int C = (nThreads * 64 < 512) ? 512 : nThreads * 64;
  static std::vector<std::vector<double> > slots;   // master-owned; reused
  static std::vector<int> dims;
  if ((int)slots.size() < C)
    slots.resize(C);
  if ((int)dims.size() < 2 * C)
    dims.resize(2 * C);

  Element *active = ops_TheActiveElement;
  result = 0;
  int i = 0;
  while (i < nFE) {
    if (!safe[i]) {
      FE_Element *fe = fes[i];
      if (kind == LADRUNO_FE_TANGENT) {
        const Matrix &K = fe->getTangent(theIntegrator);
        if (theSOE->addA(K, fe->getID()) < 0) {
          opserr << who << " failed in addA for ID " << fe->getID();
          result = failCode;
        }
      } else {
        const Vector &R = fe->getResidual(theIntegrator);
        if (theSOE->addB(R, fe->getID()) < 0) {
          opserr << who << " failed in addB for ID " << fe->getID();
          result = failCode;
        }
      }
      i++;
      continue;
    }
    int j = i;
    while (j < nFE && safe[j] && j - i < C)
      j++;

    LadrunoFECtx ctx;
    ctx.fes = &fes[i];
    ctx.integ = theIntegrator;
    ctx.kind = kind;
    ctx.active = active;
    ctx.slots = &slots;
    ctx.dims = &dims;
    ladruno_parallelFor(j - i, nThreads, 4, &ladrunoFEBody, &ctx, 0);

    // The ordered serial scatter: the same calls, in the same order, as serial.
    for (int k = 0; k < j - i; k++) {
      FE_Element *fe = fes[i + k];
      double *data = slots[k].data();
      if (kind == LADRUNO_FE_TANGENT) {
        Matrix K(data, dims[2 * k], dims[2 * k + 1]);
        if (theSOE->addA(K, fe->getID()) < 0) {
          opserr << who << " failed in addA for ID " << fe->getID();
          result = failCode;
        }
      } else {
        Vector R(data, dims[2 * k]);
        if (theSOE->addB(R, fe->getID()) < 0) {
          opserr << who << " failed in addB for ID " << fe->getID();
          result = failCode;
        }
      }
    }
    i = j;
  }
  ops_TheActiveElement = active;
  return true;
}
