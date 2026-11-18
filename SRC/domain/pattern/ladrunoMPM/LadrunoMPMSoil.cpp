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

// Authors: Nicolas Mora Bowen, Guppi (Ladruno)
// Created: 11/2026
//
// LadrunoMPMSoil — soil particles at the MPM host (WP-232k, ADR 232 §9 phase 5; the
// frozen interface tests/mpm_soil_interface.md, D1-D23). This file holds:
//
//  - the one stage walk (D7) and its three routes (D9): `mpm stage`, the
//    `parameter ... loadPattern` route (LadrunoMPMPattern::setParameter registers the
//    host) and mode A's `updateMaterialStage` (a particle element hands the host the
//    MaterialStageParameter registration);
//  - the stage guard (D10): LadrunoSANISAND's stage flag is a class-wide STATIC
//    (ManzariDafalias.h:557); a flip by a route that did not reach the particles would
//    let them flip lazily without Elastic2Plastic (LadrunoSANISAND.cpp:3397-3406), so
//    `mpm analyze` is refused until `mpm stage` walks them;
//  - the soil queries (D16): sigma_zz, material responses per particle, `info materials`;
//  - the unsymmetric-tangent check for the symmetric-SOE refusal (D15);
//  - the first refusing particle's name (D13);
//  - mode F's prescribed grid velocities, `mpm gridvel` (D17, D18; #1109: before this
//    the explicit engine accepted a gridvel and ignored it).

#include "LadrunoMPMPattern.h"

#include <Domain.h>
#include <NDMaterial.h>
#include <LadrunoSANISAND.h>
#include <LadrunoHS.h>
#include <Parameter.h>
#include <Information.h>
#include <Response.h>
#include <DummyStream.h>
#include <Vector.h>
#include <classTags.h>
#include <OPS_Globals.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <memory>

// WP-232k (review item 7): the one definition of the soil class tags
bool
LadrunoMPMPattern::isSanisandTag(int ct)
{
  return ct == ND_TAG_LadrunoSANISAND3D || ct == ND_TAG_LadrunoSANISANDPlaneStrain;
}

bool
LadrunoMPMPattern::isHSTag(int ct)
{
  return ct == ND_TAG_LadrunoHS;
}

bool
LadrunoMPMPattern::isSoilTag(int ct)
{
  return isSanisandTag(ct) || isHSTag(ct);
}

int
LadrunoMPMPattern::stageRegId(int matTag, bool intPath)
{
  for (size_t k = 0; k < stageRegs.size(); k++)
    if (stageRegs[k].matTag == matTag && stageRegs[k].intPath == intPath)
      return STAGE_PARAM_BASE + (int)k;
  stageRegs.push_back(StageReg{matTag, intPath});
  return STAGE_PARAM_BASE + (int)stageRegs.size() - 1;
}

namespace {

bool parseInt(const char *s, int &v)
{
  if (s == 0 || *s == '\0') return false;
  char *end = 0;
  const long x = std::strtol(s, &end, 10);
  if (end == s || *end != '\0') return false;
  v = (int)x;
  return true;
}

}  // namespace

// D9: `updateMaterialStage <tag>` / `materialState <tag>` naming a SANISAND or HS
// particle material. Anything else keeps the 232c forwarding (DruckerPrager's
// `materialState` route is unchanged bit for bit).
bool
LadrunoMPMPattern::soilStageRequest(const char **argv, int argc, int &matTag) const
{
  if (argc < 2 || argv == 0 || argv[0] == 0) return false;
  if (std::strcmp(argv[0], "updateMaterialStage") != 0 && std::strcmp(argv[0], "materialState") != 0)
    return false;
  if (!parseInt(argv[1], matTag)) return false;
  for (NDMaterial *m : mats)
    if (m != 0 && m->getTag() == matTag && isSoilTag(m->getClassTag())) return true;
  return false;
}

// D9 (c): MaterialStageParameter reached a particle element. The command passes an
// int (Domain::updateParameter(tag, int), OpenSeesNDMaterialCommands.cpp:371 and
// TclUpdateMaterialStageCommand.cpp), so this registration reads theInt.
// Review item 5: that Parameter is transient -- the command removes and deletes it right after
// the update (OpenSeesNDMaterialCommands.cpp:375-376, TclUpdateMaterialStageCommand.cpp:80) --
// so it is NOT recorded in paramRefs (detachParameters would hold a dangling address, and a
// later Parameter reusing tag and address would be clean()ed by mistake).
int
LadrunoMPMPattern::mpmStageParameter(const char **argv, int argc, Parameter &param)
{
  int matTag = 0;
  if (!soilStageRequest(argv, argc, matTag)) return -2;
  return param.addObject(stageRegId(matTag, true), this);
}

// D7: the one stage walk. Each material registers through its OWN
// setParameter(`materialState`, <its tag>) on a host-local Parameter (never added to
// the Domain), then update(double(stage)): SANISAND id 5 (ManzariDafalias.cpp:977-981,
// Elastic2Plastic on a flip to 1; then LadrunoSANISAND's fast path,
// ladrunoRunStageFlipOnce, LadrunoSANISAND.cpp:5309-5313), HS id 13 (LadrunoHS.cpp:
// 1225-1232), DruckerPrager id 5 -- FE's directly reached Gauss-point path, for every
// particle. mats[p] is the material in both modes (mode A keeps it as a view).
int
LadrunoMPMPattern::stageWalk(int matTag, int stage, std::string &err)
{
  // The value goes through unchanged, as update(double(stage)): each material reads its
  // own states (SANISAND 0/1, HS 0/1, UW DruckerPrager 0/1/2 -- its 1 is still elastic,
  // DruckerPrager.cpp:158-160, :446). No mapping here (orchestrator ruling at the freeze).
  bool any = false, sanisand = false, hs = false;
  for (NDMaterial *m : mats)
    if (m != 0 && m->getTag() == matTag) {
      any = true;
      if (isSanisandTag(m->getClassTag())) sanisand = true;
      if (isHSTag(m->getClassTag())) hs = true;
    }
  if (!any) {
    err = "no particle material has tag " + std::to_string(matTag);
    return -1;
  }
  // AMEND A2 (232k): SANISAND and HS know stages 0 and 1 only (SANISAND's Elastic2Plastic runs
  // on 1 alone); DruckerPrager keeps the pass-through (0/1/2)
  if ((sanisand || hs) && stage != 0 && stage != 1) {
    err = std::string(sanisand ? "LadrunoSANISAND" : "LadrunoHS") + " has stages 0 and 1 only"
          " (got " + std::to_string(stage) + ")";
    return -1;
  }
  // SANISAND's flag is class-wide: the walk reaches every SANISAND particle, of any tag,
  // since all of them flip with it.
  std::vector<int> reach;
  for (int p = 0; p < (int)mats.size(); p++) {
    NDMaterial *m = mats[p];
    if (m == 0) continue;
    if (m->getTag() == matTag || (sanisand && isSanisandTag(m->getClassTag()))) reach.push_back(p);
  }
  const int nv = eng.nv();
  std::vector<double> before(reach.size() * nv, 0.0);
  for (size_t k = 0; k < reach.size(); k++) {
    const Vector &s = mats[reach[k]]->getStress();
    for (int c = 0; c < nv && c < s.Size(); c++) before[k * nv + c] = s(c);
  }
  std::map<int, std::unique_ptr<Parameter>> P;   // one per material tag (argv[1] is checked)
  int answered = 0, answeredTag = 0;
  for (int p : reach) {
    NDMaterial *m = mats[p];
    std::unique_ptr<Parameter> &pp = P[m->getTag()];
    if (!pp) pp.reset(new Parameter(0));
    const std::string ts = std::to_string(m->getTag());
    const char *argv[2] = {"materialState", ts.c_str()};
    if (m->setParameter(argv, 2, *pp) >= 0) {
      answered++;
      if (m->getTag() == matTag) answeredTag++;
    }
  }
  if (answeredTag == 0) {        // never a silent no-op
    err = "no particle material with tag " + std::to_string(matTag) + " answers `materialState`"
          " (it has no stage: LadrunoSANISAND 0/1, LadrunoHS 0/1, DruckerPrager 0/1/2 do)";
    return -1;
  }
  int bad = 0;
  for (auto &kv : P) bad += kv.second->update((double)stage);
  if (bad != 0) {
    err = "a particle material refused the stage update";
    return -1;
  }
  // D8: Elastic2Plastic writes mSigma = mSigma_n = Pmin I where p < Pmin. The host's
  // sigHat / sig are rewritten from the material for every particle whose material stress
  // changed bitwise (mode F's next P2G reads them); nothing else is touched (S8).
  ladruno_mpm::Particles &Pt = eng.P;
  const int nd = eng.ndm(), nn = nd * nd;
  for (size_t k = 0; k < reach.size(); k++) {
    const int p = reach[k];
    const Vector &s = mats[p]->getStress();
    bool changed = false;
    for (int c = 0; c < nv && c < s.Size(); c++) changed = changed || s(c) != before[k * nv + c];
    if (!changed || p >= Pt.size()) continue;
    double sh[6];
    for (int c = 0; c < nv; c++) {
      sh[c] = Pt.sig0[nv * p + c] + s(c);
      Pt.sigHat[nv * p + c] = sh[c];
    }
    if (nd == 3) ladruno_hypo::pushStress6(&Pt.R[nn * p], sh, &Pt.sig[nv * p]);
    else ladruno_hypo::pushStressPlaneStrain(&Pt.R[nn * p], sh, &Pt.sig[nv * p]);
  }
  if (sanisand) sanisandStage = stage;
  if (hs) hsWalked[matTag] = stage;                 // AMEND A1
  return answered;
}

// D10
std::string
LadrunoMPMPattern::stageGuard(void) const
{
  bool any = false;
  for (NDMaterial *m : mats)
    if (m != 0 && isSanisandTag(m->getClassTag())) { any = true; break; }
  if (!any) return "";
  const int f = LadrunoSANISAND::ladrunoStageFlag();
  if (f == sanisandStage) return "";
  const std::string head =
      "LadrunoSANISAND's stage flag (class-wide, every SANISAND tag in the process) is " +
      std::to_string(f) + " but the particles are at stage " + std::to_string(sanisandStage) +
      ": a route that did not reach them changed it (`updateMaterialStage` answered by an FE"
      " element first, mode F's `updateMaterialStage` (W12), another SANISAND tag's flip), or"
      " `mpm reset` returned the particles to stage 0 under a flipped flag (WP-232k D10). ";
  // AMEND A3: never advise a flip the particles did not take; re-stage from their own stage
  if (sanisandStage == 0)
    return head + "Run `mpm stage <tag> 0` (the particles' stage), then stage them again from 0"
                  " as the run requires. The flag is class-wide: this also returns every FE"
                  " LadrunoSANISAND element in the process to stage 0.";
  return head + "Put the flag back by the route that moved it, or `mpm reset` and re-stage the"
                " particles from 0.";
}

// AMEND A4
bool
LadrunoMPMPattern::anySoil(void) const
{
  for (NDMaterial *m : mats)
    if (m != 0 && isSoilTag(m->getClassTag())) return true;
  return false;
}

// AMEND A5 (232k): the WP-232i B10 rule at time t, for a box that may move
std::string
LadrunoMPMPattern::gridvelGhostConflict(double t, bool staticAllTimes) const
{
  if (!ladruno_mpm::isBSpline(eng.par.basis)) return "";
  const ladruno_mpm::Grid &G = eng.grid;
  for (const GridVel &gv : gridVels) {
    // review item 4: a box that never moves is checked up front, also before its `start`
    const bool moves = gv.track >= 0 || gv.translate >= 0;
    if (!(staticAllTimes && !moves) && !gv.presentAt(t)) continue;
    double lo[3], hi[3];
    gv.boxAt(t, lo, hi);
    for (int d = 0; d < G.ndm; d++) {
      const double a = G.x0[d], b = G.x0[d] + G.n[d] * G.h, tol = 1e-9 * G.h;
      if ((eng.fc.fixMask[2 * d] != 0u && lo[d] < a - tol) ||
          (eng.fc.fixMask[2 * d + 1] != 0u && hi[d] > b + tol))
        return "gridvel '" + gv.name + "' reaches the ghost layer behind a `fix` face at t = " +
               std::to_string(t) + ": refused on a B-spline grid (those ghosts are images;"
               " WP-232i B10, checked where a moving box is)";
    }
  }
  return "";
}

// D15
std::string
LadrunoMPMPattern::unsymmetricSoil(void) const
{
  for (NDMaterial *m : mats) {
    if (m == 0) continue;
    const int ct = m->getClassTag();
    if (isSanisandTag(ct)) {
      if (!static_cast<LadrunoSANISAND *>(m)->getLadrunoImplexOptions().enabled)
        return "LadrunoSANISAND " + std::to_string(m->getTag()) + " (without -implex)";
    } else if (isHSTag(ct)) {
      if (!static_cast<LadrunoHS *>(m)->usesImplex())
        return "LadrunoHS " + std::to_string(m->getTag()) + " (without -implex)";
    }
  }
  return "";
}

// D16: sigma_zz. 3D: the particle's sig[2]. Plane strain: the material's getStressZZ()
// (an in-plane rotation leaves it invariant); the Q2 offset has no zz component.
int
LadrunoMPMPattern::particleSigZZ(std::vector<double> &out, std::string &err)
{
  const int N = eng.P.size(), nd = eng.ndm();
  out.assign(N, std::numeric_limits<double>::quiet_NaN());
  if (nd == 3) {
    for (int p = 0; p < N; p++) out[p] = eng.P.sig[6 * p + 2];
    return 0;
  }
  if ((int)mats.size() != N) {
    err = "the particle materials are not available";
    return -1;
  }
  for (int p = 0; p < N; p++)
    if (mats[p] != 0) out[p] = mats[p]->getStressZZ();
  return 0;
}

// D16: a material response per particle, flat in storage order (sizeOnly: the length
// of each particle's response). A particle whose material does not answer gives -1 in
// sizeOnly and contributes nothing.
int
LadrunoMPMPattern::particleResponse(const std::vector<std::string> &args, bool sizeOnly,
                                    std::vector<double> &out, std::string &err)
{
  out.clear();
  if (args.empty()) {
    err = "particleResponse [-size] <response args..>";
    return -1;
  }
  std::vector<const char *> argv;
  for (const std::string &a : args) argv.push_back(a.c_str());
  DummyStream ds;
  for (NDMaterial *m : mats) {
    Response *r = (m != 0) ? m->setResponse(argv.data(), (int)argv.size(), ds) : 0;
    if (r == 0) {
      if (sizeOnly) out.push_back(-1.0);
      continue;
    }
    r->getResponse();
    const Vector &d = r->getInformation().getData();
    if (sizeOnly) out.push_back(d.Size());
    else for (int k = 0; k < d.Size(); k++) out.push_back(d(k));
    delete r;
  }
  return 0;
}

// D16: `mpm info materials` -> [tag, classTag, scheme, implex, stage, count] per tag.
// scheme: SANISAND's IntScheme, -1 otherwise; stage: SANISAND's class-wide flag, -1
// otherwise (HS keeps its stage per tag inside the material).
void
LadrunoMPMPattern::materialsInfo(std::vector<double> &out)
{
  out.clear();
  std::map<int, int> count;
  std::map<int, NDMaterial *> first;
  for (NDMaterial *m : mats)
    if (m != 0) {
      count[m->getTag()]++;
      if (!first.count(m->getTag())) first[m->getTag()] = m;
    }
  for (const auto &kv : first) {
    NDMaterial *m = kv.second;
    const int ct = m->getClassTag();
    double scheme = -1.0, implex = 0.0, stage = -1.0;
    if (isSanisandTag(ct)) {
      LadrunoSANISAND *s = static_cast<LadrunoSANISAND *>(m);
      scheme = s->ladrunoIntegrationScheme();
      implex = s->getLadrunoImplexOptions().enabled ? 1.0 : 0.0;
      stage = LadrunoSANISAND::ladrunoStageFlag();
    } else if (isHSTag(ct)) {
      implex = static_cast<LadrunoHS *>(m)->usesImplex() ? 1.0 : 0.0;
    }
    out.push_back(kv.first);
    out.push_back(ct);
    out.push_back(scheme);
    out.push_back(implex);
    out.push_back(stage);
    out.push_back(count[kv.first]);
  }
}

// D13
std::string
LadrunoMPMPattern::refuserText(int p) const
{
  if (p < 0 || p >= eng.P.size())
    return "not identified (an injected or stale WP-99 count)";
  std::string s = "particle id " + std::to_string((long long)eng.P.id[p]);
  NDMaterial *m = (p < (int)mats.size()) ? mats[p] : 0;
  if (m != 0)
    s += ", nDMaterial " + std::string(m->getClassType()) + " (class tag " +
         std::to_string(m->getClassTag()) + ") tag " + std::to_string(m->getTag());
  return s;
}

// D17, D18 (#1109): mode F's prescribed grid velocities for the step [t, t + h], on the
// node set at the step start -- every extended-lattice node (real or ghost) inside the
// closed box, tolerance 1e-9 h, the rule mode A's prepareStep uses. A `fix` face that
// fixes a component the gridvel prescribes non-zero, or two gridvels with different
// values on one node component, is a malformed case and refused, as in mode A.
int
LadrunoMPMPattern::buildPrescribed(double t, double h, std::string &err)
{
  eng.presc.clear();
  if (gridVels.empty()) return 0;
  {
    const std::string why = gridvelGhostConflict(t);   // AMEND A5
    if (!why.empty()) { err = why; return -1; }
  }
  const ladruno_mpm::Grid &G = eng.grid;
  const int nd = G.ndm;
  const double eps = 1e-9;
  for (const GridVel &gv : gridVels) {
    if (!gv.presentAt(t)) continue;
    double lo[3], hi[3];
    gv.boxAt(t, lo, hi);
    const double fac = gv.stepFactor(t, h);
    int i0[3] = {0, 0, 0}, i1[3] = {0, 0, 0};
    bool empty = false;
    for (int d = 0; d < nd; d++) {
      const double a = std::max((double)-G.ghost, std::ceil((lo[d] - G.x0[d]) / G.h - eps));
      const double b = std::min((double)(G.n[d] + G.ghost), std::floor((hi[d] - G.x0[d]) / G.h + eps));
      if (a > b) { empty = true; break; }
      i0[d] = (int)a;
      i1[d] = (int)b;
    }
    if (empty) continue;
    int ijk[3] = {0, 0, 0};
    for (ijk[2] = i0[2]; ijk[2] <= i1[2]; ijk[2]++)
      for (ijk[1] = i0[1]; ijk[1] <= i1[1]; ijk[1]++)
        for (ijk[0] = i0[0]; ijk[0] <= i1[0]; ijk[0]++) {
          const long ext = G.extId(ijk);
          for (int a = 0; a < nd; a++) {
            if (gv.freeMask & (1u << a)) continue;
            double cx[3] = {0, 0, 0};
            for (int d = 0; d < nd; d++) cx[d] = G.x0[d] + ijk[d] * G.h;
            const double val = gv.valueAt(a, cx, nd, fac);
            for (int f = 0; f < 2 * nd; f++) {
              const int d = f / 2;
              const bool on = (f % 2 == 0) ? ijk[d] <= 0 : ijk[d] >= G.n[d];
              if (on && (eng.fc.fixMask[f] & (1u << a)) && val != 0.0) {
                err = "gridvel '" + gv.name + "' and a `fix` face prescribe different values on"
                      " one node component (a malformed case)";
                eng.presc.clear();
                return -1;
              }
            }
            ladruno_mpm::Engine::Presc &q = eng.presc[ext];
            if (q.mask == 0u) { q.val[0] = q.val[1] = q.val[2] = 0.0; }
            if ((q.mask & (1u << a)) && q.val[a] != val) {
              err = "gridvel '" + gv.name + "' and another gridvel prescribe different values on"
                    " one node component (a malformed case)";
              eng.presc.clear();
              return -1;
            }
            q.mask |= 1u << a;
            q.val[a] = val;
          }
        }
  }
  return 0;
}
