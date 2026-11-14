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
// LadrunoMPMIO — see the header (WP-232c phase c2, 232a contracts 1 and 2).

#include "LadrunoMPMIO.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <random>
#include <set>
#include <sstream>

#ifdef _HDF5
#include <Ladruno_Hdf5.h>
#endif

namespace ladruno_mpm_io {

static const long long MAX_SAFE_ID = (1LL << 53) - 1;

std::string
ParticleWriter::components(const std::string &f, int ndm)
{
  const bool two = (ndm == 2);
  if (f == "POSITION") return two ? "x,y" : "x,y,z";
  if (f == "VELOCITY") return two ? "Vx,Vy" : "Vx,Vy,Vz";
  if (f == "STRESS") return two ? "sigma11,sigma22,sigma12"
                                : "sigma11,sigma22,sigma33,sigma12,sigma23,sigma13";
  if (f == "STRAIN") return two ? "eps11,eps22,gamma12"
                                : "eps11,eps22,eps33,gamma12,gamma23,gamma13";
  if (f == "VOLUME") return "V";
  if (f == "GIMP_HALF_LENGTH") return two ? "lx,ly" : "lx,ly,lz";
  if (f == "DEFORMATION_GRADIENT") return two ? "F11,F12,F21,F22"
                                              : "F11,F12,F13,F21,F22,F23,F31,F32,F33";
  return "";
}

#ifdef _HDF5

bool hdf5Available(void) { return true; }

namespace {

// Silence the HDF5 error stack while probing (a missing attribute is a rule, not
// a library error to print), restoring the caller's handler on exit.
struct QuietH5 {
  H5E_auto2_t func;
  void *data;
  QuietH5() { H5Eget_auto2(H5E_DEFAULT, &func, &data); H5Eset_auto2(H5E_DEFAULT, NULL, NULL); }
  ~QuietH5() { H5Eset_auto2(H5E_DEFAULT, func, data); }
};

bool typeIs(hid_t t, char kind, size_t size)
{
  if (t < 0 || H5Tget_size(t) != size) return false;
  const H5T_class_t c = H5Tget_class(t);
  if (kind == 'f') return c == H5T_FLOAT && H5Tget_order(t) == H5T_ORDER_LE;
  if (c != H5T_INTEGER || H5Tget_order(t) != H5T_ORDER_LE) return false;
  return (kind == 'i') ? (H5Tget_sign(t) == H5T_SGN_2) : (H5Tget_sign(t) == H5T_SGN_NONE);
}

bool hasAttr(hid_t obj, const char *name) { return H5Aexists(obj, name) > 0; }

// Accepts a scalar or a 1-D attribute. Returns false when absent / not numeric.
bool readAttrNum(hid_t obj, const char *name, std::vector<double> &out, bool &isInt, size_t &width)
{
  out.clear();
  if (!hasAttr(obj, name)) return false;
  hid_t a = H5Aopen(obj, name, H5P_DEFAULT);
  if (a < 0) return false;
  hid_t t = H5Aget_type(a), s = H5Aget_space(a);
  const H5T_class_t c = H5Tget_class(t);
  width = H5Tget_size(t);
  isInt = (c == H5T_INTEGER);
  bool ok = (c == H5T_INTEGER || c == H5T_FLOAT);
  const hssize_t n = H5Sget_simple_extent_npoints(s);
  if (ok && n >= 1) {
    if (isInt) {
      std::vector<long long> v((size_t)n);
      ok = H5Aread(a, H5T_NATIVE_LLONG, v.data()) >= 0;
      for (long long x : v) out.push_back((double)x);
    } else {
      out.resize((size_t)n);
      ok = H5Aread(a, H5T_NATIVE_DOUBLE, out.data()) >= 0;
    }
  } else {
    ok = false;
  }
  H5Sclose(s); H5Tclose(t); H5Aclose(a);
  return ok;
}

bool readAttrInt64(hid_t obj, const char *name, long long &out, size_t &width)
{
  if (!hasAttr(obj, name)) return false;
  hid_t a = H5Aopen(obj, name, H5P_DEFAULT);
  hid_t t = H5Aget_type(a), s = H5Aget_space(a);
  width = H5Tget_size(t);
  bool ok = H5Tget_class(t) == H5T_INTEGER && H5Sget_simple_extent_npoints(s) >= 1;
  if (ok) {
    std::vector<long long> v((size_t)H5Sget_simple_extent_npoints(s));
    ok = H5Aread(a, H5T_NATIVE_LLONG, v.data()) >= 0;
    if (ok) out = v[0];
  }
  H5Sclose(s); H5Tclose(t); H5Aclose(a);
  return ok;
}

std::string trimFixed(const char *p, size_t n)
{
  size_t k = 0;
  while (k < n && p[k] != '\0') k++;
  std::string s(p, k);
  while (!s.empty() && s.back() == ' ') s.pop_back();   // SPACEPAD
  return s;
}

// A string attribute, scalar or 1-element, fixed or variable length.
bool readAttrStr(hid_t obj, const char *name, std::string &out)
{
  out.clear();
  if (!hasAttr(obj, name)) return false;
  hid_t a = H5Aopen(obj, name, H5P_DEFAULT);
  hid_t t = H5Aget_type(a), s = H5Aget_space(a);
  bool ok = H5Tget_class(t) == H5T_STRING && H5Sget_simple_extent_npoints(s) >= 1;
  const hssize_t n = ok ? H5Sget_simple_extent_npoints(s) : 0;
  if (ok && H5Tis_variable_str(t) > 0) {
    hid_t mt = H5Tcopy(H5T_C_S1);
    H5Tset_size(mt, H5T_VARIABLE);
    H5Tset_cset(mt, H5Tget_cset(t));
    std::vector<char *> buf((size_t)n, (char *)0);
    ok = H5Aread(a, mt, buf.data()) >= 0;
    if (ok && buf[0]) out = buf[0];
    for (char *p : buf) if (p) H5free_memory(p);
    H5Tclose(mt);
  } else if (ok) {
    const size_t w = H5Tget_size(t);
    std::vector<char> buf(w * (size_t)n + 1, '\0');
    hid_t mt = H5Tcopy(t);
    ok = H5Aread(a, mt, buf.data()) >= 0;
    if (ok) out = trimFixed(buf.data(), w);
    H5Tclose(mt);
  }
  H5Sclose(s); H5Tclose(t); H5Aclose(a);
  return ok;
}

// A dataset: its type, dims and (when wanted) its contents.
struct DS {
  bool exists = false, isDataset = false;
  hid_t type = -1;
  std::vector<hsize_t> dims;
  std::vector<int> filters;
};

DS probe(hid_t grp, const char *name)
{
  DS d;
  if (H5Lexists(grp, name, H5P_DEFAULT) <= 0) return d;
  d.exists = true;
  hid_t ds = H5Dopen2(grp, name, H5P_DEFAULT);
  if (ds < 0) return d;
  d.isDataset = true;
  d.type = H5Dget_type(ds);
  hid_t sp = H5Dget_space(ds);
  const int r = H5Sget_simple_extent_ndims(sp);
  if (r > 0) {
    d.dims.resize((size_t)r);
    H5Sget_simple_extent_dims(sp, d.dims.data(), NULL);
  }
  H5Sclose(sp);
  hid_t pl = H5Dget_create_plist(ds);
  const int nf = H5Pget_nfilters(pl);
  for (int i = 0; i < nf; i++) {
    unsigned flags = 0, cfg = 0;
    size_t nel = 0;
    d.filters.push_back((int)H5Pget_filter2(pl, (unsigned)i, &flags, &nel, NULL, 0, NULL, &cfg));
  }
  H5Pclose(pl);
  H5Dclose(ds);
  return d;
}

template <class T>
bool readAll(hid_t grp, const char *name, hid_t memType, std::vector<T> &out, size_t count)
{
  out.assign(count, T());
  if (count == 0) return true;
  hid_t ds = H5Dopen2(grp, name, H5P_DEFAULT);
  if (ds < 0) return false;
  const bool ok = H5Dread(ds, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()) >= 0;
  H5Dclose(ds);
  return ok;
}

bool readStrings(hid_t grp, const char *name, std::vector<std::string> &out)
{
  out.clear();
  hid_t ds = H5Dopen2(grp, name, H5P_DEFAULT);
  if (ds < 0) return false;
  hid_t t = H5Dget_type(ds), sp = H5Dget_space(ds);
  const hssize_t n = H5Sget_simple_extent_npoints(sp);
  bool ok = H5Tget_class(t) == H5T_STRING && n >= 0;
  if (ok && H5Tis_variable_str(t) > 0) {
    hid_t mt = H5Tcopy(H5T_C_S1);
    H5Tset_size(mt, H5T_VARIABLE);
    H5Tset_cset(mt, H5Tget_cset(t));
    std::vector<char *> buf((size_t)n, (char *)0);
    ok = n == 0 || H5Dread(ds, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()) >= 0;
    for (char *p : buf) { out.push_back(p ? std::string(p) : std::string()); if (p) H5free_memory(p); }
    H5Tclose(mt);
  } else if (ok) {
    const size_t w = H5Tget_size(t);
    std::vector<char> buf(w * (size_t)n + 1, '\0');
    hid_t mt = H5Tcopy(t);
    ok = n == 0 || H5Dread(ds, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()) >= 0;
    for (hssize_t i = 0; i < n; i++) out.push_back(trimFixed(buf.data() + w * (size_t)i, w));
    H5Tclose(mt);
  }
  H5Sclose(sp); H5Tclose(t); H5Dclose(ds);
  return ok;
}

std::vector<std::string> linkNames(hid_t grp)
{
  std::vector<std::string> names;
  H5G_info_t gi;
  if (H5Gget_info(grp, &gi) < 0) return names;
  for (hsize_t i = 0; i < gi.nlinks; i++) {
    char buf[512];
    const ssize_t len = H5Lget_name_by_idx(grp, ".", H5_INDEX_NAME, H5_ITER_INC, i, buf, sizeof(buf), H5P_DEFAULT);
    if (len > 0) names.push_back(std::string(buf, (size_t)std::min<ssize_t>(len, (ssize_t)sizeof(buf) - 1)));
  }
  return names;
}

bool isUuid4(const std::string &s)
{
  if (s.size() != 36) return false;
  for (size_t i = 0; i < 36; i++) {
    const char c = s[i];
    if (i == 8 || i == 13 || i == 18 || i == 23) { if (c != '-') return false; continue; }
    if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return false;
  }
  if (s[14] != '4') return false;
  return s[19] == '8' || s[19] == '9' || s[19] == 'a' || s[19] == 'b';
}

bool allFinite(const std::vector<double> &v)
{
  for (double x : v) if (!std::isfinite(x)) return false;
  return true;
}

}  // namespace

// ---------------------------------------------------------------------------
//  Contract 1
// ---------------------------------------------------------------------------
int
readParticleSet(const std::string &path, ParticleSet &P, std::vector<std::string> &problems)
{
  problems.clear();
  P = ParticleSet();
  P.path = path;
  QuietH5 quiet;
  auto err = [&](const std::string &m) { problems.push_back(m); };

  hid_t f = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (f < 0) {
    err("C1-REQ: cannot open '" + path + "' as an HDF5 file");
    return -1;
  }

  // ---- root attributes (§1.1) -----------------------------------------------
  const char *rootStrings[] = {"FORMAT", "GENERATOR", "SET_ID", "UNIT_LENGTH", "UNIT_FORCE",
                               "UNIT_TIME", "MASS_SOURCE", "STRESS_SIGN", "VOIGT_ORDER"};
  std::string fmt, msrc, sign, voigt;
  for (const char *k : rootStrings)
    if (!hasAttr(f, k)) err(std::string("C1-REQ: root attribute ") + k + " missing");
  for (const char *k : {"FORMAT_VERSION", "NDM", "NUM_PARTICLES"})
    if (!hasAttr(f, k)) err(std::string("C1-REQ: root attribute ") + k + " missing");
  if (readAttrStr(f, "FORMAT", fmt) && fmt != "ladruno-mpm-particles")
    err("C1-FMT: FORMAT '" + fmt + "' != 'ladruno-mpm-particles'");
  std::vector<double> nv;
  bool isInt = false;
  size_t width = 0;
  if (readAttrNum(f, "FORMAT_VERSION", nv, isInt, width)) {
    if (!isInt || width != 4) err("C1-DTYPE: root attribute FORMAT_VERSION must be int32");
    if (nv[0] != 1.0) err("C1-VER: FORMAT_VERSION is not 1");
  }
  if (readAttrStr(f, "GENERATOR", P.generator) && P.generator.empty())
    err("C1-GEN: GENERATOR is empty");
  if (readAttrStr(f, "SET_ID", P.setId) && !isUuid4(P.setId))
    err("C1-SETID: SET_ID '" + P.setId + "' is not a canonical uuid4");
  readAttrStr(f, "UNIT_LENGTH", P.unitLength);
  readAttrStr(f, "UNIT_FORCE", P.unitForce);
  readAttrStr(f, "UNIT_TIME", P.unitTime);
  readAttrStr(f, "UNIT_SYSTEM", P.unitSystem);
  for (const char *u : {"UNIT_LENGTH", "UNIT_FORCE", "UNIT_TIME"}) {
    std::string val;
    if (hasAttr(f, u) && readAttrStr(f, u, val) && val.empty())
      err(std::string("C1-UNIT: ") + u + " is empty");
  }
  if (readAttrStr(f, "STRESS_SIGN", sign) && sign != "tension-positive")
    err("C1-SIGN: STRESS_SIGN '" + sign + "' != 'tension-positive'");
  if (readAttrNum(f, "NDM", nv, isInt, width)) {
    if (!isInt || width != 4) err("C1-DTYPE: root attribute NDM must be int32");
    P.ndm = (int)nv[0];
  }
  if (hasAttr(f, "NUM_PARTICLES")) {
    if (!readAttrInt64(f, "NUM_PARTICLES", P.n, width) || width != 8)
      err("C1-DTYPE: root attribute NUM_PARTICLES must be int64");
  }
  const int ndm = P.ndm;
  if (ndm != 2 && ndm != 3) {
    err("C1-NDM: NDM must be 2 or 3");
    H5Fclose(f);
    return -1;
  }
  const int nst = (ndm == 3) ? 6 : 3;
  if (readAttrStr(f, "VOIGT_ORDER", voigt) && voigt != (ndm == 3 ? "11,22,33,12,23,13" : "11,22,12"))
    err("C1-VOIGT: VOIGT_ORDER '" + voigt + "' does not match NDM");
  const bool hasThk = readAttrNum(f, "PLANE_STRAIN_THICKNESS", nv, isInt, width);
  if (ndm == 2 && (!hasThk || !std::isfinite(nv[0]) || nv[0] <= 0.0))
    err("C1-THK: PLANE_STRAIN_THICKNESS must be finite and > 0 when NDM = 2");
  if (ndm == 2 && hasThk) P.thickness = nv[0];
  if (ndm == 3 && hasAttr(f, "PLANE_STRAIN_THICKNESS"))
    err("C1-THK: PLANE_STRAIN_THICKNESS present when NDM = 3");
  if (P.generator == "apeGmsh") {
    hid_t pr = (H5Lexists(f, "PROVENANCE", H5P_DEFAULT) > 0) ? H5Gopen2(f, "PROVENANCE", H5P_DEFAULT) : -1;
    for (const char *k : {"apegmsh_version", "apegmsh_commit", "session_id"})
      if (pr < 0 || !hasAttr(pr, k))
        err(std::string("C1-PROV: PROVENANCE lacks ") + k + " (required when GENERATOR = apeGmsh)");
    if (pr >= 0) H5Gclose(pr);
  }

  // ---- GRID (§1.2) -----------------------------------------------------------
  bool gridOk = false;
  double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
  if (H5Lexists(f, "GRID", H5P_DEFAULT) <= 0) {
    err("C1-REQ: group GRID missing");
  } else {
    hid_t g = H5Gopen2(f, "GRID", H5P_DEFAULT);
    std::vector<double> org, h, nc;
    bool i1, i2, i3;
    size_t w1, w2, w3;
    const bool a = readAttrNum(g, "ORIGIN", org, i1, w1);
    const bool b = readAttrNum(g, "CELL_SIZE", h, i2, w2);
    const bool c = readAttrNum(g, "NUM_CELLS", nc, i3, w3);
    if (!a || !b || !c) {
      err("C1-GRID: GRID needs ORIGIN, CELL_SIZE and NUM_CELLS");
    } else if ((int)org.size() != ndm || (int)h.size() != ndm || (int)nc.size() != ndm) {
      err("C1-GRID: ORIGIN/CELL_SIZE/NUM_CELLS must each be [NDM]");
    } else if (!i3) {
      err("C1-GRID: NUM_CELLS must be integer");
    } else {
      bool good = true;
      for (int d = 0; d < ndm; d++)
        good = good && std::isfinite(org[d]) && std::isfinite(h[d]) && h[d] > 0 && nc[d] >= 1;
      if (!good) err("C1-GRID: non-finite origin, CELL_SIZE <= 0 or NUM_CELLS < 1");
      if (w3 != 4) err("C1-GRID: GRID/NUM_CELLS must be int32");
      if (good) {
        gridOk = true;
        for (int d = 0; d < ndm; d++) {
          P.origin[d] = org[d]; P.cellSize[d] = h[d]; P.numCells[d] = (int)nc[d];
          lo[d] = org[d]; hi[d] = org[d] + nc[d] * h[d];
        }
      }
    }
    if (!readAttrStr(g, "BASIS", P.basis) || (P.basis != "gimp" && P.basis != "linear"))
      err("C1-BASIS: GRID/BASIS must be gimp|linear");
    H5Gclose(g);
  }

  // ---- BODIES (§1.3) ---------------------------------------------------------
  std::set<int> bodySet;
  if (H5Lexists(f, "BODIES", H5P_DEFAULT) <= 0) {
    err("C1-REQ: BODIES/ID and BODIES/NAME are required");
  } else {
    hid_t g = H5Gopen2(f, "BODIES", H5P_DEFAULT);
    DS di = probe(g, "ID"), dn = probe(g, "NAME");
    if (!di.isDataset || !dn.isDataset) {
      err("C1-REQ: BODIES/ID and BODIES/NAME are required");
    } else {
      if (!typeIs(di.type, 'i', 4)) err("C1-BODY: BODIES/ID dtype != int32");
      const size_t nb = di.dims.empty() ? 0 : (size_t)di.dims[0];
      readAll(g, "ID", H5T_NATIVE_INT, P.bodyId, nb);
      readStrings(g, "NAME", P.bodyName);
      std::set<std::string> nameSet(P.bodyName.begin(), P.bodyName.end());
      bodySet = std::set<int>(P.bodyId.begin(), P.bodyId.end());
      bool neg = false, emptyName = false;
      for (int b : P.bodyId) neg = neg || b < 0;
      for (const std::string &s : P.bodyName) emptyName = emptyName || s.empty();
      if (P.bodyId.size() != P.bodyName.size() || P.bodyId.empty())
        err("C1-BODY: BODIES ID and NAME lengths differ or are empty");
      if (bodySet.size() != P.bodyId.size() || nameSet.size() != P.bodyName.size() || neg || emptyName)
        err("C1-BODY: body ids (>= 0) and names must be unique, names non-empty");
    }
    if (di.type >= 0) H5Tclose(di.type);
    if (dn.type >= 0) H5Tclose(dn.type);
    H5Gclose(g);
  }

  // ---- PARTICLES (§1.4) ------------------------------------------------------
  if (H5Lexists(f, "PARTICLES", H5P_DEFAULT) <= 0) {
    err("C1-REQ: group PARTICLES missing");
    H5Fclose(f);
    return -1;
  }
  hid_t pg = H5Gopen2(f, "PARTICLES", H5P_DEFAULT);
  const long long n = P.n;
  if (n < 1) err("C1-N: NUM_PARTICLES must be >= 1");
  struct Spec { const char *name; char kind; size_t size; int second; bool required; };
  // second: 0 none, 1 ndm, 2 stress components
  const Spec specs[] = {
    {"ID", 'i', 8, 0, true}, {"POSITION", 'f', 8, 1, true}, {"VOLUME", 'f', 8, 0, true},
    {"MASS", 'f', 8, 0, false}, {"DENSITY", 'f', 8, 0, false},
    {"GIMP_HALF_LENGTH", 'f', 8, 1, true}, {"VELOCITY", 'f', 8, 1, false},
    {"STRESS", 'f', 8, 2, false}, {"STRESS33", 'f', 8, 0, false},
    {"BODY_ID", 'i', 4, 0, true}, {"MATERIAL_TAG", 'i', 4, 0, true}, {"FLAGS", 'u', 4, 0, false}};
  for (const std::string &name : linkNames(pg)) {
    bool known = false;
    for (const Spec &s : specs) known = known || name == s.name;
    if (!known) err("C1-UNKNOWN: PARTICLES/" + name + " is not a contract-1 field");
  }
  if (ndm == 3 && H5Lexists(pg, "STRESS33", H5P_DEFAULT) > 0)
    err("C1-STRESS33: STRESS33 is plane-strain only; in 3D sigma_33 is in STRESS");
  std::set<std::string> present;
  bool shapesOk = n >= 1;
  for (const Spec &s : specs) {
    DS d = probe(pg, s.name);
    if (!d.exists) {
      if (s.required) { err(std::string("C1-REQ: PARTICLES/") + s.name + " missing"); shapesOk = false; }
      continue;
    }
    if (!d.isDataset) { err(std::string("C1-SHAPE: PARTICLES/") + s.name + " is not a dataset"); shapesOk = false; continue; }
    present.insert(s.name);
    if (!typeIs(d.type, s.kind, s.size)) err(std::string("C1-DTYPE: PARTICLES/") + s.name + " has the wrong dtype");
    const size_t second = (s.second == 0) ? 0 : (size_t)((s.second == 1) ? ndm : nst);
    const long long rows = d.dims.empty() ? 0 : (long long)d.dims[0];
    if (rows != n) {
      err(std::string("C1-N: PARTICLES/") + s.name + " length differs from NUM_PARTICLES");
      shapesOk = false;
    } else if ((second == 0 && d.dims.size() != 1) ||
               (second > 0 && (d.dims.size() != 2 || d.dims[1] != second))) {
      err(std::string("C1-SHAPE: PARTICLES/") + s.name + " has the wrong shape");
      shapesOk = false;
    }
    for (int fid : d.filters)
      if (fid != 1 && fid != 2 && fid != 3)
        err(std::string("C1-FILTER: PARTICLES/") + s.name + " uses a filter outside shuffle/deflate/fletcher32");
    H5Tclose(d.type);
  }
  std::string msrcVal;
  readAttrStr(f, "MASS_SOURCE", msrcVal);
  if (msrcVal != "mass" && msrcVal != "density") {
    err("C1-MASS: MASS_SOURCE must be 'mass' or 'density'");
  } else {
    const bool hm = present.count("MASS") > 0, hd = present.count("DENSITY") > 0;
    if ((msrcVal == "mass" && (!hm || hd)) || (msrcVal == "density" && (!hd || hm)))
      err("C1-MASS: MASS_SOURCE=" + msrcVal + " needs exactly PARTICLES/" +
          (msrcVal == "mass" ? std::string("MASS") : std::string("DENSITY")));
  }

  if (shapesOk) {
    const size_t N = (size_t)n;
    std::vector<double> dens;
    std::vector<unsigned> flags;
    readAll(pg, "ID", H5T_NATIVE_LLONG, P.id, N);
    readAll(pg, "POSITION", H5T_NATIVE_DOUBLE, P.x, N * ndm);
    readAll(pg, "VOLUME", H5T_NATIVE_DOUBLE, P.volume, N);
    readAll(pg, "GIMP_HALF_LENGTH", H5T_NATIVE_DOUBLE, P.lp, N * ndm);
    readAll(pg, "BODY_ID", H5T_NATIVE_INT, P.body, N);
    readAll(pg, "MATERIAL_TAG", H5T_NATIVE_INT, P.mat, N);
    if (present.count("MASS")) readAll(pg, "MASS", H5T_NATIVE_DOUBLE, P.mass, N);
    if (present.count("DENSITY")) readAll(pg, "DENSITY", H5T_NATIVE_DOUBLE, dens, N);
    if ((P.hasVelocity = present.count("VELOCITY") > 0)) readAll(pg, "VELOCITY", H5T_NATIVE_DOUBLE, P.v, N * ndm);
    else P.v.assign(N * ndm, 0.0);
    if ((P.hasStress = present.count("STRESS") > 0)) readAll(pg, "STRESS", H5T_NATIVE_DOUBLE, P.stress, N * nst);
    else P.stress.assign(N * nst, 0.0);
    if ((P.hasStress33 = present.count("STRESS33") > 0)) readAll(pg, "STRESS33", H5T_NATIVE_DOUBLE, P.stress33, N);
    else P.stress33.assign(N, 0.0);
    if (present.count("FLAGS")) readAll(pg, "FLAGS", H5T_NATIVE_UINT, flags, N);

    for (size_t i = 1; i < N; i++)
      if (P.id[i] <= P.id[i - 1]) { err("C1-ID: PARTICLES/ID not strictly increasing (a repeated or unsorted id)"); break; }
    for (long long id : P.id)
      if (id < 0 || id > MAX_SAFE_ID) { err("C1-ID: PARTICLES/ID outside [0, 2^53 - 1]"); break; }
    if (!allFinite(P.x)) err("C1-FINITE: PARTICLES/POSITION has non-finite values");
    if (!allFinite(P.volume)) err("C1-FINITE: PARTICLES/VOLUME has non-finite values");
    if (!allFinite(P.mass)) err("C1-FINITE: PARTICLES/MASS has non-finite values");
    if (!allFinite(dens)) err("C1-FINITE: PARTICLES/DENSITY has non-finite values");
    if (!allFinite(P.lp)) err("C1-FINITE: PARTICLES/GIMP_HALF_LENGTH has non-finite values");
    if (!allFinite(P.v)) err("C1-FINITE: PARTICLES/VELOCITY has non-finite values");
    if (!allFinite(P.stress)) err("C1-FINITE: PARTICLES/STRESS has non-finite values");
    if (!allFinite(P.stress33)) err("C1-FINITE: PARTICLES/STRESS33 has non-finite values");
    for (double v : P.volume) if (!(v > 0.0)) { err("C1-VOL: PARTICLES/VOLUME must be > 0 (zero or negative is refused)"); break; }
    for (double v : P.mass) if (!(v > 0.0)) { err("C1-MASS: PARTICLES/MASS must be > 0"); break; }
    for (double v : dens) if (!(v > 0.0)) { err("C1-MASS: PARTICLES/DENSITY must be > 0"); break; }
    if (gridOk) {
      bool gimpBad = false, posBad = false;
      double ext = 0.0;
      for (int d = 0; d < ndm; d++) ext = std::max(ext, hi[d] - lo[d]);
      const double tol = 1e-9 * ext;
      for (size_t p = 0; p < N; p++)
        for (int d = 0; d < ndm; d++) {
          const double l = P.lp[ndm * p + d], x = P.x[ndm * p + d];
          if (!(l > 0.0) || l > P.cellSize[d] / 2 * (1 + 1e-12)) gimpBad = true;
          if (x - l < lo[d] - tol || x + l > hi[d] + tol) posBad = true;
        }
      if (gimpBad) err("C1-GIMP: GIMP_HALF_LENGTH must satisfy 0 < l <= CELL_SIZE/2 per axis");
      if (posBad) err("C1-POS: a particle or its GIMP domain lies outside the grid box");
    }
    if (!bodySet.empty())
      for (int b : P.body) if (!bodySet.count(b)) { err("C1-BODYREF: PARTICLES/BODY_ID names a body not in BODIES/ID"); break; }
    for (int m : P.mat) if (m <= 0) { err("C1-MAT: PARTICLES/MATERIAL_TAG must be > 0"); break; }
    for (unsigned fl : flags) if (fl != 0u) { err("C1-FLAGS: FLAGS has a set bit; contract-1 v1 defines none"); break; }
    if (!dens.empty() && P.mass.empty()) {
      P.mass.resize(N);
      for (size_t p = 0; p < N; p++) P.mass[p] = dens[p] * P.volume[p];   // once, at load
    }
  }
  H5Gclose(pg);
  H5Fclose(f);
  return problems.empty() ? 0 : -1;
}

// ---------------------------------------------------------------------------
//  Contract 2 — the writer
// ---------------------------------------------------------------------------
namespace h5w = ladruno::h5;

struct ParticleWriter::Impl {
  WriterConfig cfg;
  int ndm = 0;
  hid_t file = -1, gplist = -1;
  int nextStage = 0, epoch = -1;
  bool stage = false, configured = false;
  std::chrono::steady_clock::time_point lastFlush;
  hid_t dsStep = -1, dsTime = -1, dsOff = -1, dsId = -1;
  std::vector<hid_t> dsData;
  std::vector<int> ncomp;
  long long S = 0;
  int T = 0;

  void closeStage()
  {
    for (hid_t d : {dsStep, dsTime, dsOff, dsId}) if (d >= 0) H5Dclose(d);
    for (hid_t d : dsData) if (d >= 0) H5Dclose(d);
    dsStep = dsTime = dsOff = dsId = -1;
    dsData.clear();
    ncomp.clear();
    stage = false;
    S = 0;
    T = 0;
  }
};

static void
ladrunoMPMRunId(std::string &id, std::string &scope)   // LadrunoRecorder's WP-165 rule
{
  const char *v = std::getenv("LADRUNO_RUN_ID");
  if (v && *v) { id = v; scope = "user"; return; }
  v = std::getenv("SLURM_JOB_ID");
  if (v && *v) {
    id = std::string("slurm-") + v;
    const char *st = std::getenv("SLURM_STEP_ID");
    if (st && *st) id += std::string(".") + st;
    scope = "launcher";
    return;
  }
  v = std::getenv("OMPI_MCA_ess_base_jobid");
  if (v && *v) { id = std::string("ompi-") + v; scope = "launcher"; return; }
  std::stringstream ss;
  ss << "proc-" << std::hex
     << (unsigned long long)std::chrono::system_clock::now().time_since_epoch().count()
     << "-" << (unsigned long long)std::random_device{}();
  id = ss.str();
  scope = "process";
}

ParticleWriter::ParticleWriter() : impl(new Impl) {}
ParticleWriter::~ParticleWriter() { close(); delete impl; }
bool ParticleWriter::isOpen(void) const { return impl->file >= 0; }
bool ParticleWriter::isConfigured(void) const { return impl->configured; }
bool ParticleWriter::stageOpen(void) const { return impl->stage; }
int  ParticleWriter::stageEpoch(void) const { return impl->epoch; }
int  ParticleWriter::numOutputs(void) const { return impl->T; }
const WriterConfig &ParticleWriter::config(void) const { return impl->cfg; }

int
ParticleWriter::open(const WriterConfig &cfg, int ndm, std::string &err)
{
  close();
  if (cfg.path.empty()) { err = "an output path is required"; return -1; }
  impl->cfg = cfg;
  impl->ndm = ndm;
  impl->configured = true;
  return 0;
}

// The file and its INFO group, created at the first stage (see open()).
static int
createFile(ParticleWriter::Impl *impl, std::string &err)
{
  const WriterConfig &cfg = impl->cfg;
  const int ndm = impl->ndm;
  hid_t fcpl = H5Pcreate(H5P_FILE_CREATE);
  H5Pset_link_creation_order(fcpl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
  impl->gplist = H5Pcreate(H5P_GROUP_CREATE);
  H5Pset_link_creation_order(impl->gplist, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
  {
    QuietH5 quiet;
    impl->file = H5Fcreate(cfg.path.c_str(), H5F_ACC_TRUNC, fcpl, H5P_DEFAULT);
  }
  H5Pclose(fcpl);
  if (impl->file < 0) {
    H5Pclose(impl->gplist);
    impl->gplist = -1;
    err = "cannot create '" + cfg.path + "'";
    return -1;
  }
  // INFO, as LadrunoRecorder::initialize writes it, with minor 2: this file
  // carries particles (232a §2.1, ruling 1: only a file with particles is 2).
  hid_t info = H5Gcreate2(impl->file, "INFO", H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  h5w::attribute::write(info, "GENERATOR", std::string("Ladruno"));
  h5w::attribute::write(info, "FORMAT_VERSION", (int)1);
  h5w::attribute::write(info, "FORMAT_MINOR_VERSION", (int)2);
  h5w::attribute::write(info, "SOLVER_NAME", std::string("OpenSees"));
  std::vector<int> sv = {3, 5, 1};
  h5w::attribute::write(info, "SOLVER_VERSION", sv);
  h5w::attribute::write(info, "SPATIAL_DIM", ndm);
  h5w::attribute::write(info, "STORED_PRECISION", std::string("f64"));
  h5w::attribute::write(info, "PARTITION_ID", (int)0);
  h5w::attribute::write(info, "PARTITIONED", (int)0);
  h5w::attribute::write(info, "NUM_PARTITIONS", (int)1);
  std::string rid, rscope;
  ladrunoMPMRunId(rid, rscope);
  h5w::attribute::write(info, "RUN_ID", rid);
  h5w::attribute::write(info, "RUN_ID_SCOPE", rscope);
  hid_t prov = H5Gcreate2(info, "PROVENANCE", H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  H5Gclose(prov);
  H5Gclose(info);
  impl->lastFlush = std::chrono::steady_clock::now();
  return 0;
}

int
ParticleWriter::openStage(const StageMeta &m, int step, double time, std::string &err)
{
  if (!impl->configured) { err = "no output configured"; return -1; }
  if (impl->file < 0 && createFile(impl, err) < 0) return -1;
  impl->closeStage();
  const int ndm = impl->ndm;
  std::stringstream ss;
  ss << "MODEL_STAGE[" << impl->nextStage << "]";
  const std::string sname = ss.str();
  hid_t st = H5Gcreate2(impl->file, sname.c_str(), H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  if (st < 0) { err = "cannot create " + sname; return -1; }
  impl->nextStage++;
  impl->epoch = m.epoch;
  // the recorder's stage skeleton (LadrunoRecorder::writeModel), KIND pinned
  h5w::attribute::write(st, "STEP", step);
  h5w::attribute::write(st, "TIME", time);
  h5w::attribute::write(st, "KIND", std::string(m.mode == "quasi-static" ? "static" : "transient"));
  hid_t model = H5Gcreate2(st, "MODEL", H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  hid_t res = H5Gcreate2(st, "RESULTS", H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  for (const char *g : {"ON_NODES", "ON_ELEMENTS", "ON_DOMAIN", "ON_REGIONS"})
    H5Gclose(H5Gcreate2(res, g, H5P_DEFAULT, impl->gplist, H5P_DEFAULT));
  // ... and its empty-partition branch (writeModelNodes, WP-165 MP-8)
  h5w::attribute::write(st, "EMPTY_PARTITION", (int)1);
  hid_t nodes = H5Gcreate2(model, "NODES", H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  const hsize_t dId[2] = {0, 1}, dCrd[2] = {0, (hsize_t)ndm};
  H5Dclose(h5w::dataset::createEmpty(nodes, "ID", H5T_STD_I32LE, 2, dId));
  H5Dclose(h5w::dataset::createEmpty(nodes, "COORDINATES", H5T_IEEE_F64LE, 2, dCrd));
  H5Gclose(nodes);
  H5Gclose(H5Gcreate2(model, "ELEMENTS", H5P_DEFAULT, impl->gplist, H5P_DEFAULT));
  h5w::attribute::write(st, "PARTICLE_EPOCH", m.epoch);

  // ---- MODEL/PARTICLES (232a §2.4) --------------------------------------------
  hid_t mp = H5Gcreate2(model, "PARTICLES", H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  if (!m.sourceSetIds.empty())
    h5w::attribute::write(mp, "SOURCE_SET_ID", m.sourceSetIds[0]);
  if (m.sourceSetIds.size() > 1) {
    std::string js = "[";
    for (size_t i = 0; i < m.sourceSetIds.size(); i++)
      js += (i ? ",\"" : "\"") + m.sourceSetIds[i] + "\"";
    js += "]";
    h5w::attribute::write(mp, "SOURCE_SET_IDS", js);
  }
  h5w::attribute::write(mp, "PATTERN_TAG", m.patternTag);
  h5w::attribute::write(mp, "MODE", m.mode);
  if (ndm == 2) h5w::attribute::write(mp, "PLANE_STRAIN_THICKNESS", m.thickness);
  h5w::attribute::write(mp, "UNIT_LENGTH", m.unitLength);
  h5w::attribute::write(mp, "UNIT_FORCE", m.unitForce);
  h5w::attribute::write(mp, "UNIT_TIME", m.unitTime);
  h5w::attribute::write(mp, "UNIT_SYSTEM", m.unitSystem);   // absent when "" (the wrapper's rule)
  hid_t grid = H5Gcreate2(mp, "GRID", H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  std::vector<double> org(m.origin, m.origin + ndm), cs(m.cellSize, m.cellSize + ndm);
  std::vector<int> nc(m.numCells, m.numCells + ndm);
  h5w::attribute::write(grid, "ORIGIN", org);
  h5w::attribute::write(grid, "CELL_SIZE", cs);
  h5w::attribute::write(grid, "NUM_CELLS", nc);
  h5w::attribute::write(grid, "BASIS", m.basis);
  H5Gclose(grid);
  hid_t bodies = H5Gcreate2(mp, "BODIES", H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  {
    const hsize_t nb[1] = {(hsize_t)m.bodyId.size()};
    hid_t d = h5w::dataset::createEmpty(bodies, "ID", H5T_STD_I32LE, 1, nb);
    if (!m.bodyId.empty()) H5Dwrite(d, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m.bodyId.data());
    H5Dclose(d);
    H5Dclose(h5w::dataset::createAndWriteStrings(bodies, "NAME", m.bodyName));
  }
  H5Gclose(bodies);
  const int cmp = impl->cfg.compress;
  const hsize_t M = (hsize_t)m.id.size();
  const hsize_t regChunk = std::max<hsize_t>(1, std::min<hsize_t>(1 << 17, M));
  hid_t d;
  bool regOk = true;
  d = h5w::dataset::createRows(mp, "ID", H5T_STD_I64LE, 0, regChunk, cmp);
  regOk = regOk && d >= 0 && h5w::dataset::appendRows(d, H5T_NATIVE_LLONG, m.id.data(), M, 0) >= 0; if (d >= 0) H5Dclose(d);
  d = h5w::dataset::createRows(mp, "BODY_ID", H5T_STD_I32LE, 0, regChunk, cmp);
  regOk = regOk && d >= 0 && h5w::dataset::appendRows(d, H5T_NATIVE_INT, m.body.data(), M, 0) >= 0; if (d >= 0) H5Dclose(d);
  d = h5w::dataset::createRows(mp, "MATERIAL_TAG", H5T_STD_I32LE, 0, regChunk, cmp);
  regOk = regOk && d >= 0 && h5w::dataset::appendRows(d, H5T_NATIVE_INT, m.mat.data(), M, 0) >= 0; if (d >= 0) H5Dclose(d);
  d = h5w::dataset::createRows(mp, "MASS", H5T_IEEE_F64LE, 0, regChunk, cmp);
  regOk = regOk && d >= 0 && h5w::dataset::appendRows(d, H5T_NATIVE_DOUBLE, m.mass.data(), M, 0) >= 0; if (d >= 0) H5Dclose(d);
  H5Gclose(mp);

  // ---- RESULTS/ON_PARTICLES (232a §2.5) -----------------------------------------
  hid_t op = H5Gcreate2(res, "ON_PARTICLES", H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
  h5w::attribute::write(op, "PARTITION_REDUCTION", std::string("NONE"));
  // Crash layout (c2 review 2; quirk WP-232-hdf5-unflushed-file-unreadable-after-kill):
  // the writer flushes after every output by default, and lays the data out so a
  // flushed output is never rewritten. The three index datasets are tiny and stay
  // UNcompressed (an in-place rewrite of their last chunk keeps the earlier entries
  // valid, where a half-rewritten deflated chunk fails its filter); ID and every
  // DATA get one chunk per output (M rows) while it fits 16 MiB, so an output
  // appends new chunks and never touches an old one.
  impl->dsStep = h5w::dataset::createRows(op, "STEP", H5T_STD_I32LE, 0, 1024, 0);
  impl->dsTime = h5w::dataset::createRows(op, "TIME", H5T_IEEE_F64LE, 0, 1024, 0);
  impl->dsOff = h5w::dataset::createRows(op, "OFFSETS", H5T_STD_I64LE, 0, 1024, 0);
  const long long zero = 0;
  const bool offOk = impl->dsOff >= 0 && h5w::dataset::appendRows(impl->dsOff, H5T_NATIVE_LLONG, &zero, 1, 0) >= 0;
  auto perOutput = [&](int nc) -> hsize_t {
    const hsize_t one = std::max<hsize_t>(1, M);
    return (one * (hsize_t)nc * 8 <= (hsize_t)(16 << 20)) ? one
                                                           : std::max<hsize_t>(1, (hsize_t)((1 << 20) / (8 * nc)));
  };
  impl->dsId = h5w::dataset::createRows(op, "ID", H5T_STD_I64LE, 0, perOutput(1), cmp);
  struct Meta { const char *name, *display, *dim; int dataType; };
  static const Meta metas[] = {
    {"POSITION", "Position", "L", 1}, {"VELOCITY", "Velocity", "L/t", 1},
    {"STRESS", "Stress", "F/L^2", 2}, {"STRAIN", "Strain (unrotated, engineering shear)", "", 2},
    {"VOLUME", "Volume", "L^3", 0}, {"GIMP_HALF_LENGTH", "GIMP half-length", "L", 1},
    {"DEFORMATION_GRADIENT", "Deformation gradient", "", 2}};
  for (const std::string &fname : impl->cfg.fields) {
    const Meta *mt = 0;
    for (const Meta &x : metas) if (fname == x.name) mt = &x;
    const std::string comps = components(fname, ndm);
    int nc2 = 1;
    for (char c : comps) if (c == ',') nc2++;
    hid_t fg = H5Gcreate2(op, fname.c_str(), H5P_DEFAULT, impl->gplist, H5P_DEFAULT);
    h5w::attribute::write(fg, "DISPLAY_NAME", std::string(mt ? mt->display : fname.c_str()));
    h5w::attribute::write(fg, "COMPONENTS", comps);
    h5w::attribute::write(fg, "DIMENSION", std::string(mt ? mt->dim : ""));
    h5w::attribute::write(fg, "DESCRIPTION", std::string(""));
    h5w::attribute::write(fg, "TYPE", (int)0);
    h5w::attribute::write(fg, "DATA_TYPE", mt ? mt->dataType : 0);
    h5w::attribute::write(fg, "PARTITION_REDUCTION", std::string("NONE"));
    impl->dsData.push_back(h5w::dataset::createRows(fg, "DATA", H5T_IEEE_F64LE, (hsize_t)nc2, perOutput(nc2), cmp));
    impl->ncomp.push_back(nc2);
    H5Gclose(fg);
  }
  H5Gclose(op);
  H5Gclose(res);
  H5Gclose(model);
  H5Gclose(st);
  if (!regOk || !offOk) {
    err = "cannot write the particle registry or OFFSETS of " + sname;
    impl->closeStage();
    return -1;
  }
  impl->stage = true;
  return 0;
}

int
ParticleWriter::writeOutput(int step, double time, const std::vector<long long> &ids,
                            const std::vector<std::vector<double> > &rows, std::string &err)
{
  if (!impl->stage) { err = "no stage open"; return -1; }
  if (rows.size() != impl->dsData.size()) { err = "field count mismatch"; return -1; }
  const hsize_t n = (hsize_t)ids.size();
  // append order (232a §2.5): DATA, ID, STEP and TIME, OFFSETS last
  for (size_t k = 0; k < rows.size(); k++) {
    if (rows[k].size() != n * (hsize_t)impl->ncomp[k]) { err = "field rows mismatch"; return -1; }
    if (h5w::dataset::appendRows(impl->dsData[k], H5T_NATIVE_DOUBLE, rows[k].data(), n,
                                 (hsize_t)impl->ncomp[k]) < 0) { err = "cannot append DATA"; return -1; }
  }
  if (h5w::dataset::appendRows(impl->dsId, H5T_NATIVE_LLONG, ids.data(), n, 0) < 0) { err = "cannot append ID"; return -1; }
  if (h5w::dataset::appendRows(impl->dsStep, H5T_NATIVE_INT, &step, 1, 0) < 0 ||
      h5w::dataset::appendRows(impl->dsTime, H5T_NATIVE_DOUBLE, &time, 1, 0) < 0) {
    err = "cannot append STEP / TIME";
    return -1;
  }
  impl->S += (long long)n;
  const long long off = impl->S;
  // WRITE BARRIER (232a §2.5). Without SWMR, HDF5 does not reach the disk in API
  // order: a chunk or metadata eviction can put this output's OFFSETS entry on disk
  // before its ID / DATA chunks, and a crash reader would then count an output whose
  // rows are fill. Flushing here, before the OFFSETS append, means a committed
  // OFFSETS entry always has its rows on disk; the caller's flush after the output
  // then commits the entry itself.
  if (H5Fflush(impl->file, H5F_SCOPE_GLOBAL) < 0) { err = "cannot flush before the OFFSETS commit"; return -1; }
  if (h5w::dataset::appendRows(impl->dsOff, H5T_NATIVE_LLONG, &off, 1, 0) < 0) { err = "cannot append OFFSETS"; return -1; }
  impl->T++;
  return 0;
}

void
ParticleWriter::flush(void)
{
  if (impl->file >= 0) H5Fflush(impl->file, H5F_SCOPE_GLOBAL);
  impl->lastFlush = std::chrono::steady_clock::now();
}

void
ParticleWriter::flushIfDue(void)
{
  if (impl->file < 0) return;
  const double since = std::chrono::duration<double>(std::chrono::steady_clock::now() - impl->lastFlush).count();
  if (impl->cfg.flushSeconds <= 0.0 || since >= impl->cfg.flushSeconds) flush();
}

void
ParticleWriter::close(void)
{
  impl->closeStage();
  if (impl->file >= 0) { H5Fflush(impl->file, H5F_SCOPE_GLOBAL); H5Fclose(impl->file); }
  if (impl->gplist >= 0) H5Pclose(impl->gplist);
  impl->file = impl->gplist = -1;
  impl->configured = false;
}

#else  // no HDF5 in this build: refuse loudly

bool hdf5Available(void) { return false; }

int
readParticleSet(const std::string &, ParticleSet &, std::vector<std::string> &problems)
{
  problems.assign(1, "this build has no HDF5 (_HDF5): `mpm particles` is unavailable");
  return -1;
}

struct ParticleWriter::Impl { WriterConfig cfg; };
ParticleWriter::ParticleWriter() : impl(new Impl) {}
ParticleWriter::~ParticleWriter() { delete impl; }
bool ParticleWriter::isOpen(void) const { return false; }
bool ParticleWriter::isConfigured(void) const { return false; }
void ParticleWriter::flushIfDue(void) {}
bool ParticleWriter::stageOpen(void) const { return false; }
int  ParticleWriter::stageEpoch(void) const { return -1; }
int  ParticleWriter::numOutputs(void) const { return 0; }
const WriterConfig &ParticleWriter::config(void) const { return impl->cfg; }
int  ParticleWriter::open(const WriterConfig &, int, std::string &err)
{ err = "this build has no HDF5 (_HDF5): `mpm output` is unavailable"; return -1; }
int  ParticleWriter::openStage(const StageMeta &, int, double, std::string &err)
{ err = "no HDF5"; return -1; }
int  ParticleWriter::writeOutput(int, double, const std::vector<long long> &,
                                 const std::vector<std::vector<double> > &, std::string &err)
{ err = "no HDF5"; return -1; }
void ParticleWriter::flush(void) {}
void ParticleWriter::close(void) {}

#endif

}  // namespace ladruno_mpm_io
