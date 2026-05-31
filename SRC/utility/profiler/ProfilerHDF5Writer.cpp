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
** ****************************************************************** */

// Ladruno Stack Profiler — Phase P1 (HDF5 writer implementation).
//
// File: SRC/utility/profiler/ProfilerHDF5Writer.cpp
//
// Written: Ladruno / nmora
// Created: 2026-05
//
// Description: The ONLY profiler translation unit that #include <hdf5.h>.
//   Translates the profiler's plain in-memory result structs (ProfileNode tree,
//   RunMeta, Series, MemorySnapshot) into the on-disk HDF5 contract defined
//   authoritatively by Ladruno_tools/profiler_viewer/profiler_schema.py.
//
//   Every group name, attribute name, and compound-dtype field below is quoted
//   against the matching line in profiler_schema.py so the two writers stay in
//   lockstep. NANOSECOND integer time units; immutable runs; H5Fflush before
//   close; NEVER exit() — failure returns false and leaves the in-memory tree
//   intact (v1 = write-on-report).
//
//   HDF5 ISOLATION: the public header (ProfilerHDF5Writer.h) is HDF5-free and
//   stores the file handle as `long long`. This .cpp recovers it via
//   `hid_t fid = (hid_t)file_;`. hid_t is int64_t on every supported build.

#include "ProfilerHDF5Writer.h"

#include <hdf5.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ops_profiler {

// schema_version string written at root (profiler_schema.py:47, :84).
static const char* const PROFILER_SCHEMA_VERSION = "0.1.0";

// ===========================================================================
// Low-level attribute helpers. These mirror MPCORecorder.cpp:934-1001
// (h5::attribute::writei/writed/writes) but with two deliberate differences
// for this contract:
//   (1) integer run/node/memory attrs are written as 64-bit (H5T_STD_I64LE),
//       because the schema's int attrs are conceptually i8 and the reader does
//       not pin width (profiler_schema.py:103 note);
//   (2) string attrs are written UNCONDITIONALLY, even when empty: the
//       reference writer always sets all 7 RUN_ATTRS_STR (profiler_schema.py:90-91
//       `g.attrs[k] = str(...)`), so an empty value becomes a 1-byte null string
//       rather than a skipped attribute.
// All return a herr_t (<0 => error); callers fold these into a running status.
// ===========================================================================
namespace {

// SCALAR variable-length UTF-8 string attr — byte-for-byte the same on-disk
// form the Python reference writer produces (h5py default: `g.attrs[k] = str`
// emits a rank-0 variable-length H5T_C_S1/UTF8 attribute). This is REQUIRED for
// parity: a fixed-length, rank-1 (dim-1) attr (the MPCO `writes` shape) reads
// back through h5py as a 1-element numpy array of bytes (e.g. [b'0.1.0']) rather
// than the scalar str the reader (profiler_results.py) expects. An empty string
// is still written (always-present semantics, matching profiler_schema.py:90-91).
herr_t writeStrAttr(hid_t obj, const char* name, const std::string& val)
{
    hid_t space = H5Screate(H5S_SCALAR);            // rank-0 => scalar, not [1]
    if (space < 0) return -1;

    hid_t atype = H5Tcopy(H5T_C_S1);
    if (atype < 0) { H5Sclose(space); return -1; }

    herr_t st = H5Tset_size(atype, H5T_VARIABLE);   // variable-length
    if (st >= 0) st = H5Tset_cset(atype, H5T_CSET_UTF8);

    hid_t attr = -1;
    if (st >= 0) {
        attr = H5Acreate(obj, name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
        if (attr < 0) st = -1;
    }
    if (st >= 0) {
        // For a variable-length string attr, H5Awrite takes a pointer to the
        // `char*` element (one element here, the scalar).
        const char* cstr = val.c_str();
        st = H5Awrite(attr, atype, &cstr);
    }

    if (attr  >= 0) H5Aclose(attr);
    H5Tclose(atype);
    H5Sclose(space);
    return st;
}

// 64-bit integer scalar attr (H5T_STD_I64LE on disk, native LLONG in memory).
herr_t writeI64Attr(hid_t obj, const char* name, int64_t val)
{
    hid_t space = H5Screate(H5S_SCALAR);   // scalar (matches h5py g.attrs[k]=int)
    if (space < 0) return -1;
    hid_t attr = H5Acreate(obj, name, H5T_STD_I64LE, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) { H5Sclose(space); return -1; }
    long long tmp = static_cast<long long>(val);
    herr_t st = H5Awrite(attr, H5T_NATIVE_LLONG, &tmp);
    H5Aclose(attr);
    H5Sclose(space);
    return st;
}

// double scalar attr (H5T_IEEE_F64LE on disk, native double in memory).
herr_t writeF64Attr(hid_t obj, const char* name, double val)
{
    hid_t space = H5Screate(H5S_SCALAR);   // scalar (matches h5py g.attrs[k]=float)
    if (space < 0) return -1;
    hid_t attr = H5Acreate(obj, name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) { H5Sclose(space); return -1; }
    herr_t st = H5Awrite(attr, H5T_NATIVE_DOUBLE, &val);
    H5Aclose(attr);
    H5Sclose(space);
    return st;
}

// -------------------------------------------------------------------------
// Compound dtype factories. Field ORDER and NAMES are load-bearing — they back
// the numpy dtypes in profiler_schema.py and must be inserted in schema order
// with LE std types at HOFFSET. mem_t == file_t for both (round-trip exact).
// -------------------------------------------------------------------------

// ELEM_BY_TYPE_DT (profiler_schema.py:50-56):
//   ("classTag" i4) ("count" i8) ("wall_ns" i8) ("wall_ns_per_elem" f8) ("fb_coupled" u1)
hid_t mkElemByTypeDtype()
{
    hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(ElemByType));
    if (t < 0) return -1;
    herr_t st = 0;
    if (st >= 0) st = H5Tinsert(t, "classTag",         HOFFSET(ElemByType, classTag),         H5T_STD_I32LE);
    if (st >= 0) st = H5Tinsert(t, "count",            HOFFSET(ElemByType, count),            H5T_STD_I64LE);
    if (st >= 0) st = H5Tinsert(t, "wall_ns",          HOFFSET(ElemByType, wall_ns),          H5T_STD_I64LE);
    if (st >= 0) st = H5Tinsert(t, "wall_ns_per_elem", HOFFSET(ElemByType, wall_ns_per_elem), H5T_IEEE_F64LE);
    if (st >= 0) st = H5Tinsert(t, "fb_coupled",       HOFFSET(ElemByType, fb_coupled),       H5T_STD_U8LE);
    if (st < 0) { H5Tclose(t); return -1; }
    return t;
}

// COMPONENTS_LIVE_DT (profiler_schema.py:58-61):
//   ("classTag" i4) ("count" i8)
hid_t mkComponentsLiveDtype()
{
    hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(ComponentLive));
    if (t < 0) return -1;
    herr_t st = 0;
    if (st >= 0) st = H5Tinsert(t, "classTag", HOFFSET(ComponentLive, classTag), H5T_STD_I32LE);
    if (st >= 0) st = H5Tinsert(t, "count",    HOFFSET(ComponentLive, count),    H5T_STD_I64LE);
    if (st < 0) { H5Tclose(t); return -1; }
    return t;
}

// -------------------------------------------------------------------------
// 1-D contiguous-or-chunked numeric dataset helper (series scalars). Chunk =
// min(n,4096) when n>0, else contiguous (profiler_schema.py:102-104). file_t
// and mem_t are the caller-chosen LE/native pair. data may be nullptr only when
// count==0.
// -------------------------------------------------------------------------
herr_t writeDataset1D(hid_t loc, const char* name,
                      hid_t file_t, hid_t mem_t,
                      const void* data, hsize_t n)
{
    hsize_t dims[1] = { n };
    hid_t   space   = H5Screate_simple(1, dims, NULL);
    if (space < 0) return -1;

    hid_t dcpl = H5P_DEFAULT;
    if (n > 0) {
        dcpl = H5Pcreate(H5P_DATASET_CREATE);
        if (dcpl < 0) { H5Sclose(space); return -1; }
        hsize_t chunk[1] = { std::min<hsize_t>(n, 4096) };
        if (H5Pset_chunk(dcpl, 1, chunk) < 0) {
            H5Pclose(dcpl); H5Sclose(space); return -1;
        }
    }

    hid_t ds = H5Dcreate(loc, name, file_t, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    herr_t st = (ds < 0) ? -1 : 0;
    if (st >= 0 && n > 0)
        st = H5Dwrite(ds, mem_t, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    if (ds   >= 0)            H5Dclose(ds);
    if (dcpl != H5P_DEFAULT)  H5Pclose(dcpl);
    H5Sclose(space);
    return st;
}

// -------------------------------------------------------------------------
// Recursively write one rollup node => one group named node.name, carrying the
// 5 int64 NODE_ATTRS_INT and the optional elem_by_type compound dataset, then
// recursing over children in VECTOR ORDER (the stable diff key, P0#6).
// profiler_schema.py:118-126 (_write_node).
// -------------------------------------------------------------------------
bool writeNode(hid_t parent, const ProfileNode& node)
{
    hid_t g = H5Gcreate(parent, node.name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0) return false;

    herr_t st = 0;
    // NODE_ATTRS_INT (profiler_schema.py:68): calls, wall_ns, wall_ns_min, wall_ns_max, cpu_ns
    if (st >= 0) st = writeI64Attr(g, "calls",       node.calls);
    if (st >= 0) st = writeI64Attr(g, "wall_ns",     node.wall_ns);
    if (st >= 0) st = writeI64Attr(g, "wall_ns_min", node.wall_ns_min);
    if (st >= 0) st = writeI64Attr(g, "wall_ns_max", node.wall_ns_max);
    if (st >= 0) st = writeI64Attr(g, "cpu_ns",      node.cpu_ns);
    if (st < 0) { H5Gclose(g); return false; }

    // optional elem_by_type dataset (profiler_schema.py:122-124) — only when non-empty.
    if (!node.elem_by_type.empty()) {
        hid_t dt = mkElemByTypeDtype();
        if (dt < 0) { H5Gclose(g); return false; }
        hsize_t dims[1] = { static_cast<hsize_t>(node.elem_by_type.size()) };
        hid_t   space   = H5Screate_simple(1, dims, NULL);
        if (space < 0) { H5Tclose(dt); H5Gclose(g); return false; }
        hid_t ds = H5Dcreate(g, "elem_by_type", dt, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        herr_t dst = (ds < 0) ? -1 : 0;
        if (dst >= 0)
            dst = H5Dwrite(ds, dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, node.elem_by_type.data());
        if (ds >= 0) H5Dclose(ds);
        H5Sclose(space);
        H5Tclose(dt);
        if (dst < 0) { H5Gclose(g); return false; }
    }

    // recurse children in INSERTION ORDER (vector order == diff path order).
    for (const ProfileNode& child : node.children) {
        if (!writeNode(g, child)) { H5Gclose(g); return false; }
    }

    H5Gclose(g);
    return true;
}

// -------------------------------------------------------------------------
// /series group (profiler_schema.py:99-108). Six chunked 1-D scalar datasets
// + a 2-D float32 wall_ms[nSteps,nPhase] carrying a variable-length UTF-8
// string-array attr "phases". Caller guarantees series && series->nSteps()>0.
// -------------------------------------------------------------------------
bool writeSeries(hid_t runGroup, const Series& s)
{
    hid_t sg = H5Gcreate(runGroup, "series", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (sg < 0) return false;

    const hsize_t n = static_cast<hsize_t>(s.nSteps());

    herr_t st = 0;
    // SERIES_SCALAR (profiler_schema.py:70-73), in that key order.
    if (st >= 0) st = writeDataset1D(sg, "step",           H5T_STD_I64LE,  H5T_NATIVE_LLONG,  s.step.data(),           n);
    if (st >= 0) st = writeDataset1D(sg, "t",              H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, s.t.data(),              n);
    if (st >= 0) st = writeDataset1D(sg, "dt",             H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, s.dt.data(),             n);
    if (st >= 0) st = writeDataset1D(sg, "iters",          H5T_STD_I32LE,  H5T_NATIVE_INT,    s.iters.data(),          n);
    if (st >= 0) st = writeDataset1D(sg, "mem_live_bytes", H5T_STD_I64LE,  H5T_NATIVE_LLONG,  s.mem_live_bytes.data(), n);
    if (st >= 0) st = writeDataset1D(sg, "mem_peak_bytes", H5T_STD_I64LE,  H5T_NATIVE_LLONG,  s.mem_peak_bytes.data(), n);
    if (st < 0) { H5Gclose(sg); return false; }

    // wall_ms[nSteps,nPhase] float32, chunk {min(n,4096),nPhase} (profiler_schema.py:105-107).
    const hsize_t nPhase = static_cast<hsize_t>(s.nPhase());
    {
        hsize_t dims[2]  = { n, nPhase };
        hid_t   space    = H5Screate_simple(2, dims, NULL);
        if (space < 0) { H5Gclose(sg); return false; }

        // Chunking requires every chunk dim <= the fixed dataspace extent. When
        // nPhase==0 (the P1 reality until phase seams land in P2) or n==0, fall back
        // to contiguous storage — extending the Python writer's "chunks=... if n else
        // None" intent to the nPhase==0 case (correctness review C1).
        const bool chunked = (n > 0 && nPhase > 0);
        hid_t dcpl = H5P_DEFAULT;
        if (chunked) {
            dcpl = H5Pcreate(H5P_DATASET_CREATE);
            if (dcpl < 0) { H5Sclose(space); H5Gclose(sg); return false; }
            hsize_t chunk[2] = { std::min<hsize_t>(n, 4096), nPhase };
            if (H5Pset_chunk(dcpl, 2, chunk) < 0) {
                H5Pclose(dcpl); H5Sclose(space); H5Gclose(sg); return false;
            }
        }

        hid_t ds = H5Dcreate(sg, "wall_ms", H5T_IEEE_F32LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
        herr_t dst = (ds < 0) ? -1 : 0;
        if (dst >= 0 && n > 0 && nPhase > 0)
            dst = H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, s.wall_ms.data());

        // attr "phases" = variable-length UTF-8 string array (profiler_schema.py:108).
        if (dst >= 0 && ds >= 0) {
            hid_t vstr = H5Tcopy(H5T_C_S1);
            if (vstr < 0) dst = -1;
            if (dst >= 0) dst = H5Tset_size(vstr, H5T_VARIABLE);
            if (dst >= 0) dst = H5Tset_cset(vstr, H5T_CSET_UTF8);

            hsize_t adim[1] = { nPhase };
            hid_t   aspace  = -1;
            if (dst >= 0) {
                aspace = H5Screate_simple(1, adim, NULL);
                if (aspace < 0) dst = -1;
            }
            hid_t attr = -1;
            if (dst >= 0) {
                attr = H5Acreate(ds, "phases", vstr, aspace, H5P_DEFAULT, H5P_DEFAULT);
                if (attr < 0) dst = -1;
            }
            if (dst >= 0 && !s.phases.empty()) {
                // HDF5 1.14 H5Awrite rejects a NULL buf even for a 0-element space
                // ("buf parameter can't be NULL"). When nPhase==0 (the default P1
                // state) the attribute is created with 0 extent and simply left
                // unwritten — the reader's attrs.get("phases", []) yields [].
                std::vector<const char*> ptrs;
                ptrs.reserve(s.phases.size());
                for (const std::string& p : s.phases) ptrs.push_back(p.c_str());
                dst = H5Awrite(attr, vstr, ptrs.data());
            }
            if (attr   >= 0) H5Aclose(attr);
            if (aspace >= 0) H5Sclose(aspace);
            if (vstr   >= 0) H5Tclose(vstr);
        }

        if (ds >= 0) H5Dclose(ds);
        if (dcpl != H5P_DEFAULT) H5Pclose(dcpl);
        H5Sclose(space);
        if (dst < 0) { H5Gclose(sg); return false; }
    }

    H5Gclose(sg);
    return true;
}

// -------------------------------------------------------------------------
// /memory group (profiler_schema.py:110-115). ALWAYS written: 4 int64 attrs
// plus the components_live compound dataset (a zero-length dataset is legal;
// H5Dwrite is skipped when empty).
// -------------------------------------------------------------------------
bool writeMemory(hid_t runGroup, const MemorySnapshot& mem)
{
    hid_t mg = H5Gcreate(runGroup, "memory", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (mg < 0) return false;

    herr_t st = 0;
    // attrs (profiler_schema.py:112): matrix_live, vector_live, id_live, peak_bytes
    if (st >= 0) st = writeI64Attr(mg, "matrix_live", mem.matrix_live);
    if (st >= 0) st = writeI64Attr(mg, "vector_live", mem.vector_live);
    if (st >= 0) st = writeI64Attr(mg, "id_live",     mem.id_live);
    if (st >= 0) st = writeI64Attr(mg, "peak_bytes",  mem.peak_bytes);
    if (st < 0) { H5Gclose(mg); return false; }

    // components_live compound dataset (profiler_schema.py:114-115), always created.
    {
        hid_t dt = mkComponentsLiveDtype();
        if (dt < 0) { H5Gclose(mg); return false; }
        hsize_t dims[1] = { static_cast<hsize_t>(mem.components_live.size()) };
        hid_t   space   = H5Screate_simple(1, dims, NULL);
        if (space < 0) { H5Tclose(dt); H5Gclose(mg); return false; }
        hid_t ds = H5Dcreate(mg, "components_live", dt, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        herr_t dst = (ds < 0) ? -1 : 0;
        if (dst >= 0 && !mem.components_live.empty())
            dst = H5Dwrite(ds, dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem.components_live.data());
        if (ds >= 0) H5Dclose(ds);
        H5Sclose(space);
        H5Tclose(dt);
        if (dst < 0) { H5Gclose(mg); return false; }
    }

    H5Gclose(mg);
    return true;
}

// Filesystem existence check (decides H5Fopen RDWR vs H5Fcreate EXCL).
bool fileExistsOnDisk(const std::string& path)
{
    struct stat st;
    return ::stat(path.c_str(), &st) == 0;
}

} // anonymous namespace

// ===========================================================================
// ProfilerHDF5Writer public surface
// ===========================================================================

ProfilerHDF5Writer::~ProfilerHDF5Writer()
{
    close();
}

bool ProfilerHDF5Writer::open(const std::string& path)
{
    if (isOpen()) return true;

    // Silence the default HDF5 error stack printer once; we report via return
    // codes and never let HDF5 abort/exit the process (heeds the MPCO exit(-1)
    // incident: the profiler must never kill a running kernel).
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    hid_t fid;
    if (fileExistsOnDisk(path)) {
        fid = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    } else {
        fid = H5Fcreate(path.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (fid < 0) return false;

    // schema_version setdefault: write only if absent (profiler_schema.py:84).
    htri_t has = H5Aexists(fid, "schema_version");
    if (has == 0) {
        if (writeStrAttr(fid, "schema_version", PROFILER_SCHEMA_VERSION) < 0) {
            H5Fclose(fid);
            return false;
        }
    } else if (has < 0) {
        H5Fclose(fid);
        return false;
    }

    file_ = static_cast<long long>(fid);
    return true;
}

void ProfilerHDF5Writer::close()
{
    if (!isOpen()) return;
    hid_t fid = static_cast<hid_t>(file_);
    H5Fflush(fid, H5F_SCOPE_GLOBAL);
    H5Fclose(fid);
    file_ = -1;
}

bool ProfilerHDF5Writer::writeRun(const std::string&    run_id,
                                  const ProfileNode&    rollup,
                                  const RunMeta&        meta,
                                  const Series*         series,
                                  const MemorySnapshot& memory)
{
    if (!isOpen()) return false;
    hid_t fid = static_cast<hid_t>(file_);

    // /runs group (require: open if present else create). profiler_schema.py:85.
    hid_t runs;
    if (H5Lexists(fid, "runs", H5P_DEFAULT) > 0) {
        runs = H5Gopen(fid, "runs", H5P_DEFAULT);
    } else {
        runs = H5Gcreate(fid, "runs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (runs < 0) return false;

    // Immutability: fail if /runs/<run_id> already exists (profiler_schema.py:86-87).
    htri_t exists = H5Lexists(runs, run_id.c_str(), H5P_DEFAULT);
    if (exists != 0) {            // >0 = present (immutable); <0 = error
        H5Gclose(runs);
        return false;
    }

    hid_t g = H5Gcreate(runs, run_id.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0) { H5Gclose(runs); return false; }

    herr_t st = 0;

    // RUN_ATTRS_STR (profiler_schema.py:64, :90-91) — all 7, always written.
    if (st >= 0) st = writeStrAttr(g, "model",      meta.model);
    if (st >= 0) st = writeStrAttr(g, "engine_sha", meta.engine_sha);
    if (st >= 0) st = writeStrAttr(g, "integrator", meta.integrator);
    if (st >= 0) st = writeStrAttr(g, "algorithm",  meta.algorithm);
    if (st >= 0) st = writeStrAttr(g, "solver",     meta.solver);
    if (st >= 0) st = writeStrAttr(g, "units",      meta.units);
    if (st >= 0) st = writeStrAttr(g, "timestamp",  meta.timestamp);

    // RUN_ATTRS_INT (profiler_schema.py:65, :92-93) — i8.
    if (st >= 0) st = writeI64Attr(g, "nSteps",  meta.nSteps);
    if (st >= 0) st = writeI64Attr(g, "nDOF",    meta.nDOF);
    if (st >= 0) st = writeI64Attr(g, "nElem",   meta.nElem);
    if (st >= 0) st = writeI64Attr(g, "nNode",   meta.nNode);
    if (st >= 0) st = writeI64Attr(g, "nnz",     meta.nnz);
    if (st >= 0) st = writeI64Attr(g, "threads", meta.threads);

    // RUN_ATTRS_FLOAT (profiler_schema.py:66, :94-95) — f8.
    if (st >= 0) st = writeF64Attr(g, "dt_min",           meta.dt_min);
    if (st >= 0) st = writeF64Attr(g, "dt_max",           meta.dt_max);
    if (st >= 0) st = writeF64Attr(g, "dt_cr",            meta.dt_cr);
    if (st >= 0) st = writeF64Attr(g, "oversample_ratio", meta.oversample_ratio);
    if (st >= 0) st = writeF64Attr(g, "wall_ms_total",    meta.wall_ms_total);
    if (st >= 0) st = writeF64Attr(g, "cpu_ms_total",     meta.cpu_ms_total);

    if (st < 0) { H5Gclose(g); H5Gclose(runs); return false; }

    // /runs/<id>/rollup/<root>/... (profiler_schema.py:97). The rollup group
    // wraps the named root node, which writeNode emits as a subgroup.
    bool ok = true;
    {
        hid_t rollupG = H5Gcreate(g, "rollup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (rollupG < 0) { ok = false; }
        else {
            ok = writeNode(rollupG, rollup);
            H5Gclose(rollupG);
        }
    }

    // /series (optional): skip when null or empty (profiler_schema.py:99 `if series:`).
    if (ok && series != nullptr && series->nSteps() > 0)
        ok = writeSeries(g, *series);

    // /memory (always).
    if (ok)
        ok = writeMemory(g, memory);

    H5Gclose(g);
    H5Gclose(runs);

    if (!ok) return false;

    // single global flush so the run is durable across multi-run accumulation.
    if (H5Fflush(fid, H5F_SCOPE_GLOBAL) < 0) return false;
    return true;
}

} // namespace ops_profiler
