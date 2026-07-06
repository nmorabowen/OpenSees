"""ZFP re-encode benchmark for the Ladruno recorder (ADR 08, "Benchmark (do
this first)"). Zero recorder code, zero build risk: reads the f64 DATA arrays
out of existing .ladruno files and re-encodes them under candidate filter
configs, measuring size / error / speed against the f64 oracle.

Design points (post adversarial review, 2026-07-06):

* CHUNKING IS RECORDER-TRUE, NOT SOURCE-COPIED. The recorder sizes chunks by
  the ON-DISK element size (Ladruno_Hdf5.h chunkSlabsPerBlock: ~256 KiB budget,
  so an f32 chunk holds 2x the time-slabs of an f64 chunk). Re-using the f64
  source chunk shape for the f32 re-encode would under-measure the real f32
  baseline — the exact number the ADR gate divides by. Every 3-D config gets
  chunkSlabsPerBlock(dtype) chunks.

* FULL 2x2 CODEC x LAYOUT MATRIX. The component-separated / time-major
  prototype exists for gzip TOO (f32_gzip4_tm / f64_gzip4_tm), so a config-5
  win can be attributed to layout vs codec — otherwise the FORMAT_VERSION-v2
  question is unanswerable.

* THE GATE IS PER-FAMILY, NOT AGGREGATE. One huge smooth nodal field must not
  mask poor element-history compression: the marginal-vs-f32 ratio is reported
  per source dataset, and the verdict line quotes the WORST family.

* TIMING CAVEAT: enc/dec are end-to-end Python wall-clock (h5py + astype +
  transpose overhead included). They bound end-to-end behavior; they are NOT
  codec-throughput measurements and cannot confirm the ADR's literature
  "20-60x faster than gzip" claim.

Error metrics per lossy config: worst ABSOLUTE error vs the f64 oracle, worst
RANGE-RELATIVE error (abs err / per-dataset max|value|) — the headline number,
comparable across codecs whose bounds are absolute (ZFP accuracy) vs relative
(f32) — and masked pointwise relative error (only over |oracle| > 1e-6*range;
the absolute-error trap on near-zeros would otherwise dominate).

Run with the venv python (h5py + hdf5plugin):
  python reencode_bench.py <file.ladruno> [more.ladruno ...] [--out results.json]
"""
from __future__ import annotations

import json
import os
import sys
import time

import h5py
import numpy as np

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import hdf5plugin  # noqa: E402  (registers the filters on import)

ZFP_TOLS = (1e-4, 1e-6, 1e-8, 1e-10)
SZ3_TOLS = (1e-6,)
TM_TOLS = (1e-6, 1e-8)          # tolerances for the layout prototype (config 5)
BASELINE = "f32_gzip4"          # the thing to beat (ADR decision rule)


def selfcheck_filters():
    """Hard-fail if ZFP/SZ3 aren't actually usable — otherwise every lossy
    config silently becomes an 'error' row and the gate has no data."""
    for label, filt in (("ZFP", hdf5plugin.Zfp(accuracy=1e-6)),
                        ("SZ3", hdf5plugin.SZ3(absolute=1e-6))):
        with h5py.File("selfcheck.h5", "w", driver="core",
                       backing_store=False) as f:
            data = np.linspace(0.0, 1.0, 4096).reshape(16, 16, 16)
            d = f.create_dataset("x", data=data, chunks=(16, 16, 16), **dict(filt))
            back = d[...]
            if not np.allclose(back, data, atol=1e-5):
                raise RuntimeError(f"{label} round-trip mismatch")
    print("filter self-check OK (ZFP, SZ3)")


def recorder_chunks(shape, elem_bytes):
    """Mirror Ladruno_Hdf5.h chunkSlabsPerBlock/createTimeSeries3d: ~256 KiB
    chunk budget in ON-DISK bytes -> f32 chunks span 2x the time-slabs of f64.
    Clamped to the (fixed-size) re-encode dataset extent."""
    T, nid_, ncomp = shape
    per = max(1, nid_ * ncomp) * elem_bytes
    n = 262144 // per
    n = max(1, min(1024, n))
    return (min(n, max(T, 1)), max(nid_, 1), max(ncomp, 1))


def collect_data_datasets(path: str):
    """(name, ndarray f64) for every chunked 3-D [T x nIds x nComp] RESULTS
    DATA dataset in a .ladruno file."""
    out = []
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if (isinstance(obj, h5py.Dataset) and name.endswith("/DATA")
                    and obj.ndim == 3 and obj.chunks is not None):
                out.append((name, obj[...].astype(np.float64)))
        f.visititems(visit)
    return out


def encode_config(datasets, out_path, cfg):
    """Write all datasets under one config; return result row dict."""
    if os.path.exists(out_path):
        os.remove(out_path)
    dtype = cfg.get("dtype", np.float64)
    elem = 4 if dtype == np.float32 else 8

    def filt_kwargs(arr=None):
        if cfg["codec"] == "gzip":
            return dict(shuffle=True, compression="gzip", compression_opts=4)
        if cfg["codec"] == "zfp_rr" and arr is not None:
            # per-dataset RANGE-SCALED absolute tolerance — what a real
            # `-compression zfp` with per-quantity tolerance plumbing would
            # do (ADR work item (c)); rr=5.5e-8 matches f32's fidelity
            tol = max(cfg["rr"] * float(np.abs(arr).max()), 1e-300)
            return dict(hdf5plugin.Zfp(accuracy=tol))
        return dict(cfg["filt"])

    t0 = time.perf_counter()
    with h5py.File(out_path, "w") as fh:
        for name, arr in datasets:
            key = name.replace("/", "|")
            T, nid_, ncomp = arr.shape
            if cfg["layout"] == "3d":
                ck = recorder_chunks(arr.shape, elem)
                fh.create_dataset(key, data=arr.astype(dtype), dtype=dtype,
                                  chunks=ck, **filt_kwargs(arr))
            else:  # 'tm' — component-separated, time on the fast axis
                for c in range(ncomp):
                    comp = np.ascontiguousarray(arr[:, :, c].T)  # [nIds x T]
                    ck = (min(64, nid_), min(1024, T))
                    fh.create_dataset(f"{key}|c{c}", data=comp.astype(dtype),
                                      dtype=dtype, chunks=ck, **filt_kwargs(comp))
    enc = time.perf_counter() - t0

    per_ds = {}
    with h5py.File(out_path, "r") as fh:
        for name, arr in datasets:
            key = name.replace("/", "|")
            if cfg["layout"] == "tm":
                stored = sum(fh[f"{key}|c{c}"].id.get_storage_size()
                             for c in range(arr.shape[2]))
            else:
                stored = fh[key].id.get_storage_size()
            per_ds[name] = stored
    size = sum(per_ds.values())
    fsize = os.path.getsize(out_path)  # incl. metadata — the honest file cost

    # decode + error vs oracle
    t0 = time.perf_counter()
    max_abs = max_rangerel = max_rel = 0.0
    with h5py.File(out_path, "r") as fh:
        for name, arr in datasets:
            key = name.replace("/", "|")
            if cfg["layout"] == "tm":
                back = np.empty_like(arr)
                for c in range(arr.shape[2]):
                    back[:, :, c] = fh[f"{key}|c{c}"][...].astype(np.float64).T
            else:
                back = fh[key][...].astype(np.float64)
            d = np.abs(back - arr)
            rng = np.abs(arr).max()
            if rng == 0.0:
                continue
            max_abs = max(max_abs, float(d.max()))
            max_rangerel = max(max_rangerel, float(d.max() / rng))
            mask = np.abs(arr) > 1e-6 * rng
            if mask.any():
                max_rel = max(max_rel, float((d[mask] / np.abs(arr[mask])).max()))
    dec = time.perf_counter() - t0
    return {"bytes": size, "file_bytes": fsize, "enc_s": enc, "dec_s": dec,
            "max_abs": max_abs, "max_range_rel": max_rangerel,
            "max_rel": max_rel, "per_dataset": per_ds}


def build_configs():
    cfgs = [
        ("f64_gzip4", {"layout": "3d", "codec": "gzip", "dtype": np.float64}),
        (BASELINE,    {"layout": "3d", "codec": "gzip", "dtype": np.float32}),
    ]
    for tol in ZFP_TOLS:
        cfgs.append((f"zfp_acc_{tol:g}",
                     {"layout": "3d", "codec": "filt",
                      "filt": hdf5plugin.Zfp(accuracy=tol)}))
    cfgs.append(("zfp_reversible",
                 {"layout": "3d", "codec": "filt",
                  "filt": hdf5plugin.Zfp(reversible=True)}))
    for tol in SZ3_TOLS:
        cfgs.append((f"sz3_abs_{tol:g}",
                     {"layout": "3d", "codec": "filt",
                      "filt": hdf5plugin.SZ3(absolute=tol)}))
    # per-dataset range-scaled tolerance: the fair "what -compression zfp
    # would really do" config; 5.5e-8*range == f32-matched fidelity
    for rr in (5.5e-8, 1e-6):
        cfgs.append((f"zfp_rr_{rr:g}",
                     {"layout": "3d", "codec": "zfp_rr", "rr": rr}))
        cfgs.append((f"zfp_tm_rr_{rr:g}",
                     {"layout": "tm", "codec": "zfp_rr", "rr": rr}))
    # layout prototype — full codec x layout matrix so layout-vs-codec is
    # attributable (BLOCKER #2 of the methodology review)
    for tol in TM_TOLS:
        cfgs.append((f"zfp_tm_acc_{tol:g}",
                     {"layout": "tm", "codec": "filt",
                      "filt": hdf5plugin.Zfp(accuracy=tol)}))
    cfgs.append(("f32_gzip4_tm", {"layout": "tm", "codec": "gzip",
                                  "dtype": np.float32}))
    cfgs.append(("f64_gzip4_tm", {"layout": "tm", "codec": "gzip",
                                  "dtype": np.float64}))
    return cfgs


def bench_file(path: str):
    datasets = collect_data_datasets(path)
    raw = sum(a.nbytes for _, a in datasets)
    print(f"\n=== {os.path.basename(path)}: {len(datasets)} DATA datasets, "
          f"raw f64 {raw/1e6:.1f} MB ===")

    scratch = os.path.join(os.path.dirname(os.path.abspath(path)), "_reencode")
    os.makedirs(scratch, exist_ok=True)
    rows = {}
    for cname, cfg in build_configs():
        out = os.path.join(scratch, f"{os.path.basename(path)}.{cname}.h5")
        try:
            r = encode_config(datasets, out, cfg)
        except Exception as e:  # a filter refusing one dataset is a finding
            print(f"  {cname:<18} FAILED: {e}")
            rows[cname] = {"error": str(e)}
            continue
        rows[cname] = r
        print(f"  {cname:<18} {r['bytes']/1e6:9.2f} MB  {r['bytes']/raw:6.3f}x raw"
              f"  rangerel {r['max_range_rel']:.2e}  rel {r['max_rel']:.2e}"
              f"  abs {r['max_abs']:.2e}"
              f"  enc {r['enc_s']:6.2f}s dec {r['dec_s']:6.2f}s")

    base = rows.get(BASELINE)
    if base and "bytes" in base:
        print(f"  -- marginal vs {BASELINE} (the ADR gate; worst family in [];"
              " ratios >1 = smaller than baseline) --")
        for cname, r in rows.items():
            if "bytes" not in r or cname == BASELINE:
                continue
            agg = base["bytes"] / r["bytes"]
            worst_name, worst = None, float("inf")
            for dsname, b in r["per_dataset"].items():
                bb = base["per_dataset"].get(dsname)
                if bb and b > 0:
                    ratio = bb / b
                    if ratio < worst:
                        worst_name, worst = dsname, ratio
            fam = (worst_name or "?").split("RESULTS/")[-1].split("/DATA")[0]
            print(f"  {cname:<18} agg {agg:5.2f}x | worst-family {worst:5.2f}x"
                  f"  [{fam}]")
            r["marginal_agg"] = agg
            r["marginal_worst"] = worst
            r["marginal_worst_dataset"] = worst_name
    print("  NOTE: enc/dec are end-to-end Python wall-clock (astype/transpose"
          " included) — NOT codec throughput.")
    return {"file": path, "raw_bytes": raw, "n_datasets": len(datasets),
            "configs": rows}


def main():
    selfcheck_filters()
    out_json = None
    argv = sys.argv[1:]
    if "--out" in argv:
        i = argv.index("--out")
        out_json = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    args = [a for a in argv if not a.startswith("--")]
    results = [bench_file(p) for p in args]
    if out_json:
        with open(out_json, "w") as f:
            json.dump(results, f, indent=1)
        print(f"\nwrote {out_json}")


if __name__ == "__main__":
    main()
