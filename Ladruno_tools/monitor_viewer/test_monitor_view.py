"""
test_monitor_view.py  —  binary-free tests for the monitor viewer.

Writes a synthetic sink that mirrors ladruno::MonitorSink's on-disk layout
(no OpenSees build needed), then exercises MonitorReader (incremental reads) and
the FastAPI endpoints via TestClient.  Run:  python test_monitor_view.py
"""
import os
import tempfile

import numpy as np
import h5py

from monitor_reader import MonitorReader


def write_sink(path, columns, nframes, units=""):
    with h5py.File(path, "w", libver="latest") as f:
        f.attrs["FORMAT"] = "ladruno-monitor"
        f.attrs["FORMAT_VERSION"] = 1
        f.attrs["GENERATOR"] = "test"
        if units:
            f.attrs["units"] = units
        vstr = h5py.special_dtype(vlen=str)
        f.create_dataset("COLUMNS", data=np.array(columns, dtype=object), dtype=vstr)
        step = np.arange(1, nframes + 1, dtype="i4")
        t = np.arange(1, nframes + 1, dtype="f8") * 0.01
        frames = np.column_stack([
            np.sin(t * (j + 1)) for j in range(len(columns))
        ]).astype("f8") if columns else np.zeros((nframes, 0))
        f.create_dataset("STEP", data=step, maxshape=(None,), chunks=(256,))
        f.create_dataset("TIME", data=t, maxshape=(None,), chunks=(256,))
        f.create_dataset("FRAMES", data=frames,
                         maxshape=(None, len(columns)),
                         chunks=(256, max(1, len(columns))))


def test_reader_meta_and_incremental():
    p = os.path.join(tempfile.mkdtemp(), "s.h5")
    cols = ["node2.disp.dof1", "node5.disp.dof1"]
    write_sink(p, cols, 300, units="m")
    r = MonitorReader(p)

    m = r.meta()
    assert m["columns"] == cols
    assert m["nframes"] == 300 and m["nColumns"] == 2
    assert m["format"] == "ladruno-monitor" and m["units"] == "m"

    # incremental cursor semantics
    d0 = r.frames_since(0)
    assert d0["n"] == 300 and len(d0["step"]) == 300
    assert len(d0["rows"]) == 300 and len(d0["rows"][0]) == 2
    d1 = r.frames_since(290)
    assert d1["since"] == 290 and len(d1["step"]) == 10
    assert d1["step"][0] == 291                          # 1-based STEP, resumed
    assert r.frames_since(300)["step"] == []             # nothing new
    assert r.frames_since(999)["step"] == []             # cursor past end


def test_reader_rejects_foreign_h5():
    p = os.path.join(tempfile.mkdtemp(), "x.h5")
    with h5py.File(p, "w") as f:
        f.create_dataset("nope", data=[1, 2, 3])
    try:
        MonitorReader(p).meta()
        raise AssertionError("expected ValueError on non-monitor h5")
    except (ValueError, KeyError):
        pass


def test_api_endpoints():
    from fastapi.testclient import TestClient
    from monitor_server import create_app

    p = os.path.join(tempfile.mkdtemp(), "s.h5")
    write_sink(p, ["node2.disp.dof1"], 120)
    client = TestClient(create_app(p))

    h = client.get("/health").json()
    assert h["status"] == "ok" and h["nframes"] == 120

    m = client.get("/api/meta").json()
    assert m["columns"] == ["node2.disp.dof1"] and m["nframes"] == 120

    d = client.get("/api/frames", params={"since": 100}).json()
    assert d["n"] == 120 and len(d["step"]) == 20 and d["step"][0] == 101

    assert "since" in client.get("/api/frames").json()
    assert client.get("/").status_code == 200        # serves the page
    assert client.get("/api/frames", params={"since": -1}).status_code == 422  # ge=0


def test_missing_sink_503():
    from fastapi.testclient import TestClient
    from monitor_server import create_app
    client = TestClient(create_app("does_not_exist.h5"))
    assert client.get("/health").status_code == 503
    assert client.get("/api/meta").status_code == 503


def main():
    test_reader_meta_and_incremental()
    test_reader_rejects_foreign_h5()
    test_api_endpoints()
    test_missing_sink_503()
    print("OK — monitor viewer (reader incremental + api endpoints + 503 + foreign-h5 reject)")


if __name__ == "__main__":
    main()
