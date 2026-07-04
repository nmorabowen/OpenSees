"""
test_contract.py  —  guards the profiler HDF5 contract. Run: python test_contract.py
Plain asserts (no pytest dependency) so it runs anywhere h5py is present.
"""
import os
import tempfile
import h5py
from profiler_schema import write_run, SCHEMA_VERSION
from profiler_results import ProfilerResults


def _tiny_rollup(ft_ms):
    w = lambda ms: int(ms * 1e6)
    leaf = lambda n, ms, calls, **kw: {
        "name": n, "calls": calls, "wall_ns": w(ms),
        "wall_ns_min": 0, "wall_ns_max": 0, "cpu_ns": w(ms), **kw}
    return {**leaf("step", ft_ms + 10, 5), "children": [
        {**leaf("solveCurrentStep", ft_ms + 10, 5), "children": [
            leaf("formTangent", ft_ms, 15,
                 elem_by_type=[(74, 2, w(ft_ms * 0.9), w(ft_ms * 0.9) / 2, 1),
                               (12, 50, w(ft_ms * 0.1), w(ft_ms * 0.1) / 50, 0)]),
            leaf("linearSolve", 10.0, 15),
        ]},
    ]}


def _untimed_root_rollup(ft_ms):
    # Mirror the REAL C++ engine output: an untimed 'root' anchor (wall_ns=0,
    # calls=0) wrapping the first real phase. The synthetic _tiny_rollup above
    # uses 'step' as the top node; the live engine does not.
    return {"name": "root", "calls": 0, "wall_ns": 0,
            "wall_ns_min": 0, "wall_ns_max": 0, "cpu_ns": 0,
            "children": [_tiny_rollup(ft_ms)]}


def main():
    path = os.path.join(tempfile.mkdtemp(), "t.h5")
    with h5py.File(path, "w") as f:
        write_run(f, "A", meta=dict(integrator="Newmark", algorithm="Newton",
                                    nDOF=100, nSteps=5, wall_ms_total=110.0),
                  rollup=_tiny_rollup(100.0))
        write_run(f, "B", meta=dict(integrator="Newmark", algorithm="ModifiedNewton",
                                    nDOF=100, nSteps=5, wall_ms_total=60.0),
                  rollup=_tiny_rollup(50.0))
        # C mirrors the real engine: an untimed 'root' anchor above 'step'.
        write_run(f, "C", meta=dict(integrator="LoadControl", algorithm="Newton",
                                    nDOF=20, nSteps=5, wall_ms_total=110.0),
                  rollup=_untimed_root_rollup(100.0))

    # 1. immutability: re-writing an existing run id must fail (do this before any reader
    #    holds the file open — Windows HDF5 forbids concurrent write + read-only handles)
    try:
        with h5py.File(path, "a") as f:
            write_run(f, "A", meta=dict(nSteps=5), rollup=_tiny_rollup(1.0))
        raise AssertionError("expected ValueError on duplicate run id")
    except ValueError:
        pass

    pr = ProfilerResults(path)

    # 2. schema + runs present
    assert pr.schema_version == SCHEMA_VERSION
    assert set(pr.runs()) == {"A", "B", "C"}

    # 3. rollup normalization: share of root == 1, per-step = wall/nSteps
    rollA = pr.rollup("A")
    assert abs(rollA["share"] - 1.0) < 1e-9
    assert abs(rollA["wall_ms_per_step"] - rollA["wall_ms"] / 5) < 1e-6

    # 4. elem_by_type carried through with fb_coupled flag
    ft = rollA["children"][0]["children"][0]
    assert ft["name"] == "formTangent"
    tags = {e["classTag"]: e for e in ft["elem_by_type"]}
    assert tags[74]["fb_coupled"] is True and tags[12]["fb_coupled"] is False

    # 5. diff by stable path: formTangent halved -> -50%
    d = {r["path"]: r for r in pr.diff("A", "B")}
    ft_diff = d["step/solveCurrentStep/formTangent"]
    assert abs(ft_diff["base_ms"] - 100.0) < 1e-6
    assert abs(ft_diff["cand_ms"] - 50.0) < 1e-6
    assert abs(ft_diff["d_pct"] - (-50.0)) < 1e-6
    assert ft_diff["status"] == "changed"

    # 6. untimed 'root' anchor (real engine shape): the backend must backfill the
    #    root's aggregate wall from its children so `share` is a real fraction and
    #    the flame graph has a non-zero parent to lay out. Regression for the
    #    `wall_ns or 1` (== 1 ns denominator) bug that made every real-run share
    #    astronomically large and left the icicle empty.
    rollC = pr.rollup("C")
    assert rollC["name"] == "root"
    stepC = rollC["children"][0]
    assert rollC["wall_ms"] > 0, "untimed root wall must be backfilled from children"
    assert abs(rollC["wall_ms"] - stepC["wall_ms"]) < 1e-6, "root == sum(children) == step"
    assert abs(rollC["share"] - 1.0) < 1e-9, "root share must be 1.0"
    # every descendant share is a sane fraction, NOT the 1-ns-denominator blowup
    def _all_shares_sane(n):
        assert 0.0 <= n["share"] <= 1.0 + 1e-9, f"{n['name']} share out of range: {n['share']}"
        for c in n["children"]:
            _all_shares_sane(c)
    _all_shares_sane(rollC)
    assert abs(rollC["wall_ms_per_step"] - rollC["wall_ms"] / 5) < 1e-6

    pr.close()
    print("OK — contract holds (schema, immutability, normalizers, elem_by_type, "
          "diff, untimed-root backfill)")


if __name__ == "__main__":
    main()
