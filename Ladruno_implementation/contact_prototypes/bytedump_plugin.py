"""pytest plugin: after every ops.analyze()/ops.analyze_augmented-like call, record repr() of every
node's displacement (and the load factor/time) per test id; write JSON at session end
(BYTEDUMP_OUT). Used to compare two binaries bit for bit."""
import json, os, hashlib
import pytest
import opensees as _o

_rec = {}
_cur = [None]
_orig = _o.analyze


def _snap():
    h = hashlib.sha256()
    try:
        tags = _o.getNodeTags()
    except Exception:
        tags = []
    for t in tags:
        try:
            h.update(repr((t, _o.nodeDisp(t))).encode())
        except Exception:
            pass
    try:
        h.update(repr(_o.getTime()).encode())
    except Exception:
        pass
    return h.hexdigest()[:16]


def _analyze(*a, **k):
    rc = _orig(*a, **k)
    if _cur[0] is not None:
        _rec.setdefault(_cur[0], []).append((rc, _snap()))
    return rc


_o.analyze = _analyze


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    _cur[0] = item.nodeid
    yield
    _cur[0] = None


def pytest_sessionfinish(session, exitstatus):
    out = os.environ.get("BYTEDUMP_OUT")
    if out:
        json.dump(_rec, open(out, "w"), indent=0)
