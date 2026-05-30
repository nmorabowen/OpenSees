"""Shared params table (D1) — ONE file read by both pytest and Tcl, so the two
halves of a feature's test cannot drift.

Format: tests/params/<feature>.params, flat `key value` lines, `#` comments.
Values are coerced to float when they parse as a number, else kept as str.

    from _testbed.params import load_params
    P = load_params("forceBeamColumn_cantilever")
    L, E = P["L"], P["E"]

The matching 6-line Tcl reader lives in tests/tcl/_params.tcl (`read_params`).
Both point at the same file → single source of truth, no generator step.
"""
import pathlib

_PARAMS_DIR = pathlib.Path(__file__).resolve().parent.parent / "params"


def _coerce(tok: str):
    try:
        return int(tok)
    except ValueError:
        pass
    try:
        return float(tok)
    except ValueError:
        return tok


def load_params(name: str) -> dict:
    path = _PARAMS_DIR / f"{name}.params"
    if not path.exists():
        raise FileNotFoundError(f"no shared params table: {path}")
    out: dict = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        key, _, rest = line.partition(" ")
        out[key] = _coerce(rest.strip())
    return out
