# ADR-76 R4 acceptance test — `algorithm ModifiedNewton` option parsing

Verifies acceptance test 4 of [[76_ladruno_tangent_reuse_adr]]: that
`-initial` and `-factoronce` **both** take effect, in either argument order and
in every accepted spelling.

## Running it

Two steps, two different interpreters. This trips people up, so it is spelled out.

**1. Produce the profile** — needs a **fork** `OpenSees.exe` (the `profiler`
command is Ladruno-only; a vanilla build will fail at `profiler start`):

```bash
OpenSees.exe adr76_factoronce_model.tcl
```

Writes `adr76_factoronce.h5` and `adr76_factoronce_disp.txt` next to the deck.
Both are gitignored — see `.gitignore` for why that matters.

**2. Check it** — needs `h5py`, which is **not** installed in the build's Python
and is not a fork build dependency. Use a venv:

```bash
python3.12 -m venv .venv
.venv/Scripts/python -m pip install -r ../../Ladruno_tools/profiler_viewer/requirements.txt
.venv/Scripts/python adr76_factoronce_check.py .
```

Exit status is 0 on success, 1 on any failed assertion.

## Reading the result

Expect **11 PASS lines** and the terminal marker:

```
ADR-76 R4 smoke: all checks passed
```

> **Gate on that marker, not on the exit code.** `OpenSees.exe` exits **0** after
> a Tcl *parse* error — the deck aborts, prints nothing, and the process still
> reports success. That trap is recorded in `LEDGER_quirks.md`; it bit this very
> deck during development. Any CI wiring must grep for the marker line.

## What it does and does not prove

**Does:** that the option loop honours every recognised token regardless of
position or spelling. Pre-fix, `-initial -factoronce` collapsed onto `-initial`
and `-factoronce -initial` onto `-factoronce`; cases B/C/E/F all fail on the old
parser.

**Does not:** anything about tangent staleness. The profiled window is
deliberately **linear** (all ten bars sit in the same `Steel01` hardening branch
under monotonic load), which is what makes the assertions integer-crisp — and
also means this deck cannot detect a stale latched tangent, a missed
`domainChanged` re-arm, or any `factorOnce` invalidation bug. Do not extend it to
cover those; write a separate deck with a real domain change.
