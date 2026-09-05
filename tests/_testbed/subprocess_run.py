"""Shared child-process runner for tests that shell out to a fresh
interpreter -- the pattern used whenever a defect's failure mode is a
process exit (a vanilla "subclass responsibility" `exit(-1)`, an
`MPI_Abort`, ...) that no in-process `assert` could ever observe, because the
interpreter carrying the assertion would be gone with it.

Copied from (and should stay consistent with) `test_soe_zero_free_equations.
py::_run_child`, the file that paid for both Windows traps baked in here; see
`test_adr85_contact2d_t0_refusals.py`'s module docstring for the fuller
writeup of why each flag is load-bearing:

  * `stdin=subprocess.DEVNULL` -- with the inherited stdin, `Popen` tries to
    `DuplicateHandle` whatever pytest's fd-level capture left there and
    raises `OSError: [WinError 50] The request is not supported`,
    intermittently, per shell session and only once some OTHER module in the
    full suite has already run -- so the file passes standalone and fails
    under `pytest tests/`, a failure indistinguishable from the crash under
    test if this is left unguarded.
  * `encoding="utf-8", errors="replace"` -- a banner or `opserr` line with a
    non-ASCII byte otherwise raises `UnicodeDecodeError` inside the capture
    on a cp1252 console, which again looks exactly like the child crashing.
"""
import subprocess
import sys


def run_python_script(script, cwd=None, timeout=300, argv=(), merge_stderr=False):
    """Run `script` as `python -u -c script *argv` in a child process.

    Returns `(returncode, combined_output)`.

    `-u` (unbuffered) is always passed: Python's stdout is FULLY buffered
    (not line-buffered) once it is a pipe rather than a console, while C++'s
    `opserr` is unbuffered, so without `-u` the relative order between a
    script's own `print()` markers and any `opserr` output is an artifact of
    buffer-fill timing, not of when either one actually happened -- eagerly
    flushed `opserr` text can appear to precede `print()` calls that
    executed first. `-u` costs nothing for a caller that only checks whether
    some text is present ANYWHERE in the output.

    `merge_stderr=False` (default) captures stdout and stderr as two
    SEPARATE streams and returns `stdout + stderr` concatenated -- matching
    the long-standing `capture_output=True` convention elsewhere in this
    suite. That concatenation does NOT preserve chronological order BETWEEN
    the two streams (all of stdout precedes all of stderr, regardless of
    when each was written), which is fine for a caller that only checks text
    presence but wrong for one that needs to know whether some stderr text
    appeared before or after a particular stdout marker.

    `merge_stderr=True` redirects stderr into the SAME OS pipe as stdout
    (`stderr=subprocess.STDOUT`), so the OS interleaves the two in the true
    order they were written (combined with `-u`, this is also the true
    chronological order the child produced them in) -- use this whenever a
    test's assertion depends on ordering, e.g. "did the warning appear
    before or after this restore completed".
    """
    kwargs = dict(
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        text=True,
        timeout=timeout,
        encoding="utf-8",
        errors="replace",
    )
    if cwd is not None:
        kwargs["cwd"] = cwd
    if merge_stderr:
        kwargs["stderr"] = subprocess.STDOUT
    else:
        kwargs["stderr"] = subprocess.PIPE

    p = subprocess.run([sys.executable, "-u", "-c", script, *argv], **kwargs)

    if merge_stderr:
        return p.returncode, p.stdout
    return p.returncode, p.stdout + (p.stderr or "")
