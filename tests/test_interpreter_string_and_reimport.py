"""Interpreter fixes: string arguments read into a caller buffer under Tcl.

  tcl_string_from_all   OPS_GetStringFromAll in the classic Tcl interpreter
                        fills the buffer the caller passes in

KinematicCoupling reads its -k value and its -dof components through
OPS_GetStringFromAll into a local buffer. Reference values, 2D model with
the reference node fixed and two free slaves tied in x and y:

  -k 1000.0        eleResponse kt = 1000.0; a unit load in x on a slave,
                   whose only stiffness is the tie, gives u = P/Kt = 1.0e-3
  -dof 1 2         slaves with 3 DOFs, rotation left free: 2 x 2 = 4 tied
                   DOFs (6 with the default component list)

Before the fix the buffer is never written under Tcl, so the element parser
reads uninitialised memory and rejects the option; the element command fails.

The Tcl case needs the OpenSees executable named by the OPENSEES_EXE
environment variable and is skipped when it is unset or missing.
"""
import os
import subprocess

import pytest


TCL_DECK = r'''
set nfail 0
proc check {name ok detail} {
    global nfail
    if {$ok} {
        puts "PASS $name"
    } else {
        incr nfail
        puts "FAIL $name -- $detail"
    }
}

# -k value read through OPS_GetStringFromAll
model basic -ndm 2 -ndf 2
node 1 0.0 0.0 -ndf 3
node 2 1.0 0.0
node 3 0.0 1.0
fix 1 1 1 1
set rc [catch {element KinematicCoupling 1 1 2 2 3 -k 1000.0} msg]
check k_element [expr {$rc == 0}] $msg
if {$rc == 0} {
    set kt [eleResponse 1 kt]
    check k_value [expr {abs($kt - 1000.0) <= 1.0e-9}] "kt = $kt"
    timeSeries Linear 1
    pattern Plain 1 1 { load 2 1.0 0.0 }
    constraints Plain
    numberer Plain
    system FullGeneral
    test NormDispIncr 1.0e-12 10 0
    algorithm Newton
    integrator LoadControl 1.0
    analysis Static
    set arc [analyze 1]
    set ux [nodeDisp 2 1]
    check k_disp [expr {$arc == 0 && abs($ux - 1.0e-3) <= 1.0e-12}] "analyze $arc, ux = $ux"
} else {
    check k_value 0 "no element"
    check k_disp 0 "no element"
}

# -dof components read through OPS_GetStringFromAll
wipe
model basic -ndm 2 -ndf 3
node 1 0.0 0.0
node 2 1.0 0.0
node 3 0.0 1.0
fix 1 1 1 1
set rc [catch {element KinematicCoupling 1 1 2 2 3 -dof 1 2 -k 1000.0} msg]
check dof_element [expr {$rc == 0}] $msg
if {$rc == 0} {
    set nt [eleResponse 1 tiedDOFs]
    check dof_count [expr {abs($nt - 4.0) <= 1.0e-12}] "tiedDOFs = $nt"
} else {
    check dof_count 0 "no element"
}

puts "SELF-TEST: $nfail failure(s)"
'''


def _tcl_exe():
    exe = os.environ.get("OPENSEES_EXE")
    return exe if exe and os.path.isfile(exe) else None


@pytest.mark.skipif(_tcl_exe() is None,
                    reason="OPENSEES_EXE is unset or does not name a file")
def test_tcl_string_from_all_fills_buffer(tmp_path):
    deck = tmp_path / "string_from_all.tcl"
    deck.write_text(TCL_DECK)
    proc = subprocess.run(
        [_tcl_exe(), str(deck)], cwd=str(tmp_path),
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    out = proc.stdout + proc.stderr
    assert "SELF-TEST:" in out, out
    assert out.count("PASS ") + out.count("FAIL ") == 5, out
    assert "SELF-TEST: 0 failure(s)" in out, out
