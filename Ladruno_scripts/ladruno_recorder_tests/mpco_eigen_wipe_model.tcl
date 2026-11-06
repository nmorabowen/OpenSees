# MPCO modes recorder after `eigen; wipe` (apeGmsh #1983), classic Tcl exe.
#
# MPCORecorder used to gate its modal preamble on the interpreter's numEigen,
# which `wipe` never resets, and then called Domain::getEigenvalues(), which
# exit(-1)s when the domain has no spectrum. So `eigen` -> `wipe` -> new model
# -> `recorder mpco ... -N modesOfVibration` -> `analyze` killed the process
# before the deck's last line. The fix gates on Domain::getNumEigenvalues()
# (the port of LadrunoRecorder's WP-163 R1 gate).
#
# Phase A (no wipe between eigen and the recorder) is the unchanged path:
#   eig_nowipe.mpco must hold MODE_0 with LAMBDA = 3EI/(m L^3).
# Phase B is the repro: eigen, wipe, new model, MPCO modes recorder, analyze.
#   eig_wipe.mpco must hold no modal step, and the deck must reach its marker.
#
# Usage:  OpenSees.exe mpco_eigen_wipe_model.tcl <out_dir>
# Prints MPCO_EIGEN_WIPE_MODEL OK on success. OpenSees.exe exits 0 even when
# a Tcl error ends the script, so the marker is the proof the deck completed.
# Checked by mpco_eigen_wipe_check.py, which runs this deck.

if {$argc < 1} { puts "usage: mpco_eigen_wipe_model.tcl <out_dir>"; exit 1 }
set out_dir [lindex $argv 0]
set out_nowipe [file join $out_dir eig_nowipe.mpco]
set out_wipe [file join $out_dir eig_wipe.mpco]
foreach p [list $out_nowipe $out_wipe] {
    if {[file exists $p]} { file delete $p }
}

# vertical cantilever, lumped translational mass at the tip, no rotational
# mass: the lateral mode condenses the rotation, k = 3EI/L^3 = 6.0e7
proc build_model {} {
    wipe
    model basic -ndm 2 -ndf 3
    node 1 0.0 0.0
    node 2 0.0 1.0
    fix 1 1 1 1
    mass 2 1.0 1.0 0.0
    geomTransf Linear 1
    element elasticBeamColumn 1 1 2 0.01 200.0e9 1.0e-4 1
}

proc static_step {} {
    # `eigen` leaves a transient analysis behind; drop it (not the spectrum)
    wipeAnalysis
    timeSeries Linear 1
    pattern Plain 1 1 { load 2 1.0 0.0 0.0 }
    constraints Plain
    numberer Plain
    system BandGeneral
    test NormDispIncr 1.0e-9 10
    algorithm Linear
    integrator LoadControl 1.0
    analysis Static
    analyze 1
}

# ---- Phase A: eigen, then the MPCO modes recorder, no wipe in between ----
build_model
recorder mpco $out_nowipe -N modesOfVibration
set lam [eigen 1]
puts "PHASE_A_EIGEN $lam"
static_step
wipe
puts "PHASE_A OK"

# ---- Phase B: the #1983 repro ---------------------------------------------
build_model
eigen 1
build_model
recorder mpco $out_wipe -N modesOfVibration
static_step
wipe
puts "PHASE_B OK"

foreach p [list $out_nowipe $out_wipe] {
    if {![file exists $p]} { puts "FAIL: $p not written"; exit 1 }
}
puts "MPCO_EIGEN_WIPE_MODEL OK"
