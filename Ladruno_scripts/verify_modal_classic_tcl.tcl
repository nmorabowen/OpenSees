# Classic-Tcl smoke for the ADR-44 modal-response family (wired into
# SRC/tcl/commands.cpp). Proves the three commands dispatch from a sourced .tcl
# (OpenSees.exe model.tcl) — before this wiring they hard-errored with
# "invalid command name". Any failure aborts before the final ALL...PASS line.
#
# Also probes whether the frequency-domain commands return their {f,...} table as
# a Tcl list result (OPS_SetDoubleListsOutput -> Tcl_SetObjResult), in addition to
# writing -out.
#
# Run: OpenSees.exe verify_modal_classic_tcl.tcl

# --- 2-DOF axial chain (ndf=1), mirror the pytest battery ---------------------
wipe
model basic -ndm 1 -ndf 1
node 1 0.0;  fix 1 1
node 2 1.0;  mass 2 2.0
node 3 2.0;  mass 3 1.5
uniaxialMaterial Elastic 1 800.0
uniaxialMaterial Elastic 2 500.0
element Truss 1 1 2 1.0 1
element Truss 2 2 3 1.0 2

eigen -fullGenLapack 2
modalProperties
puts "OK eigen + modalProperties"

# --- frequencyResponse: dispatches, writes -out, and returns a list -----------
set frf [frequencyResponse -freq 0.1 4.0 20 -baseAccel -dir 1 -damp 0.04 \
             -node 3 -dof 1 -out _frf_classic.out]
if {![file exists _frf_classic.out]} {
    puts "FAIL frequencyResponse did not write -out file"
    exit 1
}
set nrows [llength $frf]
if {$nrows > 0} {
    set row0 [lindex $frf 0]
    puts "OK frequencyResponse dispatched; -out written; RETURN LIST works ($nrows rows, row0={$row0})"
} else {
    puts "OK frequencyResponse dispatched; -out written; RETURN LIST empty (use -out)"
}

# --- steadyStateDynamics ------------------------------------------------------
set ssd [steadyStateDynamics -freq 0.1 4.0 20 -baseAccel -dir 1 -damp 0.04 \
             -node 3 -dof 1 -out _ssd_classic.out]
if {![file exists _ssd_classic.out]} {
    puts "FAIL steadyStateDynamics did not write -out file"
    exit 1
}
puts "OK steadyStateDynamics dispatched; -out written; return rows=[llength $ssd]"

# --- modalResponseHistory (P1a transient) -------------------------------------
timeSeries Path 1 -dt 0.01 -values {0.0 1.0 0.0 -1.0 0.0} -factor 1.0
modalResponseHistory -dt 0.01 -nsteps 4 -baseAccel 1 -dir 1 -damp 0.05
puts "OK modalResponseHistory dispatched"

puts "ALL MODAL CLASSIC-TCL CHECKS PASS"
