# ADR-75b G-L3 (FORK VARIANT) — LadrunoBrick + LadrunoJ2, the fork's own
# element and material, which the Tcl nDMaterial parity fix made reachable.
# Removes the 'lower bound' caveat the stdBrick/J2Plasticity proxy carried.
#
# ADR-75b §12 parked Lane 3 pending ONE measurement: on a production-scale deck
# (>=500k DOF), what fraction of step does each INDIVIDUAL assembly loop hold?
#   loop A = Domain::update   (elem.update)
#   loop B = formTangent      (element tangent -> addA)
#   loop C = formUnbalance    (element residual -> addB)
# Gate: largest SINGLE loop >=40% => re-authorize Lane 3 at that loop.
#       no loop reaches 40%      => close Lane 3 as Amdahl-irrelevant.
# NOTE the gate is on a SINGLE LOOP, not the aggregate element-kernel fraction
# (Lane D reads kernel 85.30% but gate FAIL at loop A 38.95%).
#
# COARSE instrumentation on purpose (`profiler start -perStep`, NOT -deep):
# enabled() gates coarse phase timing while deep() additionally gates the
# per-element buckets (Profiler.h:327-334), so -perStep yields the per-loop
# fractions untaxed. It CANNOT give the kernel-vs-scatter split inside a loop
# (elemByType_ is deep-mode only) -- that is a separate, later question.
#
# MATERIAL CAVEAT, stated up front: the Tcl nDMaterial command is a hand-written
# strcmp ladder that does NOT carry the fork's LadrunoJ2 (Python-only), so this
# deck uses vanilla stdBrick + J2Plasticity. Both are CHEAPER per element than
# LadrunoBrick + LadrunoJ2 (which runs 16 function-scope statics per update), so
# the element fractions measured here are a LOWER BOUND on the fork's own decks.
# A PASS is therefore decisive; a FAIL is suggestive, not conclusive.

set pid [getPID]
set np  [getNP]

set NX 44 ; set NY 44 ; set NZ 88
if {[info exists env(GL3_NX)]} { set NX $env(GL3_NX) }
if {[info exists env(GL3_NY)]} { set NY $env(GL3_NY) }
if {[info exists env(GL3_NZ)]} { set NZ $env(GL3_NZ) }
set NSTEP 2
if {[info exists env(GL3_NSTEP)]} { set NSTEP $env(GL3_NSTEP) }
set LOADF 0.02
if {[info exists env(GL3_LOADF)]} { set LOADF $env(GL3_LOADF) }
set PROF 1
if {[info exists env(GL3_PROF)]} { set PROF $env(GL3_PROF) }

set H 100.0
set E 200.0e3
set NU 0.3
set K [expr {$E/(3.0*(1.0-2.0*$NU))}]
set G [expr {$E/(2.0*(1.0+$NU))}]

set NXP [expr {$NX + 1}]
set NYP [expr {$NY + 1}]
proc nid {i j k} { global NXP NYP; return [expr {1 + $i + $NXP*$j + $NXP*$NYP*$k}] }

wipe
model basic -ndm 3 -ndf 3
# J2Plasticity tag K G sig0 sigInf delta H   -- real plasticity, so elem.update
# does genuine constitutive work rather than a trivial elastic evaluation.
# Ladruno fork material -- reachable from Tcl only since the nDMaterial parity
# fix; before it this line failed with 'nDMaterial LadrunoJ2 not defined'.
nDMaterial LadrunoJ2 1 $K $G -iso voce 250.0 0.0 0.0 2000.0

set k0 [expr {$pid * $NZ / $np}]
set k1 [expr {($pid + 1) * $NZ / $np}]
for {set k $k0} {$k <= $k1} {incr k} {
    for {set j 0} {$j < $NYP} {incr j} {
        for {set i 0} {$i < $NXP} {incr i} {
            node [nid $i $j $k] [expr {$i*$H}] [expr {$j*$H}] [expr {$k*$H}]
        }
    }
}
for {set k $k0} {$k < $k1} {incr k} {
    for {set j 0} {$j < $NY} {incr j} {
        for {set i 0} {$i < $NX} {incr i} {
            set e [expr {1 + $i + $NX*$j + $NX*$NY*$k}]
            element LadrunoBrick $e \
                [nid $i $j $k] [nid [expr {$i+1}] $j $k] \
                [nid [expr {$i+1}] [expr {$j+1}] $k] [nid $i [expr {$j+1}] $k] \
                [nid $i $j [expr {$k+1}]] [nid [expr {$i+1}] $j [expr {$k+1}]] \
                [nid [expr {$i+1}] [expr {$j+1}] [expr {$k+1}]] [nid $i [expr {$j+1}] [expr {$k+1}]] \
                1
        }
    }
}
if {$k0 == 0} {
    for {set j 0} {$j < $NYP} {incr j} {
        for {set i 0} {$i < $NXP} {incr i} { fix [nid $i $j 0] 1 1 1 }
    }
}
set tip [nid [expr {$NX/2}] [expr {$NY/2}] $NZ]
if {$k1 == $NZ} {
    timeSeries Linear 1
    pattern Plain 1 1 {
        for {set j 0} {$j < $NYP} {incr j} {
            for {set i 0} {$i < $NXP} {incr i} {
                load [nid $i $j $NZ] 4.0e5 0.0 -4.0e5
            }
        }
    }
}

constraints Plain
numberer LadrunoParallelRCM
system Mumps -ICNTL14 200 -matrixType 2 -stats
test NormDispIncr 1.0e-6 30 0
algorithm Newton
integrator LoadControl $LOADF
analysis Static

set NDOF [expr {3*$NXP*$NYP*($NZ+1)}]
if {$pid == 0} { puts "GL3F_MODEL nx=$NX ny=$NY nz=$NZ dof~$NDOF np=$np loadf=$LOADF nstep=$NSTEP prof=$PROF" }

if {$PROF} { profiler start -perStep }
set t0 [clock milliseconds]
set ok [analyze $NSTEP]
set t1 [clock milliseconds]
if {$PROF} { profiler stop }

if {$pid == 0} { puts "GL3F_RESULT ok=$ok wall_ms=[expr {$t1-$t0}] dof~$NDOF np=$np" }
if {$k1 == $NZ} { puts [format "GL3F_TIP pid=%d ux=%.17e" $pid [nodeDisp $tip 1]] }

if {$PROF} {
    set out "gl3_prof"
    if {[info exists env(GL3_OUT)]} { set out $env(GL3_OUT) }
    profiler report "${out}.rank${pid}.h5"
    if {$pid == 0} { puts "GL3F_REPORT ${out}.rank${pid}.h5" }
}
wipe
