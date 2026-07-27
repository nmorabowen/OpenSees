# ADR-75 P2h -- Tcl `system Mumps` option-ladder smoke under OpenSeesMP.
#
# Proves by EXECUTION what the compile check could not: that the Tcl ladder now
# PARSES -BLR/-stats (and warns on a typo) instead of silently discarding them.
# Before P2h every option below was dropped with no trace.
#
# Run:  mpirun -n 2 ~/ladruno_build_test/openseesmp.sh p2h_smoke.tcl
#
# Env knobs: NX NY NZ (grid), MODE (which system line to use)

set pid [getPID]
set np  [getNP]

set NX 8 ; set NY 8 ; set NZ 16
if {[info exists env(SMOKE_NX)]} { set NX $env(SMOKE_NX) }
if {[info exists env(SMOKE_NY)]} { set NY $env(SMOKE_NY) }
if {[info exists env(SMOKE_NZ)]} { set NZ $env(SMOKE_NZ) }
set MODE "full"
if {[info exists env(SMOKE_MODE)]} { set MODE $env(SMOKE_MODE) }

set H 100.0
set E 200.0e3
set NU 0.3

set NXP [expr {$NX + 1}]
set NYP [expr {$NY + 1}]

proc nid {i j k} {
    global NXP NYP
    return [expr {1 + $i + $NXP*$j + $NXP*$NYP*$k}]
}

wipe
model basic -ndm 3 -ndf 3
nDMaterial ElasticIsotropic 1 $E $NU

# k-slab partition: rank owns [k0,k1]; the shared plane is built by both.
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
            element stdBrick $e \
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

# Ladruno ADR-75 P2h: the numberer is a CONFOUND CONTROL, not a detail.
# The first sweep used vanilla `ParallelRCM`; the fork's own mp_blr_smoke.py
# prefers `LadrunoParallelRCM` (ADR-74), whose whole reason to exist is that the
# vanilla parallel numberer degrades at scale. Numbering cost is additive and
# identical across solver modes, so it DILUTES every wall ratio toward 1.0 --
# it cannot reverse a verdict, but it can hide the size of one.
set NUMB "ParallelRCM"
if {[info exists env(SMOKE_NUMBERER)]} { set NUMB $env(SMOKE_NUMBERER) }
set t_numb0 [clock milliseconds]
numberer $NUMB
set t_numb1 [clock milliseconds]
if {$pid == 0} { puts "P2H_NUMBERER $NUMB decl_ms=[expr {$t_numb1-$t_numb0}]" }

# ---- the thing under test -------------------------------------------------
switch -- $MODE {
    full     { system Mumps -ICNTL14 200 -stats }
    blr      { system Mumps -ICNTL14 200 -BLR 1e-8 -stats }
    blrloose { system Mumps -ICNTL14 200 -BLR 1e-4 -stats }
    sym      { system Mumps -ICNTL14 200 -matrixType 2 -stats }
    symblr   { system Mumps -ICNTL14 200 -matrixType 2 -BLR 1e-8 -stats }
    typo     { system Mumps -ICNTL14 200 -BLRR 1e-8 -stats }
    trailing { system Mumps -ICNTL14 200 -stats -ICNTL7 }
    bare     { system Mumps }
    commsplit { system Mumps -ICNTL14 200 -commSplit 0 -stats }
    default  { system Mumps -ICNTL14 200 -stats }
}
if {$pid == 0} { puts "P2H_MODE $MODE" }

test NormDispIncr 1.0e-7 25 0
algorithm Newton
integrator LoadControl 0.5
analysis Static

set t0 [clock milliseconds]
set ok1 [analyze 1]
set tmid [clock milliseconds]
set ok2 [analyze 1]
set ok [expr {$ok1 + $ok2}]
set t1 [clock milliseconds]

if {$pid == 0} { puts "P2H_RESULT mode=$MODE numberer=$NUMB ok=$ok wall_ms=[expr {$t1-$t0}] step1_ms=[expr {$tmid-$t0}] step2_ms=[expr {$t1-$tmid}] np=$np" }
if {$k1 == $NZ} {
    puts [format "P2H_TIP pid=%d mode=%s ux=%.17e" $pid $MODE [nodeDisp $tip 1]]
}
wipe
