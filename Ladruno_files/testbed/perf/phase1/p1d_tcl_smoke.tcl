# ADR-75 P1d — Tcl-side `system Pardiso` smoke.
#
# P1b registered the verb only in SRC/interpreter/OpenSeesCommands.cpp, so
# `system Pardiso` worked from OpenSeesPy but NOT from OpenSees.exe (the Tcl
# `system` chain in SRC/tcl/commands.cpp is a separate if-ladder). P1d adds it
# there; this script is the proof it resolves and solves.
#
# A cantilever of elastic bricks: small, symmetric-tangent, and with a closed
# enough answer that all four spellings must agree bit-for-bit.
#
# Run:  dist\bin\OpenSees.exe p1d_tcl_smoke.tcl

set NX 6
set L  100.0

proc nid {i j k} {
    global NX
    return [expr {1 + $i + ($NX+1)*($j + ($NX+1)*$k)}]
}

proc runcase {label args} {
    global NX L
    wipe
    model basic -ndm 3 -ndf 3
    nDMaterial ElasticIsotropic 1 200000.0 0.3
    for {set k 0} {$k <= $NX} {incr k} {
      for {set j 0} {$j <= $NX} {incr j} {
        for {set i 0} {$i <= $NX} {incr i} {
          node [nid $i $j $k] [expr {$i*$L}] [expr {$j*$L}] [expr {$k*$L}]
        }
      }
    }
    for {set j 0} {$j <= $NX} {incr j} {
      for {set i 0} {$i <= $NX} {incr i} { fix [nid $i $j 0] 1 1 1 }
    }
    set e 1
    for {set k 0} {$k < $NX} {incr k} {
      for {set j 0} {$j < $NX} {incr j} {
        for {set i 0} {$i < $NX} {incr i} {
          element stdBrick $e \
            [nid $i $j $k] [nid [expr {$i+1}] $j $k] \
            [nid [expr {$i+1}] [expr {$j+1}] $k] [nid $i [expr {$j+1}] $k] \
            [nid $i $j [expr {$k+1}]] [nid [expr {$i+1}] $j [expr {$k+1}]] \
            [nid [expr {$i+1}] [expr {$j+1}] [expr {$k+1}]] \
            [nid $i [expr {$j+1}] [expr {$k+1}]] 1
          incr e
        }
      }
    }
    timeSeries Linear 1
    pattern Plain 1 1 {
      for {set j 0} {$j <= $NX} {incr j} {
        for {set i 0} {$i <= $NX} {incr i} { load [nid $i $j $NX] 4.0e5 0.0 -4.0e5 }
      }
    }
    constraints Plain
    numberer RCM
    eval system $args
    test NormDispIncr 1e-8 20
    algorithm Newton
    integrator LoadControl 1.0
    analysis Static
    if {[analyze 1] != 0} { puts [format "%-26s FAILED" $label] ; return }
    set ux [nodeDisp [nid [expr {$NX/2}] [expr {$NX/2}] $NX] 1]
    puts [format "%-26s ux = %.12e" $label $ux]
}

runcase "UmfPack"                 UmfPack
runcase "Pardiso (unsym)"         Pardiso
runcase "Pardiso -matrixType 2"   Pardiso -matrixType 2
runcase "Pardiso -matrixType 1"   Pardiso -matrixType 1
runcase "Pardiso -symmetric"      Pardiso -symmetric
runcase "Pardiso -spd"            Pardiso -spd
puts "TCL-SMOKE-DONE"
wipe
exit
