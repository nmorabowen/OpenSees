# ADR-75 P1i: the Tcl twin of the run-header fix. `profiler report` exists TWICE
# in SRC/tcl/commands.cpp and twice in OpenSeesCommands.cpp; wiring one has
# historically not wired the other (banked trap), so the Tcl path gets its own probe.
wipe
model basic -ndm 3 -ndf 3
nDMaterial ElasticIsotropic 1 200000.0 0.3
set n 0
for {set k 0} {$k <= 2} {incr k} {
  for {set j 0} {$j <= 2} {incr j} {
    for {set i 0} {$i <= 2} {incr i} {
      incr n; node $n [expr $i*100.0] [expr $j*100.0] [expr $k*100.0]
    } } }
for {set i 1} {$i <= 9} {incr i} { fix $i 1 1 1 }
set e 0
for {set k 0} {$k < 2} {incr k} {
  for {set j 0} {$j < 2} {incr j} {
    for {set i 0} {$i < 2} {incr i} {
      set b [expr 1 + $i + 3*($j + 3*$k)]
      incr e
      element stdBrick $e $b [expr $b+1] [expr $b+4] [expr $b+3] \
        [expr $b+9] [expr $b+10] [expr $b+13] [expr $b+12] 1
    } } }
timeSeries Linear 1
pattern Plain 1 1 { load 27 1000.0 0.0 -1000.0 }
constraints Plain
numberer RCM
system Pardiso -matrixType 2
test NormDispIncr 1e-8 25
algorithm Newton
integrator LoadControl 0.5
analysis Static
profiler start -perStep
analyze 2
profiler stop
profiler report p1i_tcl_probe.h5
puts "TCL-PROBE-DONE nElem=[getNumElements]"
