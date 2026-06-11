# Flag-order regression deck (classic Tcl exe) — `-G energy` in a NON-final
# position. Guards the elementAPI_TCL cursor contract: a failed OPS_GetIntInput
# does NOT consume the arg there (unlike openseespy), so the -G region-tag
# scan must not blind-rewind when it hits the next flag. Before the fix this
# deck died with a parse error on the recorder line; `-T nsteps 1 -G energy 1`
# (energy last) always worked. openseespy coverage of the same ordering lives
# in energy_model.py.
#
# Usage:  OpenSees.exe flag_order_model.tcl <out_dir>
# Prints FLAG_ORDER_MODEL OK on success.

if {$argc < 1} { puts "usage: flag_order_model.tcl <out_dir>"; exit 1 }
set out_dir [lindex $argv 0]
set out_ladruno [file join $out_dir flag_order.ladruno]
if {[file exists $out_ladruno]} { file delete $out_ladruno }

# undamped SDOF axial spring-mass, suddenly-applied load (same as energy_model.py)
set E 1.0e3
set A 1.0
set L 1.0
set m 10.0
set P 1.0
set pi 3.14159265358979323846
set k [expr {$E*$A/$L}]
set w [expr {sqrt($k/$m)}]
set T [expr {2.0*$pi/$w}]
set dt [expr {$T/100.0}]
set nsteps 20

wipe
model basic -ndm 2 -ndf 2
node 1 0.0 0.0
node 2 $L 0.0
fix 1 1 1
fix 2 0 1
uniaxialMaterial Elastic 1 $E
element truss 1 1 2 $A 1
mass 2 $m $m
region 1 -ele 1

# THE point of this deck: -G energy <tag> followed by another option (-T).
recorder ladruno $out_ladruno -N displacement -G energy 1 -T nsteps 1

timeSeries Constant 1
pattern Plain 1 1 { load 2 $P 0.0 }

constraints Transformation
numberer RCM
system BandSPD
test NormDispIncr 1.0e-12 10
algorithm Linear
integrator Newmark 0.5 0.25
analysis Transient
analyze $nsteps $dt
wipe

if {![file exists $out_ladruno]} { puts "FAIL: no .ladruno written"; exit 1 }
puts "FLAG_ORDER_MODEL OK"
