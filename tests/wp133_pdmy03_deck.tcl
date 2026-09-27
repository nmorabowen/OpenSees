# WP-133 (TIMs F23a) -- Tcl-ladder twin of tests/wp133_pdmy03_deck.py
# (case "default_tail"). Prints the per-step stress and strain with %.17g,
# which round-trips a double exactly, so two runs can be diffed byte-for-byte.
#   OpenSees.exe tests/wp133_pdmy03_deck.tcl ?extra PDMY03 args...?
# Extra args (e.g. -cs1 0.62) are appended after the positional tail.

set extra $argv
wipe
model basic -ndm 2 -ndf 2
eval nDMaterial PressureDependMultiYield03 1 2 2.0 1.3e5 2.6e5 40.0 0.1 101.0 0.5 26.0 \
    0 0.013 0.0 0.3 0.0 0.0 0.3 3.0 0.0 20 1.0 0.0 101.0 0.1 $extra
node 1 0.0 0.0
node 2 1.0 0.0
node 3 1.0 1.0
node 4 0.0 1.0
fix 1 1 1
fix 2 0 1
fix 4 1 0
equalDOF 3 4 2
element quad 1 1 2 3 4 1.0 PlaneStrain 1
timeSeries Linear 1
pattern Plain 1 1 {
    load 2 -50.0 0.0
    load 3 -50.0 -50.0
    load 4 0.0 -50.0
}
constraints Transformation
numberer Plain
system FullGeneral
test NormDispIncr 1e-9 100 0
algorithm KrylovNewton
integrator LoadControl 0.1
analysis Static
if {[analyze 10] != 0} { puts "FAIL stage0"; exit 1 }
integrator LoadControl 0.0
loadConst -time 0.0
updateMaterialStage -material 1 -stage 1
if {[analyze 1] != 0} { puts "FAIL stage1"; exit 1 }
loadConst -time 0.0
timeSeries Linear 2
pattern Plain 2 2 { sp 3 2 1.0 }
integrator LoadControl -5.0e-4
analysis Static
for {set i 0} {$i < 100} {incr i} {
    if {[analyze 1] != 0} { puts "STOP $i"; break }
    set line "$i"
    foreach v [eleResponse 1 stress] { append line " " [format %.17g $v] }
    foreach v [eleResponse 1 strain] { append line " " [format %.17g $v] }
    puts $line
}
