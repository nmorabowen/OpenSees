# Ladruno WP-221: classic-Tcl gate for `ladrunoSymmetryGuard` and the
# `-strictSymmetry` flag of `system ProfileSPD` / `system Pardiso`.
#
# The Python suite (tests/test_wp221_symmetry_guard.py) covers the behaviour; this
# deck covers what only OpenSees.exe can show: the classic engine registers the
# command, hands it ITS theSOE (not the DL `cmds` singleton, which would answer
# nothing -- see ladruno_solver_queries.tcl), and its own `system` parser accepts
# the new flag.
#
# Run standalone:  dist\bin\OpenSees.exe tests\tcl\wp221_symmetry_guard.tcl
# Driven by:       tests/test_wp221_symmetry_guard_tcl.py

set failures 0
proc check {name cond {detail ""}} {
    global failures
    if {[uplevel 1 [list expr $cond]]} {
        puts "PASS $name"
    } else {
        incr failures
        puts "FAIL $name $detail"
    }
}

# 2x2 elastic plane-strain-ish brick column: symmetric, with mass, for eigen.
proc build {} {
    wipe
    model basic -ndm 3 -ndf 3
    nDMaterial ElasticIsotropic 1 25000.0 0.35
    set n 0
    foreach z {0.0 1.0 2.0} {
        foreach {x y} {0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0} {
            incr n
            node $n $x $y $z
            mass $n 0.1 0.1 0.1
        }
    }
    foreach t {1 2 3 4} { fix $t 1 1 1 }
    element stdBrick 1 1 2 3 4 5 6 7 8 1
    element stdBrick 2 5 6 7 8 9 10 11 12 1
}

check "exists" {[llength [info commands ladrunoSymmetryGuard]] == 1}

# no system yet -> a Tcl error, never an empty success
build
check "no-system.errors" {[catch {ladrunoSymmetryGuard}] == 1}

# a full-storage system: no half-storage guard, but the eigen-side one answers
build
system UmfPack
eigen 2
check "umfpack.no-linear-guard" {[catch {ladrunoSymmetryGuard}] == 1}
check "umfpack.eigen.checks" {[ladrunoSymmetryGuard -eigen checks] == 2}

set systems [list [list ProfileSPD -strictSymmetry] [list ProfileSPD]]
if {![catch {system Pardiso -matrixType 2}]} {
    lappend systems [list Pardiso -matrixType 2 -strictSymmetry] [list Pardiso -matrixType 2]
}

foreach sys $systems {
    build
    eval system $sys
    set lam [eigen 2]
    check "$sys.eigen" {[llength $lam] == 2 && [lindex $lam 0] > 0} "lam=$lam"
    set rec [ladrunoSymmetryGuard]
    check "$sys.record" {[llength $rec] == 10} "rec='$rec'"
    check "$sys.checks" {[ladrunoSymmetryGuard checks] > 0} "rec='$rec'"
    check "$sys.latched" {[ladrunoSymmetryGuard latched] == 0} "rec='$rec'"
    set want [expr {[lsearch $sys -strictSymmetry] >= 0 ? 1 : 0}]
    check "$sys.strict" {[ladrunoSymmetryGuard strict] == $want} "rec='$rec'"
    # the eigen-side guard (ArpackSOE): every element K of the eigen
    check "$sys.eigen.checks" {[ladrunoSymmetryGuard -eigen checks] == 2}
    check "$sys.eigen.latched" {[ladrunoSymmetryGuard -eigen latched] == 0}
}

if {$failures == 0} {
    puts "SELF-TEST: PASS"
    exit 0
}
puts "SELF-TEST: FAIL ($failures)"
exit 1
