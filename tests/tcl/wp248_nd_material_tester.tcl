# Ladruno WP-248: classic-Tcl gate for `model testNDMaterial` and its `ndTest` command.
#
# Oracles are closed forms: Hooke's law in 3D, plane stress and a triaxial
# (confined) state; the uniaxial yield stress of perfectly plastic J2; elastic
# unloading after yield; commit and revert; and the refusals.
#
# Run standalone:  dist\bin\OpenSees.exe tests\tcl\wp248_nd_material_tester.tcl
# Driven by:       tests/test_wp248_nd_material_tester_tcl.py

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
proc near {a b {rel 1e-9}} {
    expr {abs($a - $b) <= $rel * max(1.0, abs($b))}
}
proc refused {script pattern} {
    if {![catch {uplevel 1 $script} msg]} { return 0 }
    return [string match $pattern $msg]
}

set E 200000.0
set nu 0.3
set K [expr {$E / (3.0 * (1.0 - 2.0 * $nu))}]
set G [expr {$E / (2.0 * (1.0 + $nu))}]
set sy 420.0

wipe
model testNDMaterial
nDMaterial ElasticIsotropic 1 $E $nu
nDMaterial J2Plasticity 2 $K $G $sy $sy 0.0 0.0

# --- before a material is chosen, every query is refused
check "query before material is refused" {[refused {ndTest getStress} "ndTest: no material under test*"]}

# --- 3D elastic: order 6, Hooke under strain control, the tangent is C
check "ThreeDimensional order" {[ndTest material 1 ThreeDimensional] == 6}
ndTest strain 1e-3 0 0 0 0 0
set s [ndTest getStress]
set lam [expr {$E * $nu / ((1.0 + $nu) * (1.0 - 2.0 * $nu))}]
check "3D strain control s11" {[near [lindex $s 0] [expr {($lam + 2.0 * $G) * 1e-3}]]} $s
check "3D strain control s22" {[near [lindex $s 1] [expr {$lam * 1e-3}]]} $s
set C [ndTest getTangent]
check "3D tangent has 36 entries" {[llength $C] == 36}
check "3D tangent C11" {[near [lindex $C 0] [expr {$lam + 2.0 * $G}]]} $C
check "3D tangent C44 (engineering shear)" {[near [lindex $C 21] $G]} $C

# --- mixed control: uniaxial stress in 3D
set it [ndTest mixed esssss 1e-3 0 0 0 0 0]
set s [ndTest getStress]
set e [ndTest getStrain]
check "uniaxial stress s11 = E e11" {[near [lindex $s 0] [expr {$E * 1e-3}]]} $s
check "uniaxial stress lateral strain = -nu e11" {[near [lindex $e 1] [expr {-$nu * 1e-3}]] && [near [lindex $e 2] [expr {-$nu * 1e-3}]]} $e
check "uniaxial stress lateral stress zero" {abs([lindex $s 1]) < 1e-9 && abs([lindex $s 2]) < 1e-9} $s
check "elastic mixed control takes one Newton step" {$it == 1} $it

# --- mixed control: confinement (triaxial), s2 = s3 = -p, shear stresses zero
set p 10.0
ndTest mixed essess 1e-3 [expr {-$p}] [expr {-$p}] 0 0 0
set s [ndTest getStress]
check "triaxial s1 = E e1 - 2 nu p" {[near [lindex $s 0] [expr {$E * 1e-3 - 2.0 * $nu * $p}]]} $s
check "triaxial confinement held" {[near [lindex $s 1] [expr {-$p}]] && [near [lindex $s 2] [expr {-$p}]]} $s

# --- plane stress view of the same material
check "PlaneStress order" {[ndTest material 1 PlaneStress] == 3}
ndTest mixed ess 1e-3 0 0
set s [ndTest getStress]
set e [ndTest getStrain]
check "plane stress uniaxial s11" {[near [lindex $s 0] [expr {$E * 1e-3}]]} $s
check "plane stress uniaxial e22" {[near [lindex $e 1] [expr {-$nu * 1e-3}]]} $e

# --- J2, perfectly plastic: uniaxial stress plateaus at the yield stress
ndTest material 2 ThreeDimensional
set ey [expr {$sy / $E}]
foreach f {0.5 1.0 1.5 2.0 3.0} {
    ndTest mixed esssss [expr {$f * $ey}] 0 0 0 0 0
    ndTest commit
}
set s [ndTest getStress]
# 1e-7: J2Plasticity's own return mapping stops at about 1e-8 relative
check "J2 uniaxial yield plateau" {[near [lindex $s 0] $sy 1e-7]} $s
set plateau [lindex $s 0]
set eCommitted [ndTest getStrain]

# --- a trial past the committed state, then revert: the committed state is back
ndTest mixed esssss [expr {5.0 * $ey}] 0 0 0 0 0
ndTest revert
set eBack [ndTest getStrain]
check "revert restores the committed strain" {[near [lindex $eBack 0] [lindex $eCommitted 0]] && [near [lindex $eBack 1] [lindex $eCommitted 1]]} "$eBack vs $eCommitted"

# --- unloading after yield is elastic, slope E under uniaxial stress
ndTest mixed esssss [expr {3.0 * $ey - 1e-4}] 0 0 0 0 0
set s [ndTest getStress]
check "J2 unloads at E" {[near [expr {$plateau - [lindex $s 0]}] [expr {$E * 1e-4}] 1e-8]} $s

# --- trials do not commit: the same trial twice gives the same stress
ndTest revert
ndTest mixed esssss [expr {4.0 * $ey}] 0 0 0 0 0
set s1 [lindex [ndTest getStress] 0]
ndTest mixed esssss [expr {4.0 * $ey}] 0 0 0 0 0
set s2 [lindex [ndTest getStress] 0]
check "repeated trial is repeatable" {$s1 == $s2} "$s1 $s2"

# --- refusals, each naming what was wrong
check "unknown view refused" {[refused {ndTest material 1 Bogus} "ndTest: nDMaterial 1 (*) has no Bogus view"]}
check "unknown tag refused" {[refused {ndTest material 99 ThreeDimensional} "ndTest: no nDMaterial with tag 99"]}
check "wrong strain count refused" {[refused {ndTest strain 1 2 3} "ndTest: want: ndTest strain e1 ... e6"]}
check "non-finite strain refused" {[refused {ndTest strain nan 0 0 0 0 0} "ndTest: strain \"nan\" is not a finite number"]}
check "short kinds refused" {[refused {ndTest mixed ess 1 0 0 0 0 0} "ndTest: kinds \"ess\" must have 6 letters*"]}
check "bad kind letter refused" {[refused {ndTest mixed exssss 1 0 0 0 0 0} "*each letter is e (strain) or s (stress)"]}
check "unknown option refused" {[refused {ndTest mixed esssss 1e-3 0 0 0 0 0 -bogus 1} "ndTest: unknown option -bogus"]}
check "bad -maxIter refused" {[refused {ndTest mixed esssss 1e-3 0 0 0 0 0 -maxIter 0} "ndTest: -maxIter must be*"]}
check "unknown subcommand refused" {[refused {ndTest frobnicate} "ndTest: unknown subcommand frobnicate"]}

# --- Newton that does not reach the tolerance is refused, and revert recovers
ndTest revert
check "no convergence is refused" {[refused {ndTest mixed esssss [expr {4.0 * $ey}] 0 0 0 0 0 -maxIter 1 -tol 1e-300} "ndTest: mixed control did not converge in 1 iterations*"]}
ndTest revert
check "revert after a refusal is clean" {[near [lindex [ndTest getStrain] 0] [lindex $eCommitted 0]]}

# --- model testNDMaterial takes no count: commit is explicit
check "a count is refused" {[catch {model testNDMaterial 5}]}

if {$failures == 0} {
    puts "SELF-TEST: PASS"
} else {
    puts "SELF-TEST: FAIL ($failures)"
    exit 1
}
