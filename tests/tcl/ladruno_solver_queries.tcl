# Ladruno (#729): classic-Tcl gate for the two solver-state query commands,
# `ladrunoDR` (ADR-31 rung 5) and `ladrunoArcLength` (ADR-20 Layer-B / rung 4).
#
# WHY A TCL DECK AND NOT A PYTHON TEST: the defect is classic-Tcl-only. Both
# commands answered correctly from openseespy and from the TclWrapper engine and
# were simply absent from SRC/tcl/commands.cpp, so no Python test could ever have
# caught it. This deck runs under OpenSees.exe -- the engine that had the gap.
#
# WHY IT CHECKS FOR EMPTY AND NOT JUST FOR ERRORS: the tempting fix (bridge to the
# no-arg OPS_Ladruno*Cmd(), as the contact family in #726 does) COMPILES, LINKS and
# runs -- and is silently useless. Those forms read the `cmds` singleton, which only
# the DL command engine constructs, so under OpenSees.exe they hit `cmds == 0`,
# return 0 (SUCCESS) and write NOTHING. The command would exist, return TCL_OK, and
# hand back an empty string. Every check below therefore asserts the result is
# non-empty and numeric BEFORE it asserts anything about the value.
#
# Run standalone:  dist\bin\OpenSees.exe tests\tcl\ladruno_solver_queries.tcl
# Driven by:       tests/test_ladruno_solver_queries_tcl.py
# Exits 0 on pass, 1 on any failure (and prints FAIL lines naming the check).

set failures 0

# `cond` arrives as an unevaluated expression string, so it must be run through
# expr in the CALLER's scope — the conditions below reference the caller's locals.
proc check {name cond {detail ""}} {
    global failures
    if {[uplevel 1 [list expr $cond]]} {
        puts "PASS $name"
    } else {
        incr failures
        puts "FAIL $name $detail"
    }
}

# A result that is present, numeric and finite. The empty string is the signature
# of the cmds==0 silent path, so it is called out by name.
proc check_numeric {name value} {
    global failures
    if {[string trim $value] eq ""} {
        incr failures
        puts "FAIL $name EMPTY RESULT (the command answered nothing -- this is the\
 cmds==0 silent path, not a missing command)"
        return 0
    }
    if {![string is double -strict [string trim $value]]} {
        incr failures
        puts "FAIL $name non-numeric result '$value'"
        return 0
    }
    puts "PASS $name numeric ($value)"
    return 1
}

# --------------------------------------------------------------------------
# 0. the commands must EXIST at all (the reported symptom was
#    `invalid command name "ladrunoDR"`)
# --------------------------------------------------------------------------
check "exists.ladrunoDR"        {[llength [info commands ladrunoDR]] == 1}
check "exists.ladrunoArcLength" {[llength [info commands ladrunoArcLength]] == 1}

# --------------------------------------------------------------------------
# 1. NEGATIVE FIRST: with no DR integrator active the query must ERROR, not
#    quietly succeed. This is what separates a real bridge from the silent path
#    -- the cmds==0 version returns 0 here too, so a deck that only tested the
#    happy path would pass on a broken build.
# --------------------------------------------------------------------------
wipe
model basic -ndm 2 -ndf 2
node 1 0.0 0.0
node 2 1.0 0.0
fix 1 1 1
fix 2 0 1
uniaxialMaterial Elastic 1 1000.0
element truss 1 1 2 1.0 1
timeSeries Constant 1
pattern Plain 1 1 { load 2 10.0 0.0 }
constraints Plain
numberer Plain
system BandGeneral
test NormDispIncr 1.0e-10 25 0
algorithm Newton
integrator LoadControl 1.0
analysis Static

# NB the catch runs on its OWN line: `check`'s detail argument is substituted at
# CALL time, so a $msg written inline would be read before catch ever set it.
set rc [catch {ladrunoDR residualNorm} msg]
check "guard.wrongIntegrator.DR" {$rc != 0} \
    "ladrunoDR answered '$msg' with a LoadControl integrator active"
set rc [catch {ladrunoArcLength arcLength} msg]
check "guard.wrongIntegrator.ArcLength" {$rc != 0} \
    "ladrunoArcLength answered '$msg' with a LoadControl integrator active"

# --------------------------------------------------------------------------
# 2. ladrunoDR against a live LadrunoDynamicRelaxation — the von Mises snap-through
#    truss, i.e. the DR-1 case of the Python battery
#    (tests/test_ladrunoDynamicRelaxation_integrator.py), reproduced here in the
#    engine that had the gap. The reference is computed IN THIS DECK by stock
#    LoadControl, so there is no hardcoded oracle to drift: a sub-critical constant
#    load relaxed by DR must land on the equilibrium LoadControl traces.
#
#    A LINEAR one-DOF oscillator is deliberately NOT used here. The stock
#    Gershgorin mass puts it at omega*dt = 2 exactly, where an undamped
#    single mode flip-flops instead of relaxing and parks at the 2x
#    dynamic-amplification point — real behaviour (it is the defect #728 fixes
#    with -massSafety), but nothing to do with command registration, so it would
#    make this gate fail for an unrelated reason.
# --------------------------------------------------------------------------
proc build_arch {} {
    wipe
    model basic -ndm 2 -ndf 2
    node 1 -1.0 0.0
    node 2  0.0 0.1
    node 3  1.0 0.0
    fix 1 1 1
    fix 3 1 1
    fix 2 1 0
    uniaxialMaterial Elastic 1 1.0e4
    element corotTruss 1 1 2 1.0 1
    element corotTruss 2 2 3 1.0 1
}

# reference: stock LoadControl to lambda = 3.0 (sub-critical; limit load ~3.80)
build_arch
timeSeries Linear 1
pattern Plain 1 1 { load 2 0.0 -1.0 }
system BandGeneral
numberer Plain
constraints Plain
test NormDispIncr 1.0e-10 50 0
algorithm Newton
integrator LoadControl 0.05
analysis Static
while {[getTime] < 3.0 - 1.0e-9} {
    if {[analyze 1] != 0} { break }
}
set refUy [nodeDisp 2 2]
check "DR.reference.traced" {abs([getTime] - 3.0) < 1.0e-6} \
    "LoadControl reference stopped at lambda [getTime], expected 3.0"

# the same load, relaxed by DR
build_arch
timeSeries Constant 1
pattern Plain 1 1 { load 2 0.0 -3.0 }
constraints Plain
numberer Plain
system Diagonal
test NormUnbalance 1.0e30 1 0
algorithm Linear
integrator LadrunoDynamicRelaxation
analysis Transient

analyze 1 1.0
set r0 [ladrunoDR residualNorm]
check_numeric "DR.residualNorm.first" $r0

analyze 2000 1.0
set r1 [ladrunoDR residualNorm]
set ke [ladrunoDR kineticEnergy]
set uy [nodeDisp 2 2]
check_numeric "DR.residualNorm.relaxed" $r1
check_numeric "DR.kineticEnergy"        $ke

# physics, not just plumbing. The query has to be reading the LIVE integrator, so
# the value must move with the model and land where the other engine lands.
check "DR.residual.started.loaded" {$r0 > 1.0} \
    "first-step residual $r0 should be order of the applied load (3.0)"
check "DR.residual.collapsed" {$r1 < 1.0e-9 * $r0} \
    "residual $r0 -> $r1 did not collapse"
check "DR.matches.loadcontrol" {abs($uy - $refUy) < 1.0e-6} \
    "DR uy = $uy vs LoadControl $refUy"
check "DR.kineticEnergy.nonNegative" {$ke >= 0.0} "KE = $ke"

set rc [catch {ladrunoDR notASubcommand} msg]
check "DR.unknownSubcommand.errors" {$rc != 0} \
    "an unknown subcommand answered '$msg'"
set rc [catch {ladrunoDR} msg]
check "DR.noSubcommand.errors" {$rc != 0} \
    "a bare ladrunoDR answered '$msg'"

# --------------------------------------------------------------------------
# 3. ladrunoArcLength against a live LadrunoArcLength: read a query, then drive
#    a MUTATOR, which proves the value argument reaches the elementAPI cursor
#    through the classic-Tcl argv path and not only the subcommand string.
#    Snap-through truss (tests/test_ladrunoArcLength_integrator.py's model).
# --------------------------------------------------------------------------
wipe
model basic -ndm 2 -ndf 2
node 1 -1.0 0.0
node 2  0.0 0.1
node 3  1.0 0.0
fix 1 1 1
fix 3 1 1
fix 2 1 0
uniaxialMaterial Elastic 1 1.0e4
element corotTruss 1 1 2 1.0 1
element corotTruss 2 2 3 1.0 1
timeSeries Linear 1
pattern Plain 1 1 { load 2 0.0 -1.0 }
constraints Plain
numberer Plain
system BandGeneral
test NormDispIncr 1.0e-10 50 0
algorithm Newton
integrator LadrunoArcLength 0.02 1
analysis Static

analyze 1
set al [ladrunoArcLength arcLength]
check_numeric "AL.arcLength" $al
check_numeric "AL.deltaLambdaStep" [ladrunoArcLength deltaLambdaStep]
check_numeric "AL.sign"            [ladrunoArcLength sign]
check_numeric "AL.deltaUstepNorm"  [ladrunoArcLength deltaUstepNorm]

# no subcommand at all => the current arc length (documented behaviour)
set alBare [ladrunoArcLength]
check_numeric "AL.bare" $alBare
check "AL.bare.matches" {abs($alBare - $al) < 1.0e-12} "$alBare vs $al"

# mutator: halve the radius and check the echoed value. A bridge that dropped the
# VALUE argument would fail here while passing every query above.
set alHalf [ladrunoArcLength reduceStep 0.5]
check_numeric "AL.reduceStep.echo" $alHalf
check "AL.reduceStep.halved" {abs($alHalf - 0.5*$al) < 1.0e-9 * $al} \
    "reduceStep 0.5: $al -> $alHalf"

set rc [catch {ladrunoArcLength notASubcommand} msg]
check "AL.unknownSubcommand.errors" {$rc != 0} \
    "an unknown subcommand answered '$msg'"

# --------------------------------------------------------------------------
puts ""
if {$failures == 0} {
    puts "SELF-TEST: PASS - the classic Tcl engine answers ladrunoDR and ladrunoArcLength"
    exit 0
} else {
    puts "SELF-TEST: FAIL - $failures check(s) failed"
    exit 1
}
