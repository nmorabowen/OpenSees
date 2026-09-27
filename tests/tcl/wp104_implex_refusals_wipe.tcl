# Ladruno (WP-104): classic-Tcl gate -- `wipe` zeroes LadrunoSANISAND's
# process-wide IMPL-EX ledger (`implexRefusals` slots 0-3 and 5,
# `implexGuards`, `avgImplexError`).
#
# THE DEFECT (apeGmsh live test, 2026-09-15, fork tip 634824e1f): a FRESH
# `nDMaterial LadrunoSANISAND` in a NEW model, after `wipe` and with a
# different tag, reported implexRefusals = {9 0 0 9 0 9} -- the nine companion
# refusals an EARLIER model in the same process had latched.  The counters live
# in one anonymous-namespace singleton in LadrunoSANISAND.cpp; nothing in
# `wipe` reached it.
#
# WHY A TCL DECK AS WELL AS THE PYTEST: `wipe` in SRC/tcl/commands.cpp and
# `wipe` in SRC/interpreter/OpenSeesCommands.cpp are two code paths.  The fix
# is hooked into OPS_clearAllNDMaterial(), which both call -- this deck proves
# the classic-Tcl one actually does.  (Same shape as
# tests/tcl/wp103_getstringfromall.tcl.)
#
# DECK: the WP-99 starved-cap LadrunoQuad (`-maxSubsteps 2`) from
# tests/test_ladrunoQuad_sanisand_implex_commit_refusal.py, transcribed.
# Model A latches a genuine commit-time companion refusal; `wipe`; model B
# (new tag, adequate cap) must read all-zero on every process-wide response
# BEFORE any analysis.  Then a second `wipe` + model A again must END with
# the same ledger model A ended with the first time (order independence).
#
# Prints SELF-TEST: PASS on success, SELF-TEST: FAIL <why> otherwise; the
# runner (tests/test_wp104_implex_refusals_wipe_tcl.py) greps for it.

set G0 264.32; set nu 0.3129; set e_init 0.6944; set Mc 1.33090; set c 0.71
set lambda_c 0.027; set e0 0.83; set ksi 0.45; set P_atm 101.0; set m 0.005
set h0 1.3; set ch 0.968; set nb 3.5; set A0 0.05; set nd 5.75; set z_max 12.5
set cz 1100.0; set Rho 2.0
set PARAMS [list $G0 $nu $e_init $Mc $c $lambda_c $e0 $ksi $P_atm $m $h0 $ch \
                 $nb $A0 $nd $z_max $cz $Rho]
set Pmin [expr {1.0e-4 * $P_atm}]

set P0 100.0
set N_CONF 5
set TOL_REL 1.0e-3
set MAXITER 60
set DQ_BIG 60.0

set XY {{0. 0.} {1. 0.} {1. 1.} {0. 1.}}

proc fail {why} {
    puts "SELF-TEST: FAIL $why"
    exit 1
}

# The same free-DOF quad as the pytest: rollers on the negative edges, the
# positive edges LOADED (a zero-free-DOF deck cannot show a refusal at all).
proc build_quad {tag maxSubsteps} {
    global PARAMS Pmin P0 N_CONF TOL_REL MAXITER XY
    wipe
    model basic -ndm 2 -ndf 2
    set j 1
    foreach pt $XY { node $j [lindex $pt 0] [lindex $pt 1]; incr j }
    eval nDMaterial LadrunoSANISAND $tag $PARAMS 1 2 1 1.0e-7 1.0e-7 \
        -Presidual 0.0 -Pmin $Pmin -implex -maxSubsteps $maxSubsteps
    element LadrunoQuad 1 1 2 3 4 $tag -thick 1.0 -type PlaneStrain -formulation bbar
    set j 1
    foreach pt $XY {
        fix $j [expr {[lindex $pt 0] == 0. ? 1 : 0}] [expr {[lindex $pt 1] == 0. ? 1 : 0}]
        incr j
    }
    set q [expr {$P0 / 2.0}]
    timeSeries Linear 1
    pattern Plain 1 1 {
        set j 1
        foreach pt $XY {
            if {[lindex $pt 0] == 1.} { load $j [expr {-$q}] 0.0 }
            if {[lindex $pt 1] == 1.} { load $j 0.0 [expr {-$q}] }
            incr j
        }
    }
    constraints Transformation
    numberer Plain
    system FullGeneral
    test NormUnbalance [expr {$TOL_REL * $P0}] $MAXITER 0
    algorithm Newton
    integrator LoadControl [expr {1.0 / $N_CONF}]
    analysis Static
}

proc confine_and_flip {tag} {
    global N_CONF
    updateMaterialStage -material $tag -stage 0
    for {set s 0} {$s < $N_CONF} {incr s} {
        if {[analyze 1] != 0} { fail "confinement step [expr {$s + 1}] failed (tag $tag)" }
    }
    loadConst -time 0.0
    updateMaterialStage -material $tag -stage 1
}

proc add_deviatoric_pattern {} {
    global XY DQ_BIG
    timeSeries Linear 2
    pattern Plain 2 2 {
        set j 1
        foreach pt $XY {
            if {[lindex $pt 1] == 1.} { load $j 0.0 [expr {-$DQ_BIG / 2.0}] }
            incr j
        }
    }
    integrator LoadControl 1.0
}

proc refusals {} { return [eleResponse 1 material 1 implexRefusals] }
proc guards   {} { return [eleResponse 1 material 1 implexGuards] }
proc avgerr   {} { return [lindex [eleResponse 1 material 1 avgImplexError] 0] }

proc is_zero_list {lst n} {
    if {[llength $lst] != $n} { return 0 }
    foreach v $lst { if {$v != 0} { return 0 } }
    return 1
}

# Model A: latch, then one more (refused) step so slot 5 moves too.
proc run_latching_model {tag} {
    build_quad $tag 2
    confine_and_flip $tag
    add_deviatoric_pattern
    set latched 0
    for {set n 1} {$n <= 6} {incr n} {
        analyze 1
        if {[lindex [refusals] 4] == 1} { set latched 1; break }
    }
    if {!$latched} { fail "the starved deck never latched in 6 steps (deck, not contract): [refusals]" }
    if {[analyze 1] == 0} { fail "the step after a latched commit converged (WP-99 contract)" }
    set r [refusals]
    if {[llength $r] != 6} { fail "implexRefusals is not the 6-slot vector: $r" }
    if {[lindex $r 3] <= 0 || [lindex $r 4] != 1 || [lindex $r 5] <= 0} {
        fail "model A did not populate companion/commitLatched/latched: $r"
    }
    return $r
}

# --- A: latch -------------------------------------------------------------
set endA1 [run_latching_model 9950]
puts "model A (first run) ended with implexRefusals = $endA1"

# --- wipe; B: fresh tag, adequate cap, NO analysis yet ----------------------
build_quad 9951 20000
set r0 [refusals]
if {![is_zero_list $r0 6]} {
    fail "a FRESH LadrunoSANISAND after wipe inherited the previous model's refusal ledger: $r0 (model A ended at $endA1)"
}
if {![is_zero_list [guards] 7]} { fail "implexGuards survived wipe: [guards]" }
if {[avgerr] != 0} { fail "avgImplexError survived wipe: [avgerr]" }
confine_and_flip 9951
if {![is_zero_list [refusals] 6]} { fail "the adequate-cap model refused during confinement: [refusals]" }
puts "model B after wipe: implexRefusals = $r0, implexGuards = [guards], avgImplexError = [avgerr]"

# --- wipe; A again: must END where the first A ended (order independence) --
set endA2 [run_latching_model 9952]
puts "model A (second run) ended with implexRefusals = $endA2"
if {$endA1 ne $endA2} {
    fail "two identical models in one process ended with different ledgers: $endA1 vs $endA2"
}

puts "SELF-TEST: PASS"
