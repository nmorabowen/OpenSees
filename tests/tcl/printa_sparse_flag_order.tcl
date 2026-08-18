# printA -sparse flag-order deck (classic Tcl engine).
#
# `-sparse <baseIndex>` is OPTIONAL, but the read used to be unconditional in
# BOTH engines: `printA -sparse -ret` consumed `-ret` as the base index and
# failed with "failed to read -sparse <baseIndex>". The only working order was
# `printA -ret -sparse 0`, which makes the flags look positional when they are
# not.
#
# The classic-Tcl path had it slightly worse than the Python one: `currentArg++`
# ran BEFORE the bounds test, so a trailing `-sparse` with nothing after it only
# survived by luck (the cursor landed on argc and the loop ended).
#
# Prints PASS/FAIL lines and a SELF-TEST verdict; the pytest wrapper asserts on
# the count, so a deck that dies early cannot masquerade as success.

set nfail 0
proc check {name ok detail} {
    global nfail
    # `ok` arrives as an already-evaluated 0/1 -- the caller uses
    # `[expr {...}]`, because a braced expression passed as a string would be
    # evaluated HERE, where the caller's locals are out of scope.
    if {$ok} {
        puts "PASS $name"
    } else {
        incr nfail
        puts "FAIL $name -- $detail"
    }
}

# --- a small 2D deck with a known equation count -------------------------
model basic -ndm 2 -ndf 2
node 1 0.0 0.0
node 2 1.0 0.0
node 3 1.0 1.0
node 4 0.0 1.0
fix 1 1 1
fix 2 1 1
nDMaterial ElasticIsotropic 1 1000.0 0.25
element quad 1 1 2 3 4 1.0 PlaneStrain 1
timeSeries Linear 1
pattern Plain 1 1 { load 3 0.0 -1.0 }
constraints LadrunoContact
numberer Plain
system FullGeneral
test NormDispIncr 1.0e-10 20 0
algorithm Newton
integrator LoadControl 1.0
analysis Static
analyze 1

# --- the orderings -------------------------------------------------------
# The pre-fix engine returns TCL_ERROR on the `-sparse -ret` form, so each call
# is wrapped: `catch` turns the refusal into a testable outcome instead of
# aborting the deck (which would look identical to "the deck never ran").

set rc0 [catch {printA -ret -sparse 0} outA]      ;# the old working order
check "ret_then_sparse0" [expr {$rc0 == 0}] "printA -ret -sparse 0 failed: $outA"

set rc1 [catch {printA -sparse -ret} outB]        ;# natural order; was broken
check "sparse_then_ret" [expr {$rc1 == 0}] "printA -sparse -ret failed: $outB"

set rc2 [catch {printA -ret -sparse 1} outC]      ;# explicit 1-based index
check "ret_sparse_base1" [expr {$rc2 == 0}] "printA -ret -sparse 1 failed: $outC"

# A bare trailing `-sparse` (index omitted, nothing after it) must still work.
set rc3 [catch {printA -sparse} outD]
check "bare_trailing_sparse" [expr {$rc3 == 0}] "printA -sparse failed: $outD"

puts "SELF-TEST: $nfail failure(s)"
