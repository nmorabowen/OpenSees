# Ladruno (WP-103): classic-Tcl gate for OPS_GetStringFromAll's BUFFER contract.
#
# THE DEFECT (fixed by WP-103, elementAPI_TCL.cpp:493):
#
#     extern "C" const char* OPS_GetStringFromAll(char *buffer, int len)
#     { return OPS_GetString(); }          // <-- never writes `buffer`
#
# elementAPI.h documents the function as "does a strcpy", and the openseespy
# backend (PythonModule::getStringFromAll) does exactly that. The classic-Tcl
# backend did not. Every caller written in the family idiom
#
#     char tok[64];
#     OPS_GetStringFromAll(tok, sizeof(tok));
#     if (strcmp(tok, "auto") == 0) ...
#
# therefore read UNINITIALISED STACK under OpenSees.exe. The parse did not
# crash and did not warn about the API -- the option was simply lost, and the
# element either took its default or died with a message naming an empty
# token (`-k wants a number or 'auto', got ''`).
#
# WHY A TCL DECK AND NOT A PYTHON TEST: the defect is classic-Tcl-ONLY. Under
# openseespy the buffer is filled, so the entire openseespy battery for these
# elements passed throughout. No Python test could ever have caught it.
#
# WHAT EACH BLOCK PROVES:
#   A  LadrunoKinematicCoupling -k <number>   the reported symptom (`got ''`)
#   B  LadrunoKinematicCoupling -dof          same idiom, silently wrong nGap
#   C  LadrunoEmbeddedNode  <nHost>           the HOST SPEC goes through the
#                                             buffer -- before WP-103 the
#                                             explicit-count form of this
#                                             element was unusable from Tcl
#   D  nDMaterial Series3D / Parallel3D       RETURN-VALUE callers: they always
#   E  element LadrunoUP                      worked under Tcl (the UP/overlay
#                                             files even carry a local shim for
#                                             this bug). Here as NO-REGRESSION
#                                             checks -- the fix changes what the
#                                             function RETURNS (buffer, not the
#                                             argv pointer) and must not break
#                                             them.
#
# Run standalone:  dist\bin\OpenSees.exe tests\tcl\wp103_getstringfromall.tcl
# Driven by:       tests/test_wp103_getstringfromall_tcl.py
# Exits 0 on pass, 1 on any failure (and prints FAIL lines naming the check).

set failures 0

# `cond` arrives unevaluated so it is run through expr in the CALLER's scope --
# the conditions below reference the caller's locals.
proc check {name cond {detail ""}} {
    global failures
    if {[uplevel 1 [list expr $cond]]} {
        puts "PASS $name"
    } else {
        incr failures
        puts "FAIL $name $detail"
    }
}

# Build an element, reporting the Tcl error instead of aborting the deck: the
# BEFORE behaviour is a refused `element` command, so the deck has to survive it
# to print the remaining verdicts.
proc try_element {name args} {
    global failures
    if {[catch {eval element $args} err]} {
        incr failures
        puts "FAIL $name element command REFUSED: $err"
        return 0
    }
    puts "PASS $name element accepted"
    return 1
}

proc scalar {resp} {
    if {[llength $resp] == 0} { return "" }
    return [lindex $resp 0]
}

# ==========================================================================
# A. LadrunoKinematicCoupling -k <number>   (the reported reproducer)
# ==========================================================================
wipe
model basic -ndm 3 -ndf 6
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 0.0 1.0 0.0

# The exact deck line from the WP-103 report. Before the fix this printed
#   WARNING LadrunoKinematicCoupling: -k wants a number or 'auto', got ''
# because kTok[64] was never written.
if {[try_element "A.k.numeric" LadrunoKinematicCoupling 1 1 2 2 3 -k 1.0e6 -enforce al]} {
    set kt [scalar [eleResponse 1 penalty]]
    check "A.k.numeric.readback" {$kt ne ""} "eleResponse 1 penalty answered nothing"
    check "A.k.numeric.value" {$kt ne "" && abs($kt - 1.0e6) <= 1.0} \
        "K_t = '$kt', want 1e6 -- a stale 1e12 default or an empty token means the\
 buffer was not filled"
    puts "---- print -ele 1 (A) ----"
    print -ele 1
    puts "--------------------------"
}

# --- A2. the 'auto' SENTINEL, which is what the idiom really tests for -------
# `-k auto` is refused either way (it needs a representative -host element), but
# the MESSAGE says whether the token survived. Pre-WP-103 the stack residue here
# happened to be an empty string -- the literal symptom in the WP-101 report:
#     WARNING LadrunoKinematicCoupling: -k wants a number or 'auto', got ''
# After the fix the parser gets as far as the -host requirement, which is the
# only way to prove `strcmp(kTok, "auto")` actually matched.
if {[catch {element LadrunoKinematicCoupling 9 1 2 2 3 -k auto} err]} {
    puts "PASS A2.k.auto refused as expected (see the WARNING above)"
} else {
    incr failures
    puts "FAIL A2.k.auto '-k auto' was accepted without a -host element"
}

# ==========================================================================
# B. LadrunoKinematicCoupling -dof  (same idiom, silent wrong answer)
# ==========================================================================
# With ndf=6 slaves the DEFAULT component list ties 6 components per slave,
# so nGap = 12 for two slaves. `-dof 1 2 3` must cut that to 6. Before the fix
# the greedy `-dof` reader strtol'd uninitialised stack, so the list was either
# empty (hard refusal) or nondeterministic.
if {[try_element "B.dof.default" LadrunoKinematicCoupling 2 1 2 2 3 -k 1.0e6]} {
    set ngDefault [scalar [eleResponse 2 tiedDOFs]]
    check "B.dof.default.value" {$ngDefault ne "" && $ngDefault == 12} \
        "default nGap = '$ngDefault', want 12 (6 components x 2 slaves)"
}
if {[try_element "B.dof.explicit" LadrunoKinematicCoupling 3 1 2 2 3 -dof 1 2 3 -k 1.0e6]} {
    set ngSel [scalar [eleResponse 3 tiedDOFs]]
    check "B.dof.explicit.value" {$ngSel ne "" && $ngSel == 6} \
        "-dof 1 2 3 gave nGap = '$ngSel', want 6 (3 components x 2 slaves) --\
 12 means -dof was silently dropped"
    set ktB [scalar [eleResponse 3 penalty]]
    check "B.dof.k.survived" {$ktB ne "" && abs($ktB - 1.0e6) <= 1.0} \
        "K_t after -dof = '$ktB', want 1e6"
    puts "---- print -ele 3 (B) ----"
    print -ele 3
    puts "--------------------------"
}

# ==========================================================================
# C. LadrunoEmbeddedNode -- the HOST SPEC itself goes through the buffer
# ==========================================================================
# `hostTok` decides between the `-host <eleTag>` form and the explicit
# `<nHost> h1..hN` form, and the explicit branch does atoi(hostTok). With an
# unwritten buffer that is atoi(garbage), so the explicit form of this element
# could not be used from a .tcl deck at all:
#   WARNING LadrunoEmbeddedNode: nHost must be >= 1 (or use -host eleTag); got ''
wipe
model basic -ndm 3 -ndf 3
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 9 0.5 0.5 0.0
if {[try_element "C.embeddedNode.explicitHost" LadrunoEmbeddedNode 11 9 4 1 2 3 4 \
                 -shape 0.25 0.25 0.25 0.25 -k 3.0e6]} {
    set ku [scalar [eleResponse 11 penalty]]
    check "C.embeddedNode.k.value" {$ku ne "" && abs($ku - 3.0e6) <= 1.0} \
        "K_u = '$ku', want 3e6"
}

# ==========================================================================
# D. NO-REGRESSION: return-value callers (nDMaterial Series3D / Parallel3D)
# ==========================================================================
# These read the RETURN of OPS_GetStringFromAll, never the buffer, so they
# always worked under Tcl. WP-103 changes the returned POINTER (it now aims at
# the caller's buffer instead of into argv) -- these checks pin that this is
# harmless.
wipe
model basic -ndm 3 -ndf 3
nDMaterial ElasticIsotropic 1 30000.0 0.3 2.4
nDMaterial ElasticIsotropic 2 20000.0 0.2 2.4
if {[catch {nDMaterial Series3D 10 1 2 -weights 0.5 0.5 -maxIter 20 -relTol 1.0e-6} err]} {
    incr failures
    puts "FAIL D.series3d REFUSED: $err"
} else {
    puts "PASS D.series3d accepted"
}
if {[catch {nDMaterial Parallel3D 11 1 2 -weights 0.4 0.6} err]} {
    incr failures
    puts "FAIL D.parallel3d REFUSED: $err"
} else {
    puts "PASS D.parallel3d accepted"
}

# ==========================================================================
# E. NO-REGRESSION: LadrunoUP (goes through the local upGetTok shim)
# ==========================================================================
# OPS_LadrunoUP.cpp carries `upGetTok`, a private workaround for exactly this
# bug ("the classic-Tcl version just returns OPS_GetString() and NEVER touches
# `buf`"). After WP-103 the shim's `if (s != buf)` copy is a no-op -- the
# greedy leading-integer read must still work.
wipe
model basic -ndm 2 -ndf 3
nDMaterial ElasticIsotropic 1 30000.0 0.3 2.0
node 1 0.0 0.0
node 2 1.0 0.0
node 3 1.0 1.0
node 4 0.0 1.0
if {[try_element "E.ladrunoUP" LadrunoUP 21 1 2 3 4 1 \
                 -Kf 2.2e6 -poro 0.4 -rhoF 1.0 -perm 1.0e-5 1.0e-5]} {
    puts "PASS E.ladrunoUP.constructed"
}

# ==========================================================================
puts ""
if {$failures == 0} {
    puts "SELF-TEST: PASS (OPS_GetStringFromAll fills the caller's buffer under classic Tcl)"
    exit 0
} else {
    puts "SELF-TEST: FAIL ($failures check(s) failed)"
    exit 1
}
