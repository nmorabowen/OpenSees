# ASDConcrete1D per-tag hardening-law registry must not survive `wipe` (apeGmsh #2044,
# the uniaxial sibling of asdconcrete3d_wipe_storage.tcl).
#
# Model 1 defines ASDConcrete1D tag 1 with E = 3e6 and is wiped. Model 2 defines
# tag 1 again with the real E = 30e6 in one truss (L = 2, A = 1) under an axial
# load P = 7200 in compression (LoadControl, one elastic step). Closed form:
#   sigma = -P/A = -7200,  u = -P L/(E A) = -4.8e-4.
# Before the fix HardeningLawStorage kept model 1's laws for tag 1 and
# regularize() copied them over model 2's: a singular truss, u = 0, stress 0.
# -autoRegularization with lch_ref != lch (= L) is what reaches that copy.
# Terminal marker: "ASDC1D_WIPE_2044: PASS" or "ASDC1D_WIPE_2044: FAIL ...".

set nFail 0
proc check {name got want rtol} {
    global nFail
    if {abs($got - $want) <= $rtol * abs($want)} {
        puts "PASS $name got=$got want=$want"
    } else {
        puts "FAIL $name got=$got want=$want"
        incr nFail
    }
}

# ---- model 1: same tag, a different (soft) concrete, then wipe ----------------
wipe
model BasicBuilder -ndm 2 -ndf 2
uniaxialMaterial ASDConcrete1D 1 3.0e6 -fc 3.0e3 -ft 3.0e2
wipe

# ---- model 2: the real material under the same tag ---------------------------
set E 30.0e6
set L 2.0
set P 7200.0
model BasicBuilder -ndm 2 -ndf 2
node 1 0.0 0.0
node 2 $L 0.0
fix 1 1 1
fix 2 0 1
uniaxialMaterial ASDConcrete1D 1 $E -fc 30.0e3 -ft 3.0e3 -autoRegularization 0.03112529146138255
element truss 1 1 2 1.0 1
timeSeries Linear 1
pattern Plain 1 1 { load 2 [expr {-$P}] 0.0 }
constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1e-10 20 0
algorithm Newton
integrator LoadControl 1.0
analysis Static

set ok [analyze 1]
if {$ok == 0} { puts "PASS converge" } else { puts "FAIL converge ok=$ok"; incr nFail }
check "u(2)" [nodeDisp 2 1] [expr {-$P * $L / $E}] 1e-5
check "stress" [lindex [eleResponse 1 material stress] 0] [expr {-$P}] 1e-5
reactions
check "R(1)" [nodeReaction 1 1] $P 1e-5

if {$nFail == 0} {
    puts "ASDC1D_WIPE_2044: PASS"
} else {
    puts "ASDC1D_WIPE_2044: FAIL ($nFail checks)"
}
