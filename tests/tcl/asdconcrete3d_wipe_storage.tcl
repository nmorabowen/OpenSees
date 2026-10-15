# ASDConcrete3D per-tag hardening-law registry must not survive `wipe` (apeGmsh #2044).
#
# Model 1 defines ASDConcrete3D tag 1 with E = 3e6 and is wiped. Model 2 defines
# tag 1 again with the real E = 30e6 on a 1x1x2 LadrunoBrick on rollers under a
# top load P = 7200 (LoadControl, one elastic step). Closed form (uniaxial):
#   sigma_zz = -P/A = -7200,  uz(top) = -P L/(E A) = -4.8e-4,  ux(x=1) = nu P/(E A).
# Before the fix HardeningLawStorage kept model 1's laws for tag 1, regularize()
# copied them over model 2's (E included): sigma_zz = -720 at iteration 1 and
# Newton diverged. -autoRegularization with lch_ref != lch is what makes
# regularize() reach that copy, so the deck keeps it.
# Terminal marker: "ASDC3D_WIPE_2044: PASS" or "ASDC3D_WIPE_2044: FAIL ...".

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
model BasicBuilder -ndm 3 -ndf 3
nDMaterial ASDConcrete3D 1 3.0e6 0.2 -fc 3.0e3 -ft 3.0e2
wipe

# ---- model 2: the real material under the same tag ---------------------------
set E 30.0e6
set nu 0.2
model BasicBuilder -ndm 3 -ndf 3
node 1 0.0 0.0 2.0
node 2 0.0 0.0 0.0
node 3 0.0 1.0 2.0
node 4 0.0 1.0 0.0
node 5 1.0 0.0 2.0
node 6 1.0 0.0 0.0
node 7 1.0 1.0 2.0
node 8 1.0 1.0 0.0
nDMaterial ASDConcrete3D 1 $E $nu \
    -Te 0.0 9e-05 0.00015 0.00507 0.0250501 0.250501 \
    -Ts 0.0 2700.0 3000.0 600.0 3.0 3.0 \
    -Td 0.0 0.0 0.0 0.960552268244576 0.9999800399998403 0.9999995660869531 \
    -Ce 0.0 0.0005 0.001 0.0015 0.002 0.18327272727272728 0.18377272727272728 \
    -Cs 0.0 15000.0 24852.8137423857 28923.048454132637 30000.0 3000.0 3000.0 \
    -Cd 0.0 0.0 0.0 0.0 0.0 0.9981618744961699 0.9981786141748574 \
    -autoRegularization 0.03112529146138255
element LadrunoBrick 2 3 1 2 4 7 5 6 8 1
# rollers: uz=0 on the base, ux=0 on x=0, uy=0 on y=0
fix 2 1 1 1
fix 4 1 0 1
fix 6 0 1 1
fix 8 0 0 1
fix 1 1 1 0
fix 3 1 0 0
fix 5 0 1 0
timeSeries Linear 1
pattern Plain 1 1 {
    foreach n {1 3 5 7} { load $n 0.0 0.0 -1800.0 }
}
constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1e-10 20 0
algorithm Newton
integrator LoadControl 1.0
analysis Static

set P 7200.0
set ok [analyze 1]
if {$ok == 0} { puts "PASS converge" } else { puts "FAIL converge ok=$ok"; incr nFail }
for {set gp 1} {$gp <= 8} {incr gp} {
    check "sigma_zz(gp$gp)" [lindex [eleResponse 2 material $gp stress] 2] [expr {-$P}] 1e-6
}
foreach n {1 3 5 7} { check "uz($n)" [nodeDisp $n 3] [expr {-$P * 2.0 / $E}] 1e-6 }
foreach n {5 7} { check "ux($n)" [nodeDisp $n 1] [expr {$nu * $P / $E}] 1e-6 }
reactions
set rz 0.0
foreach n {2 4 6 8} { set rz [expr {$rz + [nodeReaction $n 3]}] }
check "sum Rz(base)" $rz $P 1e-6

if {$nFail == 0} {
    puts "ASDC3D_WIPE_2044: PASS"
} else {
    puts "ASDC3D_WIPE_2044: FAIL ($nFail checks)"
}
