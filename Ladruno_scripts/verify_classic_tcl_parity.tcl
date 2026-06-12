# Classic-Tcl parity verification: newly wired Ladruno elements + profiler.
# Grammars mirror apeGmsh's emitters (source-verified against the OPS_ parsers).
# Any failed command aborts before the final ALL...PASS line.

# --- profiler command exists and brackets a real solve -----------------------
profiler start -deep -perStep -memory

# --- LadrunoBrick (std + ssp + eas formulations) ------------------------------
wipe
model basic -ndm 3 -ndf 3
nDMaterial ElasticIsotropic 1 30000.0 0.3 2.4
node 1 0 0 0; node 2 1 0 0; node 3 1 1 0; node 4 0 1 0
node 5 0 0 1; node 6 1 0 1; node 7 1 1 1; node 8 0 1 1
element LadrunoBrick 1 1 2 3 4 5 6 7 8 1
element LadrunoBrick 2 1 2 3 4 5 6 7 8 1 -formulation ssp
element LadrunoBrick 3 1 2 3 4 5 6 7 8 1 -formulation eas
puts "OK LadrunoBrick (std/ssp/eas)"

# tiny solve so the profiler has an analyze tree to report
fix 1 1 1 1; fix 2 1 1 1; fix 3 1 1 1; fix 4 1 1 1
timeSeries Linear 1
pattern Plain 1 1 { load 7 1.0 0.0 0.0 }
constraints Transformation; numberer RCM; system UmfPack
test NormDispIncr 1.0e-8 20; algorithm Newton; integrator LoadControl 0.5
analysis Static
analyze 2
puts "OK analyze"

# --- LadrunoQuad / LadrunoCST (2D, -thick flag) --------------------------------
wipe
model basic -ndm 2 -ndf 2
nDMaterial ElasticIsotropic 1 30000.0 0.3 2.4
node 1 0 0; node 2 1 0; node 3 1 1; node 4 0 1
element LadrunoQuad 1 1 2 3 4 1 -thick 1.0
element LadrunoQuad 2 1 2 3 4 1 -formulation ssp -thick 1.0
element LadrunoCST  3 1 2 3 1 -thick 1.0
puts "OK LadrunoQuad + LadrunoCST"

# --- LadrunoKinematicCoupling / LadrunoDistributingCoupling --------------------
wipe
model basic -ndm 3 -ndf 6
node 1 0 0 0; node 2 1 0 0; node 3 0 1 0; node 4 0.5 0.5 0
element LadrunoKinematicCoupling    1 4 2 1 2
element LadrunoDistributingCoupling 2 4 3 1 2 3
puts "OK couplings (RBE2 + RBE3)"

# --- profiler report ------------------------------------------------------------
profiler stop
profiler report _verify_profile.h5 -run classic_parity
puts "OK profiler report"

puts "ALL CLASSIC-PARITY CHECKS PASS"
