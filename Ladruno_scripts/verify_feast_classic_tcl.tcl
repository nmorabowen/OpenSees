# ADR-43 classic-Tcl -feast re-smoke: stale-band fix verification.
# 2-DOF truss chain, k=100, m=1 -> lambda = 100*(3+-sqrt(5))/2 = {38.1966, 261.8034}
# f1 = 0.9836 Hz, f2 = 2.5753 Hz.
wipe
model basic -ndm 1 -ndf 1
node 1 0.0
node 2 1.0
node 3 2.0
fix 1 1
uniaxialMaterial Elastic 1 100.0
element truss 1 1 2 1.0 1
element truss 2 2 3 1.0 1
mass 2 1.0
mass 3 1.0
system UmfPack

set fp [open feast_resmoke.out w]
# 1. full band -> expect BOTH modes
puts $fp "full_band=[eigen -feast 0.0 5.0]"
# 2. SECOND call, NARROW band -> expect ONLY mode 1 (38.1966). Stale-band
#    bug would return both.
puts $fp "narrow_low=[eigen -feast 0.0 2.0]"
# 3. THIRD call, HIGH band -> expect ONLY mode 2 (261.8034). Stale-band
#    bug would return mode 1 / empty.
puts $fp "high_band=[eigen -feast 2.0 10.0]"
# 4. EMPTY band -> expect nothing.
puts $fp "empty_band=[eigen -feast 3.0 10.0]"
# 5. FEAST -> ARPACK switch (classTag mismatch path) -> expect 1 value.
puts $fp "arpack_after=[eigen 1]"
# 6. ARPACK -> FEAST switch back -> expect BOTH modes again.
puts $fp "feast_after=[eigen -feast 0.0 5.0]"
close $fp
puts "RESMOKE_DONE"
