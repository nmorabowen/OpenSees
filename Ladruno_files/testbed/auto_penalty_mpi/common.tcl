# Shared model pieces for the `constraints Auto` cross-rank reduction test.
# Sourced by every deck so the serial references and the 2-rank deck cannot
# drift apart.
#
# The whole design exists to make ONE number observable and unambiguous: the
# `PVAL` that `constraints Auto -verbose` prints. AutoConstraintHandler sizes it
# as 10^(round(log10(KAVG)) + oom), where KAVG is the mean element-diagonal
# stiffness over the elements THIS rank can see.
#
# The element populations below are chosen so that the three possible answers
# are four orders of magnitude apart and cannot be confused:
#
#   soft part alone   : 100 trusses, k = 1e2  -> KAVG = 1e2 -> PVAL = 1e5
#   stiff part alone  :   1 truss,   k = 1e6  -> KAVG = 1e6 -> PVAL = 1e9
#   both together     : KAVG = (100*2*1e2 + 1*2*1e6)/202 = 1e4 -> PVAL = 1e7
#
# (Each truss contributes exactly 2 diagonal entries above the handler's
# 1e-12*kmax tolerance, so the counts are 200 and 2.)
#
# A 2-rank run must therefore print PVAL = 1e7 on BOTH ranks. If the cross-rank
# reduction is missing, it prints 1e5 on one and 1e9 on the other -- each rank
# sizing its penalty from its own elements. That is a four-orders-of-magnitude
# spread, not a rounding difference.

# --- the soft part: 100 unit trusses in a chain, k = E*A/L = 1e2 ------------
proc buildSoft {} {
    uniaxialMaterial Elastic 1 100.0
    for {set i 0} {$i <= 100} {incr i} {
        node [expr {1000 + $i}] [expr {1.0 * $i}] 0.0
    }
    for {set i 0} {$i < 100} {incr i} {
        element truss [expr {2000 + $i}] [expr {1000 + $i}] [expr {1001 + $i}] 1.0 1
    }
    fix 1000 1 1
    for {set i 1} {$i <= 100} {incr i} { fix [expr {1000 + $i}] 0 1 }
}

# --- the stiff part: a single truss, k = 1e6 --------------------------------
proc buildStiff {} {
    uniaxialMaterial Elastic 2 1.0e6
    node 3000 0.0 -1.0
    node 3001 1.0 -1.0
    element truss 3100 3000 3001 1.0 2
    fix 3000 1 1
    fix 3001 0 1
}

# --- the MP constraint whose penalty comes from the GLOBAL statistic --------
#
# Both nodes are deliberately ELEMENTLESS. AutoConstraintHandler accumulates a
# per-node penalty from the elements attached to that node; when a constrained
# node has none, getPenaltyValue() falls back to m_global_penalty -- the value
# the cross-rank reduction produces. That fallback is exactly the partitioned
# interface case (a replicated interface node whose elements live on the other
# rank), which is why this pair, and not a node with elements, is the probe.
proc buildFallbackMP {base} {
    node [expr {$base + 0}] 0.0 -2.0
    node [expr {$base + 1}] 1.0 -2.0
    fix  [expr {$base + 0}] 1 1
    equalDOF [expr {$base + 0}] [expr {$base + 1}] 1 2
}
