# The 2-rank deck. Same total element population as `serial.tcl whole`, split
# across two ranks in the manual-domain-decomposition shape this fork is driven
# with (one plain Domain per rank inside a getPID guard).
#
#   rank 0 : the 100 soft trusses   (rank-local KAVG = 1e2)
#   rank 1 : the 1 stiff truss      (rank-local KAVG = 1e6)
#
# Each rank also declares its own elementless MP pair, so on each rank the
# penalty actually used is m_global_penalty.
#
# PASS: both ranks print PVAL = 1.000000e+07 -- the whole-model value.
# FAIL: rank 0 prints 1.000000e+05 and rank 1 prints 1.000000e+09, i.e. each
#       rank sized its penalty from its own elements. That is the behaviour of
#       every binary built before the reduction was moved to a per-target TU,
#       because the `#ifdef _PARALLEL_*` block holding it was compiled out of
#       the OPS_Analysis object library.
#
# The two ranks' verbose reports interleave on stdout in an arbitrary order, so
# the checker does not try to attribute a line to a rank: it asserts that BOTH
# PVAL lines read 1.000000e+07. That is strictly stronger than checking either
# one, and it is order-independent.
wipe
model basic -ndm 2 -ndf 2
source [file join [file dirname [info script]] common.tcl]

set pid [getPID]

if {$pid == 0} {
    buildSoft
    buildFallbackMP 9000
}
if {$pid == 1} {
    buildStiff
    buildFallbackMP 9100
}

constraints Auto -verbose
numberer ParallelPlain
system Mumps
test NormDispIncr 1.0e-10 50 0
algorithm Linear
integrator LoadControl 1.0
analysis Static

set ok [analyze 1]
puts [format "AUTOPEN mp pid=%d analyze=%d" $pid $ok]
wipe
