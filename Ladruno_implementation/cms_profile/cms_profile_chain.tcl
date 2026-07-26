# ADR-1000 P4 section 1 -- profiling deck.
#
# The physical smoke (tests/ladruno_cms_physical_smoke.tcl) is 4 elements per
# rank: correct, but far too small for the phase timings to mean anything. This
# is the same model scaled to a size where the hierarchy phases are resolvable,
# so `-verbose` can answer "where does the time actually go" instead of us
# guessing from the source listing.
#
# Run:  mpiexec -n 4 OpenSeesMP.exe cms_profile_chain.tcl
# Size: set the env var LADRUNO_CMS_PROFILE_ELEMENTS (default 2000 per rank).
set pid [getPID]
set np [getNP]
if {$np != 4} {
    error "CMS profile deck requires exactly four MPI ranks"
}

set elementsPerRank 2000
if {[info exists env(LADRUNO_CMS_PROFILE_ELEMENTS)]} {
    set elementsPerRank $env(LADRUNO_CMS_PROFILE_ELEMENTS)
}

wipe
model BasicBuilder -ndm 1 -ndf 1

set firstElement [expr {$pid * $elementsPerRank + 1}]
set lastElement [expr {$firstElement + $elementsPerRank - 1}]
set firstNode $firstElement
set lastNode [expr {$lastElement + 1}]

for {set node $firstNode} {$node <= $lastNode} {incr node} {
    node $node [expr {double($node - 1)}]
}
if {$pid == 0} {
    fix 1 1
}

uniaxialMaterial Elastic 1 1.0
for {set element $firstElement} {$element <= $lastElement} {incr element} {
    element truss $element $element [expr {$element + 1}] 1.0 1
}

# Every free-node mass has one primary rank; interface nodes use the lower rank.
for {set node [expr {$firstElement + 1}]} {$node <= $lastNode} {incr node} {
    mass $node 1.0
}

constraints Plain
numberer ParallelRCM
system Mumps
test NormUnbalance 1.0e-12 2
algorithm Linear
integrator LoadControl 1.0
analysis Static

if {$pid == 0} {
    puts "CMS profile: $np ranks x $elementsPerRank elements = [expr {$np * $elementsPerRank}] elements"
}

set values [eigen -ladrunoCMS \
    -domainMode physical \
    -hierarchy logical -level1 2 -level2 2 \
    -modesL2 12 -modesL1 24 \
    -tol 1.0e-8 -maxEnrich 2 -maxIter 6000 \
    -maxRefineIter 160 \
    -denseMax 2000 -verbose 6]

if {$pid == 0} {
    puts "CMS profile: first eigenvalue [lindex $values 0]"
    # Analytic chain spectrum, as a correctness anchor so a profile can never be
    # taken from a run that silently produced garbage.
    set freeDOFs [expr {$np * $elementsPerRank}]
    set reference [expr {2.0 - 2.0*cos(acos(-1.0)/(2*$freeDOFs + 1))}]
    set relative [expr {abs([lindex $values 0] - $reference)/max(abs($reference), 1.0e-30)}]
    puts "CMS profile: mode 1 relative error $relative"
    if {$relative > 1.0e-6} {
        error "CMS profile deck did not reproduce the analytic chain spectrum"
    }
}
