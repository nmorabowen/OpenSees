# ADR 1000 P3-D04/P3-D08: one physical fine subdomain per MPI rank.
set pid [getPID]
set np [getNP]
if {$np != 4} {
    error "ladruno CMS physical smoke requires exactly four MPI ranks"
}

wipe
model BasicBuilder -ndm 1 -ndf 1

set elementsPerRank 4
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

set values [eigen -ladrunoCMS \
    -domainMode physical \
    -hierarchy logical -level1 2 -level2 2 \
    -modesL2 3 -modesL1 6 \
    -tol 1.0e-8 -maxEnrich 2 -maxIter 300 \
    -iterationVectors 6 -maxRefineIter 80 \
    -denseMax 128 -verbose 4]

for {set mode 1} {$mode <= 4} {incr mode} {
    set expected [expr {2.0 - 2.0*cos((2*$mode-1)*acos(-1.0)/33.0)}]
    set actual [lindex $values [expr {$mode - 1}]]
    set relativeError [expr {abs($actual-$expected)/max(abs($expected),1.0e-30)}]
    if {$relativeError > 1.0e-8} {
        error "physical CMS mode $mode relative error $relativeError"
    }
}

puts "LadrunoCMS physical rank=$pid nodes=[expr {$lastNode-$firstNode+1}] elements=$elementsPerRank"
barrier
if {$pid == 0} {
    puts "ladruno CMS physical-domain four-rank smoke passed"
}
wipe
