# ADR-78 P4 -- mortar MESH-TIE across ranks (implicit + ALM).
#
# The ADR-41 C4.1 split column, partitioned: a column [0,1]^2 x [0,2] split at
# z=1 into a bottom block (ONE LadrunoBrick, its top face = the MASTER facet)
# and a top block (TWO bricks split in x at 0.5 -> a NON-MATCHING slave face
# with a shared mid-edge node). The tie welds the interface; the column is
# loaded in TENSION (an equality bond holds -- contact would separate), so a
# silently-missing tie cannot masquerade as a pass.
#
#   np=1: everything in one domain (the serial twin).
#   np=2: rank0 = bottom block + BOTH contactSurface defs + the -mortar -tie
#         verb + GHOSTS of the slave interface nodes 9..14 (node + replayed
#         fix, no mass -- INV-6); rank1 = top block + the end load.
#
# ALM: after the loaded solve, a HELD-LOAD Uzawa augmentation sweep (ADR-41 D1,
# ladrunoBeginAugment/ladrunoEndAugment) drives the tie residual to ~0 at
# finite epsTie. The sweep runs a FIXED count of augmentations on EVERY rank:
# analyze() is collective (Mumps), so a data-dependent break keyed on the
# rank-local residual would desynchronise the ranks and deadlock. The residual
# gate is applied AFTER the sweep, on the owner rank.
#
# Analytic (E=2e4, nu=0, A=1, PTOT=20 tension): eps = 1e-3,
#   u(z=1) = +1.0e-3 (interface), u(z=2) = +2.0e-3 (tip).
#
# env: P4_TAG output tag (default np$np)

wipe
model basic -ndm 3 -ndf 3
set pid [getPID]
set np  [getNP]
if {[info exists env(P4_TAG)]} { set TAG $env(P4_TAG) } else { set TAG "np$np" }

set E      2.0e4
set PTOT   20.0
set NAUG   15
set EPSTIE 1.0e8

nDMaterial ElasticIsotropic 1 $E 0.0
timeSeries Linear 1

set haveBot [expr {$np == 1 || $pid == 0}]
set haveTop [expr {$np == 1 || $pid == 1}]

# slave-interface coordinates (top block, z=1): 9..14 = x{0,0.5,1} x y{0,1}
set SLAVEXY {9 {0.0 0.0} 10 {0.5 0.0} 11 {1.0 0.0} 12 {0.0 1.0} 13 {0.5 1.0} 14 {1.0 1.0}}

if {$haveBot} {
    # bottom block (MASTER): 1..4 z=0 fixed, 5..8 z=1 (master facet)
    node 1 0.0 0.0 0.0
    node 2 1.0 0.0 0.0
    node 3 1.0 1.0 0.0
    node 4 0.0 1.0 0.0
    node 5 0.0 0.0 1.0
    node 6 1.0 0.0 1.0
    node 7 1.0 1.0 1.0
    node 8 0.0 1.0 1.0
    foreach n {1 2 3 4} { fix $n 1 1 1 }
    foreach n {5 6 7 8} { fix $n 1 1 0 }
    element LadrunoBrick 1 1 2 3 4 5 6 7 8 1
}

if {$np > 1 && $pid == 0} {
    # GHOSTS of the slave interface (owned by rank 1): node + replayed fix,
    # no mass, no elements (INV-6).
    foreach {n xy} $SLAVEXY {
        lassign $xy x y
        node $n $x $y 1.0
        fix $n 1 1 0
    }
}

if {$haveTop} {
    # top block (SLAVE): interface z=1 nodes 9..14, top z=2 nodes 15..20
    foreach {n xy} $SLAVEXY {
        lassign $xy x y
        node $n $x $y 1.0
        fix $n 1 1 0
    }
    foreach {n xy} {15 {0.0 0.0} 16 {0.5 0.0} 17 {1.0 0.0} 18 {0.0 1.0} 19 {0.5 1.0} 20 {1.0 1.0}} {
        lassign $xy x y
        node $n $x $y 2.0
        fix $n 1 1 0
    }
    # brick A [0,0.5] x [0,1] x [1,2]; brick B [0.5,1] x [0,1] x [1,2]
    element LadrunoBrick 2  9 10 13 12 15 16 19 18 1
    element LadrunoBrick 3 10 11 14 13 16 17 20 19 1
}

if {$haveBot} {
    # MASTER facet (nps=4, face 5 6 7 8); SLAVE = the two non-matching facets.
    contactSurface 1 -master 4 5 6 7 8
    contactSurface 2 -slave-segments 4 9 10 13 12 10 11 14 13
    # -tie: a permanent bond; the interface is coincident at z=1 so the auto
    # orientation degenerates -- pass -outward explicitly (as ADR-41 C2/C3 do).
    contact 1 1 2 -mortar -tie -epsTie $EPSTIE -outward 0.0 0.0 1.0
}

pattern Plain 1 1 { }
if {$haveTop} {
    # consistent nodal loads for a uniform TENSION traction over the 2-quad
    # top face: corners t/8, shared mid-edge nodes t/4.
    pattern Plain 2 1 {
        load 15 0.0 0.0 [expr {$PTOT / 8.0}]
        load 17 0.0 0.0 [expr {$PTOT / 8.0}]
        load 18 0.0 0.0 [expr {$PTOT / 8.0}]
        load 20 0.0 0.0 [expr {$PTOT / 8.0}]
        load 16 0.0 0.0 [expr {$PTOT / 4.0}]
        load 19 0.0 0.0 [expr {$PTOT / 4.0}]
    }
}

constraints LadrunoContact
if {$np > 1} {
    numberer ParallelPlain
    system Mumps
} else {
    numberer RCM
    system UmfPack
}
test NormDispIncr 1.0e-12 60 0
algorithm Newton
integrator LoadControl 1.0
analysis Static

set ok [analyze 1]

# ---- held-load ALM augmentation (ADR-41 D1), fixed count on every rank -----
set okA 0
integrator LoadControl 0.0
ladrunoBeginAugment
for {set k 0} {$k < $NAUG} {incr k} {
    set rc [analyze 1]
    if {$rc != 0} { set okA $rc }
}
ladrunoEndAugment

if {$haveBot} {
    puts [format "P4TIEB %s pid=%d ok=%d okA=%d res=%.6e u5=%.16e u9g=%.16e" \
          $TAG $pid $ok $okA [ladrunoMortarTieResidual] [nodeDisp 5 3] [nodeDisp 9 3]]
}
if {$haveTop} {
    puts [format "P4TIET %s pid=%d ok=%d okA=%d u9=%.16e u15=%.16e" \
          $TAG $pid $ok $okA [nodeDisp 9 3] [nodeDisp 15 3]]
}
wipe
