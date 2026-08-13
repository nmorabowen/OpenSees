# ADR-78 P4 -- two-body POUNDING deck (explicit, NTS, -kn auto).
#
# One deck, three partitionings, selected by [getNP] -- np=1, 2, 4 all run the
# SAME physics, in the manual-DD shape apeGmsh's partitioned emit produces
# (one plain Domain per rank inside getPID guards, boundary nodes replicated
# by tag, ghosts = node + replayed fix, NEVER mass -- ADR-78 INV-6).
#
# Model: two 1x1xNB brick columns stacked with an initial gap. The bottom
# body is fixed at its base (MASTER, top face = master segment); the top
# body (SLAVE, bottom nodes = slave set) carries a constant downward load and
# falls onto it -> repeated impact / bounce / separation (pounding). Lateral
# DOFs fixed everywhere => a clean 1-D column.
#
# Partition map ("contact owner" = the rank holding the master face's
# backing solid, uncut -- ADR-78 D1):
#   np=1: rank0 = everything.
#   np=2: rank0 = bottom body + contact + ghost slaves; rank1 = top body + load.
#   np=4: rank0 = lower half bottom body (base);
#         rank1 = upper half bottom body + contact + ghost slaves;
#         rank2 = lower half top body; rank3 = upper half top body + load.
#         Cut-level nodes are replicated on both neighbour ranks; their nodal
#         MASS is declared by the HIGHER rank only (the distributed diagonal
#         solver sums A over shared DOFs, so declaring it twice would
#         double-count -- same mechanism as INV-6).
#
# env knobs (all optional):
#   P4_TAG       output-file tag                      (default np$np)
#   P4_NB        bricks per body, even                (default 2)
#   P4_NSTEP     explicit steps                       (default 3000)
#   P4_SOFT      add "-soft 0.5" to the contact verb  -> P2 partition-refusal control
#   P4_NOGHOST   owner omits the ghost slave nodes    -> P1 FATAL control
#   P4_GHOSTMASS ghosts get nodal mass (INV-6 mutation) -> silent wrong answer
#                the serial-comparison instrument must catch

wipe
model basic -ndm 3 -ndf 3
set pid [getPID]
set np  [getNP]

set NB    2
set NSTEP 3000
if {[info exists env(P4_NB)]}    { set NB    [expr {int($env(P4_NB))}] }
if {[info exists env(P4_NSTEP)]} { set NSTEP [expr {int($env(P4_NSTEP))}] }
if {[info exists env(P4_TAG)]}   { set TAG $env(P4_TAG) } else { set TAG "np$np" }
set DT    2.0e-5
set GAP   5.0e-3
set MNODE 1.0
set PNODE -2500.0
set half  [expr {$NB / 2}]
if {$np == 4 && 2 * $half != $NB} { puts "P4POUND ERROR NB must be even for np=4"; exit 1 }

proc seq {a b} { set r {}; for {set i $a} {$i <= $b} {incr i} { lappend r $i }; return $r }
proc bnode {i j} { expr {100 + 4 * $i + $j} }   ;# bottom body, level i (z=i), corner j=1..4
proc tnode {i j} { expr {200 + 4 * $i + $j} }   ;# top body, level i (z=NB+GAP+i)
set XY {{0.0 0.0} {1.0 0.0} {1.0 1.0} {0.0 1.0}}

nDMaterial ElasticIsotropic 1 2.0e7 0.0
timeSeries Constant 1

# ---- per-rank declaration sets --------------------------------------------
set botLevs {};  set botMassLevs {};  set botBricks {}
set topLevs {};  set topMassLevs {};  set topBricks {}
set hasContact 0; set hasGhost 0; set hasLoad 0
set hasM1 0; set hasS1 0; set hasT1 0; set hasBase 0

if {$np == 1} {
    set botLevs [seq 0 $NB]; set botMassLevs $botLevs; set botBricks [seq 0 [expr {$NB-1}]]
    set topLevs [seq 0 $NB]; set topMassLevs $topLevs; set topBricks [seq 0 [expr {$NB-1}]]
    set hasContact 1; set hasLoad 1
    set hasM1 1; set hasS1 1; set hasT1 1; set hasBase 1
} elseif {$np == 2} {
    if {$pid == 0} {
        set botLevs [seq 0 $NB]; set botMassLevs $botLevs; set botBricks [seq 0 [expr {$NB-1}]]
        set hasContact 1; set hasGhost 1; set hasM1 1; set hasBase 1
    } else {
        set topLevs [seq 0 $NB]; set topMassLevs $topLevs; set topBricks [seq 0 [expr {$NB-1}]]
        set hasLoad 1; set hasS1 1; set hasT1 1
    }
} elseif {$np == 4} {
    if {$pid == 0} {
        set botLevs [seq 0 $half]; set botMassLevs [seq 0 [expr {$half-1}]]
        set botBricks [seq 0 [expr {$half-1}]]; set hasBase 1
    } elseif {$pid == 1} {
        set botLevs [seq $half $NB]; set botMassLevs $botLevs
        set botBricks [seq $half [expr {$NB-1}]]
        set hasContact 1; set hasGhost 1; set hasM1 1
    } elseif {$pid == 2} {
        set topLevs [seq 0 $half]; set topMassLevs [seq 0 [expr {$half-1}]]
        set topBricks [seq 0 [expr {$half-1}]]; set hasS1 1
    } else {
        set topLevs [seq $half $NB]; set topMassLevs $topLevs
        set topBricks [seq $half [expr {$NB-1}]]
        set hasLoad 1; set hasT1 1
    }
} else {
    puts "P4POUND ERROR unsupported np=$np"; exit 1
}

# ---- nodes / fixes / masses ------------------------------------------------
foreach i $botLevs {
    set z [expr {double($i)}]
    foreach j {1 2 3 4} xy $XY {
        lassign $xy x y
        node [bnode $i $j] $x $y $z
        if {$i == 0} { fix [bnode $i $j] 1 1 1 } else { fix [bnode $i $j] 1 1 0 }
        if {[lsearch -exact $botMassLevs $i] >= 0} { mass [bnode $i $j] $MNODE $MNODE $MNODE }
    }
}
foreach i $topLevs {
    set z [expr {$NB + $GAP + double($i)}]
    foreach j {1 2 3 4} xy $XY {
        lassign $xy x y
        node [tnode $i $j] $x $y $z
        fix [tnode $i $j] 1 1 0
        if {[lsearch -exact $topMassLevs $i] >= 0} { mass [tnode $i $j] $MNODE $MNODE $MNODE }
    }
}
foreach i $botBricks {
    element stdBrick [expr {100 + $i}] \
        [bnode $i 1] [bnode $i 2] [bnode $i 3] [bnode $i 4] \
        [bnode [expr {$i+1}] 1] [bnode [expr {$i+1}] 2] [bnode [expr {$i+1}] 3] [bnode [expr {$i+1}] 4] 1
}
foreach i $topBricks {
    element stdBrick [expr {200 + $i}] \
        [tnode $i 1] [tnode $i 2] [tnode $i 3] [tnode $i 4] \
        [tnode [expr {$i+1}] 1] [tnode [expr {$i+1}] 2] [tnode [expr {$i+1}] 3] [tnode [expr {$i+1}] 4] 1
}

# ---- ghosts of the slave interface on the contact owner (INV-6: no mass) ---
if {$hasGhost && ![info exists env(P4_NOGHOST)]} {
    set z [expr {$NB + $GAP}]
    foreach j {1 2 3 4} xy $XY {
        lassign $xy x y
        node [tnode 0 $j] $x $y $z
        fix [tnode 0 $j] 1 1 0
        if {[info exists env(P4_GHOSTMASS)]} { mass [tnode 0 $j] $MNODE $MNODE $MNODE }
    }
}

# ---- contact (owner rank only -- INV-1: emitted exactly once) --------------
if {$hasContact} {
    contactSurface 1 -master 4 [bnode $NB 1] [bnode $NB 2] [bnode $NB 3] [bnode $NB 4]
    contactSurface 2 -slave [tnode 0 1] [tnode 0 2] [tnode 0 3] [tnode 0 4]
    if {[info exists env(P4_SOFT)]} {
        contact 1 1 2 auto -outward 0.0 0.0 1.0 -soft 0.5
    } else {
        contact 1 1 2 auto -outward 0.0 0.0 1.0
    }
}

# ---- load (constant "gravity" pull on the top face of the top body) --------
pattern Plain 1 1 { }
if {$hasLoad} {
    pattern Plain 2 1 {
        load [tnode $NB 1] 0.0 0.0 $PNODE
        load [tnode $NB 2] 0.0 0.0 $PNODE
        load [tnode $NB 3] 0.0 0.0 $PNODE
        load [tnode $NB 4] 0.0 0.0 $PNODE
    }
}

# ---- analysis chain (P0 explicit-lane shape) --------------------------------
constraints LadrunoContact
if {$np > 1} {
    numberer ParallelPlain
    system MPIDiagonal
} else {
    numberer RCM
    system Diagonal
}
integrator CentralDifferenceLadruno
algorithm Linear

# ---- recorders (per rank; sample nodes are recorded by their NATIVE rank) --
if {$hasM1} {
    recorder Node -file disp_${TAG}_m1.txt -precision 17 -time -node [bnode $NB 1] -dof 3 disp
    recorder Node -file vel_${TAG}_m1.txt -precision 17 -time -node [bnode $NB 1] -dof 3 vel
}
if {$hasS1} {
    recorder Node -file disp_${TAG}_s1.txt -precision 17 -time -node [tnode 0 1] -dof 3 disp
    recorder Node -file vel_${TAG}_s1.txt -precision 17 -time -node [tnode 0 1] -dof 3 vel
}
if {$hasT1} {
    recorder Node -file disp_${TAG}_t1.txt -precision 17 -time -node [tnode $NB 1] -dof 3 disp
    recorder Node -file vel_${TAG}_t1.txt -precision 17 -time -node [tnode $NB 1] -dof 3 vel
}
if {$hasGhost && ![info exists env(P4_NOGHOST)]} {
    recorder Node -file disp_${TAG}_s1g.txt -precision 17 -time -node [tnode 0 1] -dof 3 disp
}
if {$hasBase} {
    recorder Node -file reac_${TAG}.txt -precision 17 -time -node 101 102 103 104 -dof 3 reaction
}
recorder EnergyBalance -file energy_${TAG}_r${pid}.txt -precision 17 -time

analysis Transient
set ok [analyze $NSTEP $DT]

# ---- finals ------------------------------------------------------------------
if {$hasM1} {
    puts [format "P4VAL %s pid=%d m1=%.16e" $TAG $pid [nodeDisp [bnode $NB 1] 3]]
    if {$hasGhost && ![info exists env(P4_NOGHOST)]} {
        puts [format "P4VAL %s pid=%d s1g=%.16e" $TAG $pid [nodeDisp [tnode 0 1] 3]]
    }
}
if {$hasS1} {
    puts [format "P4VAL %s pid=%d s1=%.16e vs1=%.16e" $TAG $pid [nodeDisp [tnode 0 1] 3] [nodeVel [tnode 0 1] 3]]
}
if {$hasT1} {
    puts [format "P4VAL %s pid=%d t1=%.16e vt1=%.16e" $TAG $pid [nodeDisp [tnode $NB 1] 3] [nodeVel [tnode $NB 1] 3]]
}
puts [format "P4POUND %s pid=%d np=%d analyze=%d" $TAG $pid $np $ok]
wipe
