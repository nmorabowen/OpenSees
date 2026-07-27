# Ladruno ADR-1000 §31.5 smoke: plain `eigen` (ARPACK shift-invert) over a
# PARTITIONED deck + MumpsParallelSOE under OpenSeesMP — the composition
# apeGmsh ADR 0077 F1 refuted in vanilla (ArpackSOE left at processID -1,
# M*v never merged across ranks). PASS requires the wired build.
#
# Fixed-free chain, nMass masses, k=100, m=1:
#   lambda_j = 4*(k/m)*sin^2((2j-1)*pi/(2*(2*nMass+1)))
#
# Partition discipline: element e -> rank (e-1)*np/nel (contiguous blocks);
# boundary nodes are defined on both touching ranks; a node's MASS is owned
# by the rank of its left element only (double-defined nodal mass would be
# summed twice by the M*v merge — same owner rule every MP static deck uses
# for nodal loads).
#
#   serial oracle:  OpenSees   verify_arpack_mp_mumps.tcl
#   distributed:    mpiexec -n 2 OpenSeesMP verify_arpack_mp_mumps.tcl

if {[llength [info commands getPID]] == 0} {
    proc getPID {} { return 0 }
    proc getNP  {} { return 1 }
}
set pid [getPID]
set np  [getNP]

set nMass 8
set nel   $nMass
set k 100.0
set m 1.0
set nEig 4

proc eleRank {e np nel} { expr {($e - 1) * $np / $nel} }

wipe
model basic -ndm 1 -ndf 1
uniaxialMaterial Elastic 1 $k

# contiguous element range of this rank
set eLo 0
set eHi -1
for {set e 1} {$e <= $nel} {incr e} {
    if {[eleRank $e $np $nel] == $pid} {
        if {$eLo == 0} { set eLo $e }
        set eHi $e
    }
}
if {$eLo == 0} { puts "ARPACK_MP_SMOKE_FAIL rank=$pid: no elements (np > nel?)" ; exit 1 }

for {set j $eLo} {$j <= [expr {$eHi + 1}]} {incr j} {
    node $j [expr {double($j - 1)}]
}
for {set e $eLo} {$e <= $eHi} {incr e} {
    element truss $e $e [expr {$e + 1}] 1.0 1
}
if {$eLo == 1} { fix 1 1 }

# nodal mass of node j owned by the rank of element j-1
for {set j 2} {$j <= [expr {$nMass + 1}]} {incr j} {
    if {[eleRank [expr {$j - 1}] $np $nel] == $pid} {
        mass $j $m
    }
}

constraints Transformation
if {$np > 1} {
    numberer ParallelPlain
    system Mumps
} else {
    numberer RCM
    system UmfPack
}

set lam [eigen $nEig]

# analytic oracle
set pi [expr {acos(-1.0)}]
set maxRelErr -1.0
for {set j 1} {$j <= $nEig} {incr j} {
    set lamA [expr {4.0 * $k / $m * pow(sin((2*$j - 1) * $pi / (2.0 * (2*$nMass + 1))), 2)}]
    set lamN [lindex $lam [expr {$j - 1}]]
    if {$lamN == ""} { set rel 1.0e99 } else {
        set rel [expr {abs($lamN - $lamA) / $lamA}]
    }
    if {$rel > $maxRelErr} { set maxRelErr $rel }
}

# mode-1 shape check: phi_1(node i) ~ sin((2*1-1)*pi*(i-1)/(2*nMass+1)).
# Ratios to a local reference node cancel the arbitrary normalization; a
# rank-local (unmerged-M) Lanczos would give garbage ratios here even if
# some eigenvalue accidentally matched.
set refNode [expr {$eLo == 1 ? 2 : $eLo}]
set refA [expr {sin($pi * ($refNode - 1) / (2.0 * $nMass + 1))}]
set refN [lindex [nodeEigenvector $refNode 1] 0]
set maxShapeErr -1.0
for {set j $eLo} {$j <= [expr {$eHi + 1}]} {incr j} {
    if {$j == 1} { continue }
    set phiA [expr {sin($pi * ($j - 1) / (2.0 * $nMass + 1)) / $refA}]
    set phiN [expr {[lindex [nodeEigenvector $j 1] 0] / $refN}]
    set err [expr {abs($phiN - $phiA)}]
    if {$err > $maxShapeErr} { set maxShapeErr $err }
}

set fp [open arpack_mp_mumps_rank$pid.out w]
puts $fp "np=$np lambda=$lam"
puts $fp "maxRelErr=$maxRelErr"
puts $fp "maxShapeErr=$maxShapeErr"
close $fp

if {$maxRelErr > 1.0e-8 || $maxRelErr < 0.0 || $maxShapeErr > 1.0e-6 || $maxShapeErr < 0.0} {
    puts "ARPACK_MP_SMOKE_FAIL rank=$pid maxRelErr=$maxRelErr maxShapeErr=$maxShapeErr lambda=$lam"
    exit 1
}
puts "ARPACK_MP_SMOKE_PASS rank=$pid np=$np maxRelErr=$maxRelErr maxShapeErr=$maxShapeErr"
