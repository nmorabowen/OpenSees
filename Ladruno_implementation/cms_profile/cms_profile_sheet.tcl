# ADR-1000 P4 section 4 -- profiling deck with a REAL interface.
#
# cms_profile_chain.tcl is a 1-D chain: every rank boundary is ONE dof, so b~1
# and the b-dependent costs inside T2 -- the constraint-mode solve (O(b) right-
# hand sides) and the congruence (O(n*(k+b)^2)) -- are invisible. Those are
# exactly the costs that dominate on a real partitioned model, where an
# interface is a whole cut surface.
#
# This is a 2-D sheet of truss members cut into four vertical strips, so each
# rank boundary is a COLUMN of nodes: b = 2 * height, tunable independently of
# the subdomain size.
#
# Run:  mpiexec -n 4 OpenSeesMP.exe cms_profile_sheet.tcl
#   LADRUNO_CMS_PROFILE_WIDTH   columns of nodes per rank (default 40)
#   LADRUNO_CMS_PROFILE_HEIGHT  rows of nodes; interface b = 2 * height (default 40)
#   LADRUNO_CMS_PROFILE_RESTARTS / _MAXITER  T2 budgets (see section 27.4)
set pid [getPID]
set np [getNP]
if {$np != 4} {
    error "CMS sheet profile requires exactly four MPI ranks"
}

set width 40
if {[info exists env(LADRUNO_CMS_PROFILE_WIDTH)]} {
    set width $env(LADRUNO_CMS_PROFILE_WIDTH)
}
set height 40
if {[info exists env(LADRUNO_CMS_PROFILE_HEIGHT)]} {
    set height $env(LADRUNO_CMS_PROFILE_HEIGHT)
}
set maxRestarts 200
if {[info exists env(LADRUNO_CMS_PROFILE_RESTARTS)]} {
    set maxRestarts $env(LADRUNO_CMS_PROFILE_RESTARTS)
}
set maxIter 20000
if {[info exists env(LADRUNO_CMS_PROFILE_MAXITER)]} {
    set maxIter $env(LADRUNO_CMS_PROFILE_MAXITER)
}

wipe
model BasicBuilder -ndm 2 -ndf 2

# Global node numbering: column-major, node(i,j) = i*height + j + 1, with
# i in [0, np*width] (one extra column so the strips share a boundary column).
set columns [expr {$np * $width + 1}]
proc nodeTag {i j height} { return [expr {$i * $height + $j + 1}] }

# This rank owns columns [pid*width, (pid+1)*width]; the last column is shared
# with the next rank -- that shared COLUMN is the interface, b = height.
set firstColumn [expr {$pid * $width}]
set lastColumn [expr {$firstColumn + $width}]

for {set i $firstColumn} {$i <= $lastColumn} {incr i} {
    for {set j 0} {$j < $height} {incr j} {
        node [nodeTag $i $j $height] [expr {double($i)}] [expr {double($j)}]
    }
}

# Fix the first column of the whole sheet, on the rank that owns it.
if {$pid == 0} {
    for {set j 0} {$j < $height} {incr j} { fix [nodeTag 0 $j $height] 1 1 }
}

# Continuum quads, not trusses. A truss grid carries axial force only, so a
# rectangular one is a MECHANISM (zero-energy shear modes) and K_II comes out
# singular -- MUMPS INFOG(1) = -10. Bracing it helps but leaves the last column
# of the last rank under-connected. A quad sheet has no such failure mode.
nDMaterial ElasticIsotropic 1 1000.0 0.25
set element 1
set elementBase [expr {$pid * 1000000}]
for {set i $firstColumn} {$i < $lastColumn} {incr i} {
    for {set j 0} {$j < $height - 1} {incr j} {
        element quad [expr {$elementBase + $element}]             [nodeTag $i $j $height] [nodeTag [expr {$i + 1}] $j $height]             [nodeTag [expr {$i + 1}] [expr {$j + 1}] $height]             [nodeTag $i [expr {$j + 1}] $height]             1.0 "PlaneStress" 1
        incr element
    }
}

# Mass on every node this rank owns PRIMARILY: all its columns except the shared
# boundary column, which belongs to the lower-numbered rank... which is this one,
# so the shared column is ours and the FIRST column is the previous rank's.
set massFirst [expr {$pid == 0 ? $firstColumn : $firstColumn + 1}]
for {set i $massFirst} {$i <= $lastColumn} {incr i} {
    for {set j 0} {$j < $height} {incr j} {
        if {$pid == 0 && $i == 0} { continue }
        mass [nodeTag $i $j $height] 1.0 1.0
    }
}

constraints Plain
numberer ParallelRCM
system Mumps
test NormUnbalance 1.0e-12 2
algorithm Linear
integrator LoadControl 1.0
analysis Static

if {$pid == 0} {
    puts "CMS sheet: $np ranks, ${width}x${height} per rank, interface b=[expr {2*$height}]"
}

set values [eigen -ladrunoCMS \
    -domainMode physical \
    -hierarchy logical -level1 2 -level2 2 \
    -modesL2 12 -modesL1 24 \
    -tol 1.0e-8 -maxEnrich 2 -maxIter $maxIter \
    -maxRestarts $maxRestarts \
    -maxRefineIter 160 \
    -denseMax 4000 -verbose 6]

if {$pid == 0} {
    puts "CMS sheet: first eigenvalue [lindex $values 0]"
}
