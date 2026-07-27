# ADR-1000 P4 section 2 -- CMS against the standard solver, same model, same box.
#
# THIS IS NOT THE BUILDING-1A COMPARISON. That deck is not in the repository, so
# the 7x-slower-than-ARPACK verdict recorded in the P4 plan cannot be re-measured
# here. This is a like-for-like comparison on the sheet model instead: useful for
# tracking whether optimisation work moves the ratio, not a substitute for the
# real gate.
#
# The SAME model is built either way -- a quad sheet of `4*width` columns:
#   mpiexec -n 1 OpenSeesMP.exe cms_compare_arpack.tcl   -> monolithic + standard
#                                                           `eigen` (Arpack+UmfPack)
#   mpiexec -n 4 OpenSeesMP.exe cms_compare_arpack.tcl   -> partitioned + CMS
# so any difference is the solver, not the mesh.
#
#   LADRUNO_CMS_PROFILE_WIDTH   columns per strip (default 20); the sheet always
#                               has 4 strips regardless of rank count
#   LADRUNO_CMS_PROFILE_HEIGHT  rows of nodes (default 20)
#   LADRUNO_CMS_COMPARE_MODES   modes requested (default 6)
#   LADRUNO_CMS_COMPARE_MODESL2 fixed-interface modes k2 kept per subdomain
#                               (default 12 -- the section 30 setting)
#   LADRUNO_CMS_COMPARE_MODESL1 level-1 modes k1 (default 2*k2)
#
# section 31: k2 was held at 12 at every mesh size in the section 30 ladder, so
# that ladder could not tell a tuning artefact from an algorithmic one. These
# two knobs exist to scale k2 with the subdomain (k2 ~ sqrt(m)) and re-run.
set pid [getPID]
set np [getNP]
if {$np != 1 && $np != 4} {
    error "CMS/ARPACK comparison runs at 1 rank (standard) or 4 ranks (CMS)"
}

set width 20
if {[info exists env(LADRUNO_CMS_PROFILE_WIDTH)]} {
    set width $env(LADRUNO_CMS_PROFILE_WIDTH)
}
set height 20
if {[info exists env(LADRUNO_CMS_PROFILE_HEIGHT)]} {
    set height $env(LADRUNO_CMS_PROFILE_HEIGHT)
}
set numModes 6
if {[info exists env(LADRUNO_CMS_COMPARE_MODES)]} {
    set numModes $env(LADRUNO_CMS_COMPARE_MODES)
}
set modesL2 12
if {[info exists env(LADRUNO_CMS_COMPARE_MODESL2)]} {
    set modesL2 $env(LADRUNO_CMS_COMPARE_MODESL2)
}
set modesL1 [expr {2 * $modesL2}]
if {[info exists env(LADRUNO_CMS_COMPARE_MODESL1)]} {
    set modesL1 $env(LADRUNO_CMS_COMPARE_MODESL1)
}

# The sheet is ALWAYS 4 strips wide, so the two runs solve the same problem.
set strips 4
wipe
model BasicBuilder -ndm 2 -ndf 2
proc nodeTag {i j height} { return [expr {$i * $height + $j + 1}] }

if {$np == 1} {
    set firstColumn 0
    set lastColumn [expr {$strips * $width}]
} else {
    set firstColumn [expr {$pid * $width}]
    set lastColumn [expr {$firstColumn + $width}]
}

for {set i $firstColumn} {$i <= $lastColumn} {incr i} {
    for {set j 0} {$j < $height} {incr j} {
        node [nodeTag $i $j $height] [expr {double($i)}] [expr {double($j)}]
    }
}
if {$pid == 0} {
    for {set j 0} {$j < $height} {incr j} { fix [nodeTag 0 $j $height] 1 1 }
}

nDMaterial ElasticIsotropic 1 1000.0 0.25
set element 1
set elementBase [expr {$pid * 1000000}]
for {set i $firstColumn} {$i < $lastColumn} {incr i} {
    for {set j 0} {$j < $height - 1} {incr j} {
        element quad [expr {$elementBase + $element}] \
            [nodeTag $i $j $height] [nodeTag [expr {$i + 1}] $j $height] \
            [nodeTag [expr {$i + 1}] [expr {$j + 1}] $height] \
            [nodeTag $i [expr {$j + 1}] $height] \
            1.0 "PlaneStress" 1
        incr element
    }
}

set massFirst [expr {$pid == 0 ? $firstColumn : $firstColumn + 1}]
for {set i $massFirst} {$i <= $lastColumn} {incr i} {
    for {set j 0} {$j < $height} {incr j} {
        if {$pid == 0 && $i == 0} { continue }
        mass [nodeTag $i $j $height] 1.0 1.0
    }
}

constraints Plain
test NormUnbalance 1.0e-12 2
algorithm Linear
integrator LoadControl 1.0

if {$np == 1} {
    # The ADR's bar: "the standard OpenSees solver".
    numberer RCM
    system UmfPack
    analysis Static
    set started [clock microseconds]
    set values [eigen $numModes]
    set elapsed [expr {([clock microseconds] - $started) / 1.0e6}]
    puts "COMPARE standard np=1 modes=$numModes seconds=$elapsed"
} else {
    numberer ParallelRCM
    system Mumps
    analysis Static
    set started [clock microseconds]
    set values [eigen -ladrunoCMS \
        -domainMode physical \
        -hierarchy logical -level1 2 -level2 2 \
        -modesL2 $modesL2 -modesL1 $modesL1 \
        -tol 1.0e-8 -maxEnrich 2 -maxIter 200000 \
        -maxRestarts 4000 -maxRefineIter 160 \
        -denseMax 4000 $numModes]
    set elapsed [expr {([clock microseconds] - $started) / 1.0e6}]
    if {$pid == 0} {
        puts "COMPARE cms np=4 modes=$numModes k2=$modesL2 k1=$modesL1 seconds=$elapsed"
    }
}

if {$pid == 0} {
    for {set mode 0} {$mode < $numModes} {incr mode} {
        puts [format "COMPARE eigenvalue %d %.14e" $mode [lindex $values $mode]]
    }
}
