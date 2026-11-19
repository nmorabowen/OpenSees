# WP-240: LadrunoDynamicRelaxation in OpenSeesMP -- a bar of 8 LadrunoBrick
# elements along x, split by x over the processes (element e on process
# int(e*np/8)); the interface nodes exist on both sides with the same tag.
# Run by OpenSees.exe (serial reference) and by OpenSeesMP.exe on 1, 2 and 4
# processes; driven by tests/test_wp240_dr_mp.py, which compares the files.
#
# Environment:
#   WP240_OUT      output prefix (required)
#   WP240_MODE     hist  : 400 steps past yield, residual / KE / max|u| per step
#                  relax : `ladrunoDR relax` to rest below yield
#   WP240_SYSTEM   Diagonal | MPIDiagonal   (default: by engine)
#   WP240_REFRESH  1 = leave the KE-peak refresh of M* on (default 0: -noAutoRefresh)
#
# Files: <out>.hist<pid>  "step residualNorm kineticEnergy maxDisplacement"
#        <out>.node<pid>  "tag ux uy uz" of the nodes this process holds
#        <out>.rc<pid>    "rc <analyze result>" and, in relax mode, the result list

set out $::env(WP240_OUT)
set mode hist
if {[info exists ::env(WP240_MODE)]} { set mode $::env(WP240_MODE) }
set refresh 0
if {[info exists ::env(WP240_REFRESH)]} { set refresh $::env(WP240_REFRESH) }

set parallel [expr {[llength [info commands getNP]] > 0}]
if {$parallel} { set pid [getPID]; set np [getNP] } else { set pid 0; set np 1 }
set sys [expr {$parallel ? "MPIDiagonal" : "Diagonal"}]
if {[info exists ::env(WP240_SYSTEM)]} { set sys $::env(WP240_SYSTEM) }

set NEL 8
set E 2.0e5; set nu 0.3; set sig0 250.0
set K [expr {$E / (3.0 * (1.0 - 2.0 * $nu))}]
set G [expr {$E / (2.0 * (1.0 + $nu))}]
# hist: 1.1 x the yield force of the 1 x 1 section; relax: half of it
set P [expr {$mode eq "relax" ? 0.5 * $sig0 : 1.1 * $sig0}]

proc ntag {i j k} { return [expr {1 + 4 * $i + 2 * $j + $k}] }

model basic -ndm 3 -ndf 3
nDMaterial LadrunoJ2 1 $K $G -iso voce $sig0 0.0 0.0 200.0

set mine {}
for {set e 0} {$e < $NEL} {incr e} {
    if {int($e * $np / $NEL) == $pid} { lappend mine $e }
}
set held {}
foreach e $mine {
    foreach i [list $e [expr {$e + 1}]] {
        if {[lsearch -exact $held $i] >= 0} { continue }
        lappend held $i
        foreach j {0 1} { foreach k {0 1} {
            node [ntag $i $j $k] [expr {double($i)}] [expr {double($j)}] [expr {double($k)}]
        } }
    }
}
foreach e $mine {
    set a $e; set b [expr {$e + 1}]
    element LadrunoBrick [expr {$e + 1}] \
        [ntag $a 0 0] [ntag $b 0 0] [ntag $b 1 0] [ntag $a 1 0] \
        [ntag $a 0 1] [ntag $b 0 1] [ntag $b 1 1] [ntag $a 1 1] 1
}
if {[lsearch -exact $held 0] >= 0} {
    foreach j {0 1} { foreach k {0 1} { fix [ntag 0 $j $k] 1 1 1 } }
}

# A load on a node is given by ONE process: the one holding the element on its
# left (the tip by the last element, the mid-span interface by element 3).
timeSeries Constant 1
pattern Plain 1 1 {
    if {[lsearch -exact $::mine 7] >= 0} {
        foreach j {0 1} { foreach k {0 1} { load [ntag 8 $j $k] [expr {$::P / 4.0}] 0.0 0.0 } }
    }
    if {[lsearch -exact $::mine 3] >= 0} {
        foreach j {0 1} { foreach k {0 1} { load [ntag 4 $j $k] 0.0 0.0 [expr {-0.01 * $::P / 4.0}] } }
    }
}

constraints Transformation
if {$parallel} { numberer ParallelPlain } else { numberer Plain }
system $sys
test NormUnbalance 1.0e30 1 0
algorithm Linear
if {$refresh} {
    integrator LadrunoDynamicRelaxation
} else {
    integrator LadrunoDynamicRelaxation -noAutoRefresh
}
analysis Transient

set rc 0
set result {}
if {$mode eq "relax"} {
    if {[catch {ladrunoDR relax -tolResidual [expr {1.0e-6 * $P}] -maxSteps 20000} result]} {
        set rc -1
    }
} else {
    set fh [open "${out}.hist${pid}" w]
    for {set s 1} {$s <= 400} {incr s} {
        set rc [analyze 1 1.0]
        if {$rc != 0} { break }
        puts $fh [format "%d %.17g %.17g %.17g" $s [ladrunoDR residualNorm] \
                      [ladrunoDR kineticEnergy] [ladrunoDR maxDisplacement]]
    }
    close $fh
}

set fh [open "${out}.node${pid}" w]
foreach i [lsort -integer $held] {
    foreach j {0 1} { foreach k {0 1} {
        set t [ntag $i $j $k]
        puts $fh [format "%d %.17g %.17g %.17g" $t [nodeDisp $t 1] [nodeDisp $t 2] [nodeDisp $t 3]]
    } }
}
close $fh

set fh [open "${out}.rc${pid}" w]
puts $fh "rc $rc"
puts $fh "result $result"
close $fh
