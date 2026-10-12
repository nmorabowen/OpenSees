# WP-188 test deck: a 3-story, 2-bay RC frame (fiber forceBeamColumn, Concrete01 /
# Steel02) or its elastic twin under a synthetic ground motion, with
# `algorithm KrylovNewton` with or without `-reuseFactor`.
#
# Driven by environment variables (classic OpenSees.exe and OpenSeesMP.exe alike):
#   WP188_OUT       output prefix (required). Rank 0 writes <out>.edp (one line of
#                   EDPs) and <out>.roof (roof displacement history); every rank
#                   writes <out>.algo<rank> (`print -algorithm`).
#   WP188_MODEL     fiber | elastic                       (default fiber)
#   WP188_REUSE     -reuseFactor m; 0 = option absent      (default 0)
#   WP188_ITER      -reuseIter k   ("" = not given)
#   WP188_RATE      -reuseRate r   ("" = not given)
#   WP188_TOL       NormDispIncr tolerance                 (default 1e-8)
#   WP188_MAXIT     transient test max iterations          (default 50; gravity always 50)
#   WP188_SYSTEM    serial system: BandGeneral | UmfPack    (default BandGeneral)
#   WP188_PGA       peak ground acceleration in g          (default 0.6)
#   WP188_NSTEP     number of transient steps (dt 0.01)    (default 1000)
#   WP188_HALVE     step index at which dt is halved for 20 steps (default -1: never)
#   WP188_TPRINT    test print flag of the transient test  (default 0; 2 = one line per step)
#   WP188_MAXDIM    KrylovNewton -maxDim                   (default 3, the vanilla default)
#
# Under OpenSeesMP each rank builds the elements e with e % np == rank; rank 0
# holds every node (and all the mass), the others only the nodes they use.

proc envOr {name def} {
    global env
    if {[info exists env($name)] && $env($name) ne ""} { return $env($name) }
    return $def
}

set out    [envOr WP188_OUT ""]
set mtype  [envOr WP188_MODEL fiber]
set reuse  [envOr WP188_REUSE 0]
set rIter  [envOr WP188_ITER ""]
set rRate  [envOr WP188_RATE ""]
set tol    [envOr WP188_TOL 1.0e-8]
set maxIt  [envOr WP188_MAXIT 50]
set sys    [envOr WP188_SYSTEM BandGeneral]
set pga    [envOr WP188_PGA 0.6]
set nStep  [envOr WP188_NSTEP 1000]
set halve  [envOr WP188_HALVE -1]
set maxDim [envOr WP188_MAXDIM 3]
set tprint [envOr WP188_TPRINT 0]
if {$out eq ""} { puts "WP188_OUT not set"; exit 2 }

set pid 0
set np 1
if {[info commands getNP] ne ""} { set np [getNP]; set pid [getPID] }
set parallel [expr {$np > 1}]

model BasicBuilder -ndm 2 -ndf 3

# ---------------- geometry ---------------------------------------------------------
set nStory 3
set nBay 2
set H 3.0
set L 6.0
proc nodeTag {i j} { return [expr {100*$j + $i + 1}] }   ;# i = column line, j = floor

# elements: columns then beams, numbered from 1
set eleList {}
set e 0
for {set j 0} {$j < $nStory} {incr j} {
    for {set i 0} {$i <= $nBay} {incr i} {
        incr e
        lappend eleList [list $e col [nodeTag $i $j] [nodeTag $i [expr {$j+1}]]]
    }
}
for {set j 1} {$j <= $nStory} {incr j} {
    for {set i 0} {$i < $nBay} {incr i} {
        incr e
        lappend eleList [list $e beam [nodeTag $i $j] [nodeTag [expr {$i+1}] $j]]
    }
}

set mine {}
foreach row $eleList {
    if {[expr {[lindex $row 0] % $np}] == $pid} { lappend mine $row }
}
set myNodes {}
if {$pid == 0} {
    for {set j 0} {$j <= $nStory} {incr j} { for {set i 0} {$i <= $nBay} {incr i} { lappend myNodes [nodeTag $i $j] } }
} else {
    foreach row $mine { lappend myNodes [lindex $row 2] [lindex $row 3] }
    set myNodes [lsort -integer -unique $myNodes]
}
foreach n $myNodes {
    set i [expr {($n % 100) - 1}]
    set j [expr {$n / 100}]
    node $n [expr {$i*$L}] [expr {$j*$H}]
    if {$j == 0} { fix $n 1 1 1 }
}
# lumped floor mass (t), rank 0 only
set mFloor 40.0
if {$pid == 0} {
    for {set j 1} {$j <= $nStory} {incr j} {
        for {set i 0} {$i <= $nBay} {incr i} {
            set m [expr {($i == 0 || $i == $nBay) ? 0.5*$mFloor : $mFloor}]
            mass [nodeTag $i $j] $m $m 0.0
        }
    }
}

# ---------------- sections -----------------------------------------------------------
proc rcSection {tag b h cover nBar aBar} {
    set y [expr {$h/2.0}]
    set z [expr {$b/2.0}]
    set yc [expr {$y - $cover}]
    set zc [expr {$z - $cover}]
    section Fiber $tag {
        patch rect 1 10 1 [expr {-$yc}] [expr {-$zc}] $yc $zc
        patch rect 2 2 1 $yc [expr {-$z}] $y $z
        patch rect 2 2 1 [expr {-$y}] [expr {-$z}] [expr {-$yc}] $z
        patch rect 2 10 1 [expr {-$yc}] $zc $yc $z
        patch rect 2 10 1 [expr {-$yc}] [expr {-$z}] $yc [expr {-$zc}]
        layer straight 3 $nBar $aBar $yc $zc $yc [expr {-$zc}]
        layer straight 3 $nBar $aBar [expr {-$yc}] $zc [expr {-$yc}] [expr {-$zc}]
    }
}

geomTransf PDelta 1
geomTransf Linear 2
set Ec 2.6e7
if {$mtype eq "fiber"} {
    uniaxialMaterial Concrete01 1 -36000.0 -0.0025 -26000.0 -0.012
    uniaxialMaterial Concrete01 2 -30000.0 -0.0020 -6000.0 -0.005
    uniaxialMaterial Steel02 3 420000.0 2.0e8 0.01 18.0 0.925 0.15
    rcSection 1 0.45 0.45 0.04 3 4.9e-4
    rcSection 2 0.35 0.55 0.04 3 3.8e-4
    foreach row $mine {
        lassign $row t kind ni nj
        if {$kind eq "col"} {
            element forceBeamColumn $t $ni $nj 1 Lobatto 1 5
        } else {
            element forceBeamColumn $t $ni $nj 2 Lobatto 2 5
        }
    }
} else {
    foreach row $mine {
        lassign $row t kind ni nj
        if {$kind eq "col"} {
            element elasticBeamColumn $t $ni $nj 0.2025 $Ec 0.00342 1
        } else {
            element elasticBeamColumn $t $ni $nj 0.1925 $Ec 0.00485 2
        }
    }
}

# ---------------- analysis objects ----------------------------------------------------
proc setSystem {} {
    global parallel sys
    if {$parallel} {
        constraints Transformation
        numberer LadrunoParallelPlain
        system Mumps -ICNTL14 200
    } else {
        constraints Transformation
        numberer RCM
        system $sys
    }
}

proc reuseArgs {} {
    global reuse rIter rRate maxDim
    set a [list -maxDim $maxDim]
    if {$reuse > 0} {
        lappend a -reuseFactor $reuse
        if {$rIter ne ""} { lappend a -reuseIter $rIter }
        if {$rRate ne ""} { lappend a -reuseRate $rRate }
    }
    return $a
}

# gravity (static, reuse on as well: LoadControl is invertible, dt constant)
timeSeries Linear 1
pattern Plain 1 1 {
    if {$pid == 0} {
        for {set j 1} {$j <= $nStory} {incr j} {
            for {set i 0} {$i <= $nBay} {incr i} {
                set p [expr {($i == 0 || $i == $nBay) ? -200.0 : -400.0}]
                load [nodeTag $i $j] 0.0 $p 0.0
            }
        }
    }
}
setSystem
test NormDispIncr $tol 50 0
eval algorithm KrylovNewton [reuseArgs]
integrator LoadControl 0.1
analysis Static
set ok [analyze 10]
if {$ok != 0} { puts "WP188 gravity failed"; exit 3 }
loadConst -time 0.0
wipeAnalysis

# ---------------- ground motion --------------------------------------------------------
set dt 0.01
set g 9.81
set nGm [expr {$nStep + 200}]
set vals {}
set pi [expr {acos(-1.0)}]
for {set k 0} {$k <= $nGm} {incr k} {
    set t [expr {$k*$dt}]
    if {$t < 1.0} { set env_ $t } elseif {$t < 5.0} { set env_ 1.0 } else { set env_ [expr {exp(-0.6*($t-5.0))}] }
    set a [expr {$env_*(0.55*sin(2*$pi*1.1*$t) + 0.30*sin(2*$pi*2.7*$t + 0.4) + 0.15*sin(2*$pi*5.3*$t + 1.1))}]
    lappend vals $a
}
timeSeries Path 2 -dt $dt -values $vals -factor [expr {$pga*$g}]
pattern UniformExcitation 2 1 -accel 2

# 5 % Rayleigh on modes 1 and 3 (fixed coefficients, no eigen needed)
rayleigh 0.35 0.0 0.0035 0.0

setSystem
test NormDispIncr $tol $maxIt $tprint
eval algorithm KrylovNewton [reuseArgs]
integrator Newmark 0.5 0.25
analysis Transient

set roof [nodeTag 0 $nStory]
set roofHist {}
set maxRoof 0.0
set maxIdr 0.0
set maxV 0.0
set maxM 0.0
set fail 0
set step 0
set nSub 0
while {$step < $nStep} {
    set h $dt
    if {$halve >= 0 && $step >= $halve && $step < $halve + 20} { set h [expr {$dt/2.0}] }
    set n [expr {($h < $dt) ? 2 : 1}]
    for {set s 0} {$s < $n} {incr s} {
        set ok [analyze 1 $h]
        incr nSub
        if {$ok != 0} { set fail 1; break }
    }
    if {$fail} { break }
    incr step
    if {$pid == 0} {
        set u [nodeDisp $roof 1]
        lappend roofHist [format %.17g $u]
        if {abs($u) > $maxRoof} { set maxRoof [expr {abs($u)}] }
        set uPrev 0.0
        for {set j 1} {$j <= $nStory} {incr j} {
            set uj [nodeDisp [nodeTag 0 $j] 1]
            set idr [expr {abs($uj - $uPrev)/$H}]
            if {$idr > $maxIdr} { set maxIdr $idr }
            set uPrev $uj
        }
        if {!$parallel} {
            # base shear and base moment from the first-story columns (elements 1..3)
            set V 0.0
            for {set c 1} {$c <= $nBay+1} {incr c} {
                set f [eleForce $c]
                set V [expr {$V + [lindex $f 0]}]
                set M [expr {abs([lindex $f 2])}]
                if {$M > $maxM} { set maxM $M }
            }
            if {abs($V) > $maxV} { set maxV [expr {abs($V)}] }
        }
    }
}

if {$pid == 0} {
    set resid [nodeDisp $roof 1]
    set fh [open "$out.edp" w]
    puts $fh [format "fail %d steps %d substeps %d maxRoof %.17g maxIdr %.17g resid %.17g maxV %.17g maxM %.17g" \
        $fail $step $nSub $maxRoof $maxIdr $resid $maxV $maxM]
    close $fh
    set fh [open "$out.roof" w]
    puts $fh [join $roofHist "\n"]
    close $fh
}
print -file "$out.algo$pid" -algorithm
wipe
puts "WP188 DONE pid $pid fail $fail"
