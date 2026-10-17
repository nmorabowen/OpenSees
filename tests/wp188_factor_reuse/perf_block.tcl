# WP-188 measurement deck (not a pytest case): a 3D brick block, linear-elastic
# below and J2-plastic in its top layers, under a base acceleration, for the
# OpenSeesMP + MUMPS wall-time comparison in ADR-188 §6. The point is the 4D
# regime in miniature: one large factorization per step, a nonlinear zone that
# keeps the tangent changing, a few iterations per step.
#
# Env: WP188_OUT (prefix), WP188_N (elements per side, default 24), WP188_NZ
# (layers, default 12), WP188_NPL (plastic top layers, default 3), WP188_REUSE
# (-reuseFactor, 0 = off), WP188_TOL (1e-6), WP188_NSTEP (60), WP188_PGA (0.5),
# WP188_REUSEXTRA (extra -reuse* options, WP-193).
# Run: mpiexec -localonly -n 4 OpenSeesMP.exe perf_block.tcl

proc envOr {name def} {
    global env
    if {[info exists env($name)] && $env($name) ne ""} { return $env($name) }
    return $def
}
set out   [envOr WP188_OUT block]
set N     [envOr WP188_N 24]
set NZ    [envOr WP188_NZ 12]
set NPL   [envOr WP188_NPL 3]
set reuse [envOr WP188_REUSE 0]
set tol   [envOr WP188_TOL 1.0e-6]
set nStep [envOr WP188_NSTEP 60]
set pga   [envOr WP188_PGA 0.5]

set pid [getPID]
set np  [getNP]
model BasicBuilder -ndm 3 -ndf 3

set h 1.0
proc nid {i j k} { global N; return [expr {1 + $i + ($N+1)*($j + ($N+1)*$k)}] }

# x-slab partition
set i0 [expr {($pid*$N)/$np}]
set i1 [expr {(($pid+1)*$N)/$np}]
for {set k 0} {$k <= $NZ} {incr k} {
    for {set j 0} {$j <= $N} {incr j} {
        for {set i $i0} {$i <= $i1} {incr i} {
            set n [nid $i $j $k]
            node $n [expr {$i*$h}] [expr {$j*$h}] [expr {$k*$h}]
            if {$k == 0} { fix $n 1 1 1 }
        }
    }
}
nDMaterial ElasticIsotropic 1 2.0e5 0.3 2.0
# J2Plasticity tag K G sig0 sigInf delta H
nDMaterial J2Plasticity 2 1.6667e5 7.6923e4 60.0 90.0 5.0 500.0
set e 0
for {set k 0} {$k < $NZ} {incr k} {
    set mat [expr {($k >= $NZ - $NPL) ? 2 : 1}]
    for {set j 0} {$j < $N} {incr j} {
        for {set i 0} {$i < $N} {incr i} {
            incr e
            if {$i < $i0 || $i >= $i1} continue
            element stdBrick $e [nid $i $j $k] [nid [expr $i+1] $j $k] [nid [expr $i+1] [expr $j+1] $k] [nid $i [expr $j+1] $k] \
                [nid $i $j [expr $k+1]] [nid [expr $i+1] $j [expr $k+1]] [nid [expr $i+1] [expr $j+1] [expr $k+1]] [nid $i [expr $j+1] [expr $k+1]] $mat 0.0 0.0 0.0
        }
    }
}
# brick density is not an option of stdBrick's short form: lumped nodal mass instead,
# rank-owned columns i in [i0, i1) (+ the last column on the last rank)
for {set k 1} {$k <= $NZ} {incr k} {
    for {set j 0} {$j <= $N} {incr j} {
        for {set i $i0} {$i <= $i1} {incr i} {
            if {$i == $i1 && $i1 != $N} continue
            mass [nid $i $j $k] 2.0 2.0 2.0
        }
    }
}

set dt 0.01
set vals {}
set pi [expr {acos(-1.0)}]
for {set s 0} {$s <= $nStep + 10} {incr s} {
    set t [expr {$s*$dt}]
    lappend vals [expr {min(1.0, $t/0.2) * sin(2*$pi*3.0*$t)}]
}
timeSeries Path 2 -dt $dt -values $vals -factor [expr {$pga*9.81}]
pattern UniformExcitation 2 1 -accel 2
rayleigh 0.3 0.0 0.002 0.0

constraints Transformation
numberer LadrunoParallelPlain
system Mumps -ICNTL14 200
test NormDispIncr $tol 100 0
if {$reuse > 0} {
    # WP-193: WP188_REUSEXTRA appends -reuse* options (e.g. "-reuseCost on")
    eval algorithm KrylovNewton -maxDim 20 -reuseFactor $reuse [envOr WP188_REUSEXTRA ""]
} else {
    algorithm KrylovNewton -maxDim 20
}
integrator Newmark 0.5 0.25
analysis Transient

set corner [nid 0 0 $NZ]
set hist {}
set t0 [clock milliseconds]
set fail 0
for {set s 0} {$s < $nStep} {incr s} {
    if {[analyze 1 $dt] != 0} { set fail 1; break }
    if {$pid == 0} { lappend hist [format %.17g [nodeDisp $corner 1]] }
}
set wall [expr {([clock milliseconds] - $t0)/1000.0}]
if {$pid == 0} {
    set fh [open "$out.txt" w]
    puts $fh "fail $fail steps $s wall $wall"
    puts $fh [join $hist " "]
    close $fh
}
print -file "$out.algo$pid" -algorithm
wipe
