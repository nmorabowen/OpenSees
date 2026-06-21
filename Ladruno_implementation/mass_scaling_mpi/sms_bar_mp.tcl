# ---------------------------------------------------------------------------
# Parallel lumped-SMS shared-node validation (OpenSeesMP, Tcl).
#
# A 1-D fixed-free truss bar of 20 elements with a FINE 4-element zone in the
# middle (elements 9..12). Under a 2-rank manual partition the bar splits at
# the central node 11, which is SHARED. The two fine elements straddling that
# node -- element 10 (rank 0) and element 11 (rank 1) -- are BOTH below dtTarget,
# so CentralDifferenceSMS injects fictitious mass (s-1)*m into node 11 from BOTH
# ranks. Correct behaviour requires the DistributedDiagonal/MPIDiagonal solver
# to SUM those per-rank contributions across the partition boundary.
#
# Validation: run the SAME global model at np=1 (whole bar on one rank -> trusted
# full-model reference, all injection local) and np=2 (split at the shared node).
# If the tip-displacement history matches, the shared-node DeltaM is reduced
# correctly across ranks -> the lumped SMS path is parallel-correct.
#
# Run:  mpiexec -n 1 OpenSeesMP.exe sms_bar_mp.tcl
#       mpiexec -n 2 OpenSeesMP.exe sms_bar_mp.tcl
# Writes tip_np<NP>.out (time  tipDisp) in the cwd.
# ---------------------------------------------------------------------------

set pid [getPID]
set np  [getNP]

# ---- physical / discretisation parameters --------------------------------
set E      1.0e4
set A      1.0
set rho    1.0
set N      20          ;# total elements (must be divisible by np for the split)
set hB     1.0         ;# bulk element length
set hF     0.1         ;# fine  element length  (-> dt_e 10x smaller)
set dtTar  0.004       ;# SMS target step  (between fine dt_e=1e-3 and bulk 1e-2)
set dtRun  0.0032      ;# 0.8*dtTar, stability margin
set nSteps 150
set Fapp   1.0         ;# step load at the free tip

# fine zone = the middle 4 elements (9,10,11,12 for N=20)
set fLo [expr {$N/2 - 1}]   ;# 9
set fHi [expr {$N/2 + 2}]   ;# 12

# ---- global node x-coordinates (identical on every rank) ------------------
# node i sits at the running sum of element lengths to its left.
set x(1) 0.0
for {set e 1} {$e <= $N} {incr e} {
    set h [expr {($e >= $fLo && $e <= $fHi) ? $hF : $hB}]
    set x([expr {$e+1}]) [expr {$x($e) + $h}]
}
set tipNode [expr {$N+1}]

model basic -ndm 1 -ndf 1
uniaxialMaterial Elastic 1 $E

# ---- ownership: rank r owns elements [r*N/np+1 .. (r+1)*N/np] -------------
set e0 [expr {$pid*$N/$np + 1}]
set e1 [expr {($pid+1)*$N/$np}]
set nLo $e0
set nHi [expr {$e1 + 1}]   ;# node tags this rank touches: nLo..nHi (inclusive)

for {set n $nLo} {$n <= $nHi} {incr n} {
    node $n $x($n)
}
for {set e $e0} {$e <= $e1} {incr e} {
    element Truss $e $e [expr {$e+1}] $A 1 -rho $rho
}

# fix the left end (node 1) on whichever rank owns it
if {$nLo <= 1 && 1 <= $nHi} { fix 1 1 }

# step load at the free tip on whichever rank owns it
if {$nLo <= $tipNode && $tipNode <= $nHi} {
    timeSeries Constant 1 -factor 1.0
    pattern Plain 1 1 { load $tipNode $Fapp }
    recorder Node -file tip_np${np}.out -time -node $tipNode -dof 1 disp
}

constraints Transformation
numberer ParallelPlain
algorithm Linear
system MPIDiagonal
integrator CentralDifferenceSMS $dtTar
analysis Transient

analyze $nSteps $dtRun

if {$pid == 0} { puts "DONE np=$np  pid=$pid  tipNode=$tipNode" }
wipe
