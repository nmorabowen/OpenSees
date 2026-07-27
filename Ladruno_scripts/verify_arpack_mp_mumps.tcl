# Ladruno ADR-1000 §31.5 / §34 smoke: plain `eigen` (ARPACK shift-invert) under
# OpenSeesMP — the composition apeGmsh ADR 0077 F1 refuted in vanilla
# (ArpackSOE left at processID -1, M*v never merged across ranks).
#
# Fixed-free chain, nMass masses, k=100, m=1:
#   lambda_j = 4*(k/m)*sin^2((2j-1)*pi/(2*(2*nMass+1)))
#
# THREE CONFIGURATIONS, FOUR RECORDED CHECKS (PR #668 review: configuration
# A alone cannot detect a lost gate or a broken re-gate):
#
#   A  partitioned deck + `system Mumps`   -> the gate's TRUE side. Wiring
#      must be ACTIVE: without it, rank 0 returns an empty spectrum and
#      rank 1 deadlocks (the observed pre-fix control).
#   B  REPLICATED deck (every rank builds the whole model) + serial system
#      -> the gate's FALSE side. Wiring must stay DORMANT: a regression
#      that wires unconditionally sums np copies of M*v, scaling every
#      eigenvalue by ~1/np. Case A passes unchanged under that regression,
#      so this case is what pins the gate.
#   C  replicated deck, `system Mumps` + eigen then `system UmfPack` + eigen
#      on the SAME reused SOE -> the REUSE re-gate (commands.cpp
#      `else if typeSolver == Arpack`), BOTH directions. Reached only by a
#      SECOND consecutive plain eigen, which no other deck in the repo does.
#      Replicated because a partitioned deck cannot switch to a serial
#      system at all (ParallelPlain numbering fails LinearSOE::setSize).
#
# Partition discipline (case A): element e -> rank (e-1)*np/nel (contiguous
# blocks); boundary nodes are defined on both touching ranks; a node's MASS
# is owned by the rank of its left element only — double-defined nodal mass
# is summed twice by the M*v merge (same owner rule every MP static deck
# uses for nodal loads).
#
#   serial oracle:  OpenSees   verify_arpack_mp_mumps.tcl
#   distributed:    mpiexec -n 2 OpenSeesMP verify_arpack_mp_mumps.tcl
#
# GATE USAGE — do NOT key on the exit code. Under _PARALLEL_INTERPRETERS a
# Tcl error exits 0 (g3TclMain drops its exitCode; mpiParameterMain returns
# 0 unconditionally), and a rank that `exit`s while a peer is deadlocked
# blocks in MPI_Finalize forever. Run under a timeout and require
#   ARPACK_MP_SMOKE_PASS  ...  cases=4
# on EVERY rank; treat a timeout, a missing line, or fewer than 4 checks as
# FAIL. Everything below is catch-wrapped so a Tcl error still prints
# ARPACK_MP_SMOKE_FAIL rather than dying silently.

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
set pi [expr {acos(-1.0)}]

proc eleRank {e np nel} { expr {($e - 1) * $np / $nel} }

# --- model builders -------------------------------------------------------
# mode: "partitioned" (this rank's element block only) or "replicated"
# (the whole chain on every rank). Returns {eLo eHi}.
proc buildChain {mode pid np nel nMass k m} {
    wipe
    model basic -ndm 1 -ndf 1
    uniaxialMaterial Elastic 1 $k

    if {$mode eq "replicated"} {
        set eLo 1
        set eHi $nel
    } else {
        set eLo 0
        set eHi -1
        for {set e 1} {$e <= $nel} {incr e} {
            if {[eleRank $e $np $nel] == $pid} {
                if {$eLo == 0} { set eLo $e }
                set eHi $e
            }
        }
        if {$eLo == 0} { error "rank $pid drew no elements (np > nel?)" }
    }

    for {set j $eLo} {$j <= [expr {$eHi + 1}]} {incr j} {
        node $j [expr {double($j - 1)}]
    }
    for {set e $eLo} {$e <= $eHi} {incr e} {
        element truss $e $e [expr {$e + 1}] 1.0 1
    }
    if {$eLo == 1} { fix 1 1 }

    for {set j 2} {$j <= [expr {$nMass + 1}]} {incr j} {
        if {$j < $eLo + 1 || $j > $eHi + 1} { continue }
        # replicated: every rank owns every mass (each holds the FULL model,
        # so the merge must NOT fire). partitioned: left-element owner only.
        if {$mode eq "replicated" || [eleRank [expr {$j - 1}] $np $nel] == $pid} {
            mass $j $m
        }
    }
    return [list $eLo $eHi]
}

# --- checks ---------------------------------------------------------------
proc checkSpectrum {lam nEig nMass k m} {
    global pi
    set worst -1.0
    for {set j 1} {$j <= $nEig} {incr j} {
        set lamA [expr {4.0 * $k / $m * pow(sin((2*$j - 1) * $pi / (2.0 * (2*$nMass + 1))), 2)}]
        set lamN [lindex $lam [expr {$j - 1}]]
        if {$lamN eq ""} { return 1.0e99 }
        set rel [expr {abs($lamN - $lamA) / $lamA}]
        if {$rel > $worst} { set worst $rel }
    }
    return $worst
}

# mode-1 shape ratios against sin((2*1-1)*pi*(i-1)/(2*nMass+1)); ratios to a
# local reference node cancel the arbitrary normalization. A rank-local
# (unmerged-M) Lanczos gives garbage ratios even if an eigenvalue matched.
proc checkShape {eLo eHi nMass} {
    global pi
    set refNode [expr {$eLo == 1 ? 2 : $eLo}]
    set refA [expr {sin($pi * ($refNode - 1) / (2.0 * $nMass + 1))}]
    set refN [lindex [nodeEigenvector $refNode 1] 0]
    if {$refN == 0.0} { return 1.0e99 }
    set worst -1.0
    for {set j $eLo} {$j <= [expr {$eHi + 1}]} {incr j} {
        if {$j == 1} { continue }
        set phiA [expr {sin($pi * ($j - 1) / (2.0 * $nMass + 1)) / $refA}]
        set phiN [expr {[lindex [nodeEigenvector $j 1] 0] / $refN}]
        set err [expr {abs($phiN - $phiA)}]
        if {$err > $worst} { set worst $err }
    }
    return $worst
}

# --- the cases ------------------------------------------------------
set results {}
set failures {}

proc record {name relErr shapeErr} {
    global results failures
    lappend results [format "%s(rel=%.3g,shape=%.3g)" $name $relErr $shapeErr]
    if {$relErr > 1.0e-8 || $relErr < 0.0 || $shapeErr > 1.0e-6 || $shapeErr < 0.0} {
        lappend failures "$name rel=$relErr shape=$shapeErr"
    }
}

set ok [catch {

    # ---- case A: partitioned + Mumps (gate TRUE, wiring must be ACTIVE) ----
    set r [buildChain partitioned $pid $np $nel $nMass $k $m]
    set eLo [lindex $r 0]
    set eHi [lindex $r 1]
    constraints Transformation
    if {$np > 1} { numberer ParallelPlain ; system Mumps } \
                 else { numberer RCM ; system UmfPack }
    set lamA_ [eigen $nEig]
    record A [checkSpectrum $lamA_ $nEig $nMass $k $m] [checkShape $eLo $eHi $nMass]

    # ---- case B: replicated deck + serial system (gate FALSE, DORMANT) ----
    # Every rank builds the FULL chain. If a regression wires the merge
    # unconditionally, np copies of M*v are summed while the solve stays
    # local -> every eigenvalue scaled by ~1/np, rel err ~ (1 - 1/np), far
    # above the 1e-8 gate. This is the case that pins the gate's existence.
    #
    # (A partitioned deck CANNOT be used here: with `numberer ParallelPlain`
    # a serial system fails LinearSOE::setSize() in domainChanged before
    # eigen runs. Replicated is the only legal both-systems configuration,
    # which is also why case C below is replicated.)
    set r [buildChain replicated $pid $np $nel $nMass $k $m]
    set bLo [lindex $r 0]
    set bHi [lindex $r 1]
    constraints Transformation
    numberer RCM
    system UmfPack
    set lamB [eigen $nEig]
    record B [checkSpectrum $lamB $nEig $nMass $k $m] [checkShape $bLo $bHi $nMass]

    # ---- case C: reuse re-gate, BOTH directions, on the SAME SOE ---------
    # No wipe: theEigenSOE from case B survives, so every eigen below takes
    # the reuse-on-classTag-match branch, not the creation path.
    #
    # C1 serial -> Mumps: re-wiring must FIRE. On a replicated deck the
    #    MumpsParallelSOE sums duplicate contributions, giving np*K and
    #    np*M; the merged M*v is likewise np*(M*v), so the np cancels in the
    #    generalized problem and the correct spectrum must come back. If the
    #    re-wire does NOT fire, M*v stays local while the solve carries
    #    np*K -> eigenvalues scaled by np.
    if {$np > 1} {
        system Mumps
        set lamC1 [eigen $nEig]
        record C1 [checkSpectrum $lamC1 $nEig $nMass $k $m] [checkShape $bLo $bHi $nMass]

        # C2 Mumps -> serial on the SAME SOE: the else-branch must fire
        #    (setProcessID(-1)). If the merge stays wired, np copies of M*v
        #    are summed against a local solve -> 1/np scaling, as in case B.
        system UmfPack
        set lamC2 [eigen $nEig]
        record C2 [checkSpectrum $lamC2 $nEig $nMass $k $m] [checkShape $bLo $bHi $nMass]
    } else {
        # serial binary: the reuse branch is not compiled, but run the same
        # two calls so the case count matches and the SOE reuse is exercised.
        system UmfPack
        set lamC1 [eigen $nEig]
        record C1 [checkSpectrum $lamC1 $nEig $nMass $k $m] [checkShape $bLo $bHi $nMass]
        set lamC2 [eigen $nEig]
        record C2 [checkSpectrum $lamC2 $nEig $nMass $k $m] [checkShape $bLo $bHi $nMass]
    }

} err]

set fp [open arpack_mp_mumps_rank$pid.out w]
puts $fp "np=$np cases=[llength $results]"
puts $fp "results=$results"
if {$ok} { puts $fp "error=$err" }
if {[llength $failures]} { puts $fp "failures=$failures" }
close $fp

if {$ok} {
    puts "ARPACK_MP_SMOKE_FAIL rank=$pid tcl-error: $err"
    exit 1
}
if {[llength $failures]} {
    puts "ARPACK_MP_SMOKE_FAIL rank=$pid $failures"
    exit 1
}
puts "ARPACK_MP_SMOKE_PASS rank=$pid np=$np cases=[llength $results] $results"
