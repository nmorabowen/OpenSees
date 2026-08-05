# ADR-80 gates G1/G2/G3 — synthetic re-derivation of the Cerro Lindo mechanism.
#
# The ADR allows "a synthetic re-derivation of the same geometry" in place of the
# Cerro Lindo fuse decks. This is that: a 1D chain of hex8 bricks, ux/uy fixed
# everywhere and nu = 0, so sigma = E*eps exactly and every number below is
# hand-checkable.
#
#   base face uz fixed;  top face uz driven by a non-homogeneous `sp` (disp)
#                        or by an equivalent traction (trac).
#
# SIZING IS THE WHOLE POINT. With L = 100, delta = 0.15 over nsteps = 10:
#   converged strain  = delta/L      = 0.0015 -> sigma = 300 MPa  < fy = 379.5
#   driven-layer step = ddelta/h     = 0.015/h
# so for N = 20 (h = 5): trial sigma = 600 MPa >> fy. The CONVERGED answer never
# yields while the driven layer yields SPURIOUSLY inside every newStep --
# precisely the Cerro Lindo signature (converged 295 MPa vs fy 379.5).
#
#   OpenSees.exe g123_predictor_gates.tcl        -> ./g123_predictor_gates.json
#   set ::G123_OUT <path> before sourcing to redirect.

if {![info exists ::G123_OUT]} { set ::G123_OUT "g123_predictor_gates.json" }

set L      100.0
set AX      10.0        ;# cross-section side
set EMOD 200000.0
set KMOD [expr $EMOD/3.0]    ;# nu = 0  =>  K = E/3
set GMOD [expr $EMOD/2.0]    ;#            G = E/2
set FY     379.5
set HISO  2000.0        ;# H/E = 1 %
set DELTA    0.15       ;# total prescribed displacement (converged stays elastic)

# build a 1x1xN hex8 column; returns nothing, leaves the model current
proc build {N mat} {
    global L AX EMOD KMOD GMOD FY HISO
    wipe
    model basic -ndm 3 -ndf 3
    set h [expr $L/double($N)]
    for {set k 0} {$k <= $N} {incr k} {
        set z [expr $k*$h]
        set b [expr 4*$k]
        node [expr $b+1] 0.0   0.0   $z
        node [expr $b+2] $AX   0.0   $z
        node [expr $b+3] $AX   $AX   $z
        node [expr $b+4] 0.0   $AX   $z
    }
    if {$mat eq "j2"} {
        nDMaterial LadrunoJ2 1 $KMOD $GMOD -iso voce $FY 0.0 0.0 $HISO
    } else {
        nDMaterial ElasticIsotropic 1 $EMOD 0.0
    }
    for {set k 0} {$k < $N} {incr k} {
        set b [expr 4*$k]
        element stdBrick [expr $k+1] \
            [expr $b+1] [expr $b+2] [expr $b+3] [expr $b+4] \
            [expr $b+5] [expr $b+6] [expr $b+7] [expr $b+8] 1
    }
    # ux, uy fixed everywhere -> clean 1D chain; base uz fixed; top uz free for the drive
    for {set k 0} {$k <= $N} {incr k} {
        set b [expr 4*$k]
        for {set j 1} {$j <= 4} {incr j} {
            if {$k == 0} { fix [expr $b+$j] 1 1 1 } else { fix [expr $b+$j] 1 1 0 }
        }
    }
}

# drive: disp | trac ; mat: elastic | j2 ; algo: krylov | initial
proc run {N nsteps drive mat algo} {
    global L AX EMOD DELTA
    build $N $mat
    set top [expr 4*$N]
    timeSeries Linear 1
    pattern Plain 1 1 {
        if {$drive eq "disp"} {
            for {set j 1} {$j <= 4} {incr j} { sp [expr $top+$j] 3 $DELTA }
        } else {
            # equivalent traction: F = E*A*delta/L, split over the 4 corner nodes
            set F [expr $EMOD*$AX*$AX*$DELTA/$L/4.0]
            for {set j 1} {$j <= 4} {incr j} { load [expr $top+$j] 0.0 0.0 $F }
        }
    }
    constraints Transformation
    numberer RCM
    system BandGeneral
    test NormDispIncr 1.0e-10 60 0
    if {$algo eq "initial"} { algorithm Newton -initial } else { algorithm KrylovNewton }
    integrator LoadControl [expr 1.0/$nsteps]
    analysis Static

    set iters 0
    set fails 0
    for {set s 0} {$s < $nsteps} {incr s} {
        set rc [analyze 1]
        if {$rc != 0} { incr fails }
        set iters [expr $iters + [testIter]]
    }
    set utop [nodeDisp [expr $top+1] 3]
    return [list $iters $fails $utop]
}

set rows {}
proc emit {tag N nsteps drive mat algo} {
    global rows DELTA
    if {[catch {run $N $nsteps $drive $mat $algo} res]} {
        lappend rows "  {\"gate\": \"$tag\", \"N\": $N, \"nsteps\": $nsteps, \"drive\": \"$drive\", \"mat\": \"$mat\", \"algo\": \"$algo\", \"status\": \"error\"}"
        puts "ROW $tag N=$N ns=$nsteps $drive $mat $algo ERROR"
        return
    }
    lassign $res iters fails utop
    set err [expr abs($utop - $DELTA)]
    lappend rows "  {\"gate\": \"$tag\", \"N\": $N, \"nsteps\": $nsteps, \"drive\": \"$drive\", \"mat\": \"$mat\", \"algo\": \"$algo\", \"iters\": $iters, \"failed_steps\": $fails, \"u_top\": [format %.10g $utop]}"
    puts "ROW $tag N=$N ns=$nsteps [format %-5s $drive] [format %-7s $mat] [format %-7s $algo] iters=$iters fails=$fails utop=[format %.6f $utop]"
}

# ---- G0: the mechanism itself (drive x material), at the reference mesh -----
foreach drive {disp trac} {
    foreach mat {elastic j2} {
        emit G0 20 10 $drive $mat krylov
    }
}

# ---- G1: penalty must scale with L/h  (refine => worse) ---------------------
foreach N {5 10 20 40 80} {
    foreach mat {elastic j2} {
        emit G1 $N 10 disp $mat krylov
    }
}

# ---- G2: penalty vs step size, per unit of progress -------------------------
foreach ns {5 10 20 40} {
    foreach mat {elastic j2} {
        emit G2 20 $ns disp $mat krylov
    }
}

# ---- G3: which half dominates -- elastic tangent vs collapsed plastic one ---
# LadrunoJ2 DOES override getInitialTangent() with a true elastic tangent
# (LadrunoJ2.cpp:431), so `Newton -initial` is a real substitution here and not
# silently full Newton (the ADR-76 trap).
foreach algo {krylov initial} {
    foreach mat {elastic j2} {
        emit G3 20 10 disp $mat $algo
    }
}

set fh [open $::G123_OUT w]
puts $fh "{"
puts $fh "  \"gates\": \"ADR-80 G1/G2/G3 — synthetic Cerro Lindo re-derivation\","
puts $fh "  \"L\": $L, \"E\": $EMOD, \"nu\": 0.0, \"fy\": $FY, \"Hiso\": $HISO, \"delta\": $DELTA,"
puts $fh "  \"converged_sigma\": [expr $EMOD*$DELTA/$L],"
puts $fh "  \"rows\": \["
puts $fh [join $rows ",\n"]
puts $fh "  \]"
puts $fh "}"
close $fh
puts "WROTE $::G123_OUT"
