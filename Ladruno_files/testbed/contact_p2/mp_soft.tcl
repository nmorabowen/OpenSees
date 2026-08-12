# ADR-78 P2 refusal deck: the P0 2-rank manual-DD deck with `-soft 0.1` added
# to the contact verb. EXPECTED OUTCOME: a FATAL named refusal at handle() on
# rank 0 ("-soft ... not supported under a partitioned host") followed by
# MPI_Abort tearing the job down -- the P2SOFTMP line below must NEVER print.
# k_soft needs BOTH surfaces' assembled mass; rank 0's cache has zero mass for
# the ghosted slave nodes 11..14, so SOFT would be silently too soft (D4).
wipe
model basic -ndm 3 -ndf 3
set pid [getPID]

# rho > 0 ON PURPOSE: with a massless model the P1 slave-mass abort would fire
# too, and this deck could "pass" via the WRONG abort. With mass, the ONLY
# abort available is the P2 partitioning refusal -- and a pre-P2 build runs
# this deck to completion, silently soft (the mutation signature).
nDMaterial ElasticIsotropic 1 2.0e7 0.0 2000.0
timeSeries Linear 1

if {$pid == 0} {
    node 1 0.0 0.0 0.0
    node 2 1.0 0.0 0.0
    node 3 1.0 1.0 0.0
    node 4 0.0 1.0 0.0
    node 5 0.0 0.0 1.0
    node 6 1.0 0.0 1.0
    node 7 1.0 1.0 1.0
    node 8 0.0 1.0 1.0
    element stdBrick 101 1 2 3 4 5 6 7 8 1
    foreach n {1 2 3 4} { fix $n 1 1 1 }
    foreach n {5 6 7 8} { fix $n 1 1 0 }

    node 11 0.0 0.0 0.999
    node 12 1.0 0.0 0.999
    node 13 1.0 1.0 0.999
    node 14 0.0 1.0 0.999
    foreach n {11 12 13 14} { fix $n 1 1 0 }

    contactSurface 1 -master 4 5 6 7 8
    contactSurface 2 -slave 11 12 13 14
    contact 1 1 2 auto -outward 0.0 0.0 1.0 -soft 0.1

    pattern Plain 1 1 { }
}

if {$pid == 1} {
    node 11 0.0 0.0 0.999
    node 12 1.0 0.0 0.999
    node 13 1.0 1.0 0.999
    node 14 0.0 1.0 0.999
    node 15 0.0 0.0 1.999
    node 16 1.0 0.0 1.999
    node 17 1.0 1.0 1.999
    node 18 0.0 1.0 1.999
    element stdBrick 102 11 12 13 14 15 16 17 18 1
    foreach n {11 12 13 14 15 16 17 18} { fix $n 1 1 0 }

    pattern Plain 1 1 {
        load 15 0.0 0.0 -2500.0
        load 16 0.0 0.0 -2500.0
        load 17 0.0 0.0 -2500.0
        load 18 0.0 0.0 -2500.0
    }
}

constraints LadrunoContact
numberer ParallelPlain
system Mumps
test NormDispIncr 1.0e-10 50 0
algorithm Newton
integrator LoadControl 1.0
analysis Static

set ok [analyze 1]
puts [format "P2SOFTMP pid=%d analyze=%d SHOULD-NOT-PRINT" $pid $ok]
wipe
