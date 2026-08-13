# ADR-78 P2 review F1: a RETIRED soft interaction must not trip the
# partitioning refusal. The P0 2-rank deck plus a sacrificial NTS contact
# (tag 5) declared with -soft on bare, fully-fixed slave nodes 21..24 that
# the deck REMOVES before the first analyze. At handle():
#   pruneRemovedNodes retires contact 5 (named notice)  ->  the anySoft scan
#   must SKIP it  ->  no refusal, no mass-cache pass for it  ->  the live
#   -kn auto pair analyzes to the P0 values.
# Pre-fix signature: anySoft picked up the retired record, the per-contact
# message loops (which DO skip retired) printed nothing, and the job was
# MPI_Abort-ed behind a subject-less "...not supported under a partitioned
# host" tail. Belt and braces: had retirement itself failed, the P1
# massless-slave guard would abort this deck too -- ANY abort is a FAIL here.
wipe
model basic -ndm 3 -ndf 3
set pid [getPID]

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
    contact 1 1 2 auto -outward 0.0 0.0 1.0

    # sacrificial SOFT pair: bare fully-fixed slaves, removed before analyze
    node 21 0.0 0.0 1.5
    node 22 1.0 0.0 1.5
    node 23 1.0 1.0 1.5
    node 24 0.0 1.0 1.5
    foreach n {21 22 23 24} { fix $n 1 1 1 }
    contactSurface 8 -slave 21 22 23 24
    contact 5 1 8 1.0e6 1.0e5 0.0 -outward 0.0 0.0 1.0 -soft 0.1
    foreach n {21 22 23 24} { remove node $n }

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
if {$pid == 0} {
    puts [format "P2RETIRED pid=0 analyze=%d w5=%.15e" $ok [nodeDisp 5 3]]
} else {
    puts [format "P2RETIRED pid=1 analyze=%d w15=%.15e" $ok [nodeDisp 15 3]]
}
wipe
