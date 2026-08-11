# P1 mutation: mp.tcl with rank 0's GHOST declarations of slave nodes 11..14
# REMOVED. Before ADR-78 P1 this warned ("contact SKIPPED") and limped on to a
# non-converged analyze=-3. After P1 it must ABORT at handle().
wipe
model basic -ndm 3 -ndf 3
set pid [getPID]

nDMaterial ElasticIsotropic 1 2.0e7 0.0
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

    # NO GHOSTS -- the mutation.

    contactSurface 1 -master 4 5 6 7 8
    contactSurface 2 -slave 11 12 13 14
    contact 1 1 2 auto -outward 0.0 0.0 1.0
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
puts "MUTNOGHOST pid=$pid analyze=$ok"
wipe
