# ADR-78 P2 control: -soft on a SERIAL host must still run (the partitioning
# refusal is keyed on np>1 / Subdomain, and must NOT fire at np=1).
# Same two-block column as ../contact_parallel/serial.tcl, plus -soft 0.1.
wipe
model basic -ndm 3 -ndf 3

node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0

node 11 0.0 0.0 0.999
node 12 1.0 0.0 0.999
node 13 1.0 1.0 0.999
node 14 0.0 1.0 0.999
node 15 0.0 0.0 1.999
node 16 1.0 0.0 1.999
node 17 1.0 1.0 1.999
node 18 0.0 1.0 1.999

# rho > 0: -soft sizes k_soft from the ASSEMBLED mass, and the P1 guard
# (correctly) aborts a massless SLAVE surface even in serial. The refusal
# under test here is the np>1 one, which must NOT fire; give the model mass
# so the P1 guard stays quiet and the np=1 run demonstrably proceeds.
nDMaterial ElasticIsotropic 1 2.0e7 0.0 2000.0
element stdBrick 101 1 2 3 4 5 6 7 8 1
element stdBrick 102 11 12 13 14 15 16 17 18 1

foreach n {1 2 3 4} { fix $n 1 1 1 }
foreach n {5 6 7 8 11 12 13 14 15 16 17 18} { fix $n 1 1 0 }

contactSurface 1 -master 4 5 6 7 8
contactSurface 2 -slave 11 12 13 14
contact 1 1 2 auto -outward 0.0 0.0 1.0 -soft 0.1

timeSeries Linear 1
pattern Plain 1 1 {
    load 15 0.0 0.0 -2500.0
    load 16 0.0 0.0 -2500.0
    load 17 0.0 0.0 -2500.0
    load 18 0.0 0.0 -2500.0
}

constraints LadrunoContact
numberer RCM
system UmfPack
test NormDispIncr 1.0e-10 50 0
algorithm Newton
integrator LoadControl 1.0
analysis Static

set ok [analyze 1]
puts [format "P2SOFTSERIAL analyze=%d w15=%.15e" $ok [nodeDisp 15 3]]
wipe
