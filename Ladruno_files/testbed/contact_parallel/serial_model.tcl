# Model-only half of serial.tcl (no analysis), so query.tcl can drive a domain
# that actually has a contact defined.
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

nDMaterial ElasticIsotropic 1 2.0e7 0.0
element stdBrick 101 1 2 3 4 5 6 7 8 1
element stdBrick 102 11 12 13 14 15 16 17 18 1

foreach n {1 2 3 4} { fix $n 1 1 1 }
foreach n {5 6 7 8 11 12 13 14 15 16 17 18} { fix $n 1 1 0 }

contactSurface 1 -master 4 5 6 7 8
contactSurface 2 -slave 11 12 13 14
contact 1 1 2 auto -outward 0.0 0.0 1.0
