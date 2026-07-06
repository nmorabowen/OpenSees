model basic -ndm 3 -ndf 3
nDMaterial ElasticIsotropic 1 2.00000000e+08 0.25 2000.0
# soil column 1x1x10 m, Vs=200.0, Vp=346.4102
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0
node 9 0.0 0.0 2.0
node 10 1.0 0.0 2.0
node 11 1.0 1.0 2.0
node 12 0.0 1.0 2.0
node 13 0.0 0.0 3.0
node 14 1.0 0.0 3.0
node 15 1.0 1.0 3.0
node 16 0.0 1.0 3.0
node 17 0.0 0.0 4.0
node 18 1.0 0.0 4.0
node 19 1.0 1.0 4.0
node 20 0.0 1.0 4.0
node 21 0.0 0.0 5.0
node 22 1.0 0.0 5.0
node 23 1.0 1.0 5.0
node 24 0.0 1.0 5.0
node 25 0.0 0.0 6.0
node 26 1.0 0.0 6.0
node 27 1.0 1.0 6.0
node 28 0.0 1.0 6.0
node 29 0.0 0.0 7.0
node 30 1.0 0.0 7.0
node 31 1.0 1.0 7.0
node 32 0.0 1.0 7.0
node 33 0.0 0.0 8.0
node 34 1.0 0.0 8.0
node 35 1.0 1.0 8.0
node 36 0.0 1.0 8.0
node 37 0.0 0.0 9.0
node 38 1.0 0.0 9.0
node 39 1.0 1.0 9.0
node 40 0.0 1.0 9.0
node 41 0.0 0.0 10.0
node 42 1.0 0.0 10.0
node 43 1.0 1.0 10.0
node 44 0.0 1.0 10.0
fix 1 0 1 1
fix 2 0 1 1
fix 3 0 1 1
fix 4 0 1 1
fix 5 0 1 1
fix 6 0 1 1
fix 7 0 1 1
fix 8 0 1 1
fix 9 0 1 1
fix 10 0 1 1
fix 11 0 1 1
fix 12 0 1 1
fix 13 0 1 1
fix 14 0 1 1
fix 15 0 1 1
fix 16 0 1 1
fix 17 0 1 1
fix 18 0 1 1
fix 19 0 1 1
fix 20 0 1 1
fix 21 0 1 1
fix 22 0 1 1
fix 23 0 1 1
fix 24 0 1 1
fix 25 0 1 1
fix 26 0 1 1
fix 27 0 1 1
fix 28 0 1 1
fix 29 0 1 1
fix 30 0 1 1
fix 31 0 1 1
fix 32 0 1 1
fix 33 0 1 1
fix 34 0 1 1
fix 35 0 1 1
fix 36 0 1 1
fix 37 0 1 1
fix 38 0 1 1
fix 39 0 1 1
fix 40 0 1 1
fix 41 0 1 1
fix 42 0 1 1
fix 43 0 1 1
fix 44 0 1 1
element stdBrick 1 1 2 3 4 5 6 7 8 1
element stdBrick 2 5 6 7 8 9 10 11 12 1
element stdBrick 3 9 10 11 12 13 14 15 16 1
element stdBrick 4 13 14 15 16 17 18 19 20 1
element stdBrick 5 17 18 19 20 21 22 23 24 1
element stdBrick 6 21 22 23 24 25 26 27 28 1
element stdBrick 7 25 26 27 28 29 30 31 32 1
element stdBrick 8 29 30 31 32 33 34 35 36 1
element stdBrick 9 33 34 35 36 37 38 39 40 1
element stdBrick 10 37 38 39 40 41 42 43 44 1
element LysmerTriangle 101 1 2 3 2000.0 346.41016151 200.0
element LysmerTriangle 102 1 3 4 2000.0 346.41016151 200.0
setNodeVel 41 1 0.1 -commit
setNodeVel 42 1 0.1 -commit
setNodeVel 43 1 0.1 -commit
setNodeVel 44 1 0.1 -commit
recorder EnergyBalance -file "A1_energy.txt" -time
recorder Node -file "A1_vels.txt" -time -nodeRange 1 44 -dof 1 vel
constraints Plain
numberer RCM
if {[catch {system UmfPack} m]} { if {[catch {system ProfileSPD} m2]} { system FullGeneral } }
integrator Newmark 0.5 0.25
algorithm Linear
analysis Transient
set ok [analyze 600 0.001]
puts "ANALYZE_RETURN $ok"
wipe
puts "DONE"
exit
