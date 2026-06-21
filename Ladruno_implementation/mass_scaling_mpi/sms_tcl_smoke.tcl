# Smoke test: the 6 selective-mass-scaling integrators are reachable from the LEGACY
# Tcl `integrator` parser (SRC/tcl/commands.cpp). Before the parser wiring these all
# failed with "No Integrator type exists in specifyIntegrator". Run:
#   OpenSees.exe sms_tcl_smoke.tcl
# Each integrator builds a tiny fixed-free truss bar, runs 5 explicit steps, and the
# script prints "<name>: OK" iff analyze returns 0 (i.e. the integrator was created
# and stepped). Consistent variants take the same $dtTarget first arg as the lumped.

set lumped   {CentralDifferenceSMS ExplicitBatheSMS ExplicitBatheLNVDSMS}
set consist  {CentralDifferenceSMSConsistent ExplicitBatheSMSConsistent ExplicitBatheLNVDSMSConsistent}
# ExplicitBathe family takes $p (rho-inf) first; LNVD takes $p $alpha first; the SMS
# $dtTarget follows. Per-type arg prefixes:
array set pre {
  CentralDifferenceSMS               "0.004"
  CentralDifferenceSMSConsistent     "0.004"
  ExplicitBatheSMS                   "0.5 0.004"
  ExplicitBatheSMSConsistent         "0.5 0.004"
  ExplicitBatheLNVDSMS               "0.5 0.0 0.004"
  ExplicitBatheLNVDSMSConsistent     "0.5 0.0 0.004"
}

proc run_one {name argline} {
    wipe
    model basic -ndm 1 -ndf 1
    uniaxialMaterial Elastic 1 1.0e4
    # 6-element bar, one short (fine) element so SMS actually scales something
    set x 0.0
    node 1 $x
    for {set e 1} {$e <= 6} {incr e} {
        set h [expr {$e == 3 ? 0.1 : 1.0}]
        set x [expr {$x + $h}]
        node [expr {$e+1}] $x
        element Truss $e $e [expr {$e+1}] 1.0 1 -rho 1.0
    }
    fix 1 1
    timeSeries Constant 1 -factor 1.0
    pattern Plain 1 1 { load 7 1.0 }
    constraints Plain
    numberer Plain
    algorithm Linear
    system Diagonal
    if {[catch {eval integrator $name $argline} err]} {
        puts "$name: FAIL (integrator create: $err)"
        return
    }
    analysis Transient
    set ok [analyze 5 0.003]
    if {$ok == 0} {
        puts "$name: OK"
    } else {
        puts "$name: FAIL (analyze returned $ok)"
    }
}

foreach name [concat $lumped $consist] {
    run_one $name $pre($name)
}
puts "DONE"
