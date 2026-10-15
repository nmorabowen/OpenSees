# Per-stage KIND + stage-boundary regression deck (apeGmsh #2034), classic Tcl exe.
#
# Usage:  OpenSees.exe stage_kind_model.tcl <out.ladruno> <case>
# Prints "STAGE_KIND_MODEL OK <case>" as its last line on success (OpenSees.exe
# exits 0 even when a deck dies, so the checker keys on this marker).
#
# One 2-D cantilever (elastic beam, tip mass) and `recorder ladruno` with NO
# -kind unless the case says so. Every analysis records every step (-T dt 0).
# The cases and the MODEL_STAGE list stage_kind_check.py expects of them:
#
#   dc_static_transient  static 4 | loadConst, wipeAnalysis, domainChange | transient 5
#                        -> [static 4, transient 5]   (option B barrier + A flip)
#   flip_static_transient  the same with NO domainChange
#                        -> [static 4, transient 5]   (option A flip alone)
#   dc_static_static     static 4 | loadConst, new load pattern, domainChange | static 3
#                        -> [static 4, static 3]      (option B barrier alone)
#   repeat_static        static 4 | wipeAnalysis, a new Static analysis | static 3
#                        -> [static 7]                (a repeated kind is no boundary)
#   internal_restamp     static 2 | SP pattern | static 2 | eleLoad pattern | static 2
#                        -> [static 6]                (WP-165 folds internal re-stamps)
#   dc_first             domainChange before the first analyze | static 4
#                        -> [static 4]                (no stage for the first record)
#   transient_only       transient 5 -> [transient 5] (the G8 #1948 symptom)
#   record_before_analysis  `record` with no analysis | transient 5
#                        -> [transient 6]             (provisional KIND rewritten)
#   override_kind        -kind transient; static 4 | wipeAnalysis | transient 5 (no dc)
#                        -> [transient 9]             (-kind override: no flip split)
#   override_kind_dc     -kind static; static 4 | domainChange | transient 5
#                        -> [static 4, static 5]      (the barrier still splits)
#   flip_then_topology   static 4 | wipeAnalysis | transient 5 | add a node + element | transient 2
#                        -> [static 4, transient 5, transient 2], ids strictly increasing
#                           (a flip that did not move the stamp must not collide with
#                            the stamp of the later topology change)

if {$argc < 2} { puts "usage: stage_kind_model.tcl <out.ladruno> <case>"; exit 1 }
set out [lindex $argv 0]
set case [lindex $argv 1]
if {[file exists $out]} { file delete $out }

wipe
model basic -ndm 2 -ndf 3
node 1 0.0 0.0
node 2 0.0 3.0
fix 1 1 1 1
mass 2 10.0 10.0 0.0
geomTransf Linear 1
element elasticBeamColumn 1 1 2 0.1 2.0e8 1.0e-3 1

set kind_opt {}
if {$case == "override_kind"} { set kind_opt {-kind transient} }
if {$case == "override_kind_dc"} { set kind_opt {-kind static} }
recorder ladruno $out -N displacement {*}$kind_opt -T dt 0.0

proc static_chain {} {
    constraints Plain
    numberer RCM
    system BandGen
    test NormDispIncr 1.0e-10 10
    algorithm Newton
    integrator LoadControl 0.25
    analysis Static
}
proc transient_chain {} {
    constraints Plain
    numberer RCM
    system BandGen
    test NormDispIncr 1.0e-10 10
    algorithm Newton
    integrator Newmark 0.5 0.25
    analysis Transient
}
proc gravity {} {
    timeSeries Linear 1
    pattern Plain 1 1 { load 2 0.0 -100.0 0.0 }
}
proc lateral {} {
    timeSeries Sine 2 0.0 10.0 0.5
    pattern Plain 2 2 { load 2 100.0 0.0 0.0 }
}

switch -- $case {
    dc_static_transient {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        wipeAnalysis
        domainChange
        lateral; transient_chain; analyze 5 0.01
    }
    flip_static_transient {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        wipeAnalysis
        lateral; transient_chain; analyze 5 0.01
    }
    dc_static_static {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        timeSeries Linear 3
        pattern Plain 3 3 { load 2 50.0 0.0 0.0 }
        domainChange
        analyze 3
    }
    repeat_static {
        gravity; static_chain; analyze 4
        wipeAnalysis
        static_chain; analyze 3
    }
    internal_restamp {
        gravity; static_chain; analyze 2
        # an SP pattern moves the domain-change stamp (Domain::addLoadPattern)
        pattern Plain 5 1 { sp 2 3 0.0 }
        analyze 2
        # an element load (eleLoad) on the beam
        pattern Plain 6 1 { eleLoad -ele 1 -type -beamUniform -1.0 }
        analyze 2
    }
    dc_first {
        domainChange
        gravity; static_chain; analyze 4
    }
    transient_only {
        lateral; transient_chain; analyze 5 0.01
    }
    record_before_analysis {
        record
        lateral; transient_chain; analyze 5 0.01
    }
    override_kind {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        wipeAnalysis
        lateral; transient_chain; analyze 5 0.01
    }
    override_kind_dc {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        wipeAnalysis
        domainChange
        lateral; transient_chain; analyze 5 0.01
    }
    flip_then_topology {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        wipeAnalysis
        lateral; transient_chain; analyze 5 0.01
        node 3 1.0 3.0
        fix 3 1 1 1
        element elasticBeamColumn 2 2 3 0.1 2.0e8 1.0e-3 1
        analyze 2 0.01
    }
    default {
        puts "unknown case $case"
        exit 1
    }
}
remove recorders
puts "STAGE_KIND_MODEL OK $case"
