# ladrunoSchema capability `mdStageFlagPreserved` keeps its promise (apeGmsh #2517).
#
# Classic-Tcl twin of tests/test_md_stage_flag_capability.py. The build must
# list the capability in `ladrunoSchema`'s "capabilities", AND a ManzariDafalias /
# LadrunoSANISAND constructed after `updateMaterialStage ... -stage 1` (a new
# nDMaterial, full constructor, attached to no element) must leave the class-wide
# stage flag plastic: uz(5) equals the flip-only leg (rel 1e-9) and is far from
# the elastic (no-flip) leg. Before 43b235bbb the full constructors reset
# mElastFlag = 0 and uz(5) came back elastic.
# Stock 3.8.0 on this deck (MD): -3.4067686516e-3 (flip), -2.8970896097e-3 (no flip).
# Terminal marker: "MD_STAGE_FLAG_CAP_2517: PASS" or "MD_STAGE_FLAG_CAP_2517: FAIL ...".
# OpenSees.exe exits 0 on a Tcl error, so only the marker proves the deck ran.

set nFail 0
set MD {125 0.05 0.8 1.25 0.712 0.019 0.934 0.7 100 0.01 7.05 0.968 1.1 0.704 3.5 4 600 1.42}
set XY {{0. 0.} {1. 0.} {1. 1.} {0. 1.}}

# --- the capability is listed --------------------------------------------------
set schema [ladrunoSchema]
if {![regexp {"capabilities":\[([^\]]*)\]} $schema -> capList]} {
    puts "FAIL ladrunoSchema has no capabilities list"
    incr nFail
} elseif {[lsearch -exact [split $capList ,] {"mdStageFlagPreserved"}] < 0} {
    puts "FAIL mdStageFlagPreserved not listed: \[$capList\]"
    incr nFail
} else {
    puts "PASS mdStageFlagPreserved listed"
}

proc run {mat flip construct} {
    global MD XY
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    for {set k 0} {$k < 2} {incr k} {
        for {set j 0} {$j < 4} {incr j} {
            set xy [lindex $XY $j]
            node [expr {4*$k + $j + 1}] [lindex $xy 0] [lindex $xy 1] [expr {double($k)}]
        }
    }
    nDMaterial $mat 1 {*}$MD
    element stdBrick 1 1 2 3 4 5 6 7 8 1
    for {set n 1} {$n <= 8} {incr n} { fix $n 1 1 [expr {$n <= 4 ? 1 : 0}] }
    updateMaterialStage -material 1 -stage 0
    timeSeries Linear 1
    pattern Plain 1 1 { for {set n 5} {$n <= 8} {incr n} { load $n 0. 0. -25. } }
    constraints Plain
    numberer RCM
    system BandGeneral
    test NormDispIncr 1e-8 50 0
    algorithm Newton
    integrator LoadControl 0.25
    analysis Static
    if {[analyze 4] != 0} { return "FAILED_STAGE0" }
    loadConst -time 0.0
    if {$flip} { updateMaterialStage -material 1 -stage 1 }
    # a second instance through the full constructor, attached to no element
    if {$construct} { nDMaterial $mat 2 {*}$MD }
    timeSeries Linear 2
    pattern Plain 2 2 { for {set n 5} {$n <= 8} {incr n} { load $n 0. 0. -25. } }
    integrator LoadControl 0.05
    if {[analyze 20] != 0} { return "FAILED_STAGE1" }
    return [nodeDisp 5 3]
}

foreach mat {ManzariDafalias LadrunoSANISAND} {
    set plastic [run $mat 1 0]
    set plasticNew [run $mat 1 1]
    set elastic [run $mat 0 0]
    puts "$mat uz5: flip+construct=$plasticNew flip=$plastic elastic=$elastic"
    if {![string is double -strict $plastic] || ![string is double -strict $plasticNew] ||
        ![string is double -strict $elastic]} {
        puts "FAIL $mat: a leg did not converge"
        incr nFail
        continue
    }
    if {abs($plasticNew - $plastic) <= 1e-9 * abs($plastic)} {
        puts "PASS $mat: a new instance keeps the plastic stage"
    } else {
        puts "FAIL $mat: constructing after the flip changed uz5 (stage reset, #2517/#2444)"
        incr nFail
    }
    if {abs($plastic - $elastic) > 0.1 * abs($elastic)} {
        puts "PASS $mat: plastic and elastic legs are distinguishable"
    } else {
        puts "FAIL $mat: plastic and elastic legs coincide"
        incr nFail
    }
}

if {$nFail == 0} {
    puts "MD_STAGE_FLAG_CAP_2517: PASS"
} else {
    puts "MD_STAGE_FLAG_CAP_2517: FAIL ($nFail checks)"
}
