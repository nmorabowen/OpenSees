# WP-194 test deck: the WP-188 RC frame (tests/wp188_factor_reuse/frame_gm.tcl,
# fiber forceBeamColumn or elastic, KrylovNewton, optional -reuseFactor) run
# unchanged, with two wrappers around it:
#   * every `system Mumps ...` gets $WP194_MUMPSARGS appended
#     (e.g. "-precision single -refine 1 -stats");
#   * every `analyze` adds the Newton iterations of its last step (testIter) to a
#     counter, written by rank 0 to <WP188_OUT>.iter as
#     "iters <sum of testIter> analyze <calls>".
# All WP188_* variables of frame_gm.tcl apply. Run under OpenSeesMP with >= 2
# ranks (frame_gm.tcl only selects `system Mumps` when np > 1).
# WP194_DECK picks another WP-188 deck from tests/wp188_factor_reuse/ (e.g.
# perf_block.tcl, the 3D J2 block used for the ADR-194 measurement).

proc wp194EnvOr {name def} {
    global env
    if {[info exists env($name)] && $env($name) ne ""} { return $env($name) }
    return $def
}
set wp194Args [wp194EnvOr WP194_MUMPSARGS ""]

rename system wp194_system
proc system {args} {
    global wp194Args
    if {[lindex $args 0] eq "Mumps"} {
        return [uplevel 1 [list wp194_system {*}$args {*}$wp194Args]]
    }
    return [uplevel 1 [list wp194_system {*}$args]]
}

set wp194Iter 0
set wp194Calls 0
rename analyze wp194_analyze
proc analyze {args} {
    global wp194Iter wp194Calls
    set r [uplevel 1 [list wp194_analyze {*}$args]]
    incr wp194Iter [testIter]
    incr wp194Calls
    return $r
}

source [file join [file dirname [file normalize [info script]]] .. wp188_factor_reuse [wp194EnvOr WP194_DECK frame_gm.tcl]]

if {$pid == 0} {
    set fh [open "$out.iter" w]
    puts $fh "iters $wp194Iter analyze $wp194Calls"
    close $fh
}
