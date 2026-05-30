# Shared-params reader for the Tcl half (D1). Reads the SAME flat key/value
# file the pytest side reads via _testbed.params.load_params — single source of
# truth, no generator, impossible to diverge.
#
#   source _params.tcl
#   read_params [file join .. params forceBeamColumn_cantilever.params] P
#   set L $P(L)
proc read_params {path arrName} {
    upvar 1 $arrName P
    set f [open $path r]
    while {[gets $f line] >= 0} {
        set line [string trim [lindex [split $line "#"] 0]]
        if {$line eq ""} continue
        set key [lindex $line 0]
        set val [lrange $line 1 end]
        set P($key) [string trim $val]
    }
    close $f
}
