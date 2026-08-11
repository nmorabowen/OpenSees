# Runtime smoke for the EIGHT contact commands the P0 deck never called.
# Each is invoked bare. Reaching the shared OPS_* parser is proven either by a
# real result or by that parser's OWN usage/WARNING text -- both mean the bridge
# worked. "invalid command name" would mean the registration did not take.
source serial_model.tcl

foreach cmd {ladrunoContactInfo ladrunoContactForce ladrunoMortarPenetration
             ladrunoMortarTieResidual ladrunoEdgePenetration
             ladrunoBeginAugment ladrunoEndAugment} {
    if {[catch {$cmd} res]} {
        if {[string match "invalid command name*" $res]} {
            puts "QUERY $cmd : NOT-REGISTERED -> $res"
        } else {
            puts "QUERY $cmd : reached parser (error: [string range $res 0 60])"
        }
    } else {
        puts "QUERY $cmd : reached parser (returned: [string range $res 0 60])"
    }
}
# contactPlane is a definition verb; prove it parses too (bare = usage warning).
if {[catch {contactPlane} res]} {
    if {[string match "invalid command name*" $res]} {
        puts "QUERY contactPlane : NOT-REGISTERED -> $res"
    } else {
        puts "QUERY contactPlane : reached parser"
    }
} else {
    puts "QUERY contactPlane : reached parser"
}
wipe
