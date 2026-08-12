# Serial reference. Which element population is built is chosen by argv:
#
#   OpenSees.exe serial.tcl soft    -> expect PVAL = 1.000000e+05
#   OpenSees.exe serial.tcl stiff   -> expect PVAL = 1.000000e+09
#   OpenSees.exe serial.tcl whole   -> expect PVAL = 1.000000e+07
#
# `whole` is the reference the 2-rank run must reproduce on every rank: it is
# the same set of elements, just not partitioned.
wipe
model basic -ndm 2 -ndf 2
source [file join [file dirname [info script]] common.tcl]

set which [lindex $argv 0]
if {$which eq ""} { set which whole }

switch -- $which {
    soft   { buildSoft }
    stiff  { buildStiff }
    whole  { buildSoft ; buildStiff }
    default { puts "AUTOPEN ERROR unknown part '$which'" ; exit 1 }
}
buildFallbackMP 9000

constraints Auto -verbose
numberer Plain
system BandGeneral
test NormDispIncr 1.0e-10 50 0
algorithm Linear
integrator LoadControl 1.0
analysis Static

set ok [analyze 1]
puts [format "AUTOPEN serial part=%s analyze=%d" $which $ok]
wipe
