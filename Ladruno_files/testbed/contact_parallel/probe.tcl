model basic -ndm 3 -ndf 3
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
if {[catch {contactSurface 1 -master 4 1 2 3 4} e]} {
    puts "CONTACT-ABSENT: $e"
} else {
    puts "CONTACT-PRESENT"
}
wipe
