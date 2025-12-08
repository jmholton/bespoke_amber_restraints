#! /bin/csh -f
#
#
#
set start = "$1"
if("$start" == "") set start = 0

cat << EOF | tclsh
puts -nonewline "[clock format [clock seconds]  -format "%a %b %d %T %Z %Y"] [format "%.3f" [expr ([clock clicks -milliseconds])/1000.0 - $start]] "
EOF

