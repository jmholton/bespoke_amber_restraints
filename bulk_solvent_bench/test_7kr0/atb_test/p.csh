#! /bin/tcsh -f
cd /programs/amber22 && source amber.csh
cd /tmp
printf 'help\nquit\n' | parmed |& sed -n '/Reading input/,$p' | tr -s ' \n' ' ' | fold -w 200
echo
printf 'help changeBox\nhelp setbox\nhelp chbox\nquit\n' | parmed |& sed -n '/Reading input/,$p'
