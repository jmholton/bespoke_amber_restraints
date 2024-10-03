#! /bin/tcsh -f
#
# phased sum of a stack of MTZ files
#
#
set mtzs = ( $* )

set outfile = sum.mtz

# cannot migrate hosts because of temp files
set thishost = `hostname -s`
set test = `sinfo -h -n $thishost |& egrep -v "drain|n/a" | wc -l`
if ( $test ) then
  echo "using slurm"
  set srun = "srun -w $thishost"
else
  set srun = ""
endif

set t = /dev/shm/jamesh/temp_$$_mtzsum/
mkdir -p ${t}

cat << EOF >! ${t}mergemtz.csh
#! /bin/tcsh -f
sftools << eof
read \$1
read \$2
set labels
F1
P1
F2
P2
calc ( COL Fsum PHIsum ) = ( COL F1 P1 ) ( COL F2 P2 ) +
write \$3 col Fsum PHIsum
quit
y
eof
EOF
chmod a+x ${t}mergemtz.csh

unset done
set r = 0
while ( ! $?done ) 
set nextmtzs = ()
set done
@ r = ( $r + 1 )
set q = 0
set i = 1
while ( $i <= $#mtzs )
  @ q = ( $q + 1 )
  @ j = ( $i + 1 )
  if ( $j > $#mtzs ) set mtzs = ( $mtzs "0.mtz" )
  set mtz1 = $mtzs[$i]
  set mtz2 = $mtzs[$j]
  set a = `basename $mtz1 .mtz`
  set b = `basename "$mtz2" .mtz`
  set newmtz = ${t}/_${r}-${q}.mtz
  if(! -e "$mtz2") then
    echo "cp $mtz1 $newmtz"
    rm -f "${newmtz}" > /dev/null
sftools << EOF > /dev/null
read $mtz1
set labels
Fsum
PHIsum
write $newmtz col Fsum PHIsum
quit
y
EOF
  else
    unset done
    echo "$mtz1 + $mtz2 = $newmtz"
    $srun ${t}mergemtz.csh $mtz1 $mtz2 $newmtz >&! ${t}/merge.${a}-${b}.log &
  endif
  set nextmtzs = ( $nextmtzs $newmtz )
  @ i = ( $i + 2 )
end
wait
set mtzs = ( $nextmtzs )
end

if( $#mtzs != 1 ) then
  set BAD = WTF
  goto exit
endif

cp $mtzs $outfile
ls -l $outfile

echo "cleaning up..."
rm -rf ${t}

exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9
endif

exit


