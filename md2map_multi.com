#! /bin/tcsh -f
#
#	c-shell script for converting a molecular dynamics trajectory into structure factors
#	using the CCP4 Suite of programs
#
#
#
set SG = P1
set B  = 5
set opts = ""

set MD_mult = ( 1 1 1 )
if(-e cell.txt) set CELL = `cat cell.txt`
set GRID = `cat grid.txt`
if(-e MD_mult.txt) set MD_mult = `cat MD_mult.txt`

alias rsh /usr/bin/rsh
alias rsh ssh -x

set scriptdir = `dirname $0`
set scriptdir = `cd $scriptdir ; pwd`
if(! -x ${scriptdir}/md2map.com) set scriptdir = ../

echo "args: $*"

set pdbfiles = pdbfiles$$.txt
echo -n "" >! $pdbfiles
foreach arg ( $* )
    # change space group of comparison
    if( "$arg" =~ [PpCcIiFfRrHh][1-6]* ) then
        set SG = "$arg"
	continue
    endif
    # change B factor
    if( "$arg" =~ B=* ) then
        set B = `echo "$arg" | awk -F "=" '{print $2+0}'`
	continue
    endif

    if(("$arg" =~ *.txt)) then
        # warn about probable mispellings
        if(! -e "$arg") then
            echo "WARNING: $arg does not exist"
            continue
        endif
        
        cat $arg >> $pdbfiles
        continue
    endif

    if(("$arg" =~ *.pdb)) then
        # warn about probable mispellings
        if(! -e "$arg") then
            echo "WARNING: $arg does not exist"
            continue
        endif
        
        echo "$arg" >> $pdbfiles
        continue
    endif

    # pass-thru other options
    set opts = ( $opts $arg )
end

set maxCPUs = `cat $pdbfiles | wc -l | awk '$1>50{$1=50} {print}'`

set machines  = ( localhost )
set cpu_count = ( $maxCPUs )

foreach machine ( $machines )
    if("$machine" == "localhost") then
        killall -g md2map.com >& /dev/null
    else
        rsh -n $machine "killall -g md2map.com" >& /dev/null
    endif
end
echo "old jobs killed..."
sleep 1

cat << EOF
SG = $SG
B  = $B
EOF

set test = `cat $pdbfiles | wc -l`
if(! $test) then
    echo "defaulting to all trajectory/*.pdb files"
    ls -1 trajectory/*.pdb | \
    awk '/.pdb$/{file=$0;while(gsub("[0-9]$|.pdb$",""));printf("%05d %s\n",substr(file,length($0)+1),file) }' |\
    sort -g |\
    awk '{print $2}' >! $pdbfiles
endif
set pdbs = `cat $pdbfiles | wc -l`
echo "$pdbs pdb files to process"

set cpus = `echo $cpu_count | awk '{for(i=1;i<=NF;++i)sum+=$i;print sum}'`

echo "$cpu_count $machines" |\
 awk '{for(i=1;i<=NF/2;++i)for(c=1;c<=$i;++c){++n;printf("%03d %s\n",n,$(i+NF/2))}}' |\
cat >! cpu_map.txt

foreach cpu ( `awk '{print $1}' cpu_map.txt` )
   if("$machine" == "localhost") then
     mkdir -p /dev/shm/${USER}/cpu${cpu}_$$_/
     rm -f ./cpu$cpu
     if(-d ./cpu$cpu) rm -rf ./cpu$cpu
     ln -sf   /dev/shm/${USER}/cpu${cpu}_$$_/ ./cpu$cpu
     ln -sf `pwd` ./cpu${cpu}/back
   else
     mkdir -p cpu$cpu
     ln -sf ../ back
   endif
   awk -v cpu=$cpu -v cpus=$cpus 'NR%cpus==cpu-1{print "back/"$1}' $pdbfiles >! cpu${cpu}/pdbfiles.txt
   if("$MD_mult" != "") then
       echo "$MD_mult" >! cpu${cpu}/MD_mult.txt
   endif
   if("$GRID" != "") then
       echo "$GRID" >! cpu${cpu}/grid.txt
   endif
end

rm cpu*/*map >& /dev/null

set pwd = `pwd`
onintr killall
foreach cpu ( `awk '{print $1}' cpu_map.txt` )
  cd $pwd
  set machine = `awk -v cpu=$cpu '$1==cpu{print $2}' cpu_map.txt`
  sleep 0.1
    if("$machine" == "localhost") then
        ( cd ${pwd}/cpu${cpu} ; ${scriptdir}/md2map.com pdbfiles.txt $SG B=$B $opts |& tee md2map.log ) | grep expect &
    else
      rsh -n $machine "cd ${pwd}/cpu${cpu} ; ${scriptdir}/md2map.com pdbfiles.txt $SG B=$B $opts |& tee md2map.log" | grep expect &
  endif
end
wait

if( $?GOSLOW ) then
echo "moving final results into ./maps/"
mkdir -p ./maps/
mv cpu*/maps/sfall*.map ./maps/
endif

echo "summing sums..."
rm -f sum.map
foreach summap ( cpu*/sum.map )
  if(-e sum.map) then
     rm -f new.map
     echo "maps add" |\
     mapmask mapin1 sum.map mapin2 $summap mapout new.map > /dev/null
     mv new.map sum.map
  else
     cp $summap sum.map
  endif
end

# get number of ASUs in a unit cell
set symops = `awk -v SG=$SG '$4==SG{print $2;exit}' ${CLIBD}/symop.lib`
echo "$symops ASU/cell in $SG"

echo "computing average"
set ASUs = `echo $pdbs $MD_mult $symops | awk '{print $1*$2*$3*$4*$5}'`
echo "$pdbs x $MD_mult x $symops = $ASUs ASUs"
set scale = `echo $ASUs | awk '{print 1/$1}'`
echo "scale factor $scale 0"
echo "scale factor $scale 0" |\
mapmask mapin sum.map mapout avg.map > /dev/null
echo | mapdump mapin avg.map | egrep dens

echo "sfall"
rm -f md_map_avg.mtz
mapmask mapin avg.map mapout sfallme.map << EOF > /dev/null
xyzlim cell
AXIS Z X Y
EOF
sfall mapin sfallme.map hklout sfall.mtz << EOF | tee sfall.log > /dev/null
MODE SFCALC MAPIN
SFSG 1
EOF
sftools << EOF > /dev/null
read sfall.mtz
set labels
FP
PHI
calc Q COL SIGFP = 0.1
calc W COL FOM = 0.9
write md_map_avg.mtz col FP SIGFP PHI FOM
y
stop
EOF
rm -f sfallme.map sfall.mtz
rm -f pdbfiles$$.txt

cp avg.map avg_${pdbs}_${SG}.map
cp md_map_avg.mtz md_map_avg_${pdbs}_${SG}.mtz 

cad hklin1 md_map_avg.mtz hklout md_map_Bnorm.mtz << EOF > /dev/null
scale file 1 1 -$B
labin file 1 all
EOF
ls -l md_map_Bnorm.mtz

exit

./md2map.com $SG $pdbfiles | tee md2map.log | awk '/adding/{printf "."} END{print ""}'
    
cp avg.map avg_${realnum}_${SG}.map
cp md_map_avg.mtz md_map_avg_${realnum}_${SG}.mtz 

cad hklin1 md_map_avg.mtz hklout md_map_Bnorm.mtz << EOF
scale file 1 1 -$B
labin file 1 all
EOF

rm -f $pdbfiles

exit

killall:
foreach machine ( $machines )
   rsh -n $machine "killall -g md2map.com" >& /dev/null
end
