#! /bin/tcsh -f
#
#  find atoms that are most challenging to restraints
#  delete them, now and forever
#
#
set tootight  = 9
set toofar = 1.0
set maxE = 50
set maxfirstE = 20

set maxbad    = 10

set challengefile = restraints_challenge_vs_weight.txt
set allpos       = all_possible_refpoints.pdb
set restraintpdb = current_restraints.pdb
set outfile      = new_restraints.pdb

set outprefix = "new_"

set debug = 0

set tempfile = /dev/shm/${USER}/tempfile_rwr_$$_

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if( $assign ) then
      # re-set any existing variables
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "restpdb") set restraintpdb = "$Val"
      if("$key" == "refpoints") set restraintpdb = "$Val"

      if("$key" == "multsfile") set multsfile = "$Val"
      if("$key" == "output") set outfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set restraintpdb = $Arg
    endif
    if("$key" == "debug") set debug = "$Val"
end

if(! -e "$restraintpdb") then
    set BAD = "no coordinates provided"
    goto exit
endif

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_rwr_
set tmpdir = `dirname $tempfile`

set t = "$tempfile"


sort -gr $challengefile >! ${t}challenges.txt

echo "$tootight $toofar $maxE $maxfirstE" |\
cat - ${t}challenges.txt |\
awk 'NR==1{tootight=$1;toofar=$2;maxE=$3;maxfirstE=$4;next}\
  $3>tootight && $2>toofar || $1>maxE || NR==1 && $1>maxfirstE' |\
head -n $maxbad |\
tee ${t}challenges.txt

cat ${t}challenges.txt |\
awk '{l=length($0);print substr($0,l-16)"   BAD"}' |\
cat >! ${t}bad_refpoints.txt


foreach restraintfile ( $allpos $restraintpdb )

cat ${t}bad_refpoints.txt $restraintfile |\
awk '$NF=="BAD"{id=substr($0,1,15);++bad[id];next}\
  {id=substr($0,12,15)}\
  ! bad[id]{print}' |\
cat >! ${t}culled.pdb
set changes = `diff $restraintfile ${t}culled.pdb | egrep "^<" | wc -l`
echo "deleted $changes refpoints from $restraintfile"
mv ${t}culled.pdb ${outprefix}${restraintfile}
ls -l ${outprefix}${restraintfile}
end



