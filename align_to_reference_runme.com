#! /bin/tcsh -f
#
#  use SG and super_mult to expand a PDB and mtz into a supercell
#

#defaults
set outprefix = "aligned"
set rstfiles = ( )
set parmfile = xtal.prmtop
set pdbfile = "align_ref.pdb"
set threshold = "median"

set logfile = align_details.log

set pwd = `pwd`
set pdir = `dirname $0`
set path = ( $pdir $path )

set quiet = 0
set debug = 0
set tempfile = /dev/shm/${USER}/tempfile_a2r_$$_

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
          if( ! $quiet ) echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "reference") set pdbfile = "$Val"
      if("$key" == "output") set outprefix = "$Val"
      if("$key" == "outfile") set outprefix = `basename $Val .pdb`
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
      if("$Arg" =~ *.parm7 ) set parmfile = $Arg
      if("$Arg" =~ *.prmtop ) set parmfile = $Arg
      if("$Arg" =~ *.rst7 ) set rstfiles = ( $rstfiles $Arg )
      if("$Arg" =~ *.nc ) set rstfiles = ( $rstfiles $Arg )
    endif
    if("$arg" == "debug") set debug = "1"
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_a2r_
if( $tempfile =~ /dev/shm/$USER/* ) mkdir -p /dev/shm/$USER/

set t = ${tempfile}


foreach file ( $pdbfile $rstfiles $parmfile )
  if(! -e "$file" && "$file" != "" ) then
     echo "WARNING: $file does not exist"
  endif
end
if(! -e "$pdbfile" || $#rstfiles == 0 ) then
   set BAD = "usage: $0 amber.nc refpoints.pdb threshold=$threshold"
   goto exit
endif


if( ! $quiet ) then
cat << EOF 
outprefix = $outprefix
pdbfile = $pdbfile
rstfiles = $rstfiles

tempfile = $tempfile
debug = $debug
EOF
endif

echo "generating align_ref.crd"
update_centroid_positions_runme.com \
        topfile=$parmfile \
        pdbfile=$pdbfile \
        outfile=${t}align_ref.crd >> $logfile

echo "extracting"
rst2pdb_runme.com ${t}align_ref.crd ${t}this.pdb >> $logfile


echo "combining"
cat $pdbfile |\
awk 'substr($0,61,6)+0>0{print substr($0,1,80),"     SEL"}' |\
cat >! ${t}labeled.pdb
combine_pdbs_runme.com ${t}labeled.pdb ${t}this.pdb printref=1 outfile=${t}new.pdb >> $logfile
egrep "^ATOM|^HETAT" ${t}new.pdb |\
awk '{++n} $NF=="SEL"{print n,substr($0,61,6)}' |\
cat >! ${t}align_atom_weights.txt

set w = $threshold
if( $threshold == median ) then
    set h = `cat ${t}align_atom_weights.txt | wc -l | awk '{print int($1/2)}'`
    set w = `awk '{print $2}' ${t}align_atom_weights.txt | sort -g | head -n $h | tail -n 1`
    echo "median threshold: $w"
endif
awk -v w=$w '$2>w' ${t}align_atom_weights.txt |\
awk 'NR==1{s=e=$1;next}\
      $1==e+1{e=$1;next} \
      {print s"-"e;s=e=$1}\
      END{print s"-"e}' |\
awk -F "-" '$1==$2{print $1;next} {print}' >! ${t}align_ranges.txt 
set align_ranges = `cat ${t}align_ranges.txt `
set align_mask = `echo $align_ranges | awk '{gsub(" ",",");print}'`

if( "$align_mask" == "" ) then
  set BAD = "unable to assing alignment mask"
  goto exit
endif

foreach rstfile ( $rstfiles )

echo "aligning $rstfile to $pdbfile atoms with B > $w"

set ext = `echo $rstfile | awk -F "." '{print $NF}'`

cpptraj -p $parmfile -y $rstfile -c ${t}align_ref.crd << EOF >> $logfile
rmsd rmsd reference norotate @$align_mask out rmsd.txt savevectors combined vecsout vecsout.txt
trajout ${outprefix}.$ext
EOF
#cat rmsd.txt vecsout.txt >> align_${itr}.log

echo -n "final shift: "
tail -n 1 vecsout.txt | awk '{print $2,$3,$4,"(",sqrt($2*$2+$3*$3+$4*$4),")"}'
echo -n "biggest shift: "
cat vecsout.txt | awk '{print $2,$3,$4,"(",sqrt($2*$2+$3*$3+$4*$4),")"}' | sort -k5g | tail -n 1

if(! $quiet) ls -l ${outprefix}.$ext

end

exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
   rm -f ${t}* > /dev/null
endif

exit

