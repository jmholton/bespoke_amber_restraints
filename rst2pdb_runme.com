#! /bin/tcsh -f
#
# quick converter from rst7 file to a pdb with original names
#
#
set rstfile  = ""
set frame = 1
set parmfile  = xtal.prmtop
set paddedparm  = padded.parm7
set orignames = orignames.pdb
set Bfactors = Bfac.pdb
set outprefix = ""
set pdbfile = ""

# will be created from paddedparm and orignames if needed
set newparm = resized.parm7
set neworig = new_orignames.pdb

set make_new_orig = 0

set quiet = 0
set debug = 0
set tempfile = /dev/shm/${USER}/tempfile_r2p_$$_

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
      if("$key" == "output") set outprefix = "$Val"
      if("$key" == "outfile") set pdbfile = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
      if("$Arg" =~ *.rst7 ) set rstfile = $Arg
      if("$Arg" =~ *.nc ) set rstfile = $Arg
      if("$Arg" =~ *.rst ) set rstfile = $Arg
      if("$Arg" =~ *.crd ) set rstfile = $Arg
      if("$Arg" =~ *.parmtop ) set parmfile = $Arg
      if("$Arg" =~ *.prmtop ) set parmfile = $Arg
      if("$Arg" =~ *.parm7 ) set parmfile = $Arg
    endif
    if("$arg" == "debug") set debug = "1"
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_r2p_
if( $tempfile =~ /dev/shm/$USER/* ) mkdir -p /dev/shm/$USER/

if(! -e "$rstfile") then
   set BAD = "usage: $0 amber.rst7 [new.pdb] [$parmfile] [$orignames]"
   goto exit
endif

set dirname = `dirname $rstfile`
if(! -e "$parmfile") set parmfile = ${dirname}/$parmfile
if(! -e "$orignames") set orignames = ${dirname}/$orignames



set prefix = `basename $rstfile .rst7`
set prefix = `basename $prefix .rst`
set prefix = `basename $prefix .crd`
set prefix = `basename $prefix .nc`

if("$outprefix" == "") then
  set outprefix = ${prefix}
  if( $frame != 1 ) then
    set outprefix = ${outprefix}_f${frame}
  endif
endif
if("$pdbfile" != "") set outprefix = `echo $pdbfile | awk '{gsub(".pdb$","");print}'`

if(! $quiet) then
cat << EOF
rstfile = $rstfile
parmfile = $parmfile
orignames = $orignames

tempfile = $tempfile
outfile = ${outprefix}.pdb
EOF
endif

set t = $tempfile

# just in case xtal.prmtop does not exist
if(! -e "$parmfile" ) then
  echo "WARNING: missing $parmfile"
  if( -e "$newparm") cp "$newparm" "$parmfile"
  if( -e "$paddedparm") cp "$paddedparm" "$parmfile"
endif

echo "reading $rstfile using $parmfile"
cpptraj -p $parmfile << EOF >&! ${t}cpptraj.log
trajin $rstfile $frame $frame
outtraj ${t}.pdb include_ep sg "P 1"
go
EOF
if( $status ) then
  set BAD = "cpptraj failed and no paddedparm: $paddedparm"
endif
if( $?BAD && -e "$paddedparm") then
  unset BAD
  set  RESIZE
  echo "making new parm file from $paddedparm"
  mv ${t}cpptraj.log ${t}cpptraj_error1.log
  set rstatoms = `awk '/Error: Number of atoms in /{gsub("[)(]","");for(i=NF;i>3;--i)if($i+0>0)print $i;exit}' ${t}cpptraj_error1.log | head -n 1`
  if( "$rstatoms" == "" ) then
   set rstatoms = `awk '/Error: Number of atoms in /{gsub("[)(]","");print $(NF-2);exit}' ${t}cpptraj_error1.log`
  endif
  set maxatoms = `echo list | cpptraj -p $paddedparm | awk '$4=="atoms,"{print $3}' | head -n 1`
  set stripmask = `echo $rstatoms $maxatoms | awk '{print $1+1"-"$2}'`
  set parmfile = $newparm
  echo "new parmfile: $newparm"
  cpptraj -p $paddedparm << EOF >&! ${t}strip1.log
  parmstrip @$stripmask
  parmwrite out $newparm
EOF
  cpptraj -p $parmfile << EOF >&! ${t}cpptraj.log
  trajin $rstfile $frame $frame
  outtraj ${t}.pdb include_ep sg "P 1"
  go
EOF
  if( $status ) then
    set BAD = "conversion failed"
    goto exit
  endif
  set newres = `awk '$6=="res," && $5+0>0{print $5}' ${t}cpptraj.log | head -n 1`
  set origres = `awk '/^ATOM|^HETAT/{print substr($0,17,12)}' $orignames | sort -u | wc -l`
  echo "$newres residues in $parmfile"
  echo "$origres residues in $orignames"
  # we need to pad/strip orignames
  if( $newres > $origres ) then
    echo "padding and re-counting $orignames -> $neworig"
    egrep "^CRYST|^ATOM|^HETAT|^SSBO|^LINK" $orignames >! ${t}_orig.pdb
    egrep "^ATOM|^HETAT" ${t}.pdb |\
     convert_pdb.awk -v renumber=ordinal,w4,watS,chain,chainrestart -v fixEe=1 \
      -v append=ordresnum |\
     awk -v l=$origres '$NF>l' >> ${t}_orig.pdb
    awk '{print substr($0,1,80)}' ${t}_orig.pdb |\
      convert_pdb.awk -v renumber=ordinal,w4,watS,chain,chainrestart -v fixEe=1 \
      -v append=ordresnum >! $neworig
    set orignames = $neworig
    set origres = `awk '/^ATOM|^HETAT/{print substr($0,17,12)}' $orignames | sort -u | wc -l`
    echo "$origres residues in $orignames"
  endif
  if( $newres < $origres && $make_new_orig ) then
    echo "stripping and re-counting $orignames -> $neworig"
    egrep "^CRYST|^ATOM|^HETAT|^SSBO|^LINK" $orignames |\
    awk '{print substr($0,1,80)}' |\
      convert_pdb.awk -v renumber=ordinal,w4,watS,chain,chainrestart -v fixEe=1 \
      -v append=ordresnum |\
    awk -v l=$newres '/^ATOM|^HETAT/ && $NF>l{exit} {print}' >! $neworig
    set orignames = $neworig
    set origres = `awk '/^ATOM|^HETAT/{print substr($0,17,12)}' $orignames | sort -u | wc -l`
    echo "$origres residues in $orignames"
  endif
endif

egrep "^SSBOND|^LINK|^CRYST" $orignames >! ${t}out.pdb
# take coordinates only, using names from starting point
awk '/^ATOM|^HETAT/{print $0,"ORIG"}' $orignames |\
cat - ${t}.pdb |\
awk '$NF=="ORIG"{++o;pre[o]=substr($0,1,30);post[o]=substr($0,55,length($0)-55-4);next}\
  ! /^ATOM|^HETAT/{next}\
      {++n}\
      pre[n]==""{print "REMARK WARNING atom",n,"missing from orignames.pdb";\
       pre[n]=substr($0,1,30);post[n]=substr($0,55)}\
      {printf("%s%s%s\n",pre[n],substr($0,31,24),post[n])}' |\
cat >> ${t}out.pdb

grep "WARNING" ${t}out.pdb
if( ! $status ) set MISSING

if(-e "$Bfactors" ) then
    echo "applying B factors from $Bfactors"
    egrep "^SSBOND|^LINK|^CISP|^CRYST" $Bfactors >! ${t}Bfac.pdb

    set test = `egrep "^CRYST" ${t}Bfac.pdb | wc -l`
    if( ! $test ) then
       echo "WARNING: cell missing from $Bfactors "
       egrep "^CRYST" ${t}out.pdb | head -n 1 >> ${t}Bfac.pdb
    endif

    # take occ/B only, using names from starting point
    awk '/^ATOM|^HETAT/{print $0,"BFAC"}' $Bfactors |\
    cat - ${t}out.pdb |\
    awk '$NF=="BFAC"{++o;occB[o]=substr($0,55,12);next}\
      ! /^ATOM|^HETAT/{next}\
          {++n;pre=substr($0,1,54);post=substr($0,67)}\
          occB[n]==""{print "REMARK WARNING atom",n,"missing from Bfac, max="o;\
           occB[n]=last_occB}\
          {printf("%s%s%s\n",pre,occB[n],post);last_occB=occB[n]}' |\
    cat >> ${t}Bfac.pdb

    grep "WARNING" ${t}Bfac.pdb
    if( ! $status ) set MISSING

    cp ${t}Bfac.pdb ${outprefix}.pdb
else
    cp ${t}out.pdb ${outprefix}.pdb
endif


if( $?RESIZE || $?MISSING ) then
echo "suggestions:"
if( $?MISSING ) echo "use padded_orignames.pdb"
if(-e "$neworig" && "$neworig" != "orignames.pdb") echo "mv $neworig orignames.pdb"
if(-e "$newparm" && "$newparm" != "xtal.prmtop") echo "mv $newparm xtal.prmtop"
endif

ls -l ${outprefix}.pdb

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



foreach rst7 ( `ls -1rt *.rst7` )

set basename = `basename $rst7 .rst7`
if(-e ${basename}.pdb) continue
srun rst2pdb_runme.com $rst7 > /dev/null &

end



set parmfile = padded.parm7

set goalatoms = `echo list | cpptraj -p $parmfile | awk '$4=="atoms,"{print $3}' | head -n 1`

egrep "^CRYST1" refme.pdb | head -n 1 >! zero.pdb

echo $goalatoms |\
awk '{goal=$1;for(i=1;i<=goal;++i){\
  print "ATOM      1  N   GLY A   1        0.000   0.000   0.000  1.00  0.00";\
  }}' >> zero.pdb

cpptraj -p $parmfile -y zero.pdb << EOF
outtraj renumberme.pdb include_ep sg "P 1"
EOF
cat renumberme.pdb |\
convert_pdb.awk -v renumber=ordinal,w4,watS,chain,chainrestart -v fixEe=1 \
  -v append=ordresnum >! padded_orignames.pdb

append_file_date.com orignames.pdb

# cp padded_orignames.pdb orignames.pdb
