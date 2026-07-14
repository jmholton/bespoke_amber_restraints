#! /bin/tcsh -f
#
#   quick chirality / geometry / omega sanity check on an amber rst7 (or nc) file
#
#   if xtal.prmtop does not match the provided rst7 (eg. supercell vs ASU),
#   rst2pdb_runme.com is run to build resized.parm7 from padded.parm7, and that
#   is used instead.
#
#
set rstfile   = ""
set parmfile  = xtal.prmtop
set paddedparm = padded.parm7
set orignames = orignames.pdb

# built by rst2pdb_runme.com if xtal.prmtop does not match the rst7
set newparm  = resized.parm7
set neworig  = new_orignames.pdb

set quiet = 0
set debug = 0
set tempfile = /dev/shm/${USER}/tempfile_qgc$$_

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
      if("$key" == "rst") set rstfile = "$Val"
      if("$key" == "parm") set parmfile = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set orignames = $Arg
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

if( $tempfile =~ /dev/shm/$USER/* ) mkdir -p /dev/shm/$USER/

if(! -e "$rstfile") then
   set BAD = "usage: $0 amber.rst7 [$parmfile] [$orignames]"
   goto exit
endif

# fall back to the rst7's own directory for the support files
set dirname = `dirname $rstfile`
if(! -e "$parmfile") set parmfile = ${dirname}/$parmfile
if(! -e "$orignames") set orignames = ${dirname}/$orignames
if(! -e "$paddedparm") set paddedparm = ${dirname}/$paddedparm

if(! $quiet) then
cat << EOF
rstfile = $rstfile
parmfile = $parmfile
orignames = $orignames
EOF
endif

set t = $tempfile
set resized = 0

runchecks:
rm -f chircheck.txt geocheck.txt omega.txt omegaline.txt
cpptraj -p $parmfile << EOF >&! checks.log
trajin $rstfile lastframe
checkchirality chir out chircheck.txt
strip :WAT,HOH@Y1,EPW
check reportfile geocheck.txt
multidihedral omega omega out omegaline.txt range360
EOF
if( $status && ! $resized ) then
  set resized = 1
  echo "$parmfile does not match $rstfile"
  echo "running rst2pdb_runme.com to build $newparm from $paddedparm ..."
  rst2pdb_runme.com $rstfile parmfile=$parmfile paddedparm=$paddedparm >&! ${t}rst2pdb.log
  if(! -e "$newparm") then
    set BAD = "rst2pdb_runme.com did not produce $newparm ( see ${t}rst2pdb.log )"
    goto exit
  endif
  set parmfile = $newparm
  echo "using resized parm: $parmfile"
  if(-e "$neworig") set orignames = $neworig
  goto runchecks
endif

set ninv = `awk '$2!=1 && $1+0>0' chircheck.txt | wc -l`
set worstchir = `awk '$2!=1 && $1+0>0{print $1+0;exit}' chircheck.txt `
cat omegaline.txt |\
 awk 'NR==1{for(i=2;i<=NF;++i){split($i,w,":");oresnum[i]=w[2]};next}\
   {for(i=2;i<=NF;++i){dev=sqrt(($i-180)^2);print oresnum[i],dev,$i}}' >! omega.txt
set worstomega = `sort -k2gr omega.txt | head -n 1`
set ncis = `awk '$2>90{print}' omega.txt | wc -l`
set geoproblems = `cat geocheck.txt | wc -l`
set worstgeo = `awk '{print;exit}' geocheck.txt`
echo "$ninv inverted chiral centers, $ncis cis peptides ( $worstomega[2] deg) and $geoproblems geometry problems"


echo $worstomega |\
cat - $orignames |\
awk 'NR==1{worstres=$1;next}\
  {o=$NF} printed[o]{next}\
  o==worstres{printf("%s - ", substr($0,18,9));++printed[o]}\
  o==worstres+1{print substr($0,18,11),"   omega    ",o;++printed[o]}\
  o>worstres+1{exit}'
echo $worstchir |\
cat - $orignames |\
awk 'NR==1{worst=$1;next}\
  {o=$NF} printed[o]{next}\
  o==worst{print substr($0,18,11),"  chiral  ",o;++printed[o]}\
  o>worst{exit}'
echo $worstgeo |\
cat - $orignames |\
awk 'NR==1{geo=$0;++bad[$4+0];++bad[$6+0];next}\
  /^ATOM|^HETAT/{++a;o=$NF} printed[o]{next}\
  bad[a]{print substr($0,18,11),"    atom:",a,"res:",o;++printed[o]}\
  printed[o]>2{exit}'
echo $worstgeo

if( "$rstfile" =~ *.nc ) then

echo "extracting omegas..."
cpptraj -p $parmfile -y $rstfile << EOF >! cpptraj.txt
strip :WAT
strip @H=
multidihedral omega omega out omega.dat range360
EOF
echo "frame resnum dev omega"
cat omega.dat |\
 awk 'NR==1{for(i=2;i<=NF;++i){split($i,w,":");oresnum[i]=w[2]};next}\
   {++f;for(i=2;i<=NF;++i){dev=sqrt(($i-180)^2);if(dev>90)print f,oresnum[i],dev,$i;\
    if(dev>worst){worst=dev;worstf=f;worstres=oresnum[i]}}}\
   END{print worstf,worstres,worst,"worstomega"}' |\
head | tee first.txt
set oresnum = `awk '{print $2;exit}' first.txt`
awk -v oresnum=$oresnum 'NR==1{for(i=2;i<=NF;++i){split($i,w,":");idx[w[2]]=i};i=idx[oresnum];next}\
   {++f;{dev=sqrt(($i-180)^2);print f,dev,$i}}' omega.dat |\
cat >! worstomega_plot.txt
echo "plot of worst omega ($oresnum) vs frame in worstomega_plot.txt"

endif


exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
   rm -f ${t}* >& /dev/null
endif

exit

awk '/^ATOM|^HETAT/{print substr($0,22,1),substr($0,23,5),$NF}' orignames.pdb |\
 sort -u | sort -k3g >! ordresnums.txt

awk '/^TORS/ && $8=="CA" && $26=="CA"' tempfile_fullgeo.txt >! raw_omegas.txt

cat ordresnums.txt raw_omegas.txt |\
awk -v modulo=$modulo 'NF==3{ordresnum[$1,$2]=$3;next}\
  {o=ordresnum[$11,$12];\
    monomer=int(o/modulo);cpp=o-monomer;\
    print o,cpp,$10,$11,$12,monomer,$4}' |\
cat >! omega_count.txt

cat omega.txt omega_count.txt |\
awk 'NF==1{cppomega[NR]=$1;next}\
  {print $0,cppomega[$2]}' |\
tee compareme.txt


