#! /bin/csh -f
#
#        Script for calculating PHI-PSI-CHI values from a PDB file
#        (using geomcalc)
#

set pdbin = pdbfile.pdb
set tempfile = phipsichi_$$temp

set debug = 0

if($#argv == 0) goto Help
goto Setup
Help:
################################################################################
cat << EOF
usage: phipsichi.com protein.pdb

where:
 protein.pdb   - is the structure you want phi, psi, and chi angles from

EOF
################################################################################
exit 9
Return_from_Setup:


# simplify
awk 'substr($0,78,1)!="H" && ! /HOH/' $pdbin >! ${tempfile}pdbin.pdb

# get the atom name list
cat ${tempfile}pdbin.pdb |\
nawk '/^ATOM|^HETAT/{chain=substr($0, 22, 1); if(chain==" ") chain = "_";\
    print chain, substr($0, 23, 4), substr($0, 13, 4)}' |\
cat >! ${tempfile}atoms
# format: chain resum atomname

# store residue types
cat ${tempfile}pdbin.pdb |\
nawk '/^ATOM|^HETAT/{chain=substr($0, 22, 1); if(chain==" ") chain = "_";\
    print chain, substr($0, 23, 4), substr($0, 18, 3)}' |\
awk '! seen[$0]{print;++seen[$0]}' >! ${tempfile}types

# get distances?
cat ${tempfile}atoms |\
nawk '$1=="_"{$1=" "}\
      $3~/^[CN]$/ || $3=="CA"{print $1 "|" $2 ":" $3}' |\
nawk 'last{print "DIST", last, $0} {last=$0}' |\
geomcalc xyzin ${tempfile}pdbin.pdb |\
nawk '/Distance between atoms/{print $4, $6, $NF}' |\
cat >! ${tempfile}dists
# format: atom1 atom2 dist
# use this to identify chain breaks?


# get phi-psi dihedrals
cat ${tempfile}atoms |\
nawk '$1=="_"{$1=" "}\
      $3~/^[CN]$/ || $3=="CA"{print $1 "|" $2 ":" $3}' |\
nawk 'old && older && oldest{print "TORS", oldest, older, old, $0}\
     {oldest=older;older=old;old=$0}' |\
cat ${tempfile}dists - |\
awk '/^TORS/ && conn[$2,$3] && conn[$3,$4] && conn[$4,$5]{print;next}\
   $3<=2.1{++conn[$1,$2]}' |\
geomcalc xyzin ${tempfile}pdbin.pdb |\
nawk '/Torsion angle between/{gsub(","," "); print $5,$6,$7,$9, $NF}' |\
nawk '$3~/CA$/{print $3, "|phi", $NF}\
      $2~/CA$/{print $2, "|psi", $NF}\
      $1~/CA$/{print $1, "|omega", $NF}' |\
nawk 'BEGIN{FS="|"} {print $1, $2+0, $3}' |\
cat >! ${tempfile}phipsi


# find bonded atoms
cat ${tempfile}atoms |\
nawk '$1=="_"{$1=" "} {print $1 "|" $2 ":" $3}' |\
nawk -v maxprev=14 '{for(i=1;i<=maxprev;++i){\
    if(prev[i]!="") print "DIST", prev[i], $0};\
       for(i=maxprev;i>1;--i){prev[i]=prev[i-1]};\
       prev[1]=$0}' |\
geomcalc xyzin ${tempfile}pdbin.pdb |\
nawk '/Distance between atoms/{print $4, $6, $NF}' |\
nawk '$NF<2{print "",$1,$2}' |\
nawk '{gsub("[|:]"," "); print $1,$2,$3,"-",$4,$5,$6}' |\
cat >! ${tempfile}bonds

# get chi dihedrals
cat ${tempfile}bonds |\
nawk '$1==$5 && $2==$6 && $3!~/^[OC]$/ && $6!~/^[OC]$/{\
      print $1"|"$2":"$3, $5"|"$6":"$7}' |\
nawk '{print bound[bound[$1]], bound[$1], $1, $2} {bound[$2]=$1}' |\
nawk 'NF==4{print "TORS", $0}' |\
geomcalc xyzin ${tempfile}pdbin.pdb |\
nawk '/Torsion angle between/{gsub(","," "); print $5,$6,$7,$9, $NF}' |\
nawk '{opt=""} $4~/[1-9]$/{opt="." substr($4,length($4))}\
    {angle="";for(i=1;i<=4;++i){greek=substr($i,length($i)); \
     if(greek~/[1-9]/){greek=substr($i,length($i)-1,1)};\
     angle=angle greek;}} \
    angle=="NABG"{print $1, "|chi1" opt, $NF}\
    angle=="ABGD"{print $1, "|chi2" opt, $NF}\
    angle=="BGDE"{print $1, "|chi3" opt, $NF}\
    angle=="GDEZ"{print $1, "|chi4" opt, $NF}\
    angle=="DEZH"{print $1, "|chi5" opt, $NF}' |\
nawk 'BEGIN{FS="|"} {print $1, $2+0, $3}' |\
cat >! ${tempfile}chi

# get all bonded dihedrals
cat ${tempfile}bonds |\
nawk '$1==$5 && $2==$6' |\
nawk '{print $1"|"$2":"$3, $5"|"$6":"$7}' |\
nawk '{print bound[bound[$1]], bound[$1], $1, $2} {bound[$2]=$1}' |\
nawk 'NF==4{print "TORS", $0}' |\
geomcalc xyzin ${tempfile}pdbin.pdb |\
nawk '/Torsion angle between/{gsub(","," "); print $5,$6,$7,$9, $NF}' |\
nawk '{opt=""} $4~/[1-9]$/{opt="." substr($4,length($4))}\
    {angle="";for(i=1;i<=4;++i){greek=substr($i,length($i)); \
     if(greek~/[1-9]/){greek=substr($i,length($i)-1,1)};\
     angle=angle greek;} } \
    angle=="CNAC"{print $1, "|phi" opt, $NF}\
    angle=="TNAC"{print $1, "|phi" opt, $NF}\
    angle=="NACN"{print $1, "|psi" opt, $NF}\
    angle=="TTCN"{print $1, "|psi" opt, $NF}\
    angle=="NACT"{print $1, "|psi" opt, $NF}\
    angle=="ACNA"{print $1, "|omega" opt, $NF}\
    angle=="TTNA"{print $1, "|omega" opt, $NF}\
    angle=="NABG"{print $1, "|chi1" opt, $NF}\
    angle=="ABGD"{print $1, "|chi2" opt, $NF}\
    angle=="BGDE"{print $1, "|chi3" opt, $NF}\
    angle=="GDEZ"{print $1, "|chi4" opt, $NF}\
    angle=="DEZH"{print $1, "|chi5" opt, $NF}' |\
nawk 'BEGIN{FS="|"} {print $1, $2+0, $3}' |\
cat >! ${tempfile}chi


cat ${tempfile}phipsi ${tempfile}chi ${tempfile}types |\
nawk '{ID=$1 $2} \
      /phi/{phi[ID]=sprintf("%7.2f", $NF)}\
      /psi/{psi[ID]=sprintf("%7.2f", $NF)}\
     /omeg/{omega[ID]=sprintf("%7.2f", $NF)}\
     /chi1/&& ! chi1[ID]{chi1[ID]=sprintf("chi1= %7.2f", $NF)}\
     /chi2/&& ! chi2[ID]{chi2[ID]=sprintf("chi2= %7.2f", $NF)}\
     /chi3/&& ! chi3[ID]{chi3[ID]=sprintf("chi3= %7.2f", $NF)}\
     /chi4/&& ! chi4[ID]{chi4[ID]=sprintf("chi4= %7.2f", $NF)}\
     /chi5/&& ! chi5[ID]{chi5[ID]=sprintf("chi5= %7.2f", $NF)}\
 NF==3{if(phi[ID]=="") phi[ID]="-"; if(psi[ID]=="") psi[ID]="-"; \
       if(omega[ID]=="") omega[ID]="-";\
       if($3~/HIS|PHE|TYR|TRP/) chi3[ID]=chi4[ID]=chi5[ID]="";\
    printf "%s %s %4d phi= %7s psi= %7s omega= %7s %s %s %s %s %s\n",\
    $3, $1, $2, phi[ID], psi[ID], omega[ID], chi1[ID], chi2[ID], chi3[ID], chi4[ID], chi5[ID]}'


# clean up
if( $debug ) exit

rm -f ${tempfile}atoms >& /dev/null
rm -f ${tempfile}types >& /dev/null
rm -f ${tempfile}dists >& /dev/null
rm -f ${tempfile}phipsi >& /dev/null
rm -f ${tempfile}bonds >& /dev/null
rm -f ${tempfile}chi >& /dev/null
rm -f ${tempfile}pdbin.pdb >& /dev/null

exit


cat  ${tempfile}types - |\
nawk 'NF==3{type[$1 $2]=$3;next}\
      {printf "%s %s %4d %6s = %7.2f\n", type[$1 $2], $1, $2, $3, $4}' |\

exit

# get chi dihedrals
cat ${tempfile}atoms |\
nawk '$3!~/^[CO]$/{print $1 "|" $2 ":" $3}' |\
nawk 'old && older && oldest{print "TORS", oldest, older, old, $0}\
     {oldest=older;older=old;old=$0}' |\
geomcalc xyzin ${tempfile}pdbin.pdb |\
nawk '/Torsion angle between/{gsub(","," "); print $5,$6,$7,$9, $NF}' |\
nawk '$1~/N$/{print $1, "|chi1", $NF}\
      $1~/CA$/{print $1, "|chi2", $NF}\
      $1~/CB$/{print $1, "|chi3", $NF}\
      $1~/G$/{print $1, "|chi4", $NF}\
      $1~/D$/{print $1, "|chi5", $NF}\
      $1~/E$/{print $1, "|chi6", $NF}\
      $1~/Z$/{print $1, "|chi7", $NF}' |\
nawk 'BEGIN{FS="|"} {print $1, $2+0, $3}' |\
cat >! ${tempfile}chi



################################################################################
Setup:
nawk 'BEGIN{exit}' >& /dev/null
if($status) alias nawk awk

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
      if("$key" == "refpoints") set fullref = "$Val"
      if("$key" == "current") set notthese = "$Val"
      if("$key" == "maxmove") set maxmoves = "$Val"

      if("$key" == "output") set outfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ||("$arg" =~ *.brk)) then
        if(-e "$Arg") then
            set pdbin = "$Arg"
        else
            echo "WARNING: $Arg does not exist"
        endif
    endif
    if("$key" == "debug") set debug = 1
    if("$key" == "pdbfile") set pdbin = "$Val"
end
if(! -e "$pdbin") goto Help

goto Return_from_Setup
