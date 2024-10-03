#! /bin/tcsh -f
#
#  average atomic B factors with immediate neighbors in space                     -James Holton 2-25-24
#
#
set pdbfile    = current_restraints.pdb
set outfile    = thrubondBavg_restraints.pdb

# maximum distance to still be restrained
set maxdist = 3.2
set mindist = 0.1
# fraction of particular weight to transport into neighbors
set avgfac = 0.5
# if oldB>newB, fraction to keep
set spread = 0
# flag to debug things
set debug = 0

set path = ( `dirname $0` $path )

set tempfile = tempfile_$$_
set tempfile = /dev/shm/${USER}/tempfile_tbB_$$_


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
      if("$key" == "output") set outfile = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
      if("$key" == "maxdist") set maxdist = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
    endif
    if("$arg" == "debug") set debug = "1"
end

#if( $debug ) set tempfile = ./tempfile
mkdir -p `dirname $tempfile` > /dev/null

cat << EOF
pdbfile  $pdbfile
outfile  $outfile

maxdist  $maxdist
mindist  $mindist
avgfac   $avgfac

tempfile  $tempfile
EOF

if(! -r "$pdbfile") then
  set BAD = "cannot read $pdbfile"
  goto exit
endif

set t = "$tempfile"

echo "converting atom names"
# make sure there is space group?
egrep "^CRYST" $pdbfile        | head -n 1 >! ${t}temp.pdb
echo "MODEL        1"                      >> ${t}temp.pdb
awk '{print substr($0,1,80)}' $pdbfile |\
convert_pdb.awk -v only=atoms -v skip=H,EP -v fixEe=1 -v append=origid |\
cat >> ${t}temp.pdb
echo "MODEL        2"                      >> ${t}temp.pdb
awk '{print substr($0,1,80)}' $pdbfile |\
convert_pdb.awk -v only=atoms -v skip=H,EP -v fixEe=1 -v append=origid |\
cat >> ${t}temp.pdb

convert_pdb.awk -v renumber=w4,ordinal,chain ${t}temp.pdb >! ${t}compme.pdb

echo "measuring distances..."
ncont xyzin ${t}compme.pdb << EOF >! ${t}ncont.log
source /1/*/*/*:*
target /2/*/*/*:*
cells 0
mindist $mindist
maxdist $maxdist
sort dist
EOF

echo "formatting name pairs..."
cat ${t}ncont.log |\
awk -F ":" -v debug=1 -v maxdist=$maxdist 'NF!=3 || $NF!~/[XYZ]/{next}\
   {dist=substr($0,57,7)}\
   dist<=maxdist{if(debug)print;\
   chain1=substr($0,5,1);\
   resnum1=substr($0,7,4) substr($0,17,1);\
   typ1=substr($0,12,3);\
   atom1=substr($0,19,4);\
   Ee1=substr($0,24,2);\
   conf1=substr($0,28,1);\
   chain2=substr($0,33,1);\
   resnum2=substr($0,35,4) substr($0,45,1);\
   typ2=substr($0,40,3);\
   atom2=substr($0,47,4);\
   Ee2=substr($0,52,2);\
   conf2=substr($0,56,1);\
   if(debug)print " /1/"chain1"/"resnum1"("typ1")./"atom1"["Ee1"]:"conf1;\
   if(debug)print "atom1=",atom1,"typ1=",typ1,"chain1=",chain1,"resnum1=",resnum1;\
   # left justify, like refmac \
   gsub(" ","",typ1);typ1=sprintf("%-3s",typ1);\
   gsub(" ","",typ2);typ2=sprintf("%-3s",typ2);\
   id1= atom1 conf1 typ1" "chain1 resnum1;\
   id2= atom2 conf2 typ2" "chain2 resnum2;\
   print id1 "|" id2 "|  " dist,"DIST";}' |\
cat >! ${t}encoded_distances.txt

# change back to original names
cat ${t}compme.pdb ${t}encoded_distances.txt |\
awk -v debug=$debug '/^MODEL/{model=$NF}\
   /^ATOM|^HETAT/{id=substr($0,13,15);\
   origid[id]=substr($0,index($0,"|")+2,15)}\
   $NF=="DIST"{id1=substr($0,1,15);id2=substr($0,17,15);\
      dist=sprintf("%6.2f",substr($0,34)+0);\
    print origid[id1]" |"origid[id2]" |"dist,"DIST"}' |\
cat >! ${t}decoded_distances.txt


# now determine which atoms to restrain to what:
# protein must have same name, solvent does not
echo $avgfac $spread $debug |\
cat - ${t}decoded_distances.txt $pdbfile |\
awk 'NR==1{avgfac=$1;spread=$2;debug=$3;next}\
      $NF=="DIST"{\
      id1=substr($0,1,15);id2=substr($0,17,15);\
      atom1=substr(id1,1,4);atom2=substr(id2,1,4);\
      typ1=substr(id1,6,3);typ2=substr(id2,6,3);\
      resid1=substr(id1,10);resid2=substr(id2,10);\
      d=sprintf("%6.2f",substr($0,32)+0);\
      if(! bonded[id1,id2]){\
        ++bonded[id1,id2];++bonded[id2,id1];\
        partners[id1]=partners[id1]"|"id2;\
        partners[id2]=partners[id2]"|"id1;\
      }\
      if(debug)print "REMARK",id1 "|" partners[id1]}\
  ! /^ATOM|^HETAT/{next}\
  {id=substr($0,13,15);\
   Bfac[id]=substr($0,61,6)+0;\
   ++n;line[n]=$0;\
  }\
  END{\
    for(i=1;i<=n;++i){\
      id=substr(line[i],13,15);\
      pre=substr(line[i],1,60);post=substr(line[i],67);\
      avgB=0;\
      np=split(partners[id],w,"|");\
      if(debug)print "REMARK",id "| i=",i,"np=",np,"B=",Bfac[id],partners[id];\
      for(j=2;j<=np;++j){\
        id2=w[j];\
        avgB+=Bfac[id2]/np;\
        if(debug)print "REMARK id1=",id,"id2=",id2,"B=",Bfac[id2];\
      }\
      newB=Bfac[id];\
      oldwt=(1-avgfac);\
      if(oldwt<0)oldwt=0;\
      if(np) newB=(oldwt*newB+avgfac*avgB)/(oldwt+avgfac);\
      if(spread>0 && newB<Bfac[id]) newB=spread*Bfac[id]+(1-spread)*newB;\
      printf("%s%6.2f%s\n",pre,newB,post);\
    }\
  }' |\
cat >! ${t}avged.pdb


# make sure they are in the same order as input pdb file
combine_pdbs_runme.com ${t}avged.pdb $pdbfile \
  outfile=$outfile >! ${t}combine.log


wc $outfile



exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
rm ${t}* > /dev/null

exit

