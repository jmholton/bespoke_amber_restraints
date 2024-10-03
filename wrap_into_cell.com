#! /bin/tcsh -f
#
#   wrap-in waters and other small things that have wandered into neighboring cells
#
#

set pdbfile = ""
set outfile = wrapped.pdb

set doprotein = 0
set breakbonds = 0
set bychain = 0
set byter = 0

set logfile = wrap_details.log
set tempfile = /dev/shm/${USER}/suckin_$$_

foreach arg ( $* )
    set Key = `echo $arg | awk -F "=" '{print $1}'`
    set Val = `echo $arg | awk -F "=" '{print $2}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set val = `echo $Val | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if("$Key" =~ *.pdb && "$val" == "") then
        set pdbfile = $Key
    endif
    if("$key" == "pdbfile") set pdbfile = "$Val"
    if("$key" == "outfile") set outfile = "$Val"
    if("$key" == "output") set outfile = "$Val"
    if("$key" == "doprotein" && "$val" != "") set doprotein = "$Val"
    if("$key" == "notprotein") then
        set doprotein = 0
    endif
    if("$key" == "breakbonds" && "$val" != "") set breakbonds = "$Val"
    if("$key" == "protein") set doprotein = 1
    if("$key" == "doprotein") set doprotein = 1
    if("$key" == "bychain") set bychain = 1
    if("$key" == "byter") set byter = 1
    if("$key" == "logfile") set logfile = "$Val"
    if("$key" == "tempfile") then
        set tempfile = "$Val"
        setenv DEBUG
    endif
    if("$key" == "debug") then
        setenv DEBUG
    endif
end

touch $logfile
mkdir -p `dirname $tempfile`

cat << EOF 
pdbfile= $pdbfile
doprotein= $doprotein
breakbonds= $breakbonds
bychain= $bychain
byter= $byter
outfile= $outfile

logfile= $logfile
tempfile= $tempfile
debug= $?DEBUG
EOF

set CELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7;exit}' $pdbfile`

awk '/^CRYST1|^SSBOND/' $pdbfile | tee ${tempfile}.pdb >! ${tempfile}out.pdb
if( $doprotein ) then
    cp $pdbfile ${tempfile}.pdb
else
    egrep "^ATOM|^HETAT" $pdbfile |\
    awk '{typ=substr($0,18,3);protein=0}\
    typ~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein}\
    typ~/HID|HIE|HIS|HIP|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein}\
    protein{print;next}' |\
    cat >> ${tempfile}out.pdb
    egrep "^ATOM|^HETAT" $pdbfile |\
    awk '{typ=substr($0,18,3);protein=0}\
    typ~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein}\
    typ~/HID|HIE|HIS|HIP|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein}\
    ! protein{print;next}' |\
    cat >> ${tempfile}.pdb
endif

coordconv xyzin ${tempfile}.pdb xyzout ${tempfile}.xyz << EOF >> $logfile
CELL $CELL
INPUT PDB
OUTPUT FRAC
END
EOF

# rescue long residue numbers
echo "$doprotein $breakbonds $bychain $byter $?DEBUG" |\
cat - ${tempfile}.pdb |\
awk 'NR==1{print;next}\
     /^ATOM|^HETAT/{print ++n,substr($0,23,8)+0,"ORIG"}\
     /^TER/{print "TER",++terchain}' |\
cat - ${tempfile}.xyz |\
awk 'NR==1{doprotein=$1;breakbonds=$2;bychain=$3;byter=$4;debug=$NF;next}\
     $NF=="ORIG"{i=$1;resn[i]=$2;tc[i]=terchain+0;next}\
     /^TER/{terchain=$2;next}\
     {++n;atn[n]=$1;x[n]=$2;y[n]=$3;z[n]=$4;B[n]=$5;occ[n]=$6;Z[n]=$7;\
     oresn[n]=substr($0,56,10);\
     atom[n]=substr($0,66,5);\
     typ[n]=substr($0,71,3);chain[n]=substr($0,75,1);\
     cr=chain[n] resn[n] typ[n];\
     if(bychain)cr=chain[n];\
     if(byter)cr=tc[n];\
     if(breakbonds)cr=n;\
     dx=dy=dz=0;\
     while(x[n]+dx<0)++dx;\
     while(x[n]+dx>1)--dx;\
     while(y[n]+dy<0)++dy;\
     while(y[n]+dy>1)--dy;\
     while(z[n]+dz<0)++dz;\
     while(z[n]+dz>1)--dz;\
     if(Z[n]>1)++votes[cr,dx,dy,dz];\
     if(dx==0 && dy==0 && dz==0)votes[cr,dx,dy,dz]+=1e-6;\
     if(debug)print "DEBUG",atom[n],"voting for",cr,dx,dy,dz,"->",votes[cr,dx,dy,dz]}\
  END{for(str in votes){split(str,w,"\034");cr=w[1];dx=w[2];dy=w[3];dz=w[4];\
   if(debug)print "DEBUG:",cr,dx,dy,dz,"got",votes[str],"votes";\
   if(votes[str]>maxvote[cr]){maxvote[cr]=votes[str];\
    bestdxyz[cr]=dx" "dy" "dz};\
    if(debug)print "DEBUG: best dxyz of",cr,"->",bestdxyz[cr],"votes:",votes[str]}\
   for(i=1;i<=n;++i){\
    cr=chain[i] resn[i] typ[i];\
    if(bychain)cr=chain[i];\
    if(byter)cr=tc[i];\
    if(breakbonds)cr=i;\
    split(bestdxyz[cr],w);dx=w[1];dy=w[2];dz=w[3];\
printf("%5d%10.5f%10.5f%10.5f%10.5f%5.2f%5d%10d%-5s%3s %1s\n", \
       atn[i],x[i]+dx,y[i]+dy,z[i]+dz,B[i],occ[i],Z[i],\
      resn[i],atom[i],typ[i],chain[i])}}' |\
tee ${tempfile}new.xyz.raw |\
egrep -v "^DEBUG" >! ${tempfile}new.xyz

coordconv xyzin ${tempfile}new.xyz xyzout ${tempfile}new.pdb << EOF >> $logfile
CELL $CELL
INPUT FRAC
OUTPUT PDB
END
EOF

awk '/^ATOM|^HETAT/{print}' ${tempfile}new.pdb >> ${tempfile}out.pdb

# now inherit previous names
cat $pdbfile |\
awk '/^ATOM|^HETAT/{printf("%-100s %d %10s\n",$0,++n,"ORIG")}' |\
cat - ${tempfile}out.pdb |\
awk '$NF=="ORIG"{n=$(NF-1);\
  pref[n]=substr($0,1,30);post[n]=substr($0,55,45);next}\
 /^CRYST/ && ! crys{print;++crys;next}\
 ! /^ATOM|^HETAT/{print;next}\
 {++m;xyz=substr($0,31,24);\
  print pref[m] xyz post[m]}' |\
cat >! $outfile



if($?DEBUG || "$tempfile" == "") exit

rm -f ${tempfile}* > /dev/null

