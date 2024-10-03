#! /bin/tcsh -f
#
#   Script for expanding a PDB/MTZ file to a super cell, and back
#

set mono_map = monomer_rot_trans.txt
set asu_centers = asu_centers.pdb
set pdir = `dirname $0`
if(! -e "$mono_map") set mono_map = monomer_rot_trans.txt
if(! -e "$mono_map") set mono_map = ../alignment/$mono_map
if(! -e "$mono_map") set mono_map = ../$mono_map
if(! -e "$mono_map") set mono_map = ../$mono_map
if(! -e "$mono_map") set mono_map = ../$mono_map
if(! -e "$mono_map") set mono_map = ${pdir}/alignment/$mono_map
if(! -e "$mono_map") then
    set BAD = "cannot find monomer_rot_trans.txt"
    goto exit
endif
if(! -e "$asu_centers") set asu_centers = `dirname $mono_map`/asu_centers.pdb
set mono_step = `awk 'NR==2{print $NF-1;exit}' $mono_map`
set super_mult  = ( 2 2 3 )
set smallSG = P212121
set smallSGnum = `awk -v sg=$smallSG '$4==sg{print $1;exit}' ${CLIBD}/symop.lib`
set nsymops = `awk -v n=$smallSGnum '$1==n{print $3;exit}' ${CLIBD}/symop.lib`

set tempfile = /dev/shm/${USER}/supercellme$$
mkdir -p /dev/shm/${USER}/
set logfile = /dev/null

set ABC = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_.,;=<>+-%&*$!(){}[]^#|~?'
set confseq = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuv"

set pdbfile = refmacout.pdb
set outfile = default
set mtzfile = ""
set outputmtz = ""
set keepconf = 0
set keepchain = 0
set keepresnum = 0
set keepocc = 1
set newocc = "no"
set scale4map = 1
set blabel = 0
set debug = 0

echo "command-line arguments: $* "

foreach arg ( $* )
    set key = `echo $arg | awk -F "=" '{print $1}'`
    set val = `echo $arg | awk -F "=" '{print $2}'`
    set csv = `echo $val | awk -F "," '{gsub(","," ");print}'`

    if("$arg" =~ *.pdb && "$val" == "") set pdbfile = "$arg"
    if("$arg" =~ *.mtz && "$val" == "") set mtzfile = "$arg"
    if("$key" == "pdbfile") set pdbfile = "$val"
    if("$key" == "mtzfile") set mtzfile = "$val"
    if("$key" == "outfile" || "$key" == "output" || "$key" == "outpdb") set outfile = "$val"
    if("$key" == "outmtz" || "$key" == "mtzout") set outputmtz = "$val"
    if("$key" == "confseq") set confseq = "$val"
    if("$key" == "mono_step") set mono_step = "$val"
    if("$key" == "super_mult") set super_mult = ( $csv )
    if("$key" == "smallSG") set smallSG = "$val"
    if("$key" == "keepconf") set keepconf = "$val"
    if("$key" == "keepchain") set keepchain = "$val"
    if("$key" == "keepresnum") set keepresnum = "$val"
    if("$key" == "keepocc") set keepocc = "$val"
    if("$key" == "occ") set newocc = "$val"
    if("$key" == "scale4map") set scale4map = "$val"
    if("$key" == "blabel") set blabel = "$val"
    if("$key" == "tempfile") set tempfile = "$val"
    if("$key" == "debug") set debug = "$val"
end

if(! -e "$pdbfile") then
    set BAD = "pdbfile $pdbfile does not exist."
    goto exit
endif

echo "confseq=$confseq"

# detect if this is supercell or not
set unitcell = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7;exit}' $pdbfile`
set chains = `awk '/^ATOM/{c=substr($0,22,1)} ! seen[c]{print c;++seen[c]}' $pdbfile`
set cellsize = `awk '/^CRYST/{print int(($2*$3*$4)**(1./3));exit}' $pdbfile`

set supercell = 1
if( $cellsize < 50 || "$pdbfile" =~ *_small.pdb ) set supercell = 0


if( $supercell ) goto collapse

set ncells = `echo $super_mult | awk '{printf("%.0f\n",$1*$2*$3)}'`

echo "expanding to supercell"
if("$outfile" == "default") set outfile = supercell.pdb

# calculate the supercell
echo $unitcell $super_mult |\
 awk '{a=$1*$7;b=$2*$8;c=$3*$9;al=$4;be=$5;ga=$6;\
   printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n",a,b,c,al,be,ga);}' |\
cat >! $tempfile

echo $confseq $mono_step $keepconf $keepchain $keepresnum $keepocc $newocc $blabel |\
cat - $mono_map $pdbfile |\
awk 'NR==1{confseq=$1;modulo=$2;\
        keepconf=$3;keepchain=$4;keepresnum=$5;keepocc=$6;newocc=$7;blabel=$8;\
     next}\
  NR<=49{a=$1;c[a]=substr(confseq,a,1);asu[c[a]]=a;\
    xx[a]=$2;xy[a]=$3;xz[a]=$4;\
    yx[a]=$5;yy[a]=$6;yz[a]=$7;\
    zx[a]=$8;zy[a]=$9;zz[a]=$10;\
    tx[a]=$11;ty[a]=$12;tz[a]=$13;\
    dresnum[a]=$14-1;next}\
  /^TER/{print "TER"}\
  ! /^ATOM|^HETAT/{next}\
  /^ATOM|^HETAT/{pre=substr($0,1,16);post=substr($0,67);typ=substr($0,18,4);\
          conf=substr($0,17,1);chain=substr($0,22,1);\
          resnum=substr($0,23,4)+0;\
          occ=substr($0,55,6)+0;B=substr($0,61,6);\
          solv=( /^HETAT/ || typ~/HOH|ACY|NH4/ );\
          ins=substr($0,27,4);\
          X=substr($0,31,8)+0;\
          Y=substr($0,39,8)+0;\
          Z=substr($0,47,8)+0;\
          a=asu[chain];\
          if(a+0==0){a=1;post=post"    UNK"}\
          newresnum=resnum+dresnum[a];\
          if(! keepchain) {chain="A";if(/HOH/) chain="S"};\
          if(! keepconf ) conf=" ";\
          if(! keepocc ) occ=1;\
          if( newocc != "no" ) occ=newocc;\
          if(keepresnum || solv ) newresnum=resnum;\
          if(blabel ) B=a;\
          newX=xx[a]*X+xy[a]*Y+xz[a]*Z+tx[a];\
          newY=yx[a]*X+yy[a]*Y+yz[a]*Z+ty[a];\
          newZ=zx[a]*X+zy[a]*Y+zz[a]*Z+tz[a];\
         printf("%s%s%s%s%4d%s%8.3f%8.3f%8.3f%6.2f%6.2f%s\n",pre,conf,typ,chain,newresnum,\
            ins,newX,newY,newZ,occ,B,post)}' |\
cat >> $tempfile

mv $tempfile $outfile

if(! -e "$mtzfile") goto exit

if("$outputmtz" == "") then
    set outputmtz = `dirname $mtzfile`/`basename $mtzfile .mtz`_super.mtz
endif

echo "expanding $mtzfile to supercell in $outputmtz"
cad hklin1 $mtzfile hklout ${tempfile}.mtz << EOF >> $logfile
labin file 1 all
outlim space 1
EOF
# this operation effectively multiplies by number of symops, unless you run file through cad agains
cad hklin1 ${tempfile}.mtz hklout ${tempfile}P1.mtz << EOF  >> $logfile
labin file 1 all
scale file 1 $ncells 0
symm 1
EOF
echo $super_mult |\
awk '{print "reindex "$1"h,"$2"k,"$3"l"}' |\
 reindex hklin ${tempfile}P1.mtz hklout $outputmtz >> $logfile


if( "${tempfile}" != "" && "${tempfile}" != "./" ) then
  rm -f ${tempfile}* >& /dev/null
endif

goto exit

collapse:

echo "collapsing $super_mult supercell to one asu"

if("$outfile" == "default") set outfile = aligned.pdb

echo $unitcell $super_mult |\
 awk '{a=$1/$7;b=$2/$8;c=$3/$9;al=$4;be=$5;ga=$6;\
   printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 21 21 21\n",a,b,c,al,be,ga);}' |\
cat >! ${tempfile}.pdb

echo space $smallSG | pdbset xyzin ${tempfile}.pdb xyzout ${tempfile}sg.pdb >> $logfile
egrep "^CRYST" ${tempfile}sg.pdb >! ${tempfile}
rm -f ${tempfile}.pdb ${tempfile}sg.pdb >> $logfile

echo $confseq $mono_step $keepconf $keepchain $keepresnum $keepocc $newocc $blabel $debug |\
cat - $mono_map $asu_centers $pdbfile |\
awk 'NR==1{confseq=$1;modulo=$2;\
        keepconf=$3;keepchain=$4;keepresnum=$5;keepocc=$6;newocc=$7;blabel=$8;\
        debug=$NF;\
     next}\
  NR<=49{a=$1;c[a]=substr(confseq,a,1);asu[c[a]]=a;\
    xx[a]=$2;xy[a]=$3;xz[a]=$4;\
    yx[a]=$5;yy[a]=$6;yz[a]=$7;\
    zx[a]=$8;zy[a]=$9;zz[a]=$10;\
    tx[a]=$11;ty[a]=$12;tz[a]=$13;\
    dresnum[a]=$14-1;asu[dresnum[a]]=a;next}\
  /^TER/{print "TER"}\
  ! /^ATOM|^HETAT/{next}\
  /^ATOM|^HETAT/{pre=substr($0,1,16);post=substr($0,61);typ=substr($0,18,4);\
          conf=substr($0,17,1);chain=substr($0,22,1);\
          resnum=substr($0,23,4)+0;\
          occ=substr($0,55,6)+0;B=substr($0,61,6);\
          solv=( /^HETAT/ || typ~/HOH|ACY|NH4/ );\
          ins=substr($0,27,4);\
          X=substr($0,31,8)+0;\
          Y=substr($0,39,8)+0;\
          Z=substr($0,47,8)+0;\
          if($NF=="ASUCENTER"){++asuc;\
              Xc[asuc]=X;Yc[asuc]=Y;Zc[asuc]=Z;next}\
          newresnum=((resnum-1)%modulo)+1;\
          a=asu[resnum-newresnum];\
          if(a+0<=0 || solv){newresnum=resnum;\
            min=-1;for(i=1;i<=asuc;++i){\
            dsq=(X-Xc[i])^2+(Y-Yc[i])^2+(Z-Zc[i])^2;\
            if(dsq<min || min<0){a=i;min=dsq};\
          };}\
          if(solv){uid=typ conf chain resnum;\
            if(aof[uid]=="")aof[uid]=a;\
               a=aof[uid]};\
          if(a+0<=0){a=1;post=post"     UNK"}\
          newchain=chain;newconf=conf;\
          if(! keepchain) newchain=c[a];\
          if(! keepconf ) newconf=c[a];\
          if(! keepocc ) occ=0.02;\
          if( newocc != "no" ) occ=newocc;\
          if(keepresnum ) newresnum=resnum;\
          if(blabel)a=int(B);\
          X=X-tx[a];\
          Y=Y-ty[a];\
          Z=Z-tz[a];\
          newX=xx[a]*X+yx[a]*Y+zx[a]*Z;\
          newY=xy[a]*X+yy[a]*Y+zy[a]*Z;\
          newZ=xz[a]*X+yz[a]*Y+zz[a]*Z;\
         printf("%s%s%s%s%4d%s%8.3f%8.3f%8.3f%6.2f%s\n",\
            pre,newconf,typ,newchain,newresnum,ins,newX,newY,newZ,occ,post)}' |\
cat >> $tempfile

mv $tempfile $outfile

if(! -e "$mtzfile") goto exit

if("$outputmtz" == "") then
    set outputmtz = `dirname $mtzfile`/`basename $mtzfile .mtz`_small.mtz
endif

# re-index supercell, preserving overcompleteness
echo "collapsing $mtzfile to ASU"
set reidx = `echo $super_mult | awk '{print "reindex h/"$1",k/"$2",l/"$3}'`
set cellscale = `echo $super_mult | awk '{printf("%.20f\n",1./($1*$2*$3))}'`
set scale = `echo $cellscale $nsymops $scale4map | awk '$NF==0{print $1;exit} {printf(".20f",$1/$2)}'`
cad hklin1 $mtzfile hklout ${tempfile}.mtz << EOF >> reindex.log
labin file 1 all
scale file 1 $scale 0
EOF
reindex hklin ${tempfile}.mtz hklout $outputmtz << EOF >! reindex.log
$reidx
symm $smallSG
EOF

# note that no scaling by sympos is done here because the supercell data are in P1
# resulting small.mtz is "overcomplete" with different data in each ASU. Coot and FFT
# however, will make a normal map out of it.
# it would be good practice to divide by number of symops duing the FFT run



exit:

if( "${tempfile}" != "" && "${tempfile}" != "./" ) then
  rm -f ${tempfile}* >& /dev/null
endif

if($?BAD) then
    echo "ERROR $BAD"
    exit 9
endif


ls -l $outfile

exit






grep CRYST1 refmacout_small.pdb >! center.pdb
foreach i ( `seq 1 48` )

echo "$i $confseq" |\
awk 'NR==1{i=$1;confseq=$2;c=substr(confseq,i,1);\
  printf("ATOM %6d  O  %sHOH %s%4d       3.088   3.729   0.773  1.00  5.00           O     ASUCENTER\n",i,c,c,1);}' |\
 tee -a center.pdb

end

superxform.com center.pdb 

awk '{gsub("ATOM  ","HETATM");gsub("HOH A","HOH S");print}' supercell.pdb >! ../../alignment/asu_centers.pdb

set CELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`
egrep "^CRYST|^ATOM" refmacout.pdb >! condensed.pdb
egrep "^CRYST|^HETAT" refmacout.pdb >! hetat.pdb
coordconv xyzin hetat.pdb xyzout ${tempfile}.xyz << EOF 
CELL $CELL
INPUT PDB
OUTPUT FRAC
END
EOF

cat ${tempfile}.xyz |\
awk '{n=$1;x=$2;y=$3;z=$4;B=$5;occ=$6;Z=$7;\
     rn=substr($0,56,10);atom=substr($0,66,5);\
     TYP=substr($0,71,3);chain=$NF;\
     x-=int(x);if(x<0)x+=1;\
     y-=int(y);if(y<0)y+=1;\
     z-=int(z);if(z<0)z+=1;\
printf "%5d%10.5f%10.5f%10.5f%10.5f%5.2f%5d%10d%-5s%3s %1s\n", \
       n,x,y,z,B,occ,Z,rn,atom,TYP,chain}' |\
cat >! new.xyz

coordconv xyzin new.xyz xyzout new.pdb << EOF 
CELL $CELL
INPUT FRAC
OUTPUT PDB
END
EOF

awk '/^ATOM/{print "HETATM" substr($0,7)}' new.pdb >> condensed.pdb

superxform.pdb condensed.pdb


