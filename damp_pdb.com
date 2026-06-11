#! /bin/tcsh -f
#
#   move a refined PDB back toward the starting point                        -James Holton 9-2-17
#
#   needs:
#   rmsd
#
#   usage: damp_pdb.com  start.pdb refined_001.pdb
#
#   creates: new.pdb
#
#
set firstpdb  = ""
set toofarpdb = ""
set outpdb    = new.pdb
set logfile = /dev/null
set tempfile = tempfile$$damp

# maximum allowed changes
set maxdXYZ = 0.1
set maxdocc = 0.05
set maxdB   = 5

# maximum rms changes
set maxRMSdXYZ = 99
set maxRMSdocc = 99
set maxRMSdB   = 99

# max/min absolute values
set minocc = 0
set maxocc = 10
set minB   = 5
set maxB   = 100

# apply damping to enforce these limits
set damp   = 1
set scale_XYZ = 1
set scale_occ = 1
set scale_B   = 1

# apply overall scale and B shifts
set oscale = 1
set Badd   = 0

# put debugging messages into PDB file
set debug = 0

foreach arg ( $* )
    if("$arg" =~ *.pdb) then
        if(-e "$firstpdb" && -e "$toofarpdb") then
            set outpdb = $arg
            continue
        endif
        if(! -e "$arg") then
            echo "ERROR: $arg does not exist"
            exit 9
        endif
        if("$firstpdb" == "") then
            set firstpdb = "$arg"
        else
            set toofarpdb = "$arg"
        endif
        continue
    endif
    if("$arg" == debug) then
        set debug = 1
        set logfile = /dev/tty
    endif
    if("$arg" =~ maxdXYZ=*) set maxdXYZ = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ maxdocc=*) set maxdocc = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ maxdB=*)   set maxdB   = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ maxRMSdXYZ=*) set maxRMSdXYZ = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ maxRMSdocc=*) set maxRMSdocc = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ maxRMSdB=*)   set maxRMSdB   = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ minocc=*)  set minocc = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ maxocc=*)  set maxocc = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ minB=*)    set minB   = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ maxB=*)    set maxB   = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" == damp) set damp = 1
    if("$arg" == nodamp) set damp = 0
    if("$arg" =~ scale_XYZ=*)  set USER_scale_XYZ = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ oscale=*)  set oscale = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ Badd=*)    set Badd   = `echo $arg | awk -F "=" '{print $2+0}'`
end

if(! -e "$firstpdb" || ! -e "$toofarpdb") then
    cat << EOF
usage: $0 reference.pdb toofar.pdb [damp] [maxdXYZ=0.01]

where:
reference.pdb   is the starting point of the move to be reduced
toofar.pdb      is the end point of the move to be reduced
maxdXYZ         maximum change in any XYZ position allowed
maxdocc         maximum change in any occupancy allowed
maxdB           maximum change in any B factor allowed
maxRMSdXYZ      maximum RMS change over all XYZ positions allowed
maxRMSdocc      maximum RMS change over all occupancies allowed
maxRMSdB        maximum RMS change over all B factors allowed
damp            implement move reduction with an overall scale-back
nodamp          implement move reduction as a cut-off for each atom
scale_XYZ       manual scale
oscale          specify scale factor on all occupancy values
Badd            add a constant to all B factors
minocc          minimum final absolute value of occupancy allowed
maxocc          maximum final absolute value of occupancy allowed
minB            minimum final absolute value of B factor allowed
maxB            maximum final absolute value of B factor allowed
EOF
    exit 1
endif


if($damp) then
    echo "gathering stats"
    set movestats = `rmsd -v xlog=1 $firstpdb $toofarpdb | head -n 1`

    set RMSmove_XYZ = `echo $movestats | awk '{print $2}'`
    set RMSmove_occ = `echo $movestats | awk '{print $3}'`
    set RMSmove_B   = `echo $movestats | awk '{print $4}'`

    set maxmove_XYZ = `echo $movestats | awk '{print $5}'`
    set maxmove_occ = `echo $movestats | awk '{print sqrt($6^2)}'`
    set maxmove_B   = `echo $movestats | awk '{print sqrt($7^2)}'`

    # scale-back due to MAX limit
    set scale_XYZ_max = `echo $maxmove_XYZ $maxdXYZ | awk '$1+0<$2{$1=$2} $1==0{print 1;exit} {print $2/$1}'`
    set scale_occ_max = `echo $maxmove_occ $maxdocc | awk '$1+0<$2{$1=$2} $1==0{print 1;exit} {print $2/$1}'`
    set scale_B_max   = `echo $maxmove_B   $maxdB   | awk '$1+0<$2{$1=$2} $1==0{print 1;exit} {print $2/$1}'`

    # scale-back due to RMS limits
    set scale_XYZ_rms = `echo $RMSmove_XYZ $maxRMSdXYZ | awk '$1+0<$2{$1=$2} $1==0{print 1;exit} {print $2/$1}'`
    set scale_occ_rms = `echo $RMSmove_occ $maxRMSdocc | awk '$1+0<$2{$1=$2} $1==0{print 1;exit} {print $2/$1}'`
    set scale_B_rms   = `echo $RMSmove_B   $maxRMSdB   | awk '$1+0<$2{$1=$2} $1==0{print 1;exit} {print $2/$1}'`

    # take most conservative of the two
    set scale_XYZ = `echo $scale_XYZ_max $scale_XYZ_rms | awk '{print ($1>$2?$2:$1)}'`
    set scale_occ = `echo $scale_occ_max $scale_occ_rms | awk '{print ($1>$2?$2:$1)}'`
    set scale_B   = `echo $scale_B_max $scale_B_rms   | awk '{print ($1>$2?$2:$1)}'`
endif

if($?USER_scale_XYZ) then
    set scale_XYZ = $USER_scale_XYZ
endif

cat << EOF
oscale = $oscale
Badd   = $Badd

scale_XYZ = $scale_XYZ
scale_occ = $scale_occ
scale_B   = $scale_B

maxdXYZ = $maxdXYZ
maxdocc = $maxdocc
maxdB   = $maxdB

maxRMSdXYZ = $maxRMSdXYZ
maxRMSdocc = $maxRMSdocc
maxRMSdB   = $maxRMSdB

minocc = $minocc
maxocc = $maxocc
minB   = $minB
maxB   = $maxB
EOF


# get headers from the first one
awk '/^ATOM|^HETAT/{exit} {print}' $firstpdb >! ${tempfile}outpdb

# apply filters
echo "$maxdXYZ $maxdocc $maxdB  $minocc $maxocc  $minB $maxB  $scale_XYZ $scale_occ $scale_B   $oscale $Badd" |\
cat - $firstpdb $toofarpdb |\
awk 'NR==1{maxdX=$1;maxdocc=$2;maxdB=$3;minocc=$4;maxocc=$5;minB=$6;maxB=$7;\
           scale_XYZ=$8;scale_occ=$9;scale_B=$10;oscale=$11;Badd=$12;next}\
   ! /^ATOM|^HETAT/{next}\
   {id=substr($0,12,15);++seen[id];X=substr($0,31,8)+0;Y=substr($0,39,8)+0;Z=substr($0,47,8)+0;\
    occ=substr($0,55,6)+0;B=substr($0,61,6)+0;\
    rest=substr($0,67)}\
   seen[id]>2{print "WARNING: more than two copies of",id}\
   seen[id]==1{++n;X0[n]=X;Y0[n]=Y;Z0[n]=Z;occ0[n]=occ;B0[n]=B;onum[id]=n;next}\
   seen[id]==2{m=onum[id];occ*=oscale;B+=Badd;\
      dx=X-X0[m];dy=Y-Y0[m];dz=Z-Z0[m];docc=occ-occ0[m];dB=B-B0[m];\
      dx*=scale_XYZ;dy*=scale_XYZ;dz*=scale_XYZ;docc*=scale_occ;dB*=scale_B;\
      absdxyz=sqrt(dx^2+dy^2+dz^2);absdO=sqrt(dO^2);absdB=sqrt(dB^2);\
      if(absdxyz>maxdX){scale=(maxdX/absdxyz);dx*=scale;dy*=scale;dz*=scale}\
      if(absdocc>maxdO){docc=maxdO*docc/absdocc}\
      if(absdB>maxdB){dB=maxdB*dB/absdB}\
      X=X0[m]+dx;\
      Y=Y0[m]+dy;\
      Z=Z0[m]+dz;\
      occ=occ0[m]+docc;\
      B=B0[m]+dB;\
      if(occ<minocc){occ=minocc};if(occ>maxocc){occ=maxocc};\
      if(B<minB){B=minB};if(B>maxB){B=maxB};\
      printf("ATOM  %5d%s    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",m%100000,id,X,Y,Z,occ,B,rest);\
   }' |\
cat >> ${tempfile}outpdb
mv ${tempfile}outpdb $outpdb

egrep "^ATOM|^HETAT" $outpdb | wc -l | awk '{print $1,"atoms processed."}'



