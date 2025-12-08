#! /bin/tcsh -f
#
#        Script for extracting atoms from one PDB that are close      -James Holton 3-3-21
#         to the atoms in another PDB (symmetry-aware)
#
#

set subject = ""
set probe   = ""
set outfile = "neighbors.pdb"
set radius  = "3"
set SG      = ""

set nnearest = ""
set nfarthest = ""

set tempfile = ${CCP4_SCR}/neighbor_temp_$$_
if($#argv == 0) goto Help
goto Setup
Help:
################################################################################
cat << EOF
usage:

neighbors.com subject.pdb probe.pdb 3.0 P212121 -gather -label

where:
subject.pdb   - is the pdb you want to extract atoms from
probe.pdb     - is the pdb whose atoms mark "interesting" positions in space
3.0           - is the radius from the atoms in probe.pdb that will be
                extracted from subject.pdb (default=${radius})
P212121       - is the space group 
-include      - include precise overlaps (identical coordinates)
-nearest n    - limit output to closest n atoms within the radius
-farthest n   - limit output to most distant n atoms within the radius
-gather       - move subject atoms to be within the probe radius
-label        - label output file with nearest-neighbor distances
-outfile x    - specify output file other than $outfile 

This script extracts any atoms in subject.pdb that are within a given distance
of the atoms in probe.pdb.  Symmetry is considered.  A copy of the subject.pdb,
atoms within the specified radius of probe.pdb's atoms is output to "$outfile".

Using a negative value for the radius will extract all of the atoms in
subject.pdb EXCEPT those that are withing the specified distance of the atoms
in probe.pdb.

EOF
################################################################################
exit 9
Return_from_Setup:
if(! -e "$subject" || ! -e "$probe") goto Help
if($#CELL != 6) then
    set BAD = "need a unit cell in at least one PDB file"
    goto exit
endif

if($?debug) set tempfile = neighbors_temp
set outprefix = `basename $outfile .pdb`

set temp = ""
if("$logic" != "") set temp = " not"
if( "$nnearest" != "" ) set temp = " the $nnearest nearest, out to"
if( "$nfarthest" != "" ) set temp = " the $nfarthest farthest away, but still"

echo "extracting atoms from $subject in $SG"
echo "that are$temp within $radius A of $probe"

# extract header of the subject file
# convert to easy-to-pass labels: atom# = (resid + 1000*int(B))
cat $subject |\
awk '/^CRYST|^SCALE/{print;next}\
  /^ATOM/ || /^HETATM/{++n; \
 resid=n%10000;B=int(n/10000);\
 printf "ATOM %6d  O%-2d XXX A%4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",\
 n,B,resid,substr($0, 31, 8), substr($0, 39, 8), substr($0, 47, 8),B}' |\
cat >! ${tempfile}subj.pdb
# count number of atoms in first (subject) file
set subject_atoms  = `awk '/^ATOM/ || /^HETATM/{print} END{print "extra"}' ${tempfile}subj.pdb | wc -l`

# if more than 60000, need multiple runs
@ passes = ( ( $subject_atoms / 50000 ) + 1 )

foreach pass ( `seq 1 $passes` )

echo $pass |\
cat - ${tempfile}subj.pdb |\
awk 'NR==1{pass=$1;next}\
  /^CRYST/{print} ! /^ATOM|^HETAT/{next}\
  {++n} n>(pass-1)*50000 && n<=pass*50000{print}' |\
cat >! ${tempfile}both.pdb

set these_atoms = `egrep "^ATOM|^HETAT" ${tempfile}both.pdb | wc -l`
@ these_atoms = ( $these_atoms + 1 )

# convert probe atoms to easily-recognized atom labels
cat $probe |\
awk '/^ATOM/ || /^HETATM/{++n;\
 resid=n%10000;B=int(n/10000);\
 printf "ATOM   %4d  X%-2d YYY X%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",\
 n,B,1,substr($0, 31, 8), substr($0, 39, 8), substr($0, 47, 8)}' |\
cat >> ${tempfile}both.pdb

# renumber the atoms so that distang won't get confused
cat ${tempfile}both.pdb |\
awk '/^ATOM/{++n; $0 = sprintf("ATOM%7d%s",n,substr($0,12))} {print $0}' |\
cat >! ${tempfile}.pdb
mv ${tempfile}.pdb ${tempfile}both.pdb >& /dev/null

set half_radius = `echo "$radius" | awk '{print $1/2}'`
set dmin = 0.0001
if($?INCLUDE_PROBE) set dmin = 0
distang xyzin ${tempfile}both.pdb << EOF >! ${tempfile}distang
SYMM $SG
DIST ALL
RADII O $half_radius
RADII X $half_radius
DMIN $dmin
FROM ATOM 1 to $these_atoms
TO   ATOM $these_atoms to 99999
END
EOF

# convert back to ordinal atom numbers
cat ${tempfile}distang |\
awk '/Atom I  Atom J   Dij/,/Symmetry matrix/' |\
awk -v offset=0 '$1=="Z" && $3!~/^X/{\
    hibits=substr($3,2)*10000;\
    atom=$2%10000+hibits+offset;\
    dist=$9;symop=$10;\
    split($0,w,"[][]");\
    split(w[2],x,"[ +-]");\
    print "CLOSE", atom,dist,symop,x[1]+0,x[2]+0,x[3]+0}' |\
sort -k3n |\
cat >! ${tempfile}close_atoms_${pass}
end
sort -k3n ${tempfile}close_atoms_* >! ${tempfile}close_atoms

set close_atoms = `awk '{print $2}' ${tempfile}close_atoms | sort -u | wc -l`


if("$nnearest" != "") then
    echo "$close_atoms atoms within $radius A " 
    set near_atoms = $nnearest
    if( "$logic" != "" ) then
        # want nnearest atoms out of these
        set logic = ""
#        @ near_atoms = ( $subject_atoms - $nnearest - 1 )
#        echo "really want $near_atoms , right? "
    endif
    set distrange = `awk 'NR==1{mind=$3} END{print mind,"-",$3}' ${tempfile}close_atoms`
    echo "distance range: $distrange A away"
    awk '! seen[$2]{print;++seen[$2]}' ${tempfile}close_atoms |\
    head -n $near_atoms >! ${tempfile}nnearest
    mv ${tempfile}close_atoms ${tempfile}close_atoms1
    mv ${tempfile}nnearest ${tempfile}close_atoms
    set distrange = `awk 'NR==1{mind=$3} END{print mind,"-",$3}' ${tempfile}close_atoms`
    echo "taking nearest $nnearest"
    set close_atoms = `awk '{print $2}' ${tempfile}close_atoms | sort -u | wc -l`
endif

if("$nfarthest" != "") then
    echo "$close_atoms atoms within $radius A "
    set far_atoms = $nfarthest
    if( "$logic" != "" ) then
        # want nfarthest atoms that are still within $radius
        set logic = ""
#        @ far_atoms = ( $subject_atoms - $nfarthest - 1 )
#        echo "selecting from $far_atoms , right? "
    endif
    set distrange = `awk 'NR==1{mind=$3} END{print mind,"-",$3}' ${tempfile}close_atoms`
    echo "distance range: $distrange A away"
    awk '! seen[$2]{print;++seen[$2]}' ${tempfile}close_atoms |\
    tail -n $far_atoms >! ${tempfile}nfarthest
    mv ${tempfile}close_atoms ${tempfile}close_atoms1
    mv ${tempfile}nfarthest ${tempfile}close_atoms
    set distrange = `awk 'NR==1{mind=$3} END{print mind,"-",$3}' ${tempfile}close_atoms`
    echo "taking farthest $nfarthest"
    set close_atoms = `awk '{print $2}' ${tempfile}close_atoms | sort -u | wc -l`
endif

# now extract atoms from subject that are within the specified distance
# from to probe atoms
echo $?LABELSUFFIX |\
cat - ${tempfile}close_atoms $subject |\
awk 'NR==1{distlabel=$1;next}\
     /^CLOSE/{Close[$2]=1;if(distlabel) label[$2]="          "$3; next} \
     /^ATOM/ || /^HETATM/{++n; if('"$logic"' Close[n]) print $0 label[n]; next}\
      {print}' |\
cat >! ${outfile}

set n = `egrep "^ATOM|^HETAT" ${outfile} | wc -l`
echo "$n atoms is ready in $outfile"

if($?GATHER_ATOMS && "$logic" == "") then
    echo -n "gathering atoms together "
    set xl = 1.5
    cp -p ${outfile} ${tempfile}sparse.pdb
    echo -n "" >! ${tempfile}gathered.pdb
    foreach atom ( `awk '/^ATOM|^HETAT/{print ++n}' ${tempfile}sparse.pdb` )
        echo -n "."
        awk -v atom=$atom '! /^ATOM|^HETAT/{print;next} {++n} n==atom{print;exit}' ${tempfile}sparse.pdb |\
        cat >! ${tempfile}atom.pdb
gensymagain:
        gensym XYZIN ${tempfile}atom.pdb XYZOUT ${tempfile}test.pdb << EOF >! ${tempfile}gensym.log
BROOK
SYMM $SG
CELL $CELL
XYZLIM -$xl $xl -$xl $xl -$xl $xl
read ${tempfile}atom.pdb
EOF
        awk '{print $0,"PROBE"}' $probe |\
        cat - ${tempfile}test.pdb |\
        awk -v minD=$radius '! /^ATOM|^HETAT/{next}\
             {X=substr($0,31,8)+0;Y=substr($0,39,8)+0;Z=substr($0,47,8)+0}\
             $NF=="PROBE"{++n;X0[n]=X;Y0[n]=Y;Z0[n]=Z;next}\
             {for(i=1;i<=n;++i){\
               d=sqrt((X-X0[i])^2+(Y-Y0[i])^2+(Z-Z0[i])^2);\
               if(d<=minD) print $0,d;\
             }\
           }' |\
        cat ${tempfile}atom.pdb - |\
        awk '! /^ATOM|^HETAT/{next}\
           id==""{id=substr($0,1,30);rest=substr($0,55);next}\
           {print id substr($0,31,24) rest,$NF}' |\
        tee ${tempfile}check.pdb |\
        cat >> ${tempfile}gathered.pdb
        set test = `cat ${tempfile}check.pdb | wc -l`
        if($test == 0) then
            echo ""
            echo "ERROR: atom disappeared! "
            set xl = `echo $xl | awk '{print $1*1.5}'`
            echo "increasing range to $xl unit cells"
            cat ${tempfile}atom.pdb
            goto gensymagain
        endif
    end
    echo ""
    set n = `egrep "^ATOM|^HETAT" ${tempfile}gathered.pdb | wc -l`
    echo "$n atoms processed"
    

    # stitch probe together with other neighbors, eliminating any duplicates
    set include = $probe
    if(! $?INCLUDE_PROBE) set include = ""
    egrep "^CRYST1" ${tempfile}sparse.pdb | head -n 1 >! ${outprefix}_gathered.pdb
    cat $include ${tempfile}gathered.pdb |\
    awk '! /^ATOM|^HETAT/{next}\
       {++n;line[n]=$0;X[n]=substr($0,31,8)+0;Y[n]=substr($0,39,8)+0;Z[n]=substr($0,47,8)+0}\
       {for(i=1;i<n;++i){\
           d=sqrt((X[i]-X[n])^2+(Y[i]-Y[n])^2+(Z[i]-Z[n])^2);\
           if(d<0.001){\
              print "REMARK OVERLAP",n,$0;\
              print "REMARK OVERLAP",i,line[i];\
              next};\
         }\
        print;\
       }' |\
    cat >> ${outprefix}_gathered.pdb

    set n = `egrep "^ATOM|^HETAT" ${outprefix}_gathered.pdb | wc -l`
    echo "$n atoms is ready in ${outprefix}_gathered.pdb"

#    rm -f ${tempfile}sparse.pdb >& /dev/null
#    rm -f ${tempfile}gathered.pdb >& /dev/null
endif

exit:
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

if($?debug) exit
rm -f ${tempfile}close_atoms >& /dev/null
rm -f ${tempfile}distang >& /dev/null
rm -f ${tempfile}both.pdb >& /dev/null

exit

###################################



Setup:

set logic = ""
set CELL

# scan command line
set i = 0
while( $i < $#argv )
    @ i = ( $i + 1 )
    @ nexti = ( $i + 1 )
    @ lasti = ( $i - 1 )
    if($nexti > $#argv) set nexti = $#argv
    if($lasti < 0) set lasti = 0
    set arg = "$argv[$i]"
    set nextarg = ""
    if($nexti <= $#argv) set nextarg = $argv[$nexti]


    # recognize pdb files
    if(("$arg" =~ *.pdb)||("$arg" =~ *.brk)) then
        if(! -e "$arg") then
            echo "WARNING: $arg does not exist! "
            continue
        endif
        set test = `awk '/^CRYST/{print $2, $3, $4, $5, $6, $7}' $arg | tail -1`
        if($#test == 6) set CELL = ( $test )

        # maybe contains space group?
        if("$SG" == "") then
            set pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $arg | head -1`
            if("$pdbSG" != "") then
                if("$pdbSG" == "R 32") set pdbSG = "R 3 2"
                if("$pdbSG" == "P 21") set pdbSG = "P 1 21 1"
                if("$pdbSG" == "R 3 2" && $CELL[6] == 120.00) set pdbSG = "H 3 2"
                set SG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
                if("$SG" == R3 && $CELL[6] == 120.00) set SG = H3
            endif
        endif
        
        if(! -e "$subject") then
            # first PDB is the subject
            set subject = "$arg"

            # probably won't run
            if("$CELL" == "") echo "WARNING: no unit cell in $arg"
            continue
        endif
        
        if(! -e "$probe") then
            # second PDB is probe
            set probe = "$arg"
            
            # see if a radius is indicated in the probe file
            set temp = `awk 'BEGIN{RS=" "} $1=="RADIUS"{getline; if($1+0>0) print $1+0}' $probe`
            if("$temp" != "") set radius = "$temp"
            continue
        endif
        
        continue
    endif
    
    # recognize radius
    set temp = `echo "$arg" | awk '$1+0 > 0{print $1+0}'`
    if("$temp" != "") then
        set radius = "$temp"
        continue
    endif
    set temp = `echo "$arg" | awk '$1+0 < 0{print 0-$1}'`
    if("$temp" != "") then
        set radius = "$temp"
        # negative radius inverts selection
        set logic  = " ! "
        continue
    endif

    # recognise SG
    if ("$arg" =~ [PpCcIiFfHhRr][1-6]*) then
        set temp = `echo $arg | awk '{print toupper($1)}'`
        set temp = `awk -v SG=$temp '$4 == SG {print $4}' $CLIBD/symop.lib | head -1`
        if("$temp" != "") then
            # add this SG to the space group list
            set SG = "$temp"
            continue
        endif
    endif

    # options with sub-arguments
    if ("$arg" =~ "-nearest"*) then
        set nnearest = $nextarg
        @ i = ( $i + 1 )
        continue
    endif
    if ("$arg" =~ "-farthest"*) then
        set nfarthest = $nextarg
        @ i = ( $i + 1 )
        continue
    endif
    if ("$arg" =~ "-outfile"*) then
        set outfile = $nextarg
        @ i = ( $i + 1 )
        continue
    endif

    # recognise switches
    if ("$arg" =~ "-gath"*) then
        set GATHER_ATOMS
        continue
    endif
    if ("$arg" == "-label") then
        set LABELSUFFIX
        continue
    endif
    if ("$arg" =~ "-includ"*) then
        set INCLUDE_PROBE
        continue
    endif
    if ("$arg" =~ "-exclud"*) then
        unset INCLUDE_PROBE
        continue
    endif
    if ("$arg" == "-debug") then
        set debug
        continue
    endif
end



goto Return_from_Setup

