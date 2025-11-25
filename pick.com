#! /bin/tcsh -f
#
#        pick.com                - James Holton        1-3-21
#
#        Pick unique peaks in a map, 
#        avoiding "map-edge" false peaks
#        optionally avoiding a list of "boring" positions
#        in a "symmetry-aware" fashion
#
# defaults
set mapfile = "maps/FH_Four.map"
set pdbfile = ""
set logfile = "pick_details.log"
set outfile = "pick.pdb"

set sigma   = 6
set top_peaks = ""
set bottom_peaks = ""

set tempfile = ${CCP4_SCR}/pick_temp$$
set tempfile = pick_temp
rm -f ${tempfile}* >& /dev/null

# set this to wherever your awk program is
alias nawk /usr/bin/nawk
nawk 'BEGIN{print}' >& /dev/null
if($status) alias nawk awk

if("$1" == "") goto Help
echo -n "" >! $logfile
################################################################################
goto Setup
# set/reset $mapfile $pdbfile $sigma from command-line
Help:
cat << EOF

usage: $0 $mapfile [$sigma] [boring.pdb]

where: $mapfile is the map you want to pick
       $sigma (optional) is the minimum peak height (sigma units)
       boring.pdb (optional) sites to avoid in peak-picking

EOF
exit 9
Return_from_Setup:
################################################################################

if(! -e "$mapfile") goto Help

set sign = `echo "$sigma" | nawk '$1+0<0{print "+/-" 0-$1} $1+0>0{print $1}'`
set print_sigma = "${sign}*sigma peaks"
if($?top_only) set print_sigma = "highest peak"
if($?bottom_only) set print_sigma = "lowest peak"
echo -n "looking for $print_sigma in $mapfile "
if(-e "$sitefile") then
    echo "not already within ${CLOSE_peaks}A"
    echo -n "of the $boring_sites atoms listed in $sitefile"
endif
echo ""

# extract a single ASU from the input map
mapmask mapin $mapfile mapout ${tempfile}asu.map << EOF-xtend | tee ${tempfile}xtend >> $logfile
#scale sigma
xyzlim ASU
# re-axis to X,Y,Z
AXIS X Y Z
# fill blank spaces with zero
pad 0
EOF-xtend

# make sigma=1 and mean=0
echo "scale sigma" | mapmask mapin ${tempfile}asu.map mapout ${tempfile}sigma.map >> $logfile
mv ${tempfile}sigma.map ${tempfile}asu.map

# get size of the ASU
cat ${tempfile}xtend |\
nawk '/Grid sampling on x, y, z/{print $(NF-2), $(NF-1), $NF;} \
      /Start and stop points on x, y, z/{print $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1), $NF}' |\
nawk 'NF==3{gx=$1;gy=$2;gz=$3}\
      NF==6{print "ASU", $1/gx, $2/gx, $3/gy, $4/gy, $5/gz, $6/gz}' |\
cat >! ${tempfile}asu
if(! $?DEBUG) rm -f ${tempfile}xtend >& /dev/null

# calculate a 10% "edge pad"
set xyzlim = `nawk '/^ASU/{print $2-0.1, $3+0.1, $4-0.1, $5+0.1, $6-0.1, $7+0.1}' ${tempfile}asu`
set asu = `nawk '/^ASU/{print $2, $3, $4, $5, $6, $7}' ${tempfile}asu`

# re-extend the map by a 10% pad in every direction
echo "xyzlim $xyzlim" |\
mapmask mapin ${tempfile}asu.map mapout ${tempfile}pick.map | tee ${tempfile}stats >> $logfile
if(! $?DEBUG) rm -f ${tempfile}asu.map >& /dev/null

set max = `awk '/Maximum density/{print $NF}' ${tempfile}stats`
set min = `awk '/Minimum density/{print -$NF}' ${tempfile}stats`
if(! $?DEBUG) rm -f ${tempfile}stats

repeat:
set NUMPEAKS 
if("$top_peaks" != "") then
    set NUMPEAKS = "NUMPEAKS $top_peaks"
    echo "looking for $top_peaks peaks only"
endif
if("$bottom_peaks" != "") then
    set NUMPEAKS = "NUMPEAKS $bottom_peaks"
    echo "looking for bottom $bottom_peaks peaks only"
    set sigma = `echo $sigma | awk '{print -$1}'`
endif
if($?top_only) set sigma = `echo $max | awk '{print $1-0.01}'`
if($?bottom_only) set sigma = `echo $min $max | awk '$1>$2{$1=$2} {print -$1+0.01}'`

# reformat to peakmax vernacular
set sigma = `echo $sigma | awk '$1+0>0{print $1} $1+0<0{print -$1,"NEGATIVES"}'`


# do the actual peak-pick
peakmax MAPIN ${tempfile}pick.map PEAKS ${tempfile}.xyz << eof-pick >> $logfile
THRESHOLD $sigma
OUTPUT PEAKS
#$NUMPEAKS
END
eof-pick
if($status) then
    grep "Threshold too high" $logfile >& /dev/null
    if(! $status) then
        set sigma = `awk '/is greater than maximum density/{print $NF*0.99}' $logfile`
        if("$sigma" == "") then
            set sigma = `echo "$sigma" | nawk '$1+0>0.5{print  $1*2/3}'`
        endif
        if("$sigma" == "") then
            echo "no peaks."
            set BAD
            goto cleanup
        endif
        if("$origsigma" =~ -*) set sigma = `echo $sigma | awk '{print -sqrt($1*$1)}'`
        echo "reducing sigma to $sigma"
        goto repeat
    endif
endif
if(! $?DEBUG) rm -f ${tempfile}pick.map >& /dev/null

# re-format the peaks list (with no stuck-together numbers)
cat ${tempfile}.xyz |\
nawk 'NF>6 && /[^1-9.-]/{ \
print substr($0,23,8), substr($0,31,8), substr($0,39,8),\
      substr($0,49,8), substr($0,57,8), substr($0,65,8), substr($0,6,8)}' |\
cat >! ${tempfile}peaks.pick
if(! $?DEBUG) rm -f ${tempfile}.xyz >& /dev/null

# filter with sigma cutoff (in case it was changed above)
echo $origsigma |\
cat - ${tempfile}peaks.pick |\
awk 'NR==1{sigma0=sqrt($1*$1);next}\
  {sigma=sqrt($NF*$NF)}\
  sigma>sigma0{print}' |\
cat >! ${tempfile}filtered.pick
set test = `cat ${tempfile}filtered.pick | wc -l`
if( $test == 0 ) then
    cp ${tempfile}peaks.pick ${tempfile}filtered.pick
endif

# now filter out out-of-bounds (map-edge) peaks
################################################################
# trim off peaks outside the CCP4 ASU limits
# (they should either be redundant, or map-edge peaks)
cat ${tempfile}asu ${tempfile}filtered.pick |\
nawk '/^ASU/{xmin=$2-0;ymin=$4-0;zmin=$6-0;\
             xmax=$3+0;ymax=$5+0;zmax=$7+0;next;} \
      {x=$1+0;y=$2+0;z=$3+0}\
      x>=xmin && x<=xmax && y>=ymin && y<=ymax && z>=zmin && z<=zmax {print}' |\
sort -nr --key=7 >! ${tempfile}peaks.trimmed

# add another 5% pad (just in case we lost a few)
cat ${tempfile}asu ${tempfile}filtered.pick |\
nawk -v pad=0.05 '/^ASU/{xmino=$2-pad;ymino=$4-pad;zmino=$6-pad;\
                         xmaxo=$3+pad;ymaxo=$5+pad;zmaxo=$7+pad;\
                         xmini=$2-0;ymini=$4-0;zmini=$6-0;\
                         xmaxi=$3+0;ymaxi=$5+0;zmaxi=$7+0;\
                         next;} \
      {x=$1+0;y=$2+0;z=$3+0}\
      x>=xmini && x<=xmaxi && y>=ymini && y<=ymaxi && z>=zmini && z<=zmaxi {next}\
      x>xmino && x<xmaxo && y>ymino && y<ymaxo && z>zmino && z<zmaxo {print}' |\
sort -nr --key=7 >> ${tempfile}peaks.trimmed
if(! $?DEBUG) rm -f ${tempfile}filtered.pick ${tempfile}asu >& /dev/null

# these peaks are sorted by "priority" of their ASU convention
# ${tempfile}peaks.trimmed had
# format: xf yf zy X Y Z height

if("$top_peaks" != "") then
    # make sure there are not more than the top specified number
    sort -k7gr ${tempfile}peaks.trimmed |\
    awk -v n=$top_peaks '! seen[$NF]{++m;++seen[$NF]}\
       m>n{exit} {print}' >! ${tempfile}
    mv ${tempfile} ${tempfile}peaks.trimmed
endif
if("$bottom_peaks" != "") then
    sort -k7g ${tempfile}peaks.trimmed |\
    awk -v n=$bottom_peaks '! seen[$NF]{++m;++seen[$NF]}\
       m>n{exit} {print}' >! ${tempfile}
    mv ${tempfile} ${tempfile}peaks.trimmed
endif

################################################################
# near-edge peaks have probably been counted twice, so we need
# to filter them out

# generate ALL symmetry-equivalent positions for the trimmed peaks
cat << EOF >! ${tempfile}gensym.in
SYMM $SGnum
CELL $CELL
XYZLIM $xyzlim
EOF
cat ${tempfile}peaks.trimmed |\
nawk '{++n; print "RESIDUE",n; print "ATOM X", $1, $2, $3}' |\
cat >> ${tempfile}gensym.in
cat ${tempfile}gensym.in | gensym |&\
nawk '/List of sites/,/Normal termination/' |\
nawk '$2 ~ /[01].[0-9][0-9][0-9]/{print $2, $3, $4, $5, $6, $7, $(NF-1), "sym"}' |\
cat >! ${tempfile}peaks.symm
if(! $?DEBUG) rm -f ${tempfile}gensym.in >& /dev/null

#  ${tempfile}peaks.symm is now an indexed list of all 
# symmetry-related peak positions within $xyzlim
#format: xf yf zf X Y Z peak# "sym" 

# to preserve the ASU coordinates:
# sort the trimmed coordinates into the list so 
# they will be the "first" symmetry mate considered 
# for each "site"
cat ${tempfile}peaks.trimmed |\
nawk '{++n; print $1, $2, $3, $4, $5, $6, n, $NF*$NF, $NF}' |\
cat - ${tempfile}peaks.symm |\
sort --key=7n --key=8nr >! ${tempfile}peaks.xpanded
# sort +6n -7 +7nr -8 
if(! $?DEBUG) rm -f ${tempfile}peaks.trimmed >& /dev/null
if(! $?DEBUG) rm -f ${tempfile}peaks.symm >& /dev/null


# now filter the symmetry-expanded list for the
# unique list of non-symmetry related peaks
cat ${tempfile}peaks.xpanded |\
nawk '! seen[$1 " " $2 " " $3] {print} {seen[$1 " " $2 " " $3]=1}' |\
nawk -v cut=$CLOSE_map '$NF!="sym"{height[$7]=$NF}\
        NF>3{++n; X[n]=$4; Y[n]=$5; Z[n]=$6; site[n]=$7;\
        # compare this peak to all sites seen so far \
        for(i=1;i<n;++i){\
            dist=sqrt(($4-X[i])^2 +($5-Y[i])^2 +($6-Z[i])^2);\
            # see if an equivalent peak has already been printed \
            if(dist < cut){ ++taken[site[n]]; break}}; \
        if(! taken[site[n]]) print $1, $2, $3, $4, $5, $6, height[site[n]];\
        # register all the symm mates as taken too \
        ++taken[site[n]]}' |\
cat >! ${tempfile}peaks.reduced
if(! $?DEBUG) rm -f ${tempfile}peaks.xpanded >& /dev/null

# ${tempfile}peaks.reduced should now contain only unique peaks from $mapin
# format: xf yf zf X Y Z height


################################################################
# now look for special positions
cat ${tempfile}peaks.reduced |\
nawk 'NF>0{++n; print "RESIDUE",n; print "ATOM X", $1, $2, $3}' |\
cat >! ${tempfile}gensym.in
gensym << EOF >! ${tempfile}.log
SYMM $SGnum
CELL $CELL
XYZLIM 0 0.999999 0 0.999999 0 0.999999
@${tempfile}gensym.in
EOF
# count the number of times each site is "seen"
# more than once implies a special position
cat ${tempfile}.log |\
nawk '/List of sites/,/Normal termination/' |\
nawk '$2 ~ /[01].[0-9][0-9][0-9]/{print $2, $3, $4, $(NF-1), $NF}' |\
nawk '{++seen[$1 " " $2 " " $3 " " $4]}\
   END{for(site in seen) print site, seen[site]}' |\
sort -un --key=4 |\
nawk '$5+0>0{print $4, 1/$5}' >! ${tempfile}occs
# sort -un +3
if(! $?DEBUG) rm -f  ${tempfile}.log >& /dev/null


# now add these "occupancies" to the master list
cat ${tempfile}occs ${tempfile}peaks.reduced |\
nawk 'NF==2{occ[$1]=$2} \
       NF>2{++n;print $1,$2,$3,$4,$5,$6,occ[n],$7}' |\
cat >! ${tempfile}peaks.final
if(! $?DEBUG) rm -f ${tempfile}peaks.reduced >& /dev/null
# ${tempfile}peaks.final now contains the "final" list of output peaks
# and should faithfully represent the unique peaks in the map
# format: xf yf zf X Y Z 1/mult height

set peaks = `cat ${tempfile}peaks.final | wc -l`
echo -n "$peaks found "
if(! $?DEBUG) rm -f  ${tempfile}occs >& /dev/null

################################################################
# filter out "boring" sites
if(! -e ${tempfile}boring_sites) touch ${tempfile}boring_sites
# format: xf yf zf

# count number of "boring" sites
set boring_sites = `cat ${tempfile}boring_sites | wc -l`

# symmetry-expand the "boring" sites
cat << EOF >! ${tempfile}gensym_head.in
SYMM $SGnum
CELL $CELL
XYZLIM $xyzlim
EOF
# gensym can only handle 5000 atoms at a time!!! ?
echo -n "" >! ${tempfile}boring_sites.symm
set chunksize = 501
set chunks = `cat ${tempfile}boring_sites | wc -l | awk -v c=$chunksize '{print int($1/c)+1}'`
foreach chunk ( `seq 1 $chunks` )
  echo -n "."
  set head = `echo $chunk $chunksize | awk '{print ($1-1)*$2}'`
  set tail = `echo $chunk $chunksize | awk '{print $1*$2-1}'`
  echo $head $tail |\
    cat - ${tempfile}boring_sites |\
  awk 'NR==1{head=$1;tail=$2;next} \
     NR>=head{print} NR>=tail{exit}' |\
  nawk 'NF>2{print "ATOM X", $1, $2, $3}' |\
  cat ${tempfile}gensym_head.in - >! ${tempfile}gensym.in

  cat ${tempfile}gensym.in | gensym |&\
  nawk '/List of sites/,/Normal termination/' |\
  nawk '$2 ~ /[01].[0-9][0-9][0-9]/{print $2, $3, $4, $5, $6, $7, "boring"}' |\
  cat >> ${tempfile}boring_sites.symm
  #rm -f ${tempfile}gensym.in >& /dev/null
end

# ${tempfile}all_boring_sites now contains ALL symmetry-equivalent
# positions to the sites the user entered on the command-line
# format: xf yf zf X Y Z "boring"

# now remove peaks that were too close to "boring" sites
cat ${tempfile}boring_sites.symm ${tempfile}peaks.final |\
nawk -v cut=$CLOSE_peaks ' \
      $NF=="boring"{++n; X[n]=$4; Y[n]=$5; Z[n]=$6} \
      $NF!="boring"{minD=999999; \
         # find nearest "boring" site \
         for(i=1;i<=n;++i){\
             dist=sqrt(($4-X[i])^2 +($5-Y[i])^2 +($6-Z[i])^2);\
             if(dist < minD){\
             minD=dist;}}\
         # now see if it is too close \
         if(minD > cut) {print}}' |\
cat >! ${tempfile}peaks.interesting
if(! $?DEBUG) rm -f  ${tempfile}boring_sites.symm >& /dev/null
if(! $?DEBUG) rm -f  ${tempfile}peaks.final >& /dev/null

# sort the picked peaks by height
sort -nr -k8 ${tempfile}peaks.interesting >! ${tempfile}
mv ${tempfile} ${tempfile}peaks.interesting

# ${tempfile}peaks.interesting now contains only "interesting" peaks
# format: xf yf zf X Y Z 1/mult height

# see how many are left
set interesting = `cat ${tempfile}peaks.interesting | wc -l`
if($boring_sites) echo -n " ($interesting) new"
echo ""

################################################################
# make a pdb output file

cat ${tempfile}peaks.interesting |\
awk '{++i; printf "ATOM   %4d  OW  WAT X%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",\
         i, i, $4, $5, $6, $7, $8}' |\
cat >> ${tempfile}atoms.pdb
echo "END" >> ${tempfile}atoms.pdb
pdbset xyzin ${tempfile}atoms.pdb xyzout ${tempfile}cryst.pdb << EOF > /dev/null
CELL $CELL
SPACE $SGnum
EOF
echo "REMARK B-factors are peak heights in $mapfile" >! $outfile
echo "REMARK Occs are 1/multiplicity (for special positions)" >> $outfile
egrep "^CRYST1" ${tempfile}cryst.pdb >> $outfile
cat ${tempfile}atoms.pdb >> $outfile
if(! $?DEBUG) rm -f ${tempfile}atoms.pdb ${tempfile}cryst.pdb

################################################################
# calculate distance to nearest "old" site

# by default, use the "boring" list as "old" sites
cat ${tempfile}boring_sites >! ${tempfile}old_sites
if(! -e "$sitefile") then
    # no sites input, so just get inter-peak distances
    cat ${tempfile}peaks.interesting >! ${tempfile}old_sites
endif
if($?scriptfile) then
    # input file was an mlphare script, so only interested in real atoms
    cat "$scriptfile" |\
    nawk '$1~/^ATOM/{print $3, $4, $5}' |\
    cat >! ${tempfile}old_sites
endif
# count number of "old" sites
set old_sites = `cat ${tempfile}old_sites | wc -l`

# calculate the largest expected inter-atom distance (cell center to origin)
echo "0.5 0.5 0.5" |\
nawk 'NF>2{++n; printf "%5d%10.5f%10.5f%10.5f%10.5f%5.2f%5d%10d%2s%3s%3s %1s\n", \
       n, $1, $2, $3, 80, 1, "38", n, "H", "", "IUM", " "}' |\
cat >! ${tempfile}.frac
coordconv XYZIN ${tempfile}.frac \
         XYZOUT ${tempfile}.pdb << EOF-conv >& /dev/null
CELL $CELL
INPUT FRAC
OUTPUT PDB ORTH 1
END
EOF-conv
cat ${tempfile}.pdb |\
nawk '/^ATOM/{print substr($0, 31, 8), substr($0, 39, 8), substr($0, 47, 8)}' |\
cat >! ${tempfile}center
if(! $?max_dist) set max_dist = `nawk '{print sqrt($1*$1 + $2*$2 + $3*$3)+3}' ${tempfile}center`
if(! $?DEBUG) rm -f ${tempfile}.frac >& /dev/null
if(! $?DEBUG) rm -f ${tempfile}.pdb >& /dev/null
if(! $?DEBUG) rm -f ${tempfile}center >& /dev/null

if($?FAST) then
    set max_dist = `echo $CLOSE_peaks 3 | awk '$2>$1{$1=$2} {print $1}'`
    set CLOSE_peaks = 0.0001
endif
if($?max_dist) then
    echo "not measuring distances > $max_dist A"
endif

# convert "old" sites to a pdb
cat ${tempfile}old_sites |\
nawk 'NF>2{++n; rn=n%10000; un=(n-rn)/10000\
    printf "%5d%10.5f%10.5f%10.5f%10.5f%5.2f%5d%10d%2s%3s%3s %1s\n", \
       n, $1, $2, $3, 80, un/100+1, "38", rn, "C", "", "IUM", "A"}' |\
cat >! ${tempfile}old.frac
coordconv XYZIN ${tempfile}old.frac \
         XYZOUT ${tempfile}old.pdb << EOF-conv | tee ${tempfile}coordconv.log >& /dev/null
CELL $CELL
INPUT FRAC
OUTPUT PDB ORTH 1
END
EOF-conv
if(! $?DEBUG) rm -f ${tempfile}old.frac >& /dev/null
if(! $?DEBUG) rm -f ${tempfile}coordconv.log >& /dev/null

# keep distang from crashing?
if($old_sites > 9999) then
    echo "WARNING: $old_sites is a lot of atoms! "
    cat ${tempfile}old.pdb |\
    awk '! /^ATOM|^HETAT/{print;next} {++n}\
      substr($0,23,4)~/[A-Z]/{\
       $0=substr($0,1,22) sprintf("%4d",n%10000) substr($0,27)}\
      {print}' |\
    cat >! ${tempfile}.pdb
    mv ${tempfile}.pdb ${tempfile}old.pdb
endif

# strip off unneeded cards
nawk '/^ATOM|^HETAT/ || /^CRYS/ || /^SCALE/' ${tempfile}old.pdb |\
cat >! ${tempfile}both.pdb
if(! $?DEBUG) rm -f ${tempfile}old.pdb >& /dev/null

# append peak list to the combined PDB file
nawk '/^ATOM|^HETAT/{print}' $outfile >> ${tempfile}both.pdb
@ start_peaks = ( $old_sites + 1 )

# renumber the atoms so that distang won't get confused
cat ${tempfile}both.pdb |\
nawk '/^ATOM|^HETAT/{++n; $0 = sprintf("ATOM  %5d%s",n,substr($0,12))} {print $0}' |\
cat >! ${tempfile}.pdb
mv ${tempfile}.pdb ${tempfile}both.pdb >& /dev/null

# use distang to calculate all inter-atom distances (and then sort them)
distang xyzin ${tempfile}both.pdb << EOF | tee ${tempfile}distang.log |\
    nawk '$1=="Z"{print "dist", $6, $2+($(NF-1)-1)*1000000, $9}' | sort -n --key=4 |\
    nawk '! seen[$2]{seen[$2]=1;print}' | sort -n --key=2 >! ${tempfile}dists
SYMM $SGnum
DIST ALL
RADII C 1
RADII OW $max_dist
DMIN $CLOSE_peaks
FROM ATOM 1 to $old_sites
TO   ATOM $start_peaks to 99999
END
EOF
if(! $?DEBUG) rm -f ${tempfile}both.pdb >& /dev/null
#${tempfile}dists now contains the minimum 
# format: peak# old# min_dist
# if -fast was set, missing distances are bigger than specified limit

# gaurentee a label file for the "old" sites
cat ${tempfile}old_sites |\
nawk 'NF>2{++n; printf "label %5d atom %d\n", n, n}' |\
cat >! ${tempfile}labels

# get a descriptive label from the "old site" source file
if(-e "$pdbfile") then
    # input file was a PDB file
    cat $sitefile |\
    nawk '/^ATOM/ || /^HETATM/{++n; printf "label %5d %s\n", n, substr($0,12,15)}' |\
    cat >! ${tempfile}labels
endif
if($?scriptfile) then
    # input file was an mlphare script
    cat $sitefile |\
    nawk '$1~/^DERIV/{deriv=$0}\
          $1~/^ATOM/{++n; printf "label %5d %6s in %s\n", n, $1, deriv}' |\
    cat >! ${tempfile}labels
endif
if(! -e "$sitefile") then
    # peaks list was used for distance self-caclulation
    cat ${tempfile}peaks.interesting |\
    nawk 'NF>2{++n;printf "label %5d peak %d\n", n, n}' |\
    cat >! ${tempfile}labels
endif

# add these descriptive labels to the peaks list
cat ${tempfile}dists ${tempfile}labels ${tempfile}peaks.interesting |\
nawk 'BEGIN{label[""]="unk"}\
      /^dist/{dist[$2]=$4; neighbor[$2]=$3; next}\
      /^label/{label[$2]=substr($0,13); next}\
      {++n} dist[n]==""{dist[n]="?"}\
      {print $0, dist[n], "label:", label[neighbor[n]]}' |\
cat >! ${tempfile}peaks.distlabel
#
# format: xf yf zf X Y Z 1/mult height   dist neighbor name ...

# additionally sort by distance and sigma?  For peaks with identical sigmas?
if($?FAR_FIRST) then
    sort -k8gr -k9gr ${tempfile}peaks.distlabel >! ${tempfile}sorted.distlabel
    head -n 3 $outfile >! ${tempfile}new.pdb
    cat ${tempfile}sorted.distlabel |\
    awk '{++i; printf "ATOM   %4d  OW  WAT X%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",\
             i, i, $4, $5, $6, $7, $8}' |\
    cat >> ${tempfile}new.pdb
    echo "END" >> ${tempfile}new.pdb
    mv ${tempfile}sorted.distlabel ${tempfile}peaks.distlabel
    mv ${tempfile}new.pdb $outfile
endif

if($?FAST) then
    # filter out too-close peaks now.
    mv $outfile ${tempfile}.pdb
    egrep "^CRYST" ${tempfile}.pdb >! $outfile
    mv ${tempfile}peaks.distlabel ${tempfile}
    cat ${tempfile} ${tempfile}.pdb |\
    awk -v mindist=$CLOSE_peaks '{id=substr($0,12,15)}\
      /^ATOM|^HETAT/ && ! sel[id]{print;next}\
      $9+0>mindist || $9=="?"{++sel[id]}' |\
    cat >> $outfile
    # now remove these labels from the output
    cat ${tempfile} |\
    awk -v mindist=$CLOSE_peaks '$9+0>=mindist || $9=="?"' |\
    cat >! ${tempfile}peaks.distlabel
endif

################################################################

# print surviving peaks out to screen
echo ""
echo "unique peaks:"
set xyz = "x        y        z"
if($?PATT) set xyz = "u        v        w"
echo "  $xyz        mult  height/sigma   dist  from nearest neighbor"
cat ${tempfile}peaks.distlabel |\
nawk -v max=$max_dist '{pdist=sprintf("%10.1f",$9)} $9=="?"{pdist=">"max}\
  {printf "%8.5f %8.5f %8.5f   %4d %8.2f %10sA  %s\n", $1, $2, $3, 1/$7, $8,\
       pdist, substr($0, index($0,"label:")+7)}'


echo "written to $outfile"
#exit

cleanup:
# clean up 
#exit
if($?BAD) exit 9
if(! $?DEBUG) rm -f  ${tempfile}* >& /dev/null

exit


################################################################
################################################################
################################################################
Setup:
set siteCELL
set sitefile

# scan command line
foreach arg ( $* )
    # recognize map files
    if(("$arg" =~ *.map)||("$arg" =~ *.ext)) then
        if(! -e "$arg") then
            echo "WARNING: $arg does not exist! "
            continue
        endif
        set mapfile = "$arg"

        continue
    endif
    
    # recognize pdb files
    if(("$arg" =~ *.pdb)||("$arg" =~ *.brk)) then
        if(! -e "$arg") then
            echo "WARNING: $arg does not exist! "
            continue
        endif
        set pdbfile = "$arg"
        set siteCELL = `nawk '/^CRYST/{print $2, $3, $4, $5, $6, $7}' $pdbfile | tail -1`
        
        continue
    endif
    
    # recognize mlphare scripts? 
    if(-e "$arg") then
        cat "$arg" |\
        nawk '$1~/^ATOM/ || $1~/^BADATOM/ || $1~/^OLDATOM/{print $3, $4, $5}' |\
        cat >! ${tempfile}boring_sites
        # format: xf yf zf
        set sitefile = "$arg"
        set scriptfile = "$arg"
    endif
    
    # recognize sigma-cutoff
    set temp = `echo "$arg" | awk '$1+0 != 0 && ! /A$/{print $1+0}'`
    if("$temp" != "") then
        set sigma = "$temp"
        continue
    endif

    # recognize closeness-cutoff
    set temp = `echo "$arg" | awk '$1+0 != 0 && /A$/{print $1+0}'`
    if("$temp" != "") then
        set CLOSE_peaks = "$temp"
        continue
    endif

    # recognize closeness-cutoff
    set temp = `echo "$arg" | awk -F "=" '/^max_dist/ && /=/ && /A$/{print $2+0}'`
    if("$temp" != "") then
        set max_dist = "$temp"
        continue
    endif

    # recognize "top peak only" flag
    if("$arg" =~ "-top"*) then
        set top_peaks = `echo $arg | awk -F "=" '{print $2+0}'`
        if("$top_peaks" == "") set top_peaks = 1
        if("$top_peaks" == "0") set top_peaks = 1
        continue
    endif

    # recognize "bottom peak only" flag
    if("$arg" =~ "-bottom"*) then
        set bottom_peaks = `echo $arg | awk -F "=" '{print $2+0}'`
        if("$bottom_peaks" == "") set bottom_peaks = 1
        if("$bottom_peaks" == "0") set bottom_peaks = 1
        continue
    endif

    # recognize "short distances only" flag
    if("$arg" =~ "-fast") then
        set FAST
        continue
    endif
end

if("$top_peaks" == "1") set top_only
if("$bottom_peaks" == "1") set bottom_only
set origsigma = "$sigma"

# get map parameters
echo "go" | mapdump mapin $mapfile >! ${tempfile}.mapdump
if($status) goto Help

set CELL = `nawk '/Cell dimensions/{print $4, $5, $6, $7, $8, $9; exit}' ${tempfile}.mapdump`
set SGnum = `nawk '/Space-group/{print $3; exit}' ${tempfile}.mapdump`
set SG = `nawk -v SGnum=$SGnum '$1==SGnum {print $4}' $CLIBD/symop.lib | head -1`
set GRID = `nawk '/Grid sampling on x, y, z/{print $(NF-2), $(NF-1), $NF; exit}' ${tempfile}.mapdump`
if(! $?DEBUG) rm -f ${tempfile}.mapdump >& /dev/null

# see if this is a Patterson map
set temp = `nawk -v SGnum="$SGnum" '$1==SGnum{print $4}' $CLIBD/symop.lib | nawk '/[abcdmn-]/{print "PATT"}'`
if("$temp" != "") set PATT


# convert input coordinate file formats to fractional
if($#siteCELL != 6) set siteCELL = ( $CELL )

if(-e "$pdbfile") then
    # convert orthogonal PDB coordinates to fractional
    coordconv xyzin $pdbfile xyzout ${tempfile}.xyz << EOF >> $logfile
CELL $siteCELL 
INPUT PDB
OUTPUT FRAC
END
EOF
    # all we need are fractional coordinates
    cat ${tempfile}.xyz |\
    nawk '{print $2, $3, $4}' |\
    cat >! ${tempfile}boring_sites
    if(! $?DEBUG) rm -f ${tempfile}.xyz >& /dev/null
    set sitefile = "$pdbfile"
endif
if(! -e "${tempfile}boring_sites") touch ${tempfile}boring_sites
set boring_sites = `cat ${tempfile}boring_sites | wc -l`
# format: xf yf zf 

# set the "close" criteria to be one grid unit
echo "$GRID $CELL" |\
nawk '$1+0>0{print $4/$1}\
      $2+0>0{print $5/$2}\
      $3+0>0{print $6/$3}' |\
sort -n >! ${tempfile}close
set CLOSE_map = `nawk 'NR==1{printf "%.2f", $1}' ${tempfile}close`
if(! $?DEBUG) rm -f ${tempfile}close >& /dev/null
if("$CLOSE_map" == "") set CLOSE_map = 0.5

# decide on a "closeness" cutoff for two peaks being the same
if(! $?CLOSE_peaks) set CLOSE_peaks
if("$CLOSE_peaks" == "") set CLOSE_peaks = $CLOSE_map
# guess? 
if("$CLOSE_peaks" == "") set CLOSE_peaks = 0.5


goto Return_from_Setup

exit
#################################
# the future? 
- support other coordinate file formats
