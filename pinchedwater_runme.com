#! /bin/tcsh -f
#
#   look for waters with bad non-bonds  these will probably explode amber
#   also check for vacuum-filled voids
#
set pdbfile = ""
set restraints = ""

set nonbonds = ""
set outprefix = ""

set water_radius  = 2.4
set energy_thresh = 30
set max_rejects   = 1000

set tempfile = /dev/shm/${USER}/temp_up_$$_
mkdir -p /dev/shm/${USER}
setenv CCP4_SCR ${tempfile}_dir/
mkdir -p $CCP4_SCR

# scan command line for settings
foreach Arg ( $* )
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk -F "=" '{print $2}'`
    set assign = `echo $Arg | awk '{print ( /=/ )}'`
    set arg = `echo $Arg | awk '{print tolower($0)}'`

    if( $assign ) then
        set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
        if ( $test ) then
           set $Key = $Val
           echo "$Key = $Val"
        endif
    else
        if( "$Arg" =~ *.pdb ) then
          if( "$pdbfile" == "" ) then
            set pdbfile = $Arg
          else
            set restraints = $Arg
          endif
        endif
        if( "$arg" == "debug" ) set debug
    endif
end

if( $?DEBUG ) set debug
if( $?debug ) set tempfile = tempfile


cat << EOF
pdbfile = $pdbfile
restraints = $restraints
outprefix  = $outprefix
water_radius = $water_radius

nonbonds       = $nonbonds
energy_thresh  = $energy_thresh
max_rejects    = $max_rejects

tempfile       = $tempfile
EOF


set t = ${tempfile}

convert_pdb.awk -v skip=H $pdbfile |\
 awk '/^CRYST/ && ! c{print;++c;next}\
 ! /^ATOM|^HETAT/{next} \
 {print "ATOM  ",substr($0,8,73)}' >! ${t}_noH.pdb


echo "looking for voids"
sfall xyzin ${t}_noH.pdb mapout ${t}voidme.map << EOF >! ${t}sfall.log
MODE ATMMAP SOLVMAP
SYMM 1
grid 128 128 128
vdwrad $water_radius
EOF
if($status) then
   mv ${t}_noH.pdb sfallme.pdb
   mv ${t}sfall.log sfall.log
   goto geom
endif
float_func -func segment -header 1104 -xsize 128 -ysize 128 \
  ${t}voidme.map -outfile ${t}segments.map >! ${t}segment.log

float_add -histogram -header 1104 ${t}segments.map -outfile ${t}.bin |\
 tee ${t}segment.log |\
 awk 'NF==2 && $1+0!=0 && $2+0!=0{print $2}' |\
 sort -gr | head -n 1 >! ${t}biggest_void.txt

set void = `awk '{print $1;exit}' ${t}biggest_void.txt`
if("$void" == "") set void = 0
echo "biggest void: $void voxels"


map_scaleB_runme.com ${t}voidme.map B=50 outfile=${t}scaled.map  >> /dev/null
echo xyzlim -5 133 -5 133 -5 133 |\
 mapmask mapin ${t}scaled.map mapout ${outprefix}extended.map >> /dev/null

geom:
if(! -e "$nonbonds") then
echo "phenix geo"
  phenix.geometry_minimization macro_cycles=0 ${t}_noH.pdb \
    output_file_name_prefix=${t}phenix >! ${t}phenix_geom.log 

  cat ${t}phenix.geo |\
  awk '/nonbonded pdb=/{key="NONBOND";split($0,w,"\"");id1=w[2];\
       getline;split($0,w,"\"");id2=w[2];\
       getline;getline;\
       obs=$1;ideal=$2;sigma=1;energy=lj(obs,ideal);}\
    id1!=""{\
       a1=substr(id1,1,4);a2=substr(id2,1,4);a3=substr(id3,1,4);a4=substr(id4,1,4);\
       f1=substr(id1,5,1);f2=substr(id2,5,1);f3=substr(id3,5,1);f4=substr(id4,5,1);\
       t1=substr(id1,6,4);t2=substr(id2,6,4);t3=substr(id3,6,4);t4=substr(id4,6,4);\
       c1=substr(id1,10,1);c2=substr(id2,10,1);c3=substr(id3,10,1);c4=substr(id4,10,1);\
       r1=substr(id1,11,5);r2=substr(id2,11,5);r3=substr(id3,11,5);r4=substr(id4,11,5);\
       gsub(" ","_",f1);gsub(" ","_",f2);gsub(" ","_",f3);gsub(" ","_",f4);\
       gsub(" ","_",c1);gsub(" ","_",c2);gsub(" ","_",c3);gsub(" ","_",c4);\
       nid1=sprintf("%4s %s %4s %s %5s",a1,f1,t1,c1,r1);\
       nid2=sprintf("%4s %s %4s %s %5s",a2,f2,t2,c2,r2);\
       nid3=sprintf("%4s %s %4s %s %5s",a3,f3,t3,c3,r3);\
       nid4=sprintf("%4s %s %4s %s %5s",a4,f4,t4,c4,r4);\
       ids=nid1; if(nid2~/[0-9]/)ids=ids" - "nid2; if(nid3~/[0-9]/)ids=ids" - "nid3;if(nid4~/[0-9]/)ids=ids" - "nid4;\
       print key,energy+0,ideal-obs,obs,ideal,sigma+0,"|",ids;\
       id1=id2=id3=id4=""} \
     function lj0(r,r0) {if(r==0)return 1e40;return 4*((r0*2**(-1./6)/r)**12-(r0*2**(-1./6)/r)**6)}\
   function lj(r,r0) {return lj0(r,r0)-lj0(6,r0)}' |\
  sort -k2gr |\
  awk '$2>0.1' >! ${outprefix}nonbond_sorted.txt

  set nonbonds = ${outprefix}nonbond_sorted.txt
endif

echo "using nonbonds in: $nonbonds"
#grep HOH $nonbonds | head

if(-e "$restraints") then
  echo "keeping all atoms in $restraints regardless of pinching"
  cat "$restraints" |\
  awk '/^ATOM|^HETAT/{print substr($0,18,3),substr($0,22,1),substr($0,23,4)}' |\
  cat >! ${t}restraints.txt
else
  touch ${t}restraints.txt
endif

cat ${t}restraints.txt $nonbonds |\
awk 'NF==3{++restr[$1,$2,$3];next}\
 $8~/^H/ || $14~/^H/{next}\
 $2<0{exit}\
 ! restr[$16,$17,$18] && $16=="HOH" && ! rejected[$11,$12]{print;++rejected[$17,$18]}\
 ! restr[$10,$11,$12] && $10=="HOH" && ! rejected[$17,$18]{print;++rejected[$11,$12]}' |\
head


cat $nonbonds |\
awk -v thresh=$energy_thresh 'NF==3{++restr[$1,$2,$3];next}\
 $2<thresh{exit}\
 ! restr[$16,$17,$18] && $16=="HOH" && ! rejected[$11,$12]{print $2,$17,$18;++rejected[$17,$18]}\
 ! restr[$10,$11,$12] && $10=="HOH" && ! rejected[$17,$18]{print $2,$11,$12;++rejected[$11,$12]}' |\
awk '{sum[$2" "$3]+=$1}\
  END{for(cr in sum)print sum[cr],cr}' |\
sort -gr >! ${t}badwaters.txt
set badwaters = `cat ${t}badwaters.txt | wc -l`
echo "$badwaters waters are pinched with energy threshold > $energy_thresh"

cat ${t}restraints.txt $nonbonds |\
awk -v thresh=$energy_thresh 'NF==3{++restr[$1,$2,$3];next}\
 $2<thresh{exit}\
 restr[$16,$17,$18] && $16=="HOH"{print $17,$18}\
 restr[$10,$11,$12] && $10=="HOH"{print $11,$12}' |\
sort -u >! ${outprefix}restraints_and_pinched.txt
set pinrest = `cat ${outprefix}restraints_and_pinched.txt | wc -l`
echo "$pinrest restrained waters are pinched"

cat ${t}restraints.txt $nonbonds |\
awk -v thresh=$energy_thresh 'NF==3{++restr[$1,$2,$3];next}\
 $2<thresh{exit}\
 ! restr[$16,$17,$18] && $16=="HOH" && ! rejected[$11,$12]{print $2,$17,$18;++rejected[$17,$18]}\
 ! restr[$10,$11,$12] && $10=="HOH" && ! rejected[$17,$18]{print $2,$11,$12;++rejected[$11,$12]}' |\
awk '{sum[$2" "$3]+=$1}\
  END{for(cr in sum)print sum[cr],cr}' |\
sort -gr >! ${outprefix}badwaters.txt
set badwaters = `cat ${outprefix}badwaters.txt | wc -l`
echo "$badwaters unrestrained waters are pinched with energy threshold > $energy_thresh"

cat ${outprefix}badwaters.txt ${t}_noH.pdb |\
awk 'NF==3{++sel[$2,$3];next}\
  ! /^ATOM|^HETAT/{print;next}\
 {chain=substr($0,22,1);resnum=substr($0,23,8)+0}\
 sel[chain,resnum]{print}' |\
cat >! ${outprefix}badwaters.pdb

set badwaters = `egrep "O   HOH" ${outprefix}badwaters.pdb | head -n $max_rejects | wc -l`

head -n $max_rejects ${outprefix}badwaters.txt |\
cat - $pdbfile |\
awk 'NF==3{++elim[$2,$3];next}\
  ! /^ATOM|^HETAT/{print;next}\
 {chain=substr($0,22,1);resnum=substr($0,23,8)+0}\
 ! elim[chain,resnum]{print}' |\
cat >! ${outprefix}unpinched.pdb

wc -l ${outprefix}badwaters.txt 
echo -n "badwaters: "
egrep "O   HOH" ${outprefix}badwaters.pdb | wc -l
set waterin = `egrep "O   HOH" $pdbfile | wc -l`
set waterout = `egrep "O   HOH" ${outprefix}unpinched.pdb | wc -l`

@ removed = ( $waterin - $waterout )

echo "removed $removed out of $waterin waters, leaving $waterout"

exit:
if( $?BAD ) then
  echo "ERROR: $BAD"
  exit 9
endif

if(! $?debug) then
  echo "clearng temp files"
  rm -rf ${t}*
endif

exit




foreach test ( test8 test9 test10_verywet test11 test12 )
cd ../$test
foreach itr ( `seq 1 300` )
if(! -e wrapped_${itr}.pdb) continue

#append_file_date.com pinchedwater_${itr}.log 

set water_radius = `awk '/^water_radius/{print $NF}' pinchedwater_${itr}.log | tail -n 1`

$sruncpu ../pinchedwater_runme.com wrapped_${itr}.pdb \
     restraints=restraints_for_${itr}.pdb \
     water_radius=$water_radius energy_thresh=10 max_rejects=1000 \
     nonbonds=itr${itr}_nonbond_sorted.txt \
     outprefix=itr${itr}_ |\
  cat >! pinchedwater_${itr}.log &

sleep 0.2

end
end


