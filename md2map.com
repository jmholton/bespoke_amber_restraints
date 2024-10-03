#! /bin/tcsh -f
#
#	c-shell script outline for converting a molecular dynamics trajectory into structure factors
#	using the CCP4 Suite of programs
#
#
#
set CELL = ( 43.000   52.610   89.120  90.00  90.00  90.00 )
set CELL
set MD_mult = ( 1 1 1 )
set SG = P1
set reso = 1.0
#set GRID = "GRID 108  128  216"
set GRID
if(-e cell.txt) set CELL = `cat cell.txt`
if(-e grid.txt) set GRID = `cat grid.txt`
if(-e MD_mult.txt) set MD_mult = `cat MD_mult.txt`
set B = 5
set tempfile = /tmp/${USER}/tempfile$$.txt
mkdir -p /dev/shm/${USER}
if(-w /dev/shm/${USER}) set tempfile = /dev/shm/${USER}/tempfile$$_
echo -n "" >! ${tempfile}pdbfiles.txt

mkdir -p maps

echo "args: $*"
foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set ARG = `echo $Arg | awk '{print toupper($0)}'`
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
    endif

    if(("$arg" =~ *.txt)) then
        # warn about probable mispellings
        if(! -e "$arg") then
            echo "WARNING: $arg does not exist"
            continue
        endif
        
        cat $arg >> ${tempfile}pdbfiles.txt
        continue
    endif

    if(("$arg" =~ *.pdb)) then
        # warn about probable mispellings
        if(! -e "$arg") then
            echo "WARNING: $arg does not exist"
            continue
        endif
        
        echo "$arg" >> ${tempfile}pdbfiles.txt
        continue
    endif

    # change B factor
    if( "$arg" =~ B=* ) then
        set B = `echo "$arg" | awk -F "=" '{print $2+0}'`
    endif

    # change space group of comparison
    if( "$ARG" =~ [PCIFRH][1-6]* ) then
        set SG = "$ARG"
    endif

    # specify resolution
    if( "$arg" =~ *A ) then
        set reso = "$num"
    endif

    # allow atom type selections
    if( "$arg" == "-wateronly" || "$arg" == "-solventonly" ) then
        set WATERONLY
    endif
    if( "$arg" == "-proteinonly" ) then
        set PROTEINONLY
    endif
    if( "$arg" == "-nocache" ) then
        set NOCACHE
    endif
    if( "$arg" == "-allatoms" || "$arg" == "-allatom" ) then
        unset PROTEINONLY
        unset WATERONLY
    endif
end
if($?WATERONLY && $?PROTEINONLY) then
    unset WATERONLY
    unset PROTEINONLY
endif
if( ! $?WATERONLY && ! $?PROTEINONLY ) then
    echo "all atoms selected"
endif
if($?WATERONLY) echo "solvent atoms selected"
if($?PROTEINONLY) echo "non-solvent atoms selected"


# get number of ASUs in a unit cell
set symops = `awk -v SG=$SG '$4==SG{print $2;exit}' ${CLIBD}/symop.lib`
echo "$symops ASU/cell in $SG"

if(! -s "${tempfile}pdbfiles.txt") then
    set BAD = "please specify a pdb file"
    goto exit
#    ls -1 md[0-9]*.pdb | sort -k1.3g >! ${tempfile}pdbfiles.txt
endif
set frames = `cat ${tempfile}pdbfiles.txt | wc -l`
set frame = 0
set avg_step = 99999



# run through all the "frames" in the simulation 
set start_time = `echo "puts [clock clicks -milliseconds]" | tclsh`
rm -f sum.map
while ( $frame < $frames )
  @ frame = ( $frame + 1 )
  set file = `head -n $frame ${tempfile}pdbfiles.txt | tail -1`
  set pframe = `echo $file | awk '{file=$0;while(gsub("[0-9]$|.pdb$",""));printf("%05d",substr(file,length($0)+1)) }'`

  set weight = `egrep "^${pframe} " pdb_weights.txt |& awk '/ weight /{w+=$3} END{print w}'`
  if("$weight" == "") then
      set weight = `egrep "^${pframe} " ../pdb_weights.txt |& awk '/ weight /{w+=$3} END{print w}'`
  endif
  if("$weight" == "") set weight = 1

  echo "adding ${weight}x $file"

  if( $?NOCACHE ) then
     rm -f maps/sfalled_${pframe}.map >& /dev/null
  endif
  if(! -e maps/sfalled_${pframe}.map) then

  set shift = `egrep "^${pframe} " pdb_shifts.txt |& awk '/ shift /{x+=$3;y+=$4;z+=$5} END{print x,y,z}'`
  if("$shift" == "") then
      set shift = `egrep "^${pframe} " ../pdb_shifts.txt |& awk '/ shift /{x+=$3;y+=$4;z+=$5} END{print x,y,z}'`
  endif
  set scale = `egrep "^${pframe} " pdb_scales.txt |& awk '/ scale /{x+=$3;y+=$4;z+=$5} END{print x,y,z}'`
  if("$scale" == "") then
      set scale = `egrep "^${pframe} " ../pdb_scales.txt |& awk '/ scale /{x+=$3;y+=$4;z+=$5} END{print x,y,z}'`
  endif
  set cell = `egrep "^${pframe} " pdb_cells.txt |& awk '/ cell /{print $3,$4,$5,$6,$7,$8}'`
  if("$cell" == "") then
      set cell = `egrep "^${pframe} " ../pdb_cells.txt |& awk '/ cell /{print $3,$4,$5,$6,$7,$8}'`
  endif

  # awk program to convert file format
  # note that "HG11" will look like mercury to SFALL if it is not aligned properly
  # for "safety", we will convert all atoms into just their element names
cat << EOF | tee ${tempfile}params.txt
MDMULT $MD_mult
SHIFT $shift
SCALE $scale
CELL  $cell
BFAC $B
EOF
  cat ${tempfile}params.txt $file |\
  awk 'BEGIN{sx=sy=sz=1;B=20} \
       /^MDMULT/{na=$2;nb=$3;nc=$4;next}\
       /^SHIFT/{dx=$2;dy=$3;dz=$4;next}\
       /^SCALE/ && $2*$3*$4>0{sx=$2;sy=$3;sz=$4;next}\
       /^BFAC/{B=$2;next}\
       /^CELL/{a=$2;b=$3;c=$4;al=$5;be=$6;ga=$7;next}\
      /^CRYST/{if(! a){a=$2/na*sx;b=$3/nb*sy;c=$4/nc*sz;al=$5;be=$6;ga=$7}\
        printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",a,b,c,al,be,ga}\
      /^ATOM|^HETAT/{\
        RES= substr($0, 18, 3);\
        X = (substr($0, 31, 8)+dx)*sx;\
        Y = (substr($0, 39, 8)+dy)*sy;\
        Z = (substr($0, 47, 8)+dz)*sz;\
#	while(X>a){X-=a};while(X<0){X+=a};\
#	while(Y>b){Y-=b};while(Y<0){Y+=b};\
#	while(Z>c){Z-=c};while(Z<0){Z+=c};\
       atom=substr($0,12,5);gsub(" ","",atom)\
       Ee = $NF;\
       if(atom=="SE")Ee=atom;\
    printf("ATOM   %4d  %-2s  %3s     1    %8.3f%8.3f%8.3f  1.00%6.2f%12s\n",++n%10000,Ee,RES,X,Y,Z,B,Ee);}\
    END{print "END"}' |\
  cat >! sfallme.pdb

  if($?WATERONLY) then
    egrep "^CRYST| HOH | WAT | AMM | ACT | CL | LI | EDO " sfallme.pdb >! new.pdb
    mv new.pdb sfallme.pdb
  endif
  if($?PROTEINONLY) then
    egrep -v " WAT | HOH | AMM | ACT | CL | LI | EDO " sfallme.pdb >! new.pdb
    mv new.pdb sfallme.pdb
  endif

  # sfall needs "confirmation" of cell
  set cell = `awk '/^CRY/{print $2,$3,$4,$5,$6,$7;exit}' sfallme.pdb`

if( $?GOSLOW ) then
  pdbset xyzin sfallme.pdb xyzout new.pdb << EOF > /dev/null
  SPACE $SG
  CELL $cell
EOF
  mv new.pdb sfallme.pdb

  coordconv xyzin sfallme.pdb xyzout new.xyz << EOF > /dev/null
INPUT PDB
OUTPUT FRAC
END
EOF
  cat new.xyz |\
  awk '{fx=$2;while(fx<0)++fx;while(fx>1)--fx;\
        fy=$3;while(fy<0)++fy;while(fy>1)--fy;\
        fz=$4;while(fz<0)++fz;while(fz>1)--fz;\
   printf("%5d%10.5f%10.5f%10.5f%s\n",\
      $1,fx,fy,fz,substr($0,36))}' |\
  cat >! sfallme.xyz
  coordconv xyzin sfallme.xyz xyzout sfallme.pdb << EOF > /dev/null
CELL $cell
INPUT FRAC
OUTPUT PDB
END
EOF
endif

  setenv MEMSIZE `echo $cell $reso | awk '{print int($1*$2*$3/($NF**3)*100)}'`

  # actual run of sfall, perhaps use a slightly finer grid that default?
  echo "sfall"
  sfall xyzin sfallme.pdb mapout maps/sfalled_${pframe}.map << EOF | tee sfall.log > /dev/null
  mode atmmap
  CELL $cell
  SYMM $SG
  FORMFAC NGAUSS 5
  $GRID
  RESO $reso
EOF
  endif

  rm -f ${tempfile}weighted.map
  
  echo "scale factor $weight 0" |\
  mapmask mapin maps/sfalled_${pframe}.map mapout ${tempfile}weighted.map > /dev/null

  # sum the electron density
  echo "mapmask" 
  if(-e sum.map) then
     rm -f new.map
     echo "maps add" |\
     mapmask mapin1 sum.map mapin2 ${tempfile}weighted.map mapout new.map > /dev/null
     mv new.map sum.map
     set norm = `echo $norm $weight | awk '{print $1+$2}'`
  else
     cp ${tempfile}weighted.map sum.map
     set norm = $weight
     if("$GRID" == "") then
        set GRID = `echo "" | mapdump mapin sum.map | awk '/Grid sampling on x, y, z/{print "GRID",$(NF-2), $(NF-1), $NF; exit}'`
        echo "selected $GRID"
     endif
  endif

  if($frame % $avg_step == 0) then
      @ avg_step = ( 2 * $avg_step )
      set ASUs = `echo $frame $MD_mult $symops | awk '{print $1*$2*$3*$4*$5}'`
      set scale = `echo $ASUs | awk '{print 1/$1}'`
      echo "scale factor $scale 0" |\
      mapmask mapin sum.map mapout avg.map > /dev/null
      echo | mapdump mapin avg.map | egrep dens
  endif

  set now_time = `echo "puts [clock clicks -milliseconds]" | tclsh`
  set sofar = `echo $now_time $start_time | awk '{print ($1-$2)/1000}'`
  set togo  = `echo $sofar $frame $frames | awk '{print int($1/$2*$3-$1)}'`
  set finish = `echo "puts [clock format [expr [clock seconds] + $togo]]" | tclsh`
  echo "expect to finish at $finish"
end



echo "computing average"
set ASUs = `echo $frame $MD_mult $symops | awk '{print $1*$2*$3*$4*$5}'`
echo "$frame x $MD_mult x $symops = $ASUs ASUs"
set scale = `echo $ASUs $frame $norm | awk '{print 1/$1*($2/$3)}'`
echo "scale = $scale after weighting"
echo "scale factor $scale 0" |\
mapmask mapin sum.map mapout avg.map > /dev/null
echo | mapdump mapin avg.map | egrep dens

echo scale sigma |\
mapmask mapin avg.map mapout sigscale.map > /dev/null
echo | mapdump mapin sigscale.map | egrep dens


rm -f md_map_avg.mtz
sfall mapin avg.map hklout sfall.mtz << EOF | tee lastsfall.log > /dev/null
MODE SFCALC MAPIN
EOF
sftools << EOF > /dev/null
read sfall.mtz
set labels
FP
PHI
calc Q COL SIGFP = 0.1
calc W COL FOM = 0.9
write md_map_avg.mtz col FP SIGFP PHI FOM
y
stop
EOF

exit:

rm -f ${tempfile}pdbfiles.txt
rm -f ${tempfile}params.txt
rm -f ${tempfile}*

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif


exit

ls -1 md*.pdb | sort -k1.3g | awk '{print "../"$1}' >! pdbfiles.txt

set cpus = 24
foreach cpu ( `seq -f%03.0f 1 $cpus` )
mkdir -p cpu$cpu
cd cpu$cpu
awk -v cpu=$cpu -v cpus=$cpus 'NR%cpus==cpu-1' ../pdbfiles.txt >! pdbfiles.txt
cd ..
end
wc cpu*/pdbfiles.txt
rm cpu*/busy
rm cpu*/*map

set jobs = `grep processor /proc/cpuinfo | wc -l | awk '{print $1-1}'`
set pwd = `pwd`
foreach cpu ( `seq -f%03.0f 1 30` )
cd $pwd
if(-e cpu${cpu}/busy) continue
ls cpu${cpu}/ > /dev/null
if(-e cpu${cpu}/busy) continue
cd cpu$cpu
echo "`hostname -s` $$" >! busy
../md2map.com pdbfiles.txt |& log_timestamp.tcl >! md2map.log &
@ jobs = ( $jobs - 1 )
if($jobs == 0) break
end
pwd
tail -f md2map.log


mv cpu*/sfall*.map .



foreach num ( 500 200 100 )

ls -1 md*.pdb | sort -k1.3g | tail -n $num >! pdbfiles.txt

./md2map.com pdbfiles.txt |& log_timestamp.tcl
mv avg.map avg_${num}.map
mv md_map_avg.mtz md_map_avg_${num}.mtz 

end


./md2mtz.com pdbfiles.txt |& log_timestamp.tcl | tee md2mtz_${num}.log
mv F_Favg.mtz F_Favg_${num}.mtz
mv F_Iavg.mtz F_Iavg_${num}.mtz

end



