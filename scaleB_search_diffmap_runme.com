#! /bin/tcsh -f
#
#  hone in on optimal scale and B with direct map comparison and outlier rejection
#
#
set mtz2 = "$1"

set Brange = 100
set Bsteps = 10
set scalerange = 2
set scalesteps = 10
set max_rounds = 10

set mtz1 = reference.mtz
if(! -e "$mtz2") then
   set BAD = "need two mtzs"
   goto exit
endif

set pwd = `pwd`
set tempfile = /dev/shm/${USER}/tmpsBs_$$_dir/
mkdir -p $tempfile
set t = $tempfile

cp $mtz1 ${t}ref.mtz
cp $mtz2 ${t}test.mtz
cd $t

diff.com ref.mtz test.mtz | tee diff.log
set scale0 = `awk '/scale=/{s=$2;B=$4;print $2}' diff.log`
set B0 = `awk '/scale=/{s=$2;B=$4;print $4}' diff.log`

cad hklin1 ref.mtz hklin2 test.mtz hklout cadded.mtz << EOF
labin file 1 E1=Fref E2=PHIref
labin file 2 E1=FCavg E2=PHICavg
labou file 1 E1=Fref E2=PHIref
labou file 2 E1=Ftest E2=PHItest
EOF

fft hklin cadded.mtz mapout ref.map << EOF >! fft.log
labin F1=Fref PHI=PHIref
EOF

echo go | mapdump mapin ref.map | tee mapdump.log |\
awk '/Grid sampling on x, y, z/{gx=$8;gy=$9;gz=$10}\
     /Maximum density /{max=$NF}\
     /Cell dimensions /{xs=$4/gx;ys=$5/gy;zs=$6/gz}\
     /Number of columns, rows, sections/{nc=$7;nr=$8;ns=$9}\
 END{print xs,ys,zs,nc,nr,ns,max}' >! mapstuff.txt

set voxels = `awk '{print $4*$5*$6}' mapstuff.txt`
set size = `echo $voxels | awk '{print 4*$1}'`
set head = `ls -l ref.map | awk -v size=$size '{print $5-size}'`
set skip = `echo $head | awk '{print $1+1}'`
set offset = `awk '/Mean density/{print $NF}' mapdump.log`
set GRID = `awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}' mapdump.log`

set round = 0
again:
@ round = ( $round + 1 )

echo "scalerange = $scalerange in $scalesteps steps"
echo "Brange = $Brange in $Bsteps steps"

set dBs = `echo $Brange $Bsteps | awk '{r=$1;d=2*r/$2;for(B=-r;B<=r;B+=d)print B}'`
set dscales = `echo $scalerange $scalesteps | awk '{r=$1;d=exp(log($1)/($2-1)*2);for(s=1/r;s<=r;s*=d)print s}'`
echo dBs: $dBs
echo dscales: $dscales 
foreach dB ( $dBs )

set B = `echo $dB $B0 | awk '{print $1+$2}'`

fft hklin cadded.mtz mapout test.map << EOF >! fft.log
labin F1=Ftest PHI=PHItest
scale F1 1 $B
GRID $GRID
EOF

echo 
foreach dscale ( $dscales )

set scale = `echo $dscale $scale0 | awk '{print $1*$2}'`

float_add -header $skip -scale1 -1 -scale2 $scale ref.map test.map -outfile /dev/null -reject |\
  tee float_add_${B}_${scale}.log | egrep "^after" &

end
end
wait

awk '/pixels rejected/{rej=$1}\
  /^after/{split(FILENAME,w,"_");print w[4]+0,w[3]+0,$9,rej}' float_add* |\
 sort -k3g |\
 tee fitme_raw.txt |\
 awk 'NF<4{next}\
   ! seen[$1]{++seen[$1];++uniqx}\
   ! seen[$2]{++seen[$2];++uniqy}\
   {print}\
   uniqx>6 && uniqy>6 && NR>20{exit}' |\
 tee  fitme
 
set r0 = `head -n 1 fitme | awk '{print $3;exit}'`
set sstats = `awk '{print $1}' fitme | avg.awk | awk '{print $1,$3*20}' | tail -n 1`
set Bstats = `awk '{print $2}' fitme | avg.awk | awk '{print $1,$3*20}' | tail -n 1`

cat << EOF >! fitparams.txt
r0 = $r0
s0 = $sstats[1]
sw = $sstats[2]
B0 = $Bstats[1]
Bw = $Bstats[2]
EOF

cat << EOF >! gnuplot.in
r(s,B) = r0+((B-B0)/Bw)**2+((s-s0)/sw)**2
fit r(x,y) 'fitme' using 1:2:3:(1) via r0
fit r(x,y) 'fitme' using 1:2:3:(1) via r0,Bw,sw
fit r(x,y) 'fitme' using 1:2:3:(1) via r0,Bw,sw,B0,s0
update 'fitparams.txt'
EOF

set scale = `awk '/^s0/{print $3}' fitparams.txt`
set B = `awk '/^B0/{print $3}' fitparams.txt`
echo "best-fit scale,B = $scale $B"

set scale0 = $scale
set B0 = $B

echo $scale $B |\
cat - fitme |\
awk 'NR==1{S0=$1;B0=$2;next}\
  {++n;gtS=gtB=0} $1>S0{++hiS;++gtS} $2>B0{++hiB;++gtB}\
  ! seenS[$1]{++seenS[$1];++hiloS[gtS]}\
  ! seenB[$2]{++seenB[$2];++hiloB[gtB]}\
  END{print hiS/n,hiB/n,hiloS[0]+0,hiloS[1]+0,hiloB[0]+0,hiloB[1]+0}' |\
tee fracs.txt
  
set tests = `awk '{print ( $1>0.25 && $1<0.75 && $2>0.25 && $1<0.75 && $3>2 && $4>2 && $5>2 && $6>2 ) }' fracs.txt`
if ( "$tests" != "0" && $round < $max_rounds ) then
  echo "another round..."
#  set Bstep = `echo $Bstep 2 | awk '{print $1/$2}'`
  set test = `awk '{print ( $5>2 && $6>2 ) ; exit}' fracs.txt`
  if( $test ) set Brange = `echo $Brange | awk '{print $1/1.5}'`
  set test = `awk '{print ( $3>2 && $4>2 ) ; exit}' fracs.txt`
  if( $test ) set scalerange = `echo $scalerange | awk '$1>1.1{$1=$1/1.1} {print $1}'`
  echo "new ranges: s: $scalerange B: $Brange"
  goto again
endif

cad hklin1 cadded.mtz hklout scaled.mtz << EOF
labin file 1 E1=Ftest E2=PHItest
scale file 1 $scale $B
EOF

fft hklin cadded.mtz mapout diff.map << EOF >! fft.log
labin F1=Fref PHI=PHIref F2=Ftest PHI2=PHItest
scale F2 $scale $B
GRID $GRID
EOF

rm -f cootme-grid.mtz
sftools << EOF
read ref.mtz
read scaled.mtz
map correl Fref PHIref Ftest PHItest
calc ( COL DELFWT PHDELWT ) = ( COL Fref PHIref ) ( COL Ftest PHItest ) -
calc ( COL FWT PHWT ) = ( COL Fref PHIref )
calc ( COL FC PHIC ) = ( COL Ftest PHItest )
write cootme-grid.mtz col FWT PHWT DELFWT PHDELWT FC PHIC
quit
EOF


diff.com ref.mtz scaled.mtz -noscale

echo "did $round rounds"
set n = `cat fitme | wc -l`
echo -n "fraction of top-$n points above optima: "
cat fracs.txt
echo "best-fit scale,B = $scale $B"

cp diff.map cootme-grid.mtz fitme* fitparams.txt $pwd

cd $pwd

if( "$t" == "" || "$t" == ".") exit

rm -rf ${t}

exit





fft hklin cootme.mtz mapout ref.map << EOF 
labin F1=FWT PHI=PHWT
GRID $GRID
EOF

fft hklin cootme.mtz mapout scaled.map << EOF 
labin F1=FC PHI=PHIC
GRID $GRID
EOF




