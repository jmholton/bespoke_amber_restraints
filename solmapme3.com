#! /bin/tcsh -f
#
#   make solvent map from amber trajectories
#
#
set trajs = ()
set lastn = "99999"
set topfile = ""
set tag = auto

set MD_mult = ( 1 1 1 )
set smallSG = P1
set GRID = ""

set scratch = "/dev/shm/${USER}"

set opts = "-nocache"

set pdir = `dirname $0`
set path = ( $path $pdir )

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($1)}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk -F "=" '{print $2}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set val = `echo $Val | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if("$Key" =~ *.nc && "$val" == "") then
        set trajs = ( $trajs $Key )
    endif
    if("$key" == "topfile") set topfile = "$Val"
    if("$key" == "lastn") set lastn = "$num"
    if("$key" == "tag") set tag = "$Val"

    if( "$arg" == "allatom") set ALLATOM
    if( "$arg" == "allatomonly") set ALLATOMONLY

    set opts = ( $opts $Arg )
end


if("$trajs" == "") set trajs = Prod.nc

set traj = ""
if( $#trajs > 0 ) then
  set traj = $trajs[1]
endif
if(! -e "$traj") then
    set BAD = "cannot read $traj"
    goto exit
endif

if(! -e "$topfile") set topfile = start.top
if(! -e "$topfile") set topfile = `ls -1rt *top | tail -n 1`
if(! -e "$topfile") set topfile = `ls -1rt *parm7 | tail -n 1`
if(! -e "$topfile") then
   set BAD = "no topfile"
   goto exit
endif
echo "using parmtop: $topfile"

set base = `echo $traj | awk '{while(gsub("[^.]$",""));gsub(".$","");print}'`
set num = `basename $traj | awk '{n=substr($0,match($0,/[0-9]/))+0} n>0{print n}'`
if( "$tag" == "auto" ) then
  set tag = "$num"
  if("$base" != "Prod" && "$base" !~ Prod[0-9]*) set tag = "_${base}"
  set tag = "_$base"
endif
echo "tag = $tag"

cat ${base}.out |\
awk '/ nstlim /{nstlim+=$3} \
     / dt /{dt=$6} \
     / ntwx /{ntwx=$6} \
     / begin time read / && $8!="ps"{t0=$8}\
     /NSTEP/{step=($6-t0)/dt+1;if(t0=="")t0=$6}\
     END{print dt*ntwx/1000,int(step/ntwx),step/nstlim}' |\
tee info.txt
set info = `cat info.txt`
if($#info != 3) then
    set BAD = "reading ${base}.out "
    goto exit
endif

set nsstep = $info[1]
@ nframes = ( $info[2] * $#trajs )
set fracdone = $info[3]

if($nframes < 50) then
    echo "WARNING: not enough frames, need 50"
endif

rm -rf ./trajectory
mkdir -p ${scratch}/trajectory_$$_/
ln -sf ${scratch}/trajectory_$$_/ trajectory
set pwd = `pwd`
echo "cpptraj ..."
echo "parm $topfile" >&! cpptraj.in
foreach traj ( $trajs )
  echo "trajin $traj" >> cpptraj.in
end
cat << EOF >> cpptraj.in
outtraj trajectory/md.pdb pdb multi pdbv3 keepext sg "P 1"
go
EOF
cpptraj < cpptraj.in >&! cpptraj.log
echo "status: $status"

echo "$MD_mult" >! MD_mult.txt
echo "$GRID" >! grid.txt

rm -rf cpu* maps/ >& /dev/null
ls -1 trajectory/md*.pdb |\
 sort -k1.15g |\
tail -n $lastn >! pdbfiles.txt
set npdbs = `cat pdbfiles.txt | wc -l`
if($npdbs < 50) then
    echo "WARNING: not enough pdb files, want 50"
endif

if( $?ALLATOMONLY) goto allatom

echo "md2map solvent on $npdbs pdb files ..."
md2map_multi.com $smallSG pdbfiles.txt -solventonly  $opts >&! md2map_solvent.log

if(! -e md_map_Bnorm.mtz) then
    cat md2map.log
    sleep 3600
    goto exit
endif
mv avg.map solvent${tag}.map
mv md_map_Bnorm.mtz solvent${tag}.mtz
ln -sf solvent${tag}.mtz solvent.mtz

if(! $?ALLATOM) goto exit

echo "clearing old scratch"
set cpudirs = `ls -l | awk '/^l/ && NF>7 && $(NF-2)~/^cpu/{print $NF}'`
rm -rf $cpudirs

allatom:
echo "md2map all-atom on $npdbs pdb files ..."
md2map_multi.com $smallSG pdbfiles.txt $opts -allatom >&! md2map_allatom.log

if(! -e md_map_Bnorm.mtz) then
    cat md2map.log
    sleep 3600
    goto exit
endif
mv avg.map allatom${tag}.map
mv md_map_Bnorm.mtz allatom${tag}.mtz


goto exit



set refmtz = ${pdir}/1aho.mtz


mkdir -p refmac${tag}
cad hklin1 $refmtz hklin2 solvent${tag}.mtz hklout refmac${tag}/refme.mtz << EOF >! refmacsetup.log
labin file 1 all
labin file 2 E1=FP E2=PHI
labou file 2 E1=Fpart E2=PHIpart
EOF

mkdir -p erefmac${tag}


cad hklin1 refmac${tag}/refme.mtz hklout temp.mtz << EOF >> refmacsetup.log
labin file 1 all
outlim space 1
EOF
cad hklin1 temp.mtz hklout P1.mtz << EOF >> refmacsetup.log
labin file 1 all
scale file 1 12 0
symm 1
EOF
echo reindex 2h,2k,3l |\
 reindex hklin P1.mtz hklout erefmac${tag}/refme.mtz >> refmacsetup.log
rm -f temp.mtz P1.mtz


ls -l refmac${tag}/refme.mtz
ls -l erefmac${tag}/refme.mtz


exit:
echo "cleaning up..."

set cpudirs = `ls -l | awk '/^l/ && NF>7 && $(NF-2)~/^cpu/{print $NF}'`
rm -rf $cpudirs
rm -rf cpu* maps/ trajectory header.map >& /dev/null
rm -rf ${scratch}/trajectory_$$_/ >& /dev/null

exit


mapmask mapin solvent${tag}.map mapout sfallme.map << EOF
xyzlim cell
axis Z X Y
EOF

sfall mapin sfallme.map hklout sfalled.mtz << EOF
mode sfcalc mapin
resolution $reso
sfsg 1
EOF


echo reindex h/2,k/2,l/3 |\
 reindex hklin erefmac_Prod/refme.mtz hklout temp.mtz
cad hklin1 temp.mtz hklout this.mtz << EOF
labin file 1 E1=Fpart E2=PHIpart
labou file 1 E1=Fthis E2=PHIthis
symm P212121
EOF
cad hklin1 refmac_Prod/refme.mtz hklout that.mtz << EOF
labin file 1 E1=Fpart E2=PHIpart
labou file 1 E1=Fthat E2=PHIthat
symm P212121
EOF
cad hklin1 this.mtz hklin2 sfalled.mtz hklout phases.mtz << EOF
labin file 1 E1=PHIthis
labin file 2 E1=PHIC
EOF
mtz2txt phases.mtz 

awk '$4~/[0-9]/ && $5~/[0-9]/{print sqrt(($4-$5)^2),$0}' phases.csh  | sort -g




echo "setting up refmac${tag} "
cd refmac${tag}
ln -sf ../../../1aho/vs_real/best_so_far.pdb start.pdb
cat << EOF >! refmac_opts.txt
LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag  FPART1=Fpart PHIP1=PHIpart
SCPART 1
SOLVENT NO
EOF
refmac_occupancy_setup_new.com start.pdb >> refmac_opts.txt

mkdir -p ../refmac${tag}/dry
cd ../refmac${tag}/dry
ln -sf ../refme.mtz .
grep -v HOH ../start.pdb >! start.pdb
cat << EOF >! refmac_opts.txt
LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag  FPART1=Fpart PHIP1=PHIpart
SCPART 1
SOLVENT NO
EOF
refmac_occupancy_setup_new.com start.pdb >> refmac_opts.txt


ls -l refme.mtz
# ssh graphics1 
# converge_refmac.com start.pdb refme.mtz nosalvage >&! converge1.log &

exit:
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

exit


./runme_allmaps.com


awk '/correct F:/{print $7,FILENAME}' */diff*.log | sort -gr
see diff_results_runme.com for diff_results.txt




diff.com ../1aho.mtz allatom${tag}.mtz
cad_phases.com scaleited.mtz allatom${tag}.mtz
fft hklin phase_cadded.mtz mapout ffted.map << EOF
labin F1=Fref PHI=PHI F2=Ftest
GRID 192 144 256
EOF
pick.com ffted.map -10 0.0001A ../best_so_far.pdb | tee pick${tag}.log


set pwd = `pwd`
cd $pwd
foreach mtz ( */allatom_Ramp*.mtz )
set dir = `dirname $mtz`
set n = `echo $mtz | awk '{gsub(".mtz","");n=substr($0,index($0,"allatom")+7);print n}'`
cd $pwd
cd $dir
echo "$dir $n $mtz"
if(-s diff${n}_3A.log) continue
pwd
diff.com ../1aho.mtz ../$mtz 3A | tee diff${n}_3A.log
diff.com ../1aho.mtz ../$mtz 2A | tee diff${n}_2A.log
diff.com ../1aho.mtz ../$mtz | tee diff${n}.log
mv diff_details.log diff_details${n}.log
cd $pwd
end


foreach log ( */diff*.log )
    if("$log" =~ *details*) continue
    set test = `grep "correct F:" $log | wc -l`
    if("$test" == "0") then
        echo $log
        rm -f $log
    endif
end



