#! /bin/tcsh -f
#
#

#set pdbid = 6c2r

#phenix.fetch_pdb ${pdbid}


set prefix = Dirk

phenix.ready_set ${prefix}.pdb

set ligcifs = ( ???.cif )

echo 1 | mmtbx.quantum_interface iterate_NQH=HIS ${prefix}.updated.pdb | tee options_${prefix}.log

set opts = `awk '/resname HIS/ && $2==":"{print $1}' options_${prefix}.log`

set confs = `awk '/^ATOM|^HETAT/{print substr($0,17,1)}' ${prefix}.pdb | sort -u`

foreach conf ( $confs )
cat ${prefix}.pdb |\
awk -v conf=${conf} '/^CRYST/{print;next}\
  ! /^ATOM|^HETAT/{next}\
  {f=substr($0,17,1)}\
  conf==f || f==" "{print substr($0,1,16),substr($0,18)}' |\
cat >! ${prefix}${conf}.pdb

phenix.ready_set ${prefix}${conf}.pdb
end

foreach conf ( $confs )

echo 1 | mmtbx.quantum_interface iterate_NQH=HIS ${prefix}${conf}.updated.pdb | tee options_${prefix}${conf}.log

end


foreach i ( $opts )
set conf = ""
echo $i | mmtbx.quantum_interface iterate_NQH=HIS ${prefix}${conf}.updated.pdb | tee iterate_${i}.log
end

set phils = `ls -1rt ${prefix}*.phil`

foreach i ( `seq 1 $#phils` )
set phil = $phils[$i]
set conf = `grep altloc $phil | awk '{gsub(/[^A-Za-z]/," ");print}' | awk '$NF~/^[A-Za-z]$/{print $NF}'`
if( "$conf" == "" ) continue
mmtbx.quantum_interface ${prefix}${conf}.updated.pdb iterate_NQH=HIS $phil run_qmr=True qi.nproc=3 $ligcifs |& tee phil${i}.log 
end


rm HIS_votes.txt
foreach i ( `seq 1 $#phils` )
set phil = $phils[$i]
cat phil${i}.log |\
awk '/ resid /{for(i=1;i<=NF;++i)if($i=="resid")resid=$(i+1)}\
 /\!\!\!|></{print $0,resid}' |\
awk '/HD1, HE2/{print $NF,++seen[$NF],"HIP"}\
  /HD1 only/{print $NF,++seen[$NF],"HID"}\
  /HE2 only/{print $NF,++seen[$NF],"HIE"}\
' |\
tee -a HIS_votes.txt

awk '/\!\!\!/{i=$2+0;print file[i]} /phenix.pymol/{gsub(".pml",".pdb");++i;file[i]=$2}' phil${i}.log | tee file.txt
set file = `cat  file.txt`
cp $file frag${i}.pdb

end

awk '! seen[$1,$NF]{string[$1]=string[$1]" "$NF;++seen[$1,$NF]} END{for(r in string)print r,string[r]}' HIS_votes.txt  | sort -g

cat HIS_votes.txt |\
awk '! seen[$1,$NF]{string[$1]=string[$1]" "$NF;++seen[$1,$NF]}\
  END{for(r in string)print r,string[r]}' |\
sort -g |\
awk '{print $2,"A"$1}' |\
tee HIS_settings.txt

