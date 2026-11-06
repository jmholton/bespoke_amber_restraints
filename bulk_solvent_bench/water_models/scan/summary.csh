#! /bin/tcsh -f
# one line per scan run: restraints, waters, pressures, R-free of both refinement protocols, difference-map blobs
set W = /home2/jamesh/projects/amber/bulk_solvent/water_models
cd $W/scan
printf "%-8s %-38s %6s %8s %8s %8s | %-15s | %-15s | %s\n" tag settings nwater P_hyd P_prod P_unrst "example Rfree" "deposited Rfree" "diff rms solv, +core score/+total"
foreach tag ( `awk '{print $1}' settings.txt` opc_r2 )
  set set = "`awk -v t=$tag '"'$1==t{$1="";$2="";print}'"' settings.txt`"
  set nw = "-" ; set ph = "-" ; set pp = "-"
  if( -e $tag/pressure_vs_seg.txt ) then
    set nw = `awk 'NF==5{n=$5} END{print n}' $tag/pressure_vs_seg.txt`
    set ph = `awk 'NF==5{p=$2} END{print p}' $tag/pressure_vs_seg.txt`
    set pp = `awk '/^prod_/{s+=$2;n++} END{if(n) printf("%.0f",s/n); else print "-"}' $tag/pressure_vs_seg.txt`
  endif
  set pu = "-"
  if( -e $tag/barometer.log ) set pu = `awk '/^barometer/{print $(NF-4)}' $tag/barometer.log`
  if( "$pu" == "" ) set pu = "-"
  set ex = `cat $W/example_check/exocc_starthere_scan_${tag}_j?.out |& awk '/Final R-work/{gsub(",","");s+=$NF;n++} END{if(n) printf("%.4f n=%d",s/n,n)}'`
  set dp = `cat $W/phenix/jig_scan_${tag}_j?.out |& awk '/Final R-work/{gsub(",","");s+=$NF;n++} END{if(n) printf("%.4f n=%d",s/n,n)}'`
  set db = "-"
  if( -e diffblobs_${tag}.txt ) set db = `awk '/map rms/{for(i=1;i<NF;i++)if($i=="region"){r=$(i+1);gsub(",","",r);break}} /positive blobs/ && ! p{++p;gsub(/[,;:]/," ");for(i=1;i<=NF;i++){if($i=="total")t=$(i+2);if($i=="(<1")c=$(i+4)}} END{print r, c "/" t}' diffblobs_${tag}.txt`
  printf "%-8s %-38s %6s %8s %8s %8s | %-15s | %-15s | %s %s\n" $tag "$set" $nw $ph $pp $pu "$ex" "$dp" "$db" ""
end
