#! /bin/tcsh -f
#
#  delete the per-iteration debugging files optimize_weights_runme.com leaves behind
#                                                                -James Holton 11-5-26
#
#  For every iteration optimize_weights keeps copies of the B-factor field and the
#  per-atom B, the restraint and water-remap bookkeeping PDBs, and the restart files
#  of each MD sub-stage.  They are how you debug a run, but production does not need
#  them, and they add up: Bfac_N.map alone is ~2 MB per iteration for 1aho and ~160 MB
#  for 6c2r 2x2x2.
#
#  Always kept: the last keep= iterations whole (a restart, or the run itself if it
#  is still going, reads them), every keepevery=-th iteration whole (places to hand
#  a stage on from), and in every iteration the products: amber_N.nc/.rst7/.in/.out,
#  barometer_N.out, avg_N.mtz and all the logs.
#
#  usage: cleanup_optimize_runme.com [rundir] [keep=2] [keepevery=0] [dryrun=1]
#
set rundir = .
# keep the last this many iterations whole
set keep = 2
# also keep every iteration divisible by this whole (0: none)
set keepevery = 0
# 1: only list what would go
set dryrun = 0
# what goes (1) or stays (0)
#  maps       Bfac_N.map premod_Bfac_N.map
#  Bfiles     Bfac_N.pdb and the pre-step copies of it (premod_, Bfac_preadjust_, presmooth_, predry_ ...),
#             sorted_Bmods_N.txt
#  waters     preremap_N/remapped_N/teleported_waters_N/new_positions_N.pdb, remapped_pairs_N.txt, wraps_N.txt
#  restraints restraints_for_N.pdb and its pre-step copies, actual_restraints_N.pdb, full_refpoints_N.pdb,
#             sorted_mults_N.txt, restraints_challenge_vs_weight_N.txt
#  restarts   start_N.rst7/.pdb, amber_N_equi*.rst7, amber_N_settle*.rst7, amber_N_unwrapped.rst7
set maps = 1
set Bfiles = 1
set waters = 1
set restraints = 1
set restarts = 1

foreach Arg ( $* )
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    if( "$Arg" !~ *=* ) then
      if( -d "$Arg" ) then
        set rundir = "$Arg"
        continue
      endif
      set BAD = "unrecognized argument: $Arg"
      goto exit
    endif
    set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
    if( ! $test ) then
      set BAD = "unrecognized option: $Key"
      goto exit
    endif
    set $Key = "$Val"
end

cd $rundir
if( $status ) then
    set BAD = "cannot cd to $rundir"
    goto exit
endif

# the last finished iteration, as optimize_weights finds it on a restart
set last = `ls -1 | awk -F "[_.]" '/^amber_[0-9]+\.rst7$/{n=$2+0; if(n>m)m=n} END{print m+0}'`
if( "$last" == "0" ) then
    set BAD = "no amber_N.rst7 in ${rundir}: not an optimize_weights run directory?"
    goto exit
endif
set lastdel = `echo $last $keep | awk '{print $1-$2}'`
set everymsg = ""
if( $keepevery > 0 ) set everymsg = "and every ${keepevery}th"
echo "last iteration ${last}: keeping iterations after $lastdel whole $everymsg"

ls -1 |\
awk -v last=$last -v keep=$keep -v every=$keepevery -v maps=$maps -v Bfiles=$Bfiles \
    -v waters=$waters -v restraints=$restraints -v restarts=$restarts '\
  BEGIN{N="_([0-9]+)";\
    if(maps)      re["^(premod_)?Bfac" N "\\.map$"]=1;\
    if(Bfiles)  { re["^(premod_|prepad_|prephenix_|prerefmac_|pretrunc_|predry_)?Bfac" N "\\.pdb$"]=1;\
                  re["^(Bfac_preadjust|presmooth_Bfac_for)" N "\\.pdb$"]=1;\
                  re["^sorted_Bmods" N "\\.txt$"]=1}\
    if(waters)  { re["^(preremap|remapped|teleported_waters|old_teleported_waters|new_positions)" N "\\.pdb$"]=1;\
                  re["^(remapped_pairs|wraps)" N "\\.txt$"]=1}\
    if(restraints){re["^(restraints_for|prerelease_restraints_for|predel_restraints_for|predelalts_restraints_for)" N "\\.pdb$"]=1;\
                  re["^(actual_restraints|randelled_restraints|full_refpoints)" N "\\.pdb$"]=1;\
                  re["^(sorted_mults|restraints_challenge_vs_weight)" N "\\.txt$"]=1}\
    if(restarts){ re["^start" N "\\.(rst7|pdb)$"]=1;\
                  re["^amber" N "_(equi|settle)(_[0-9]+)?\\.rst7$"]=1;\
                  re["^amber" N "_unwrapped\\.rst7$"]=1}}\
  {for(r in re) if(match($0,r)){s=$0;sub(/^[^0-9]*_/,"",s);n=s+0;\
      if(n>last-keep) break; if(every>0 && n%every==0) break; print; break}}' |\
cat >! cleanup_list.txt

set n = `cat cleanup_list.txt | wc -l`
set MB = `cat cleanup_list.txt | xargs -r stat -c %s | awk '{s+=$1} END{printf("%.0f",s/1e6)}'`
if( $dryrun ) then
    echo "would delete $n files ($MB MB), listed in $cwd/cleanup_list.txt"
    exit
endif
cat cleanup_list.txt | xargs -r rm -f
echo "deleted $n files ($MB MB)"
rm -f cleanup_list.txt

exit:
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif
exit
