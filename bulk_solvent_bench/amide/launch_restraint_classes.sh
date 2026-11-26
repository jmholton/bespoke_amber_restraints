#!/bin/bash
# Restraint-class tests (2026-11-24; bulk_solvent/MD_BULK_SOLVENT_JOINT.md, "Open question"): one build
# (build_2s2r_amide), code_391d233, uniform B 5, 200 itr, min_align_weight 0 in every arm (so the alignment floor
# cannot bring back a restraint class an arm leaves out).  Each run: seed (maxitr 0), edit the start list
# (restraints_for_0.pdb + current_restraints.pdb), then the wrapper to itr 200.
#   rc_dw_w  density-weighted protein heavy atoms + density-picked waters        b5sd1 (reference)
#   rc_dw    the same, no waters: start list without HOH, waters never re-picked  b5sd1 + repick_hohscale=1000
#   rc_bb_w  backbone N CA C O + waters (+ salt), density-weighted                b5sd1 + repick_itr=0 release_itr=0
#   rc_ca1   CA only at fixed k 1, no waters                                      fixed
#   rc_k10   every restrained atom at fixed k 10 (9.9999)                         fixed
# Run on octamus1:  ssh octamus1 'bash /scratch/jamesh/projects/amber/amide/launch_restraint_classes.sh'
S=/scratch/jamesh/projects/amber; H=$S/1aho_hurry_test; code=code_391d233
runs="rc_dw_w:b5sd1:none: rc_dw:b5sd1:nohoh:repick_hohscale=1000 rc_bb_w:b5sd1:backbone:repick_itr=0,release_itr=0
      rc_ca1:fixed:ca1: rc_k10:fixed:k10:"
[ -n "$1" ] && runs="$*"
for r in $runs; do
  IFS=: read tag variant edit extra <<< "$r"
  name=ff19_opc_$tag; d=$H/$name
  extra="min_align_weight=0 ${extra//,/ }"
  if pgrep -f "run_local.csh $name " >/dev/null; then echo "$name wrapper running, skipping"; continue; fi
  if [ ! -e $d/amber_0.rst7 ]; then
    ( cd $H; BUILD=build_2s2r_amide /bin/tcsh -f $S/run_local.csh $name 0 $variant $code </dev/null >> $S/$name.log 2>&1 )
    [ -e $d/amber_0.rst7 ] || { echo "$name: seeding failed, see $S/$name.log"; continue; }
    for f in restraints_for_0.pdb current_restraints.pdb; do
      awk -v e=$edit '/^ATOM|^HETAT/{n=substr($0,13,4);gsub(/ /,"",n);rn=substr($0,18,3);\
            w=(rn=="HOH"||rn=="WAT"); s=(rn=="ACY"||rn=="NH4"); bb=(n=="N"||n=="CA"||n=="C"||n=="O");\
            if(e=="nohoh" && w) next;\
            if(e=="backbone" && !w && !s && !bb) next;\
            if(e=="ca1"){ if(n!="CA" || w || s) next; $0=substr($0,1,60) "100.00" substr($0,67)}\
            if(e=="k10") $0=substr($0,1,60) "999.99" substr($0,67);\
            ++c[(w?"water":(s?"salt":"protein"))]}\
          {print} END{printf("protein %d water %d salt %d\n",c["protein"],c["water"],c["salt"]) > "/dev/stderr"}' \
        $d/$f > $d/$f.rc 2> $d/$f.count
      mv $d/$f.rc $d/$f
    done
    echo "$tag: $variant, start list $edit ($(cat $d/restraints_for_0.pdb.count)), extra: $extra" > $d/restraint_class.txt
    rm -f $d/*.count
  fi
  ( cd $H; BUILD=build_2s2r_amide setsid nohup /bin/tcsh -f $S/run_local.csh $name 200 $variant $code $extra \
      </dev/null >> $S/$name.log 2>&1 & )
  echo "launched $name: $(cat $d/restraint_class.txt)"
  sleep 20
done
