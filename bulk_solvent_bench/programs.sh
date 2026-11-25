# programs.sh: bash version of programs.sourceme (same tree, same phenix link dir), for interactive use and bash
# scripts: source /home2/jamesh/projects/amber/bulk_solvent/programs.sh   (the MD engine settings are tcsh only)
export PROGRAMS=${PROGRAMS:-/data/home/programs}
programs_bench=/home2/jamesh/projects/amber/bulk_solvent
source $PROGRAMS/amber22/amber.sh
export PHENIX=$PROGRAMS/phenix-2.2-6143 PHENIX_PREFIX=$PROGRAMS/phenix-2.2-6143 PHENIX_VERSION=2.2-6143
programs_pbin=$programs_bench/programs_bin/$(echo $PROGRAMS | sed 's#^/##; s#/#_#g')/phenix-2.2-6143
if [ ! -d $programs_pbin ]; then
  mkdir -p $programs_pbin
  for f in $PHENIX/phenix_bin/*; do ln -sf $PHENIX/bin/${f##*/} $programs_pbin/${f##*/}; done
fi
export PATH=$programs_pbin:$PATH
source <(sed -e 's#"/programs"#"'$PROGRAMS'"#' -e "s#'/programs/#'$PROGRAMS/#" $PROGRAMS/ccp4-9/bin/ccp4.setup-sh)
