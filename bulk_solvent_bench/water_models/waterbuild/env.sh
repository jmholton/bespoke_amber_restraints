# sourced by queue.sh and check_nohome.sh: an environment that names only /scratch and the local system.  CCP4_SCR
# is the original /tmp/jamesh (local disk): CCP4 programs given longer temporary file names behave differently
P=/scratch/jamesh/programs
W=/scratch/jamesh/projects/amber/waterbuild
export HOME=$W/home OPENSSL_CONF=$P/phenix-2.2-6143/ssl/openssl.cnf
export PATH=$W/bin:$P/ccp4-9/bin:/usr/local/bin:/usr/bin:/bin
export CCP4_MASTER=$P CCP4=$P/ccp4-9 CBIN=$P/ccp4-9/bin CLIB=$P/ccp4-9/lib CLIBD=$P/ccp4-9/lib/data
export CINCL=$P/ccp4-9/include CETC=$P/ccp4-9/etc CCP4_SCR=/tmp/jamesh CCP4_OPEN=UNKNOWN
export PHENIX_OVERWRITE_ALL=true
mkdir -p $CCP4_SCR
unset TMPDIR PHENIX PHENIX_PREFIX PHENIX_VERSION PHENIX_ROSETTA_PATH PYTHONPATH PYTHONHOME LD_LIBRARY_PATH MMCIFDIC CLIBD_MON
unset CCP4I_TOP CCP4I_TCLTK CRANK CHTML CEXAM CCP4_HELPDIR AMBERHOME
