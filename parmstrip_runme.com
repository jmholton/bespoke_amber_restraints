#! /bin/tcsh -f
#
#  strip atoms out of an Amber topology without losing its CMAP terms
#
#  cpptraj before V7.00.0 silently drops all CMAP terms whenever it strips a topology that did not come from
#  CHARMM (chamber).  For ff19SB that is the entire backbone phi/psi energy: frcmod.ff19SB zeroes the phi and psi
#  torsions and puts them in the CMAP grids.  So: use a cpptraj that is new enough (cpptraj_strip, else the one on
#  the path), fall back to ParmEd if there is none, and refuse to write a topology that lost its CMAP terms.
#
#  usage: parmstrip_runme.com padded.parm7 mask=:816-2115 outfile=xtal.prmtop
#
set parmfile = ""
set mask = ""
set outfile = ""
# cpptraj to use for strips if it is V7.00.0 or later (e.g. a standalone build); "" = only the cpptraj on the path.
# Usually set in compute_settings.sourceme
set cpptraj_strip = ""
set debug = 0

set tempfile = /dev/shm/${USER}/tempfile_parmstrip_$$_

foreach sourceme ( compute_settings.sourceme )
   if(-e $sourceme ) then
      source $sourceme
   endif
end

foreach Arg ( $* )
    set assign = `echo $Arg | awk '{print ( /=/ )}'`
    set key = `echo $Arg | awk -F "=" '{print tolower($1)}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    if( $assign ) then
      if("$key" == "mask") set mask = "$Val"
      if("$key" == "outfile") set outfile = "$Val"
      if("$key" == "parmfile") set parmfile = "$Val"
      if("$key" == "cpptraj_strip") set cpptraj_strip = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
      if("$key" == "debug") set debug = "$Val"
    else
      if("$parmfile" == "") set parmfile = "$Arg"
    endif
end

if(! -e "$parmfile") then
  set BAD = "cannot read topology $parmfile"
  goto exit
endif
if("$mask" == "" || "$outfile" == "") then
  set BAD = "usage: $0 topology.parm7 mask=<amber mask> outfile=<new topology>"
  goto exit
endif

set t = "$tempfile"
mkdir -p `dirname $t` >& /dev/null

# CMAP terms in the original
set cmap_in = `awk '/^%FLAG (CHARMM_)?CMAP_COUNT/{getline;getline;print $1+0;exit}' $parmfile`
if("$cmap_in" == "") set cmap_in = 0

# first cpptraj that is new enough: cpptraj_strip, then the one on the path
set tool = ""
set oldcpptraj = ""
foreach exe ( $cpptraj_strip cpptraj )
  if("$exe" != "cpptraj" && ! -x "$exe") continue
  set version = `$exe --version |& awk '{for(i=1;i<=NF;++i) if($i~/^V[0-9]/){print substr($i,2);exit}}'`
  if("$version" == "") continue
  set new = `echo $version | awk -F "." '{print ( $1+0 >= 7 )}'`
  if( $new ) then
    set tool = "$exe"
    set toolname = "cpptraj V$version"
    break
  endif
  if("$oldcpptraj" == "") set oldcpptraj = "$exe V$version"
end
if("$tool" == "" && "$oldcpptraj" != "" && $cmap_in == 0) then
  # no CMAP terms to lose: any cpptraj will do
  set tool = `echo $oldcpptraj | awk '{print $1}'`
  set toolname = `echo $oldcpptraj | awk '{print "cpptraj",$2}'`
endif
if("$tool" == "") then
  set parmed = `which parmed |& egrep -v "not found|no parmed" | head -n 1`
  if("$parmed" == "") then
    set BAD = "no cpptraj V7.00.0 or later and no ParmEd: cannot strip $parmfile without losing its $cmap_in CMAP terms"
    goto exit
  endif
  if("$oldcpptraj" != "") then
    echo "WARNING: $oldcpptraj drops CMAP terms when it strips a topology (fixed in cpptraj V7.00.0)."
    echo "WARNING: stripping with ParmEd instead (slower).  Set cpptraj_strip in compute_settings.sourceme"
    echo "WARNING: to a cpptraj V7.00.0 or later to make this fast."
  endif
  set toolname = "ParmEd"
endif

# write to a temporary file first: the output may replace the input
rm -f ${t}out.parm7
if("$toolname" == "ParmEd") then
  $parmed -O -p $parmfile << EOF >&! ${t}strip.log
strip $mask
outparm ${t}out.parm7
quit
EOF
  set err = $status
else
  $tool -p $parmfile << EOF >&! ${t}strip.log
parmstrip $mask
parmwrite out ${t}out.parm7
EOF
  set err = $status
endif
if( $err || ! -s ${t}out.parm7 ) then
  tail -n 5 ${t}strip.log
  set BAD = "$toolname failed to strip $mask from $parmfile"
  goto exit
endif

set natom_in  = `awk '/^%FLAG POINTERS/{getline;getline;print $1;exit}' $parmfile`
set natom_out = `awk '/^%FLAG POINTERS/{getline;getline;print $1;exit}' ${t}out.parm7`
set cmap_out = `awk '/^%FLAG (CHARMM_)?CMAP_COUNT/{getline;getline;print $1+0;exit}' ${t}out.parm7`
if("$cmap_out" == "") set cmap_out = 0
if( $cmap_in > 0 && $cmap_out == 0 ) then
  set BAD = "$toolname dropped all $cmap_in CMAP terms of $parmfile"
  goto exit
endif

mv ${t}out.parm7 $outfile
echo "parmstrip: $parmfile -> $outfile by ${toolname}: $natom_in -> $natom_out atoms, CMAP terms $cmap_in -> $cmap_out"

exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   if(! $debug) rm -f ${t}* >& /dev/null
   exit 9
endif

if( $debug ) exit
rm -f ${t}* >& /dev/null

exit
