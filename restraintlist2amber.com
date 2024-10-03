#! /bin/tcsh -f
#
#   Convert a list of desired restraint weights per atom to AMBER group restraint format      -James Holton 11-5-24
#
#   format:  atomname ordinal_resnum any other things  weight
#
#   can also extract B factors from a second PDB file as the weights
#

set tempfile = /dev/shm/${USER}/temp_r2a$$_
mkdir -p /dev/shm/${USER}

set logfile = debuglog.log


set infile = msander.in
set suffixfile = restraint_suffix.in
set outfile = new.in
set restraint_list = all_restraint_weights.txt
set out_restraint_list = current_restraints.txt
set orignames  = ""
set scale = 1
set pdbscale = 0.01
set minwt = 1e-6
set allatom_weight = 0
set binthresh = 1.1


echo "command-line arguments: $* "

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print tolower($1)}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set key = `echo $Key | awk '{print tolower($0)}'`
    set val = `echo $Val | awk '{print tolower($0)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`

    if( $assign ) then
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "infile") set infile = "$Val"
      if("$key" == "list") set restraint_list = "$Val"
      if("$key" == "scale") set scale = "$Val"
      if("$key" == "restraint_minwt") set minwt = "$Val"
      if("$key" == "thresh") set binthresh = "$Val"
    else
      # no equal sign
      if("$key" =~ *.in ) set infile = "$Arg"
      if("$key" =~ *.txt ) set restraint_list = "$Arg"
    endif
    if("$arg" == "debug") set debug = 1
end

if( $?debug ) then
    set tempfile = tempfile
endif
touch $logfile

cat << EOF
list       $restraint_list
infile     $infile
binthresh  $binthresh
scale      $scale
minwt      $minwt
allatom    $allatom_weight
outfile    $outfile
outlist    $out_restraint_list
suffixfile $suffixfile

tempfile   $tempfile
debug      $?debug
EOF

set t = "$tempfile"

# if we already have master list, just generate amber restraints
if(-e "$restraint_list") goto pdbextract

echo "extracting restraints from $infile"
cat $infile |\
awk '/Specific/{atom=$2;getline;w=$1}\
  /^RES/{for(n=$2;n<=$3;++n){print atom,n,w}}' |\
sort -u | sort -k3gr >! "$out_restraint_list"

set restraint_list = "$out_restraint_list"


pdbextract:

if("$restraint_list" =~ *.pdb) then

  if(! -e "$orignames") then
    set BAD = "need orignames.pdb to interpret $restraint_list "
    goto exit
  endif
  
  # append ordinal residue number at end of each line ( the ordinal residue that amber uses )
  echo "counting residues in $orignames"
  awk '/^ATOM|HETAT/{print substr($0,1,80)}' $orignames |\
  convert_pdb.awk -v append=ordresnum -v skip=H |\
  cat >! ${t}ordresnum.pdb

  # combine partial-list B factors with ordinal residue numbers as new partial list
  # extracts only atom IDs named in refpointspdb
  echo "extracting atoms named in $restraint_list"
  combine_pdbs_runme.com saveBfac=1 ${t}ordresnum.pdb \
     $restraint_list outfile=${t}restme.pdb >> $logfile

  echo "using $pdbscale x B factors in $restraint_list as restraint weights"
  # use B factors in reference points PDB as restraint weights
  cat ${t}restme.pdb |\
  awk -v pdbscale=$pdbscale '/^ATOM|^HETAT/{\
    print substr($0,12,5),$NF,substr($0,18,5),substr($0,23,8),pdbscale*substr($0,61,6)}' |\
  cat >! ${out_restraint_list}
  # atom ordresnum typ chain resnum  weight

  touch ${t}allother.txt
  if ( "$allatom_weight" != "0" ) then
    echo "extracting every non-H atom not mentioned in restraint list"
    combine_pdbs_runme.com xor=1 $restraint_list ${t}ordresnum.pdb \
      outfile=${t}unrestrained.pdb >> $logfile
    cat ${t}unrestrained.pdb |\
    awk '/^ATOM|^HETAT/{\
      print substr($0,12,5),$NF,substr($0,18,5),substr($0,23,8),"w"}' |\
    cat >! ${t}allother.txt
    set n = `cat ${t}allother.txt | wc -l`
    echo "found $n"
  endif

  # no longer a pdb file
  set restraint_list = ${out_restraint_list}
endif



generate:

# allow variable number of columns, last one is the weight
set NF = `awk '{print NF;exit}' $restraint_list`

echo $scale $minwt |\
cat - $restraint_list |\
awk 'NR==1{scale=$1;minwt=$2;next}\
  {$NF*=scale}\
  $NF>=minwt{print $1,$2,$NF}' |\
sort -k3gr -k1,1 >! ${t}anw.txt

if ( "$allatom_weight" != "0" ) then
    # apply it to all atoms not otherwise restrained
    echo $scale $minwt |\
    cat - $restraint_list |\
    awk 'NR==1{scale=$1;minwt=$2;next}\
      {$NF*=scale}\
      $NF<minwt{print $1,$2}' |\
    cat - ${t}allother.txt |\
    awk -v w=$allatom_weight '{print $1,$2,w}' |\
    cat >> ${t}anw.txt
endif

# assign all weight ranges to bins
sort -k3gr -k1,1 ${t}anw.txt |\
awk -v bt=$binthresh 'NR==1{print $NF;hi=$NF}\
  hi/$NF>bt{print $NF;hi=$NF}\
 END{print $NF}' |\
sort -u | sort -gr |\
awk '{printf("BREAK %.12g\n",$1-0.000001)}' |\
cat - ${t}anw.txt |\
awk '{print $NF,$0}' |\
sort -gr |\
awk '/BREAK/{++bin;next} {$1="";print $0,bin+1}' |\
cat >! bin_assignments.txt
# atom oresnum weight bin

set bins = `awk '{print $NF}' bin_assignments.txt | sort -u | sort -g`
echo "$#bins bins created"

if("$suffixfile" == "" ) set suffixfile = new_restraints.in

echo "writing restraints to $suffixfile"
echo -n "" >! "$suffixfile"
# restrain high-weight atoms first
foreach bin ( $bins )

  # select which atoms we will restrain
  awk -v bin=$bin '$NF==bin{$NF="";print}' bin_assignments.txt >! selected.txt
  set median_weight = `awk '{print $NF}' selected.txt | median.awk | awk '{print $1}'`
  set weight = $median_weight

  # check for crazy weights
  set test = `echo $weight | awk '{print ( $1 <= 0 )}'`
  if( $test ) then
    echo "WARNIING: non-positive weight $weight in $restraint_list "
  endif

  set nselected = `cat selected.txt | wc -l`

  set atomtypes = `awk '{print $1}' selected.txt | sort -u`

  # loop over atom types
  foreach atom ( $atomtypes )

    awk -v atom=$atom '$1==atom' selected.txt |\
    cat >! ${t}atoms.txt
    set natoms = `cat ${t}atoms.txt | wc -l`
    set wrange = `awk '{print $NF}' ${t}atoms.txt | sort -g | awk 'NR==1{min=$NF} END{print min"-"$NF}'`

    set median_weight = `awk '{print $NF}' ${t}atoms.txt | median.awk | awk '{print $1}'`
    set weight = $median_weight

    cat << EOF >> "$suffixfile"
     Specific ${atom} restraints on $natoms atoms in range: $wrange
     $weight
FIND
$atom * * *
SEARCH
EOF
    # list as ranges if they are there - more compact
    awk '{print $2}' ${t}atoms.txt |\
    sort -g |\
    awk 'NR==1{f=$1;ex=f+1;next}\
           $1!=ex{print "RES",f,ex-1;f=$1}\
           {ex=$1+1} \
           END{print "RES",f,$1}' |\
    cat >> "$suffixfile"

    # end of RES lines for this atom
    echo "END" >> "$suffixfile"
  end
end
# end of all restraints
echo "END" >> "$suffixfile"

echo "amber restraints written to $suffixfile"

if( -e "$infile" && $outfile != "" ) then
#     / ntr=1,/ && allatom_weight!=0{\
# using restraintmask will cancel all group restraints
#        print "  restraintmask=\047! @H= & ! @EP=\047,";\
#        print "  restraint_wt="w",";\
#     }\
    echo "and updated in $outfile "
    echo $allatom_weight |\
    cat - $infile |\
    awk 'NR==1{allatom_weight=$1;next}\
     /Specific|^END/{exit} \
     / ntr=0,/{gsub("ntr=0","ntr=1")}\
     /restraint/{next}\
     {print}' >! ${t}.in
    cat $suffixfile >> ${t}.in
    mv ${t}.in $outfile
endif

exit:
if($?BAD) then
   echo "ERROR: $BAD"
   exit 9
endif

# clean up
if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $?debug && ! ( "$tempdir" == "." || "$tempdir" == "" ) ) then
    rm -f ${tempfile}*
endif


exit

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################




