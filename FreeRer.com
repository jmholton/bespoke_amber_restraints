#! /bin/csh -f
#
#   FreeRer.com                                         -James Holton  8-2-04
#
#
#  Automatically generated script for setting up consistent FREE R flags
#
#  Unlike "uniqueify", FreeRer.com can "inherit" free R flags from another file 
#  (mtz or X-plor)
#
#
# set this to wherever your awk program is
set nawk = nawk
$nawk 'BEGIN{print}' >& /dev/null
if($status) set nawk = awk
alias nawk $nawk
#
set tempfile = $CCP4_SCR/FreeRer$$

# defaults
set unfree = ""
set flagfile = ""
set outfile = FreeRed.mtz
set XPLORfile = XPLOR.cv
set FRAC = 0.05
set hires = ""

# process command-line input
foreach arg ( $* )
    # user may specify free-R fraction
    if("$arg" =~ [1-9]*%) set FRAC = `echo "$arg" | nawk '{print ($1 + 0)/100}'`
    if("$arg" =~ 0.*) set FRAC = `echo "$arg" | nawk 'if(($1+0)>0) {print $1 + 0}'`

    # resolution
    if("$arg" =~ [1-9]*A) set hires = `echo "$arg" | nawk '{print $1+0}'`

    # look for MTZ files
    if(("$arg" =~ *.mtz)&&(-e "$arg")) then
	if(! -e "$unfree") then
	    # file to recieve free-R flags is first MTZ encountered
	    set unfree = "$arg"	    
	else
	    # make sure that potential free-R source contains free-R flags
	    mtzdump HKLIN "$arg" << EOF-dump >! ${tempfile}
HEAD
go
EOF-dump
	    grep "FreeR_flag" ${tempfile} >& /dev/null
	    if(! $status) set flagfile = "$arg"
	endif
	rm -f ${tempfile} >& /dev/null
    else
	# free-R source can also be an X-plor hkl file
	if(-e "$arg") then
	    grep "TEST" "$arg" >& /dev/null
	    if(! $status) set flagfile = "$arg"
	endif
    endif
end

if(! -e "$unfree") then
    cat << EOF
usage: $0 mtzfile.mtz [free-R source] [fraction[%]]

Free R flags can be inherited from "free-R source" (mtz or X-plor)
or a new set can be defined as "fraction" of the data in mtzfile.mtz
and the resulting file will be "polished" to fill in missing flags.

EOF
    exit 1
endif

# print out what we are going to do
if(-e "$flagfile") then
    set temp = "${flagfile}'s"
else
    set temp = `echo $FRAC | nawk '{print 100*$1"%"}'`
endif
echo "\n\nadding $temp FreeR_flag to $unfree in output files: $outfile and $XPLORfile"
echo ""
echo ""
# get variables from the input MTZ
mtzdump HKLIN $unfree << EOF-dump >! ${tempfile}
HEAD
go
EOF-dump
if($status) goto bad
if("$hires" == "") then
    set hires = ` nawk '/Resolution Range/{getline; getline; print $6}' ${tempfile} `
endif
set CELL = ` nawk '/Cell Dimensions/ {getline; getline; print}' ${tempfile} `
set SGnum = ` nawk '/Space group/{print $NF+0}' ${tempfile} | tail -1`
set SG = ` nawk -F "[\047]" '/Space group/{print $2}' ${tempfile} `
set SG = ` nawk -v num=$SGnum '$1==num && NF>5{print $4}' ${CLIBD}/symop.lib `

grep FreeR_flag ${tempfile} >& /dev/null
if($status) goto flags_removed
################################################################################
#
# if input file has pre-existing FreeR_flag, I assume you don't want it anymore
# otherwise, you should just use the COMPLETE option in a run of "freerflags"
# or use the same file as the freer source
#
remove_flags:
# purge old Free R flags from $unfree file
mtzutils hklin $unfree hklout ${tempfile}.unfree.mtz << eof-purge
EXCLUDE FreeR_flag
eof-purge
if($status) goto bad
set unfree = ${tempfile}.unfree.mtz
echo "FreeR_flag removed from $1"
flags_removed:

# if flag file has already been given, use it instead
if(-e "$flagfile") then
    mtzdump hklin $flagfile << EOF >&! ${tempfile}.dump2
    HEAD
    go
EOF
    grep FreeR_flag ${tempfile}.dump2  >& /dev/null
    if(! $status) then
	# suggested file is a good mtz, so just use it
	goto add_flags
    else
	grep INDE $flagfile >& /dev/null
	if(! $status) then
#	    echo not done yet...
	    # this is an X-plor file.  We need to import it
	    cat  $flagfile |\
	    nawk '{print toupper($0)}' |\
	    nawk '{for(i=1;i<NF;++i){ \
	           if($i ~ /INDE/){printf "%4d %4d %4d ", $(i+1), $(i+2), $(i+3)}\
	           if($i ~ /TEST/){printf "%4d\n", ($(i+1)+1)%2}\
		   }}' |\
	    cat >! ${tempfile}.import
	    
	    # import the flags into an MTZ
	    f2mtz HKLIN ${tempfile}.import HKLOUT ${tempfile}.mtz << EOF-import
TITLE FREE-R flags Imported from $flagfile
CELL $CELL
SYMM $SG
LABOUT H K L FreeR_flag
CTYPOUT H H H I
EOF-import
	    if($status) goto bad
	    if(! $?debug) rm -f ${tempfile}.import
	    sortmtz HKLOUT ${tempfile}.sorted.mtz << EOF-sort
H K L
${tempfile}.mtz
EOF-sort
	    if($status) goto bad
	    if(! $?debug) rm -f ${tempfile}.mtz
	    # reduce to CCP4 asymmetric unit
	    reindex HKLIN ${tempfile}.sorted.mtz HKLOUT ${tempfile}.imported.mtz << EOF-reindex
EOF-reindex
	    if($status) goto bad
	    set flagfile = ${tempfile}.imported.mtz
	    goto add_flags
	else
	    # can't use this file
	    echo "ERROR: $flagfile is no good!!! "
	endif
    endif
endif


make_flags:
################################################################################
# use unique to generate all possible HKLs
# if you are ever so fortunate to extend your data to 1.5x the RS volume
# this FREE-R set will still be appropriate! 
#
# if you want your free-R to be a different %-age, change FREERFAC
#
set higheres = `echo $hires | awk '{print (1.5*$1^(-3))^(-1./3)}'`
unique HKLOUT ${tempfile}.unique.mtz << EOF-unique
TITLE  Unique data for $SG
LABOUT  F=FP SIGF=SIGFP
# get RESOL, SYMM and CELL
RESO $higheres
SYMM $SGnum
CELL $CELL
EOF-unique
if($status) goto bad
#
freerflag HKLIN ${tempfile}.unique.mtz HKLOUT freeR_flag.mtz << EOF-FreeR
FREERFAC $FRAC
END
EOF-FreeR
if($status) goto bad
set flagfile = freeR_flag.mtz

# the file freeRflags.mtz now contains the Free-R flags


add_flags:
################################################################################
#
#	Combine these FREE R flags with the merged data
mtzutils        \
HKLIN1 $unfree \
HKLIN2 $flagfile \
HKLOUT ${tempfile}.freed.mtz \
<< EOF-lastcad

RESOLUTION 1000 $hires

INCLUDE 1 ALL
INCLUDE 2 FreeR_flag

END
EOF-lastcad
if($status) goto bad
if($status) exit 9

# polish flags here (that is, fill in holes)
freerflag HKLIN ${tempfile}.freed.mtz HKLOUT $outfile << EOF-polish
COMPLETE FREE=FreeR_flag
END
EOF-polish
if($status) goto bad
echo ""
echo ""
echo "$outfile is now identical to $1 with Free-R flags added" `if(-e $flagfile) echo "from $flagfile"`

################################################################################
#
#	list them in X-PLOR format as well
#

# output new (polished) flags from output file
sortmtz HKLOUT ${tempfile}.sorted.mtz << EOF-sort
L K H
$outfile
EOF-sort
if($status) goto bad
mtz2various hklin ${tempfile}.sorted.mtz hklout ${tempfile}.cv << EOF
OUTPUT XPLOR
MISS 0.0
LABIN FP=FreeR_flag SIGFP=FreeR_flag FREE=FreeR_flag
END
EOF
if($status) goto bad
nawk '/TEST= 1/{print substr($0, 6, 13)}' ${tempfile}.cv >! ${tempfile}_1.cv

# transform X-plor output into something X-plor will read without problems
set NREF = `wc ${tempfile}.cv | nawk '{print $1 +1}'`

echo " NREFlections= $NREF"                                         >! $XPLORfile
echo " ANOMalous=FALSE"                                             >> $XPLORfile
echo " DECLare NAME=TEST   DOMAin=RECIprocal   TYPE=INTE END"       >> $XPLORfile
nawk '{print substr($0, 1, 18) substr($0, 47) }' ${tempfile}.cv     >> $XPLORfile
if($status) goto bad



#########################################
# output flags from original flag file
sortmtz HKLOUT ${tempfile}.sorted.mtz << EOF-sort
L K H
$flagfile
EOF-sort
if($status) goto bad
mtz2various hklin ${tempfile}.sorted.mtz hklout ${tempfile}.cv << EOF
OUTPUT XPLOR
MISS 0.0
LABIN FP=FreeR_flag SIGFP=FreeR_flag FREE=FreeR_flag
END
EOF
if($status) goto bad
nawk '/TEST= 1/{print substr($0, 6, 13)}' ${tempfile}.cv >! ${tempfile}_0.cv

# find differences between new and old free R flags
diff ${tempfile}_0.cv ${tempfile}_1.cv | grep '> ' >! ${tempfile}.diff


echo ""
echo ""
echo `cat ${tempfile}.diff | wc -l`" new Free HKLs assigned during polishing:"
if($temp == 0) then
    echo none
else
    cat ${tempfile}.diff
endif
echo "above "`cat ${tempfile}.diff | wc -l`" new Free HKLs assigned during polishing."
echo ""
echo ""


set temp = `echo $FRAC | nawk '{print 100*$0}'`
if(-e "${tempfile}.unique.mtz") echo "$flagfile was created with ${temp}% free-R flags."
if("$flagfile" == "${tempfile}.imported.mtz") set flagfile = $2
echo "$outfile is now identical to $1 with Free-R flags added from the file $flagfile"
echo "$XPLORfile contains these Free-R flags in X-PLOR format"


if($?debug) exit
rm -f ${tempfile}* >& /dev/null


exit


bad:
echo "ACK! \07"
exit 9

#########################
# future improvements
#
- use a target of 1000 free-R reflections, instead of a %age

