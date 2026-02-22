#! /bin/tcsh -f
#
# make a bond edit list that disables nonbonds that were not in the original structure
#
#
set pdbfile = ""
set ciffiles = ""
set mtzfile = refme.mtz

set outfile = phenix_opts_unbump.eff

set tempfile = tempfile_nnnb_$$_
set debug = 0

if(-e xtal_properties.sourceme) source xtal_properties.sourceme

# read the command line to update variables and other settings
foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if( $assign ) then
      # re-set any existing variables
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
      if("$Arg" =~ *.mtz ) set mtzfile = $Arg
      if("$Arg" =~ *.cif ) set ciffiles = ( $ciffiles $Arg )
    endif
    if("$key" == "debug") set debug = "1"
end

# shorthand for temporary file
set t = $tempfile


if(! -e opts.eff) then
  echo "WARNING: no opts.eff provided. making one"
cat << EOF >! opts.eff
refinement {
  refine {
    strategy = *individual_sites individual_sites_real_space rigid_body \
               *individual_adp group_adp tls *occupancies group_anomalous den
    occupancies {
      individual = water
    }
  }
  pdb_interpretation {
    restraints_library {
      cdl = False
      mcl = False
    }
    flip_symmetric_amino_acids = False
    correct_hydrogens = False
    allow_polymer_cross_special_position = True
    automatic_linking {
      link_none = True
    }
    exclude_from_automatic_linking {
      selection_1 = All
      selection_2 = All
    }
  }
  bulk_solvent_and_scale {
    apply_back_trace=False
  }
  main {
    nqh_flips=False
    max_number_of_iterations=100
  }
}
EOF
endif

if( "$pdbfile" != "${t}nnnb.pdb" ) cp $pdbfile ${t}nnnb.pdb

if( 0 ) then
cat opts.eff |\
awk '$NF=="\\"{$NF="";printf("%s",$0);getline}\
    {print}' |\
awk '$NF=="{"{++ind} $1=="}"{--ind}\
    $1=="strategy"{next}\
    /nqh_flip/{next}\
    $1=="bulk_solvent_and_scale"{skip=ind}\
    /occupanc/{skip=ind}\
   # {print "DEBUG",$0,ind,skip}\
    skip-1>ind{skip=0}\
    skip {next}\
    {print $0}' |\
cat >! geo_opts.eff

echo "geometry pre-run"
rm -f nnnb.geo
phenix.geometry_minimization ${t}nnnb.pdb $ciffiles macro_cycles=0 \
  stop_for_unknowns=false geo_opts.eff \
  allow_polymer_cross_special_position=True \
  output_file_name_prefix=nnnb >! geom_nnnb.log
if( $status || ! -e nnnb.geo ) then
  tail geom_nnnb.log
  set BAD = "initial geometry_minimization failed"
  goto exit
endif
endif

echo "zero-cycle pre-run"
phenix.refine ${t}nnnb.pdb $ciffiles $mtzfile main.number_of_macro_cycles=0 \
  opts.eff \
  allow_polymer_cross_special_position=True \
  prefix=nnnbp >! pr_nnnb.log

cp nnnbp_001.geo nnnb.geo

if(! -e start.geo ) then
  echo "WARNING: using geometry for start.geo"
  cp nnnb.geo start.geo
endif

foreach geo ( start nnnb )
cat ${geo}.geo |\
awk '/nonbonded pdb=/{key="NONBOND";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;getline;\
     obs=$1;ideal=$2;sym=$3;\
     a1=substr(id1,1,4);a2=substr(id2,1,4);\
     f1=substr(id1,5,1);f2=substr(id2,5,1);\
     t1=substr(id1,6,4);t2=substr(id2,6,4);\
     c1=substr(id1,10,1);c2=substr(id2,10,1);\
     r1=substr(id1,11,5);r2=substr(id2,11,5);\
     gsub(" ","_",f1);gsub(" ","_",f2);\
     gsub(" ","_",c1);gsub(" ","_",c2);\
     print a1,f1,t1,c1,r1,"-",a2,f2,t2,c2,r2,"|",obs,sym}' |\
cat >! nonbonds_${geo}.txt
end

# make sure we dont duplicate bonds
cat ${geo}.geo |\
awk '/bond pdb=/{key="BOND";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;getline;\
     obs=$2;ideal=$1;\
     a1=substr(id1,1,4);a2=substr(id2,1,4);\
     f1=substr(id1,5,1);f2=substr(id2,5,1);\
     t1=substr(id1,6,4);t2=substr(id2,6,4);\
     c1=substr(id1,10,1);c2=substr(id2,10,1);\
     r1=substr(id1,11,5);r2=substr(id2,11,5);\
     gsub(" ","_",f1);gsub(" ","_",f2);\
     gsub(" ","_",c1);gsub(" ","_",c2);\
     print a1,f1,t1,c1,r1,"-",a2,f2,t2,c2,r2,"|",obs,ideal,key}' |\
cat >! bonds_${geo}.txt

awk '{print $0,"PREV"}' nonbonds_start.txt |\
 cat -  bonds_nnnb.txt nonbonds_nnnb.txt |\
 awk -F "|" '/PREV/{++seen[$1];next} \
 /BOND/{++bonded[$1];++seen[$1];next} \
 ! seen[$1] || / HOH /{print}' |\
 tee ${t}not_bonds.txt |\
 awk '{\
     a1=$1;f1=$2;c1=$4;r1=$5;\
     a2=$7;f2=$8;c2=$10;r2=$11;\
     v=$13;s=99;sym=$14;\
     if(v>4)next;\
     if(v<0.1)v=0.1;\
   print "    #",c1,c2,a1,a2,r1,r2,f1,f2,v,s;\
   print "    bond {"\
   print "      action = *add";\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   if(sym!="")print "      symmetry_operation =",sym;\
   print "      distance_ideal =",v;\
   print "      sigma = 99";\
   print "      slack = 3";\
   print "    }";}' |\
awk '{gsub("and altid _","");gsub("and chain _","");print}' |\
awk 'BEGIN{print "  geometry_restraints.edits {";\
           print "    excessive_bond_distance_limit=15";\
       } {print}\
       END{print "  }"}' |\
awk 'BEGIN{print "refinement {"} {print} END{print "}"}' |\
cat >! $outfile

set count = `cat ${t}not_bonds.txt | wc -l`
echo "$count potential bumps suppressed in $outfile"


exit:

if("$tempfile" == "") set  tempfile = "./"
set tempbase = `basename $tempfile`
set tempdir = `dirname $tempfile`
if(! $debug && ! ( "$tempdir" == "." && "$tempbase" == "" ) ) then
    echo "cleaning up"
    rm -f ${tempfile}*
endif

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif


