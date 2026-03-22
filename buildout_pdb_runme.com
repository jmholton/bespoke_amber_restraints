#! /bin/tcsh -f
#
#  build all missing atoms while avoiding symmetry mates
#
set set pdbfile = ""
set sequence = ""
set firstresnum = ""
set ciffiles = ""

set tempfile = tempfile
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
      if("$Arg" =~ *.pdb ) then
        set pdbfile = $Arg
        echo "pdbfile = $pdbfile"
        continue
      endif
      if("$Arg" =~ *.seq || "$arg" =~ *.fasta || "$arg" =~ *.fa ) then
        set sequence = $Arg
        echo "sequence = $sequence"
        continue
      endif
     if("$Arg" =~ *.cif ) then
        set ciffiles = ( $ciffiles $Arg )
        echo "pdbfile = $pdbfile"
        continue
      endif
      if("$Arg" =~ *.mtz ) then
        set mtzfile = $Arg
        echo "mtzfile = $mtzfile"
        continue
      endif
    endif
    if("$key" == "debug") set debug = "1"
end

# shorthand for temporary file
set t = $tempfile

# make sure other scripts are in the $path
set path = ( $path `dirname $0` )

foreach dependency ( sequence.awk filter_pdb.awk build_n2c.awk build_c2n.awk build_side.awk )
   echo -n "using: "
   which $dependency
   if( $status ) then
       set BAD = "need $dependency in "'$'"path"
       goto exit
   endif
end

set test = `grep SEQRES $pdbfile | wc -l`
if( ! $test ) then
  set BAD = "no SEQRES in $pdbfile"
  goto exit
endif

if( "$firstresnum" == "" ) then
  echo "trying to get first residue number from SEQRES"
  set firstresnum = `sequence.awk $pdbfile | awk '/^seqres start/{print $NF}'`
  echo "it is $firstresnum"
endif
if( "$firstresnum" == "" ) set firstresnum = 1


if(! -e "$sequence") then
  # get full sequence from PDB
  grep SEQRES $pdbfile |\
  sequence.awk |\
  awk 'NR==1{$0="> "$0} {print} NF==0{exit}' |\
  cat >! seq.fasta
  set sequence = seq.fasta
  awk '/SEQRES/{print substr($0,18)}' $pdbfile |\
  awk '{for(i=1;i<=NF;++i)print $i}' >! seq_TLC.txt
endif

if( ! -e seq_TLC.txt && -e "$sequence" ) then
  sequence.awk -v tlc=1 $sequence | awk '/^[A-V][A-Y][A-Y]$/' >! seq_TLC.txt
endif
set nres = `cat seq_TLC.txt | wc -l`
set lastresnum = `echo $firstresnum $nres | awk '{print $1+$2-1}'`

# maybe sanity check here
echo "sequence: $nres residues from $firstresnum to $lastresnum "
cat $sequence

# take first protien chain
set chain = `filter_pdb.awk -v only=protein,atoms $pdbfile | awk '{print substr($0,22,1)}' | head -n 1`

# build everything as a helix
awk '{print "BUILD",$1,-60,-40}' seq_TLC.txt | build_n2c.awk >! backbone.pdb
awk '{print "BUILD",$1, ++n,"? ? ? ? ?"}' seq_TLC.txt |\
cat - backbone.pdb |\
build_side.awk >! side.pdb
echo "renumber ${firstresnum}\nchain $chain" | pdbset xyzin side.pdb xyzout helix.pdb > /dev/null
# all-helix version of the full sequence

# now add split versions of every residue that is split
egrep "^CRYST1|^LINK|^SSB" $pdbfile >! complete.pdb
cat $pdbfile |\
awk '! /^ATOM|^HETAT/{next}\
  {residue=substr($0,22,8);f=substr($0,17,1)}\
  f!=" " && ! seen[f" "residue]{\
     print f,residue,"SPLIT"\
     ++seen[f" "residue]}' |\
cat - helix.pdb |\
awk '$NF=="SPLIT"{res=substr($0,3,8);confs[res]=confs[res]" "$1;next}\
  /^CRYST|^SSB|^LINK/{print}\
  ! /^ATOM|^HETAT/{next}\
  {res=substr($0,22,8);pre=substr($0,1,16);post=substr($0,18)}\
  {print}\
  ! confs[res]{next}\
  {n=split(confs[res],f);for(i=1;i<=n;++i){\
   print pre f[i] post;\
}}' >> complete.pdb
filter_pdb.awk -v only=atoms -v skip=protein $pdbfile >> complete.pdb


 
egrep "^ATOM|^HETAT" helix.pdb |\
awk '{c=substr($0,22,1)} c==" "{c="_"}\
   {print substr($0,12,5),substr($0,18,3),c,substr($0,23,6)}' |\
awk '{while(gsub("  "," "));gsub(" $","");print $0,"EXPECT"}' |\
sort -u |\
sort -k3,3 -k4g >! expected_atoms.txt
set test = `cat expected_atoms.txt | wc -l`
echo "$test non-H protein atoms expected from sequence"


# check if anything is missing
awk '/^ATOM|^HETAT/{print substr($0,1,16),substr($0,18)}' $pdbfile |\
awk '{id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
filter_pdb.awk -v only=protein -v skip=H |\
cat >! noalt.pdb

egrep "^ATOM|^HETAT" noalt.pdb |\
awk '{c=substr($0,22,1)} c==" "{c="_"}\
   {print substr($0,12,5),substr($0,18,3),c,substr($0,23,6)}' |\
awk '{while(gsub("  "," "));gsub(" $","");print $0,"EXIST"}' |\
sort -u |\
sort -k3,3 -k4g >! existing_atoms.txt
set test = `cat existing_atoms.txt | wc -l`
echo "$test non-H protein atoms accounted for in $pdbfile"
cp existing_atoms.txt existing_atoms0.txt 



# Flory coefficients for mean and rms inter-atom distances
cat << EOF >! Flory_coeff.txt
N O   5.82 1.17 0.46 1.56 2.34 0.58 
N C   5.79 1.35 0.46 1.56 2.59 0.58 
CA O  5.82 1.42 0.46 1.49 2.47 0.60 
CA C  5.78 1.60 0.46 1.48 2.71 0.60 
N CA  5.79 1.65 0.46 1.51 2.72 0.59 
C O   5.77 1.84 0.46 1.36 2.61 0.62 
O O   5.76 1.85 0.46 1.34 2.37 0.63 
CA CA 5.76 1.88 0.46 1.43 2.82 0.61 
C C   5.73 2.01 0.46 1.34 2.80 0.63 
O C   5.71 2.01 0.46 1.31 2.57 0.63 
N N   5.72 2.07 0.46 1.37 2.76 0.62 
O CA  5.65 2.25 0.47 1.27 2.66 0.64 
C CA  5.67 2.26 0.47 1.28 2.86 0.64 
CA N  5.65 2.27 0.47 1.29 2.83 0.64 
O N   5.44 2.55 0.48 1.16 2.75 0.67 
C N   5.49 2.58 0.48 1.15 2.88 0.67 
EOF

# convert to mean and rmsd distance between each particular atom pair
foreach gapincr ( `seq 1 30` )
  cat Flory_coeff.txt |\
  awk -v iN=$gapincr '{a1=$1;a2=$2;b=$3;x0=$4;v=$5;q=$6;x0=$7;n=$8;\
    print a1,a2,b*(iN-x0)**v,"+/-",q*(iN-x1)**n}' |\
  cat >! dist_sig_1-${gapincr}.txt
end

foreach gapincr ( `seq 3 20` )
  if(! -e dist_sig_1-${gapincr}.txt) continue
  sort -k5g dist_sig_1-${gapincr}.txt |\
  awk '{a1=$1;a2=$2;dist=$3;sigma=$5;\
        print "    bond {";\
        print "      action = add";\
        print "      atom_selection_1 = \"name",a1,"and chain A and resseq i1 \"";\
        print "      atom_selection_2 = \"name",a2,"and chain A and resseq iN \"";\
        print "      distance_ideal =",dist;\
        print "      sigma =",sigma;\
        print "    }";\
    }' >! bridge_1-${gapincr}.efftemplate

end


# cp existing_atoms0.txt existing_atoms.txt
echo -n "" >! buildorder.txt
set n = 0
set missing_atoms = 1
while ( $missing_atoms )

  @ n = ( $n + 1 )

  cat existing_atoms.txt expected_atoms.txt |\
  awk '{id=$1" "$2" "$3" "$4}\
      $NF=="EXIST"{++exist[id]}\
      $NF=="EXPECT" && ! exist[id]{print id,"MISSING"}' |\
  tee missing_atoms.txt |\
  awk '{print $2,$3,$4}' |\
  sort -u | sort -k1.6g |\
  cat >! missing_residues.txt
  set nmissatom = `cat missing_atoms.txt | wc -l`
  set nmissres = `cat missing_residues.txt | wc -l`
  echo "$nmissatom atoms in $nmissres residues still unplanned"
  if ( ! $nmissatom ) break

  # find the thing that will be easiest to build
  echo $firstresnum $lastresnum |\
  cat - existing_atoms.txt missing_residues.txt |\
  awk 'NR==1{f=$1;l=$2;nres=l-f+1;next}\
      $NF=="EXIST" && $1=="CA"{++existCA[$3" "$4]}\
      $NF=="EXIST"{++exist[$3" "$4];next}\
      {chain=$2;resnum=$3;p=n=resnum;\
        while(! existCA[chain" " n] && n<=l)++n;\
        while(! existCA[chain" " p] && p>=f)--p;\
        before=resnum-p;after=n-resnum;\
        if(! existCA[chain" "p])before=nres;\
        if(! existCA[chain" "n])after=nres;\
        gotCA=existCA[chain" "resnum]+0;\
        gotatoms=exist[chain" "resnum]+0;\
        terminal=( ! ( existCA[chain" "p] && existCA[chain" "n] ) );\
        gotCAprev=existCA[chain" "(resnum-1)]+0;\
        gotCAnext=existCA[chain" "(resnum+1)]+0;\
        score=before*after*100+gotatoms/2+$3/1000;\
        print $0,score,"  ",before,after,"  ",gotCA,gotatoms,gotCAprev,gotCAnext;\
      }' |\
  sort -k4g >! priorities.txt
  # TYP A num score   nbefore nafter   gotCA gotatoms gotCAprev gotCAnext
  #  1  2  3    4        5      6        7      8        9         19

   set info = `head -n 1 priorities.txt`
   set typ = $info[1]
   set chain = $info[2]
   set resnum = $info[3]

   set gotCA = $info[7]
   set gotCAprev = $info[9]
   set gotCAnext = $info[10]


   # predict gaps after this build
   echo "CA $typ $chain $resnum" |\
   sort -k3,3 -k4g existing_atoms.txt - >! sorted.txt
   set gaps = `awk '$3==lastC && $4>lastN+1{print $3 lastN"-"$4} {lastC=$3;lastN=$4}' sorted.txt`


   echo -n "looking to build $typ $chain $resnum "
   if( $debug ) echo "gotCA=$gotCA gotCAprev=$gotCAprev gotCAnext=$gotCAnext"


   # get main chain atoms that currently exist
   echo $chain $resnum |\
     cat - existing_atoms.txt |\
     awk 'NR==1{c=$1;r=$2;next}\
        $3==c && $4==r && ( $1~/^[CNOH]$/ || $1=="CA" ){print $1,$2,$3,$4,"EXIST"}' |\
   cat >! old_main.txt
   # get side chain atoms that currently exist
   echo $chain $resnum |\
     cat - existing_atoms.txt |\
     awk 'NR==1{c=$1;r=$2;next}\
        $3==c && $4==r && ! ( $1~/^[CNOH]$/ || $1=="CA" ){print $1,$2,$3,$4,"EXIST"}' |\
   cat >! old_side.txt

   # get main chain atoms that will be built
   echo $chain $resnum |\
     cat - expected_atoms.txt |\
     awk 'NR==1{c=$1;r=$2;next}\
        $3==c && $4==r && ( $1~/^[CNOH]$/ || $1=="CA" ){print $1,$2,$3,$4,"EXIST"}' |\
   cat >! new_main.txt
   # get side chain atoms that will be built
   echo $chain $resnum |\
     cat - expected_atoms.txt |\
     awk 'NR==1{c=$1;r=$2;next}\
        $3==c && $4==r && ! ( $1~/^[CNOH]$/ || $1=="CA" ){print $1,$2,$3,$4,"EXIST"}' |\
   cat >! new_side.txt

   set got_all_main = `cat old_main.txt new_main.txt | awk '{++seen[$1]} END{got=1;for(a in seen)if(seen[a]!=2)got=0;print got}'`
   set got_all_side = `cat old_side.txt new_side.txt | awk '{++seen[$1]} END{got=1;for(a in seen)if(seen[a]!=2)got=0;print got}'`
   if( $debug ) echo "got_all_main= $got_all_main got_all_side= $got_all_side"

   if( $#gaps ) echo "$n GAP $gaps" >> buildorder.txt

   if( $got_all_main && $got_all_side ) then
     echo "already taken care of"
     cat new_main.txt new_side.txt >> existing_atoms.txt
     continue
   endif
   if( $got_all_main && ! $got_all_side ) then
     echo "just need side"
     echo "$n SIDE $info" >> buildorder.txt
     cat new_side.txt >> existing_atoms.txt
     continue
   endif
   if(   $gotCAprev && ! $gotCA &&   $gotCAnext ) then
     echo "both neighbors exist: $gotCAprev $gotCA $gotCAnext"
     echo "$n JOIN $info" >> buildorder.txt
     cat new_main.txt >> existing_atoms.txt
     continue
   endif
   if( ! $gotCAprev && ! $gotCA &&   $gotCAnext ) then
     echo "C-terminal neighbor exists"
     echo "$n C2N $info" >> buildorder.txt
     cat new_main.txt new_side.txt  >> existing_atoms.txt
     continue
   endif
   if(   $gotCAprev && ! $gotCA && ! $gotCAnext ) then
     echo "N-terminal neighbor exists"
     echo "$n N2C $info" >> buildorder.txt
     cat new_main.txt new_side.txt  >> existing_atoms.txt
     continue
   endif
   if( ! $gotCAprev && ! $gotCA && ! $gotCAnext  ) then
     set BAD = "no neighbors exist for $info"
     goto exit
   endif
   set BAD = "dont know what to do with: $info"
   goto exit
end

#phipsichi.com noalt.pdb >! phipsichi.txt

egrep "^CRYST1|^LINK|^SSBO|^ATOM|^HETAT" $pdbfile >! initial.pdb
# perhaps filter out ligands?
cp initial.pdb incomplete.pdb

# make selection mask so that only new stuff is minimized
awk '$2~/N2C|C2N|JOIN|SIDE/{print $4,$5}' buildorder.txt |\
sort -u | sort -k1,1 -k2g |\
awk 'NF==2{s = "chain "$1" and resseq "$2}\
    {printf("(%s) or ",s)}' |\
awk '{print "selection = \"" $0 " water or (chain z)\""}' >! selection.eff


foreach n ( `awk '{print $1}' buildorder.txt | sort -u | sort -g` )

  set info = `egrep "^$n " buildorder.txt | awk '$2!="GAP"{$1="";print}'`
  set gaps = `egrep "^$n GAP " buildorder.txt | awk '{$1=$2="";print}'`

  echo "building round $n : $info"
  if( $debug ) echo "gaps: $gaps"

  echo "refinement {\n  geometry_restraints.edits {" >! gaps.eff
  echo "    excessive_bond_distance_limit = 99" >> gaps.eff
  foreach gap ( $gaps )
    # keep gaps from getting too wide to join
    set chain = `echo $gap | awk '{print substr($0,1,1)}'`
    set gaprange = `echo $gap | awk '{print substr($0,2)}' | awk -F "-" '{print $1,$2}'`
    set gapincr = `echo $gaprange | awk '{print $2-$1+1}'`

    foreach gapincr ( $gapincr )
      if(-e bridge_1-${gapincr}.efftemplate) then
        echo "using random-coil restraint for i--i+${gapincr} "
        echo $chain $gaprange |\
        cat - bridge_1-${gapincr}.efftemplate |\
        awk 'NR==1{c=$1;i1=$2;iN=$3;next}\
          {gsub("chain A","chain " c)\
           gsub("resseq i1","resseq "i1);\
           gsub("resseq iN","resseq "iN);\
           print}' >> gaps.eff
      else
        # always tight restraint at extended-chain limit
        echo "wide gap: $chain $gaprange"
        echo $chain $gaprange |\
        awk '{c=$1;nc=$2;nn=$3;sigma=1;\
          gapdist=3.8*(nn-nc+1);\
          slack=(gapdist-3.8)/2;\
          center=gapdist-slack;\
          print "    bond {";\
          print "      action = add";\
          print "      atom_selection_1 = \"name CA and resseq",nc,"and chain",c,"\"";\
          print "      atom_selection_2 = \"name CA and resseq",nn,"and chain",c,"\"";\
          print "      distance_ideal =",center;\
          print "      slack =",slack;\
          print "      sigma =",sigma;\
          print "    }";}' >> gaps.eff
      endif
    end
  end
  echo "  }\n}" >> gaps.eff

  if( $#info < 4 ) then
    echo "no build this round"
    set info = ( x x x x )
    continue
  endif
 
  set action = $info[1]
  set typ = $info[2]
  set chain = $info[3]
  set resnum = $info[4]
  @ prevresnum = ( $resnum - 1 )
  @ nextresnum = ( $resnum + 1 )

  # see if we are at the end
  set noOXT = `echo $resnum $lastresnum | awk '{print ( $1 != $2 ) }' `

  # isolate any existing atoms
  echo $chain $resnum |\
  cat - incomplete.pdb |\
  awk 'NR==1{c=$1;r=$2;}\
      ! /^ATOM|^HETAT/{next}\
    {chain=substr($0,22,1);resnum=substr($0,23,8)+0}\
    c==chain && r==resnum{print}' |\
  cat >! thisres.pdb

  # isolate appropriate N term to build on in C2N direction
  echo $chain $nextresnum |\
  cat - incomplete.pdb |\
  awk 'NR==1{c=$1;r=$2;}\
      ! /^ATOM|^HETAT/{next}\
    {chain=substr($0,22,1);resnum=substr($0,23,8)+0}\
    c==chain && r==resnum{print}' |\
  cat >! nextres.pdb

  # isolate appropriate C term to build upon in N2C direction
  echo $chain $prevresnum |\
  cat - incomplete.pdb |\
  awk 'NR==1{c=$1;r=$2;}\
      ! /^ATOM|^HETAT/{next}\
    {chain=substr($0,22,1);resnum=substr($0,23,8)+0}\
    c==chain && r==resnum{print}' |\
  cat >! prevres.pdb

  cp incomplete.pdb built.pdb

  if( "$action" == "JOIN" ) then
    cp built.pdb prebuild.pdb
    set prevtyp = `awk '/^ATOM|^HETAT/{print substr($0,18,3)}' prevres.pdb | head -n 1`
    set nexttyp = `awk '/^ATOM|^HETAT/{print substr($0,18,3)}' nextres.pdb | head -n 1`
    echo "BUILD $typ - -" >! build.txt
    echo "BUILD $nexttyp - -" >> build.txt
    cat prevres.pdb build.txt | build_n2c.awk -v noOXT=1 >! frag.pdb
    echo "renumber ${prevresnum}\nchain z" | pdbset xyzin frag.pdb xyzout tripeptide.pdb > /dev/null

    cat incomplete.pdb tripeptide.pdb >! alignme.pdb

    echo "refinement {\n  geometry_restraints.edits {" >! align.eff
    echo "    excessive_bond_distance_limit = 99" >> align.eff
    cat prevres.pdb nextres.pdb |\
    awk '! /^ATOM/{next}\
          {a=substr($0,12,5);c=substr($0,22,1);n=substr($0,23,5);\
            gsub(" ","",a);\
          sigma=0.01;}\
          a~/^[CNO]$/ || a~/^C[AB]/{\
          print "    bond {";\
          print "      action = add";\
          print "      atom_selection_1 = \"name",a,"and resseq",n,"and chain",c,"\"";\
          print "      atom_selection_2 = \"name",a,"and resseq",n,"and chain z \"";\
          print "      distance_ideal =",0.001;\
          print "      sigma =",sigma;\
          print "    }";}' >> align.eff
    echo "  }\n}" >> align.eff

    phenix.geometry_minimization alignme.pdb align.eff selection.eff $ciffiles \
      cdl=false apply_all_trans=True > align.log
    if( $status ) then
      set BAD = "error closing gap"
      goto exit
    endif

    cat alignme_minimized.pdb |\
    awk -v chain=$chain '! /^ATOM|^HETAT/{next}\
     {c=substr($0,22,1);pre=substr($0,1,21);post=substr($0,23)}\
     c!="z"{print;next}\
     {print pre chain post}' |\
    rmsd2B -v averageB=1 >! average.pdb

    combine_pdbs_runme.com average.pdb complete.pdb outfile=built.pdb > /dev/null    

    combine_pdbs_runme.com prebuild.pdb built.pdb xor=1 outfile=new.pdb > /dev/null 
    set newatoms = `egrep "^ATOM|^HETAT" new.pdb | wc -l`
    echo "$newatoms new atoms joined to $chain $resnum"
  endif

  if( "$action" == "N2C" ) then
    cp built.pdb prebuild.pdb

    cat prevres.pdb thisres.pdb nextres.pdb >! tripeptide.pdb
    phipsichi.com tripeptide.pdb | tee phipsichi.txt
    set phipsi = `awk '{print $0,"- - - - - - - - - - - - -"}' phipsichi.txt | awk '{print $5,$7}'`
    echo "BUILD $typ $phipsi" |\
    tee buildcmds_n2c_${n}.txt |\
    cat prevres.pdb - |\
    build_n2c.awk -v noOXT=$noOXT >! build_n2c.pdb
    combine_pdbs_runme.com build_n2c.pdb incomplete.pdb complete.pdb outfile=built.pdb > /dev/null

    combine_pdbs_runme.com prebuild.pdb built.pdb xor=1 outfile=new.pdb > /dev/null 
    set newatoms = `egrep "^ATOM|^HETAT" new.pdb | wc -l`
    echo "$newatoms new atoms n2c to $chain $resnum"
  endif
  if( "$action" == "C2N" ) then
    cp built.pdb prebuild.pdb

    cat prevres.pdb thisres.pdb nextres.pdb >! tripeptide.pdb
    phipsichi.com tripeptide.pdb | tee phipsichi.txt
    set phipsi = `awk '{print $0,"- - - - - - - - - - - - -"}' phipsichi.txt | awk '{print $5,$7}'`
    echo "BUILD $typ $phipsi\nBUILD -" |\
    tee buildcmds_c2n_${n}.txt |\
    cat - nextres.pdb |\
    build_c2n.awk >! build_c2n.pdb
    combine_pdbs_runme.com build_c2n.pdb incomplete.pdb complete.pdb outfile=built.pdb > /dev/null

    combine_pdbs_runme.com prebuild.pdb built.pdb xor=1 outfile=new.pdb > /dev/null 
    set newatoms = `egrep "^ATOM|^HETAT" new.pdb | wc -l`
    echo "$newatoms new atoms c2n to $chain $resnum"
  endif

  # always rebuild the side chain?
#  if( "$action" == "SIDE" || "$typ" == "PRO" ) then
    # isolate any existing atoms
    echo $chain $resnum |\
    cat - built.pdb |\
    awk 'NR==1{c=$1;r=$2;}\
        ! /^ATOM|^HETAT/{next}\
      {chain=substr($0,22,1);resnum=substr($0,23,8)+0}\
      c==chain && r==resnum{print}' |\
    cat >! thisres.pdb

    cp built.pdb prebuild.pdb

    phipsichi.com thisres.pdb | tee phipsichi.txt
    set chis = `awk '{print $0,"- - - - - - - - - - - - -"}' phipsichi.txt | awk '{print $11,$13,$15,$17,$19}'`
    echo "BUILD $typ ${chain}$resnum $chis" |\
    tee buildcmds_side_${n}.txt |\
    cat - built.pdb |\
    build_side.awk >! side.pdb
    combine_pdbs_runme.com side.pdb built.pdb complete.pdb outfile=built.pdb > /dev/null
#  endif

    combine_pdbs_runme.com prebuild.pdb built.pdb xor=1 outfile=new.pdb > /dev/null 
    set newatoms = `egrep "^ATOM|^HETAT" new.pdb | wc -l`
    echo "$newatoms new atoms side to $chain $resnum"

    cp built.pdb built_${n}.pdb

   # average over any duplicates?
   #rmsd2B -v averageB=1 incomplete.pdb built.pdb >! average.pdb
  
  if( "$action" != "SIDE" ) then 

    rm built_minimized.pdb >& /dev/null
    echo "minimizing geometry:"
    phenix.geometry_minimization built.pdb gaps.eff selecton.eff $ciffiles \
      write_geo_file=False cdl=false | tee geomin.log | egrep "target:" 
    if( ! -e built_minimized.pdb ) then
      set BAD = "error minimizing geometry"
      goto exit
    endif


    cp built_minimized.pdb built.pdb

  endif

  egrep -v "^LINK" built.pdb >! incomplete.pdb

end

# check if anything is missing
awk '/^ATOM|^HETAT/{print substr($0,1,16),substr($0,18)}' built.pdb |\
awk '{id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
filter_pdb.awk -v only=protein |\
awk '{c=substr($0,22,1)} c==" "{c="_"}\
   {print substr($0,12,5),substr($0,18,3),c,substr($0,23,6)}' |\
awk '{while(gsub("  "," "));gsub(" $","");print $0,"EXIST"}' |\
sort -u |\
sort -k3,3 -k4g >! existing_atoms.txt

cat existing_atoms.txt expected_atoms.txt |\
awk '{id=$1" "$2" "$3" "$4}\
    $NF=="EXIST"{++exist[id]}\
    $NF=="EXPECT" && ! exist[id]{print id,"MISSING"}' |\
cat >! missing_atoms.txt 
set nmissatom = `cat missing_atoms.txt | wc -l`
echo "$nmissatom missing atoms in built.pdb"


exit:

if("$tempfile" == "") set  tempfile = "./"
set tempbase = `basename $tempfile`
set tempdir = `dirname $tempfile`
if(! $debug && ! ( "$tempdir" == "." && "$tempbase" == "" ) ) then
    rm -rf ${tempfile}* >& /dev/null
endif

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif



exit

##########################################################################
#
#



grep "CA CA " stats_1-* | awk -F "-" '{print $2}' | awk '{print $1+0,$3}' | sort -g | tee plotme.txt



