#! /bin/tcsh -f
#
#

set rstfile = "$1"

if(! -e "$rstfile") then
  echo "usage: $0 amber.rst7"
  exit 9
endif

rm -f chircheck.txt geocheck.txt omega.txt
cpptraj -p xtal.prmtop << EOF >&! checks.log
trajin $rstfile lastframe
checkchirality chir out chircheck.txt
strip :WAT,HOH@Y1,EPW
check reportfile geocheck.txt
multidihedral omega omega out omegaline.txt range360
EOF
set ninv = `awk '$2!=1 && $1+0>0' chircheck.txt | wc -l`
set worstchir = `awk '$2!=1 && $1+0>0{print $1+0;exit}' chircheck.txt `
cat omegaline.txt |\
 awk 'NR==1{for(i=2;i<=NF;++i){split($i,w,":");oresnum[i]=w[2]};next}\
   {for(i=2;i<=NF;++i){dev=sqrt(($i-180)^2);print oresnum[i],dev,$i}}' >! omega.txt
set worstomega = `sort -k2gr omega.txt | head -n 1`
set ncis = `awk '$2>90{print}' omega.txt | wc -l`
set geoproblems = `cat geocheck.txt | wc -l`
set worstgeo = `awk '{print;exit}' geocheck.txt`
echo "$ninv inverted chiral centers, $ncis cis peptides ( $worstomega[2] deg) and $geoproblems geometry problems"


echo $worstomega |\
cat - orignames.pdb |\
awk 'NR==1{worstres=$1;next}\
  {o=$NF} printed[o]{next}\
  o==worstres{printf("%s - ", substr($0,18,9));++printed[o]}\
  o==worstres+1{print substr($0,18,11),"   omega    ",o;++printed[o]}\
  o>worstres+1{exit}'
echo $worstchir |\
cat - orignames.pdb |\
awk 'NR==1{worst=$1;next}\
  {o=$NF} printed[o]{next}\
  o==worst{print substr($0,18,11),"  chiral  ",o;++printed[o]}\
  o>worst{exit}'
echo $worstgeo |\
cat - orignames.pdb |\
awk 'NR==1{geo=$0;++bad[$4+0];++bad[$6+0];next}\
  /^ATOM|^HETAT/{++a;o=$NF} printed[o]{next}\
  bad[a]{print substr($0,18,11),"    atom:",a,"res:",o;++printed[o]}\
  printed[o]>2{exit}'
echo $worstgeo

if( "$rstfile" =~ *.nc ) then

echo "extracting omegas..."
cpptraj -p xtal.prmtop -y $rstfile << EOF >! cpptraj.txt
strip :WAT
strip @H=
multidihedral omega omega out omega.dat range360
EOF
cat omega.dat |\
 awk 'NR==1{for(i=2;i<=NF;++i){split($i,w,":");oresnum[i]=w[2]};next}\
   {++f;for(i=2;i<=NF;++i){dev=sqrt(($i-180)^2);if(dev>90)print f,oresnum[i],dev,$i;\
    if(dev>worst){worst=dev;worstf=f;worstres=oresnum[i]}}}\
   END{print worstf,worstres,worst,"worstomega"}' |\
head | tee first.txt
set oresnum = `awk '{print $2;exit}' first.txt`
awk -v oresnum=$oresnum 'NR==1{for(i=2;i<=NF;++i){split($i,w,":");idx[w[2]]=i};i=idx[oresnum];next}\
   {++f;{dev=sqrt(($i-180)^2);print f,dev,$i}}' omega.dat |\
cat >! worstomega_plot.txt

endif


exit

awk '/^ATOM|^HETAT/{print substr($0,22,1),substr($0,23,5),$NF}' orignames.pdb |\
 sort -u | sort -k3g >! ordresnums.txt

awk '/^TORS/ && $8=="CA" && $26=="CA"' tempfile_fullgeo.txt >! raw_omegas.txt

cat ordresnums.txt raw_omegas.txt |\
awk -v modulo=$modulo 'NF==3{ordresnum[$1,$2]=$3;next}\
  {o=ordresnum[$11,$12];\
    monomer=int(o/modulo);cpp=o-monomer;\
    print o,cpp,$10,$11,$12,monomer,$4}' |\
cat >! omega_count.txt

cat omega.txt omega_count.txt |\
awk 'NF==1{cppomega[NR]=$1;next}\
  {print $0,cppomega[$2]}' |\
tee compareme.txt



