#! /bin/awk -f
#
#  jiffy script for re-numbering and filtering large PDB files in differnt conventions
#
#
BEGIN{
   # may be: pdb amber
   if(! output) output = "pdb"
   # may be: "none" chain w4 w8 h36 dedupe ordinal
   if(! renumber) renumber = "none"
   # may be: all atoms protein notprotein ligand notligand water notwater notEP zeroH
   if(! only) only = "all"
   if(! skip) skip = "none"

   # minimum residue number
   if(minresnum=="") minresnum = 1;
   if(maxchain4resnum=="") maxchain4resnum=9999;

   # default to try to remove duplicates
   if(dedupe=="") dedupe = 1;

   # for pattern matching
   only = ","only","
   skip = ","skip","
   renumber = ","renumber","
   append = ","append","

   # for renumbering
   offset = 0;

   # default
   if(only !~ /protein,|ligand|water,|EP,|H,/) only = ","only",all,"
   if(only ~ /,all,/) only = ","only",protein,ligand,water,EP,zeroH,"

   # setting variables as shortcuts
   if(! ordinal) ordinal = ( renumber ~ /,ordinal,/ )
   if(! w4) w4 = ( renumber ~ /,w4/ )
   if(! w8) w8 = ( renumber ~ /,w8/ )
   if(! watS) watS = ( renumber ~ /,watS/ )
   if(! PLW) PLW = ( renumber ~ /,PLW,/ )
   if(! hy36) hy36 = ( renumber ~ /,hy36/ )
   if(! rechain) rechain = ( renumber ~ /,chain,/ )
   if(! nolock) nolock = ( renumber ~ /,nolock,/ )
   if( renumber ~ /,none,/ ) {
      ordinal=w4=w8=watS=hy36=rechain=dedupe=0;
   }

   # always start at N terminus
   nterm=2
   protein=prevprotein=1

   # build up alphabet
   for(c=65;c<=90;++c) ABC = ABC sprintf("%c",c);
   for(c=97;c<=122;++c) ABC = ABC sprintf("%c",c);
   ABCabc=ABC;
   for(c=23;c<=128;++c) {
     ch=sprintf("%c",c);
     if(ch !~ /[[:print:]]/) continue;
     if(ch ~ /[[:cntrl:]]/) continue;
     if(ch ~ /[[:digit:]]/) continue;
     if(ch ~ /[\!:\"\047\$\&\*]/) continue;
     if(index(ABC,ch)) continue;
     if(ch !~ /[[:graph:]]/) continue;
     ABC = ABC sprintf("%c",c)
   }
   for(c=1;c<=255;++c) {
     ch=sprintf("%c",c);
     if(ch !~ /[[:print:]]/) continue;
     if(index(ABC,ch)) continue;
     if(ch !~ /[[:graph:]]/) continue;
     ABC = ABC sprintf("%c",c)
   }
}


# special input cards for influencing re-naming
/^PROTON/{
  for(i=2;i<=NF;++i){
    resid= substr($i,1,1)" "substr($i,2);
    ++proton[resid];
  }
  next
}
/^DEPROTON/{
  for(i=2;i<=NF;++i){
    resid= substr($i,1,1)" "substr($i,2);
    ++deproton[resid];
  }
  next
}
/^HID/{
  resid= substr($2,1,1)" "substr($2,2);
  ++hid[resid];
  next
}
/^HIE/{
  resid= substr($2,1,1)" "substr($2,2);
  ++hie[resid];
  next
}
/^HIP/{
  resid= substr($2,1,1)" "substr($2,2);
  ++hip[resid];
  next
}
/^CYX/{
  resid= substr($2,1,1)" "substr($2,2);
  ++cyx[resid];
  next
}

/^ONLY/{only = ","$2",";next}
/^SKIP/{skip = ","$2",";next}
/^RENUMBER/{renumber = ","$2",";
   ordinal = ( renumber ~ /,ordinal,/ )
   w4 = ( renumber ~ /,w4/ )
   w8 = ( renumber ~ /,w8/ )
   watS = ( renumber ~ /,watS/ )
   hy36 = ( renumber ~ /,hy36/ )
   rechain = ( renumber ~ /,chain,/ )
   next;
}


/^END/{++seenend}

/^TER/ && renumber ~ /,ter,/ {
  new_chain=1;
  restart_counter=1;
  if(debug > 5) print "DEBUG restart counter and chain"
}

# scrape out non-atom records if asked
! /^ATOM|HETAT/ && only ~ /,atoms,/ {next}
# otherwise, always pass through

# skip unneccesary TER records
/^TER/ && output == "pdb" && printedatoms < 10 { printedatoms=0; next }

# shorten TER records
/^TER/ && skip !~ /,ter,/{
  print "TER";
  printedatoms = 0;
  next;
}

! /^ATOM|HETAT/{print;next}


# parse the atom line
{
    if(debug) print "DEBUG" substr($0,6);
    pre=substr($0,1,11);
    aftB=substr($0,67,10);
    post=substr($0,81);
      while(gsub(/ $/,"",post));post="       " post;
    origid=substr($0,12,19);
    atom=substr($0,12,5);
      atm=atom;gsub(" ","",atm);
    resnum=substr($0,23,8);
    if(substr(resnum,1,4)~/[A-Za-z]/) {
      gsub(" ","",resnum);
      resnum=encresnum=hy36decode(resnum);
    }
    resnum+=0; resnum0=resnum;
    chain = chain0 = substr($0,22,1);
    conf = substr($0,17,1);
    typ = substr($0,18,3);
    X = substr($0, 31, 8)+0;
    Y = substr($0, 39, 8)+0;
    Z = substr($0, 47, 8)+0;
    occ = substr($0,55,6)+0;
    Bfac = substr($0,61,6)+0;
    Ee = substr($0,77,4);gsub(" ","",Ee);
    atomEe=substr($0,13,2);gsub(" ","",atomEe);
    id=atm" "chain" "resnum;
    resid=   chain" "resnum;
    prevprotein = protein;
    protein=water=ligand=EP=0;
    cidx=index(ABC,chain);
}

# fix common bugs
typ=="MSE" && atm=="SE" && Ee=="S"{
    Ee = atm;
}


# classify the residue
typ~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein}
typ~/HID|HIE|HIP|HIS|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein}
typ~/HOH|WAT/{++water}
! protein && ! water{++ligand}

debug > 5 {print "DEBUG: protein",protein,"ligand",ligand,"water",water,"prevprotein",prevprotein}

# skip if asked
protein && skip ~ /,protein,/{next}
water && skip ~ /,water,/{next}
ligand && skip ~ /,ligand/{next}
Ee=="XP" && skip ~ /,EP,|,H,/{next}
atm=="EPW" && skip ~ /,EP,|,H,/{next}
Ee=="H" && occ==0 && skip ~ /,zeroH,/{next}
Ee=="H" && skip ~ /,H,/{next}


protein   && ( only !~ /,protein,/ ) {next}
water     && ( only !~ /,water,/ ) {next}
ligand    && ( only !~ /,ligand/) {next}
protein && only ~ /notprotein|nonprotein/ {next}
ligand && only ~ /notligand|nonligand/ {next}
water && only ~ /notwater|nonwater/ {next}
Ee=="XP" && only ~ /noEP|notEP|nonEP/{next}
Ee=="H" && occ==0 && only ~ /nonzeroH|nozeroH|notzeroH/{next}

debug > 5 {print "DEBUG: not skipped"}

rechain || renumber ~ /,ter,/ && prevochain!="" {
   if(debug>5) print "DEBUG rechain=",rechain,"renumber-ter=",(renumber ~/,ter,/),"prevochain=",prevochain
   cidx = index(ABC,prevochain);
}

# default to chain A if renumbering
watS || rechain {cidx=1}
watS && ligand {cidx=index(ABC,"L")}
watS && water  {cidx=index(ABC,"S")}
#PLW || rechain {cidx=1}
#PLW && ligand {cidx=index(ABC,"C")}
#PLW && water  {cidx=index(ABC,"D")}

# option to re-set things
CONF  != ""{conf=CONF}
CHAIN != ""{
  chain=CHAIN;
  cidx=index(ABC,chain);
}
OCC   != ""{occ=OCC}
BFAC  !="" {Bfac=BFAC}

# increment chain ID for whatever reason
new_chain {
   cidx = index(ABC,prevochain);
   chain=substr(ABCabc ABC,cidx,1);
   while( seen[chain] && chain != "" ) {
      ++cidx;
      chain=substr(ABCabc ABC,cidx,1);
   }
   new_chain=0;
   if(debug>5) print "DEBUG: new chain",cidx,chain
}

# convert chain index to actual chain letter
{chain=substr(ABCabc ABC,cidx,1)}

# always count up total residues by looking at input resids
prevresid != resid{
   ++ordresnum;
   # try not to re-use same chain-resnum far apart in file
   ++taken[prevoresid];
   if(nterm) --nterm;
   if(debug>10) print "DEBUG incrementing ordinal resnum to:",ordresnum
   if(debug>10 && nterm==0) print "DEBUG no longer N-terminus"
}

# user selected re-start at 1 for each chain
debug>5{print "DEBUG chainrestart=",( renumber~/,chainrestart,/ ),"seen["chain"]=",seen[chain]}
renumber~/,chainrestart,/ && ! seen[chain] {
    ++seen[chain];
    ++restart_counter;
}

# here is where counters get re-started for whatever reason
restart_counter {
    offset=1-ordresnum;
    restart_counter=0;
    if(debug>5) print "DEBUG restarting ordinal counter"
}


# re-name from amber to pdb conventions
output=="pdb" && typ=="WAT" {typ="HOH"}
output=="pdb" && typ~/HI[DEPS]/ {typ="HIS"}
output=="pdb" && typ=="CYX" {typ="CYS"}
output=="pdb" && typ=="ASH" {typ="ASP"}
output=="pdb" && typ=="GLH" {typ="GLU"}

debug>10{print "DEBUG output=",output,"nterm=",nterm,"atm=",atm}
output=="pdb" && nterm && protein && atm=="H1" && ordresnum!=last_cterm {
    if(debug>10) print "DEBUG: changing N-terminal H1 to H "
    # only do this for N-terminal NH3
    atom="  H  ";
    nterm=0;
}
output=="amber" && nterm && protein && atm=="H" && ordresnum!=last_cterm {
    if(debug>10) print "DEBUG: changing N-terminal H to H1"
    # only do this for N-terminal NH3
    atom="  H1 ";
    nterm=0;
}

output=="pdb" && typ=="ACT" { typ="ACY" }
typ=="ACY" && atm=="HB1" {atom="  H1 "}
typ=="ACY" && atm=="HB2" {atom="  H2 "}
typ=="ACY" && atm=="HB3" {atom="  H3 "}
typ=="ACY" && atm=="CA"  {atom=" C   "}
typ=="ACY" && atm=="CB"  {atom=" CH3 "}
typ=="ACY" && atm=="OA1" {atom=" O   "}
typ=="ACY" && atm=="OA2" {atom=" OXT "}

output=="pdb" && typ=="AMM" { typ="NH4" }
typ=="NH4" && atm=="H1"  {atom=" HN1 "}
typ=="NH4" && atm=="H2"  {atom=" HN2 "}
typ=="NH4" && atm=="H3"  {atom=" HN3 "}
typ=="NH4" && atm=="H4"  {atom=" HN4 "}


output=="amber" && typ=="ACY" && ! norenamelig { typ="ACT" }
typ=="ACT" && atm=="H1" {atom=" HB1 "}
typ=="ACT" && atm=="H2" {atom=" HB2 "}
typ=="ACT" && atm=="H3" {atom=" HB3 "}
typ=="ACT" && atm=="C"  {atom=" CA  "}
typ=="ACT" && atm=="CH3"{atom=" CB  "}
typ=="ACT" && atm=="O"  {atom=" OA1 "}
typ=="ACT" && atm=="OXT"{atom=" OA2 "}

output=="amber" && typ=="NH4" && ! norenamelig { typ="AMM" }
typ=="AMM" && atm=="HN1" {atom=" H1  "}
typ=="AMM" && atm=="HN2" {atom=" H2  "}
typ=="AMM" && atm=="HN3" {atom=" H3  "}
typ=="AMM" && atm=="HN4" {atom=" H4  "}

output=="amber" && typ=="HOH"{typ="WAT"}

# user-selectable protonation
proton[resid] && typ=="ASP" {typ="ASH"}
proton[resid] && typ=="GLU" {typ="GLH"}
proton[resid] && typ~/HI[SDEP]/ {typ="HIP"}

hid[resid] && typ~/HI[SDEP]/ {typ="HID"}
hie[resid] && typ~/HI[SDEP]/ {typ="HIE"}
hip[resid] && typ~/HI[SDEP]/ {typ="HIP"}

cyx[resid] && typ~/CY[XS]/ {typ="CYX"}

# ignore special hydrogens
typ=="HID" && atm=="HE2"{next}
typ=="HIE" && atm=="HD1"{next}
typ=="KCX" && atom=="HQ2"{next}


# switch to detect when we have decided resnum
{locked=0}

# use previously determined mapping, if available
mapping[ordresnum] {
    oresid= mapping[ordresnum];

    chain = substr(oresid,1,1);
    resnum = substr(oresid,3)+0;
    id = atm" "chain" "resnum" "conf;
    if(! nolock) ++locked;
    if(debug>3) print "DEBUG locked in resid",oresid,"chain=",chain,"ordinal=",ordresnum;
}

# user selected renumber to ordinal increasing residue
! locked && ordinal{
    resnum = ordresnum-1+minresnum+offset;
    if(debug>3) print "DEBUG ordinal and offset resnum=",resnum
}

# user selected funky hex-like chain numbers if they get too big
! locked && hy36 {
    encresnum = hy36decode( resnum );
    resnum = hy36encode( resnum );
}  

# user selected renumber, but always fit resnum in 4 ltters
! locked && w4 && resnum>maxchain4resnum || dedupe && seen[id] && ! rechain {
  id= atm" "chain" "resnum;
  oresid=   chain" "resnum;
  # register this as seen, even if we didnt print it
  if(resnum>maxchain4resnum) {
    if(debug>5) print "DEBUG: too long resnum",id;
    id= atm" "chain" "resnum;
    if(debug>5) print "DEBUG: trying",id;
    ++seen[id];
  }
  # try keeping original?
  if( seen[id] && ! ordinal && resnum<=maxchain4resnum ) {
    resnum=resnum0;
    id= atm" "chain" "resnum;
    oresid=   chain" "resnum;
    if(debug>5) print "DEBUG: trying original resnum:",id;
  }
  if( seen[id] && ! ordinal && resnum<=maxchain4resnum ) {
    chain=chain0;
    resnum=resnum0;
    id= atm" "chain" "resnum;
    oresid=   chain" "resnum;
    if(debug>5) print "DEBUG: trying original chain and resnum:",id;
  }
  # try incrementing previous
  if((resnum>maxchain4resnum || seen[id] || taken[oresid]) && prevoresnum) {
     if(debug>5) print "DEBUG: seen["id"]=",seen[id],"taken["oresid"]=",taken[oresid];
     resnum=prevoresnum+1;
     id= atm" "chain" "resnum;
     oresid=   chain" "resnum;
     if(debug>5) print "DEBUG: trying incremented previous output resnum:",id;
  }
  # try incrementing previous while keeping previous chain
  if((resnum>maxchain4resnum || seen[id] || taken[oresid]) && prevochain && prevoresnum) {
     chain=prevochain;
     resnum=prevoresnum+1;
     id= atm" "chain" "resnum;
     oresid=   chain" "resnum;
     if(debug>5) print "DEBUG: trying previous chain and incremented previous resnum:",id;
  }
  # start at zero if still no good
  if(resnum>maxchain4resnum || seen[id] || taken[oresid]) {
    resnum=minresnum;
    id= atm" "chain" "resnum;
    oresid=   chain" "resnum;
    if(debug>5) print "DEBUG: rewinding to zero:",id,"taken["oresid"]=",taken[oresid];
  }
  # try to stay within chain first
  while(seen[id] || taken[oresid]){
    if(debug>5) print "DEBUG: seen["id"]=",seen[id],"taken["oresid"]=",taken[oresid]
    ++resnum;
    id= atm" "chain" "resnum;
    oresid=   chain" "resnum;
    if(resnum>maxchain4resnum) {
      # too big, so pretend we printed it
      ++seen[id];
      break;
    }
  }
  # again, try keeping original number if ordinal counter not selected
  if(resnum>maxchain4resnum && ! ordinal) resnum=resnum0;
  # if being ordinal, reset to zero
  if(resnum>maxchain4resnum) resnum=minresnum ;
  id= atm" "chain" "resnum;
  oresid=   chain" "resnum;
  if(resnum>maxchain4resnum) ++taken[oresid];
  while(seen[id] || taken[oresid]){
    if(debug>5) print "DEBUG: seen["id"]=",seen[id],"taken["oresid"]=",taken[oresid]
    # now loop over chains
    ++cidx;
    chain=substr(ABCabc ABC,cidx,1);
    if(chain=="") {
      # nothing left but hy36 ?
      print "REMARK unavoidable duplicate:", atm,chain0,resnum0;
      chain="!";
      resnum=resnum0;
      # must break out of loop
      id= atm" "chain" "resnum;
      oresid=   chain" "resnum;
      seen[id]=0;
      break;
    };
    id= atm" "chain" "resnum;
    oresid=   chain" "resnum;
    if(resnum>maxchain4resnum) ++taken[oresid];
    if(debug>5) print "DEBUG: seen["id"]=",seen[id],"taken["oresid"]=",taken[oresid]
    if(seen[id] && renumber !~ /,ordinal,/) {
      # try keeping original resnum
      resnum=resnum0;
      id= atm" "chain" "resnum;
      oresid=   chain" "resnum;
      if(resnum>maxchain4resnum) ++taken[oresid];
    }
    while(seen[id] || taken[oresid]){
      if(debug>5) print "DEBUG: seen["id"]=",seen[id],"taken["oresid"]=",taken[oresid]
      ++resnum;
      id= atm" "chain" "resnum;
      oresid=   chain" "resnum;
      if(resnum>maxchain4resnum) {
        ++seen[id];
        break;
      }
    }
    if(seen[id] || taken[oresid]) {
      # try starting at zero
      resnum=minresnum;
      id= atm" "chain" "resnum;
      oresid=    chain" "resnum;
    }
    while(seen[id] || taken[oresid]){
      if(debug>5) print "DEBUG: seen["id"]=",seen[id],"taken["oresid"]=",taken[oresid]
      ++resnum;
      id= atm" "chain" "resnum;
      oresid=    chain" "resnum;
      if(resnum>maxchain4resnum) {
        ++seen[id];
        break;
      }
    }
    if(debug && seen[id]) print "DEBUG chain",chain,"no good"
  }
  id = atm" "chain" "resnum" "conf;
  oresid=    chain" "resnum" "conf;
}

! locked && w8 && length(resnum)>8{
   ++cidx;
   chain=substr(ABCabc ABC,cidx,1);
   resnum=minresnum;
}

! locked {
  id= atm" "chain" "resnum" "conf;
  oresid=   chain" "resnum" "conf;
}

! locked && rechain && seen[id] && renumber !~ /,none,/{
  loop1=0;
  while(seen[id]){
    ++cidx;
    chain=substr(ABCabc ABC,cidx,1);
    if(chain=="" && ! loop1) {cidx=0;++loop1};
    if(chain=="") {print "REMARK unavoidable duplicate";break};
    id = atm" "chain" "resnum" "conf;
    oresid=    chain" "resnum" "conf;
  }
}

# enforce mappings
#mapping[ordresnum] {
#    oresid= mapping[ordresnum];
#
#    chain = substr(oresid,1,1);
#    resnum = substr(oresid,3)+0;
#    id = atm" "chain" "resnum" "conf;
#}

# enforce not changing anything
renumber~/,none,/ {
    chain = chain0;
    resnum = resnum0;
    id = atm" "chain" "resnum" "conf;
}

# command line override of above
CHAIN != ""{
  chain=CHAIN;
  cidx=index(ABC,chain);
}

#  check for duplicates
dedupe && seen[id]{
  while( seen[id] ){
    ++resnum;
    id = atm" "chain" "resnum" "conf;
  }
}

# final check for duplicates
seen[id] && renumber !~ /,none,/{
   print "REMARK: ERROR duplicate atom id:",id
}

# stick terminators between chains
( output == "amber" || renumber ~ /,terify,/ ) && nexTER && prevresid != resid{
   print "TER";
   printedatoms = nexTER = 0;
}
# stick terminators between waters
output == "amber" && ! prevprotein && prevresid != resid{
   if(printedatoms || nexTER) print "TER";
   printedatoms = nexTER = 0;
}


fixEe && atomEe != Ee {
  if(debug>5) print "DEBUG fixEe atomEe=",atomEe,"Ee=",Ee,"typ=",typ
  typEe = typ
  gsub(" ","",typEe);
  if(typEe == atm) atomEe=typEe;
  if(typ ~ /WAT|HOH/ && atomEe!~/^E/) Ee=atomEe;
  if(atomEe=="C" && Ee=="CL" && substr(atm,1,2)==Ee) {
    if(debug) print "stupid AddToBox |"atom"|"
    atomEe=Ee;
    atom=sprintf(" %-4s",atm);
    if(debug) print "stupid AddToBox |"atom"|"
  }
  if(substr(atm,1,1) == Ee || atm == Ee || Ee=="XP" ) {
    # this is OK
  }
  else
  {
    Ee = atomEe;
  }
}

fixEe && atm == Ee {
  typEe = typ
  gsub(" ","",typEe);
  if(typEe == Ee || length(Ee)==2 && typEe=="MSE") {
    if(debug) print "DEBUG justifying"
    # standardize justification
    atom=sprintf("%-4s",atm);
    typ =sprintf("%-3s",typEe);
    if(debug>1) print "DEBUG",substr($0,7);
  }
}


{apdx=""}

append ~ /,ordresnum,/ {
   apdx = apdx" "sprintf("%10s",ordresnum)
}
append ~ /,encresnum,/ {
   apdx = apdx" "sprintf("%10s",encresnum)
}
append ~ /,origid,/ {
   apdx = apdx"     |" origid;
}

# actual, final printing out
{
   presnum=sprintf("%4s",resnum);
   s=8-length(presnum);
   for(;s>0;--s)presnum=presnum" ";
   printf("%11s%5s%1s%-4s%1s%8s%8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s%s%s\n", pre,atom,conf,typ,chain,presnum,X,Y,Z,occ,Bfac,aftB,Ee,post,apdx);
   ++printedatoms;
   id = atm" "chain" "resnum" "conf;
   
   ++seen[id];
   ++seen[chain];
   if(! mapping[ordresnum]) mapping[ordresnum]= oresid;
}


atm == "OXT"{
  # next residue is an N terminus
  nterm=2;
  last_cterm=ordresnum;
  # maybe need a TER
  if(printedatoms && only !~ /,atoms,/) ++nexTER
  printedatoms = 0;
}



# keep track of previous for next atom
{
  previd = id;
  prevresid = resid;
  prevtyp = typ;
  prevoresnum = resnum;
  prevochain = chain;
  prevoresid = chain" "resnum;
}

END{
  if( only ~ /,atoms,/ ) exit
  if(printedatoms || nexTER) print "TER"
  if(! seenend ) print "END"
}


# these work, but gemmi cannot understand lower case
function hy36decode( code )
{
  if(! match(code,/[A-Za-z]/)) return code;
  udigits="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  n=length(code);
  value=0;
  ucode=toupper(code);
  m=match(code,/[a-z]/);
  for(i=1;i<=n;++i) {
    c = substr(ucode,i,1);
    q=index(udigits,c);
    #print "DEBUG: c=",c,"q=",q;
    value+=(q-1+m*35)*(36**(n-i));
  }
  #print "DEBUG:",value,n,m
  return value - (1+m)*10*36**(n-1) + 10**n + m;
}

function hy36encode( value )
{
  n=4;code="";
  digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  if(value < 10**n) return value;

  value = value + 10*36**(n-1) - 10**n;
  #print "DEBUG: value=",value;
  if(value>=36**n){digits=tolower(digits);value=value-36**4+10*36**(n-1)};
  for(i=1;i<=n;++i) {
    q = int( value / 36**(n-i) )
    c = substr(digits,q+1,1);
    #print "DEBUG: value=",value,"c=",c,"q=",q;
    code = code c;
    value = value - q*36**(n-i)
  }
  return code;
}

