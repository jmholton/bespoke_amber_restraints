#! /bin/tcsh -f
#
#
#  re-organize conformers in a PDB file
#
#
set pdbfiles = ""
foreach arg ( $* )
    if("$arg" =~ *.pdb) set pdbfiles = ( $pdbfiles "$arg" )
end
if(! -e "$pdbfiles[1]") then
    echo "ERROR: pdb file $pdbfiles does not exits"
    exit 9
endif

cat $pdbfiles |\
egrep "^ANISO" >! tempfile_$$_aniso.txt

egrep "^CRYST1" $pdbfiles[1] | head -1
cat $pdbfiles |\
awk '/^ATOM|^HETAT/{++n;printf("ATOM%7d%s\n",n,substr($0,12))}' |\
awk 'BEGIN{conflib="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_.,:;=<>+-%&*$!(){}[]^#|~?"}\
     {id=substr($0,12,15);while(seen[id]){\
        newconf=substr(conflib,index(conflib,substr($0,17,1))+1,1);\
        if(newconf=="")newconf="?";\
        $0=substr($0,1,16) newconf substr($0,18);\
        id=substr($0,12,15);\
       }\
       ++seen[id];\
       print;\
     }' |\
awk '{c=substr($0,17,1)} c==" "{c="\!"}\
     {o=substr($0,55,6);r=substr($0,23,4);s=substr($0,22,1);a=substr($0,5,7)}\
     o+0>0{print s,r,o,c,a,$0}' |\
sort -k1,1 -k2,2g -k3,3gr -k4,4 -k5,5g |\
awk '{print substr($0,index($0,"ATOM"))}' |\
awk '{++n;line[n]=$0;conf=substr($0,17,1);\
     id=substr($0,18,10)}\
     ! seen[conf id]{++seen[conf id];++confs[id];\
        occs[id]=occs[id]" "substr($0,55,6)}\
     END{for(i=1;i<=n;++i){\
          id=substr(line[i],18,10);\
          print line[i],"CONFORMERS:",occs[id];\
     }}' |\
awk 'BEGIN{conflib="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_.,:;=<>+-%&*$!(){}[]^#|~?"}\
     {res=substr($0,23,4);conf=substr($0,17,1);\
      atom=substr($0,12,5);occ=substr($0,55,6)+0}\
      res!=lastr{c=0}\
      conf!=lastc || res!=lastr {++c}\
     {newconf=substr(conflib,c,1)} newconf==""{newconf="?"}\
     /CONFORMERS/{nconf=split(substr($0,index($0,"CONFORMERS:")+12),occs);\
        $0=substr($0,1,index($0,"CONFORMERS")-1)}\
     nconf==1 && newconf=="A" && occ==1{newconf=" "}\
     nconf>1 && conf==" " && occ==1{\
	for(i=1;i<nconf;++i){\
           newconf=substr(conflib,i,1);\
           if(newconf=="")newconf="?";\
           print substr($0,1,16) newconf substr($0,18,39) occs[i+1] substr($0,61);\
        }\
        next;\
     }\
     {print substr($0,1,16) newconf substr($0,18,62)}\
     {lastc=conf;lastr=res}' |\
cat tempfile_$$_aniso.txt - |\
awk '{id=substr($0,12,16);id0=substr($0,12,5)" "substr($0,18,10)}\
  /^ANISO/{anis[id]=anis0[id0]=substr($0,28);next} {print}\
     anis[id]=="" && anis0[id0] != "" {anis[id]=anis0[id0]}\
     anis[id]!=""{printf("ANISOU%5d%s%s\n",$2,id,anis[id])}' |\
cat
echo "END"

rm -f tempfile_$$_aniso.txt

