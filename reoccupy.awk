#! /bin/awk -f
#
#
#  re-organize conformers in a PDB file
#
#
BEGIN{
    conflib = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_.,:;=<>+-%&*$!(){}[]^#|~?"
}

! /^ATOM|^HETAT|^ANISOU/ {
    print;
    next;
}

{
    conf=substr($0,17,1);
    atom=substr($0,12,5);
    res=substr($0,18,10);
    restyp=substr($0,18,3);

    id=atom res;
    if(ignorerestyp) id= atom substr($0,21,7)
}

/^ANISOU/{
     idc=substr($0,12,16)
     anisc[idc]=anis[id]=substr($0,28);
    next;
}

/^ATOM|^HETAT/{
    ++count[id];
    XYZ[id,count[id]]  = substr($0,31,24);
    occ[id,count[id]]  = substr($0,55,6);
    resocc[res,count[id]]  = substr($0,55,6);
    Bfac[id,count[id]] = substr($0,61);
    hetat[id]= ( /^HETAT/ );
}


! seen[id]{
    ++m;
    ido[m]=id
} 

{++seen[id]}

seen[id]+0>confs[res]+0{
    confs[res]=seen[id];
}

END{
    # loop over atom names
    for(i=1;i<=m;++i){

        id=ido[i];
        atom=substr(id,1,5);
        atm=atom;gsub(atm," ","");
        res=substr(id,6);

        if(confs[res]=="") confs[res]=1;

        # split main chain of multi-conf residues
        if(Bfac[id,confs[res]]=="" && ( atm ~ /^[CNO]$/ || atom=="  CA " )){
            for(c=1;c<=confs[res];++c){
                occ[id,c]=resocc[res,c];
                Bfac[id,c]=Bfac[id,1];
                XYZ[id,c]=XYZ[id,1];
            }
        }
        # normalize occupancy of protein atoms
        occsum=0;
        for(c=1;c<=confs[res];++c){
            occsum+=resocc[res,c];
        }
        if(occsum+0==0) occsum=1;

        for(c=1;c<=seen[id];++c){
            o=resocc[res,c]/occsum;
            if(o<=0)o=0.01;
            if(res ~ /ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/) {
                resocc[res,c]=o;
            }
            if(res ~ /HID|HIE|HIP|HIS|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/) {
                resocc[res,c]=o;
            }
        }
        for(c=1;c<=seen[id];++c){
            o=resocc[res,c];
            conf=" ";
            if(confs[res]>1 || (o+0<1) ) conf=substr(conflib,c,1);
            idc=atom conf res
            printf("ATOM%7d%s   %s%6.2f%s\n",++atnum,idc,XYZ[id,c],o,Bfac[id,c]);


            idc=atom conf res;
            if(anisc[idc]=="") anisc[idc]=anis[id];
            if(anisc[idc]!=""){
                printf("ANISOU%5d%s%s\n",atnum,idc,anisc[idc]);
            }
        }
    }
    print "END"
}

