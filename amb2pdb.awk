#! /bin/awk -f
#
#
#
BEGIN{nterm=1}

/^TER/{print "TER";next}
! /^ATOM|HETAT/{print;next}

{c=65}
/ACT|AMM/{c=65}
/WAT|HOH/{c=83}

{
    typ=substr($0,18,3);
    protein=0;
    if(typ~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL/){++protein};
    if(typ~/HID|HIE|HIS|ILE|LEU|LYS|PHE|PRO|SER|THR|TRP|TYR/){++protein};

    gsub("HID","HIS");
    gsub("HIE","HIS");
    gsub("HIP","HIS");
    gsub("ASH","ASP");
    gsub("GLH","GLU");
    gsub("CYX","CYS");
    gsub("WAT","HOH");
    if(nterm && protein && /  H1 / ){
        # only do this for N-terminal NH3
        gsub("  H1 ","  H  ");
        nterm=0;
    }

    gsub("HB1 ACT","H1  ACY");
    gsub("HB2 ACT","H2  ACY");
    gsub("HB3 ACT","H3  ACY");
    gsub("CA  ACT","C   ACY");
    gsub("CB  ACT","CH3 ACY");
    gsub("OA1 ACT","O   ACY");
    gsub("OA2 ACT","OXT ACY");

    gsub("H1  AMM","HN1 NH4");
    gsub("H2  AMM","HN2 NH4");
    gsub("H3  AMM","HN3 NH4");
    gsub("H4  AMM","HN4 NH4");
    gsub("AMM","NH4");
    atom=substr($0,13,4)
    resnum=substr($0,23,4)+0
    id = atom" "c" "resnum
}

seen[id]{
    loop1=0
    while(seen[id]){
        ++c;
        if(c==91) c=97;
        if(c==123 && ! loop1) {c=66;++loop1};
        if(c==123) {print "REMARK unavoidable duplicate";break};
        id = atom" "c" "resnum;
    }
}

id != last{
    ++seen[last]
}

{last=id}

# next residue is an N terminus
/OXT/{nterm=1}

{printf("%s%c%s\n",substr($0,1,21),c,substr($0,23))}

