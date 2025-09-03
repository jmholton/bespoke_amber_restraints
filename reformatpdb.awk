#! /usr/bin/awk -f
#
#
# Re-writes PDB into something more canonnical
# 
BEGIN {

# memorize ordering of the alphabet
abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
for(i=1;i<=length(abc);++i) {
  c=substr(abc,i,1);
  alphabet[i]=c;
  alphabet[c]=i;
}
alphabet[0]=" "; alphabet[" "]=0;

# canonical ordering of atom types in PDB
align[1]  = "N"  ; align[2]  = "A " ; align[3]  = "C"  ; align[4]  = "O"  ;
align[5]  = "B " ; align[6]  = "G " ; align[7]  = "G1" ; align[8]  = "G2" ; 
align[9]  = "G3" ; align[10] = "D " ; align[11] = "D1" ; align[12] = "D2" ; 
align[13] = "D3" ; align[14] = "E " ; align[15] = "E1" ; align[16] = "E2" ;
align[17] = "E3" ; align[18] = "Z " ; align[19] = "Z1" ; align[20] = "Z2" ;
align[21] = "Z3" ; align[22] = "H " ; align[23] = "H1" ; align[24] = "H2" ;
align[25] = "H3" ;

}

/^ATOM/ || /^HETATM/ {

    if(debug) print tolower($0)

#######################################################################################
    electrons = substr($0, 67,6)        # number of electrons in this atom (not always there)
    XPLORSegid = substr($0, 73, 4)      # XPLOR-style segment ID
    split(XPLORSegid, a)                # (remove spaces)
    XPLORSegid = a[1];
    Element = substr($0, 67)            # sometimes element is given here
    gsub("[ 0-9]","",Element)

    Atomnum= substr($0,  7, 5)+0        # atom number
    Element= substr($0, 13, 2);         # actual element number
    Greek= substr($0, 15, 2);           # "remoteness" number of this atom (i.e. "A" for C-alpha)
    split(Element Greek, a)             # (remove spaces)
    Atom   = a[1];                      # store whole atom name
    Conf   = substr($0, 17, 1)          # conformer letter
    Restyp = substr($0, 18, 3)          # residue name
    Segid  = substr($0, 22, 1)          # O/Brookhaven-style segment ID
    Resnum = substr($0, 23, 4)          # residue number
    Insert = substr($0, 27, 1)          # insertion code
    X      = substr($0, 31, 8)+0        # coordinates
    Y      = substr($0, 39, 8)+0
    Z      = substr($0, 47, 8)+0
    Occ    = substr($0, 55, 6)+0        # occupancy
    Bfac   = substr($0, 61, 6)+0        # B-factor
#   rest   = substr($0, 67)             # rest of the line after B-factor?
    ATOM   = toupper(substr($0, 1, 6))  # store given atom name
    ID     = substr($0,12,15)

#######################################################################################
#   correct for alternate formatting

    if((Segid == " ") && (substr(XPLORSegid,1,1) ~ /[A-Z]/))
    {
        Segid = substr(XPLORSegid,1,1)
    }
    if(Resnum ~ /[A-Z]/)
    {
        # incorrect residue numbers: A14, etc.
        Segid = substr(Resnum,match(Resnum,"[A-Z]"),1);
        Resnum = substr(Resnum,match(Resnum,"[A-Z]")+1);
    }
    Resnum+=0;

    # fix X-plor's non-standard MET
    if(Restyp == "MSE")
    {
#        Restyp = "MET"
        if((Greek == "E ")&&(Element == " S"))
        {
            Element = "SE"
            Greek = "  "
        }
    }

    Ee = Element
    if(Ee ~ /^H/ && Greek ~ /[1-9][1-9]/) Ee = " H"

    if(electrons+0 == 0)
    {
        if(Ee == " C") electrons = 6
        if(Ee == " O") electrons = 8
        if(Ee == " N") electrons = 7
        if(Ee == " H") electrons = 1
        if(Ee == " S") electrons = 16
        if(Ee == "SE") electrons = 34
    }

#######################################################################################
# user-directed globlal changes
    if(BFAC != "") 
    {
        if(BFAC !~ /^[+-]/) Bfac = 0
        Bfac += BFAC
    }
    if(OCC != "")
    {
        if(OCC !~ /^[+-]/) Occ = 0
        Occ += OCC
    }
    if(RESNUM != "")
    {
        Resnum = RESNUM
    }
    if(CONF != "") Conf=CONF;
    if(CHAIN != "") Segid = CHAIN
    if(renumber)
    {
        ++ATOMNUM
        Atomnum = ATOMNUM
    }
    # increment/decrement Segids?
    if(map[Segid]=="") map[Segid]=Segid
    ID = substr(ID,1,10) map[Segid] substr(ID,12)
    seen[ID]=seen[ID]+1
prev=Segid
    while((seen[ID]>1)&&(Segid != "z")&&(1))
    {
        map[Segid] = alphabet[alphabet[map[Segid]]+1]
        print "REMARK:", Segid,"->",map[Segid]
        ID = substr(ID,1,10) map[Segid] substr(ID,12)
        seen[ID]=seen[ID]+1
    }
    Segid=map[Segid]

#######################################################################################
    if(Resnum != lastResnum || Insert != lastInsert)
    {
        for(i=1;i<30;++i)
        {
            if((residue[align[i]] != "")&&(reorder))
            {
                print residue[align[i]]
            }
        }
        # clear for next time
        for(x in residue) residue[x] = ""
    }
    lastResnum = Resnum
    lastInsert = Insert
    
    order = Greek
    if(order == "  ") order = Atom
    residue[order] = sprintf("%6s%5d %2s%-2s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f %5.2f%6.2f%4s  %-4s%2s",\
        ATOM, Atomnum,Element,Greek,Conf,Restyp,Segid,Resnum,Insert,X,Y,Z,Occ,Bfac,"",XPLORSegid,Ee);
    
#######################################################################################
    # default, non-reordering mode
    if(!reorder)
    {
        printf("%6s%5d %2s%-2s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f %5.2f%6.2f%4s  %-4s%2s\n",\
       ATOM, Atomnum,Element,Greek,Conf,Restyp,Segid,Resnum,Insert,X,Y,Z,Occ,Bfac,"",XPLORSegid,Ee);        
    }
}

! /^ATOM/ && ! /^HETATM/

END{
        for(i=1;i<30;++i)
        {
            if((residue[align[i]] != "")&&(reorder))
            {
                print residue[align[i]]
            }
        }    
}
