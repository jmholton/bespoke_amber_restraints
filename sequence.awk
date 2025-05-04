#! /bin/awk -f
#
#   Process/identify protein sequences in a text/pdb file             -James Holton  2-19-25
#   as > 20 consecutive, aa letters
#
#   plus a few other goodies, such as monoisotopic mass, identifying 
#   chemically unstable sequences, and common cleavage sites (using chop=yes)
#
#
BEGIN {

    if(! minlength) minlength = 20

    # one-letter amino-acid code
    OLC["ALA"] = "A"
    OLC["CYS"] = "C"
    OLC["ASP"] = "D"
    OLC["GLU"] = "E"
    OLC["PHE"] = "F"
    OLC["GLY"] = "G"
    OLC["HIS"] = "H"
    OLC["ILE"] = "I"
    OLC["LYS"] = "K"
    OLC["LEU"] = "L"
    OLC["MET"] = "M"
    OLC["MSE"] = "M"
    OLC["ASN"] = "N"
    OLC["PRO"] = "P"
    OLC["GLN"] = "Q"
    OLC["ARG"] = "R"
    OLC["SER"] = "S"
    OLC["THR"] = "T"
    OLC["VAL"] = "V"
    OLC["TRP"] = "W"
    OLC["TYR"] = "Y"    
    

    # average mass
    H  =  1.007947
    C  = 12.0111  
    N  = 14.006747
    O  = 15.999943
    P  = 30.973762
    S  = 32.0666  
    Cl = 35.45279 
    Se = 78.963   

    # residue masses -NH-CH(R)-CO-
    aMass["G"] =  2*H+2*C+N+O+ H
    aMass["A"] =  2*H+2*C+N+O+ C+3*H
    aMass["V"] =  2*H+2*C+N+O+ C+H + 2*(C+3*H)
    aMass["I"] =  2*H+2*C+N+O+ C+H + C+2*H + 2*(C+3*H)
    aMass["L"] =  2*H+2*C+N+O+ C+2*H + C+H + 2*(C+3*H)
    aMass["D"] =  2*H+2*C+N+O+ C+2*H + C + 2*O  +H
    aMass["N"] =  2*H+2*C+N+O+ C+2*H + C + O+N+2*H
    aMass["E"] =  2*H+2*C+N+O+ C+2*H + C+2*H + C + 2*O  +H
    aMass["Q"] =  2*H+2*C+N+O+ C+2*H + C+2*H + C + O+N+2*H
    aMass["R"] =  2*H+2*C+N+O+ C+2*H + C+2*H + C+2*H + N + C + 2*(N+ 2*H)
    aMass["K"] =  2*H+2*C+N+O+ C+2*H + C+2*H + C+2*H + C+2*H + N+2*H
    aMass["M"] =  2*H+2*C+N+O+ C+2*H + C+2*H + S     + C+3*H
    aMass["C"] =  2*H+2*C+N+O+ C+2*H + S+H
    aMass["S"] =  2*H+2*C+N+O+ C+2*H + O+H
    aMass["T"] =  2*H+2*C+N+O+ C+H   + O+H   + C+3*H
    aMass["H"] =  2*H+2*C+N+O+ C+2*H + C     + C+H + N+H + C+H + N 
    aMass["W"] =  2*H+2*C+N+O+ C+2*H + 3*C + 5*(C+H) + N+H
    aMass["F"] =  2*H+2*C+N+O+ C+2*H + C +   5*(C+H)
    aMass["Y"] =  2*H+2*C+N+O+ C+2*H + C +   5*(C+H) + O
    aMass["P"] =  1*H+2*C+N+O+ 3*(C+2*H)

    # monoisotopic mass (most abundant isotope)
    H  = 1.007825 
    C  = 12.00000 
    N  = 14.003074
    O  = 15.994915
    P  = 30.973762
    S  = 31.972070
    Cl = 34.968852
    Se = 79.916520

    # residue masses -NH-CH(R)-CO-
    iMass["G"] =  2*H+2*C+N+O+ H
    iMass["A"] =  2*H+2*C+N+O+ C+3*H
    iMass["V"] =  2*H+2*C+N+O+ C+H + 2*(C+3*H)
    iMass["I"] =  2*H+2*C+N+O+ C+H + C+2*H + 2*(C+3*H)
    iMass["L"] =  2*H+2*C+N+O+ C+2*H + C+H + 2*(C+3*H)
    iMass["D"] =  2*H+2*C+N+O+ C+2*H + C + 2*O  +H
    iMass["N"] =  2*H+2*C+N+O+ C+2*H + C + O+N+2*H
    iMass["E"] =  2*H+2*C+N+O+ C+2*H + C+2*H + C + 2*O  +H
    iMass["Q"] =  2*H+2*C+N+O+ C+2*H + C+2*H + C + O+N+2*H
    iMass["R"] =  2*H+2*C+N+O+ C+2*H + C+2*H + C+2*H + N + C + 2*(N+ 2*H)
    iMass["K"] =  2*H+2*C+N+O+ C+2*H + C+2*H + C+2*H + C+2*H + N+2*H
    iMass["M"] =  2*H+2*C+N+O+ C+2*H + C+2*H + S     + C+3*H
    iMass["C"] =  2*H+2*C+N+O+ C+2*H + S+H
    iMass["S"] =  2*H+2*C+N+O+ C+2*H + O+H
    iMass["T"] =  2*H+2*C+N+O+ C+H   + O+H   + C+3*H
    iMass["H"] =  2*H+2*C+N+O+ C+2*H + C     + C+H + N+H + C+H + N 
    iMass["W"] =  2*H+2*C+N+O+ C+2*H + 3*C + 5*(C+H) + N+H
    iMass["F"] =  2*H+2*C+N+O+ C+2*H + C +   5*(C+H)
    iMass["Y"] =  2*H+2*C+N+O+ C+2*H + C +   5*(C+H) + O
    iMass["P"] =  1*H+2*C+N+O+ 3*(C+2*H)
    iMass["formyl"]=  C + O
    iMass["p"] =  P + 3*O
}

# read standard SEQRES cards from a PDB file
/SEQRES/{
    if(! seqres) seq = ""
    seqres = 1;
    for(i=4;i<=NF;++i)
    {
	seq = seq OLC[$i]
	if($i !~ /^[A-Z][A-Z].$/) continue
	if(OLC[$i]=="") seq = seq "X"
    }
}


# don't do other kinds of search in a PDB file
seqres {next}

# read sequence of a PDB
/^ATOM / && substr($0,12,5) == "  CA " {
    if(! pdb) seq = ""
    pdb = 1

    Restype = substr($0, 18, 3)
    Segid   = substr($0, 22, 1)    	# O/Brookhaven-style segment ID
    Resnum  = substr($0, 23, 4)+0
    
    if(seen[Segid Resnum]) next

    # check for breaks
    if((Segid != lastSegid)||(nextResnum != Resnum && Resnum != lastResnum)) {
	    # break in chain
	    seq = seq " "
	    lastSegid = Segid
    }
    lastResnum = Resnum
    nextResnum = Resnum +1

    # translate three-letter code to one letter
    seq = seq OLC[Restype]
    if(OLC[Restype]=="") seq = seq "X"

    ++seen[Segid Resnum]
}
/^TER/{seq = seq " "}

# don't do other kinds of search in a PDB file
pdb {next}


# recognize ENTREZ files
$1~/[0-9]/ && $1+0==$1 && $2 !~ /[^a-z]/ && $NF !~ /[^a-z]/ && NF<10{
    if(! entrez) seq=""
    entrez = 1
    $0=$2 $3 $4 $5 $6 $7 $8 $9 $10 $11
}

(length($0) > 9 || seq != "") && ! /^>/ {
#{
    # remove leading spaces
    line = ""
    for(i=1;i<NF;++i) line = line $i " " 
    line = line $NF    

    # scan for aa letters
    for(i=1;i<=length(line);++i)
    {
	c = toupper(substr(line, i, 1));
	# ignore these characters
	if(c == "\"") c = ""
	if(c == "\t") c = ""
	if(c == " ")  c = " "
	if(c == " ")  c = ""

	if(c !~ /[A,C-I,K-N,P-T,V-Y]/ && c != "")
#	if(c !~ /[A,C-I,K-N,P-T,V-Y]/)
	{
	    c = " "
	}
	seq = seq c 
    }
}

# blank lines terminate a sequence
NF==0 {
    seq = seq " "
    seqres = pdb = entrez = 0
}

END{

if(debug) print seq

# break up strings of protein letters into "words"
num = split(seq, sequence)

for(n=1;n<=num;++n)
    if(length(sequence[n]) >= minlength || seqres)
    {
	# look for all the horrible things that can happen to the peptide
	acid = ""
	base = ""
	race = ""
	pyroQ = ""
	CNBr = ""

	factorXa = ""
	chymotrypsin = ""
	endoproteinaseDN = ""
	endoproteinaseKC = ""
	thrombin = ""
	trypsin = ""
	pepsin = ""
	V8 = ""

	# weigh this chain
	weight = Met = His = Cys = A280 = "";
	weight = O + 3*H;
	mass   = O + 3*H;
	for(i=1;i<=length(sequence[n]);++i)
	{
	    c = substr(sequence[n], i, 1);
	    
	    # weigh this chain
	    weight += aMass[c];
	    mass   += aMass[c];
	
	    # count potentially derivitized residues
	    if(c == "M") ++Met;
	    if(c == "C") ++Cys;
	    if(c == "H") ++His;
	    
	    # add up (denatured) extinction coefficient
	    if(c == "W") A280 += 5600
	    if(c == "Y") A280 += 1400
	    if(c == "F") A280 += 197
	    
	    # chemical instabilities  (add up single-cleavage MWs)
	    if(c == "M" )                       CNBr = CNBr " " mass
	    if(substr(sequence[n],i,2) == "DP") acid = acid " " mass
	    if(substr(sequence[n],i,2) == "NG") base = base " " mass
	    if(substr(sequence[n],i,2) == "NG") race = race ", " i
	    if(substr(sequence[n],1,1) == "Q") pyroQ = 1
	    if((substr(sequence[n],1,1) == "M")&&(substr(sequence[n],2,1) == "Q")) pyroQ = 2

	    # proteolytic recognition sites (add up single-cleavage MWs)?
	    if(substr(sequence[n],i-3,4) == "IEGR") factorXa = factorXa  " " mass
	    if(substr(sequence[n],i+1,1) == "D") endoDN = endoDN " " mass
	    if(c ~ /[Y,F,W]/)     chymotrypsin = chymotrypsin    " " mass
#	    if(c ~ /[L,M,A,N,E]/) chymotrypsin = chymotrypsin    " " mass "*"
	    if(c ~ /[K]/)               endoKC = endoKC          " " mass
	    if(c ~ /[R]/)             thrombin = thrombin        " " mass
	    if(c ~ /[R,K]/)            trypsin = trypsin         " " mass
	    if(c ~ /[F,L]/)             pepsin = pepsin          " " mass
#	    if(c ~ /[Y,W,I,M]/)         pepsin = pepsin          " " mass "*"
	    if(c ~ /[E]/)                   V8 = V8              " " mass
	    if(c ~ /[E]/)                  V82 = V82             " " mass
	    if(c ~ /[D]/)                  V82 = V82             " " mass
	}
	
	# finish off cleavages
	acid = acid " " mass
	base = base " " mass
	CNBr = CNBr " " mass

	factorXa = factorXa         " " mass
	chymotrypsin = chymotrypsin " " mass
	endoDN = endoDN             " " mass
	endoKC = endoKC             " " mass
	thrombin = thrombin         " " mass
	trypsin = trypsin           " " mass
	pepsin = pepsin             " " mass
	V8 = V8                     " " mass
	
	
	# we have found an acceptable protein sequence
	print mass " Da chain: "
	l=length(sequence[n])
	while(length(sequence[n]) > 0)
	{
	    # actually print out sequence here
	    print substr(sequence[n], 1, 80)
	    sequence[n] = substr(sequence[n], 81)
	}
	print ""
	print l "aa"
	print Met+0 "met"
	print Cys+0 "cys"
	print His+0 "his"
	print ""
	printf "denatured A(280nm) = %.4f*l*c (c in g/L)\n", A280/weight
	printf "    SeMET MAD Rano = %.3f%%\n", 100*(Met*8^2)/(7^2 * (weight/14))
	print ""
	
	f=split(acid base, Split)
	if((f>2) || pyroE != "")
	{
	    print "Chemical Instabilities: "
	}
	f=split(acid, Split)
	if(f>1)
	{
	    printf "acid (D*P):                       "
	    for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	    print ""
	}
	f=split(base, Split)
	if(f>1)
	{
	    printf "base (N*G):                       "
	    for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	    print ""
	    print "racemization hazard at" substr(race,2)
	}
	if(pyroQ)
	{
	    print "residue " pyroQ " could form an N-cyclized glutamine "
	}

	if(chop)
	{
	    print ""
	    
	    f=split(CNBr, Split)
	    if(f>1)
	    {
		printf "CNBr (M*):                        "
		for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
		print ""
	    }
	    print ""
	    print "Common proteases: "
	    f=split(factorXa, Split)
	    if(f>1)
	    {
	        printf "factorXa (IEGR*):                 "
	        for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	        print ""
	    }
	    f=split(thrombin, Split)
	    if(f>1)
	    {
	        printf "thrombin (R*):                    "
	        for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	        print ""
	    }
	    f=split(trypsin, Split)
	    if(f>1)
	    {
	        printf "trypsin (R*, K*):                 "
	        for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	        print ""
	    }
	    f=split(endoKC, Split)
	    if(f>1)
	    {
	        printf "endoproteinase Lys-C (K*):        "
	        for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	        print ""
	    }
	    f=split(endoDN, Split)
	    if(f>1)
	    {
	        printf "endoproteinase Asp-N (*D):        "
	        for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	        print ""
	    }
	    f=split(chymotrypsin, Split)
	    if(f>1)
	    {
	        printf "chymotrypsin (W*,Y*,F*, +others): "
	        for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	        print ""
	    }
	    f=split(pepsin, Split)
	    if(f>1)
	    {
	        printf "pepsin (F*, L*, +others):         "
	        for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	        print ""
	    }
	    f=split(V8, Split)
	    if(f>1)
	    {
	        printf "V8 protease (E*):                 "
	        for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	        print ""
	    }
	    f=split(V82, Split)
	    if(f>1)
	    {
	        printf "V8 protease (E*,D*):              "
	        for(i=1;i<=f;++i) printf Split[i] - Split[i-1] " "
	        print ""
	    }
	}
    }
}

