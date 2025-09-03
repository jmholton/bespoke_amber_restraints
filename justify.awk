#! /bin/awk -f
#
#  clean up messy lists
#
#
BEGIN{\
  if(! justify) justify="left";
  j="";
  if(justify=="left") j="-";
}

{
  # count existing spaces
  w=0;p=""
  for(i=1;i<=length($0);++i) {
    c=substr($0,i,1);
    if(c==" ")++sp[w];
    if(c=="\t"){sp[w]+=4;c=" "};
    if(c!=" " && p==" "){++w};
    p=c;
  }
  for(w=1;w<=NF;++w){
    v[NR,w]=$w
    ln=length($w)
    if(l[w]<ln) l[w]=ln;
  }
}

nf<NF{nf=NF}

END{
  for(w=1;w<=nf;++w){
    sp[w]=int(sp[w]/NR-1);
    if(sp[w]<0)sp[w]=0;
  }
  # look for columns that are all same
  for(w=1;w<=nf;++w){
   first=v[1,w]
   for(nr=2;nr<=NR;++nr){
    if(v[nr,w]!=first && v[nr,w]!=""){++diff[w];break}
   }
  }
  for(nr=1;nr<=NR;++nr){
    for(w=1;w<=nf;++w){
      if(dropallsame && ! diff[w]) continue
      printf("%" j (l[w]+sp[w]) "s ",v[nr,w]);
    }
    print ""
  }
}

