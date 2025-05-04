#! /bin/awk -f
#
#  jiffy script for re-numbering large PDB files        -James Holton 1-31-24
#  code adopted from https://cci.lbl.gov/hybrid_36/
#
#
BEGIN{
   if(! ordinal) ordinal=0
}

! /^ATOM|HETAT/{print;next}


# parse the atom line
{
    if(debug) print "DEBUG" substr($0,6);
    pre=substr($0,1,22);
    post=substr($0,31);
    prevnum=resnum;
    resnum=substr($0,23,8);
    oresnum=resnum;
    if(ordinal){
      if(resnum!=prevnum){
        prevnum=resnum;
        ++rescounter;
      }
      oresnum=rescounter;
    }
    if(resnum~/[A-Za-z]/) {
      gsub(" ","",resnum);
      oresnum=hy36decode(resnum);
    }
    code=hy36encode(oresnum+0);
    printf("%s%4s    %s       %d\n",pre,code,post,oresnum);
}

# these seem to work
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

function hy36decode_gemmi( code )
{
  if(! match(code,/[A-Za-z]/)) return code;
  udigits="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  ldigits=tolower(udigits);
  n=length(code);
  value=0;
  ucode=toupper(code);
  for(i=1;i<=n;++i) {
    c = substr(ucode,i,1);
    value+=(index(udigits,c)-1)*(36**(n-i));
  }
  return value - 10*36**(n-1) + 10**n;
}

function hy36encode_gemmi( value )
{
  if(value < 10000) return value;

  n=4;code="";
  value = value + 10*36**3 - 10**n;
  udigits="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  ldigits=tolower(udigits);
  for(i=1;i<=n;++i) {
    q = int( value / 36**(n-i) )
    c = substr(udigits,q+1,1);
    code = code c;
    value = value - q*36**(n-i)
  }
  return code;
}

