#! /bin/awk -f
# Nudge any water whose O lies within <cut> A of a protein/ligand HEAVY atom
# radially outward to <push> A from that atom.
#
# Fixes the hurry / un-buildout case where tleap completes a long sidechain
# (ARG/LYS) with an atom that lands on top of a bulk water: a ~0-distance overlap
# whose 1/r^12 van der Waals term overflows to NaN and makes the GPU MD crash
# ("illegal memory access") on the first heating step.  Waters are fungible, but
# rather than delete them (which would desync the topology from orignames /
# ref.crd / refpoints) we just move the offending water out - only COORDINATES
# change, so atom count, order and topology are untouched.  A legitimate
# crystallographic water H-bonds at >=2.4 A, so <1.8 A from a solute heavy atom is
# always a clash, never a real contact.  Strict no-op when nothing clashes.
#
#   declash_waters.awk [-v cut=1.8] [-v push=3.0] model.pdb   (nudged pdb -> stdout)
BEGIN{ if(cut=="")cut=1.8; if(push=="")push=3.0 }
{
  L[NR]=$0
  if($0 ~ /^ATOM|^HETATM/){
    na++; ln[na]=NR
    resn=substr($0,18,3); gsub(/ /,"",resn)
    atn =substr($0,13,4); gsub(/ /,"",atn)
    X[NR]=substr($0,31,8)+0; Y[NR]=substr($0,39,8)+0; Z[NR]=substr($0,47,8)+0
    isw=(resn=="WAT"||resn=="HOH"||resn=="SOL"||resn=="T3P"||resn=="OPC"||resn=="TP3")
    ISW[NR]=isw
    isH=(substr(atn,1,1)=="H")
    if(isw==0 && isH==0){ nh++; HX[nh]=X[NR];HY[nh]=Y[NR];HZ[nh]=Z[NR] }  # solute heavy
    if(isw==1 && atn=="O"){ WO[NR]=1 }
  }
}
END{
  c2=cut*cut
  for(k=1;k<=na;k++){
    i=ln[k]; if((i in WO)==0) continue
    best=1e18; bx=by=bz=0
    for(h=1;h<=nh;h++){
      dx=X[i]-HX[h]; if(dx>cut||dx<-cut)continue
      dy=Y[i]-HY[h]; if(dy>cut||dy<-cut)continue
      dz=Z[i]-HZ[h]; if(dz>cut||dz<-cut)continue
      d=dx*dx+dy*dy+dz*dz; if(d<best){best=d;bx=HX[h];by=HY[h];bz=HZ[h]}
    }
    if(best<c2){
      bd=sqrt(best); if(bd<1e-4){bx-=1;bd=1}
      ux=(X[i]-bx)/bd; uy=(Y[i]-by)/bd; uz=(Z[i]-bz)/bd
      ddx=(bx+ux*push)-X[i]; ddy=(by+uy*push)-Y[i]; ddz=(bz+uz*push)-Z[i]
      MV[i]=1; DX[i]=ddx; DY[i]=ddy; DZ[i]=ddz          # the O
      m=k+1
      while(m<=na && ISW[ln[m]]==1 && (ln[m] in WO)==0){  # its H's / EPW, up to next water O
        j=ln[m]; MV[j]=1; DX[j]=ddx; DY[j]=ddy; DZ[j]=ddz; m++
      }
      nmv++
    }
  }
  for(r=1;r<=NR;r++){
    if(r in MV){
      printf("%s%8.3f%8.3f%8.3f%s\n", substr(L[r],1,30), X[r]+DX[r], Y[r]+DY[r], Z[r]+DZ[r], substr(L[r],55))
    } else print L[r]
  }
  printf("declash_waters: nudged %d clashing water(s) (cut=%.1f push=%.1f)\n", nmv+0, cut, push) > "/dev/stderr"
}
