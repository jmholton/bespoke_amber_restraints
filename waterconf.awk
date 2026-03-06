#! /bin/awk -f
#
#
#  re-organize conformers in a PDB file to avoid clashes   - James Holton 9-1-25
#  also allow pre-cached distances vs atom number as "DIST i j dist" for symmetry mates
#
#
BEGIN{
    if(! conflib) conflib = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_.,:;=<>+-%&*$!(){}[]^#|~?"
    if(! tooclose) tooclose = 0.0001
    if(! clashdist) clashdist = 2.4
    if(! renumber) renumber = 0
}

/^ANISO/{next}

/^DIST/{i=$2;j=$3;distcache[i,j]=$4;next}

! /^ATOM|^HETAT/{print;next}

{
    ++n;
    conf[n]=substr($0,17,1);
    X[n] = substr($0,31,8)+0;
    Y[n] = substr($0,39,8)+0;
    Z[n] = substr($0,47,8)+0;
    occ[n] = substr($0,55,6)+0;
    cn[n] = substr($0,22,6);
    typ=substr($0,18,3);
}

typ !~ /HOH|CL |NA /{
    print
#    next
}

typ ~ /HOH|CL |NA /{
    delete involved;
    delete worstdist;
    for(i=1;i<n;++i){
        if(! distcache[i,n]) {
           distcache[i,n]=sqrt((X[i]-X[n])^2+(Y[i]-Y[n])^2+(Z[i]-Z[n])^2);
#           print "DEBUG: cacheing",i,n,distcache[i,n]
        }
        dist=distcache[i,n];
        if(dist<worstdist[conf[i]] || worstdist[conf[i]]=="") {
            worstdist[conf[i]]=dist;
        }
        if(dist <= tooclose) next;
        if(dist < clashdist) {
            ++involved[conf[i]];
            partner=cn[i];
            print "REMARK: partner",conf[i],cn[i],"at",dist
        }
    }
    cnf = conf[n];
    j=0;
    while(involved[cnf] && cnf!="") {
        cnf = substr(conflib,++j,1);
    }
    if(cnf=="") {
        print "REMARK: warning: all confs clash with:",cn[n]
        bestdist=0;
        bestconf="";
        for(cnf in worstdist) {
           if(worstdist[cnf]>bestdist) {
               bestdist=worstdist[cnf];
               bestconf=cnf;
           }
        }
        print "REMARK: warning: all confs clash with",cn[n],"best of worst",bestconf,"at",bestdist;
        cnf=bestconf;
    }
    printf("%s%1s%s\n",substr($0,1,16),cnf,substr($0,18));
    conf[n]=cnf;
}



