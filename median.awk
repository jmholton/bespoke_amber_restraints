#! /bin/awk -f
#
#
#
#

BEGIN{min=1e99;max=-1e99}

min>$1{min=$1}
max<$1{max=$1}

{++n;v[n]=$1}

END{
    range=max-min;
    if(range==0) range=1;
    for(i=1;i<=n;++i){
	sortme[i]=sprintf("%.20f",(v[i]-min)/(range));
	idx[sortme[i]]=i;
    }
    asort(sortme);
    i=int(n/2.0)+1;
    if(n%2==0){
	median = (v[idx[sortme[i]]] + v[idx[sortme[i-1]]])/2
    }
    else
    {
        median = v[idx[sortme[i]]];
    }

    min=1e99;max=1e-99
    for(i=1;i<=n;++i){
	d[i] = sqrt((v[i]-median)^2);
	if(min>d[i])min=d[i];
	if(max<d[i])max=d[i];
    }
    range=max-min;
    if(range==0) range=1;
    for(i=1;i<=n;++i){
	sortme[i]=sprintf("%.20f",(d[i]-min)/(range));
	idx[sortme[i]]=i;
    }
    asort(sortme);
    i=int(n/2.0)+1;
    if(n%2==0){
	mad = (d[idx[sortme[i]]] + d[idx[sortme[i-1]]])/2
    }
    else
    {
        mad = d[idx[sortme[i]]];
    }

    print median,"+/-",mad;
}
