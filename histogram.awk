#! /bin/awk -f
#
#	make a histogram of the input data
#
#
BEGIN{
    if(bin) bs=bin;
    if(bw) bs=bw;
    if(bs<=0) bs=1;
    if(logscale) { }
}

{   
    if(! logscale){
	++n;
        bin=sprintf("%.0f",$1/bs)+0;
	++count[bin];
	binlist[bin]=sprintf("%+020d",bin);
    }
    if(logscale && $1+0>0){
	++n;
        bin=sprintf("%.0f",(log($1)/bs))+0;
	++count[bin];
	binlist[bin]=sprintf("%+020d",bin);
    }
} 

END{
    bins=asort(binlist);
    for(b=1;b<=bins;++b) {
#    for(bin in count) {
	bin=binlist[b]+0;
	if(logscale && bin*bs>-300 && bin*bs<300){
	    print exp(bin*bs),count[bin]/n/bs/exp(bin*bs),count[bin];
	}
	else
	{
	    print bin*bs,count[bin]/n/bs,count[bin];
	}
    }
}

