#! /bin/awk -f
#
#	compute best least-squares linear fit to data
#
#
# slope  =  (n*<x y> - <x><y>) / ( n*<x^2> - <x>^2 )
# intercept = (<x^2>*<y> - <x><x y>) / ( n*<x^2> - <x>^2 )
# 

NF==2{
    x=$1;y=$2;
    ++n;
    xv[n]=x;yv[n]=y;

    sumx+=x;
    sumy+=y;
    sumxy+=x*y;
    sumxx+=x*x;
    sumyy+=y*y;
}

END{
    denom = n*sumxx-sumx^2;
    if( denom == 0 || n < 2){
	print 0,0; exit;
    }
    slope     = (n*sumxy - sumx*sumy)/denom;
    intercept = (sumxx*sumy - sumx*sumxy)/denom;

    if(! printout){
        printf("%.12g %.12g\n",slope,intercept)
    }
    else
    {
	for(i=1;i<=n;++i){
	    printf("%s %s %.12g\n",xv[i],yv[i],slope*xv[i]+intercept);
	}
    }
}
