#! /bin/awk -f
#
#
# takes the average of a list of numbers, then prints out another list
# of deviates and SIGMAs
# 
BEGIN {

if (!start) start = 1;
if (!stop) stop = 1000;
if (col) start=stop=col;

maxcol=0;
maxnumber=0;
}

{
    for(col=start; (col<=NF)&&(col <= stop); ++col)
    {
	# count number of items in this column
	++number[col];
	# store them in memory
        line[col, number[col]] = $col;
	# compute the average
        sum[col] += $col;	
    }
}

END{
    # find number of columns
    for(col in sum)
    {
	if(maxcol+0 < col+0) maxcol = col;
        if(maxnumber+0 < number[col]+0) maxnumber = number[col];
    }
    
    for(col=start; col <= maxcol; ++col)
    {
	if(number[col]+0<=0) break
	avg[col] = sum[col]/number[col]
	sum[col] = 0;

        for (i = 1; i <= number[col]; ++i)
	{
	    dev= line[col, i]-avg[col];
	    sum[col] += dev*dev;
	}
	rms[col] = sqrt(sum[col]/number[col]);

	print avg[col] " +/- " rms[col];
    }

    # reject outliers
    rejects=9999
    while(reject && rejects>0)
    {
	rejects=0
	for(col=start; col <= maxcol; ++col)
	{
	    count[col]=0
	    sum[col]=0
	    for (i = 1; i <= number[col]; ++i)
	    {
		if(! rejected[col, i])
		{
		    # this datum is still in play
		    if(sqrt((line[col, i] - avg[col])^2) > reject*rms[col])
		    {
			# this is outside the specified reject limit
			print "rejecting column",col,"row",i,"("line[col,i]")"
			++rejects;
			rejected[col, i]=1;
		    }
		}
		
		if(! rejected[col, i])
		{
		    # this datum is still valid
		    ++count[col];
		    sum[col] += line[col, i];
		}
	    }
	    avg[col] = sum[col]/count[col];
	    sum[col]=0;
	    for (i = 1; i <= number[col]; ++i)
	    {
		if(! rejected[col, i])
		{
		    dev= line[col, i]-avg[col];
		    sum[col] += dev*dev;
		}
	    }
	    rms[col] = sqrt(sum[col]/count[col]);
	}
    }

    # now print this out
    for(col=start; col <= maxcol; ++col)
    {
	if(number[col]+0>0) print avg[col] " +/- " rms[col];
    }

    if(xlog)
    {
	# print out as xloggraph file
	print "";
	print " $TABLE : Stats:"
	printf " $GRAPHS:Deviates:A:1"
	for(col=start; col <=maxcol; ++col)
	{
	    printf ", "2*(col-start+1)
	}
	print ":"
	printf         ":Sigmas:A:1"
	for(col=start; col <=maxcol; ++col)
	{
	    printf ", "2*(col-start+1)+1
	}
	print ": $$"
	printf "line "
	for(col=start; col <=maxcol; ++col)
	{
	    printf " delta"col " sigma"col
	}	
	print " $$"
	print " $$"
        for (i = 1; i <= maxnumber; ++i)
	{
	    printf i" ";
	    for(col=start; col <= maxcol; ++col)
	    {
		if(!rms[col])
		{
		    rms[col]=10000000000000;
		}
		printf line[col, i]+0 " ";
		printf (line[col, i]+0)/rms[col]" ";
	    }
	    print "";
	}
	print "$$\n\n";
    }
}
