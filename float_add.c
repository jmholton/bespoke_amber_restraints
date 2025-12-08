
/* add two binary "float" files together - report on resulting stats                           -James Holton  3-29-24

example:

gcc -O -O -o float_add float_add.c -lm
./float_add file1.map file2.map output.map

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#ifndef NAN
#define NAN strtod("NAN",NULL)
#endif


char *infile1name = "";
char *infile2name = "";
FILE *infile1 = NULL;
FILE *infile2 = NULL;
char *outfilename = "output.bin\0";
FILE *outfile = NULL;

float fmedian(unsigned int n, float arr[]);
float fmedian_with_rejection(unsigned int n, float arr[],float sigma,float *mad,int *final_n);
float fmedian_absolute_deviation(unsigned int n, float arr[], float median_value);
float fmean_with_rejection(unsigned int starting_points, float arr[], float sigma_cutoff, double *final_rmsd, int *final_n);

int main(int argc, char** argv)
{

    int n,i,j,k,pixels;
    float *outimage;
    float *inimage1;
    float *inimage2;
    unsigned short int *invalid_pixel;
    float *medbuffer;
    float median,mad;
    double sum_arej,sumsq_arej,sumd_arej,sumdsq_arej,avg_arej,rms_arej,rmsd_arej;
    double sumd3_arej,sumd4_arej,skewness_arej,kurtosis_arej,max_arej,min_arej;
    int valid_pixels,n_arej;
    float scale1=1.0,scale2=1.0,outscale=1.0,offset=0.0,outoffset=0.0;
    double float1, float2, float3, float4;
    double result,outdouble;
    float outfloat;
    int header=0;
    int normalize=0;
    int overflows=0;
    int underflows=0;
    int ignore_values=0;
    float ignore_value[70000];
    int reject_outliers=0;
    float outlier_sigma=6.0;
    int histogramize=0;
    unsigned long int *histogram;
    float binsize=1.0,epsilon=0.0;
    int bin,bins;
    double diff,deviate;
    double sum,sumsq,sumd,sumdsq,avg,rms,rmsd,max,min;
    double sumd3,sumd4,skewness,kurtosis;
    double sumx=0,sumy=0,sumxy=0,sumxx=0,sumyy=0;
    double avgx,avgy,avgxy,avgxx,avgyy,CC;
    double denom,lsq_scale,lsq_offset;


    /* check argument list */
    for(i=1; i<argc; ++i)
    {
        if(strlen(argv[i]) > 4)
        {
            if(strstr(argv[i]+strlen(argv[i])-4,".bin") || strstr(argv[i]+strlen(argv[i])-4,".map"))
            {
                printf("filename: %s ",argv[i]);
                if(infile1 == NULL){
                    infile1 = fopen(argv[i],"r");
                    if(infile1 == NULL){
                        perror(argv[i]);
                    }else{
                        infile1name = argv[i];
                        printf(" input1");
                    }
                }
                else
                {
                    if(infile2 == NULL){
                        infile2 = fopen(argv[i],"r");
                        if(infile2 == NULL){
                            perror(argv[i]);
                        }else{
                            infile2name = argv[i];
                            printf(" input2");
                        }
                    }
                    else
                    {
                        outfilename = argv[i];
                        printf(" output");
                    }
                }
                printf("\n");
            }
        }

        if(argv[i][0] == '-')
        {
            /* option specified */
            if(strstr(argv[i], "-scale1") && (argc >= (i+1)))
            {
                scale1 = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-scale2") && (argc >= (i+1)))
            {
                scale2 = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-outscale") && (argc >= (i+1)))
            {
                outscale = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-outfile") && (argc >= (i+1)))
            {
                outfilename = argv[i+1];
                ++i;
                continue;
            }
            if(strstr(argv[i], "-output") && (argc >= (i+1)))
            {
                outfilename = argv[i+1];
                ++i;
                continue;
            }
            if(strstr(argv[i], "-normalize"))
            {
                normalize = 1;
            }
            if(strstr(argv[i], "-histogram"))
            {
                histogramize = 1;
            }
            if(strstr(argv[i], "-bin") && (argc >= (i+1)))
            {
                binsize = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-reject"))
            {
                reject_outliers = 1;
            }
            if(strstr(argv[i], "-offset") && (argc >= (i+1)))
            {
                offset = atof(argv[i+1]);
                outoffset = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-ignore") && (argc >= (i+1)))
            {
                ++ignore_values;
                ignore_value[ignore_values] = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-epsilon") && (argc >= (i+1)))
            {
                epsilon = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-outoffset") && (argc >= (i+1)))
            {
                outoffset = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-header") && (argc >= (i+1)))
            {
                header = atof(argv[i+1]);
            }
        }
    }

    if(infile1 == NULL){
        printf("usage: float_add file1.bin file2.bin [outfile.bin] -scale1 1.0 -scale2 1.0 -offset 40 -header 512\n");
        printf("options:\n");\
        printf("\t-scale1 \tscale factor for first file\n");
        printf("\t-scale2 \tscale factor for second file\n");
        printf("\t-offset \tinteger offset to be subtracted from each file before applying scale\n");
        printf("\t-outoffset\tinteger offset to be added to the output file after all scaling\n");
        printf("\t-header \tnumber of bytes to ignore in header of each file\n");
        printf("\t-normalize \toutput relative difference instead of difference\n");
        printf("\t-histogram \tprint out a histogram of output values\n");
        printf("\t-binsize \tsize of the bins to use for the histogram\n");
        printf("\t-epsilon \tmaximum difference between values to be considered equal\n");
        printf("\t-ignore \tflag particular values as bad pixels, kept in output\n");
        printf("\t-reject \treject outliers from printed statistics\n");
        exit(9);
    }


    /* load first float-image */
    fseek(infile1,0,SEEK_END);
    n = ftell(infile1);
    fseek(infile1,0,SEEK_SET);
    inimage1 = calloc(n,1);
    inimage2 = calloc(n,1);
    outimage = calloc(n,1);
    invalid_pixel = calloc(n,1);
    fread(inimage1,n,1,infile1);
    fclose(infile1);
    if(infile2 != NULL)
    {
        fseek(infile2,0,SEEK_SET);
        fread(inimage2,n,1,infile2);
        fclose(infile2);
    }

    pixels = (n-header)/sizeof(float);
    medbuffer = calloc(2*pixels+2,sizeof(float));
    outfile = fopen(outfilename,"wb");
    if(outfile == NULL)
    {
        printf("ERROR: unable to open %s\n", outfilename);
        exit(9);
    }

    printf("header = %d bytes pixel zero offset = %g (%g in output)\n",header,offset,outoffset);

    for(k=1;k<=ignore_values;++k)
    {
        printf("ignoring value %g in both files\n",ignore_value[k]);
    }

    /* copy the header of the first file */
    for(i=0;i<=header/sizeof(float);++i)
    {
        outimage[i] = inimage1[i];
    }
    sum = sumsq = sumd = sumdsq = 0.0;
    min = NAN;
    max = NAN;
    valid_pixels = 0;
    for(j=0;j<pixels;++j)
    {
        i=j+header/sizeof(float);

        float1 = (double) inimage1[i];
        float2 = (double) inimage2[i];
        if( normalize )
        {
            /* report relative difference */
            avg = ( scale1*(float1-offset) + scale2*(float2-offset) ) / 2.0 ;
            diff = scale1*(float1-offset) - scale2*(float2-offset) ;
            result = outscale * (diff / avg);
        }
        else
        {
            result = outscale * (scale1*(float1-offset) + scale2*(float2-offset));
        }

        /* convert back to float */
        outfloat = outdouble = result+outoffset;
//        i = fpclassify(outfloat);

        if(isnan(outfloat))
        {
            ++invalid_pixel[i];
        }
#ifdef isnormal
        if(! isnormal(outfloat))
        {
#endif
            if(isinf(outfloat))
            {
                ++overflows;
            }
#ifdef isnormal
            if(! isfinite(outfloat))
            {
                ++underflows;
            }
#endif
#ifdef isnormal
        }
#endif
        /* assign to output array */
        outimage[i] = outfloat;

        /* skip any invalid values, propagate to output */
        for(k=1;k<=ignore_values;++k)
        {
            if(inimage1[i]==ignore_value[k]){
                outimage[i] = ignore_value[k];
                ++invalid_pixel[i];
                /* no need to check others */
                break;
            }
            if(inimage2[i]==ignore_value[k] && infile2 != NULL){
                outimage[i] = ignore_value[k];
                ++invalid_pixel[i];
                /* no need to check others */
                break;
            }
        }
        /* skip stats for this one */
        if(invalid_pixel[i]) continue;

        /* only use valid pixels in statistics */
        ++valid_pixels;
        /* don't do median stuff if we don't have to */
        medbuffer[valid_pixels]=result;

        /* basic stats */
        sum += result;
        sumsq += result*result;
        if(result>max || isnan(max)) max=result;
        if(result<min || isnan(min)) min=result;
        /* set up for CC calculation */
        sumx+=float1;
        sumy+=float2;
        sumxy+=float1*float2;
        sumxx+=float1*float1;
        sumyy+=float2*float2;
    }
    avg = sum/valid_pixels;
    rms = sqrt(sumsq/valid_pixels);

    denom = (sumxx*valid_pixels)-sumx*sumx;
    if( denom == 0 || valid_pixels < 2 )
    {
        lsq_scale = 0.0;
        lsq_offset = 0.0;
    }
    else
    {
        lsq_scale  = (sumxy*valid_pixels - sumx*sumy)/denom;
        lsq_offset = (sumxx*sumy - sumx*sumxy)/denom;
    }
    printf("xy lsq scale= %lg offset= %lg\n",lsq_scale,lsq_offset);

    denom = (sumyy*valid_pixels)-sumy*sumy;
    if( denom == 0 || valid_pixels < 2 )
    {
        lsq_scale = 0.0;
        lsq_offset = 0.0;
    }
    else
    {
        lsq_scale  = (sumxy*valid_pixels - sumx*sumy)/denom;
        lsq_offset = (sumyy*sumx - sumy*sumxy)/denom;
    }
    printf("yx lsq scale= %lg offset= %lg\n",lsq_scale,lsq_offset);

    sumd=sumdsq=sumd3=sumd4=0.0;
    for(j=0;j<pixels;++j)
    {
        i=j+header/sizeof(float);

        if(invalid_pixel[i]) continue;
        diff = outimage[i]-outoffset - avg;
        sumd   += diff;
        sumdsq += diff*diff;
        sumd3  += diff*diff*diff;
        sumd4  += diff*diff*diff*diff;
    }
    rmsd = sqrt(sumdsq/valid_pixels);
    skewness = sumd3/valid_pixels/(rmsd*rmsd*rmsd);
    kurtosis = sumd4/valid_pixels/(rmsd*rmsd*rmsd*rmsd);


    if(histogramize)
    {
        /* default to reasonable bin size */
        if(binsize <= 0)
        {
            binsize = rmsd/3.0;
        }
        /* data are all same value? */
        if(binsize <= 0) binsize = 1.0;

        /* calculate number of bins needed for histogram, hope its not ridiculous */
        bins = ( ceil( (max - min + 2) / binsize ) );
        if(bins > 100000) bins = 100000;
        histogram = calloc(bins+1,sizeof(unsigned long int));
        if(histogram == NULL)
        {
            perror("histogram");
        }
        else
        {
            for(i=1;i<=valid_pixels;++i)
            {
                bin = ( ceil( (medbuffer[i] - min) / binsize ) - 0.5 );
                ++histogram[bin];
            }
            printf("histogram:\n");
            for(bin=0;bin<=bins;++bin)
            {
                if(histogram[bin]>0)
                {
                    printf("%g %ld\n",binsize*bin+min,histogram[bin]);
                }
            }
        }
    }

    avgx = sumx/valid_pixels;
    avgy = sumy/valid_pixels;
    avgxy = sumxy/valid_pixels;
    avgxx = sumxx/valid_pixels;
    avgyy = sumyy/valid_pixels;
    if(avgyy<avgy*avgy)avgyy=avgy*avgy;
    if(avgxx<avgx*avgx){
        CC=0.0;
    }
    else
    {
        CC=(avgxy-avgx*avgy)/(sqrt(avgxx-avgx*avgx)*sqrt(avgyy-avgy*avgy));
    }


    if(ignore_values || reject_outliers) printf("%d invalid of ",pixels-valid_pixels);
    printf("%d pixels ",pixels);
    if(ignore_values || reject_outliers) printf("( %d left)",valid_pixels);
    printf("\n");
    printf("%d overflows, %d underflows\n",overflows,underflows);
    printf("max = %g min = %g\n",max,min);
    printf("mean = %g rms = %g rmsd = %g skewness = %g kurtosis = %g\n",avg,rms,rmsd,skewness,kurtosis);
    if(reject_outliers)
    {
        median = fmedian_with_rejection(valid_pixels+1,medbuffer,outlier_sigma,&mad,&n_arej);
        printf("median = %g mad = %g\n",median,mad);
        for(j=1;j<=n_arej;++j)
        {
        //    deviate = sqrt(medbuffer[j]*medbuffer[j]);
        //    if(deviate>4*mad) printf("GOTHERE: medbuffer[%d]=%g >%g\n",j,medbuffer[j],4*mad);
            outfloat = medbuffer[j]-outoffset;
            if(j==1) max_arej=min_arej=outfloat; 
            if(outfloat>max_arej)max_arej=outfloat;
            if(outfloat<min_arej)min_arej=outfloat;
        }
        printf("%d pixels rejected after median-mad filter, new max= %g min= %g\n",valid_pixels-n_arej,max_arej,min_arej);
        avg_arej = fmean_with_rejection(n_arej,medbuffer,outlier_sigma,&rmsd_arej,&n_arej);
        sumsq_arej=sumd_arej=sumdsq_arej=sumd3_arej=sumd4_arej=0.0;
        for(j=1;j<=n_arej;++j)
        {
            outfloat = medbuffer[j]-outoffset;
            if(j==1) max_arej=min_arej=outfloat; 
            if(outfloat>max_arej)max_arej=outfloat;
            if(outfloat<min_arej)min_arej=outfloat;
            sumsq_arej   += outfloat*outfloat;
            diff = outfloat - avg_arej;
            sumd_arej   += diff;
            sumdsq_arej += diff*diff;
            sumd3_arej  += diff*diff*diff;
            sumd4_arej  += diff*diff*diff*diff;
        }
        rms_arej=sqrt(sumsq_arej/n_arej);
        rmsd_arej=sqrt(sumdsq_arej/n_arej);
        skewness_arej = sumd3_arej/n_arej/(rmsd_arej*rmsd_arej*rmsd_arej);
        kurtosis_arej = sumd4_arej/n_arej/(rmsd_arej*rmsd_arej*rmsd_arej*rmsd_arej);
        printf("%d pixels rejected after mean-rmsd filter, new max= %g min= %g\n",valid_pixels-n_arej,max_arej,min_arej);
        printf("after outlier rejection: mean= %g rms= %g rmsd= %g skewness= %g kurtosis= %g ( %d points)\n",avg_arej,rms_arej,rmsd_arej,skewness_arej,kurtosis_arej,n_arej);

    }
    printf("CC = %.12g\n",CC);

    if(outfile != NULL)
    {
        printf("writing %s as a %d-byte header with %d %ld-byte floats\n",outfilename,header,pixels,sizeof(float));
        fwrite(outimage,n,1,outfile);
        fclose(outfile);
    }

    return 0;
}


/* find the median value of an array of data - scrambles the array */
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
float fmedian(unsigned int n, float arr[])
{
    unsigned int i,j,k,l,ir,mid;
    float a,temp;

    l=1;
    ir=n;
    k=(n+1)/2;
//printf("n=%d; k=%d\n",n,k);

//for(i=1;i<=n;++i) printf("arr[%d]=%f\n",i,arr[i]);

    for(;;)
    {
        if(ir <= l+1)
        {
            if(ir == l+1 && arr[ir] < arr[l])
            {
                SWAP(arr[l],arr[ir]);
            }
//for(i=1;i<=n;++i) printf("arr[%d]=%f\n",i,arr[i]);
            return arr[k];
        } else {
            mid=(l+ir) >> 1;
            SWAP(arr[mid],arr[l+1]);
            if(arr[l+1] > arr[ir])
            {
                SWAP(arr[l+1],arr[ir]);
            }
            if(arr[l] > arr[ir])
            {
                SWAP(arr[l],arr[ir]);
            }
            if(arr[l+1] > arr[l])
            {
                SWAP(arr[l+1],arr[l]);
            }
            i=l+1;      // initialize pointers for partitioning
            j=ir;
            a=arr[l];   // partitioning element
            for(;;)     // innermost loop
            {
                do i++; while(arr[i]<a);        // scan up to find element > a
                do j--; while(arr[j]>a);        // scan down to find element < a
                if( j < i ) break;              // pointers crossed, median is in between
                SWAP(arr[i],arr[j]);
            }
            arr[l]=arr[j];                      // insert partitioning element
            arr[j]=a;
            if( j >= k ) ir=j-1;                // Keep partition that contains the median active
            if( j <= k ) l=i;
        }
    }
}


float fmedian_with_rejection(unsigned int n, float arr[],float sigma_cutoff, float *final_mad, int *final_n)
{
    float median_value;
    int i,orig_n,done;
    float min_frac,deviate,mad;

    orig_n = n;
    min_frac = 0.7;

    done = 0;
    while(! done)
    {
        /* compute the median (centroid) value */
        median_value = fmedian(n,arr);

        /* now figure out what the mean absolute deviation from this value is */
        mad = fmedian_absolute_deviation(n,arr,median_value);
        //if(flag) printf("mad = %f\n",mad);

        done = 1;
        /* reject all outliers */
        for(i=1;i<=n;++i)
        {
            /* reject positive and negative outliers */
            deviate = fabs(arr[i]-median_value);
            if(deviate > sigma_cutoff*mad)
            {
                /* needs to go */
                /* move value at the end of the array to this "reject" and then shorten the array */
                //if(flag) printf("rejecting arr[%d] = %f (%f)\n",i,arr[i],deviate);
                //arr[worst]+=10000;
                if(i != n)
                {
                    //temp=arr[worst];
                    arr[i] = arr[n];
                    //arr[n]=temp;
                }
                --n;
                done = 0;
            }
        }
    }

    /* basically two return values */
    *final_mad = mad;
    *final_n = n;
    return median_value;
}


/* note: there must be 2*n elements in this array! */
float fmedian_absolute_deviation(unsigned int n, float arr[], float median_value)
{
    int i;
    for(i=1;i<=n;++i)
    {
        arr[i+n] = fabs(arr[i]-median_value);
    }

    return fmedian(n,arr+n);
}



/* this function keeps track of outliers by swapping them to the end of the array */
/* counting starts at 0 and "points" is the number of points */
float fmean_with_rejection(unsigned int starting_points, float arr[], float sigma_cutoff, double *final_rmsd, int *final_n)
{
    int points,i;
    int rejection,worst;
    float temp,sum,avg,sumd,rmsd,deviate,worst_deviate;

    points=starting_points;
    rejection = 1;
    while ( rejection && points>starting_points/2.0 )
    {
        /* find the mean and rms deivation */
        sum = sumd = 0.0;
        for(i=0;i<points;++i)
        {
            sum+=arr[i];
        }
        avg=sum/points;
        worst=-1;
        worst_deviate=-1.0;
        for(i=0;i<points;++i)
        {
            deviate=fabs(arr[i]-avg);
            if(deviate > worst_deviate)
            {
                worst=i;
                worst_deviate=deviate;
            }
            sumd+=deviate*deviate;
        }
        rmsd=sqrt(sumd/points);

        rejection=0;
        if(worst_deviate>sigma_cutoff*rmsd)
        {
            /* we have a reject! */
            rejection=1;
//printf("GOTHERE: rejecting arr[%d] = %f = %f > %f * %f\n",worst,arr[worst],worst_deviate,sigma_cutoff,rmsd);

            /* move it to end of the array and forget about it */
            SWAP(arr[worst],arr[points]);
            --points;
        }
    }

    *final_rmsd = rmsd;
    *final_n = points;
    return avg;
}
