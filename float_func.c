
/* apply a C function to each value in a raw "float" input file                -James Holton               8-24-25

example:

gcc -O -O -o float_func float_func.c -lm
./float_func -func erf snr.bin occ.bin 

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char *infile1name = NULL;
char *infile2name = NULL;
FILE *infile1 = NULL;
FILE *infile2 = NULL;
char *outfilename = "output.bin\0";
FILE *outfile = NULL;

typedef enum { UNKNOWN, SQRT, CBRT, CEIL, FLOOR, FABS,
       SIGN,
       SIN, ASIN, SINH, ASINH, 
       COS, ACOS, COSH, ACOSH,
       TAN, ATAN, TANH, ATANH,
       ERF, ERFC, EXP, SAFEXP, LOG, LOG10,
       J0, J1, JN, Y0, Y1, YN, GAMMA, LGAMMA,
        POW, ERFPOW,
       NORM,
       URAND, GRAND, LRAND, PRAND, ERAND, TRAND,
       SET,
       ADD, SUBTRACT, MULTIPLY, DIVIDE, INVERSE, NEGATE,
       MAXIMUM, MINIMUM, NANZERO,
       THRESH,
       FFT, INVFFT, REALFFT, INVREALFFT, AB2PHIF, PHIF2AB,
       SWAB2, SWAB4,
       ODD, EVEN, ODDEVEN, EVENODD,
       REVERSE, STRIDE, STITCH, FLIP, FLOP, FLIPFLOP, SWAPXY,
       MAXPOOL, AVGBOX, GAUSSBLUR,
       MAXRADIUS, AVGRADIUS, MEDRADIUS,
       SEGMENT, EDGE, FEATHER } func_name;

/* for byte swapping */
typedef union {
    unsigned char byte[4];
    unsigned short int word[2];
    float floaty;
} FOURBYTES;
FOURBYTES swapper,swappee;

/* random deviate with Poisson distribution */
float poidev(float xm, long *idum);
/* random deviate with Gaussian distribution */
float gaussdev(long *idum);
/* random deviate with Lorentzian distribution */
float lorentzdev(long *idum);
/* random deviate with triangle-shaped distribution */
float triangledev(long *idum);
/* random deviate with exponential distribution (>0) */
float expdev(long *idum);
/* random deviate with uniform distribution */
float ran1(long *idum);
/* median value of a list */
float fmedian(unsigned long n, float arr[]);
/* FFT function */
float *Fourier(float *data, unsigned long length, int direction);

long seed;

int debug = 0;

int main(int argc, char** argv)
{
     
    long n,i,j,k,pixels,outpixels,samples;
    long xsize=0,ysize=0,zsize=0,x,y,z;
    long xosize=0,yosize=0,zosize=0,xo,yo,zo;
    long xstart=0,ystart=0,zstart=0;
    float *outimage;
    float *inimage1;
    float *inimage2;
    float *scratch;
    long *count;
    long *segments,*joins,*segmap,nsegments=0,njoins=0;
    char *headerstuff;
    long header=0,outheader=0;
    func_name func = UNKNOWN;
    float param=1.0;
    int parami;
    int user_param = 0, map_param = 0;
    float sum,sumd,sumsq,sumdsq,avg,rms,rmsd,min,max;
    int ignore_values=0,valid_pixels=0;
    float ignore_value[70000];
    unsigned short int *invalid_pixel;
    float phi,F,a,b;
    float radius=-1.0,fradsq,fx,fy,fz,fwidthsq;
    int dx,dy,dz,narr,xrad=-1,yrad=-1,zrad=-1;

    /* check argument list */
    for(i=1; i<argc; ++i)
    {
        if(strlen(argv[i]) > 4)
        {
            if(strstr(argv[i]+strlen(argv[i])-4,".bin") || strstr(argv[i]+strlen(argv[i])-4,".map"))
            {
                printf("filename: %s\n",argv[i]);
                if(infile1name == NULL){
                    infile1name = argv[i];
                }
                else
                {
                    if(infile2name == NULL){
                        infile2name = argv[i];
                    }
                    else
                    {
                        outfilename = argv[i];
                    }
                }
            }
        }

        if(argv[i][0] == '-')
        {
            /* option specified */
            if(strstr(argv[i], "-debug"))
            {
                ++debug;
                continue;
            }
            if(strstr(argv[i], "-func") && (argc >= (i+1)))
            {
                func = UNKNOWN;
                if(strstr(argv[i+1],"sqrt")) func = SQRT;
                if(strstr(argv[i+1],"cbrt")) func = CBRT;
                if(strstr(argv[i+1],"ceil")) func = CEIL;
                if(strstr(argv[i+1],"floor")) func = FLOOR;
                if(strstr(argv[i+1],"abs")) func = FABS;
                if(strstr(argv[i+1],"fabs")) func = FABS;
                if(strstr(argv[i+1],"sign")) func = SIGN;
                if(strstr(argv[i+1],"sin")) func = SIN;
                if(strstr(argv[i+1],"asin")) func = ASIN;
                if(strstr(argv[i+1],"sinh")) func = SINH;
                if(strstr(argv[i+1],"asinh")) func = ASINH;
                if(strstr(argv[i+1],"cos")) func = COS;
                if(strstr(argv[i+1],"acos")) func = ACOS;
                if(strstr(argv[i+1],"cosh")) func = COSH;
                if(strstr(argv[i+1],"acosh")) func = ACOSH;
                if(strstr(argv[i+1],"tan")) func = TAN;
                if(strstr(argv[i+1],"atan")) func = ATAN;
                if(strstr(argv[i+1],"tanh")) func = SINH;
                if(strstr(argv[i+1],"atanh")) func = ATANH;
                if(strstr(argv[i+1],"erf")) func = ERF;
                if(strstr(argv[i+1],"erfc")) func = ERFC;
                if(strstr(argv[i+1],"pow")) func = POW;
                if(strstr(argv[i+1],"erfpow")) func = ERFPOW;
                if(strstr(argv[i+1],"exp")) func = EXP;
                if(strstr(argv[i+1],"safexp")) func = SAFEXP;
                if(strstr(argv[i+1],"log")) func = LOG;
                if(strstr(argv[i+1],"log10")) func = LOG10;
                if(strstr(argv[i+1],"j0")) func = J0;
                if(strstr(argv[i+1],"j1")) func = J1;
                if(strstr(argv[i+1],"jn")) func = JN;
                if(strstr(argv[i+1],"y0")) func = Y0;
                if(strstr(argv[i+1],"y1")) func = Y1;
                if(strstr(argv[i+1],"yn")) func = YN;
                if(strstr(argv[i+1],"gamma")) func = GAMMA;
                if(strstr(argv[i+1],"lgamma")) func = LGAMMA;
                if(strstr(argv[i+1],"norm")) func = NORM;
                if(strstr(argv[i+1],"urand")) func = URAND;
                if(strstr(argv[i+1],"grand")) func = GRAND;
                if(strstr(argv[i+1],"lrand")) func = LRAND;
                if(strstr(argv[i+1],"prand")) func = PRAND;
                if(strstr(argv[i+1],"erand")) func = ERAND;
                if(strstr(argv[i+1],"trand")) func = TRAND;
                if(strstr(argv[i+1],"zero")) {   func = SET; param = 0.0; user_param=1;};
                if(strstr(argv[i+1],"one")) {    func = SET; param = 1.0; user_param=1;};
                if(strstr(argv[i+1],"unity")) {  func = SET; param = 1.0; user_param=1;};
                if(strstr(argv[i+1],"set"))      func = SET;
                if(strstr(argv[i+1],"const"))    func = SET;
                if(strstr(argv[i+1],"constant")) func = SET;
                if(strstr(argv[i+1],"add"))      func = ADD;
                if(strstr(argv[i+1],"subtract")) func = SUBTRACT;
                if(strstr(argv[i+1],"multiply")) func = MULTIPLY;
                if(strstr(argv[i+1],"mult"))     func = MULTIPLY;
                if(strstr(argv[i+1],"divide"))   func = DIVIDE;
                if(strstr(argv[i+1],"div"))      func = DIVIDE;
                if(strstr(argv[i+1],"inverse"))  func = INVERSE;
                if(strstr(argv[i+1],"inv"))      func = INVERSE;
                if(strstr(argv[i+1],"negative")) func = NEGATE;
                if(strstr(argv[i+1],"negate"))   func = NEGATE;
                if(strstr(argv[i+1],"max"))      func = MAXIMUM;
                if(strstr(argv[i+1],"maximum"))  func = MAXIMUM;
                if(strstr(argv[i+1],"min"))      func = MINIMUM;
                if(strstr(argv[i+1],"minimum"))  func = MINIMUM;
                if(strstr(argv[i+1],"nanzero"))  func = NANZERO;
                if(strstr(argv[i+1],"thresh"))   func = THRESH;
                if(strstr(argv[i+1],"fft"))      func = FFT;
                if(strstr(argv[i+1],"invfft"))   func = INVFFT;
                if(strstr(argv[i+1],"realfft"))  func = REALFFT;
                if(strstr(argv[i+1],"invrealfft")) func = INVREALFFT;
                if(strstr(argv[i+1],"ab2phif"))   func = AB2PHIF;
                if(strstr(argv[i+1],"phif2ab"))   func = PHIF2AB;
                if(strstr(argv[i+1],"swab"))     func = SWAB4;
                if(strstr(argv[i+1],"swab2"))    func = SWAB2;
                if(strstr(argv[i+1],"swap"))     func = SWAB4;
                if(strstr(argv[i+1],"swap2"))    func = SWAB2;
                if(strstr(argv[i+1],"odd"))      func = ODD;
                if(strstr(argv[i+1],"even"))     func = EVEN;
                if(strstr(argv[i+1],"oddeven"))  func = ODDEVEN;
                if(strstr(argv[i+1],"evenodd"))  func = EVENODD;
                if(strstr(argv[i+1],"reverse"))  func = REVERSE;
                if(strstr(argv[i+1],"stride"))   func = STRIDE;
                if(strstr(argv[i+1],"stitch"))   func = STITCH;
                if(strstr(argv[i+1],"flip"))     func = FLIP;
                if(strstr(argv[i+1],"flop"))     func = FLOP;
                if(strstr(argv[i+1],"flipflop")) func = FLIPFLOP;
                if(strstr(argv[i+1],"swapxy")) func = SWAPXY;
                if(strstr(argv[i+1],"maxpool"))  func = MAXPOOL;
                if(strstr(argv[i+1],"avgbox"))   func = AVGBOX;
                if(strstr(argv[i+1],"segment"))  func = SEGMENT;
                if(strstr(argv[i+1],"edge"))     func = EDGE;
                if(strstr(argv[i+1],"feather"))  func = FEATHER;
                if(strstr(argv[i+1],"gaussblur"))  func = GAUSSBLUR;
                if(strstr(argv[i+1],"maxradius"))  func = MAXRADIUS;
                if(strstr(argv[i+1],"avgradius"))  func = AVGRADIUS;
                if(strstr(argv[i+1],"medradius"))  func = MEDRADIUS;
                if(strstr(argv[i+1],"+")) func = ADD;
                if(strstr(argv[i+1],"-")) func = SUBTRACT;
                if(strstr(argv[i+1],"*")) func = MULTIPLY;
                if(strstr(argv[i+1],"/")) func = DIVIDE;
                if(strstr(argv[i+1],"1/")) {func = INVERSE ; param = 1.0; user_param=1;};
                if(strstr(argv[i+1],"1-")) {func = NEGATE ; param = 1.0; user_param=1;};
                if(func == UNKNOWN) printf("WARNING: unknown function: %s\n",argv[i+1]);
            }
            if(strstr(argv[i], "-param") && (argc >= (i+1)))
            {
                param = atof(argv[i+1]);
                 user_param=1;
            }
            if(strstr(argv[i], "-xsize") && (argc >= (i+1)))
            {
                xsize = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-ysize") && (argc >= (i+1)))
            {
                ysize = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-zsize") && (argc >= (i+1)))
            {
                zsize = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-radius") && (argc >= (i+1)))
            {
                radius = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-xrad") && (argc >= (i+1)))
            {
                xrad = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-yrad") && (argc >= (i+1)))
            {
                yrad = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-zrad") && (argc >= (i+1)))
            {
                zrad = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-ignore") && (argc >= (i+1)))
            {
                ++ignore_values;
                ignore_value[ignore_values] = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-seed") && (argc >= (i+1)))
            {
                seed = -atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-box") || strstr(argv[i], "-roi"))
            {
                if(argc >= (i+2))
                {
                    xstart  = atol(argv[i+1]);
                    xosize  = atol(argv[i+2]);
                    xosize -= xstart;
                }
                if(argc >= (i+4))
                {
                    ystart  = atol(argv[i+3]);
                    yosize  = atol(argv[i+4]);
                    yosize -= ystart;
                }
                if(argc >= (i+6))
                {
                    zstart  = atol(argv[i+5]);
                    zosize  = atol(argv[i+6]);
                    zosize -= zstart;
                }
            }
            if((strstr(argv[i], "-output") || strstr(argv[i], "-outfile")) && (argc >= (i+1)))
            {
                outfilename = argv[i+1];
                ++i;
                continue;
            }
            if((strstr(argv[i], "-input1") || strstr(argv[i], "-infile1")) && (argc >= (i+1)))
            {
                infile1name = argv[i+1];
                ++i;
                continue;
            }
            if((strstr(argv[i], "-input2") || strstr(argv[i], "-infile2")) && (argc >= (i+1)))
            {
                infile2name = argv[i+1];
                ++i;
                continue;
            }
            if(strstr(argv[i], "-header") && (argc >= (i+1)))
            {
                header = outheader = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-outheader") && (argc >= (i+1)))
            {
                outheader = atof(argv[i+1]);
            }
            ++i;
        }
    }

    if(infile1name != NULL){
        if( strstr(infile1name,"-") && strlen(infile1name) == 1 )
        {
            infile1 = stdin;
        }
        else {
            infile1 = fopen(infile1name,"rb");
        }
    }
    if(infile1 == NULL){
        printf("ERROR: unable to open %s as input 1\n", infile1name);
        perror("");
    }
    if(infile1 == NULL || func == UNKNOWN){
        printf("usage: float_func file.bin [file2.bin] [outfile.bin] -func jn -param 1.0 -ignore 0\n");
        printf("options:\n");
        printf("\tfile.bin\t binary file containing the arguments to the desired function\n");
        printf("\tfile2.bin\t binary file containing second argument to the desired function\n");
        printf("\t-func\t may be one of:\n");
        printf("\t sqrt cbrt ceil floor abs sign\n");
        printf("\t zero one set (one output value)\n");
        printf("\t add subtract multiply divide inverse negate maximum minimum \n");
        printf("\t thresh (1-or-0 threshold)\n");
        printf("\t nanzero (set all NaN and Inf values to zero)\n");
        printf("\t sin asin sinh asinh cos acos cosh acosh tan atan tanh atanh (trigonometry)\n");
        printf("\t pow erf erfpow norm erfc (power and error functions)\n");
        printf("\t exp log log10 (natural or base-10 log)\n");
        printf("\t safexp (safe exponential, wont over/underflow)\n");
        printf("\t j0 j1 jn y0 y1 yn gamma lgamma (Bessel and gamma functions)\n");
        printf("\t urand grand lrand prand erand trand (uniform gaussian lorentzian poisson exponential or triangle randomness)\n");
        printf("\t fft invfft realfft invrealfft (complex or real FFT - base 2)\n");
        printf("\t ab2phif phif2ab (cartesian to amplitude-phase conversion)\n");
        printf("\t swab4 swab2 (swap bytes)\n");
        printf("\t odd even oddeven evenodd (extract or interleave values)\n");
        printf("\t reverse (re-order 4-byte floats in file)\n");
        printf("\t stride (skip blocks of size param)\n");
        printf("\t stitch (interleave blocks from different files)\n");
        printf("\t flip (reverse data in interleaved blocks)\n");
        printf("\t flop (reverse data as interleaved blocks)\n");
        printf("\t flipflop (reverse x and y data with -xsize or -param xsize)\n");
        printf("\t swapxy (swap x and y axes with -xsize)\n");
        printf("\t maxpool (maximum value within box with -param edge pixels)\n");
        printf("\t avgbox (average value of box with -param edge pixels)\n");
        printf("\t edge (find rising/falling edges, return 0,1,or -1)\n");
        printf("\t feather (find edges and soften them)\n");
        printf("\t maxradius (maximum value within radius xrad, yrad, zrad pixels no down-sampling)\n");
        printf("\t avgradius (local average within radius xrad, yrad, zrad  pixels, no down-sampling)\n");
        printf("\t medradius (local median within radius xrad, yrad, zrad  pixels, no down-sampling)\n");
        printf("\t-param\t second value for add, subtract, set, etc. or parameter needed by the function (such as jn)\n");
        printf("\t-xsize\t fast dimension when treating data as 2D or 3D array\n");
        printf("\t-ysize\t 2nd fastest dimension when treating data as 2D or 3D array\n");
        printf("\t-zsize\t slowest dimension when treating data as 3D array\n");
        printf("\t-xrad\t max radial extent for gaussblur maxradius avgradius and medradius functions\n");
        printf("\t-yrad -zrad\t radial extent in y and directions\n");
        printf("\t-seed\t seed for random number functions\n");
        printf("\t-ignore\t value found in either input file will pass through\n");
        printf("\t-header \tnumber of bytes to ignore in header of each file\n");
        printf("\t-outheader \tnumber of bytes of the header to pass on to output file\n");
        printf("\toutput.bin\t output binary file containing the results to the desired function\n");
        exit(9);
    }

    printf("selected function: ");
    switch(func){
        case AB2PHIF: printf("AB2PHIF\n"); break;
        case PHIF2AB: printf("PHIF2AB\n"); break;
        case ACOS: printf("ACOS\n"); break;
        case ACOSH: printf("ACOSH\n"); break;
        case ADD: printf("ADD\n"); break;
        case ASIN: printf("ASIN\n"); break;
        case ASINH: printf("ASINH\n"); break;
        case ATAN: printf("ATAN\n"); break;
        case ATANH: printf("ATANH\n"); break;
        case CBRT: printf("CBRT\n"); break;
        case CEIL: printf("CEIL\n"); break;
        case COS: printf("COS\n"); break;
        case COSH: printf("COSH\n"); break;
        case DIVIDE: printf("DIVIDE\n"); break;
        case ERAND: printf("ERAND\n"); break;
        case ERF: printf("ERF\n"); break;
        case ERFC: printf("ERFC\n"); break;
        case ERFPOW: printf("ERFPOW\n"); break;
        case EVEN: printf("EVEN\n"); break;
        case EXP: printf("EXP\n"); break;
        case SAFEXP: printf("SAFEXP\n"); break;
        case FABS: printf("FABS\n"); break;
        case SIGN: printf("SIGN\n"); break;
        case FFT: printf("FFT\n"); break;
        case INVFFT: printf("INVFFT\n"); break;
        case REALFFT: printf("REALFFT\n"); break;
        case INVREALFFT: printf("INVREALFFT\n"); break;
        case FLOOR: printf("FLOOR\n"); break;
        case GAMMA: printf("GAMMA\n"); break;
        case GRAND: printf("GRAND\n"); break;
        case INVERSE: printf("INVERSE\n"); break;
        case J0: printf("J0\n"); break;
        case J1: printf("J1\n"); break;
        case JN: printf("JN\n"); break;
        case LGAMMA: printf("LGAMMA\n"); break;
        case LOG: printf("LOG\n"); break;
        case LOG10: printf("LOG10\n"); break;
        case LRAND: printf("LRAND\n"); break;
        case MAXIMUM: printf("MAXIMUM\n"); break;
        case MINIMUM: printf("MINIMUM\n"); break;
        case MULTIPLY: printf("MULTIPLY\n"); break;
        case NANZERO: printf("NANZERO\n"); break;
        case NEGATE: printf("NEGATE\n"); break;
        case NORM: printf("NORM\n"); break;
        case ODD: printf("ODD\n"); break;
        case ODDEVEN: printf("ODDEVEN\n"); break;
        case EVENODD: printf("EVENODD\n"); break;
        case REVERSE: printf("REVERSE\n"); break;
        case STRIDE: printf("STRIDE\n"); break;
        case STITCH: printf("STITCH\n"); break;
        case FLIP: printf("FLIP\n"); break;
        case FLOP: printf("FLOP\n"); break;
        case FLIPFLOP: printf("FLIPFLOP\n"); break;
        case SWAPXY: printf("SWAPXY\n"); break;
        case MAXPOOL: printf("MAXPOOL\n"); break;
        case AVGBOX: printf("AVGBOX\n"); break;
        case POW: printf("POW\n"); break;
        case SEGMENT: printf("SEGMENT %f\n",param); break;
        case EDGE: printf("EDGE %f\n",param); break;
        case FEATHER: printf("FEATHER %f\n",param); break;
        case GAUSSBLUR: printf("GAUSSBLUR %f  %f %f %f\n",param,xrad,yrad,zrad); break;
        case AVGRADIUS: printf("AVGRADIUS %f %f %f\n",xrad,yrad,zrad); break;
        case MEDRADIUS: printf("MEDRADIUS %f %f %f\n",xrad,yrad,zrad); break;
        case PRAND: printf("PRAND\n"); break;
        case SET: printf("SET\n"); break;
        case SIN: printf("SIN\n"); break;
        case SINH: printf("SINH\n"); break;
        case SQRT: printf("SQRT\n"); break;
        case SUBTRACT: printf("SUBTRACT\n"); break;
        case SWAB2: printf("SWAB2\n"); break;
        case SWAB4: printf("SWAB4\n"); break;
        case TAN: printf("TAN\n"); break;
        case TANH: printf("TANH\n"); break;
        case THRESH: printf("THRESH\n"); break;
        case TRAND: printf("TRAND\n"); break;
        case UNKNOWN: printf("UNKNOWN\n"); break;
        case URAND: printf("URAND\n"); break;
        case Y0: printf("Y0\n"); break;
        case Y1: printf("Y1\n"); break;
        case YN: printf("YN\n"); break;
    }
/*
awk '/func = / && $NF !~ /}/{print substr($NF,1,length($NF)-1);}' ~/projects/bin_stuff/float_func.c |\
 sort -u | awk '{ print "\tcase "$1": printf(\""$1"\\n\"); break;"}'
*/

    for(k=1;k<=ignore_values;++k)
    {
        printf("ignoring value %g in both files\n",ignore_value[k]);
    }

    /* load first float-image */
    printf("header = %d bytes\n",header);
    if(header)
    {
        headerstuff = calloc(header+10,1);
        rewind(infile1);
        fread(headerstuff,header,1,infile1);
    }
    fseek(infile1,0,SEEK_END);
    n = ftell(infile1)-header;
    inimage1 = calloc(n+10,1);
    invalid_pixel = calloc(n+10,1);
    printf("reading %u floats from %s\n",n/sizeof(float),infile1name);
    fseek(infile1,header,SEEK_SET);
    fread(inimage1,n,1,infile1);
    fclose(infile1);

    /* assume file is all pixels */
    pixels = samples = outpixels = n/sizeof(float);

    /* handy integer version of param */
    parami = param;
    /* apply defaults to array dimensions */
    if(xsize <= 0 && ( func == FLIP || func == FLOP || func == FLIPFLOP || func == SWAPXY ) ) xsize = abs(parami);
    /* nothing specified = all x */
    if(xsize <= 0 && ysize <= 0 && zsize <= 0) xsize = pixels;
    if(xsize <= 0) xsize = 1;
    if(ysize <= 0) ysize = pixels/xsize;
    if(ysize <= 0) ysize = 1;
    zsize = pixels/xsize/ysize;
    if(zsize <= 0) zsize = 1;
    printf("array dimensions: x %d y %d z %d\n",xsize,ysize,zsize);

    /* if the function takes two args ... */
    if( func == ADD || func == SUBTRACT || func == MULTIPLY || func == DIVIDE || func == POW ||
        func == MAXIMUM || func == MINIMUM || func == THRESH ||
        func == EVENODD || func == ODDEVEN ||
        func == STITCH ) {
        if(infile2name == NULL || user_param )
        {
            /* use the "param" as the second value */
            map_param = 0;
        }
        else
        {
            map_param = 1;
        }
    }
    else
    {
        if(infile2name != NULL)
        {
            /* second filename must be the output file */
            outfilename = infile2name;
            infile2name = NULL;
        }
    }

    /* open second file, if it was specified */
    if(infile2name != NULL){
        infile2 = fopen(infile2name,"r");
        if(infile2 == NULL){
            printf("ERROR: unable to open %s as input 2\n", infile2name);
            perror("");
            exit(9);
        }            
        inimage2 = calloc(n,1);
        fseek(infile2,header,SEEK_SET);
        fread(inimage2,n,1,infile2);
        fclose(infile2);
    }

    outfile = fopen(outfilename,"wb");
    if(outfile == NULL)
    {
        printf("ERROR: unable to open %s for output\n", outfilename);
        perror("");
        exit(9);
    }
    printf("input1 is: %s\n",infile1name);
    if(map_param)
    {
        printf("input2 is: %s\n",infile2name);
    }else{
        printf("input2 is: %g\n",param);
    }
    printf("output to: %s\n",outfilename);
    
    /* just to keep track */
    outimage = NULL;
    count = NULL;

    /* special functions */
    if( func == FFT || func == INVFFT || func == REALFFT )
    {
        /* must be a power of two */
        valid_pixels = (int) ceil(pow(2.0,floor(log(pixels)/log(2))));
        if(valid_pixels != pixels){
            pixels = outpixels = xsize = valid_pixels;
            printf("truncating data to %d floats\n",pixels);
        }
        valid_pixels = 0;
    }
    if( func == INVREALFFT )
    {
        /* must be a power of two plus 2 */
        valid_pixels = 2+(int) ceil(pow(2.0,floor(log(pixels-2)/log(2))));
        if(valid_pixels != pixels){
            pixels = outpixels = xsize = valid_pixels;
            printf("truncating data to %d floats\n",pixels);
        }
        valid_pixels = 0;
    }

    if( func == FFT ){
        /* assume input array is complex numbers (even=real, odd=imag) */
        /* fortranish, so zeroeth item will not be touched */
        samples = pixels/2;
        outimage = Fourier(inimage1-1, samples, 1);
        /* correct for fortranish */
        ++outimage;
        /* output file size in floats */
        outpixels = pixels;
        /* correct coefficients so that they back-transform on the same scale */
        for(j=0;j<outpixels;++j)
        {
            outimage[j]/=(outpixels/2);
        }
    }
    if( func == INVFFT ){
        /* assume input array is complex numbers (even=real, odd=imag) */
        /* fortranish, so zeroeth item will not be touched */
        samples = pixels/2;
        outimage = Fourier(inimage1-1, samples, -1);
        /* correct for fortranish */
        ++outimage;
        /* output file size in floats */
        outpixels = pixels;
    }
    if( func == REALFFT ){
        /* need to convert to complex numbers */
        samples = pixels+2;
        outimage = calloc(samples,2*sizeof(float));
        /* copy input file, staggering into output array */
        for(i=0;i<pixels;++i)
        {
            j = 2*i;
            outimage[j] = inimage1[i]/(pixels);
        }
        /* FFT array needs to be twice the size of input file */
        /* fortranish, so zeroeth item will not be touched */
        outimage = Fourier(outimage-1, pixels, 1);
        /* correct for fortranish */
        ++outimage;
        /* note: we are only interested in first half (plus one more for the Nyquist frequency) */
//        pixels=pixels+2;
        /* output file size in floats has two extra floats for F0 and Fnyquist */
        outpixels = samples;
    }
    if( func == INVREALFFT ){
        /* realfft has two extra floats for F0 and Fnyquist */
        samples = pixels-2;
//        pixels-=2;
        outimage = calloc(samples,2*sizeof(float));
        /* compensate for zero-value negative index coefficients by doubling non-trivial ones */
        outimage[0]=inimage1[0];
        outimage[samples]=inimage1[samples];
        for(j=1;j<samples;++j)
        {
            i=j;
            outimage[j] = 2*inimage1[i];
        }
        /* fortranish, so zeroeth item will not be touched */
        outimage = Fourier(outimage-1, samples, -1);
        /* correct for fortranish */
        ++outimage;
        /* but we are only interested in even numbered ones (real part) */
        for(j=0;j<samples;++j)
        {
            outimage[j] = outimage[2*j];
        }
        /* output file size in floats is two floats less than input file size */
        outpixels = samples;
    }

    /* other functions with inequal input and output file sizes */
    if( func == ODD || func == EVEN ){
        /* output will be half the size */
        outpixels = pixels/2;
    }
    if( func == ODDEVEN || func == EVENODD ){
        /* interleaving values, so output will be twice the size */
        outpixels = 2*pixels;
        outimage = calloc(outpixels,sizeof(float));
    }
    if( func == AB2PHIF || func == PHIF2AB ){
        /* values of interest are two floats */
        outimage = calloc(pixels,sizeof(float));
        outpixels = pixels;
        samples = pixels/2;
        /* we will restore true float-count before output */
    }
    if( func == STITCH ){
        /* output file is twice the size */
        outpixels = 2*pixels;
        outimage = calloc(outpixels,sizeof(float));
    }
    if( func == MAXPOOL || func == AVGBOX ){
        /* param is the down-sampling box size */
        if(xosize <= 0 ) xosize = (int) xsize/param;
        if(yosize <= 0 ) yosize = (int) ysize/param;
        if(zosize <= 0 ) zosize = (int) zsize/param;
        if(xosize <= 0 ) xosize = 1;
        if(yosize <= 0 ) yosize = 1;
        if(zosize <= 0 ) zosize = 1;

        /* output file is smaller */
        outpixels = xosize*yosize*zosize;
        outimage = calloc(outpixels,sizeof(float));
        count    = calloc(outpixels,sizeof(int));
    }

    /* from here on outpixels should be fixed */
    if( func == SWAPXY ){
        /* output has swapped x and y axes */
        xosize = ysize;
        yosize = xsize;
    }
    if( func == GAUSSBLUR || func == MAXRADIUS || func == AVGRADIUS || func == MEDRADIUS ){
        /* param is the neighbor radius */
        if(xosize <= 0 ) xosize = (int) xsize;
        if(yosize <= 0 ) yosize = (int) ysize;
        if(zosize <= 0 ) zosize = (int) zsize;
        if(xosize <= 0 ) xosize = 1;
        if(yosize <= 0 ) yosize = 1;
        if(zosize <= 0 ) zosize = 1;

        if(radius<0) radius = 1.0;
        if(xrad<0) xrad = (int) radius;
        if(yrad<0) yrad = (int) radius;
        if(zrad<0) zrad = (int) radius;
        if( xrad > xosize-1 ) xrad = xosize-1;
        if( yrad > yosize-1 ) yrad = yosize-1;
        if( zrad > zosize-1 ) zrad = zosize-1;
        printf("filter radii: xrad=%d yrad=%d zrad=%d\n",xrad,yrad,zrad);

        /* output file is same size */
        outpixels = xosize*yosize*zosize;
        outimage = calloc(outpixels,sizeof(float));
        count    = calloc(outpixels,sizeof(int));
    }
    if( func == GAUSSBLUR ){
        radius = xrad;
        if(yrad > radius) radius = yrad;
        if(zrad > radius) radius = zrad;
        fwidthsq = param*param/radius/radius/log(2);
        printf("fractional width = %f  %f\n",sqrt(fwidthsq),radius);
    }
    if( func == MEDRADIUS ){
        /* allocate enough to hold the sphere to be median-ed */
        scratch = calloc((2*xrad+1)*(2*yrad+1)*(2*zrad+1),sizeof(int));
    }
    if( func == SEGMENT ){
        /* allocate enough to hold all segments */
        joins = calloc(pixels,sizeof(int));
        segmap = calloc(pixels,sizeof(int));
    }
    /* make sure these get set */
    if(xosize <= 0 && yosize <= 0 && zosize <= 0 && xsize > 0 && ysize > 0 && zsize > 0 ) {
        xosize = xsize; yosize=ysize ; zosize=zsize; 
    }
    if(xosize <= 0 ) xosize = (int) xsize;
    if(yosize <= 0 ) yosize = (int) ysize;
    if(zosize <= 0 ) zosize = (int) zsize;
    if(xosize <= 0 ) xosize = 1;
    if(yosize <= 0 ) yosize = 1;
    if(zosize <= 0 ) zosize = 1;
    printf("output array dimensions: x %d y %d z %d\n",xosize,yosize,zosize);
    outpixels = xosize*yosize*zosize;

    
    /* see if we need to allocate memory for output image */
    if( outimage == NULL ) {
        printf("allocating %ld pixels\n",pixels);
        outimage = calloc(outpixels+10,sizeof(float));
    }

    /* i is input  j is output */
    for(i=0;i<pixels;++i)
    {
        j=i;
        if(debug > 2) printf("inimage1[%d] = %f\n",i,inimage1[i]);

        /* skip any invalid values, propagate to output */
        for(k=1;k<=ignore_values;++k)
        {
            if(inimage1[i]==ignore_value[k]){
                outimage[j] = ignore_value[k];
                ++invalid_pixel[j];
                /* no need to check others */
                break;
            }

            if(infile2name == NULL) continue;
            if(inimage2[i]==ignore_value[k]){
                outimage[j] = ignore_value[k];
                ++invalid_pixel[j];
                /* no need to check others */
                break;
            }
        }
        ++valid_pixels;
    }

    for(k=0;k<samples;++k)
    {
        /* remap for functions that treat data in pairs */
        if( func == AB2PHIF ){
            i=j=k*2;
            a = inimage1[i];
            b = inimage1[i+1];

            F = sqrtf(a*a+b*b);
            phi = atan2(b,a);
            outimage[j] = F;
            outimage[j+1] = phi;
        }
        if( func == PHIF2AB ){
            i=j=k*2;
            F = inimage1[i];
            phi = inimage1[i+1];
            a = F*cos(phi);
            b = F*sin(phi);
            outimage[j] = a;
            outimage[j+1] = b;
        }

        if( func == ODD ){
            /* skip even pixels, for stats */
            i = 2*k+1;
            j = k;
            outimage[j] = inimage1[i];
        }
        if( func == EVEN ){
            /* skip odd pixels, for stats */
            i = 2*k;
            j = k;
            outimage[j] = inimage1[i];
        }
        if( func == EVENODD ){
            /* even-number output from first file, odd-number from second */
            i = k;
            j = 2*k;
            if(map_param) param=inimage2[i];
            outimage[j]   = inimage1[i];
            outimage[j+1] = param;
        }
        if( func == ODDEVEN ){
            /* even-number output from second file, odd-number from first */
            i = k;
            j = 2*k;
            if(map_param) param=inimage2[i];
            outimage[j]   = param;
            outimage[j+1] = inimage1[i];
        }
    }

    for(i=0;i<pixels;++i)
    {
        j=i;

        /* up to 3D array mappings */
        x = ( i % xsize ) + xstart;
        y = (( i / xsize ) % ysize ) + ystart;
        z = i / xsize / ysize + zstart;
        xo = x - xstart;
        yo = y - ystart;
        zo = z - zstart;
        j = xo + xosize*yo + zo*xosize*yosize;

        if(j > outpixels) continue;

        /* skip most calcs for invalid pixels */
        if(invalid_pixel[j]) {
            if( ! ( func == MEDRADIUS || func == MAXRADIUS || func == AVGRADIUS || func == AVGBOX || func == MAXPOOL ) ) continue;
        }

        if( func == SQRT ){
            outimage[j] = sqrtf(inimage1[i]);
        }
        if( func == CBRT ){
            outimage[j] = cbrtf(inimage1[i]);
        }
        if( func == CEIL ){
            outimage[j] = ceilf(inimage1[i]);
        }
        if( func == FLOOR ){
            outimage[j] = floorf(inimage1[i]);
        }
        if( func == FABS ){
            outimage[j] = fabsf(inimage1[i]);
        }
        if( func == SIGN ){
            outimage[j] = 2.0*(inimage1[i] >= 0.0)-1.0;
        }
        if( func == SIN ){
            outimage[j] = sinf(inimage1[i]);
        }
        if( func == ASIN ){
            outimage[j] = asinf(inimage1[i]);
        }
        if( func == SINH ){
            outimage[j] = sinhf(inimage1[i]);
        }
        if( func == ASINH ){
            outimage[j] = asinhf(inimage1[i]);
        }
        if( func == COS ){
            outimage[j] = cosf(inimage1[i]);
        }
        if( func == ACOS ){
            outimage[j] = acosf(inimage1[i]);
        }
        if( func == COSH ){
            outimage[j] = coshf(inimage1[i]);
        }
        if( func == ACOSH ){
            outimage[j] = acoshf(inimage1[i]);
        }
        if( func == TAN ){
            outimage[j] = tanf(inimage1[i]);
        }
        if( func == ATAN ){
            outimage[j] = atanf(inimage1[i]);
        }
        if( func == TANH ){
            outimage[j] = tanhf(inimage1[i]);
        }
        if( func == ATANH ){
            outimage[j] = atanhf(inimage1[i]);
        }
        if( func == ERF ){
            outimage[j] = erff(inimage1[i]);
        }
        if( func == ERFC ){
            outimage[j] = erfcf(inimage1[i]);
        }
        if( func == EXP ){
            outimage[j] = expf(inimage1[i]);
        }
        if( func == SAFEXP ){
            float v;
            v = inimage1[i];
            if(v > 3.4e38) v = 3.4e38;
            if(v < -3.4e38) v = -3.4e38;
            outimage[j] = expf(v);
        }
        if( func == LOG ){
            outimage[j] = logf(inimage1[i]);
        }
        if( func == LOG10 ){
            outimage[j] = log10f(inimage1[i]);
        }
#ifndef __APPLE__
        if( func == J0 ){
            outimage[j] = j0f(inimage1[i]);
        }
        if( func == J1 ){
            outimage[j] = j1f(inimage1[i]);
        }
        if( func == JN ){
            if(map_param) param=inimage2[i];
            outimage[j] = jnf((int) param,inimage1[i]);
        }
        if( func == Y0 ){
            outimage[j] = y0f(inimage1[i]);
        }
        if( func == Y1 ){
            outimage[j] = y1f(inimage1[i]);
        }
        if( func == YN ){
            if(map_param) param=inimage2[i];
            outimage[j] = ynf((int) param,inimage1[i]);
        }
        if( func == GAMMA ){
            outimage[j] = gammaf(inimage1[i]);
        }
#endif
        if( func == LGAMMA ){
            outimage[j] = lgammaf(inimage1[i]);
        }
        if( func == POW ){
            if(map_param) param=inimage2[i];
            outimage[j] = pow(inimage1[i],param);
        }
        if( func == ERFPOW ){
            outimage[j] = (float) pow(erf(inimage1[i]),param);
        }
        if( func == NORM ){
            outimage[j] = erff(inimage1[i]/sqrt(2.0))/2.0+0.5;
        }
        if( func == URAND ){
            outimage[j] = ran1(&seed);
        }
        if( func == GRAND ){
            outimage[j] = gaussdev(&seed);
        }
        if( func == LRAND ){
            outimage[j] = lorentzdev(&seed);
        }
        if( func == PRAND ){
            outimage[j] = poidev(fabsf(inimage1[i]),&seed);
        }
        if( func == ERAND ){
            outimage[j] = expdev(&seed);
        }
        if( func == TRAND ){
            outimage[j] = triangledev(&seed);
        }
        if( func == SET ){
            outimage[j] = param;
        }
        if( func == ADD ){
            if(map_param) param=inimage2[i];
            outimage[j] = inimage1[i]+param;
        }
        if( func == SUBTRACT ){
            if(map_param) param=inimage2[i];
            outimage[j] = inimage1[i]-param;
        }
        if( func == MULTIPLY ){
            if(map_param) param=inimage2[i];
            outimage[j] = inimage1[i]*param;
        }
        if( func == DIVIDE ){
            if(map_param) param=inimage2[i];
            outimage[j] = inimage1[i]/param;
            if(debug) printf("divide: %d %f / %f = %f\n",i,inimage1[i],param,outimage[j]);
        }
        if( func == INVERSE ){
            if(map_param) param=inimage2[i];
            outimage[j] = param/inimage1[i];
        }
        if( func == NEGATE ){
            if(map_param) param=inimage2[i];
            outimage[j] = param-inimage1[i];
        }
        if( func == MAXIMUM ){
            if(map_param) param=inimage2[i];
            outimage[j] = inimage1[i];
            if(outimage[j]<param) outimage[j]=param;
        }
        if( func == MINIMUM ){
            if(map_param) param=inimage2[i];
            outimage[j] = inimage1[i];
            if(outimage[j]>param) outimage[j]=param;
        }
        if( func == NANZERO ){
            outimage[j] = inimage1[i];
            if(isnan(outimage[j]) || isinf(outimage[j])) outimage[j]=0.0;
        }
        if( func == THRESH ){
            if(map_param) param=inimage2[i];
            outimage[j] = 1.0;
            if(inimage1[i]<param) outimage[j]=0.0;
        }
        if( func == SWAB4 ){
             swapper.floaty = inimage1[i];
             swappee.byte[0] = swapper.byte[3];
             swappee.byte[1] = swapper.byte[2];
             swappee.byte[2] = swapper.byte[1];
             swappee.byte[3] = swapper.byte[0];
             outimage[j] = swappee.floaty;
        }
        if( func == SWAB2 ){
             swapper.floaty = inimage1[i];
             swappee.byte[0] = swapper.byte[1];
             swappee.byte[1] = swapper.byte[0];
             swappee.byte[2] = swapper.byte[3];
             swappee.byte[3] = swapper.byte[2];
             outimage[j] = swappee.floaty;
        }
        if( func == REVERSE ){
            j = pixels-i-1;
            outimage[j]   = inimage1[j];
        }
        if( func == STRIDE ){
            k = i/xsize;
            if( k % 2 == ( param < 0.0 ) )
            {
                outimage[j] = inimage1[i];
                ++j;
            }
        }
        if( func == STITCH ){
            k = i/xsize;
            j = k*2*xsize+( i % xsize );
            outimage[j]       = inimage1[i];
            outimage[j+xsize] = inimage2[i];
        }
        if( func == FLIP ){
            x = i % xsize;
            y = i / xsize;
            j = y*xsize+(xsize-x-1);
            outimage[j] = inimage1[i];
        }
        if( func == FLOP ){
            x = i % xsize;
            y = i / xsize;
            j = (ysize-y-1)*xsize+x;
            outimage[j] = inimage1[i];
        }
        if( func == FLIPFLOP ){
            x = i % xsize;
            y = i / xsize;
            j = x*ysize+y;
            if(debug) printf("flipflop: %d %d -> %d %d   %f\n",x,y,i,j,inimage1[j]);
            outimage[j] = inimage1[i];
        }
        if( func == SWAPXY ){
            x = i % xsize;
            y = i / xsize;
            z = i / xsize / ysize;
            xo = y;
            yo = x;
            zo = z;
            j = xo + xosize*yo + zo*xosize*yosize;
            if(debug )printf("swapxy: %d %d -> %d %d   %f\n",x,y,i,j,inimage1[i]);
            outimage[j] = inimage1[i];
        }
        if( func == MAXPOOL ){
            x = i % xsize;
            y = ( i / xsize ) % ysize;
            z = i / xsize / ysize;
            xo = (int) x / param;
            yo = (int) y / param;
            zo = (int) z / param;
            j = xo + xosize*yo + zo*xosize*yosize;
            if( inimage1[i] > outimage[j]) outimage[j] = inimage1[i];
            if(debug)printf("maxpool: %d %d %d -> %d %d   %f %f\n",x,y,z,i,j,inimage1[i],outimage[j]);
        }
        if( func == AVGBOX ){
            x = i % xsize;
            y = ( i / xsize ) % ysize;
            z = i / xsize / ysize;
            xo = (int) x / param;
            yo = (int) y / param;
            zo = (int) z / param;
            j = xo + xosize*yo + zo*xosize*yosize;
            outimage[j] += inimage1[i];
            ++count[j];
            if(debug)printf("avgbox: %d %d %d -> %d %d   %f %f\n",x,y,z,i,j,inimage1[i],outimage[j]);
        }
        if( func == GAUSSBLUR ){
            /* skip over zeroes for speed */
            if(inimage1[i] == 0.0) continue;
            x = i % xsize;
            y = ( i / xsize ) % ysize;
            z = i / xsize / ysize;
            for ( dz = -zrad ; dz <= zrad; ++dz ) 
            for ( dy = -yrad ; dy <= yrad; ++dy ) 
            for ( dx = -xrad ; dx <= xrad; ++dx ) {
              fx = dx*dx; fy = dy*dy; fz = dz*dz;
              fx /= (xrad+1)*(xrad+1);
              fy /= (yrad+1)*(yrad+1);
              fz /= (zrad+1)*(zrad+1);
              fradsq = fx+fy+fz;
              if( fradsq > 1.0 ) continue;
              xo = x + dx;
              yo = y + dy;
              zo = z + dz;
              if( xo < 0 || xo >= xosize ) continue;
              if( yo < 0 || yo >= yosize ) continue;
              if( zo < 0 || zo >= zosize ) continue;
              j = xo + xosize*yo + zo*xosize*yosize;
              outimage[j] += inimage1[i]*exp(-fradsq/fwidthsq);
              ++count[j];
            }
            if(debug)printf("gaussblur: %d %d %d -> %d %d   %f %f\n",x,y,z,i,j,inimage1[i],outimage[j]);
        }
        if( func == EDGE || func == FEATHER ){
            int up=0, down=0, same=0;
            xo = i % xosize;
            yo = ( i / xosize ) % yosize;
            zo = i / xosize / yosize;
            for ( dz = -1 ; dz <= 1; ++dz ) 
            for ( dy = -1 ; dy <= 1; ++dy ) 
            for ( dx = -1 ; dx <= 1; ++dx ) {
              x = xo + dx;
              y = yo + dy;
              z = zo + dz;
              if( x < 0 || x >= xsize ) continue;
              if( y < 0 || y >= ysize ) continue;
              if( z < 0 || z >= zsize ) continue;
              j = x + xsize*y + z*xsize*ysize;
              if( inimage1[i] > inimage1[j] ) {
                 /* vote for going uphill */
                 ++up;
              }
              if( inimage1[i] < inimage1[j] ) {
                 /* vote for going downhill */
                 ++down;
              }
              if( inimage1[i] == inimage1[j] ) {
                 /* vote for exact same */
                 ++same;
              }
            }
            /* three possible outputs */
            if( up > down ) outimage[i] = 1.0;
            if( up < down ) outimage[i] = -1.0;
            if( up == down ) outimage[i] = 0.0;
            
            if( func == FEATHER ) {
               outimage[i] = inimage1[i]-0.25*outimage[i];
            }

            if(debug)printf("up: %d down: %d same: %d -> %f at: %d %d %d\n",up,down,same,outimage[i],x,y,z);
        }
        if( func == SEGMENT ){
            /* anything zero is always its own segment */
            if( inimage1[i] == 0.0 ) continue;
            int segment = 0;
            int othersegment = 0;
            /* look at all neighboring pixels */
            xo = i % xosize;
            yo = ( i / xosize ) % yosize;
            zo = i / xosize / yosize;
            for ( dz = -1 ; dz <= 1; ++dz ) 
            for ( dy = -1 ; dy <= 1; ++dy ) 
            for ( dx = -1 ; dx <= 1; ++dx ) {
              x = xo + dx;
              y = yo + dy;
              z = zo + dz;
              if( x < 0 || x >= xsize ) continue;
              if( y < 0 || y >= ysize ) continue;
              if( z < 0 || z >= zsize ) continue;
              j = x + xsize*y + z*xsize*ysize;
              /* consider value of neighbor */
              othersegment = (int) outimage[j];
              if( othersegment > 0.0 ) {
                 /* must have already been assigned to a segment */
                 if ( segment == 0.0 ) {
                   /* we havent discovered a segment yet, so inherit it */
                   segment = othersegment;
                 }
                 if( segment != othersegment ) {
                   /* hmm, a conflict */
                   /* this should never happen */
                   if(k > pixels) printf("PANIC: segment %d > pixels (%d) \n",othersegment,pixels);
                   /* build a mapping of segments that are connected */
                    int lo = segment < othersegment ? segment : othersegment;
  
                    int hi = segment < othersegment ? othersegment : segment;
                    /* only record if not already pointing somewhere better */
                    if( joins[hi] == 0 || joins[hi] > lo ) joins[hi] = lo;
                    if( njoins < hi ) njoins = hi;
                    segment = lo;  
                    /* continue scan using the lower id */
                 }
              }
              if(debug > 1)printf("segment: %d %d %d -> %d %d   %d   %f %f\n",x,y,z,i,j,segment,inimage1[i],outimage[j]);
            }
            if( segment == 0.0) {
              /* we must be the first? */
              ++nsegments;
              segment = nsegments;
              if(debug)printf("new segment %d at: %d %d %d\n",segment,x,y,z);
            }
            /* label the output image with this assigned segment */
            outimage[i] = (float) segment;
            if(debug)printf("output at: %d %d %d = %f\n",x,y,z,outimage[i]);
        }
        if( func == MAXRADIUS ){
            xo = i % xosize;
            yo = ( i / xosize ) % yosize;
            zo = i / xosize / yosize;
            for ( dz = -zrad ; dz <= zrad; ++dz ) 
            for ( dy = -yrad ; dy <= yrad; ++dy ) 
            for ( dx = -xrad ; dx <= xrad; ++dx ) {
              fx = dx*dx; fy = dy*dy; fz = dz*dz;
              fx /= (xrad+1)*(xrad+1);
              fy /= (yrad+1)*(yrad+1);
              fz /= (zrad+1)*(zrad+1);
              fradsq = fx+fy+fz;
              if( fradsq > 1.0 ) continue;
              x = xo + dx;
              y = yo + dy;
              z = zo + dz;
              if( x < 0 || x >= xsize ) continue;
              if( y < 0 || y >= ysize ) continue;
              if( z < 0 || z >= zsize ) continue;
              j = x + xsize*y + z*xsize*ysize;
              if( inimage1[i] > outimage[j]) outimage[j] = inimage1[i];
              if(debug)printf("maxradius: %d %d %d -> %d %d  %d   %f %f\n",x,y,z,i,j,parami,inimage1[i],outimage[j]);
            }
        }
        if( func == AVGRADIUS ){
            xo = i % xosize;
            yo = ( i / xosize ) % yosize;
            zo = i / xosize / yosize;
            for ( dz = -zrad ; dz <= zrad; ++dz ) 
            for ( dy = -yrad ; dy <= yrad; ++dy ) 
            for ( dx = -xrad ; dx <= xrad; ++dx ) {
              fx = dx*dx; fy = dy*dy; fz = dz*dz;
              fx /= (xrad+1)*(xrad+1);
              fy /= (yrad+1)*(yrad+1);
              fz /= (zrad+1)*(zrad+1);
              fradsq = fx+fy+fz;
              if( fradsq > 1.0 ) continue;
              x = xo + dx;
              y = yo + dy;
              z = zo + dz;
              if( x < 0 || x >= xsize ) continue;
              if( y < 0 || y >= ysize ) continue;
              if( z < 0 || z >= zsize ) continue;
              j = x + xsize*y + z*xsize*ysize;
              outimage[j] += inimage1[i];
              ++count[j];
            }
            if(debug)printf("avgradius: %d %d %d -> %d %d   %f %f\n",xo,yo,zo,i,j,inimage1[i],outimage[j]);
        }
        if( func == MEDRADIUS ){
            j=i;
            // speedup, if pixel already done
            if( count[j] ) continue;
            
            xo = j % xosize;
            yo = ( j / xosize ) % yosize;
            zo = j / xosize / yosize;
            narr=0;
            for ( dz = -zrad ; dz <= zrad; ++dz ) 
            for ( dy = -yrad ; dy <= yrad; ++dy ) 
            for ( dx = -xrad ; dx <= xrad; ++dx ) {
              fx = dx*dx; fy = dy*dy; fz = dz*dz;
              fx /= (xrad+1)*(xrad+1);
              fy /= (yrad+1)*(yrad+1);
              fz /= (zrad+1)*(zrad+1);
              fradsq = fx+fy+fz;
              if( fradsq > 1.0 ) continue;
              x = xo + dx;
              y = yo + dy;
              z = zo + dz;
              if( x < 0 || x >= xsize ) continue;
              if( y < 0 || y >= ysize ) continue;
              if( z < 0 || z >= zsize ) continue;
              k = x + xsize*y + z*xsize*ysize;
              if( invalid_pixel[k] ) continue;
              ++narr;
              scratch[narr] = inimage1[k];
              //printf("debug: %d %d %d -> %d %d %d   %f %f\n",xo,yo,zo,x,y,z,inimage1[i],outimage[i]);
            }
            if(narr) {
              outimage[j] = fmedian(narr,scratch);
            } else {
              /* no better ideas? */
              outimage[j] = inimage1[i];
            }

            /* speed up if this applies to whole row */
            if( xosize && xosize <= 2*xrad ) {
              for(xo=0;xo<=xosize;++xo){
                j = xo + xosize*yo + zo*xosize*yosize;
                outimage[j]=outimage[i];
                ++count[j];
              }
            }
            if( yosize && yosize <= 2*yrad ) {
              xo = i % xosize;
              yo = ( i / xosize ) % yosize;
              zo = i / xosize / yosize;
              for(yo=0;yo<yosize;++yo){
                j = xo + xosize*yo + zo*xosize*yosize;
                outimage[j]=outimage[i];
                ++count[j];
              }
            }
            if(debug)printf("medradius: %d %d %d -> %d %d   %f %f\n",xo,yo,zo,i,j,inimage1[i],outimage[i]);
        }
    }
//    if( func == AB2PHIF || func == PHIF2AB || func == EVENODD || func == ODDEVEN || func == STITCH )
//    {
//        pixels = 2*pixels;
//    }
    if( func == STRIDE )
    {
        pixels = j;
    }

    /* any post-processing */
    if( func == AVGBOX || func == AVGRADIUS ) {
        for(j=0;j<outpixels;++j)
        {
            outimage[j] /= count[j];
        }
        free(count);
    }

    if( func == SEGMENT ) {
        int segment = 0;
        int othersegment = 0;
        int minseg = 0;
        printf("%d segment discoveries and %d joinings\n",nsegments,njoins);
        for(k=1;k<=nsegments;++k)
        {
          if( ! joins[k] ) joins[k]=k;
          segment = minseg = k;
          if(debug) printf("segment %d joins with %d\n",k,joins[k]);
          while ( ! ( joins[segment] == 0 || joins[segment] == segment ) ) {
            /* not at the end of a chain yet */
            othersegment = joins[segment];
            if( othersegment < minseg ) minseg = othersegment;
            if( segmap[segment]>0 && segmap[segment] < minseg ) minseg = segmap[segment];
            if( segmap[othersegment]>0 && segmap[othersegment] < minseg ) minseg = segmap[othersegment];
            if(debug) printf("minseg=%d \n",minseg);
            segmap[segment] = minseg;
            segmap[othersegment] = minseg;
            if(debug) printf("mapping segment %d to %d\n",segment,minseg);\
            joins[segment] = minseg;
           segment = othersegment;
           }
          if(debug) printf(" min=%d\n",minseg);
          if(debug) printf("asigning segment %d to segment %d\n",segmap[k],k);
        }
        for(k=1;k<=nsegments;++k)
        {
            segmap[k]=0;
        }
        /* now re-number */
        int count = 0;
        for(k=1;k<=nsegments;++k)
        {
            segment = joins[k];
            if(! segment) segment = k;
            if(! segmap[segment]) 
            {
                ++count;
                segmap[segment]=count;
            }
        }
        nsegments = count;
        for(j=0;j<outpixels;++j)
        {
            segment = (int) outimage[j];
            outimage[j] = (float) segmap[joins[segment]];
        }

        for(k=0;k<nsegments;++k)
        {
            if(debug) printf("debug: k=%d joins[%d]=%d segmap[%d]=%d\n",k,k,joins[k],k,segmap[k]);
        }
    }

    /* now do stats */
    sum = sumsq = sumd = sumdsq = 0.0;
    min=max=outimage[0];
    for(j=0;j<outpixels;++j)
    {
        if(outimage[j]>max) max=outimage[j];
        if(outimage[j]<min) min=outimage[j];
        sum += outimage[j];
        sumsq += outimage[j]*outimage[j];
    }
    avg = sum/outpixels;
    rms = sqrt(sumsq/outpixels);
    if(ignore_values) printf("%d invalid of ",outpixels-valid_pixels);
    printf("%d pixels ",outpixels);
    if(ignore_values) printf("( %d left)",valid_pixels);
    printf("\n");
    for(j=0;j<outpixels;++j)
    {
        sumd   += outimage[j] - avg;
        sumdsq += (outimage[j] - avg) * (outimage[j] - avg);
    }
    rmsd = sqrt(sumdsq/outpixels);
    printf("max = %g min = %g\n",max,min);
    printf("mean = %g rms = %g rmsd = %g\n",avg,rms,rmsd);

    printf("writing %s as a %d-byte header and %u %u-byte floats\n",outfilename,outheader,outpixels,sizeof(float));
    outfile = fopen(outfilename,"wb");
    if(outheader)
    {
        /* copy header from first file */
        if(header > outheader)
        {
            /* truncate the original header */
            fwrite(headerstuff,outheader,1,outfile);
        }
        else
        {
            /* write full original header */
            fwrite(headerstuff,header,1,outfile);
            /* pad the rest with zeroes */
            j = header;
            k = 0;
            while(j < outheader)
            {
                fwrite(&k,1,1,outfile);
                ++j;
            }
        }
    }
    fwrite(outimage,outpixels,sizeof(float),outfile);
    fclose(outfile);

    return 0;
}



float poidev(float xm, long *idum)
{
    float gammln(float xx);
    float ran1(long *idum);
    /* oldm is a flag for whether xm has changed since last call */
    static float sq,alxm,g,oldm=(-1.0);
    float em,t,y;
        
    if (xm < 12.0) {
        /* use direct method: simulate exponential delays between events */
        if(xm != oldm) {
            /* xm is new, compute the exponential */
            oldm=xm;
            g=exp(-xm);
        }
        /* adding exponential deviates is equivalent to multiplying uniform deviates */
        /* final comparison is to the pre-computed exponential */
        em = -1;
        t = 1.0;
        do {
            ++em;
            t *= ran1(idum);
        } while (t > g);
    } else {
        /* Use rejection method */
        if(xm != oldm) {
            /* xm has changed, pre-compute a few things... */
            oldm=xm;
            sq=sqrt(2.0*xm);
            alxm=log(xm);
            g=xm*alxm-gammln(xm+1.0);
        }
        do {
            do {
                /* y is a deviate from a lorentzian comparison function */
                y=tan(M_PI*ran1(idum));
                /* shift and scale */
                em=sq*y+xm;
            } while (em < 0.0);                /* there are no negative Poisson deviates */
            /* round off to nearest integer */
            em=floor(em);
            /* ratio of Poisson distribution to comparison function */
            /* scale it back by 0.9 to make sure t is never > 1.0 */
            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
        } while (ran1(idum) > t);
    }
        
    return em;
}


/* return gaussian deviate with rms=1 and FWHM = 2/sqrt(log(2)) */
float gaussdev(long *idum)
{
    float ran1(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
        
    if (iset == 0) {
        /* no extra deviats handy ... */
        
        /* so pick two uniform deviates on [-1:1] */
        do {
            v1=2.0*ran1(idum)-1.0;
            v2=2.0*ran1(idum)-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0);
        /* restrained to the unit circle */
        
        /* apply Box-Muller transformation to convert to a normal deviate */
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;                /* we now have a spare deviate */
        return v2*fac;
    } else {
        /* there is an extra deviate in gset */
        iset=0;
        return gset;
    }
}


/* generate Lorentzian deviate with FWHM = 2 */
float lorentzdev(long *seed) {
    float ran1(long *idum);
 
    return tan(M_PI*(ran1(seed)-0.5));
}

/* return triangular deviate with FWHM = 1 */
float triangledev(long *seed) {
    float ran1(long *idum);
    float value;

    value = ran1(seed);
    if(value > 0.5){
        value = sqrt(2*(value-0.5))-1;
    }else{
        value = 1-sqrt(2*value);
    }

    return value;
}



float expdev(long *idum)
{
    float dum;
    
    do
    dum=ran1(idum);
    while( dum == 0.0);
    return -log(dum);
}



/* ln of the gamma function */
float gammln(float xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser = 1.000000000190015;
    for(j=0;j<=5;++j) ser += cof[j]/++y;
    
    return -tmp+log(2.5066282746310005*ser/x);
}





/* returns a uniform random deviate between 0 and 1 */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <= 0 || !iy) {
        /* first time around.  don't want idum=0 */
        if(-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        
        /* load the shuffle table */
        for(j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if(*idum < 0) *idum += IM;
            if(j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    /* always start here after initializing */
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}




#define SWAP(a,b) swp_temp=(a);(a)=(b);(b)=swp_temp

float fmedian(unsigned long n, float arr[])
{
    unsigned int i,j,k,l,ir,mid;
    float a,swp_temp;

    l=1;
    ir=n;
    k=(n+1)/2;

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
          i=l+1;    // initialize pointers for partitioning
          j=ir;    
          a=arr[l];    // partitioning element
          for(;;)    // innermost loop
          {
            do i++; while(arr[i]<a);    // scan up to find element > a
            do j--; while(arr[j]>a);    // scan down to find element < a
            if( j < i ) break;        // pointers crossed, median is in between
            SWAP(arr[i],arr[j]);
          }
          arr[l]=arr[j];            // insert partitioning element
          arr[j]=a;
          if( j >= k ) ir=j-1;        // Keep partition that contains the median active
          if( j <= k ) l=i;
      }
    }
}




/********************************************************************
*
*        Fourier()        Numerical Recipes's Fast Fourier Transform
*
*********************************************************************
*
*   local variables:
*        swp_temp                - used by SWAP macro
*
*        i1, i2                        - indicies used in bit-reversal
*        data_length                - size (in floats) of the data buffer
*        new_FFT_len                - length (in cplx) of the next sub-FFT
*        last_FFT_size         - size (in floats) of the last sub-FFT
*        sub_FFT_spacing - spacing between adjacent sub-FFTs
*
*   FFT_idx                 - index along sub-FFT
*        evn_idx                        - index of even-data Fourier coefficient
*        odd_idx                        - index of odd-data Fourier coefficient
*
*        w_temp
*        w_real, w_imag        - the complex value of the basis wave
*        w_p_real        - previous ''
*   w_p_imag
*        theta                        - argument to complex exponential
*        temp_real                - temporary complex number
*        temp_imag
*
*        Description:
*                This function computes the Fast Fourier Transform of the
*        complex data represented in data[].  Odd indicies (starting
*        with 1) are the real values, and even indicies (starting with
*        2) are imaginaries.  Note the Fortranish array indexing.
*                The Fourier spectrum is returned as a series of similar
*        complex coefficients.  The first pair ([1],[2]) is the coefficient
*        of zero-frequency.  Higher and higher positive frequencies are
*        stored in the higher indexed pairs.  The Nyquist frequency
*        coefficient (which is the same for both positive and negative
*        frequencies) is stored at data[length+1]. Progressively lower
*        (more positive) negative frequencies are entered until
*        data[length+2];
*                It is not recommended to call this particular function
*        directly, but if you do, be sure to make length a power of
*        two, and pass your data pointer as data-1.
*
********************************************************************/

float *Fourier(float *data, unsigned long length, int direction)
{
        float swp_temp;

        unsigned long i1, i2, temp_int;

        unsigned long data_length;
        unsigned long new_FFT_len, last_FFT_size, sub_FFT_spacing;
        unsigned long FFT_idx, evn_idx, odd_idx;

        double w_temp, w_real, w_imag, w_p_real, w_p_imag, theta;
        double temp_real, temp_imag;

        /* data size is 2 * complex numbers */
        data_length = 2*length;
        i2 = 1;

        /* reorganize data in bit-reversed order */
        for(i1 = 1; i1 < data_length; i1 += 2)
        {
                /* swap if appropriate */
                if (i2 > i1)
                {
                        SWAP(data[i2], data[i1]);
                        SWAP(data[i2+1], data[i1+1]);
                }

                /* calculate bit-reverse of index1 in index2 */
                temp_int = data_length/2;
                while ((temp_int >= 2)&&(i2 > temp_int))
                {
                        i2 -= temp_int;
                        temp_int >>= 1;
                }
                i2 += temp_int;
        }

        /* FFTs of length 1 have been "calculated" and organized so
           that the odd and even Fourier coefficents are grouped */

        /* first sub-FFT is of length 2 */
        last_FFT_size = (new_FFT_len = 2);

        /* for as long as sub-FFT is not full FFT */
        while(data_length > new_FFT_len)
        {
                /* separation of previous sub-FFT coefficients */
                sub_FFT_spacing = 2*last_FFT_size;

                /* this is a trig recurrence relation that will yield
                   the W number in the Danielson-Lanczos Lemma */
                theta = direction*(2*M_PI/new_FFT_len);
                w_temp = sin(0.5*theta);

                w_p_real = -2.0*w_temp*w_temp;
                w_p_imag = sin(theta);

                /* initialize W number */
                w_real = 1.0;
                w_imag = 0.0;

                /* recursively combine the sub-FFTs */
                for (FFT_idx = 1; FFT_idx < last_FFT_size; FFT_idx += 2)
                {
                        for (evn_idx = FFT_idx;
                                 evn_idx <= data_length;
                                 evn_idx += sub_FFT_spacing)
                        {
                                /* FFT coefficients to combine */
                                odd_idx = evn_idx + last_FFT_size;

                                /* calculate W*F(odd) */
                                temp_real = w_real*data[odd_idx] - w_imag*data[odd_idx+1];
                                temp_imag = w_real*data[odd_idx+1] + w_imag*data[odd_idx];

                                /* complete Lemma: F = F(even) + W*F(odd) */
                                /* since F(even) and F(odd) are pereodic... */
                                data[odd_idx] = data[evn_idx] - temp_real;
                                data[odd_idx+1] = data[evn_idx+1] - temp_imag;

                                /* "regular D/L here" */
                                data[evn_idx] += temp_real;
                                data[evn_idx+1] += temp_imag;
                        }
                        /* now calculate the next W */
                        w_temp = w_real;
                        w_real = w_temp*w_p_real - w_imag*w_p_imag + w_real;
                        w_imag = w_imag*w_p_real + w_temp*w_p_imag + w_imag;
                }
                /* prepare to combine next FFT */
                new_FFT_len = (last_FFT_size = sub_FFT_spacing);
        }

        return data;
}

