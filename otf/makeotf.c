/* A new small project to radially average OTF starting from */
/* PSF dataset. Mainly for structured illumination microscopy ( */
/* 1 or 2 objectives). Also take into account the following:    */
/*  compensation for finite bead size; */
/*  the phase factor associated with the side bands; */
/*  estimate the sub-pixel position of the bead by fitting parabolas; */
/*  cleanup out-of-band noises. */


/* from Lin:
I compile it with gcc with the Priism library (-limlib) and 
single-precision FFTW2 libraries (-lsrfftw and -lsfftw). No other dependencies.

A couple of option flags need to be specified, for example:
./otf /Users/talley/Dropbox/528_150912_1515_glyc_a2_001.dv test1.otf -angle -1.855500 -ls 0.2075 -na 1.4 -nimm 1.515 -fixorigin 3 20 -leavekz 7 11 3
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <IMInclude.h>
#include <sfftw.h>
#include <srfftw.h>

#define PI 3.1415926536
#define MAXDIRS 3
#define SPOTRATIO 0.1
/* #define NIMM 1.515 */
#define MAXPHASES 10
#define DEFAULTPHASES 5
#define MAXPATH 80

#define mag2(b) ((b).re*(b).re+(b).im*(b).im)

typedef fftw_complex complex;
typedef struct { double re; double im;} zomplex;
typedef struct { float x; float y;} vector;

float cmag(complex a)
{
  return (sqrt(a.re*a.re+a.im*a.im));
}

complex cmult(complex a, complex b)
{
  complex c;
  c.re = a.re*b.re - a.im*b.im;
  c.im = a.re*b.im + a.im*b.re;
  return(c);
}


float fitparabola(float a1, float a2, float a3);
void estimate_background(float *image, int nx, int ny, int border_size, float *background);
void determine_center_and_background(float **stack5phases, float *I2Mstack, int nx, int ny, int nz, int nphases, float *xc, float *yc, float *zc, float *background, float *background_i2m, float *xc_i2m, float *yc_i2m, float *zc_i2m, int two_lens, int I2M_inc);
void apodize(int napodize,int nx, int ny, float *image);
void cosapodize(int nx, int ny, float *image);
void makematrix (int nphases, float ** sepMatrix);
void separate(int nx,int ny,int z,int nphases, float *floatimage[], float ** sepMatrix);

void combine_reim(complex **bands, int norders, int nx, int nz, int bForcedPIshift);
void shift_center(complex *bands, int nx, int ny, int nz, float xc, float yc, float zc);

void beadsize_compensate(complex **bands, float k0angle, float linespacing, float bead_diameter, int norders,
                         int nx, int ny, int nz, float dr, float dz);
double sphereFFT(float k, float radius);

void radialft(complex *bands, int nx, int ny, int nz, complex *avg);
void cleanup(complex *otfkxkz, int order, int nx, int nz, float dkr, float dkz, float linespacing,
             int lambdanm, int icleanup, int ileave1_1, int ileave1_2, int ileave2, int twolens, float NA, float NIMM);
void cleanup_I2M(complex *otfkxkz, int nx, int nz, float dkr, float dkz, int lamdanm, int icleanup, float NA, float NIMM);
void modify(complex *otfkxkz, int nx, int nz, int *ifixz, int ifixr, int order, int twolens);
void fixorigin(complex *otfkxkz, int nx, int nz, int kx1, int kx2);
void rescale(complex *otfkxkz, int order, int nx, int nz, float *scalefactor, int dorescale);

void outputdata(int ostream_no, complex **bands, IW_MRC_HEADER *header, int norders, int nx, int ny, int nz, float dkr, float dkz, int five_bands);

void getbg_and_slope(char *corrfiles, float *background, float *slope, int nx, int ny);

void mrc_file_write(float *buffer, int nx, int ny, int nz, float rlen, float zlen, int mode, int iwave, char *files);
int commandline(int argc, char *argv[], int * twolens, int *rescale, float *beaddiam, float *k0angle, float *linespacing, int *five_bands, int *nphases, int *interpkr, int *leavekz, int *do_compen, int *I2M, char *I2Mfiles, float *background, int *bBgInExtHdr, int *order0gen, char *order0files, int *conjugate, float *na, float *nimm, int *ifixz, int *ifixr, int *bUsecorr, int *bForcedPIshift, char *corrfiles, char *ifiles, char *ofiles, int * ifilein, int *ofilein);

int main(argc,argv)
     int argc;
     char **argv;
{
  char ifiles[MAXPATH], ofiles[MAXPATH], I2Mfiles[MAXPATH], order0files[MAXPATH];
  int istream_no=1, ostream_no=2;
  int ixyz[3], mxyz[3], pixeltype;      /* variables for IMRdHdr call */
  float min, max, mean;                 /* variables for IMRdHdr call */
  int nx, ny, nz, nxy, phase, nphases=DEFAULTPHASES, norders, order, napodize=10;
  int i, j, z, border_size;
  float dr, dz, dkr, dkz, *background, background_dummy=-1., background_i2m, xcofm, ycofm, zcofm, xcofm_i2m, ycofm_i2m, zcofm_i2m;
  float **floatimage, **sepMatrix, *floatsection, *I2M_image=0;
  float *buffer;
  complex **bands, *avg_output[MAXPHASES], *I2Mavg_output=0;
  float bead_diameter=0.12, scalefactor=1.0;
  float k0angleguess=0, linespacing = /* 0.1945 after 11/22/02 for Rhodamine channel. 0.17536 for FITC channel. otherwise */ 0.20365 ;
  rfftwnd_plan rfftplan3d;
  IW_MRC_HEADER header, otfheader;
  int ifixz[3], ifixr, twolens=0, dorescale=0, ifilein=0, ofilein=0, five_bands=0, do_compen=1, I2M_inc=0, Generate_band0=0, conjugate=0;
  int bBgInExtHdr = 0; /* if the background of each section is recorded in the extended header's 3rd float number (in Lin's scope) */
  int zsec;
  int icleanup, ileavekz[3]={0,0,0}, ileave1_1, ileave1_2, ileave2, interpkr[2];
  float NA=1.4, NIMM=1.515;
  int bUseCorr=0;
  int bForcedPIshift=0;
  char corrfiles[MAXPATH];
  float *background2D=0, *slope2D=0;

  IMAlPrt(0);

  interpkr[0] = 0;
  interpkr[1] = 0;
  ifixz[0] = 1;ifixz[1]=1; ifixz[2]=1; ifixr = 0;
  if (!commandline(argc, argv, &twolens, &dorescale, &bead_diameter, &k0angleguess, &linespacing, &five_bands, &nphases, interpkr, ileavekz, &do_compen, &I2M_inc, I2Mfiles, &background_dummy, &bBgInExtHdr, &Generate_band0, order0files, &conjugate, &NA, &NIMM, ifixz, &ifixr, &bUseCorr, &bForcedPIshift, corrfiles, ifiles, ofiles, &ifilein, &ofilein))
    exit(0);

  if (nphases > MAXPHASES) {
    fprintf(stderr, "nphases is larger than MAXPHASES\n");
    exit(-1);
  }

  if (!ifilein) {
    printf(" PSF dataset file name: ");
    fgets(ifiles, MAXPATH, stdin);
    ifiles[strlen(ifiles)-1]='\0';
  }
  if (IMOpen(istream_no, ifiles, "ro")) {
    fprintf(stderr, "File %s does not exist.\n", ifiles);
    exit(-1);
  }

  if (!ofilein) {
    printf(" Output OTF file name: ");
    fgets(ofiles, MAXPATH, stdin);
    ofiles[strlen(ofiles)-1]='\0';
  }
  if (IMOpen(ostream_no, ofiles, "new")) {
    fprintf(stderr, "File %s can not be created.\n", ofiles);
    exit(-1);
  }

  IMRdHdr(istream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(istream_no, &header);

  nx=header.nx;
  ny=header.ny;
  if (!I2M_inc)
    nz=header.nz/nphases;
  else
    nz=header.nz/(nphases+1);

  dr = header.ylen;
  dz = header.zlen;

  dkr = 1/(ny*dr);
  dkz = 1/(nz*dz);

  printf("nx=%d, ny=%d, nz=%d\n", nx, ny, nz);

  norders = (nphases+1)/2;

  sepMatrix = (float **) malloc(nphases * sizeof(float *));
  for (i=0;i<nphases;i++)
    sepMatrix[i] = (float *) malloc(nphases * sizeof(float));
  makematrix(nphases, sepMatrix);

  nxy=(nx+2)*ny;

  buffer =  (float *) malloc((nx*ny) * sizeof(float));
  floatimage = (float **) malloc(nphases*sizeof(float *));
  bands = (complex **)malloc(nphases*sizeof(complex *));

  for(i=0;i<nphases;i++) {
    floatimage[i] = (float *) calloc(nxy*nz, sizeof(float));
    bands[i] = (complex *)floatimage[i];  /* reuse same allocation to save memory */
  }

  if (I2M_inc)
    I2M_image = (float *) calloc(nxy*nz, sizeof(float));

  background = calloc(nz*nphases, sizeof(float));
  border_size = 20;

  if (bUseCorr) {  // flatfield correction of measured data using calibration data
    printf("loading CCD calibration file\n");
    background2D = (float *) malloc((nx*ny) * sizeof(float));
    slope2D = (float *) malloc((nx*ny) * sizeof(float));
    getbg_and_slope(corrfiles, background2D, slope2D, nx, ny);
  }

  printf("Reading data...\n\n");
  zsec = 0;
  for(z=0; z<nz; z++) {
    for( phase=0; phase<nphases; phase++) {
      IMRdSec(istream_no, buffer);

      if(buffer[(ny-1)*nx + (nx-1)] < 0)  /* fix camera error at even binnings */
        buffer[(ny-1)*nx + (nx-1)] = buffer[(ny-1)*nx + (nx-2)];

      if (bBgInExtHdr) {
        int extInts;
        float extFloats[3];
        IMRtExHdrZWT(istream_no, zsec, 0, 0, &extInts, extFloats);
        background[zsec] = extFloats[2];
      }
      else if (bUseCorr) {
        int l, k;
        for (l=0; l<ny; l++)
          for (k=0; k<nx; k++) {
            buffer[l*nx+k] -= background2D[l*nx+k];
            buffer[l*nx+k] *= slope2D[l*nx+k];
          }
      }
      else if (background_dummy >=0)
        background[zsec] = background_dummy;
      else
        estimate_background(buffer, nx, ny, border_size, background+zsec);
      printf("%.3f\n",background[zsec]);

      for(i=0;i<ny;i++)
        for(j=0;j<nx;j++)
          floatimage[phase][z*nxy + i*(nx+2) + j] = buffer[i*nx + j];
      floatsection = &floatimage[phase][z*nxy];
      if(napodize>=0)
        apodize(napodize,nx,ny,floatsection);
      else if(napodize==-1)
        cosapodize(nx,ny,floatsection);

      zsec ++;
    } /* end for(phase) */
    if (I2M_inc) {  /* Data contains one extra section of I2M image */
      IMRdSec(istream_no, buffer);
      if(buffer[(ny-1)*nx + (nx-1)] == -1.0)  /* fix camera error at even binnings */
        buffer[(ny-1)*nx + (nx-1)] = buffer[(ny-1)*nx + (nx-2)];
      for(i=0;i<ny;i++)
        for(j=0;j<nx;j++)
          I2M_image[z*nxy + i*(nx+2) + j] = buffer[i*nx + j];
    }
  } /* end for(z) */

  /* Before FFT, use center band to estimate bead center position */
  determine_center_and_background(floatimage, I2M_image, nx, ny, nz, nphases, &xcofm, &ycofm, &zcofm, &background_dummy, &background_i2m, &xcofm_i2m, &ycofm_i2m, &zcofm_i2m, twolens, I2M_inc);

  printf("Center of mass is (%.3f, %.3f, %.3f)\n\n", xcofm, ycofm, zcofm);

  if (I2M_inc) {
    printf("I2M psf's background is %.3f\n", background_i2m);
    printf("I2M psf's center of mass is (%.3f, %.3f, %.3f)\n\n", xcofm_i2m, ycofm_i2m, zcofm_i2m);
  }

  for(z=0; z<nz; z++) {
    for(i=0;i<ny;i++)
      for(j=0;j<nx;j++)
        for( phase=0; phase<nphases; phase++)
          floatimage[phase][z*nxy + i*(nx+2) + j] -= background[z*nphases+phase];

    separate(nx, ny, z, nphases, floatimage, sepMatrix);
    if (I2M_inc)
      for(i=0;i<ny;i++)
        for(j=0;j<nx;j++)
          I2M_image[z*nxy + i*(nx+2) + j] -= background_i2m;
  }

  if (Generate_band0) {
    mrc_file_write(floatimage[0], nx+2, ny, nz, dr, dz, 0, header.iwav1, order0files);
  }
  rfftplan3d = rfftw3d_create_plan(nz, ny, nx, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  printf("Before fft\n");
  for( phase=0; phase<nphases; phase++)
    rfftwnd_one_real_to_complex(rfftplan3d, floatimage[phase], NULL);
  if (I2M_inc)
    rfftwnd_one_real_to_complex(rfftplan3d, I2M_image, NULL);

  fftwnd_destroy_plan(rfftplan3d);
  printf("After fft\n\n");

  /* modify the phase of bands, so that it corresponds to FFT of a bead at origin */
  printf("Shifting center...\n");
  for (phase=0; phase<nphases; phase++)
    shift_center(bands[phase], nx, ny, nz, xcofm, ycofm, zcofm);

  if (I2M_inc)
    shift_center((complex *)I2M_image, nx, ny, nz, xcofm_i2m, ycofm_i2m, zcofm_i2m);

  /* Now compensate finite bead size for all bands */
  if (do_compen)
    beadsize_compensate(bands, k0angleguess, linespacing, bead_diameter, norders, nx, ny, nz, dkr, dkz);

  for (phase=0; phase < nphases; phase ++)
    avg_output[phase] = (complex *) calloc(nz*(nx/2+1), sizeof(complex));
  if (!avg_output[nphases-1]) {
    printf("No memory available before radialft()\n");
    exit(-1);
  }

  for (order=0; order < norders; order ++) {
    if (order==0)
      radialft(bands[0], nx, ny, nz, avg_output[0]);

    else {
      radialft(bands[2*order-1], nx, ny, nz, avg_output[2*order-1]);
      radialft(bands[2*order], nx, ny, nz, avg_output[2*order]);
    }
  }

  if (I2M_inc) {
    I2Mavg_output = (complex *) calloc(nz*(nx/2+1), sizeof(complex));
    radialft((complex *)I2M_image, nx, ny, nz, I2Mavg_output);
  }

  icleanup = nx/2+1;
  ileave1_1 = ileavekz[0]; ileave1_2 = ileavekz[1]; ileave2 = ileavekz[2];

  for (phase=0; phase < nphases; phase ++) {
    modify(avg_output[phase], nx, nz, ifixz, ifixr, (phase+1)/2, twolens);
    if (ileave1_1 > 0)
      cleanup(avg_output[phase], (phase+1)/2, nx, nz, dkr, dkz, linespacing, header.iwav1, icleanup, ileave1_1, ileave1_2, ileave2, twolens, NA, NIMM);
    if (phase==0 && interpkr[0] > 0)
      fixorigin(avg_output[phase], nx, nz, interpkr[0], interpkr[1]);
    rescale(avg_output[phase], (phase+1)/2, nx, nz, &scalefactor, dorescale);
  }

  /* For side bands, combine bandre's and bandim's into bandplus */
  /* Shouldn't this be done later on the averaged bands? */
  if (!five_bands)
    combine_reim(avg_output, norders, nx, nz, bForcedPIshift);

  if (conjugate)
    for (phase=0; phase<nphases; phase++)
      for (i=0; i<(nx/2+1)*nz; i++)
        avg_output[phase][i].im *= -1;

  otfheader = header;

  outputdata(ostream_no, avg_output, &otfheader, norders, nx, ny, nz, dkr, dkz, five_bands);

  if (I2M_inc) {
    /*     cleanup_I2M(I2Mavg_output, nx, nz, dkr, dkz, header.iwav1, icleanup, NA, NIMM); */
    fixorigin(I2Mavg_output, nx, nz, interpkr[0], interpkr[1]);
    rescale(I2Mavg_output, 0, nx, nz, &scalefactor, dorescale);
    mrc_file_write((float *)I2Mavg_output, nz, nx/2+1, 1, dkz, dkr, 2, header.iwav1, I2Mfiles);
  }

  return 0;
}

/***************************** makematrix ************************************/
/*     generates the matrix that is to be used to separate the raw indata    */
/*     into the different bands of sample information.                       */
/*****************************************************************************/
void makematrix (int nphases, float ** sepMatrix)
{
  int j,order,norders;
  float phi;

  norders = (nphases+1)/2;
  phi = 2*PI/nphases;
  for(j=0;j<nphases;j++) {
    sepMatrix[0][j] = 1.0/nphases;
    for(order=1;order<norders;order++)
      {
        sepMatrix[2*order-1][j] = cos(j*order*phi)/nphases;
        sepMatrix[2*order  ][j] = sin(j*order*phi)/nphases;
      }
  }
}

/***************************** separate **************************************/
/*     Applies image arithmetic to the image sequence floatimage[][],        */
/*     which was acquired with different phases of the illumination          */
/*     pattern, to separate the different bands of sample information.       */
/*     The coefficients are pre-stored in the matrix sepMatrix.              */
/*     The bands are returned in the same array floatimage where the         */
/*     input data was.                                                       */
/*****************************************************************************/
void separate(int nx,int ny,int z,int nphases, float *floatimage[], float ** sepMatrix)
{
  int i,j,k,l, ind;
  float output[MAXPHASES];

  if (nphases > MAXPHASES) {
    fprintf(stderr, "In separate(), nphases is larger than MAXPHASES\n");
    exit(-1);
  }

  for(l=0;l<ny;l++) {
    for(k=0;k<nx;k++) {
      ind= z*(nx+2)*ny + l*(nx+2) + k;
      for(i=0;i<nphases;i++){
        output[i]=0.0;
        for(j=0;j<nphases;j++)
          output[i] += floatimage[j][ind] * sepMatrix[i][j];
      }
      for(i=0;i<nphases;i++)
        floatimage[i][ind] = output[i];
    }
  }
}

void estimate_background(float *image, int nx, int ny, int border_size, float *background)
{
  int i, j, total=0;
  double sum=0;

  for (i=0; i<ny; i++)
    if (i<border_size || i>ny-border_size)
      for (j=0; j<nx; j++)
        if (j<border_size || j>nx-border_size) {
          sum += image[i*nx+j];
          total ++;
        }
  *background = sum/total;
}


/*  locate peak pixel to subpixel accuracy by fitting parabolas  */
void determine_center_and_background(float **stack5phases, float *I2M_image, int nx, int ny, int nz, int nphases, float *xc, float *yc, float *zc, float *background, float *background_i2m, float *xc_i2m, float *yc_i2m, float *zc_i2m, int twolens, int I2M)
{
  int i, j, k, maxi, maxj, maxk, ind, nxy2, infocus_sec;
  int iminus, iplus, jminus, jplus, kminus, kplus;
  float maxval, reval, valminus, valplus, *stack3d;
  double sum;

  printf("In determine_center_and_background()\n");
  nxy2 = (nx+2)*ny;

  stack3d = calloc(nxy2*nz, sizeof(float));
  for (j=0; j<nxy2*nz; j++) {
    for (i=0; i<nphases; i++)
      stack3d[j] += stack5phases[i][j];
    stack3d[j] /= nphases;
  }

  /* Search for the peak pixel */
  /* Be aware that stack3d is of dimension (nx+2)xnyxnz */
  maxval=0.0;
  maxi=0; maxj=0; maxk=0;
  for(k=0;k<nz;k++)
    for(i=0;i<ny;i++)
      for(j=0;j<nx;j++) {
        ind=k*nxy2+i*(nx+2)+j;
        reval=stack3d[ind];
        if( reval > maxval ) {
          maxval = reval;
          maxi=i; maxj=j;
          maxk=k;
        }
      }

  iminus = maxi-1; iplus = maxi+1;
  if( iminus<0 ) iminus+=ny;
  if( iplus>=ny ) iplus-=ny;
  jminus = maxj-1; jplus = maxj+1;
  if( jminus<0 ) jminus+=nx;
  if( jplus>=nx ) jplus-=nx;
  kminus = maxk-1; kplus = maxk+1;
  if( kminus<0 ) kminus+=nz;
  if( kplus>=nz ) kplus-=nz;

  valminus = stack3d[kminus*nxy2+maxi*(nx+2)+maxj];
  valplus  = stack3d[kplus *nxy2+maxi*(nx+2)+maxj];
  *zc = maxk + fitparabola(valminus, maxval, valplus);

#ifdef ACCURATE_PEAK
  if (!twolens)
    *zc = maxk + fitparabola(valminus, maxval, valplus);
  else {
    float zcs[5], peakval[5], minval_right, minval_left;
    int n_samples = 5, kmin, kminus2, kplus2;
    float Xs[5], Ys[5];

    kminus2 = maxk-n_samples/2; kplus2 = maxk+n_samples/2;
    for (i=0; i<n_samples; i++) {
      Xs[i] = kminus2 + i;
      Ys[i] = stack3d[(kminus2+i)*nxy2 + maxi*(nx+2) + maxj];
    }
    zcs[0] = fitparabola_Nsample_points(Xs, Ys, n_samples, peakval);

    minval_right = Ys[n_samples-1];
    minval_left = Ys[0];

    /* search forward for the next trough */
    kmin = kplus2;
    for (k=kplus2+1; k<=kplus2+5; k++)
      if (stack3d[k*nxy2+maxi*(nx+2)+maxj] < minval_right) {
        minval_right = stack3d[k*nxy2+maxi*(nx+2)+maxj];
        kmin = k;
      }
    for (i=0; i<n_samples; i++) {
      Xs[i] = kmin + 1 + i;
      Ys[i] = stack3d[(kmin + 1 + i)*nxy2 + maxi*(nx+2) + maxj];
    }
    zcs[1] = fitparabola_Nsample_points(Xs, Ys, n_samples, peakval+1);

    /* search backward for the previous trough */
    kmin = kminus2;
    for (k=kminus2-1; k>=kminus-5; k--)
      if (stack3d[k*nxy2+maxi*(nx+2)+maxj] < minval_left) {
        minval_left = stack3d[k*nxy2+maxi*(nx+2)+maxj];
        kmin = k;
      }
    for (i=0; i<n_samples; i++) {
      Xs[i] = kmin - 1 - i;
      Ys[i] = stack3d[(kmin - 1 - i)*nxy2 + maxi*(nx+2) + maxj];
    }
    zcs[2] = fitparabola_Nsample_points(Xs, Ys, n_samples, peakval+2);

    printf("zc0=%.3f, zc_plus=%.3f, zc_minus=%3f\n", zcs[0], zcs[1], zcs[2]);

    zcs[3] = 94; peakval[3] = 251;
    zcs[4] = 68; peakval[4] = 304;
    *zc = fitparabola_Nsample_points(zcs, peakval, 3, peakval);
  }
#endif

  valminus = stack3d[maxk*nxy2+iminus*(nx+2)+maxj];
  valplus  = stack3d[maxk*nxy2+iplus *(nx+2)+maxj];
  *yc = maxi + fitparabola(valminus, maxval, valplus);

  valminus = stack3d[maxk*nxy2+maxi*(nx+2)+jminus];
  valplus  = stack3d[maxk*nxy2+maxi*(nx+2)+jplus];
  *xc = maxj + fitparabola(valminus, maxval, valplus);

  free(stack3d);

  sum = 0;
  infocus_sec = floor(*zc);
  for (k=0; k<nphases; k++)
    for (i=0; i<*yc-20; i++)
      for (j=0; j<nx; j++)
        sum += stack5phases[k][infocus_sec*nxy2 + i*(nx+2) + j];
  *background = sum / (nphases*(*yc-20)*nx);

  if (I2M) {
    sum = 0;
    for (i=0; i<*yc-20; i++)
      for (j=0; j<nx; j++)
        sum += I2M_image[infocus_sec*nxy2 + i*(nx+2) + j];
    *background_i2m = sum / (nx*(*yc-20));

    /* Search for the peak pixel */
    /* Be aware that I2M_image is of dimension (nx+2)xnyxnz */
    maxval=0.0;
    for(k=0;k<nz;k++)
      for(i=0;i<ny;i++)
        for(j=0;j<nx;j++) {
          ind=k*nxy2+i*(nx+2)+j;
          reval=I2M_image[ind];
          if( reval > maxval ) {
            maxval = reval;
            maxi=i; maxj=j;
            maxk=k;
          }
        }

    iminus = maxi-1; iplus = maxi+1;
    if( iminus<0 ) iminus+=ny;
    if( iplus>=ny ) iplus-=ny;
    jminus = maxj-1; jplus = maxj+1;
    if( jminus<0 ) jminus+=nx;
    if( jplus>=nx ) jplus-=nx;
    kminus = maxk-1; kplus = maxk+1;
    if( kminus<0 ) kminus+=nz;
    if( kplus>=nz ) kplus-=nz;

    valminus = I2M_image[kminus*nxy2+maxi*(nx+2)+maxj];
    valplus  = I2M_image[kplus *nxy2+maxi*(nx+2)+maxj];
    *zc_i2m = maxk + fitparabola(valminus, maxval, valplus);

    valminus = I2M_image[maxk*nxy2+iminus*(nx+2)+maxj];
    valplus  = I2M_image[maxk*nxy2+iplus *(nx+2)+maxj];
    *yc_i2m = maxi + fitparabola(valminus, maxval, valplus);

    valminus = I2M_image[maxk*nxy2+maxi*(nx+2)+jminus];
    valplus  = I2M_image[maxk*nxy2+maxi*(nx+2)+jplus];
    *xc_i2m = maxj + fitparabola(valminus, maxval, valplus);
  }
}

/***************************** fitparabola **********************************/
/*     Fits a parabola to the three points (-1,a1), (0,a2), and (1,a3).     */
/*     Returns the x-value of the max (or min) of the parabola.             */
/****************************************************************************/

float fitparabola( float a1, float a2, float a3 )
{
  float slope,curve,peak;

  slope = 0.5* (a3-a1);         /* the slope at (x=0). */
  curve = (a3+a1) - 2*a2;       /* (a3-a2)-(a2-a1). The change in slope per unit of x. */
  if( curve == 0 )
    {
      printf("no peak: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f\n",a1,a2,a3,slope,curve);
      return( 0.0 );
    }
  peak = -slope/curve;          /* the x value where slope = 0  */
  if( peak>1.5 || peak<-1.5 )
    {
      printf("bad peak position: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f, peak=%f\n",a1,a2,a3,slope,curve,peak);
      return( 0.0 );
    }
  return( peak );
}

/*********************************   apodize   **********************************************/
/*   softens the edges of a singe xy section to reduce edge artifacts and improve the fits  */
/********************************************************************************************/
void apodize(int napodize,int nx,int ny,float *image)
{
  float diff,fact;
  int k,l;

  for(k=0;k<nx;k++)
    {
      diff = (image[ (ny-1)*(nx+2) +k ] - image[ /* 0*(nx+2)+ */ k ]) / 2;
      for(l=0;l<napodize;l++)
        {
          fact = 1 - sin((((float)l+0.5)/napodize)*PI*0.5);
          image[l*(nx+2) + k] += diff*fact;
          image[(ny-1-l)*(nx+2) + k] -= diff*fact;
        }
    }
  for(l=0;l<ny;l++)
    {
      diff = (image[ l*(nx+2) + nx-1 ] - image[ l*(nx+2) /* +0 */ ]) / 2;
      for(k=0;k<napodize;k++)
        {
          fact = 1 - sin((((float)k+0.5)/napodize)*PI*0.5);
          image[l*(nx+2) + k] += diff*fact;
          image[l*(nx+2) + ny-1-k] -= diff*fact;
        }
    }
}


/*******************************   cosapodize   ********************************/
/*   softens the edges to reduce edge artifacts and improve the fits           */
/*******************************************************************************/
void cosapodize(int nx,int ny,float *image)
{
  float xfact,yfact;
  int k,l;

  printf("in cosapodize\n");
  for(k=0;k<nx;k++)
    {
      xfact = sin(PI*((float)k+0.5)/nx);
      for(l=0;l<ny;l++)
        {
          yfact = sin(PI*((float)l+0.5)/ny);
          image[l*(nx+2) + k] *= (xfact*yfact);
        }
    }
}

/* Find out the phase of the side bands relative to the center band    */
/* and get the real OTF of the side bands. It's based on the fact that */
/* bandre=BAND*cos(phi), and bandim=BAND*sin(-phi)                     */
void combine_reim(complex **otf, int norders, int nx, int nz, int bForcedPIshift)
{
  int order, i, nxz;
  complex czero={0,0};
  double bandre_mag, bandim_mag, phi;
  complex otfval;

  printf("In combine_reim()\n");
  nxz = (nx/2+1)*nz;

  /* for (order=1; order<norders; order++) { */
  /*   for (i=0; i<nxz; i++) { */
  /*     otfval.re = otf[2*order-1][i].re - otf[2*order][i].im; */
  /*     otfval.im = otf[2*order-1][i].im + otf[2*order][i].re; */
  /*     otf[2*order-1][i] = otfval; */
  /*     otf[2*order][i] = czero; */
  /*   } */
  /* } */

  for (order=1; order<norders; order++) {
    bandre_mag = 0;
    bandim_mag = 0;
    for (i=0; i<nxz; i++) {
      bandre_mag += cmag(otf[2*order-1][i]);
      bandim_mag += cmag(otf[2*order][i]);
    }
    phi = atan(bandim_mag/bandre_mag);
    /* phi is now in the first quadrant only, since bandim_mag and bandre_mag are both positive */
    /* which quadrant phi should be in is decided by the kz=0 plane of bandre and bandim values */
    if (otf[2*order-1][nx/8*nz].re < 0 && otf[2*order][nx/8*nz].re >0)
      phi += PI;
    else if (otf[2*order-1][nx/8*nz].re < 0 && otf[2*order][nx/8*nz].re <0)
      phi = PI-phi;
    else if (otf[2*order-1][nx/8*nz].re > 0 && otf[2*order][nx/8*nz].re >0)
      phi = -phi;

    /* Sometimes the above logic still won't find the correct phase; user can specify additional Pi phase shift*/
    if (order==1 && bForcedPIshift)
      phi = PI-phi;

    printf("  phi=%f\n", phi);

    for (i=0; i<nxz; i++) {
      otfval.re = otf[2*order-1][i].re * cos(phi) + otf[2*order][i].re * sin(-phi);
      otfval.im = otf[2*order-1][i].im * cos(phi) + otf[2*order][i].im * sin(-phi);
      otf[2*order-1][i] = otfval;
      otf[2*order][i] = czero;  /* since bandim = band*sin(phi) and phi is essentially 0  */
    }
  }
}

void beadsize_compensate(complex **bands, float k0angle, float linespacing, float bead_diameter, int norders, int nx, int ny, int nz, float dkr, float dkz)
{
  int order, kin, iin, jin, kout, iout, indin, kycent, kxcent, kzcent, nxy;
  complex *bandptr=0, *bandptr1=0;
  float rho; /* the distance in Fourier space from a point to the origin */
  float limit_at_origin;   /* the limit value of the FT of a sphere at the origin of Fourier space */
  float radius;   /* the radius of the fluorescent bead, according to the number provided by vendor */
  float kz, ky, kx, k0y, k0x, ratio, k0mag;

  printf("In beadsize_compensate()\n");
  kycent = ny/2;
  kxcent = nx/2;
  kzcent = nz/2;
  nxy = (nx/2+1)*ny;

  k0mag = 1/linespacing;
  radius = bead_diameter * 0.5;
  limit_at_origin = 4*PI*radius*radius*radius/3;

  for (order=0; order<norders; order++) {
    if(order==0)
      bandptr = bands[order];
    else {
      bandptr  = bands[2*order-1];
      bandptr1 = bands[2*order];
    }

    k0x = ((float)order)/(norders-1) * k0mag * cos(k0angle);
    k0y = ((float)order)/(norders-1) * k0mag * sin(k0angle);

    printf("order=%d, k0x=%f, k0y=%f\n", order, k0x, k0y);

    for (kin=0; kin<nz; kin++) {
      kout = kin;
      if (kout>kzcent) kout -= nz;
      kz = kout * dkz;
      for (iin=0; iin<ny; iin++) {
        iout = iin;
        if (iout>kycent) iout -= ny;
        ky = iout * dkr + k0y;
        for (jin=0; jin<kxcent+1; jin++) {
          kx = jin * dkr + k0x;
          indin = kin*nxy+iin*(nx/2+1)+jin;
          if (! (order==0 && indin==0) ) {
            rho = sqrt(kz*kz + ky*ky + kx*kx);
            ratio = sphereFFT(rho,radius) / limit_at_origin;
          }
          else
            ratio = 1;
          bandptr[indin].re /= ratio;
          bandptr[indin].im /= ratio;
          if (order > 0) {
            bandptr1[indin].re /= ratio;
            bandptr1[indin].im /= ratio;
          }
        }
      }
    }
  }
}

/* To get rid of checkerboard effect in the OTF bands */
/* (xc, yc, zc) is the estimated center of the point source, which in most cases is the bead */
/* Converted from Fortran code. kz is treated differently than kx and ky. don't know why */
void shift_center(complex *bands, int nx, int ny, int nz, float xc, float yc, float zc)
{
  int kin, iin, jin, indin, nxy, kz, kx, ky, kycent, kxcent, kzcent;
  complex exp_iphi;
  float phi1, phi2, phi, dphiz, dphiy, dphix;

  kycent = ny/2;
  kxcent = nx/2;
  kzcent = nz/2;
  nxy = (nx/2+1)*ny;

  dphiz = 2*PI*zc/nz;
  dphiy = 2*PI*yc/ny;
  dphix = 2*PI*xc/nx;

  for (kin=0; kin<nz; kin++) {    /* the origin of Fourier space is at (0,0) */
    kz = kin;
    if (kz>kzcent) kz -= nz;
    phi1 = dphiz*kz;      /* first part of phi */
    for (iin=0; iin<ny; iin++) {
      ky = iin;
      if (iin>kycent) ky -= ny;
      phi2 = dphiy*ky;   /* second part of phi */
      for (jin=0; jin<kxcent+1; jin++) {
        kx = jin;
        indin = kin*nxy+iin*(nx/2+1)+jin;
        phi = phi1+phi2+dphix*kx;  /* third part of phi */
        /* kz part of Phi has a minus sign, I don't know why. */
        exp_iphi.re = cos(phi);
        exp_iphi.im = sin(phi);
        bands[indin] = cmult(bands[indin], exp_iphi);
      }
    }
  }
}

void radialft(complex *band, int nx, int ny, int nz, complex *avg_output)
{
  int kin, iin, jin, indin, indout, indout_conj, kz, kx, ky, kycent, kxcent, kzcent;
  int *count, nxz, nxy;
  float rdist;

  printf("In radialft()\n");
  kycent = ny/2;
  kxcent = nx/2;
  kzcent = nz/2;
  nxy = (nx/2+1)*ny;
  nxz = (nx/2+1)*nz;

  count = (int *) calloc(nxz, sizeof(int));

  if (!count) {
    printf("No memory availale in radialft()\n");
    exit(-1);
  }

  for (kin=0; kin<nz; kin++) {
    kz = kin;
    if (kin>kzcent) kz -= nz;
    for (iin=0; iin<ny; iin++) {
      ky = iin;
      if (iin>kycent) ky -= ny;
      for (jin=0; jin<kxcent+1; jin++) {
        kx = jin;
        rdist = sqrt(kx*kx+ky*ky);
        if (rint(rdist) < nx/2+1) {
          indin = kin*nxy+iin*(nx/2+1)+jin;
          /* indout = floor(rdist)*nz+kin; /\* rint(rdist*nz)+kin caused trouble *\/  */
          /* the above line has been used up till 2010, which is probably less right compared to the following line: */
          indout = rint(rdist)*nz+kin;
          avg_output[indout].re += band[indin].re;
          avg_output[indout].im += band[indin].im;
          count[indout] ++;
        }
      }
    }
  }

  for (indout=0; indout<nxz; indout++) {
    if (count[indout]>0) {
      avg_output[indout].re /= count[indout];
      avg_output[indout].im /= count[indout];
    }
  }

  /* Then complete the rotational averaging and scaling*/
  for (kx=0; kx<nx/2+1; kx++) {
    indout = kx*nz+0;
    avg_output[indout].im = 0;
    for (kz=1; kz<=nz/2; kz++) {
      indout = kx*nz+kz;
      indout_conj = kx*nz + (nz-kz);
      avg_output[indout].re = (avg_output[indout].re + avg_output[indout_conj].re)/2;
      avg_output[indout].im = (avg_output[indout].im - avg_output[indout_conj].im)/2;
      avg_output[indout_conj] = avg_output[indout];
      avg_output[indout_conj].im *= -1;
    }
  }
  free(count);
}

void cleanup(complex *otfkxkz, int order, int nx, int nz, float dkr, float dkz, float linespacing,
             int lamdanm, int icleanup, int ileave1_1, int ileave1_2, int ileave2, int twolens, float NA, float NIMM)
{
  int ix, iz, kzstart, kzend, jotfshape;
  float lamda, sinalpha, cosalpha, kr, krmax, beta, kzedge, k0mag;
  complex czero={0.0,0.0};
  float NA_local;

  if (order==0)
    NA_local = NA;
  else
    NA_local = NA * .92;

  lamda = lamdanm * 0.001;
  sinalpha = NA_local/NIMM;
  cosalpha = cos(asin(sinalpha));
  krmax = 2*NA_local/lamda;
  k0mag = 1.0/linespacing;

  printf("nz=%d\n", nz);

  if (!twolens) {
    for (ix=0; ix<icleanup; ix++) {
      kr = ix * dkr;
      if ( kr <= krmax ) {
        beta = asin( ( NA_local - kr*lamda ) /NIMM );
        kzedge = (NIMM/lamda) * ( cos(beta) - cosalpha );
        /* printf("kzedge=%f, dkz=%f\n", kzedge, dkz); */
        if(order==0) {
          kzstart = rint((kzedge/dkz) + 1);
          kzend = nz - kzstart;
          for (iz=kzstart; iz<=kzend; iz++)
            otfkxkz[ix*nz+iz] = czero;
        }
        else if (order==1) {
          jotfshape = rint((kzedge/dkz) + 0.999);
          kzend = ileave1_1 - jotfshape;
          if (kzend>=0)
            for (iz=0; iz<=kzend; iz++) {
              if(iz==0)
                otfkxkz[ix*nz] = czero;
              else {
                otfkxkz[ix*nz+iz] = czero;
                otfkxkz[ix*nz+nz-iz] = czero;
              }
            }
          kzstart = ileave1_2 + jotfshape;
          for (iz=kzstart+1; iz<=nz/2/*+1*/; iz++) {
            otfkxkz[ix*nz+iz] = czero;
            otfkxkz[ix*nz+nz-iz] = czero;
          }
        }
        else { /* order == 2 */
          jotfshape = rint((kzedge/dkz) + 0.999);
          kzstart = ileave2 + jotfshape;
          for (iz=kzstart; iz<=nz/2/*+1*/; iz++) {
            otfkxkz[ix*nz+iz] = czero;
            otfkxkz[ix*nz+nz-iz] = czero;
          }
        }
      }
      else {   /* outside of lateral resolution limit */
        for (iz=0; iz<nz; iz++)
          otfkxkz[ix*nz+iz] = czero;
      }
    }
  }
  else {
    float lambdaem, lambdaexc, two_over_lambdaem, two_over_lambdaexc, alpha, beta, betamin;
    int center_of_arc;

    lambdaem = lamda/NIMM;
    lambdaexc = 0.88* lambdaem;;  /* 0.88 approximates a typical lambdaexc/lambdaem  */
    two_over_lambdaem = 2/lambdaem;
    two_over_lambdaexc = 2/lambdaexc;
    alpha = asin(NA_local/NIMM);  /* aperture angle of objectives */
    beta = asin(k0mag/two_over_lambdaexc);   /* angle of side illumination beams */
    betamin = asin((k0mag/two_over_lambdaexc) - sin(alpha)*SPOTRATIO);

    for (ix=0; ix<icleanup; ix++) {
      kr = ix * dkr;
      if ( kr <= krmax ) {
        if (order==0)
          center_of_arc = ceil(two_over_lambdaexc/dkz)+3;
        else if (order == 1) {
          center_of_arc = ceil((1/lambdaexc)*(cos(beta) + 1)/dkz)+3;
        }
        else { /* order == 2 */
          center_of_arc = ceil(two_over_lambdaexc * cos(beta)/dkz)+3;
        }
        kzstart = ceil(center_of_arc + sqrt(two_over_lambdaem * two_over_lambdaem - kr * kr)/dkz);
        kzend = nz - kzstart;
        for (iz=kzstart; iz<=kzend; iz++)
          otfkxkz[ix*nz+iz] = czero;
      }
      else
        for (iz=0; iz<nz; iz++)
          otfkxkz[ix*nz+iz] = czero;
    }
  }
}

void cleanup_I2M(complex *otfkxkz, int nx, int nz, float dkr, float dkz, int lamdanm, int icleanup, float NA, float NIMM)
{
  float lamda, sinalpha, cosalpha, krmax, two_over_lambdaem, beta, lamdaem, kr;
  int ix, iz, kzstart, kzend;
  complex czero={0.0,0.0};

  lamda = lamdanm * 0.001;
  lamdaem = lamda/NIMM;
  sinalpha = NA/NIMM;
  cosalpha = cos(asin(sinalpha));
  krmax = 2*NA/lamda;
  two_over_lambdaem = 2/lamdaem;

  for (ix=0; ix<icleanup; ix++) {
    kr = ix * dkr;
    if (kr < krmax) {
      /* clean up stuff outside the side detection bands */
      kzstart = ceil(sqrt(two_over_lambdaem * two_over_lambdaem - kr * kr)/dkz)+1;
      kzend = nz - kzstart;
      for (iz=kzstart; iz<=kzend; iz++)
        otfkxkz[ix*nz+iz] = czero;

      /* Then, clean up stuff in-between center and side band */
      beta = asin( ( NA - kr*lamdaem ) /NIMM );
      kzstart = ceil((NIMM/lamdaem) * ( cos(beta) - cosalpha )/dkz);
      kzend = floor(sqrt(two_over_lambdaem * two_over_lambdaem - krmax*krmax)/dkz)-2;
      for (iz=kzstart; iz<=kzend; iz++)
        otfkxkz[ix*nz+iz] = czero;
      for (iz=nz-kzend; iz<nz-kzstart; iz++)
        otfkxkz[ix*nz+iz] = czero;
    }
    else
      for (iz=0; iz<nz; iz++)
        otfkxkz[ix*nz+iz] = czero;
  }
}

void outputdata(int ostream_no, complex **bands, IW_MRC_HEADER *header, int norders, int nx, int ny, int nz, float dkr, float dkz, int five_bands)
{
  int i;

  printf("In outputdata()\n");

  header->nx = nz;
  header->ny = nx/2+1;
  if (!five_bands)
    header->nz = norders;
  else
    header->nz = norders*2-1;
  header->mode = IW_COMPLEX;
  header->xlen = dkz;
  header->ylen = dkr;
  header->zlen = 0;
  header->amin = 0;
  header->amax = 1;

  IMPutHdr(ostream_no, header);
  for (i=0; i<norders; i++)
    if (i==0)
      IMWrSec(ostream_no, bands[i]);
    else {
      IMWrSec(ostream_no, bands[2*i-1]);
      if (five_bands)
        IMWrSec(ostream_no, bands[2*i]);
    }
  IMWrHdr(ostream_no, header->label, 0, header->amin, header->amax, header->amean);
  IMClose(ostream_no);
}

/* According to Mats's derivation, FT of a sphere with radius R is f(k)=(R/Pi*k^2)*g(2*Pi*k*R), where g(x)=sin(x)/x - cos(x) */
/* The limit at the origin of F space is 4*Pi*R^3/3 */
double sphereFFT(float k, float radius)
{
  double a, x;

  if (k>0 || k<0) {
    x = 2*PI*radius*k;
    a = radius / (PI*k*k) * (sin(x)/x - cos(x));
  }
  else
    a = 4*PI*pow(radius, 3)/3;

  return a;
}

void modify(complex *otfkxkz, int nx, int nz, int *ifixz, int ifixr, int order, int twolens)
{
  int kz, kx, ind0, ind1, ind2;

  if (ifixz[0])
    for (kz=0; kz<nz; kz++)
      if (kz !=0 || order > 0)
        otfkxkz[kz] = otfkxkz[nz+kz];    /* replace the kx=0 line with kx=1 line */

  if (ifixr) {
    for (kx=ifixr; kx<nx/2+1; kx++) {
      ind0 = kx*nz;
      ind1 = ind0+1;
      ind2 = ind0+nz-1;
      otfkxkz[ind0].re = (otfkxkz[ind1].re + otfkxkz[ind2].re)/2;
      otfkxkz[ind0].im = (otfkxkz[ind1].im + otfkxkz[ind2].im)/2;
    }
  }
}

void fixorigin(complex *otfkxkz, int nx, int nz, int kx1, int kx2)
{
  float meani, slope, avg, lineval;
  int numvals, i, j;
  double *sum, totsum=0, ysum=0, sqsum=0;

  printf("In fixorigin()\n");
  meani = 0.5*(kx1+kx2);
  numvals = kx2-kx1+1;
  sum = (double *) malloc((kx2+1)*sizeof(double));

  for (i=0; i<=kx2; i++) {  /* '<=' is correct. For a long time, '<' was used, which is a bug. */
    sum[i] = 0;

    if (i == 0)
      sum[i] = otfkxkz[0].re;  /* don't want to add up garbages on kz axis */
    else
      for (j=0; j<nz; j++)
        sum[i] += otfkxkz[i*nz+j].re; /* imaginary parts will cancel each other because of symmetry */

    if (i>=kx1) {
      totsum += sum[i];
      ysum += sum[i] * (i-meani);
      sqsum += (i-meani) * (i-meani);
    }
  }
  slope = ysum / sqsum;
  avg = totsum / numvals;

  for (i=0; i<kx1; i++) {
    lineval = avg + (i-meani)*slope;
    otfkxkz[i*nz].re -= sum[i] - lineval;
  }

  free(sum);
}

void rescale(complex *otfkxkz, int order, int nx, int nz, float *scalefactor, int dorescale)
{
  int nxz, ind;
  float valmax=0, mag;

  nxz = (nx/2+1)*nz;
  if (order==0 || dorescale) { /* need to find out scalefactor */
    for (ind=0; ind<nxz; ind++) {
      mag = cmag(otfkxkz[ind]);
      if (mag > valmax)
        valmax = mag;
    }
    *scalefactor = 1/valmax;
  }
  for (ind=0; ind<nxz; ind++) {
    otfkxkz[ind].re *= *scalefactor;
    otfkxkz[ind].im *= *scalefactor;
  }
}


complex otfinterpolate(complex * otf, int kr, float krscale, int kz, float kzscale, int nzotf)
/* (kr, kz) is Fourier space coords with origin at kr=kz=0 and going  betwen -nx(or ny,nz)/2 and +nx(or ny,nz)/2 */
{
  fftw_complex otfval;

  int irindex, izindex, indices[2][2];
  float krindex, kzindex;
  float ar, az;

  krindex = kr * krscale;
  kzindex = kz * kzscale;

  irindex = floor(krindex);
  izindex = floor(kzindex);

  ar = krindex - irindex;
  az = kzindex - izindex;  // az is always 0 for 2D case, and it'll just become a 1D interp

  if (izindex == nzotf-1) {
    indices[0][0] = irindex*nzotf+izindex;
    indices[0][1] = irindex*nzotf+0;
    indices[1][0] = (irindex+1)*nzotf+izindex;
    indices[1][1] = (irindex+1)*nzotf+0;
  }
  else {
    indices[0][0] = irindex*nzotf+izindex;
    indices[0][1] = irindex*nzotf+(izindex+1);
    indices[1][0] = (irindex+1)*nzotf+izindex;
    indices[1][1] = (irindex+1)*nzotf+(izindex+1);
  }
  otfval.re = (1-ar)*(otf[indices[0][0]].re*(1-az) + otf[indices[0][1]].re*az) +
    ar*(otf[indices[1][0]].re*(1-az) + otf[indices[1][1]].re*az);
  otfval.im = (1-ar)*(otf[indices[0][0]].im*(1-az) + otf[indices[0][1]].im*az) +
    ar*(otf[indices[1][0]].im*(1-az) + otf[indices[1][1]].im*az);

  return otfval;
}

void downsize_pixel(complex *otfra[], complex *newotf[], int norders, int nkr, int nkz, int downfactor)
/*
   Downsize the k-space pixel size by downfactor (an integer), so that the number of pixels increases but
   the same k-space extent remains.
*/
{
  int nkr_new, nkz_new;
  int order, kr, kz;

  nkr_new = nkr*downfactor;
  nkz_new = nkz*downfactor;

  for (order=0; order <norders; order ++) {
    complex * otforiginal = otfra[order];
    int index_new = 0;
    for (kr=0; kr<nkr_new; kr++)
      for(kz=0; kz<nkz_new; kz++) {
        newotf[order][index_new] = otfinterpolate(otforiginal, kr, 1.0/downfactor, kz, 1.0/downfactor, nkz);
        index_new ++;
      }
  }
}

void Usage()
{

  /* not complete */
  printf("Usage:\nradialft [[input file] [output file]] [Options]\n");
  printf("\nOptions:\n");
  printf("\t-beaddiam F -- the diameter of the bead in microns (default is 0.12)\n");
  printf("\t-angle F -- the k0 vector angle with which the PSF is taken (default is 0)\n");
  printf("\t-ls x -- the illumination pattern's line spacing (in microns, default is 0.2)\n");
  printf("\t-nocompen -- do not perform bead size compensation (default is to perform)\n");
  printf("\t-5bands -- to output 5 OTF bands (default is combining higher-order's real and imag bands into one output\n");
  printf("\t-fixorigin kr1 kr2 -- the starting and end pixel for interpolation along kr axis (default is 2 and 9)\n");
  printf("\t-na x -- the (effective) NA of the objective\n");
  printf("\t-nimm x -- the index of refraction of the immersion liquid\n");
  printf("\t-leavekz kz1_1 kz1_2 kz2 -- the pixels to be retained on kz axis\n");
  printf("\t-I2M I2m_otf_file_name -- data contains I2M PSF and supply the I2M otf file name\n");
  printf("\t-background bkgd_value -- use the supplied number as the background to subtract\n");
  printf("\t-help or -h -- print this message\n");
}

#include <ctype.h>  /* isdigit() */

int getinteger(char * buffer, int * ret)
{
  int i, len;

  len = strlen(buffer);
  for (i=0; i<len; i++)
    if (!isdigit(buffer[i])) {
      if (i==0 && buffer[i] !='+' && buffer[i] !='-' )
        return 0;
      else if (i!=0)
        return 0;
    }
  sscanf(buffer, "%d", ret);
  return 1;
}

int getintegers(char *buffer[], int *ret, int num)
{
  int i;

  for (i=0; i<num; i++) {
    if (!getinteger(buffer[i], ret+i))
      return 0;
  }
  return 1;
}

int getfloat(char *buffer, float *ret)
{
  int i, len;
  int dotcount = 0;

  len = strlen(buffer);
  for (i=0; i<len; i++)
    if (!isdigit(buffer[i])) {
      if (i==0 && buffer[i] !='+' && buffer[i] !='-' && buffer[i]!='.' )
        return 0;
      else if (i!=0 && buffer[i] != '.')
        return 0;
      if (buffer[i] == '.')
        dotcount ++;
    }
  if (dotcount > 1) return 0;
  sscanf(buffer, "%f", ret);
  return 1;

}

int getfloats(char *buffer[], float *ret, int num)
{
  int i;

  for (i=0; i<num; i++) {
    if (!getfloat(buffer[i], ret+i))
      return 0;
  }
  return 1;
}

void getbg_and_slope(char *corrfiles, float *background, float *slope, int nx, int ny)
{
  int cstream_no=10;
  int ixyz[3], mxyz[3], pixeltype;      /* variables for IMRdHdr call */
  float min, max, mean;      /* variables for IMRdHdr call */
  IW_MRC_HEADER header;

  if (IMOpen(cstream_no, corrfiles, "ro")) {
    fprintf(stderr, "File %s does not exist\n", corrfiles);
    exit(-1);
  }

  IMRdHdr(cstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(cstream_no, &header);
  if (header.nx != nx || header.ny != ny) {
    fprintf(stderr, "calibration file %s has different dimension than data file", corrfiles);
    exit(-1);
  }
  IMAlCon(cstream_no, 0);

  IMRdSec(cstream_no, background);
  if (header.nz>1)  /* correction file with only 1 section is assumed to contain background image */
    IMRdSec(cstream_no, slope);
  else {
    int i;
    for (i=0; i<nx*ny; i++)
      slope[i] = 1.;
  }

  IMClose(cstream_no);
}

int commandline(int argc, char *argv[], int * twolens, int *rescale, float *beaddiam, float *k0angle, float *linespacing, int *five_bands, int *nphases, int *interpkr, int *leavekz, int *do_compen, int *I2M_inc,  char *I2Mfiles, float * background, int *bBgInExtHdr, int *order0gen, char *order0files, int *conjugate, float *na, float *nimm, int * ifixz, int *ifixr, int *bUsecorr, int *bForcedPIshift, char *corrfiles, char *ifiles, char *ofiles, int * ifilein, int *ofilein) {
  int ncomm=1;

  while (ncomm < argc) {
    if (strcmp(argv[ncomm], "-angle") == 0) {
      if (!getfloat(argv[ncomm+1], k0angle)) {
        printf("Invalid input for switch -angle\n");
        return 0;
      }
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-ls") == 0) {
      if (!getfloat(argv[ncomm+1], linespacing)) {
        printf("Invalid input for switch -ls\n");
        return 0;
      }
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-fixorigin") == 0) {
      if (!getintegers(argv+ncomm+1, interpkr, 2)) {
        printf("Invalid input for switch -fixorigin\n");
        return 0;
      }
      ncomm +=3;
    }
    else if (strcmp(argv[ncomm], "-leavekz") == 0) {
      if (!getintegers(argv+ncomm+1, leavekz, 3)) {
        printf("Invalid input for switch -leavekz\n");
        return 0;
      }
      ncomm +=4;
    }
    else if (strcmp(argv[ncomm], "-2lenses") == 0) {
      *twolens = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-rescale") == 0) {
      *rescale = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-5bands") == 0) {
      *five_bands = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-conj") == 0) {
      *conjugate = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-nocompen") == 0) {
      *do_compen = 0;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-I2M") == 0) {
      *I2M_inc = 1;
      strcpy(I2Mfiles, argv[ncomm+1]);
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-gen_order0") == 0) {
      *order0gen = 1;
      strcpy(order0files, argv[ncomm+1]);
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-nphases") == 0) {
      if (!getinteger(argv[ncomm+1], nphases)) {
        printf("Invalid input for switch -nphases\n");
        return 0;
      }
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-beaddiam") == 0) {
      if (!getfloat(argv[ncomm+1], beaddiam)) {
        printf("Invalid input for switch -beaddiam\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-background") == 0) {
      if (!getfloat(argv[ncomm+1], background)) {
        printf("Invalid input for switch -background\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-bgInExtHdr") == 0) {
      *bBgInExtHdr = 1;
      ncomm += 1;
    }
    else if (strcmp(argv[ncomm], "-na") == 0) {
      if (!getfloat(argv[ncomm+1], na)) {
        printf("Invalid input for switch -na\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-nimm") == 0) {
      if (!getfloat(argv[ncomm+1], nimm)) {
        printf("Invalid input for switch -nimm\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-ifixkz") == 0) {
      if (!getinteger(argv[ncomm+1], ifixz)) {
        printf("Invalid input for switch -ifixkz\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-ifixkr") == 0) {
      if (!getinteger(argv[ncomm+1], ifixr)) {
        printf("Invalid input for switch -ifixkr\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-usecorr") == 0) {
      *bUsecorr = 1;
      strcpy(corrfiles, argv[ncomm+1]);
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-PIshift") == 0) {
      *bForcedPIshift = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-help") == 0 || strcmp(argv[ncomm], "-h") == 0) {
      Usage();
      return 0;
    }
    else if (ncomm < 3) {
      if (ncomm==1) {
        strcpy(ifiles, argv[ncomm]);
        *ifilein = 1;
      }
      else if (ncomm==2) {
        strcpy(ofiles, argv[ncomm]);
        *ofilein = 1;
      }
      ncomm ++;
    }
    else {
      printf("Invalid command line option %s\n", argv[ncomm]);
      Usage();
      return 0;
    }
  }
  return 1;
}

void mrc_file_write(float *buffer, int nx, int ny, int nz, float rlen, float zlen, int mode, int iwave, char *files)
{
  int ostream_no=19;
  IW_MRC_HEADER header;
  int dimx, dimy, dimz, nxy, i;
  float amin=0, amax=1, amean=0.1; /* make-shift way of doing this */

  printf("Writing output file: %s\n", files);

  if (IMOpen(ostream_no, files, "new")) {
    fprintf(stderr, "File %s can not be created.\n", files);
    exit(-1);
  }

  dimx = nx;
  dimy = ny;
  dimz = nz;
  nxy = nx*ny;

  switch (mode) {
  case 0: /* normal floating-point data */
    header.mode = IW_FLOAT;
    break;
  case 1: /* half Fourier space complex data */
    header.mode = IW_COMPLEX;
    dimx = nx/2 + 1;
    nxy = (nx+2)*ny;
    break;
  case 2: /* full Fourier space complex data */
    header.mode = IW_COMPLEX;
    nxy = nx*ny*2;
    break;
  case 3: /* half Fourier floating-point data */
    header.mode = IW_FLOAT;
    dimx = nx/2+1;
    nxy = dimx * dimy;
    break;
  default:
    fprintf(stderr, "Illegal mode in mrc_file_write()\n");
    exit(-1);
  }

  header.nx = dimx;
  header.mx = dimx;
  header.ny = dimy;
  header.my = dimy;
  header.nz = dimz;
  header.mz = dimz;
  header.ylen = rlen;
  header.xlen = rlen;
  header.zlen = zlen;
  header.nxst = 0;
  header.nyst = 0;
  header.nzst = 0;
  header.num_waves = 1;
  header.num_times = 1;
  header.iwav1 = iwave;
  header.alpha = 0;
  header.beta = 0;
  header.gamma = 0;
  header.mapc = 1;
  header.mapr = 2;
  header.maps = 3;
  header.ispg = 0;
  header.nDVID = -16244;
  header.ntst = 0;
  header.inbsym = 0;
  header.nint = 0;
  header.nreal = 0;
  header.nres = 1;
  header.nzfact = 1;
  header.file_type = 0;  /* normal image type */
  header.lens = 0;
  header.interleaved = 2;
  header.nlab = 1;
  IMPutHdr(ostream_no, &header);

  for (i=0; i<nz; i++)
    IMWrSec(ostream_no, buffer + i*nxy);

  IMWrHdr(ostream_no, "Processed", 1, amin, amax, amean);
  IMClose(ostream_no);
}
