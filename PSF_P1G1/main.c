#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nr_utils.c"

#include "fft.c"
#include "mtfc.c"

#define INV_ARCSEC 206264.806247

#define N_APER 4

#define N_SUBPIX 8

typedef struct {
  double reff,g1,g2;
  double ra,rb,pa;
} GALAXY_PAR;

int main(int argc, char **argv) {

  long i,j,ii,jj,psfdim;
  double xctr,yctr,delx,dely,factor;
#if 0
  double theta, Xrot, Yrot; int irot;
#endif
  int k,l,lref;
  double u, v;
  APERTURE_TYPE A[N_APER];
#if 0
  int n_ab, m_ab;
#endif
  double obscure = 0., *pt, *pt2, temp, Rap;
  double **PSF, **Objtot;
  double pix, cd, pix_out, lambda, weight, lambda_min, lambda_max;
  long pad = 40*N_SUBPIX, Nx=1536, Nlam, ilam;
  unsigned long flags;
  GALAXY_PAR g;
  double w, u1, v1, Dx, Dy, ow;
  double mom[5], ZA[8], jit[3];
  double alpha, alphaplus;

  PSF = NULL;

  if (argc<24) {
    fprintf(stderr, "Error -- too few arguments\nUsage: ./main.x pupil lambdamin lambdamax Nlam hexflags reff g1 g2 ZA0..ZA7 jitxx jitxy jityy Dx Dy alpha alpha+\n");
    return(1);
  }

  /* Galaxy parameters */
  sscanf(argv[6], "%lg", &(g.reff));
  sscanf(argv[7], "%lg", &(g.g1));
  sscanf(argv[8], "%lg", &(g.g2));
  /* ... derived ... */
  g.ra = g.reff*exp(sqrt(g.g1*g.g1+g.g2*g.g2));
  g.rb = g.reff*exp(-sqrt(g.g1*g.g1+g.g2*g.g2));
  g.pa = atan2(g.g2,g.g1)/2;

  /* Amplitudes of each Zernike -- nm rms input, converted to m here */
  /* Ordering: focus; astig +,x; coma x,y; trefoil x,y; spherical */
  for(i=0;i<8;i++) {
    sscanf(argv[9+i], "%lg", ZA+i);
    ZA[i] /= 1e9;
  }

  /* Jitter in mas^2 */
  for(i=0;i<3;i++) {
    sscanf(argv[17+i], "%lg", jit+i);
  }

  /* Chromatic astrometric gradient, mas across bandpass */
  sscanf(argv[20], "%lg", &Dx);
  sscanf(argv[21], "%lg", &Dy);

  /* IPC parameters */
  sscanf(argv[22], "%lg", &alpha);
  sscanf(argv[23], "%lg", &alphaplus);

  /* Legal pupils:
   * 1 = AFTA - test
   */

  /* Flags: [hexadecimal digits]
   * digits 00-03 --> code for pixel parameters
   * digits 04    --> exclude smearing terms? (0=include them)
   * digits 05    --> exclude pixel tophat? (0=include it)
   */

  /* Wavelengths in microns */

  sscanf(argv[1], "%d", &lref); l=lref;
  sscanf(argv[2], "%lg", &lambda_min);  lambda_min *= 1e-6;
  sscanf(argv[3], "%lg", &lambda_max);  lambda_max *= 1e-6;
  sscanf(argv[4], "%ld", &Nlam);
  sscanf(argv[5], "%lx", &flags);

  Objtot = dmatrix(0,2*pad,0,2*pad);

  for(ilam=0;ilam<Nlam;ilam++) {

    lambda = lambda_min + (lambda_max-lambda_min)*(0.5+ilam)/Nlam;
    weight = 1./Nlam;

    /* Set basic parameters */
    pix = 0.11; cd = 0.03234*0.03234 + 0.025*0.025*0;
    pix_out = pix/N_SUBPIX;
    A[l].lambda = lambda;
    A[l].dx = lambda/pix_out*INV_ARCSEC/(4.*Nx);

    allocate_aperture(A+l,Nx);
    Rap=1e-2; obscure = 0.0;
    obscure = 0.32;
    set_aperture_mask(A+l, l?1.161*obscure:0, Rap=1.161);

    xctr = 0.5*A[l].dx*(A[l].Nx-1);
    yctr = 0.5*A[l].dx*(A[l].Nx-1);
    for(i=0;i<A[l].Nx;i++) for(j=0;j<A[l].Nx;j++) {
      delx = i*A[l].dx - xctr;
      dely = j*A[l].dx - yctr;
      /* Spider removal would go here*/
    }

    /* Insert chromatic astrometric gradient */
    ow = (ilam-(Nlam-1)/2.)/Nlam;
    shift_phase_Z(A+l, 1, 1, Dx*Rap/INV_ARCSEC/1e3*ow, Rap);
    shift_phase_Z(A+l, 1,-1, Dy*Rap/INV_ARCSEC/1e3*ow, Rap);

    /* To get Zernike mode with rms 1, use amplitude sqrt(2j_order+|m|+1) * sqrt(m==0? 1: 2),
     * where j_order = (n-|m|)/2.   ==> amplitude is sqrt(n+1) * sqrt(m==0? 1: 2).
     */
    shift_phase_Z(A+l, 2, 0, ZA[0]*sqrt(3.), Rap);
    shift_phase_Z(A+l, 2, 2, ZA[1]*sqrt(6.), Rap);
    shift_phase_Z(A+l, 2,-2, ZA[2]*sqrt(6.), Rap);
    shift_phase_Z(A+l, 3, 1, ZA[3]*sqrt(8.), Rap);
    shift_phase_Z(A+l, 3,-1, ZA[4]*sqrt(8.), Rap);
    shift_phase_Z(A+l, 3, 3, ZA[5]*sqrt(8.), Rap);
    shift_phase_Z(A+l, 3,-3, ZA[6]*sqrt(8.), Rap);
    shift_phase_Z(A+l, 4, 0, ZA[7]*sqrt(5.), Rap);

    generate_mtf(A+l);

    /* Insert smearing */
    for(i=0;i<A[l].Nu;i++) for(j=0;j<A[l].Nu;j++) {
      u = A[l].du * (i-A[l].Nu/2);
      v = A[l].du * (j-A[l].Nu/2);
      factor  = 1.;
      if ((flags&0x20)==0) {
        factor *= fabs(u)>1? sin(M_PI*u*pix/INV_ARCSEC)/(M_PI*u*pix/INV_ARCSEC): 1;
        factor *= fabs(v)>1? sin(M_PI*v*pix/INV_ARCSEC)/(M_PI*v*pix/INV_ARCSEC): 1;
      }
      if ((flags&0x10)==0) {
        factor *= exp(-2.*M_PI*M_PI*u*u/INV_ARCSEC/INV_ARCSEC*cd);
        factor *= exp(-2.*M_PI*M_PI*v*v/INV_ARCSEC/INV_ARCSEC*cd);
      }

      /* Jitter */
      factor *= exp(-2.*M_PI*M_PI*u*u/INV_ARCSEC/INV_ARCSEC*jit[0]/1e6);
      factor *= exp(-4.*M_PI*M_PI*u*v/INV_ARCSEC/INV_ARCSEC*jit[1]/1e6);
      factor *= exp(-2.*M_PI*M_PI*v*v/INV_ARCSEC/INV_ARCSEC*jit[2]/1e6);

      /* Extended galaxy */
      u1 = u*cos(g.pa)+v*sin(g.pa);
      v1 = v*cos(g.pa)-u*sin(g.pa);
      w = (g.ra*g.ra*u1*u1+g.rb*g.rb*v1*v1)/INV_ARCSEC/INV_ARCSEC;
      factor *= pow(1. + 4.*M_PI*M_PI/1.67834/1.67834*w, -1.5);

      /* IPC */
      factor *= 1.-4.*alpha + 2.*(alpha+alphaplus)*cos(2*M_PI*u*pix/INV_ARCSEC) + 2.*(alpha-alphaplus)*cos(2*M_PI*v*pix/INV_ARCSEC);

      A[l].mtf[i][2*j  ] *= factor;
      A[l].mtf[i][2*j+1] *= factor;
    }

    /* PSF Output. First get the pupil autocorrelation, zero-padded */
    l = lref;
    psfdim = A[l].Nu * 2;
    PSF = dmatrix(0,psfdim-1,0,2*psfdim-1);
    for(ii=0;ii<A[l].Nu;ii++) for(jj=0;jj<A[l].Nu;jj++) {
      i = ii-A[l].Nu/2; if (i<0) i += psfdim;
      j = jj-A[l].Nu/2; if (j<0) j += psfdim;
      PSF[i][2*j  ] = A[l].mtf[ii][2*jj  ];
      PSF[i][2*j+1] = A[l].mtf[ii][2*jj+1];
    }
    deallocate_aperture(A+l);

    /* Flip signs to put PSF in the middle of the frame */
    for(i=0;i<psfdim;i++) for(j=0;j<psfdim;j++) if ((i+j)%2) {
        PSF[i][2*j  ] *= -1;
        PSF[i][2*j+1] *= -1;
    }

    /* Do a 2D FFT */
    for(k=0;k<2;k++) {
      for(i=0; i<psfdim; i++) fourier_trans_2(PSF[i]-1, psfdim, -1);
      /* Transposition */
      for(i=1;i<psfdim;i++) for(j=0;j<i;j++) {
        pt  = PSF[i] + 2*j;
        pt2 = PSF[j] + 2*i;
        temp=*pt; *pt=*pt2; *pt2=temp;
        pt++; pt2++;
        temp=*pt; *pt=*pt2; *pt2=temp;
      }
    }

    for(i=psfdim/2-pad;i<=psfdim/2+pad;i++) {
      for(j=psfdim/2-pad;j<=psfdim/2+pad;j++) {
        Objtot[i-psfdim/2+pad][j-psfdim/2+pad] += weight*PSF[i][2*j];
      }
    }

    free_dmatrix(PSF,0,psfdim-1,0,2*psfdim-1);
  } /* end loop over wavelengths */

  /* Get the moments! */
  mom[0] = 0; /* x */
  mom[1] = 0; /* y */
  mom[2] = 4.*N_SUBPIX*N_SUBPIX; /* Cxx */
  mom[3] = 0; /* Cxy */
  mom[4] = 4.*N_SUBPIX*N_SUBPIX; /* Cyy */
  get_adaptive_gaussmomX(Objtot, pad, mom, mom+1, mom+2, mom+3, mom+4, 1e-9, 180);

#ifndef SPOT_OUTPUT
  printf("%14.7le %14.7le %14.7le %14.7le %14.7le     %10.8lf %11.8lf %11.8lf\n",
      mom[0], mom[1], mom[2], mom[3], mom[4], (mom[2]+mom[4])/N_SUBPIX/N_SUBPIX, (mom[2]-mom[4])/(mom[2]+mom[4]),
      2*mom[3]/(mom[2]+mom[4]));
#endif

#ifdef SPOT_OUTPUT
  printf("# %14.7le %14.7le %14.7le %14.7le %14.7le     %10.8lf %11.8lf %11.8lf\n",
      mom[0], mom[1], mom[2], mom[3], mom[4], (mom[2]+mom[4])/N_SUBPIX/N_SUBPIX, (mom[2]-mom[4])/(mom[2]+mom[4]),
      2*mom[3]/(mom[2]+mom[4]));

  /* Output of the spot */
  for(i=0;i<=2*pad;i++) {
    for(j=0;j<=2*pad;j++)
      printf("%4ld %4ld %13.6le\n", i-pad, j-pad, Objtot[i][j]/psfdim/psfdim*N_SUBPIX*N_SUBPIX);
    printf("\n");
  }
#endif
  free_dmatrix(Objtot,0,2*pad,0,2*pad);

  return(0);
}
