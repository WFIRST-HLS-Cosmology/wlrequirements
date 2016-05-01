/* Definitions for MTF computer */

typedef struct {
  double dx; /* grid spacing in meters */
  long Nx; /* number of grid points in 1D */
  double lambda; /* wavelength in meters */
  double **aperture; /* array of complex numbers: (i,j) pixel is aperture[i][2*j] (re), 2*j+1 (im) */
  double du; /* grid spacing of MTF in cycles per radian */
  long Nu, ictr; /* number of grid points in Fourier plane, and index of center */
  double **mtf; /* array of complex numbers: (i,j) pixel in uv-plane is mtf[i][2*j] (re), 2*j+1 (im) */
} APERTURE_TYPE;

/* Create and destroy apertures */
void allocate_aperture(APERTURE_TYPE *A, long Nx) {
  long i,j;
  A->Nx = Nx;
  A->Nu = 2*Nx;
  A->ictr = Nx;
  A->aperture = dmatrix(0,Nx-1,0,2*Nx-1);  
  A->mtf = dmatrix(0,A->Nu-1,0,2*A->Nu-1);  
  for(i=0;i<Nx;i++) for(j=0;j<Nx;j++) A->aperture[i][2*j] = A->aperture[i][2*j+1] = 0.;
}

void deallocate_aperture(APERTURE_TYPE *A) {
  free_dmatrix(A->aperture,0,A->Nx-1,0,2*A->Nx-1);
  free_dmatrix(A->mtf,0,A->Nu-1,0,2*A->Nu-1);
}

/* Zernike polynomial in aperture Ra */
double ZernikePoly(int n, int m, double x, double y, double Ra) {

  double trig, phi, r, xj;
  int J_order, nj, absm;
  double RR, Pnext, P, Pprev;

  /* Error checking */
  if (abs(m)>n || (n+m)%2==1) {
    fprintf(stderr, "Warning: illegal Zernike mode (%d,%d)\n", n, m);
    return(0);
  }

  /* Rectangular --> polar */
  r = sqrt(x*x + y*y)/Ra;
  phi = atan2(y,x);
  trig = m>=0? cos(m*phi): sin(-m*phi);

  /* Need the Jacobi polynomials -- alpha = 0, beta = |m| */
  xj = 2.*r*r-1.;
  absm = abs(m);
  J_order = (n-absm)/2;

  RR = 1.;
  if (J_order>=1) {
    Pprev = 0.; P=1.;
    for(nj=0; nj<J_order; nj++) {
      if (absm>0 || nj>0) {
        Pnext = ( -(2*nj+absm+1.)*absm*absm + (2.*nj+absm)*(2*nj+absm+1)*(2*nj+absm+2)*xj ) * P
                - 2.*nj*(nj+absm)*(2*nj+2+absm) * Pprev;
        Pnext /= 2.*(nj+1)*(nj+1+absm)*(2*nj+absm);
      } else {
        Pnext = xj;
      }
      Pprev = P; P = Pnext;
    }
    RR = P;
  }

  for(nj=0;nj<absm;nj++) RR *= r;

  return(RR*trig);
}

/* Sets up an annular mask for the aperture, i.e. all 1's within, 0's outside */
void set_aperture_mask(APERTURE_TYPE *A, double rmin, double rmax) {

  long i,j;
  double delta_x, delta_y, r;
  double xctr,yctr;

  xctr = 0.5*A->dx*(A->Nx-1);
  yctr = 0.5*A->dx*(A->Nx-1);

  for(i=0;i<A->Nx;i++) for(j=0;j<A->Nx;j++) {
    delta_x = i*A->dx-xctr;
    delta_y = j*A->dx-yctr;
    r = sqrt(delta_x*delta_x+delta_y*delta_y);
    A->aperture[i][2*j] = r>=rmin && r<rmax? 1: 0;
    A->aperture[i][2*j+1] = 0.;
  }
}

/* Applies a Zernike mode path length shift, with amplitude in meters, with radius rmax */
void shift_phase_Z(APERTURE_TYPE *A, int n, int m, double shift, double rmax) {

  long i,j;
  double xctr,yctr,delta_x,delta_y;
  double phase, phaseamp, old_re, old_im;

  xctr = 0.5*A->dx*(A->Nx-1);
  yctr = 0.5*A->dx*(A->Nx-1);
  phaseamp = 2.*M_PI/A->lambda*shift;

  for(i=0;i<A->Nx;i++) for(j=0;j<A->Nx;j++) {
    delta_x = i*A->dx-xctr;
    delta_y = j*A->dx-yctr;
    phase = phaseamp * ZernikePoly(n, m, delta_x, delta_y, rmax);
    old_re = A->aperture[i][2*j  ];
    old_im = A->aperture[i][2*j+1];
    A->aperture[i][2*j  ] = old_re*cos(phase) -old_im*sin(phase);
    A->aperture[i][2*j+1] = old_im*cos(phase) +old_re*sin(phase);
  }
}

/* Generates the MTF for the aperture.  This is an autocorrelation,
 * so we take the Fourier transform, square, and transform back.
 */
void generate_mtf(APERTURE_TYPE *A) {

  double c;
  long i,j,ic,jc;
  double **temp;

  /* Set spacing */
  A->du = A->dx/A->lambda;

  /* Build temporary grid on which to do FFTs, and populate.
   * The grid must be 2x size to avoid artifacts due to overlap across the boundaries. */
  temp = dmatrix(0,A->Nu-1,0,2*A->Nu-1);
  for(i=0;i<A->Nu;i++) for(j=0;j<A->Nu;j++) temp[i][2*j] = temp[i][2*j+1] = 0.;
  for(i=0;i<A->Nx;i++) for(j=0;j<A->Nx;j++) {
    temp[i][2*j  ] = A->aperture[i][2*j  ];
    temp[i][2*j+1] = A->aperture[i][2*j+1];
  }

  /* To accomplish 2D FFT: FT the y-direction, transpose, FT the y-direction, transpose. */
  for(i=0;i<A->Nu;i++) fourier_trans_2(temp[i]-1, A->Nu, 1);
  for(i=0;i<A->Nu;i++) for(j=0;j<i;j++) {
    c = temp[i][2*j  ]; temp[i][2*j  ] = temp[j][2*i  ]; temp[j][2*i  ] = c;
    c = temp[i][2*j+1]; temp[i][2*j+1] = temp[j][2*i+1]; temp[j][2*i+1] = c;
  }
  for(i=0;i<A->Nu;i++) fourier_trans_2(temp[i]-1, A->Nu, 1);
  for(i=0;i<A->Nu;i++) for(j=0;j<i;j++) {
    c = temp[i][2*j  ]; temp[i][2*j  ] = temp[j][2*i  ]; temp[j][2*i  ] = c;
    c = temp[i][2*j+1]; temp[i][2*j+1] = temp[j][2*i+1]; temp[j][2*i+1] = c;
  }

  /* We now have the amplitude at position theta, take square norm */
  for(i=0;i<A->Nu;i++) for(j=0;j<A->Nu;j++) {
    c = temp[i][2*j  ]*temp[i][2*j  ] + temp[i][2*j+1]*temp[i][2*j+1];
    temp[i][2*j  ] = c;
    temp[i][2*j+1] = 0;
  }

  /* To accomplish 2D FFT: FT the y-direction, transpose, FT the y-direction, transpose. */
  for(i=0;i<A->Nu;i++) fourier_trans_2(temp[i]-1, A->Nu, -1);
  for(i=0;i<A->Nu;i++) for(j=0;j<i;j++) {
    c = temp[i][2*j  ]; temp[i][2*j  ] = temp[j][2*i  ]; temp[j][2*i  ] = c;
    c = temp[i][2*j+1]; temp[i][2*j+1] = temp[j][2*i+1]; temp[j][2*i+1] = c;
  }
  for(i=0;i<A->Nu;i++) fourier_trans_2(temp[i]-1, A->Nu, -1);
  for(i=0;i<A->Nu;i++) for(j=0;j<i;j++) {
    c = temp[i][2*j  ]; temp[i][2*j  ] = temp[j][2*i  ]; temp[j][2*i  ] = c;
    c = temp[i][2*j+1]; temp[i][2*j+1] = temp[j][2*i+1]; temp[j][2*i+1] = c;
  }

  /* Output, normalized by (0,0) */
  for(i=0;i<A->Nu;i++) for(j=0;j<A->Nu;j++) {
    ic = (i + A->Nu/2)%A->Nu;
    jc = (j + A->Nu/2)%A->Nu;
    A->mtf[i][2*j  ] = temp[ic][2*jc  ]/temp[0][0];
    A->mtf[i][2*j+1] = temp[ic][2*jc+1]/temp[0][0];
  }

  free_dmatrix(temp,0,A->Nu-1,0,2*A->Nu-1);
}

/* Extract MTF at a particular point in the (u,v) plane by interpolation.
 * The flag is either 0 (real part) or 1 (imag part).
 */
double mtf_extract(APERTURE_TYPE *A, double u, double v, int flag) {

  long i_int, j_int;
  double ie,je, i_frac, j_frac;

  /* Effective position to extract MTF */
  ie = A->Nu/2 + u/A->du;
  je = A->Nu/2 + v/A->du;

  /* Interpolation parameters */
  i_int = (long)floor(ie);  i_frac = ie-i_int;
  j_int = (long)floor(je);  j_frac = je-j_int;

  /* Off grid */
  if (i_int<0 || i_int>=A->Nu-1 || j_int<0 || j_int>=A->Nu-1) return(0.);

  return(
      (1.-i_frac)*(1.-j_frac)*A->mtf[i_int  ][2*j_int  +flag]
    + (1.-i_frac)*(   j_frac)*A->mtf[i_int  ][2*j_int+2+flag]
    + (   i_frac)*(1.-j_frac)*A->mtf[i_int+1][2*j_int  +flag]
    + (   i_frac)*(   j_frac)*A->mtf[i_int+1][2*j_int+2+flag]
  );
}

/* Gets adaptive centroid and covariance matrix for a postage stamp.
 * region is object[0..2*offset][0..2*offset]
 * 
 * Computes Gaussian weighted moments of the postage stamp:
 * out[0] = int G(x,y) exp(-r C^-1 r/2) d^2r
 * out[1] has factor of x-xctr
 * out[2] has factor of y-yctr
 * out[3] has factor of (x-xctr)^2
 * out[4] has factor of (x-xctr)(y-yctr)
 * out[5] has factor of (y-yctr)^2
 *
 * All integrals are implemented in real space, in units of the pixel in object.
 */
void get_gaussmomX(double **object, long offset, double xctr, double yctr, double Cxx, double Cxy, double Cyy,
  double *out) {

  long i, j;
  double G, Cixx, Cixy, Ciyy;
  double x,y;

  Cixx =  Cyy/(Cxx*Cyy-Cxy*Cxy);
  Cixy = -Cxy/(Cxx*Cyy-Cxy*Cxy);
  Ciyy =  Cxx/(Cxx*Cyy-Cxy*Cxy);

  out[0] = out[1] = out[2] = out[3] = out[4] = out[5] = 0.;

  for(i=-offset; i<=offset; i++) for(j=-offset; j<=offset; j++) {

    x = i-xctr; y = j-yctr;
    G = exp(-0.5*Cixx*x*x - 0.5*Ciyy*y*y - Cixy*x*y) * object[i+offset][j+offset];

    out[0] += G;
    out[1] += G*x;
    out[2] += G*y;
    out[3] += G*x*x;
    out[4] += G*x*y;
    out[5] += G*y*y;
  }

  for(i=0;i<6;i++) out[i] /= 2*M_PI*sqrt(Cxx*Cyy-Cxy*Cxy);
}

/* Gets adaptive centroid and covariance matrix for a postage stamp.
 * region is object[0..2*offset][0..2*offset]
 *
 * The tolerance on changes to the centroid and covariance is given.
 */
void get_adaptive_gaussmomX(double **object, long offset, double *xctr, double *yctr, double *Cxx, double *Cxy, double *Cyy,
  double tol, int itermax) {

  int i=0, j;
  double err=1, con;
  double out[6], sigma, errsig[5];

  while (i<itermax && err>tol) {

    get_gaussmomX(object, offset, *xctr, *yctr, *Cxx, *Cxy, *Cyy, out);
    sigma = sqrt(*Cxx>*Cyy?*Cxx:*Cyy);
    errsig[0] = out[1]/out[0]/sigma;
    errsig[1] = out[2]/out[0]/sigma;
    errsig[2] = (out[3]/out[0]-*Cxx/2.)/sigma/sigma;
    errsig[3] = (out[4]/out[0]-*Cxy/2.)/sigma/sigma;
    errsig[4] = (out[5]/out[0]-*Cyy/2.)/sigma/sigma;

    /* Compute errors */
    err = 0; for(j=0;j<5;j++) if (fabs(errsig[j])>err) err=fabs(errsig[j]);
    i++;

    /* Updates */
    con = i<8? 0.5: 1.;
    *xctr += con * errsig[0] * sigma;
    *yctr += con * errsig[1] * sigma;
    *Cxx  += con * errsig[2] * sigma * sigma;
    *Cxy  += con * errsig[3] * sigma * sigma;
    *Cyy  += con * errsig[4] * sigma * sigma;

#if 0
    fprintf(stderr, "%3d %12.5le | %12.5le %12.5le %12.5le %12.5le %12.5le\n", i, err, *xctr, *yctr, *Cxx, *Cxy, *Cyy);
#endif
  }

  if (err>tol) {fprintf(stderr, "Convergence failure: err = %12.5le\n", err);}
}
