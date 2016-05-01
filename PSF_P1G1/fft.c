#define Pi                3.1415926535897932384626433832795
#define TwoPi             6.283185307179586476925286766559

#define INTERP_XMIN (1.0e-12) /* This is used in interp_coefs_1 */

/* fourier_trans_1
 * *** FOURIER TRANSFORMS A DATA SET WITH LENGTH A POWER OF 2 ***
 *
 * This is a Fourier transform routine.  It has the same calling
 * interface as Numerical Recipes four1 but uses recursion.
 *
 * Replaces data[1..2*nn] by its discrete Fourier transform, if
 * isign is input as 1; or replaces data[1..2*nn] by nn times its
 * inverse discrete Fourier transform, if isign is input as -1.
 * data is a complex array of length nn.
 *
 */
      
void fourier_trans_1(double *data, long nn, int isign) {
  double *Dodd, *Deven, c, s, dc, ds, temp[8], *myd;
  long i, j, twonn;

  data++; /* Remove unit offset */

  /* Length 1 is trivial */
  if (nn==1) return;

  /* Length 2 is trivial as well */
  if (nn==2) {
    temp[0] = data[0];
    temp[1] = data[1];
    data[0] = temp[0] + data[2];
    data[1] = temp[1] + data[3];
    data[2] = temp[0] - data[2];
    data[3] = temp[1] - data[3];
    return;
  }

  /* Length 4 */
  if (nn==4) {
    /* Get real space data in "temp" array; reverse order for isign==-1 */
    temp[0] = data[0]; temp[1] = data[1];
    temp[4] = data[4]; temp[5] = data[5];
    if (isign==1) {
      temp[2] = data[2]; temp[3] = data[3];
      temp[6] = data[6]; temp[7] = data[7];
    } else {
      temp[2] = data[6]; temp[3] = data[7];
      temp[6] = data[2]; temp[7] = data[3];
    }
    data[0] = temp[0] +temp[2] +temp[4] +temp[6];
    data[1] = temp[1] +temp[3] +temp[5] +temp[7];
    data[2] = temp[0] -temp[3] -temp[4] +temp[7];
    data[3] = temp[1] +temp[2] -temp[5] -temp[6];
    data[4] = temp[0] -temp[2] +temp[4] -temp[6];
    data[5] = temp[1] -temp[3] +temp[5] -temp[7];
    data[6] = temp[0] +temp[3] -temp[4] -temp[7];
    data[7] = temp[1] -temp[2] -temp[5] +temp[6];
    return;
  }

  /* Recursive method for length 8 or longer */
  Deven = (double*)malloc((size_t)(2*nn*sizeof(double)));
  Dodd = Deven+nn;
  myd = data-1;
  for(j=0;j<nn;j+=2) {
    Deven[j] = *++myd;
    Deven[j+1] = *++myd;
    Dodd[j] = *++myd;
    Dodd[j+1] = *++myd;
  }
  fourier_trans_1(Deven-1,nn>>1,isign);
  fourier_trans_1(Dodd-1,nn>>1,isign);
  dc = cos(TwoPi/nn); ds = sin(TwoPi/nn)*isign;
  c = 1.; s = 0.;
  twonn=2*nn;
  for(j=0;j<twonn;j+=2) {
    i = j%nn;
    data[j  ] = Deven[i  ] + c*Dodd[i  ] - s*Dodd[i+1];
    data[j+1] = Deven[i+1] + c*Dodd[i+1] + s*Dodd[i  ];
    s = c*ds+(temp[0]=s)*dc;
    c = c*dc-temp[0]*ds;
  }
  free((char*)Deven);
}

/* fourier_trans_2
 * *** FOURIER TRANSFORMS A DATA SET OF ARBITRARY LENGTH ***
 *
 * This is a Fourier transform routine.  It has the same calling
 * interface as NR four1 and our fourier_trans_1, but it works for
 * non-power-of-two lengths.
 *
 * This code works by finding the largest factor of 2 dividing nn
 * and doing this portion via FFT.  The odd factors of nn are computed
 * via brute force.  Consequently this routine will run much faster
 * if the length of your data set is a power of 2 times a small odd
 * number.  If your data set has length exactly a power of 2, it is
 * recommended that you use fourier_trans_1 to reduce RAM usage.
 */
      
void fourier_trans_2(double *data, long nn, int isign) {
      
   double *data_ptr;
   double *data2;
   long i,j,p,b,nu;
   double oneminuscosphi, sinphi, temp;
   double cexpr, cexpi;
   double *costable, *sintable;
   unsigned long nn2, unsigned_comp_nn, mantissa, two_nn2;
   unsigned long i_in, i_out;
   long phase;
            
   data2 = dvector(1,2*nn);
            
   /* Separate nn into mantissa*nn2, where nn2 is a power of 2.
    * We want nn2 to be the largest power of 2 for which the
    * mantissa is an integer.
    */
   unsigned_comp_nn = ~ (unsigned long)nn;
   nn2=1;
   while (nn2 & unsigned_comp_nn) nn2 <<= 1;
   two_nn2 = nn2 << 1;
   mantissa = nn/nn2;
          
   /* Now re-arrange the data as follows: element mantissa*p+b of data
    * is moved to element nn2*b+p of data2.  Note that we're moving
    * complex numbers and the array indices are doubled accordingly.
    */
   i_in = 0;
   for(p=0;p<nn2;p++) {
      for(b=0;b<mantissa;b++) {
         /* We already have i_in == (mantissa*p+b) << 1 */
         i_out = (nn2*b+p) << 1;
         data2[++i_out] = data[++i_in];
         data2[++i_out] = data[++i_in];
      } /* end b loop */
   } /* end p loop */

   /* Do the FFT of each component of data2 */
   for(b=0;b<mantissa;b++) fourier_trans_1(data2+b*two_nn2,nn2,isign);
 
   /* Build the trig function tables for the order-mantissa Fourier
    * transform.  costable[i] = cos ( TwoPi/nn * i ).  Similar for 
    * sintable.
    */
   oneminuscosphi = sin(Pi/nn * isign);
   oneminuscosphi *= 2*oneminuscosphi;
   sinphi = sin(TwoPi/nn * isign);
   costable = dvector(0,nn-1);
   sintable = dvector(0,nn-1);
   cexpr = 1.;   
   cexpi = 0.;
   for(i=0;i<nn;i++) {
      costable[i] = cexpr;
      sintable[i] = cexpi;
      cexpr -= oneminuscosphi * (temp=cexpr) + sinphi * cexpi;
      cexpi += sinphi * temp - oneminuscosphi * cexpi;
   }
            
   /* It remains to do the order-mantissa Fourier transform to
    * complete the operation.  We do this by clearing the input
    * array and doing the order-mantissa Fourier transform from
    * data2 to data.
    */
   data_ptr = data;
   for(i=0;i<nn;i++) {
      data_ptr[1] = data_ptr[2] = 0.;
      nu = 2*(i % nn2);
      phase = 0; /* phase will be (i*j)%nn */
      for(j=0;j<mantissa;j++) {
         data_ptr[1] += costable[phase] * data2[nu+1] - sintable[phase] * data2[nu+2];
         data_ptr[2] += sintable[phase] * data2[nu+1] + costable[phase] * data2[nu+2];
         if ((phase += i) >= nn) phase -= nn;
         nu += two_nn2;
      } /* end j loop */
   
      data_ptr += 2;   
   } /* end i loop */
      
   /* Clean up memory */
   free_dvector(data2,1,2*nn);
   free_dvector(costable,0,nn-1);
   free_dvector(sintable,0,nn-1);
}

/* interp_coefs_1
 * *** FINDS THE COEFFICIENTS FOR A POLYNOMIAL INTERPOLATION ***
 *
 * This routine computes the coefficients c_i for determining
 * y(x) given y evaluated at -M, -M+1, ... N
 * using an interpolating polynomial of order M+N.  We require x to
 * be in the range [0,1].  The formula for y(x) is:
 *
 * y(x) = sum_i c[i] y(i), i=-M .. N.
 */
  
void interp_coefs_1(double x, int M, int N, double *c) {
   
   int i,j;
   double product;
   static unsigned short int is_initialized = 0;
   static double **factorial_products;
      
#define MAX_FACTORIAL 21
   /* factorial_array[n] = n! */
   static double factorial_array[]={1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880.,
      3628800.,39916800.,479001600.,6227020800.,1307674368000.,20922789888000.,
      355687428096000.,355687428096000.,6402373705728000.,121645100408832000.,
      2.43290200817664e+18,5.109094217170944e+19};
      
   /* If this routine hasn't been called yet, allocate the factorial products,
    * factorial_products[i][j] = i!j!.  Note that we never de-allocate this
    * memory, it is freed only when the program exits.  This is no big loss
    * since we only allocate it once and even then it's only a few kbytes.
    */
   if (!is_initialized) {
      is_initialized = 1;
      factorial_products = dmatrix(0,MAX_FACTORIAL,0,MAX_FACTORIAL);
      for(i=0;i<=MAX_FACTORIAL;i++)
         for(j=0;j<=MAX_FACTORIAL;j++)
            factorial_products[i][j] = factorial_array[i] * factorial_array[j];
   }
   
#ifdef N_CHECKVAL
   /* Error indicators */
   if (M<0 || N<1 || M+N>21) {
      fprintf(stderr,"Error: Bad order for interp_coefs_1: M=%d, N=%d\n",M,N);
      exit(1);
   }
#endif
 
   /* First get rid of the cases where x is close to 0 or to 1;
    * we need this branching because for the general case we divide
    * by x and (x-1).
    */
   if (fabs(x)<INTERP_XMIN) {
      for(i= -M;i<=N;i++) c[i]=0.;
      c[0] = 1.;
      return;
   }
   if (fabs(x-1.)<INTERP_XMIN) {
      for(i= -M;i<=N;i++) c[i]=0.;
      c[1] = 1.;
      return;
   }
   
   /* Compute product = (x+M)(x+M-1) ... (x+1)x(x-1) ... (x-N) */
   product = 1.;
   for(i= -M;i<=N;i++) product *= x-i;
      
   /* Compute interpolation coefficients */
   for(i=N;i>= -M;i--) {
      if ((N-i)%2) {
         c[i] = -product / ( (x-i) * factorial_products[N-i][M+i] );
      } else {
         c[i] =  product / ( (x-i) * factorial_products[N-i][M+i] );
      }
   }
}

