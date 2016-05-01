#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/* Inversion of symmetric matrix, C-->Cinv, of dimension N.
 * Uses Cholesky algorithm. Here C is in the 'pointer to pointer' format.
 * Can be used in place by setting Cinv==C (does not cause a memory problem).
 *
 * Returns 0 if successful, 1 if failed (non positive definite).
 */
int sym_matrix_invert(double **C, double **Cinv, int N) {
#define S(a,b) (a+b*(long)N)
  int i,j,k;
  double arg;
  double *L, *M;

  L = (double*)malloc((size_t)(N*(long)N*sizeof(double)));
  M = (double*)malloc((size_t)(N*(long)N*sizeof(double)));
  
  /* Construct the L-matrix: C = L L^T */
  for(j=N-1; j>=0; j--) {
    /* Diagonal element */
    arg = C[j][j];
    for(k=j+1; k<N; k++)
      arg -= L[S(j,k)]*L[S(j,k)];
    if (arg<=0)
      return(1);
    L[S(j,j)] = sqrt(arg);

    /* Off-diagonal elements */
    for(i=j-1; i>=0; i--) {
      arg = C[i][j];
      for(k=j+1; k<N; k++)
        arg -= L[S(i,k)]*L[S(j,k)];
      L[S(i,j)] = arg/L[S(j,j)];
    }
  }

  /* Now the M-matrix */
  for(i=0; i<N; i++) {
    /* Diagonal element */
    M[S(i,i)] = 1./L[S(i,i)];

    /* Off-diagonal elements */
    for(j=i+1; j<N; j++) {
      arg = 0.;
      for(k=i; k<j; k++)   
        arg += M[S(i,k)]*L[S(k,j)];
      M[S(i,j)] = -arg/L[S(j,j)];
    }
  }

  /* Now the C-invese */
  for(i=0; i<N; i++)
    for(j=0; j<N; j++) {
      arg = 0.;
      for(k=0; k<=i && k<=j; k++)
        arg += M[S(k,i)]*M[S(k,j)];
      Cinv[i][j] = arg;
    }
  
  free((char*)L);
  free((char*)M);
#undef S
  return(0);
}

int main(void) {
  int N, i, j;
  double **C, **Cinv;

  /* Read in data, set up arrays */
  scanf("%d", &N);
  C = (double**)malloc((size_t)(N*sizeof(double*)));
  Cinv = (double**)malloc((size_t)(N*sizeof(double*)));
  for(i=0;i<N;i++) {
    C[i] = (double*)malloc((size_t)(N*sizeof(double)));
    Cinv[i] = (double*)malloc((size_t)(N*sizeof(double)));
  }
  for(i=0;i<N;i++) for(j=0;j<N;j++) scanf("%lg", C[i]+j);

  if (sym_matrix_invert(C,Cinv,N)) {
    fprintf(stderr, "Matrix not positive definite\n");
    exit(1);
  }

  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) printf(" %19.12lE", Cinv[i][j]);
    printf("\n");
  }

  /* Cleanup */
  for(i=0;i<N;i++) {
    free((char*)C[i]);
    free((char*)Cinv[i]);
  }
  free((char*)C);
  free((char*)Cinv);

  return(0);  
}
