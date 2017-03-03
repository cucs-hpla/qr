#define _XOPEN_SOURCE 700

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define MallocA(n, ptr) posix_memalign((void**)ptr, 64, (n) * sizeof(**(ptr)));

double gettime() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return (double)tp.tv_sec + 1e-6*tp.tv_usec;
}

void Cosspace(double a, double b, int n, double **x) {
  double *y;
  MallocA(n, &y);
#pragma omp parallel for  
  for (int i=0; i<n; i++) {
    y[i] = (a+b)/2 + (b-a)/2 * cos(M_PI*i / (n-1));
  }  
  *x = y;
}

// Vandermonde matrix at m points x[0:m], Chebyshev polynomial basis
void VanderCheb(int m, int n, double *x, double **a) {
  double *b;
  MallocA(m*n, &b);
#pragma omp parallel for
  for (int i=0; i<m; i++) {
    int j = 0;
    for ( ; j<1; j++) b[i+m*j] = 1;
    for ( ; j<2; j++) b[i+m*j] = 2*x[i];
    for ( ; j<n; j++) {
      b[i+m*j] = 2*x[i]*b[i+m*(j-1)] - b[i+m*(j-2)];
    }
  }
  *a = b;
}

void FPrintMatrix(FILE *stream, const char *name, const char *fmt, int m, int n, const double *a) {
  for (int i=0; i<m; i++) {
    fprintf(stream, "%s[%d,:]", name, i);
    for (int j=0; j<n; j++) {
      fputc(' ', stream);
      fprintf(stream, fmt, a[i+m*j]);
    }
    fputc('\n', stream);
  }
}

double MaxDifference(int m, int n, const double *a, const double *b) {
  double max = 0;
#pragma omp parallel for
  for (int i=0; i<m*n; i++) {
    // Logic intended to treat NaN as big
    if (!(fabs(a[i] - b[i]) < max)) max = fabs(a[i] - b[i]);
  }
  return max;
}

// Compute the dot product of two vectors
double VecDot(int n, const double *x, const double *y) {
  double sum = 0;
#pragma omp parallel for reduction(+: sum)
  for (int i=0; i<n; i++) {
    sum += x[i]*y[i];
  }
  return sum;
}

// Scale a vector
void VecScale(int n, double *x, double scale) {
#pragma omp parallel for
  for (int i=0; i<n; i++) x[i] *= scale;
}

// Compute y += alpha * x
void VecAXPY(int n, double *y, double alpha, const double *x) {
#pragma omp parallel for
  for (int i=0; i<n; i++) {
    y[i] += alpha * x[i];
  }
}

// Apply the elementary reflection
//   A = (I - tau v v') A
// in-place.  The reflector plane is defined by the vector v = [1; x].
// That is, the first entry of v is implicitly 1 and not stored.
void Reflect1(int m, int n, double *a, int lda, const double *x, double tau) {
#pragma omp parallel for
  for (int i=0; i<n; i++) { // One column at a time
    double *ai = &a[0+lda*i];
    double dot = ai[0] + VecDot(m-1, x, ai+1);
    ai[0] -= tau*dot;
    VecAXPY(m-1, ai+1, -tau*dot, x);
  }
}

// Compute the in-place Householder QR factorization of A.
// This function is meant to be a simple version of dgeqrf.
void dgemv_(char,int,int,double,double**,int,double*,int,double,double*,int);
void QRFactor_WY(int m, int n, double *a, int lda, double *T, double *work, int lwork, int *info) {
  for (int i=0; i<n; i++) {
    double *v = &a[i+lda*i];
    double d = VecDot(m-i-1, v+1, v+1);
    double norm = sqrt(d + v[0]*v[0]);
    double Rii = -copysign(norm, v[0]);
    v[0] -= Rii;
    double tau = 2*pow(v[0],2)/(pow(v[0],2)+d);
    VecScale(m-i, v, 1/v[0]);
    //update remaining panel
    Reflect1(m-i, n-i-1, &a[i+lda*(i+1)], lda, &v[1], tau);
    T[i*n+i]=tau;
    //Add this column to T
    if(i>0){
      double **tmp = malloc((m-i)*sizeof(double*));
      for (int j=0;j<(m-i);j++){ tmp[j] = malloc((i)*sizeof(double));}

      for (int j=0;j<m-i;j++){
        for (int k=0;k<i;k++){
          tmp[j][k] = a[(j+i)+lda*k];
        }
      }
      double *tmp2 = calloc(i,sizeof(double));
      for (int j=0;j<i;j++){
        for (int l=0;l<m-i;l++){
          tmp2[j]+=v[l]*tmp[l][j];
        }
      }
      double *tmp3 = calloc(i,sizeof(double));
      for (int j=0;j<i;j++){
        for (int l=0;l<i;l++){
          tmp3[j]+=tmp2[l]*T[l*n+j];
        }
      }
      for(int j=0;j<i;j++){
        tmp3[j]*=-tau;    
      }
      for(int j=0;j<i;j++){
         T[i*n+j]=tmp3[j];
      } 
      for (int j=0;j<m-i;j++){free(tmp[j]);}
        free(tmp);free(tmp2);free(tmp3);
    }
    a[i+lda*i]=Rii;  
    }
*info = 0;
}

// Declare LAPACK dgeqrf (assuming the most common Fortran name mangling)
extern void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

int main(int argc, char **argv) {
  int verbose = 0;
  int m, n, lwork, info;
  double *x, *a, *b, *tau_a, *tau_b, *work, t;
  if (argc < 3 || argc > 4) {
    fprintf(stderr, "usage: %s M N [Verbose]\n", argv[0]);
    return 1;
  }
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  if (argc == 4) verbose = atoi(argv[3]);
  lwork = 8*n;
  MallocA(lwork, &work);
  Cosspace(-1, 1, m, &x);
  VanderCheb(m, n, x, &a);

  MallocA(n*n, &tau_a);
  t = gettime();
  QRFactor_WY(m, n, a, m, tau_a, work, lwork, &info);
  t = gettime() - t;
  printf("%10s %10.6f s\t%7.3f GF/s\n", "QRFactor", t, (2.*m*n*n - 2./3*n*n*n)*1e-9/t);
  if (verbose) {
    FPrintMatrix(stdout, "A", "% 10.6f", m, n, a);
    FPrintMatrix(stdout, "tauA", "% 10.6f", n, n, tau_a);
  }

  VanderCheb(m, n, x, &b);
  MallocA(n, &tau_b);
  t = gettime();
  dgeqrf_(&m, &n, b, &m, tau_b, work, &lwork, &info);
  t = gettime() - t;
  printf("%10s %10.6f s\t%7.3f GF/s\n", "dgeqrf", t, (2.*m*n*n - 2./3*n*n*n)*1e-9/t);
  if (verbose) {
    FPrintMatrix(stdout, "B", "% 10.6f", m, n, b);
    FPrintMatrix(stdout, "tauB", "% 10.6f", 1, n, tau_b);
  }

  printf("max difference %e, tau %e\n", MaxDifference(m, n, a, b), MaxDifference(1, n, tau_a, tau_b));

  free(x);
  free(a);
  free(b);
  free(tau_a);
  free(tau_b);
  free(work);
  return 0;
}
