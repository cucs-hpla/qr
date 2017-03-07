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

// Compute the compact WY QR factorization of A.
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

// Compute the TS recursive QR factorization of A
void QRFactor_TS(int m, int n, double *a, int lda, double *T, double *work, int lwork, int *info, double *R) {
if(m<=2*n){
  QRFactor_WY(m, n, a, m, T, work, lwork, info);
  *info = 0;
}
else{
  double *T1,*T2,*T3,*b,*a1,*a2;
  int m1,m2;
  //separate matricies
  m1 = m/2;
  m2 = m-m/2;
  a1 = calloc(m1*n,sizeof(double));
  a2 = calloc(m2*n,sizeof(double));
  T1 = calloc(n*n,sizeof(double));
  T2 = calloc(n*n,sizeof(double));
  T3 = calloc(n*n,sizeof(double));
  b = calloc(2*n*n,sizeof(double));
  for (int i=0;i<m;i++){
    if(i<m1){
      for(int j=0;j<n;j++){
        a1[i+m1*j] = a[i+m*j];
      }
    }
    else{
      for(int j=0;j<n;j++){
        a2[(i-m1)+m2*j] = a[i+m*j];
      }
    }
  }
  QRFactor_TS(m1,n,a1,lda,T1,work,lwork,info,R);
  QRFactor_TS(m2,n,a2,lda,T2,work,lwork,info,R);
  //concatenate
  for (int j=0;j<n;j++){
    for (int i=0;i<j+1;i++){
      b[i+(2*n)*j] = a1[i+m1*j];
    }
  }
  for (int j=0;j<n;j++){
    for (int i=0;i<j+1;i++){
      b[(i+n)+(2*n)*j] = a2[i+m2*j];
    }
  }

  //new qr
  QRFactor_TS(2*n,n,b,2*n,T3,work,lwork,info,R);

  //Final R is extracted from b after qr
  //Save R but not Q
  for (int j=0;j<n;j++){
    for (int i=0;i<j+1;i++){
      R[i+n*j] = b[i+(2*n)*j];
    }
  }
  free(T1);free(T2);free(b);
  *info = 0;
}
}

// Declare LAPACK dgeqrf (assuming the most common Fortran name mangling)
extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dgetri_(int*, double*, int*, int*, double*, int*, int*);
extern void dtrtri_( char*, char*, int*, double*, int*, int*);
extern void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

int main(int argc, char **argv) {
  char uplo = 'U';
  char diag = 'N';
  char transa = 'n';
  char transb = 'n';
  double alpha=1.0,beta=1.0;
  int verbose = 0;
  int m, n, lwork, info;
  int *ipiv;
  double *x, *a, *b, *tau_a, *tau_b, *work, t, *R, *AA, *QQ;
  if (argc < 3 || argc > 4) {
    fprintf(stderr, "usage: %s M N [Verbose]\n", argv[0]);
    return 1;
  }
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  if (argc == 4) verbose = atoi(argv[3]);
  lwork = n*n;
  MallocA(lwork, &work);
  MallocA(n+1,&ipiv);
  Cosspace(-1, 1, m, &x);
  VanderCheb(m, n, x, &a);
  VanderCheb(m, n, x, &AA);
  QQ = calloc(m*n,sizeof(double));
  MallocA(n*n, &tau_a);
  MallocA(n*n,&R);
  t = gettime();
  QRFactor_TS(m, n, a, m, tau_a, work, lwork, &info, R);
  t = gettime() - t;
  printf("%10s %10.6f s\t%7.3f GF/s\n", "QRFactor", t, (2.*m*n*n - 2./3*n*n*n)*1e-9/t);
  if (verbose) {
    FPrintMatrix(stdout, "R", "% 10.6f", n, n, R);
    printf("\n");
  }
//Compute Q from r^-1 and see
//how accurate A is
  dtrtri_(&uplo,&diag,&n,R,&n,&info);
  dgemm_(&transa,&transb,&m,&n,&n,&alpha,AA,&m,R,&n,&beta,QQ,&m);
  if (verbose) {
    FPrintMatrix(stdout, "Q", "% 10.6f", m, n, QQ);
    printf("\n");
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

  printf("max difference %e, tau %e\n", MaxDifference(m, n, a, b), MaxDifference(1, n,  tau_a, tau_b));

  free(x);
  free(a);
  free(b);
  free(tau_a);
  free(tau_b);
  free(work);
  return 0;
}
