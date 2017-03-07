#define _XOPEN_SOURCE 700

// from http://stackoverflow.com/questions/3437404/min-and-max-in-c
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cblas.h>
#include <lapacke.h>
#include <string.h>
#include <assert.h>

#define MallocA(n, ptr) posix_memalign((void**)ptr, 64, (n) * sizeof(**(ptr)));

double gettime() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return (double)tp.tv_sec + 1e-6*tp.tv_usec;
}

void Cosspace(double a, double b, int n, double **x) {
  double *y;
  MallocA(n, &y);
  for (int i=0; i<n; i++) {
    y[i] = (a+b)/2 + (b-a)/2 * cos(M_PI*i / (n-1));
  }
  *x = y;
}

// Vandermonde matrix at m points x[0:m], Chebyshev polynomial basis
void VanderCheb(int m, int n, double *x, double **a) {
  double *b;
  MallocA(m*n, &b);
  for (int i=0; i<m; i++) {
    int j = 0;
    for ( ; j<1; j++) b[i+m*j] = 1;
    for ( ; j<2; j++) b[i+m*j] = x[i];
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

void FPrintMatrix2(FILE *stream, const char *name, const char *fmt, int m, int n, const double *a, int lda) {
  for (int i=0; i<m; i++) {
    fprintf(stream, "%s[%d,:]", name, i);
    for (int j=0; j<n; j++) {
      fputc(' ', stream);
      fprintf(stream, fmt, a[i+lda*j]);
    }
    fputc('\n', stream);
  }
}

double MaxDifference(int m, int n, const double *a, const double *b) {
  double max = 0;
  for (int i=0; i<m*n; i++) {
    // Logic intended to treat NaN as big
    if (!(fabs(a[i] - b[i]) < max)) max = fabs(a[i] - b[i]);
  }
  return max;
}

// Compute the dot product of two vectors
double VecDot(int n, const double *x, const double *y) {
  double sum = 0;
  #pragma omp simd reduction(+:sum)
  for (int i=0; i<n; i++) {
    sum += x[i]*y[i];
  }
  return sum;
}

// Scale a vector
void VecScale(int n, double *x, double scale) {
  for (int i=0; i<n; i++) x[i] *= scale;
}

// Compute y += alpha * x
void VecAXPY(int n, double *restrict y, double alpha, const double *restrict x) {
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
void QRFactor(int m, int n, double *a, int lda, double *tau, double *work, int lwork, int *info) {
  for (int i=0; i<n; i++) {
    if (m-i-1 == 0) { // No sub-diagonal
      tau[n-1] = 0;
      break;
    }
    double *v = &a[i+lda*i];
    double d = VecDot(m-i-1, v+1, v+1);
    double norm = sqrt(d + v[0]*v[0]);
    double Rii = -copysign(norm, v[0]);
    v[0] -= Rii;
    norm = sqrt(d + v[0]*v[0]) / v[0]; // New norm of v after modification above and scaling below
    VecScale(m-i, v, 1/v[0]);
    tau[i] = 2 / (norm*norm);
    a[i + lda*i] = Rii;
    Reflect1(m-i, n-i-1, &a[i+lda*(i+1)], lda, &v[1], tau[i]);
  }
  *info = 0;
}

// tau is n x n matrix
void qr_wy(int m, int n, double *a, int lda, double *tau, int ldt) {
  // zero out tau
  for(int i = 0; i < n; i++) {
    if (m-i-1 == 0) { // No sub-diagonal
      tau[n-1] = 0;
      break;
    }
    double *v = &a[i+lda*i];
    double d = VecDot(m-i-1, v+1, v+1);
    double norm = sqrt(d + v[0]*v[0]);
    double Rii = -copysign(norm, v[0]);
    v[0] -= Rii;
    tau[i + ldt * i] = 2 * v[0]*v[0] / (v[0]*v[0] + d);
    VecScale(m-i, v, 1/v[0]);

    // apply to rest of panel
    Reflect1(m-i, n-i-1, &a[i+lda*(i+1)], lda, &v[1], tau[i + ldt * i]);

    /* FPrintMatrix2(stdout, "Tau", "% 10.6f", 4, n, tau, ldt); */
    if (i > 0) {
      // compute y = F[i:,:i].T.dot(v)
      double* y;
      MallocA(i, &y);
      cblas_dgemv(CblasColMajor, CblasTrans, m - i, i, 1, a + i, lda, v, 1, 0, y, 1);
      // compute x = T[:i,:i].dot(F[i:,:i].T.dot(v))
      double* x;
      MallocA(i, &x);
      cblas_dgemv(CblasColMajor, CblasNoTrans, i, i, 1, tau, ldt, y, 1, 0, x, 1);
      VecScale(i, x, -tau[i + ldt * i]);
      for(int j = 0; j < i; j++) {
        tau[j + ldt * i] = x[j];
      }
    }

    a[i + lda * i] = Rii;
  }
}

// Applies Q in place to submatrix
// from http://www.netlib.org/utk/papers/factor/node8.html
void apply_qr_wy(int m, int n, int split, double *a, int lda, double *tau, int ldt) {
  // A = [V A_2]
  // A: m x n
  // V : m x split
  // V = V1
  //     V2
  // V1 : split x split
  // V2 : (m - split) x (n - split)
  // A_2 : m x (n - split)
  double *w;
  MallocA(split*(n-split), &w);
  // A12 : split x (n - split)
  // W = V^T A_2 = V1^T A12 + V2^T A22
  // W: split x (n - split)
  // W = A12
  for(int j = 0; j < n - split; j++) {
    for(int i = 0; i < split; i++) {
      w[i + j * split] = a[i + (j + split) * lda];
    }
  }
  // V1^T A12
  cblas_dtrmm(CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, split, n-split, 1, a, lda, w, split);
  // V2^T A22
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, split, n - split, m - split, 1, a + split, lda, a + lda * split + split, lda, 1, w, split);

  // W = T^T W
  // T: split x split
  cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, split, n-split, 1, tau, ldt, w, split);

  // A_2 = A_2 - V*W = A12 - V1 W
  //                   A22 - V2 W
  // V2*W: (m - split) x (n - split)
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m-split, n-split, split, -1, a + split, lda, w, split, 1, a + lda * split + split, lda);
  // W = V1 * W : split x (n - split)
  cblas_dtrmm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, split, (n - split), 1, a, lda, w, split);
  // A12 = W
  for(int j = 0; j < n - split; j++) {
    for(int i = 0; i < split; i++) {
      a[i + (j + split) * lda] -= w[i + j * split];
    }
  }
  free(w);
}

// t is n x nb
// based off of https://github.com/cucs-hpla/class/LinearAlgebra.ipynb
void blocked_qr(int m, int n, int nb, double *a, int lda, double *tau, int ldt) {
  int b = ceil(n*1.0/nb);
  for(int i = 0; i < b; i++) {
    int off = i * nb;
    int n_blk;
    if(i == b-1) n_blk = n - off;
    else n_blk = nb;

    double *a_sub = a + off + off * lda;
    double *tau_sub = tau + off * ldt;

    // factor left column
    qr_wy(m - off, n_blk, a_sub, lda, tau_sub, ldt);
    if(i < b-1)
    apply_qr_wy(m - off, n - off, n_blk, a_sub, lda, tau_sub, ldt);
  }
}

struct tsqr_data {
  double *tau1;
  double *tau2;
  double *tau3;
  double *r;
};

// based off of https://github.com/cucs-hpla/class/LinearAlgebra.ipynb
struct tsqr_data tsqr(int m, int n, double *a, int lda, int nb) {
  assert(m > 2 * n);

  double * tau1;
  MallocA(nb * n, &tau1);
  blocked_qr(m/2, n, nb, a, lda, tau1, nb);

  double * tau2;
  MallocA(nb * n, &tau2);
  blocked_qr(m-m/2, n, nb, a + m/2, lda, tau2, nb);

  double * R;
  MallocA(n * n * 2, &R);
  int ldr = 2 *n;

  // pack R into a matrix
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < i+1; j++) {
      R[i * ldr + j] = a[i * lda + j];
    }
    for(int j = 0; j < i+1; j++) {
      R[i * ldr + j + n] = a[i * lda + j + m/2];
    }
  }

  // qr on R
  double * tau3;
  MallocA(nb * n, &tau3);
  blocked_qr(2*n, n, nb, R, ldr, tau3, nb);

  struct tsqr_data qr;
  qr.tau1 = tau1;
  qr.tau2 = tau2;
  qr.tau3 = tau3;
  qr.r = R;

  return qr;
}

void manifest_tsqr_q(int m, int n, double *a, int lda, struct tsqr_data qr) {
  // R: (2*n) x n
  // Manifest Q_bar
  double * Q;
  MallocA(2*n*n, &Q);
  memset(Q, 0, 2*n*n*sizeof(double));
  int ldq = 2*n;
  for(int i = 0; i < n; i++) {
    Q[i + ldq * i] = 1;
    for(int j = i; j < 2*n; j++) {
      Q[j + ldq * i] = qr.r[j + ldq * i];
    }
  }

  // multiply upper half
}

// Declare LAPACK dgeqrf (assuming the most common Fortran name mangling)
void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
/* void dgeqrt_(int *m, int *n, int *nb, double *a, int *lda, double *tau, int *ldt, int *lwork,  int *info); */
void dgeqrt2_(int *m, int *n, double *a, int *lda, double *tau, int *ldt, int *info);

int main(int argc, char **argv) {
  int verbose = 0;
  int m, n, nb, lwork, info;
  double *x, *a, *b, *tau_a, *tau_b, *work, t;
  if (argc < 4 || argc > 5) {
    fprintf(stderr, "usage: %s M N NB [Verbose]\n", argv[0]);
    return 1;
  }
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  nb = atoi(argv[3]);
  if (argc == 5) verbose = atoi(argv[4]);
  lwork = 8*n;
  MallocA(lwork, &work);
  Cosspace(-1, 1, m, &x);

  VanderCheb(m, n, x, &a);
  MallocA(n, &tau_a);
  t = gettime();
  QRFactor(m, n, a, m, tau_a, work, lwork, &info);
  t = gettime() - t;
  printf("%10s %10.6f s\t%7.3f GF/s\n", "QRFactor", t, (2.*m*n*n - 2./3*n*n*n)*1e-9/t);
  if (verbose) {
    FPrintMatrix(stdout, "A", "% 10.6f", m, n, a);
    FPrintMatrix(stdout, "tauA", "% 10.6f", 1, n, tau_a);
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

  free(tau_a);
  free(tau_b);

  // qr_wy
  VanderCheb(m, n, x, &a);
  MallocA(n*n, &tau_a);
  t = gettime();
  qr_wy(m, n, a, m, tau_a, n);
  t = gettime() - t;
  printf("%10s %10.6f s\t%7.3f GF/s\n", "qr_wy", t, (2.*m*n*n - 2./3*n*n*n)*1e-9/t);
  if (verbose) {
    FPrintMatrix(stdout, "A", "% 10.6f", m, n, a);
    FPrintMatrix(stdout, "tauA", "% 10.6f", n, n, tau_a);
  }

  VanderCheb(m, n, x, &b);
  MallocA(n*n, &tau_b);
  t = gettime();
  dgeqrt2_(&m, &n, b, &m, tau_b, &n, &info);
  t = gettime() - t;
  printf("%10s %10.6f s\t%7.3f GF/s\n", "dgeqrt2", t, (2.*m*n*n - 2./3*n*n*n)*1e-9/t);
  if (verbose) {
    FPrintMatrix(stdout, "B", "% 10.6f", m, n, b);
    FPrintMatrix(stdout, "tauB", "% 10.6f", n, n, tau_b);
  }

  printf("max difference %e, tau %e\n", MaxDifference(m, n, a, b), MaxDifference(n, n, tau_a, tau_b));

  free(tau_a);
  free(tau_b);

  // blocked qr
  VanderCheb(m, n, x, &a);
  MallocA(nb*n, &tau_a);
  for(int i = 0; i < nb *n; i++) tau_a[i] = 0;
  t = gettime();
  blocked_qr(m, n, nb, a, m, tau_a, nb);
  t = gettime() - t;
  printf("%10s %10.6f s\t%7.3f GF/s\n", "blocked_qr", t, (2.*m*n*n - 2./3*n*n*n)*1e-9/t);
  if (verbose) {
    FPrintMatrix(stdout, "A", "% 10.6f", m, n, a);
    FPrintMatrix(stdout, "tauA", "% 10.6f", nb, n, tau_a);
  }

  VanderCheb(m, n, x, &b);
  MallocA(nb*n, &tau_b);
  for(int i = 0; i < nb *n; i++) tau_b[i] = 0; // make sure lower triangular part is zero
  t = gettime();
  /* dgeqrt_(&m, &n, &nb, b, &m, tau_b, &n, &work, &info); */
  LAPACKE_dgeqrt(LAPACK_COL_MAJOR,m,n,nb,b,m,tau_b,nb);
  t = gettime() - t;
  printf("%10s %10.6f s\t%7.3f GF/s\n", "dgeqrt", t, (2.*m*n*n - 2./3*n*n*n)*1e-9/t);
  if (verbose) {
    FPrintMatrix(stdout, "B", "% 10.6f", m, n, b);
    FPrintMatrix(stdout, "tauB", "% 10.6f", nb, n, tau_b);
  }

  printf("max difference %e, tau %e\n", MaxDifference(m, n, a, b), MaxDifference(nb, n, tau_a, tau_b));

  if(m > 2*n) {
    VanderCheb(m, n, x, &a);
    t = gettime();
    tsqr(m, n, a, m, nb);
    t = gettime() - t;
    printf("%10s %10.6f s\t%7.3f GF/s\n", "tsqr", t, (2.*m*n*n - 2./3*n*n*n)*1e-9/t);
    if (verbose) {
      FPrintMatrix(stdout, "A", "% 10.6f", m, n, a);
    }
  }

  free(x);
  free(a);
  free(b);
  free(tau_a);
  free(tau_b);
  free(work);
  return 0;
}
