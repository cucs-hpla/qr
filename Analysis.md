###################################################################
Analysis.md<br />
Eric Peters<br />
cucs-hpla<br />
HW2:Due 2017-03-01<br />
###################################################################

#1.) The starter code does not actually use OpenMP.  Can you annotate loops to improve performance? For which sizes?
First I had to download a new version of gcc via homebrew so that I could use OpenMP.  Then I had to update the makefile so that I could properly use OpenMP.  Below you can see the new makefile.
```c
#matrix size input parameters m & n
#print qr: a=3, no print qr a =0
m = 1000000
n = 100
a = 0

OMP_NUM_THREADS=4

CC = gcc
LINKER = $(CC)
CFLAGS = -Wall -std=c99 -O3 -fopenmp
LDFLAGS = -L/usr/local/opt/openblas/lib -lblas -lm

UTIL := qr-omp.o

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

all : 
	make clean;
	make qr.x

qr.x: $(UTIL) 
	$(LINKER) $(UTIL) $(LDFLAGS) \
	$(BLAS_LIB) $(CFLAGS) -o $(TEST_BIN) $@ 

run:
	make all
	OMP_NUM_THREADS=1 ./qr.x $(m) $(n) $(a)

clean:
	rm -f *.o *~ core *.x

```
Additionally, to improve performance in the code I was able to add some omp parallel loops as well as reductions in order split loops up between multiple threads.  This can be seen in the examples below.
```c
#include <omp.h>

#pragma omp parallel for  
  for (int i=0; i<n; i++) {
    y[i] = (a+b)/2 + (b-a)/2 * cos(M_PI*i / (n-1));
  }  

#pragma omp parallel for
  for (int i=0; i<m; i++) {
    int j = 0;
    for ( ; j<1; j++) b[i+m*j] = 1;
    for ( ; j<2; j++) b[i+m*j] = x[i];
    for ( ; j<n; j++) {
      b[i+m*j] = 2*x[i]*b[i+m*(j-1)] - b[i+m*(j-2)];
    }
  }

#pragma omp parallel for
  for (int i=0; i<m*n; i++) {
    // Logic intended to treat NaN as big
    if (!(fabs(a[i] - b[i]) < max)) max = fabs(a[i] - b[i]);
  }

#pragma omp parallel for reduction(+: sum)
  for (int i=0; i<n; i++) {
    sum += x[i]*y[i];
  }

#pragma omp parallel for
  for (int i=0; i<n; i++) x[i] *= scale;
}

#pragma omp parallel for
  for (int i=0; i<n; i++) {
    y[i] += alpha * x[i];
  }

#pragma omp parallel for
  for (int i=0; i<n; i++) { // One column at a time
    double *ai = &a[0+lda*i];
    double dot = ai[0] + VecDot(m-1, x, ai+1);
    ai[0] -= tau*dot;
    VecAXPY(m-1, ai+1, -tau*dot, x);
  }

```
This does enable some performance enhancment pending on the size of the matrix.  If the size is too small communication becomes relatively more expensive and actually makes the program run slower.  However, for very large matricies great speedups can be seen.  Especially for tall skinny matricies!  I believe this comes from the speed up seen in the dot products when they are longer.  We can even see for certain sizes "TSQR" that the QRFactor can outperform the lapack implementation.  But in general lapack is faster.


#2.) What impediments remain for higher performance?
We did see some speed ups when implementing OpenMP, but there are still some gains to be made.  Much like we saw in the matrix matrix multiply we can take advantage of blocking so that we minimize cache instructions and misses.  It seems that being diligent with memory will always allow for performance enhancments.


#3.) A block QR algorithm (known as Tall-Skinny QR, TSQR) is based on the recurrence<br \>

A0 = [Q0 0] [R0] = [Q0 0] [U0] R<br \>
A1 = [0 Q1] [R1] = [0 Q1] [U1] <br \>

I.e., split a tall skinny matrix A into two parts.  Factor each of these parts (can do in parallel).  The result is two right triangular parts, [R0;R1] which can be factored as [U0; U1] R.   This final right triangular matrix is the R corresponding to a Householder QR.  The product of the orthogonal matrices is the Q.  Can you implement a blocked QR algorithm that orthogonalizes a few columns at a time using this method, then applies the block reflector to the remaining panel of the matrix (analagous to Reflect1 in the code)?  What observations can you make about this algorithm?

The first thing I did was implement a naive compact wy qr factorization.  This takes advantage of the block reflectors instead of just the scalar reflectors.  This implementation can be seen qr-omp-wy.c.  This implementation can be very useful if you choose the correct block size for the reflectors, which is dependent upon your original matrix size.  Next I implemented the tall skinny qr factorization in qr-omp-ts.c, which takes advantage of recursion and the block reflectors.  This algorithm is very interesting in that it breaks the problem down into smaller problems the size of the block reflectors, which works very well for parallel because you don't have to use sequential reflectors.  One can also see in the code that the last R given in the recursion formula is actually the full R for the decomposition.  One can then reconstruct Q by combining all the Q's that were created and stored in the recursion.  I ended up multipling r inverse by A on the right to recover Q in a new way, which we can see is equivelant to what we would get if we reconstructed it using the block reflectors.

All the code can be found in this directory and run using make run.  If you would like to change the code to run all you have to do is change the main file from qr-omp.c, qr-omp-wy.c, qr-omp-ts.c in the Makefile.

The thing to take away from this exercise, is that a naive implementation of qr factorization does not lend itself to parallel computing with the sequential reflectors, but if block reflectors are used it can become very parallelizable with openmp or even distributed memory architecures.
