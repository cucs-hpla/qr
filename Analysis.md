To improve performance, we apply some OpenMP multithreading to the reflection function, because that function does quite a bit of work. In particular, we annotate the main loop in the **Reflect** function as follows:

```
void Reflect1(int m, int n, double *a, int lda, const double *x, double tau) {
    #pragma omp parallel for
    for (int i=0; i<n; i++) { // One column at a time
        ...
    }
```

We specify the number of OMP threads at runtime, using **OMP_NUM_THREADS**, and test the QR factorization without OpenMP against OpenMP with 2 and 4 threads. The following plot shows the results (in GFLOPS/s) for varying matrix sizes, with the **dgeqrf** results for each case included for reference.

![Figure 1](https://github.com/seblaud/qr/blob/master/multithreading.png "OpenMP vs. no OpenMP")

We can clearly see an improvement by using OpenMP, and this is particularly apparent when using 2 threads. In addition, we notice that performance for 4 threads is fairly bad for small matrix sizes, which makes sense because of the overhead of multithreading.

Clearly, there is still a great discrepancy between the performance of our implementation and the **dgeqrf** implementation. One thing that would improve performance further is to use blocking in order to orthogonalize several columns at a time. While this was not actually implemented, we expect the blocking variant to improve performance, especially when coupled with OpenMP multithreading.
