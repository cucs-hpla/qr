Jorge Luis Barrera Cruz
03/01/2017

---------
HPLA-HW02
---------

1 & 2) Performance improvement and impediments (QR).

In order to improve performance one can add an openmp directive in the "reflect1" function at the beginning of the for loop as follows: '# pragma omp parallel for'. Using this simple modification we can get a significant improvement in performance. 

./qr-omp 900 850:

Initial results:
----------------
  QRFactor   1.479912 s	  0.602 GF/s
    dgeqrf   0.489956 s	  1.819 GF/s
max difference 3.197442e-14, tau 1.554312e-15


Results with openmp directive:
------------------------------
  QRFactor   0.038557 s	  0.380 GF/s
    dgeqrf   0.012633 s	  1.161 GF/s
max difference 8.881784e-15, tau 6.661338e-16

we can see a reduction of approximately 50 times the time in computing the results. However, there are some circumstances under which such improvement is not seen. In small matrices, for example, multithreating does not improve the code as in big matrices. Matrices with relative small number of columns (compared to number of rows) would not see a significant improvement in performance because there are not many reflections to be done. 

If we select some columns, q is represented in v, a part of tau corresponding to the selected region is filled in, and we don't have to apply the reflectors one at a time but we can build up the set of reflectors and apply these blocks to the rest of the panel.
Similar to matrix-matrix multiplications, we have a considerable amount of cache misses if we don't group the data in blocks. A second step for optimizing the code is blocking the number of operations, and applying several reflectors at the same time. We can do, for example, 4 columns at the same time (block reflector) to have better data reuse.

 
 
3) Tall-Skinny QR factorization (TSQR)

(Codes attached)

Observations:
-------------

If the matrices we are dealing with are big and square, all of the work is concentrated on applying the reflector on the rest of the matrix, not to the panel factorization. So we could do a sequential panel factorization and the paralelizing the reflection. But if the matrix doesn't have many columns (tall-skinny matrix) most of the cost goes into the panel factorization. If the number of columns is very small (as it is in this case), then practically everything goes into the panel factorization. This could be benefitial for parallelism. If we do a Q factorization and send it back we can redundantly store small pieces of Q and we only have to do local work when we want to apply Q to some vector. Since the array Q only have a few columns storing pieces of Q is relatively cheap.
It is worth mentioning that we are not restricted to 2 subdivisions at a time. Form a general perspective, the block size is such that the number of rows is defined as a multiple of the number of columns, and changing the block size can affect the performance. Also, communication costs of this approach is less than Householder-QR.
