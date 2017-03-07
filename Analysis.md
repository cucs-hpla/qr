1) Found below are the results of implementing omp for loops for different sections of of the QR decomposition code. For this analysis I decided to add omp for loops to both vecscale and the reflector functions. Optimizing the for loops at these two points seemed to create the most substantial gains in performance over the serial version. In my git commit, the file with these optimizations is qr-omp_2.c. I attempted to implement omp for loops in two other sections of the code on top of the what is contained in qr-omp_2.c, specifically VecDot and VecAXPY which resulted in a decrease in performance relative to the qr-omp_2.c. The reason I believe this decrease was seen is dude to the fact that too many processes are being forked at once in turn decreasing the performance of the implementation. Since VecDot and VecAXPY are both called inside of a process already being split by omp, this leads to double splitting and too many processes being run for my limited machine.

The sizes that saw the most performance increase were non-square matrices where m or n were much larger than the other. Both when m was much larger than n and n was much larger than m, we saw roughly a factor of 3 increase in both time to complete and GF/s. While there was only a slight increase in performance for  square matrices, the increase was negligible relative to the performance increase seen in rectangular matrices. The output of all runs can be found below the answers to these questions.

2) Multiple impediments lie in the way of greater performance increase when it comes to this implementation of parallel QR decomposition.
	As noted above, since I am working on a personal laptop, we can only fork so many processes before we actually begin to see a decrease in performance. I believe that the QR decomposition could drun faster given more processors in turn enabling us to fork at different locations in the code such as VecDot and VecAXPY. Given more processors, I believe we could then fork inside of already forked threads and see a further increase in performance. Along with that, given more processors we could also push our number of threads larger than the 4 threads we are currently implementing now and see greater performance increases, this holds especially true as our matrices grow larger and larger in size. 
	On top of throwing more processors at the problem, I believe that the memory associated with the QR decompositions could be handled better when creating separate threads. By associating specific cache lines with specific threads we could theoretically decrease the time required to grab data from memory. Although I do not know specifically how PAD works, I have scene parallel omp for loop performance increase when we associate specific chunks of cache memory to threads. 


3) uhhh.... implementing TSQR in C seems to be a bit above my skill level. I haven't personally done all that much C coding, I generally work with higher languages. Thus i apparently can not implement a blocked QR algorithm that decomposes A. I began to outline the structure of the code as given in the notebook on github.  



Serial performance
 ./qr-omp 100 100
  QRFactor   0.002946 s	  0.453 GF/s
    dgeqrf   0.036249 s	  0.037 GF/s
max difference 4.551914e-15, tau 8.881784e-16

./qr-omp 1000 1000
  QRFactor   2.492595 s	  0.535 GF/s
    dgeqrf   0.258535 s	  5.157 GF/s
max difference 3.197442e-14, tau 3.552714e-15

./qr-omp 2000 2000
  QRFactor  19.962204 s	  0.534 GF/s
    dgeqrf   2.511081 s	  4.248 GF/s
max difference 4.840572e-14, tau 4.884981e-15

After omp for on reflect -2 threads

./qr-omp_1 100 100
  QRFactor   0.101836 s	  0.013 GF/s
    dgeqrf   0.013953 s	  0.096 GF/s
max difference 4.551914e-15, tau 8.881784e-16


./qr-omp_1 1000 1000
  QRFactor   0.924971 s	  1.441 GF/s
    dgeqrf   0.281715 s	  4.733 GF/s
max difference 3.197442e-14, tau 3.552714e-15

./qr-omp_1 2000 2000
  QRFactor   6.915847 s	  1.542 GF/s
    dgeqrf   2.239114 s	  4.764 GF/s
max difference 4.840572e-14, tau 4.884981e-15

After omp for reflect -4 threads
./qr-omp_1 100 100
  QRFactor   0.103434 s	  0.013 GF/s
    dgeqrf   0.007348 s	  0.181 GF/s
max difference 4.551914e-15, tau 8.881784e-1

./qr-omp_1 1000 1000
  QRFactor   0.936439 s	  1.424 GF/s
    dgeqrf   0.277876 s	  4.798 GF/s
max difference 3.197442e-14, tau 3.552714e-15

./qr-omp_1 2000 2000
  QRFactor   6.918472 s	  1.542 GF/s
    dgeqrf   2.252109 s	  4.736 GF/s
max difference 4.840572e-14, tau 4.884981e-15


After omp for reflect and VecAXPY -2 threads
./qr-omp_2 100 100
  QRFactor   0.102003 s	  0.013 GF/s
    dgeqrf   0.006385 s	  0.209 GF/s
max difference 4.551914e-15, tau 8.881784e-16

./qr-omp_2 1000 1000
  QRFactor   1.114044 s	  1.197 GF/s
    dgeqrf   0.274278 s	  4.861 GF/s
max difference 3.197442e-14, tau 3.552714e-15

./qr-omp_2 2000 2000
  QRFactor   8.088032 s	  1.319 GF/s
    dgeqrf   2.262155 s	  4.715 GF/s
max difference 4.840572e-14, tau 4.884981e-15


After omp for reflect and VecAXPY -4 thread

./qr-omp_2 100 100
  QRFactor   0.100064 s	  0.013 GF/s
    dgeqrf   0.007563 s	  0.176 GF/s
max difference 4.551914e-15, tau 8.881784e-16


./qr-omp_2 1000 1000
  QRFactor   1.107742 s	  1.204 GF/s
    dgeqrf   0.282392 s	  4.722 GF/s
max difference 3.197442e-14, tau 3.552714e-15

./qr-omp_2 2000 2000
  QRFactor   8.090398 s	  1.318 GF/s
    dgeqrf   2.236178 s	  4.770 GF/s
max difference 4.840572e-14, tau 4.884981e-15

As we can see, adding omp for to vec axpy actually decreased performance, I believe that this decreased performance is due to the fact that we are forking too many processes when doing this. Since VecAXPY is nested inside our Reflect for loop which has already been set to run on multiple threads, when we enter VecAXPY we fork threads inside of already forked threads. Since I am running on my laptop, and only have access to 2 processors, having this extra forking is probably a hindrance. If I were to be on a system with more processors, this may not be the case and it is possible we would actually see an increase in performance.  




After omp for reflect and VecScale -2 thread

./qr-omp_2 100 100
  QRFactor   0.106065 s	  0.013 GF/s
    dgeqrf   0.006435 s	  0.207 GF/s
max difference 4.551914e-15, tau 8.881784e-16

./qr-omp_2 1000 1000
  QRFactor   0.934085 s	  1.427 GF/s
    dgeqrf   0.279365 s	  4.773 GF/s
max difference 3.197442e-14, tau 3.552714e-15

./qr-omp_2 2000 2000
  QRFactor   6.863758 s	  1.554 GF/s
    dgeqrf   2.230581 s	  4.782 GF/s
max difference 4.840572e-14, tau 4.884981e-15

After omp for reflect and VecScale -4 thread

./qr-omp_2 100 100
  QRFactor   0.099697 s	  0.013 GF/s
    dgeqrf   0.002259 s	  0.590 GF/s
max difference 4.551914e-15, tau 8.881784e-16

./qr-omp_2 1000 1000
  QRFactor   0.935355 s	  1.425 GF/s
    dgeqrf   0.275324 s	  4.843 GF/s
max difference 3.197442e-14, tau 3.552714e-15

./qr-omp_2 2000 2000
  QRFactor   6.909139 s	  1.544 GF/s
    dgeqrf   2.232362 s	  4.778 GF/s
max difference 4.840572e-14, tau 4.884981e-15

As we can see, this second optimization did not significantly enhance the performance of the qr decomposition for a square matrix. Instead of operating on a square matrix, I decided to make my input matrix A more rectangular in order to see a greater performance increase than I have been seeing to test how well omp for loop optimizations affects non-square matrices. 











Changing Matrix input size

When we make m large and n smaller, we see more significant performance increase when using openmp then we had when using a square matrix. For this section we will only compare the results of the serial version of the qr factorization with the second implementation of OMP ( with omp loops for the reflect and  VecScale.) The results can be seen below. 

Serial
./qr-omp 1000 500
  QRFactor   0.786514 s	  0.530 GF/s
    dgeqrf   0.081827 s	  5.092 GF/s
max difference 3.197442e-14, tau 6.661338e-16

 ./qr-omp 10000 500
  QRFactor   9.106410 s	  0.540 GF/s
    dgeqrf   0.976660 s	  5.034 GF/s
max difference 2.866596e-13, tau 6.661338e-16

 ./qr-omp 20000 500
  QRFactor  18.364245 s	  0.540 GF/s
    dgeqrf   2.038289 s	  4.865 GF/s
max difference 8.242296e-13, tau 6.661338e-16



Reflect and VecScale OMP implementation -2 threads

./qr-omp_2 1000 500
  QRFactor   0.327670 s	  1.272 GF/s
    dgeqrf   0.093328 s	  4.465 GF/s
max difference 3.197442e-14, tau 6.661338e-16

./qr-omp_2 10000 500
  QRFactor   3.191995 s	  1.540 GF/s
    dgeqrf   0.975683 s	  5.039 GF/s
max difference 2.866596e-13, tau 6.661338e-16

./qr-omp_2 20000 500
  QRFactor   6.410528 s	  1.547 GF/s
    dgeqrf   2.059463 s	  4.815 GF/s
max difference 8.242296e-13, tau 6.661338e-16





Reflect and VecScale OMP implementation - 4 threads
./qr-omp_2 1000 500
  QRFactor   0.351599 s	  1.185 GF/s
    dgeqrf   0.086918 s	  4.794 GF/s
max difference 3.197442e-14, tau 6.661338e-16
 ./qr-omp_2 10000 500
  QRFactor   3.226675 s	  1.524 GF/s
    dgeqrf   0.982304 s	  5.005 GF/s
max difference 2.866596e-13, tau 6.661338e-16

 ./qr-omp_2 20000 500
  QRFactor   6.407659 s	  1.548 GF/s
    dgeqrf   2.076408 s	  4.776 GF/s
max difference 8.242296e-13, tau 6.661338e-16

As we can see, our implementation of OMP for loops for both Reflect and VecScale significantly increase our performance for large m and small n matrices. The performance increase here is nearly a factor of 3 increase over the serial implementation given to us to start with. While this is true, there is practically no increase in perfomance between 2 and 4 threads, most of this is probably due to only having two physical processors as well as poor memory managment. For small matrices our performance actually decreases when  initializing 4 threads instead of 2.


Although not shown, we see approximately the same performance increase between the two implementations of QR decomp when we switch the instances of m and n so that n is now large and m is small. 
	-This kinda confuses me that we see similar increase in both situations since both fors loop over n and not m, I would think that a higher performance increase would be seen for matrices where n is large and not m.
