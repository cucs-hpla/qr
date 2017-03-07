#matrix size input parameters m & n
#print qr: a=3, no print qr a =0
m = 10
n = 4
a = 4

OMP_NUM_THREADS=4

CC = gcc
LINKER = $(CC)
CFLAGS = -Wall -std=c99 -Og -fopenmp #-g -fsanitize=address -fstack-protector 
LDFLAGS = -L/usr/local/opt/openblas/lib -lblas -lm

UTIL := qr-omp-ts.o

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
	OMP_NUM_THREADS=4 ./qr.x $(m) $(n) $(a)

clean:
	rm -f *.o *~ core *.x

